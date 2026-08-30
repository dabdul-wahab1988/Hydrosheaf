from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hydrosheaf.graph.build import infer_edges_from_coordinates
from hydrosheaf.graph.parameter_regularize import (
    parameter_graph_regularize,
    build_aicc_edge_weights,
    reconstruct_age_from_regularized_tau,
)
from hydrosheaf.nuclear import audit_graph_age_coherence
from hydrosheaf.validation import regression_metrics

import run_m3_usgs_benchmark as usgs


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
DOCS_DIR = BENCHMARK_DIR / "docs"

DEFAULT_POINTWISE_RESULTS = RESULTS_DIR / "m3_design_matrix_results.csv"
DEFAULT_SCENARIO = "tracerlpm_strict_parity"

PRIOR_STRENGTHS = {
    "none": 0.0,
    "weak": 0.25,
    "medium": 0.55,
    "strong": 0.85,
}


def _finite_float(value) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius_km = 6371.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)
    a = (
        math.sin(d_phi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2.0) ** 2
    )
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return radius_km * c


def _preferred_pointwise_results() -> Path:
    for candidate in (
        RESULTS_DIR / "m3_design_matrix_results.csv",
        RESULTS_DIR / "m3_usgs_benchmark_results.csv",
    ):
        if candidate.exists():
            return candidate
    return DEFAULT_POINTWISE_RESULTS


def log10_rmse(reference: Iterable[float], estimate: Iterable[float]) -> float:
    metrics = regression_metrics(list(reference), list(estimate), log10=True)
    return float(metrics["rmse"])


def median_abs_log10_error(reference: Iterable[float], estimate: Iterable[float]) -> float:
    values = []
    for ref, est in zip(reference, estimate):
        ref_value = _finite_float(ref)
        est_value = _finite_float(est)
        if not math.isfinite(ref_value) or not math.isfinite(est_value) or ref_value <= 0 or est_value <= 0:
            continue
        values.append(abs(math.log10(est_value) - math.log10(ref_value)))
    return float(pd.Series(values).median()) if values else float("nan")


def graph_regularize(
    single_ages: Mapping[str, float],
    edges: list[tuple[str, str]],
    strength: float,
    *,
    iterations: int = 8,
    min_increment_years: float = 0.0,
) -> dict[str, float]:
    ages = {node: float(age) for node, age in single_ages.items()}
    if strength <= 0:
        return ages
    alpha = min(max(float(strength), 0.0), 1.0)
    for _ in range(iterations):
        for upstream, downstream in edges:
            if upstream not in ages or downstream not in ages:
                continue
            required = ages[upstream] + min_increment_years
            if ages[downstream] < required:
                ages[downstream] = (1.0 - alpha) * ages[downstream] + alpha * required
    return ages


def load_pointwise_nodes(path: Path, scenario_id: str) -> pd.DataFrame:
    df_results = pd.read_csv(path, low_memory=False)
    if "scenario_id" in df_results.columns:
        available = set(df_results["scenario_id"].dropna().astype(str))
        if scenario_id not in available:
            raise ValueError(
                f"Scenario {scenario_id!r} is absent from {path}; available scenarios: "
                + ", ".join(sorted(available))
            )
        df_results = df_results[df_results["scenario_id"] == scenario_id].copy()
    if "scenario_withdrawn" in df_results.columns:
        withdrawn = df_results["scenario_withdrawn"].astype(str).str.lower().eq("true")
        df_results = df_results[~withdrawn].copy()
    df_results = df_results[
        df_results["site_id"].notna()
        & pd.to_numeric(df_results["ref_age"], errors="coerce").notna()
        & pd.to_numeric(df_results["est_age_multi"], errors="coerce").notna()
    ].copy()
    source_columns = [
        "site_id",
        *[
            column
            for column in ("StudyUnit", "lat", "lon", "depth_m")
            if column not in df_results.columns
        ],
    ]
    df_source = usgs.load_benchmark_dataset()[source_columns].copy()
    merged = pd.merge(df_results, df_source, on="site_id", how="left")
    merged["ref_age"] = pd.to_numeric(merged["ref_age"], errors="coerce")
    merged["est_age_multi"] = pd.to_numeric(merged["est_age_multi"], errors="coerce")
    merged["lat"] = pd.to_numeric(merged["lat"], errors="coerce")
    merged["lon"] = pd.to_numeric(merged["lon"], errors="coerce")
    merged["depth_m"] = pd.to_numeric(merged["depth_m"], errors="coerce")
    nodes = merged[
        merged["site_id"].notna()
        & merged["ref_age"].notna()
        & merged["est_age_multi"].notna()
        & merged["lat"].notna()
        & merged["lon"].notna()
        & merged["StudyUnit"].notna()
    ].copy()
    if nodes.empty:
        raise ValueError(
            f"No reportable coordinate-bearing nodes remain for {scenario_id!r} in {path}."
        )
    return nodes


def _samples_for_edges(df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        sample = {
            "site_id": row["site_id"],
            "lat": float(row["lat"]),
            "lon": float(row["lon"]),
        }
        depth = _finite_float(row.get("depth_m"))
        if math.isfinite(depth):
            sample["elevation"] = -depth
        rows.append(sample)
    return rows


def _depth_lookup(df: pd.DataFrame) -> dict[str, float]:
    return {
        str(row["site_id"]): _finite_float(row["depth_m"])
        for _, row in df.iterrows()
    }


def _study_unit_lookup(df: pd.DataFrame) -> dict[str, str]:
    return {
        str(row["site_id"]): str(row["StudyUnit"])
        for _, row in df.iterrows()
    }


def build_graph_families(df: pd.DataFrame, *, min_unit_size: int = 5) -> tuple[dict[str, list[tuple[str, str]]], pd.DataFrame]:
    family_edges: dict[str, list[tuple[str, str]]] = {
        "coordinate_global": [],
        "study_unit_coordinate": [],
        "depth_constrained": [],
        "hydraulic_proxy_constrained": [],
        "parameter_smooth_aicc": [],
        "wrong_direction_negative_control": [],
        "randomized_negative_control": [],
    }
    edge_rows: list[dict[str, object]] = []
    depth_lookup = _depth_lookup(df)
    study_lookup = _study_unit_lookup(df)

    global_edges = infer_edges_from_coordinates(
        _samples_for_edges(df),
        max_neighbors=1,
        allow_uphill=False,
        elevation_key="elevation",
    )
    family_edges["coordinate_global"] = [(edge.u, edge.v) for edge in global_edges]

    rng = random.Random(42)
    unit_groups = [group.copy() for _, group in df.groupby("StudyUnit") if len(group) >= min_unit_size]
    study_unit_edges: list[tuple[str, str]] = []
    depth_edges: list[tuple[str, str]] = []
    hydraulic_proxy_edges: list[tuple[str, str]] = []
    randomized_edges: list[tuple[str, str]] = []

    for unit_df in unit_groups:
        unit_name = str(unit_df["StudyUnit"].iloc[0])
        unit_edges = infer_edges_from_coordinates(
            _samples_for_edges(unit_df),
            max_neighbors=1,
            allow_uphill=False,
            elevation_key="elevation",
        )
        unit_pairs = [(edge.u, edge.v) for edge in unit_edges]
        study_unit_edges.extend(unit_pairs)

        candidates = [(str(site_id), str(other_id)) for site_id in unit_df["site_id"] for other_id in unit_df["site_id"] if site_id != other_id]
        rng.shuffle(candidates)
        randomized_edges.extend(candidates[: len(unit_pairs)])

        for upstream, downstream in unit_pairs:
            up_depth = depth_lookup.get(upstream, float("nan"))
            down_depth = depth_lookup.get(downstream, float("nan"))
            up_row = unit_df[unit_df["site_id"] == upstream].iloc[0]
            down_row = unit_df[unit_df["site_id"] == downstream].iloc[0]
            distance = _haversine_km(float(up_row["lat"]), float(up_row["lon"]), float(down_row["lat"]), float(down_row["lon"]))
            depth_diff = down_depth - up_depth if math.isfinite(up_depth) and math.isfinite(down_depth) else float("nan")
            if math.isfinite(depth_diff) and depth_diff >= 5.0 and depth_diff <= 200.0:
                depth_edges.append((upstream, downstream))
                if distance <= 25.0:
                    hydraulic_proxy_edges.append((upstream, downstream))
            edge_rows.append(
                {
                    "family": "study_unit_coordinate",
                    "study_unit": unit_name,
                    "upstream": upstream,
                    "downstream": downstream,
                    "distance_km": distance,
                    "depth_diff_m": depth_diff,
                }
            )

    family_edges["study_unit_coordinate"] = study_unit_edges
    family_edges["depth_constrained"] = depth_edges
    family_edges["hydraulic_proxy_constrained"] = hydraulic_proxy_edges
    family_edges["parameter_smooth_aicc"] = hydraulic_proxy_edges
    family_edges["wrong_direction_negative_control"] = [(downstream, upstream) for upstream, downstream in hydraulic_proxy_edges]
    family_edges["randomized_negative_control"] = randomized_edges

    def _append_family_rows(family: str, edges: list[tuple[str, str]]) -> None:
        for upstream, downstream in edges:
            up_row = df[df["site_id"] == upstream].iloc[0]
            down_row = df[df["site_id"] == downstream].iloc[0]
            edge_rows.append(
                {
                    "family": family,
                    "study_unit": study_lookup.get(upstream, ""),
                    "upstream": upstream,
                    "downstream": downstream,
                    "distance_km": _haversine_km(float(up_row["lat"]), float(up_row["lon"]), float(down_row["lat"]), float(down_row["lon"])),
                    "depth_diff_m": _finite_float(down_row["depth_m"]) - _finite_float(up_row["depth_m"]),
                }
            )

    _append_family_rows("coordinate_global", family_edges["coordinate_global"])
    _append_family_rows("depth_constrained", family_edges["depth_constrained"])
    _append_family_rows("hydraulic_proxy_constrained", family_edges["hydraulic_proxy_constrained"])
    _append_family_rows("wrong_direction_negative_control", family_edges["wrong_direction_negative_control"])
    _append_family_rows("randomized_negative_control", family_edges["randomized_negative_control"])

    return family_edges, pd.DataFrame(edge_rows)


def summarize_unit_deltas(df: pd.DataFrame) -> tuple[float, int]:
    rows = []
    for _, unit_df in df.groupby("StudyUnit"):
        if len(unit_df) < 3:
            continue
        baseline_rmse = log10_rmse(unit_df["ref_age"], unit_df["baseline_age"])
        graph_rmse = log10_rmse(unit_df["ref_age"], unit_df["graph_age"])
        rows.append(graph_rmse - baseline_rmse)
    if not rows:
        return float("nan"), 0
    return float(pd.Series(rows).median()), int(sum(value < 0 for value in rows))


def benchmark_graph_families(df: pd.DataFrame, family_edges: dict[str, list[tuple[str, str]]]) -> pd.DataFrame:
    single_age_map = {str(row["site_id"]): float(row["est_age_multi"]) for _, row in df.iterrows()}
    ref_map = {str(row["site_id"]): float(row["ref_age"]) for _, row in df.iterrows()}
    baseline_rmse = log10_rmse(df["ref_age"], df["est_age_multi"])
    baseline_median = median_abs_log10_error(df["ref_age"], df["est_age_multi"])
    baseline_factor2 = float(((df["est_age_multi"] / df["ref_age"]).clip(lower=1 / 1e12, upper=1e12).apply(lambda x: 0.5 <= x <= 2.0)).mean())
    rows: list[dict[str, object]] = []
    # Pre-extract parameter-level data for parameter_smooth families
    has_params = "fit_mean_age_years" in df.columns and "fit_model_aiccs_json" in df.columns
    single_tau_map = {}
    node_aiccs = {}
    node_families = {}
    if has_params:
        for _, row in df.iterrows():
            sid = str(row["site_id"])
            tau = row["fit_mean_age_years"]
            if math.isfinite(tau) and tau > 0:
                single_tau_map[sid] = float(tau)
                node_families[sid] = str(row.get("multi_model", ""))
                try:
                    node_aiccs[sid] = json.loads(str(row.get("fit_model_aiccs_json", "{}")))
                except Exception:
                    node_aiccs[sid] = {}

    for family, edges in family_edges.items():
        control_type = "negative_control" if "negative_control" in family else "candidate_family"
        audit_before = audit_graph_age_coherence(edges, {site_id: {"age_years": age} for site_id, age in single_age_map.items()})
        for prior_label, strength in PRIOR_STRENGTHS.items():
            if family == "parameter_smooth_aicc" and has_params and single_tau_map:
                # Filter edges to nodes with valid taus
                valid_nodes = set(single_tau_map)
                param_edges = [(u, v) for u, v in edges if u in valid_nodes and v in valid_nodes]
                edge_weights = build_aicc_edge_weights(param_edges, node_aiccs)
                reg_taus = parameter_graph_regularize(
                    single_tau_map, param_edges, edge_weights, strength, min_increment_years=0.0
                )
                graph_age_map = {}
                for sid in single_age_map:
                    if sid in reg_taus:
                        graph_age_map[sid] = reconstruct_age_from_regularized_tau(
                            reg_taus[sid], node_families.get(sid, ""), {}
                        )
                    else:
                        graph_age_map[sid] = single_age_map[sid]
            else:
                graph_age_map = graph_regularize(single_age_map, edges, strength)
            graph_df = df.copy()
            graph_df["baseline_age"] = graph_df["site_id"].map(single_age_map)
            graph_df["graph_age"] = graph_df["site_id"].map(graph_age_map)
            audit_after = audit_graph_age_coherence(edges, {site_id: {"age_years": age} for site_id, age in graph_age_map.items()})
            median_unit_delta, units_improved = summarize_unit_deltas(graph_df)
            graph_factor2 = float(((graph_df["graph_age"] / graph_df["ref_age"]).clip(lower=1 / 1e12, upper=1e12).apply(lambda x: 0.5 <= x <= 2.0)).mean())
            rows.append(
                {
                    "graph_family": family,
                    "control_type": control_type,
                    "prior_strength": prior_label,
                    "prior_weight": strength,
                    "n_nodes": len(df),
                    "n_edges": len(edges),
                    "n_study_units": int(df["StudyUnit"].nunique()),
                    "rmse_single_log10": baseline_rmse,
                    "rmse_graph_log10": log10_rmse(graph_df["ref_age"], graph_df["graph_age"]),
                    "delta_rmse_graph_minus_single": log10_rmse(graph_df["ref_age"], graph_df["graph_age"]) - baseline_rmse,
                    "median_abs_log10_error_single": baseline_median,
                    "median_abs_log10_error_graph": median_abs_log10_error(graph_df["ref_age"], graph_df["graph_age"]),
                    "within_factor_2_single": baseline_factor2,
                    "within_factor_2_graph": graph_factor2,
                    "delta_within_factor_2_graph_minus_single": graph_factor2 - baseline_factor2,
                    "n_violations_before": audit_before["n_violations"],
                    "n_violations_after": audit_after["n_violations"],
                    "n_severe_violations_before": audit_before["n_severe_violations"],
                    "n_severe_violations_after": audit_after["n_severe_violations"],
                    "median_unit_delta_rmse_log10": median_unit_delta,
                    "study_units_improved_count": units_improved,
                    "improved_vs_single": bool(log10_rmse(graph_df["ref_age"], graph_df["graph_age"]) < baseline_rmse),
                }
            )
    return pd.DataFrame(rows)


def write_qa(results: pd.DataFrame, edge_rows: pd.DataFrame, pointwise_path: Path, scenario_id: str) -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    best = results.sort_values("delta_rmse_graph_minus_single").head(10)
    negative = results[results["control_type"] == "negative_control"].sort_values("delta_rmse_graph_minus_single").head(10)
    lines = [
        "# M3 Real-USGS Graph Benchmark QA",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Pointwise input: `{pointwise_path}`",
        f"- Scenario: `{scenario_id}`",
        f"- Benchmark rows: {len(results)}",
        f"- Edge rows: {len(edge_rows)}",
        "",
        "## Best Candidate Rows",
        "",
        "```text",
        best.to_string(index=False),
        "```",
        "",
        "## Negative Controls",
        "",
        "```text",
        negative.to_string(index=False),
        "```",
        "",
        "## Guardrail",
        "",
        "This exploratory graph benchmark uses only scalar ages that passed the reference-free reportability guard. StudyUnit is the aquifer grouping and depth is a hydraulic proxy because direct head data are unavailable. USGS ages are model-derived reference values, not field-observed truth. Candidate-family changes must be compared with wrong-direction and randomized controls and are described as reference-age agreement, not general predictive accuracy.",
    ]
    (DOCS_DIR / "m3_graph_benchmark_qa.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run the real-USGS graph benchmark for M3.")
    parser.add_argument("--pointwise-results", type=Path, default=None)
    parser.add_argument("--scenario", default=DEFAULT_SCENARIO)
    parser.add_argument("--min-unit-size", type=int, default=5)
    args = parser.parse_args(argv)

    pointwise_path = args.pointwise_results or _preferred_pointwise_results()
    df = load_pointwise_nodes(pointwise_path, args.scenario)
    family_edges, edge_rows = build_graph_families(df, min_unit_size=args.min_unit_size)
    results = benchmark_graph_families(df, family_edges)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    edge_rows.to_csv(RESULTS_DIR / "m3_real_usgs_graph_edges.csv", index=False)
    results.to_csv(RESULTS_DIR / "m3_real_usgs_graph_benchmark.csv", index=False)
    write_qa(results, edge_rows, pointwise_path, args.scenario)
    manifest = {
        "schema_version": "1.0",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "pointwise_results": str(pointwise_path),
        "pointwise_results_sha256": hashlib.sha256(pointwise_path.read_bytes()).hexdigest(),
        "scenario": args.scenario,
        "scope": "reportable coordinate-bearing scalar ages only",
        "reference_interpretation": "USGS model-derived age; not field-observed truth",
        "n_nodes": int(len(df)),
        "n_benchmark_rows": int(len(results)),
        "n_edge_rows": int(len(edge_rows)),
        "python_version": platform.python_version(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (RESULTS_DIR / "m3_real_usgs_graph_benchmark_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote graph benchmark rows to {RESULTS_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
