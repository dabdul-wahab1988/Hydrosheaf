"""Run E2b MODPATH-conditioned geochemical replay validation.

E2 validates that real MODPATH endpoint/pathline outputs can be converted into
Hydrosheaf graph priors. E2b adds a controlled geochemical layer on the same
real graph: known transport and reaction signals are injected along MODPATH
edges, then Hydrosheaf edge fitting is asked to recover them.
"""
from __future__ import annotations

import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.config import Config  # noqa: E402
from hydrosheaf.inference.edge_fit import fit_edge  # noqa: E402
from hydrosheaf.models.reactions import build_reaction_dictionary  # noqa: E402

BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
MODPATH_RESULT_DIR = BENCHMARK_ROOT / "external" / "modpath" / "results"
RESULT_DIR = MODPATH_RESULT_DIR
FIGURE_DIR = BENCHMARK_ROOT / "figures"

GRAPH_PRIORS = MODPATH_RESULT_DIR / "modpath_graph_priors.csv"
SOURCE_DOI = "10.5066/F7J102FK"
ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]
ACTIVE_REACTIONS = ["calcite", "gypsum", "NO3src", "denit"]
DETECTION_THRESHOLD = 0.01


def _stable_float(text: str, modulus: int, scale: float) -> float:
    value = sum((idx + 1) * ord(char) for idx, char in enumerate(text)) % modulus
    return float(value) / float(max(modulus - 1, 1)) * scale


def _config() -> Config:
    cfg = Config(
        phreeqc_enabled=False,
        gibbs_enabled=False,
        isotope_enabled=False,
        lambda_l1=0.0,
        transport_models_enabled=["evap"],
        active_minerals=["calcite", "gypsum"],
    )
    cfg.weights = [1.0] * len(ION_ORDER)
    cfg.conservative_weights = [0.01] * len(ION_ORDER)
    cfg.conservative_weights[ION_ORDER.index("Cl")] = 1.0
    return cfg


def _reaction_vectors(config: Config) -> Dict[str, np.ndarray]:
    matrix, labels, _ = build_reaction_dictionary(config)
    return {label: np.asarray(vector, dtype=float) for label, vector in zip(labels, matrix)}


def _source_vector(node_id: str) -> np.ndarray:
    return np.asarray(
        [
            0.85 + _stable_float(node_id + "ca", 19, 0.35),
            0.25 + _stable_float(node_id + "mg", 17, 0.20),
            0.90 + _stable_float(node_id + "na", 23, 0.50),
            0.035 + _stable_float(node_id + "k", 11, 0.025),
            2.20 + _stable_float(node_id + "hco3", 29, 0.90),
            1.10 + _stable_float(node_id + "cl", 31, 0.45),
            0.25 + _stable_float(node_id + "so4", 13, 0.22),
            0.25 + _stable_float(node_id + "no3", 37, 0.35),
            0.008 + _stable_float(node_id + "f", 7, 0.010),
            0.001,
            0.001,
        ],
        dtype=float,
    )


def _truth_for_edge(edge: Mapping[str, object]) -> Dict[str, float]:
    travel_time = float(edge["travel_time_mean"])
    particle_count = float(edge["particle_count"])
    source = str(edge["source_node"])
    target = str(edge["target_node"])
    age_factor = min(max(math.log10(max(travel_time, 1.0)) / 4.0, 0.0), 1.0)
    fast_edge = travel_time < 700.0
    sink_factor = _stable_float(target, 41, 1.0)
    source_factor = _stable_float(source, 43, 1.0)
    return {
        "gamma": 1.0 + 0.015 + 0.045 * age_factor,
        "calcite": 0.04 + 0.12 * age_factor + 0.02 * source_factor,
        "gypsum": 0.02 + 0.08 * sink_factor,
        "NO3src": (0.05 + 0.10 * sink_factor) if fast_edge else 0.0,
        "denit": 0.0 if fast_edge else min(0.22, 0.035 + 0.000035 * travel_time + 0.002 * particle_count),
    }


def _downstream_vector(source: np.ndarray, truth: Mapping[str, float], reaction_vectors: Mapping[str, np.ndarray]) -> np.ndarray:
    out = source * float(truth["gamma"])
    for label in ACTIVE_REACTIONS:
        out = out + float(truth.get(label, 0.0)) * reaction_vectors[label]
    return np.maximum(out, 1.0e-8)


def _row_to_sample(site_id: str, vector: Sequence[float]) -> Dict[str, object]:
    row: Dict[str, object] = {"site_id": site_id}
    row.update({ion: float(value) for ion, value in zip(ION_ORDER, vector)})
    return row


def _fit_replay(edge: Mapping[str, object], config: Config, reaction_vectors: Mapping[str, np.ndarray]) -> Dict[str, object]:
    source_id = str(edge["source_node"])
    target_id = str(edge["target_node"])
    edge_id = str(edge["edge_id"])
    x_u = _source_vector(source_id)
    truth = _truth_for_edge(edge)
    x_v = _downstream_vector(x_u, truth, reaction_vectors)
    result = fit_edge(
        x_u.tolist(),
        x_v.tolist(),
        config,
        edge_id=edge_id,
        u=source_id,
        v=target_id,
        obs_u=_row_to_sample(source_id, x_u),
        obs_v=_row_to_sample(target_id, x_v),
        residence_time_days=float(edge["travel_time_mean"]),
    )
    estimate_by_label = {label: float(value) for label, value in zip(result.z_labels, result.z_extents)}
    row: Dict[str, object] = {
        "edge_id": edge_id,
        "source_node": source_id,
        "target_node": target_id,
        "particle_count": int(edge["particle_count"]),
        "travel_time_mean": float(edge["travel_time_mean"]),
        "gamma_true": float(truth["gamma"]),
        "gamma_estimated": float(result.gamma or float("nan")),
        "gamma_abs_error": abs(float(result.gamma or float("nan")) - float(truth["gamma"])),
        "objective_score": float(result.objective_score),
        "transport_residual_norm": float(result.transport_residual_norm),
        "anomaly_norm": float(result.anomaly_norm),
    }
    active_truth: List[str] = []
    active_estimated: List[str] = []
    extent_errors: List[float] = []
    for label in ACTIVE_REACTIONS:
        true_value = float(truth.get(label, 0.0))
        estimated = float(estimate_by_label.get(label, 0.0))
        row[f"{label}_true"] = true_value
        row[f"{label}_estimated"] = estimated
        row[f"{label}_abs_error"] = abs(estimated - true_value)
        extent_errors.append(abs(estimated - true_value))
        if abs(true_value) >= DETECTION_THRESHOLD:
            active_truth.append(label)
        if abs(estimated) >= DETECTION_THRESHOLD:
            active_estimated.append(label)
    truth_set = set(active_truth)
    estimated_set = set(active_estimated)
    row["active_truth"] = ";".join(active_truth)
    row["active_estimated"] = ";".join(active_estimated)
    row["reaction_tp"] = len(truth_set & estimated_set)
    row["reaction_fp"] = len(estimated_set - truth_set)
    row["reaction_fn"] = len(truth_set - estimated_set)
    row["mean_extent_abs_error"] = float(np.mean(extent_errors))
    row["max_extent_abs_error"] = float(np.max(extent_errors))
    return row


def _summary(results: pd.DataFrame) -> Dict[str, float]:
    tp = int(results["reaction_tp"].sum())
    fp = int(results["reaction_fp"].sum())
    fn = int(results["reaction_fn"].sum())
    return {
        "n_edges": int(len(results)),
        "gamma_mae": float(results["gamma_abs_error"].mean()),
        "mean_extent_mae": float(results["mean_extent_abs_error"].mean()),
        "max_extent_mae": float(results["max_extent_abs_error"].mean()),
        "reaction_tp": tp,
        "reaction_fp": fp,
        "reaction_fn": fn,
        "reaction_precision": float(tp / max(tp + fp, 1)),
        "reaction_recall": float(tp / max(tp + fn, 1)),
        "reaction_f1": float((2 * tp) / max(2 * tp + fp + fn, 1)),
        "median_objective_score": float(results["objective_score"].median()),
    }


def _write_figure(results: pd.DataFrame, summary: Mapping[str, float]) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = FIGURE_DIR / "figure_s2b_modpath_conditioned_geochemical_replay.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.scatter(results["gamma_true"], results["gamma_estimated"], s=24, color="#3b6f8f", alpha=0.75)
    lo = min(results["gamma_true"].min(), results["gamma_estimated"].min())
    hi = max(results["gamma_true"].max(), results["gamma_estimated"].max())
    ax.plot([lo, hi], [lo, hi], color="#333333", linewidth=1)
    ax.set_xlabel("injected transport gamma")
    ax.set_ylabel("recovered transport gamma")
    ax.set_title("Transport recovery on MODPATH edges")

    ax = axes[1]
    labels = ACTIVE_REACTIONS
    true_means = [results[f"{label}_true"].mean() for label in labels]
    est_means = [results[f"{label}_estimated"].mean() for label in labels]
    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x - width / 2, true_means, width, label="injected", color="#627c62")
    ax.bar(x + width / 2, est_means, width, label="recovered", color="#c47f3f")
    ax.set_xticks(x, labels=labels, rotation=30, ha="right")
    ax.set_ylabel("mean reaction extent")
    ax.set_title("Reaction replay recovery")
    ax.legend(fontsize=8)
    ax.text(
        0.98,
        0.95,
        f"F1={summary['reaction_f1']:.2f}\nextent MAE={summary['mean_extent_mae']:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )

    fig.suptitle("Figure S2B. MODPATH-conditioned geochemical replay", fontsize=12)
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)
    return figure_path


def _update_tables(summary: Mapping[str, float], report_path: Path) -> None:
    metric = (
        f"n={int(summary['n_edges'])}; gamma MAE={summary['gamma_mae']:.3f}; "
        f"reaction F1={summary['reaction_f1']:.2f}; extent MAE={summary['mean_extent_mae']:.3f}"
    )
    table4 = BENCHMARK_ROOT / "tables" / "table4_validation_design_and_results.csv"
    if table4.exists():
        df = pd.read_csv(table4)
        benchmark = "MODPATH-conditioned geochemical replay"
        row = {
            "benchmark": benchmark,
            "data_source": f"USGS Savage MODPATH graph with injected known-truth chemistry, DOI {SOURCE_DOI}",
            "target_variable": "transport gamma and reaction identity/extent on real MODPATH topology",
            "performance_metric": metric,
            "expected_evidence": "Hydrosheaf recovers known geochemical transformations when conditioned on real MODPATH graph priors",
            "key_reference": f"Harte USGS data release {SOURCE_DOI}; Hydrosheaf controlled replay",
            "m2_status": "completed",
            "run_note": "semi-synthetic chemistry replay on real physical topology",
            "notes": f"Report: {report_path.name}",
        }
        if benchmark in set(df["benchmark"]):
            mask = df["benchmark"] == benchmark
            for key, value in row.items():
                df.loc[mask, key] = value
        else:
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
        df.to_csv(table4, index=False)

    workplan = BENCHMARK_ROOT / "tables" / "external_validation_workplan.csv"
    if workplan.exists():
        df = pd.read_csv(workplan)
        tier = "E2b_modpath_geochemical_replay"
        row = {
            "validation_tier": tier,
            "section_or_figure": "3.4 and Figure S2B",
            "dataset": "USGS Savage MODPATH graph with injected known-truth chemistry",
            "doi": SOURCE_DOI,
            "access_url": "local MODPATH archive plus controlled replay script",
            "task": "Inject known transport/reaction signals along real MODPATH edges and test Hydrosheaf recovery",
            "outputs": "external/modpath/results/modpath_geochemical_replay.csv; figures/figure_s2b_modpath_conditioned_geochemical_replay.png",
            "status": "completed",
            "notes": metric,
        }
        if tier in set(df["validation_tier"]):
            mask = df["validation_tier"] == tier
            for key, value in row.items():
                df.loc[mask, key] = value
        else:
            df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
        df.to_csv(workplan, index=False)


def main() -> int:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    if not GRAPH_PRIORS.exists():
        raise FileNotFoundError(f"Run E2 first; missing {GRAPH_PRIORS}.")

    config = _config()
    reaction_vectors = _reaction_vectors(config)
    graph = pd.read_csv(GRAPH_PRIORS).sort_values(["particle_count", "travel_time_mean"], ascending=[False, True])
    graph = graph.head(80).reset_index(drop=True)
    results = pd.DataFrame([_fit_replay(row, config, reaction_vectors) for row in graph.to_dict("records")])
    summary = _summary(results)
    figure_path = _write_figure(results, summary)

    results_path = RESULT_DIR / "modpath_geochemical_replay.csv"
    summary_path = RESULT_DIR / "modpath_geochemical_replay_summary.csv"
    report_path = RESULT_DIR / "e2b_modpath_geochemical_replay_report.md"
    manifest_path = RESULT_DIR / "e2b_modpath_geochemical_replay_manifest.json"
    results.to_csv(results_path, index=False)
    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    report_lines = [
        "# E2b MODPATH-Conditioned Geochemical Replay Report",
        "",
        f"Run timestamp UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"Physical topology source: USGS Savage MODPATH archive, DOI `{SOURCE_DOI}`.",
        "",
        "## Scope",
        "",
        "This is a semi-synthetic integration validation. The graph topology and travel-time priors are real MODPATH outputs from E2. The geochemical observations are controlled known-truth injections placed on those real edges because the Savage archive does not provide paired chemistry/tracer samples for the same graph nodes.",
        "",
        "## Outputs",
        "",
        f"- Replay results: `{results_path}`",
        f"- Replay summary: `{summary_path}`",
        f"- Figure S2B: `{figure_path}`",
        "",
        "## Summary",
        "",
        pd.DataFrame([summary]).to_markdown(index=False),
        "",
        "## Interpretation",
        "",
        "E2b closes the integration gap left by E2: it tests whether Hydrosheaf can recover known geochemical transport and reaction signals when conditioned on a real MODPATH-derived physical graph. It does not replace a future fully paired physical-geochemical field dataset, but it prevents the M2 paper from implying that the MODPATH archive alone validates chemistry.",
        "",
    ]
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    manifest_path.write_text(
        json.dumps(
            {
                "validation_tier": "E2b_modpath_geochemical_replay",
                "source_doi": SOURCE_DOI,
                "graph_priors": str(GRAPH_PRIORS),
                "n_edges_replayed": int(summary["n_edges"]),
                "active_reactions": ACTIVE_REACTIONS,
                "outputs": {
                    "results": str(results_path),
                    "summary": str(summary_path),
                    "figure": str(figure_path),
                    "report": str(report_path),
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    _update_tables(summary, report_path)
    print(f"E2b completed with {int(summary['n_edges'])} MODPATH-conditioned replay edges.")
    print(f"Reaction F1={summary['reaction_f1']:.3f}; extent MAE={summary['mean_extent_mae']:.4f}")
    print(f"Figure S2B written to {figure_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
