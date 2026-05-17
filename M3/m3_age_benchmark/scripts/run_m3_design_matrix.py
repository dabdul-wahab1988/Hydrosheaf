"""Run the Phase-2 M3 experimental design matrix."""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parents[0]
CONFIG_DIR = BENCHMARK_DIR / "configs"
RESULT_DIR = BENCHMARK_DIR / "results"
DOCS_DIR = BENCHMARK_DIR / "docs"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_m3_usgs_benchmark as usgs


DEFAULT_CONFIG = CONFIG_DIR / "design_matrix.yaml"
DEFAULT_OUTPUT = RESULT_DIR / "m3_design_matrix_results.csv"


def companion_output_paths(output: Path) -> dict[str, Path]:
    if output == DEFAULT_OUTPUT:
        return {
            "summary": RESULT_DIR / "m3_design_matrix_summary.csv",
            "pairwise": RESULT_DIR / "m3_design_matrix_pairwise_deltas.csv",
            "manifest": RESULT_DIR / "m3_design_matrix_manifest.json",
            "qa": DOCS_DIR / "m3_design_matrix_qa.md",
        }
    stem = output.stem
    return {
        "summary": output.with_name(f"{stem}_summary.csv"),
        "pairwise": output.with_name(f"{stem}_pairwise_deltas.csv"),
        "manifest": output.with_name(f"{stem}_manifest.json"),
        "qa": DOCS_DIR / f"{stem}_qa.md",
    }


def load_design_matrix(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or "scenarios" not in config:
        raise ValueError("Design matrix config must contain a top-level 'scenarios' list.")
    for scenario in config["scenarios"]:
        if "scenario_id" not in scenario:
            raise ValueError("Every design-matrix scenario must define scenario_id.")
    return config


def resolve_age_steps(
    config: dict,
    age_steps: int | None,
    *,
    max_rows: int | None,
    allow_coarse_full_grid: bool = False,
) -> int:
    resolved = int(age_steps or config.get("defaults", {}).get("age_steps") or usgs.M3_DEFAULT_AGE_STEPS)
    if max_rows is None and resolved < usgs.M3_DEFAULT_AGE_STEPS and not allow_coarse_full_grid:
        raise ValueError(
            "Full design-matrix runs require age_steps >= "
            f"{usgs.M3_DEFAULT_AGE_STEPS}; got {resolved}. "
            "Use --allow-coarse-full-grid only for explicit runtime diagnostics."
        )
    return resolved


def _with_reference_age(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_m3_ref_age"] = out["reference_age_years"].map(usgs._parse_age)
    out["_m3_age_class"] = out["_m3_ref_age"].map(usgs._age_class)
    return out[out["_m3_ref_age"].notna()].copy()


def select_rows(df: pd.DataFrame, max_rows: int | None) -> pd.DataFrame:
    df = _with_reference_age(df)
    if max_rows is None or max_rows <= 0 or len(df) <= max_rows:
        return df.drop(columns=["_m3_ref_age", "_m3_age_class"], errors="ignore")

    per_class = max(1, math.ceil(max_rows / max(1, df["_m3_age_class"].nunique())))
    selected = (
        df.sort_values(["_m3_age_class", "site_id"])
        .groupby("_m3_age_class", group_keys=False)
        .head(per_class)
        .head(max_rows)
    )
    return selected.drop(columns=["_m3_ref_age", "_m3_age_class"], errors="ignore")


def run_design_matrix(
    df: pd.DataFrame,
    config: dict,
    *,
    scenario_ids: set[str] | None = None,
    max_rows: int | None = 80,
    age_steps: int | None = None,
    skip_selection: bool = False,
) -> pd.DataFrame:
    defaults = dict(config.get("defaults") or {})
    base_age_steps = int(age_steps or defaults.get("age_steps") or usgs.M3_DEFAULT_AGE_STEPS)
    scenarios = list(config["scenarios"])
    if scenario_ids:
        scenarios = [scenario for scenario in scenarios if scenario["scenario_id"] in scenario_ids]
    if skip_selection:
        scenarios = [scenario for scenario in scenarios if scenario.get("lpm_strategy") != "selection"]
    if not scenarios:
        raise ValueError("No design-matrix scenarios selected.")

    df_run = select_rows(df, max_rows)
    rows: list[dict] = []
    for scenario in scenarios:
        scenario_id = scenario["scenario_id"]
        model_strategy = scenario.get("lpm_strategy", "reported")
        print(f"Running scenario {scenario_id} on {len(df_run)} rows...", flush=True)
        for i, (_, row) in enumerate(df_run.iterrows(), start=1):
            if i == 1 or i % 10 == 0:
                print(f"  {scenario_id}: row {i}/{len(df_run)}", flush=True)
            result = usgs.fit_benchmark_row(
                row,
                age_steps=base_age_steps,
                model_strategy=model_strategy,
                factors=scenario,
            )
            result.update(
                {
                    "scenario_id": scenario_id,
                    "scenario_label": scenario.get("label", scenario_id),
                    "scenario_purpose": scenario.get("purpose", ""),
                    "design_age_steps": base_age_steps,
                }
            )
            rows.append(result)
    return pd.DataFrame(rows)


def summarize_results(results: pd.DataFrame) -> pd.DataFrame:
    totals = (
        results.groupby("scenario_id")
        .agg(total_rows=("site_id", "count"))
        .reset_index()
    )
    metric_rows = results[results["ref_age"].notna() & results["est_age_multi"].notna()].copy()
    if metric_rows.empty:
        return totals
    metrics = (
        metric_rows.groupby("scenario_id")
        .agg(
            metric_rows=("site_id", "count"),
            finite_estimates=("est_age_multi", "count"),
            median_abs_log10_error=("log10_error", "median"),
            log10_rmse=("log10_error", lambda s: float((s.dropna().pow(2).mean()) ** 0.5) if s.dropna().size else float("nan")),
            within_factor_2=("within_factor_2", "mean"),
            within_factor_10=("within_factor_10", "mean"),
            calibrated_he4_rows=("he4_calibrated", "sum"),
        )
        .reset_index()
    )
    return pd.merge(totals, metrics, on="scenario_id", how="left")


def summarize_pairwise_deltas(
    results: pd.DataFrame,
    *,
    baseline_scenario: str = "parity_reported_corrected",
) -> pd.DataFrame:
    """Summarize paired scenario effects against the baseline scenario."""
    baseline = results[results["scenario_id"] == baseline_scenario][
        ["site_id", "log10_error", "est_age_multi", "within_factor_2", "within_factor_10"]
    ].rename(
        columns={
            "log10_error": "baseline_log10_error",
            "est_age_multi": "baseline_est_age_multi",
            "within_factor_2": "baseline_within_factor_2",
            "within_factor_10": "baseline_within_factor_10",
        }
    )
    if baseline.empty:
        return pd.DataFrame()

    comparisons = results[results["scenario_id"] != baseline_scenario].merge(baseline, on="site_id", how="inner")
    comparisons = comparisons[
        comparisons["log10_error"].notna() & comparisons["baseline_log10_error"].notna()
    ].copy()
    if comparisons.empty:
        return pd.DataFrame()

    comparisons["delta_log10_error"] = comparisons["log10_error"] - comparisons["baseline_log10_error"]
    comparisons["improved_vs_baseline"] = comparisons["delta_log10_error"] < 0
    comparisons["lost_factor_2"] = comparisons["baseline_within_factor_2"].astype(bool) & ~comparisons["within_factor_2"].astype(bool)
    comparisons["gained_factor_2"] = ~comparisons["baseline_within_factor_2"].astype(bool) & comparisons["within_factor_2"].astype(bool)

    return (
        comparisons.groupby("scenario_id")
        .agg(
            paired_rows=("site_id", "count"),
            median_delta_log10_error=("delta_log10_error", "median"),
            mean_delta_log10_error=("delta_log10_error", "mean"),
            improved_fraction=("improved_vs_baseline", "mean"),
            gained_factor_2_rows=("gained_factor_2", "sum"),
            lost_factor_2_rows=("lost_factor_2", "sum"),
        )
        .reset_index()
    )


def write_summary(
    results: pd.DataFrame,
    summary: pd.DataFrame,
    pairwise: pd.DataFrame,
    output: Path,
    config_path: Path,
    max_rows: int | None,
    qa_path: Path,
) -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    lines = [
        "# M3 Phase-2 Design Matrix QA",
        "",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Config: `{config_path}`",
        f"- Output: `{output}`",
        f"- Output rows: {len(results)}",
        f"- Unique scenarios: {results['scenario_id'].nunique() if not results.empty else 0}",
        f"- Row limit: {'full' if max_rows is None else max_rows}",
        "",
        "## Scenario Metrics",
        "",
    ]
    if summary.empty:
        lines.append("- No finite metric rows.")
    else:
        lines.append("```text")
        lines.append(summary.to_string(index=False))
        lines.append("```")
    lines.append("")
    lines.append("## Paired Effects Versus `parity_reported_corrected`")
    lines.append("")
    if pairwise.empty:
        lines.append("- No paired effect rows.")
    else:
        lines.append("```text")
        lines.append(pairwise.to_string(index=False))
        lines.append("```")
    lines.append("")
    qa_path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Phase-2 M3 design-matrix scenarios.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-rows", type=int, default=80, help="Deterministic stratified row limit. Use --full for all rows.")
    parser.add_argument("--full", action="store_true", help="Run all eligible rows.")
    parser.add_argument("--age-steps", type=int, default=None)
    parser.add_argument(
        "--allow-coarse-full-grid",
        action="store_true",
        help="Allow a full run to use age_steps below the canonical M3 default. Intended only for runtime diagnostics.",
    )
    parser.add_argument("--scenario", action="append", help="Scenario ID to run. May be repeated.")
    parser.add_argument("--skip-selection", action="store_true", help="Skip Hydrosheaf model-selection scenarios.")
    args = parser.parse_args(argv)

    config = load_design_matrix(args.config)
    df = usgs.load_usgs_national_dataset()
    max_rows = None if args.full else args.max_rows
    resolved_age_steps = resolve_age_steps(
        config,
        args.age_steps,
        max_rows=max_rows,
        allow_coarse_full_grid=args.allow_coarse_full_grid,
    )
    results = run_design_matrix(
        df,
        config,
        scenario_ids=set(args.scenario or []) or None,
        max_rows=max_rows,
        age_steps=resolved_age_steps,
        skip_selection=args.skip_selection,
    )
    paths = companion_output_paths(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, index=False)
    primary_scenario = next(
        (
            scenario
            for scenario in ("screened_dgm_gases", "parity_reported_corrected")
            if scenario in set(results["scenario_id"])
        ),
        str(results["scenario_id"].iloc[0]) if not results.empty else "",
    )
    primary_output = RESULT_DIR / "m3_usgs_benchmark_results.csv"
    if primary_scenario:
        results[results["scenario_id"] == primary_scenario].to_csv(primary_output, index=False)
    summary = summarize_results(results)
    summary.to_csv(paths["summary"], index=False)
    pairwise = summarize_pairwise_deltas(results)
    pairwise.to_csv(paths["pairwise"], index=False)
    write_summary(results, summary, pairwise, args.output, args.config, max_rows, paths["qa"])
    manifest = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(args.config),
        "output": str(args.output),
        "summary": str(paths["summary"]),
        "pairwise": str(paths["pairwise"]),
        "qa": str(paths["qa"]),
        "primary_pointwise_output": str(primary_output) if primary_scenario else "",
        "primary_pointwise_scenario": primary_scenario,
        "max_rows": max_rows,
        "age_steps": resolved_age_steps,
        "scenario_ids": sorted(results["scenario_id"].unique().tolist()),
        "n_output_rows": int(len(results)),
        "pairwise_baseline": "parity_reported_corrected",
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote {len(results)} rows to {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
