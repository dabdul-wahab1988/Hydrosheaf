"""Create data-driven M3 manuscript tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = BENCHMARK_ROOT / "tables" / "Manuscript_Ready"
RESULT_DIR = BENCHMARK_ROOT / "results"
TABLE_DIR.mkdir(parents=True, exist_ok=True)


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _design_summary() -> pd.DataFrame:
    frames = []
    for candidate in (
        RESULT_DIR / "m3_tracerlpm_strict_parity_full_summary.csv",
        RESULT_DIR / "m3_tracerlpm_parity_hier_oldwater_full_summary.csv",
        RESULT_DIR / "m3_tracerlpm_parity_agefractions_full_summary.csv",
        RESULT_DIR / "m3_hydrosheaf_selection_corrected_full_summary.csv",
        RESULT_DIR / "m3_parity_reported_corrected_full_summary.csv",
        RESULT_DIR / "m3_design_matrix_summary.csv",
    ):
        df = _read_csv(candidate)
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    summary = pd.concat(frames, ignore_index=True)
    if "scenario_id" in summary.columns:
        summary = summary.drop_duplicates("scenario_id", keep="first")
    return summary


def _mode_results() -> pd.DataFrame:
    frames = []
    for candidate in (
        RESULT_DIR / "m3_tracerlpm_strict_parity_full.csv",
        RESULT_DIR / "m3_tracerlpm_parity_hier_oldwater_full.csv",
        RESULT_DIR / "m3_tracerlpm_parity_agefractions_full.csv",
        RESULT_DIR / "m3_hydrosheaf_selection_corrected_full.csv",
        RESULT_DIR / "m3_parity_reported_corrected_full.csv",
        RESULT_DIR / "m3_design_matrix_results.csv",
    ):
        df = _read_csv(candidate)
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    if {"scenario_id", "site_id"}.issubset(out.columns):
        out = out.drop_duplicates(["scenario_id", "site_id"], keep="first")
    return out


def _primary_pointwise_results() -> pd.DataFrame:
    for candidate in (
        RESULT_DIR / "m3_tracerlpm_strict_parity_full.csv",
        RESULT_DIR / "m3_tracerlpm_parity_hier_oldwater_full.csv",
        RESULT_DIR / "m3_tracerlpm_parity_agefractions_full.csv",
        RESULT_DIR / "m3_design_matrix_results.csv",
        RESULT_DIR / "m3_phase4_screened_full_results.csv",
        RESULT_DIR / "m3_phase4_younggas_full_results.csv",
        RESULT_DIR / "m3_phase4_younggas_results.csv",
        RESULT_DIR / "m3_usgs_benchmark_results.csv",
    ):
        df = _read_csv(candidate)
        if df.empty:
            continue
        for preferred in ("tracerlpm_strict_parity", "tracerlpm_parity_hier_oldwater", "screened_dgm_gases", "parity_reported_corrected"):
            if "scenario_id" in df.columns and preferred in set(df["scenario_id"].dropna()):
                return df[df["scenario_id"] == preferred].copy()
        return df
    return pd.DataFrame()


def _log_rmse(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(np.sqrt(np.mean(values**2))) if len(values) else np.nan


def make_table1() -> None:
    rows = [
        ("Tritium (3H)", "12.32", "0-60", "TU", "bomb-pulse young-water tracer"),
        ("3H/3He", "12.32 parent", "0-50", "TU-equivalent", "closed-system ingrowth apparent age"),
        ("Carbon-14 (14C)", "5730", "1000-30000", "pmC", "old-groundwater radiocarbon constraint"),
        ("Helium-4 (4He)", "accumulation", "1000-100000", "ccSTP/g", "radiogenic accumulation screening constraint"),
        ("SF6", "stable", "0-50", "pptv", "atmospheric-equivalent gas tracer"),
        ("CFC-11/12/113", "stable", "0-60", "pptv", "atmospheric-equivalent gas tracers"),
    ]
    df = pd.DataFrame(
        rows,
        columns=["Tracer", "Decay or process scale", "Target age range (yr)", "Unit", "Benchmark role"],
    )
    df.to_csv(TABLE_DIR / "Manuscript_Table1_Nuclear_Tracer_Physics.csv", index=False)


def make_table2() -> None:
    summary = _design_summary()
    if summary.empty:
        return
    columns = [
        "scenario_id",
        "total_rows",
        "metric_rows",
        "median_abs_log10_error",
        "log10_rmse",
        "log10_r2",
        "within_factor_2",
        "within_factor_10",
        "calibrated_he4_rows",
    ]
    summary[[col for col in columns if col in summary.columns]].to_csv(
        TABLE_DIR / "Manuscript_Table2_Design_Matrix_Performance.csv",
        index=False,
    )


def make_table3() -> None:
    pairwise = _read_csv(RESULT_DIR / "m3_design_matrix_pairwise_deltas.csv")
    if pairwise.empty:
        return
    pairwise.to_csv(TABLE_DIR / "Manuscript_Table3_Paired_Ablation_Effects.csv", index=False)


def make_table4() -> None:
    graph = _read_csv(RESULT_DIR / "m3_real_usgs_graph_benchmark.csv")
    if graph.empty:
        return
    columns = [
        "graph_family",
        "control_type",
        "prior_strength",
        "n_nodes",
        "n_edges",
        "rmse_single_log10",
        "rmse_graph_log10",
        "delta_rmse_graph_minus_single",
        "within_factor_2_single",
        "within_factor_2_graph",
        "n_violations_before",
        "n_violations_after",
        "improved_vs_single",
    ]
    graph[[col for col in columns if col in graph.columns]].to_csv(
        TABLE_DIR / "Manuscript_Table4_Real_USGS_Graph_Benchmark.csv",
        index=False,
    )


def _scenario_metrics(df: pd.DataFrame, scenario_id: str) -> dict[str, float]:
    """Compute core parity metrics for a specific scenario."""
    subset = df[df["scenario_id"] == scenario_id].copy() if "scenario_id" in df.columns else df.copy()
    subset = subset[subset["ref_age"].notna() & subset["est_age_multi"].notna()]
    subset = subset[pd.to_numeric(subset.get("log10_error"), errors="coerce").notna()]
    if subset.empty:
        return {}
    log_err = pd.to_numeric(subset.get("log10_error"), errors="coerce").dropna()
    within_raw = subset.get("within_factor_2", pd.Series(np.nan, index=subset.index))
    within_2 = pd.to_numeric(within_raw, errors="coerce")
    within_bool = within_raw.astype(str).str.lower().map({"true": 1.0, "false": 0.0})
    within_2 = within_2.where(within_2.notna(), within_bool).dropna()
    ref = pd.to_numeric(subset.get("ref_age"), errors="coerce")
    est = pd.to_numeric(subset.get("est_age_multi"), errors="coerce")
    valid = (ref > 0) & (est > 0)
    y = np.log10(ref[valid])
    residual = np.log10(est[valid]) - y
    ss_tot = float(((y - y.mean()) ** 2).sum()) if len(y) else float("nan")
    ss_res = float((residual**2).sum()) if len(residual) else float("nan")
    return {
        "N": int(len(subset)),
        "finite_metrics": int(len(log_err)),
        "median_abs_log10_error": float(log_err.median()) if len(log_err) else float("nan"),
        "log10_rmse": float(np.sqrt(np.mean(log_err**2))) if len(log_err) else float("nan"),
        "log10_r2": 1.0 - ss_res / ss_tot if ss_tot and np.isfinite(ss_tot) else float("nan"),
        "within_factor_2": float(within_2.mean()) if len(within_2) else float("nan"),
    }


def make_table5_mode_comparison() -> None:
    """Report strict parity, selection, and calibrated-emulation metrics separately."""
    design = _mode_results()
    rows: list[dict] = []
    if not design.empty and "scenario_id" in design.columns:
        for scenario_id, label in (
            ("tracerlpm_strict_parity", "Strict TracerLPM parity"),
            ("tracerlpm_parity_hier_oldwater", "Strict parity + hierarchical old-water priors"),
            ("tracerlpm_parity_agefractions", "Strict parity + age-fraction constraints"),
            ("hydrosheaf_selection_corrected", "Hydrosheaf model selection"),
            ("screened_dgm_gases", "Screened young-gas correction"),
            ("parity_reported_corrected", "Reported-model parity"),
        ):
            if scenario_id in set(design["scenario_id"].dropna()):
                metrics = _scenario_metrics(design, scenario_id)
                if metrics:
                    rows.append({"mode": label, **metrics})

    calibrated = _read_csv(RESULT_DIR / "m3_usgs_calibrated_parity.csv")
    if not calibrated.empty:
        target = pd.to_numeric(calibrated.get("log10_reported_age"), errors="coerce")
        pred = pd.to_numeric(calibrated.get("log10_calibrated_age"), errors="coerce")
        valid = target.notna() & pred.notna()
        if valid.any():
            signed_residual = pred[valid] - target[valid]
            residual = signed_residual.abs()
            y = target[valid]
            ss_tot = float(((y - y.mean()) ** 2).sum()) if len(y) else float("nan")
            ss_res = float((signed_residual**2).sum()) if len(signed_residual) else float("nan")
            rows.append({
                "mode": "USGS-calibrated benchmark emulation",
                "N": int(valid.sum()),
                "finite_metrics": int(len(residual)),
                "median_abs_log10_error": float(residual.median()),
                "log10_rmse": float(np.sqrt(np.mean(residual**2))),
                "log10_r2": 1.0 - ss_res / ss_tot if ss_tot and np.isfinite(ss_tot) else float("nan"),
                "within_factor_2": float((residual <= np.log10(2.0)).mean()),
            })

    if rows:
        pd.DataFrame(rows).to_csv(
            TABLE_DIR / "Manuscript_Table5_Mode_Comparison.csv", index=False
        )


def make_supp_tables() -> None:
    design = _primary_pointwise_results()
    if not design.empty:
        age_class = (
            design[design["log10_error"].notna()]
            .groupby(["scenario_id", "age_class"])
            .agg(
                n=("site_id", "count"),
                median_abs_log10_error=("log10_error", "median"),
                log10_rmse=("log10_error", _log_rmse),
                within_factor_2=("within_factor_2", "mean"),
                within_factor_10=("within_factor_10", "mean"),
                median_proxy_coherence=("young_gas_proxy_coherence", "median"),
            )
            .reset_index()
        )
        age_class.to_csv(TABLE_DIR / "Manuscript_Supp_TableS1_Age_Class_Performance.csv", index=False)

        old = (
            design.groupby(["scenario_id", "old_groundwater_status"], dropna=False)
            .agg(
                n=("site_id", "count"),
                median_c14_age=("old_groundwater_apparent_c14_age", "median"),
                median_he4_age=("old_groundwater_apparent_he4_age", "median"),
                median_constraint_gap_log10=("old_groundwater_constraint_gap_log10", "median"),
            )
            .reset_index()
        )
        old.to_csv(TABLE_DIR / "Manuscript_Supp_TableS2_Old_Groundwater_Diagnostics.csv", index=False)

    audit = _read_csv(RESULT_DIR / "m3_gas_correction_audit_summary.csv")
    if not audit.empty:
        audit.to_csv(TABLE_DIR / "Manuscript_Supp_TableS3_Gas_Correction_Audit.csv", index=False)


def main() -> int:
    make_table1()
    make_table2()
    make_table3()
    make_table4()
    make_table5_mode_comparison()
    make_supp_tables()
    print(f"Wrote M3 manuscript tables to {TABLE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
