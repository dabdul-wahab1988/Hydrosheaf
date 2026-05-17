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
    summary = _read_csv(RESULT_DIR / "m3_design_matrix_summary.csv")
    if summary.empty:
        return
    columns = [
        "scenario_id",
        "total_rows",
        "metric_rows",
        "median_abs_log10_error",
        "log10_rmse",
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


def make_supp_tables() -> None:
    design = _read_csv(RESULT_DIR / "m3_design_matrix_results.csv")
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
    make_supp_tables()
    print(f"Wrote M3 manuscript tables to {TABLE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
