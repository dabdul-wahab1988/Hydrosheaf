"""Reader-friendly, presentation-only reformatting of the manuscript tables.

Renames columns to plain labels and rounds values for display; does not alter
any underlying number. The canonical, full-precision CSVs remain at
M3.1/manuscript/artifacts/data/Manuscript_Table*.csv and are the evidence of
record (registered in artifact_registry.json alongside these display copies).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "manuscript" / "artifacts" / "data"
OUT = DATA / "display"
OUT.mkdir(parents=True, exist_ok=True)

SCENARIO_LABELS = {
    "tracerlpm_strict_parity": "Strict reported-configuration parity",
    "tracerlpm_parity_agefractions": "Reported-output-constrained fraction sensitivity",
    "hydrosheaf_selection_corrected": "Hydrosheaf model selection",
    "parity_reported_corrected": "Reported-model parity",
    "screened_dgm_gases": "Screened young-gas correction",
    "oldwater_he4_uncertainty": "Old-water 4He uncertainty",
    "oldwater_c14_ensemble": "Old-water 14C ensemble",
    "oldwater_ensemble_he4_uncertainty": "Old-water ensemble 4He uncertainty",
    "ablation_no_he4": "Ablation: no 4He",
    "ablation_raw_c14": "Ablation: raw 14C",
    "ablation_raw_gases": "Ablation: raw gases",
    "tracer_young_only": "Young tracers only",
    "tracer_old_only": "Old tracers only",
}

FAMILY_LABELS = {
    "coordinate_global": "Coordinate (global)",
    "study_unit_coordinate": "Coordinate (study unit)",
    "depth_constrained": "Depth constrained",
    "hydraulic_proxy_constrained": "Hydraulic proxy",
    "parameter_smooth_aicc": "Parameter smoothing",
    "wrong_direction_negative_control": "Wrong-direction control",
    "randomized_negative_control": "Randomised control",
}


def r(df: pd.DataFrame, cols: dict[str, int]) -> pd.DataFrame:
    for col, digits in cols.items():
        if col in df.columns:
            df[col] = df[col].round(digits)
    return df


def pct(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = (df[col] * 100).round(1)
    return df


def table2():
    df = pd.read_csv(DATA / "Manuscript_Table2_Design_Matrix_Performance.csv")
    df["Scenario"] = df["scenario_id"].map(SCENARIO_LABELS).fillna(df["scenario_id"])
    df = r(df, {"median_abs_log10_error": 4, "log10_rmse": 4, "log10_r2": 3})
    df = pct(df, ["within_factor_2", "within_factor_10"])
    out = df[["Scenario", "total_rows", "identifiable_rows", "median_abs_log10_error",
              "log10_rmse", "log10_r2", "within_factor_2", "within_factor_10"]]
    out.columns = ["Scenario", "Total wells", "Reportable N", "Median |log10 error|",
                   "log10 RMSE", "log10 R2", "Within factor 2 (%)", "Within factor 10 (%)"]
    out.to_csv(OUT / "Table2_display.csv", index=False)


def table3():
    df = pd.read_csv(DATA / "Manuscript_Table3_Paired_Ablation_Effects.csv")
    df["Scenario"] = df["scenario_id"].map(SCENARIO_LABELS).fillna(df["scenario_id"])
    df = r(df, {"median_delta_log10_error": 4, "mean_delta_log10_error": 4})
    df = pct(df, ["improved_fraction"])
    out = df[["Scenario", "paired_rows", "median_delta_log10_error", "mean_delta_log10_error",
              "improved_fraction", "gained_factor_2_rows", "lost_factor_2_rows"]]
    out.columns = ["Scenario", "Paired N", "Median delta |log10 error|", "Mean delta |log10 error|",
                   "Improved fraction (%)", "Gained factor-2", "Lost factor-2"]
    out.to_csv(OUT / "Table3_display.csv", index=False)


def table4():
    df = pd.read_csv(DATA / "Manuscript_Table4_Real_USGS_Graph_Benchmark.csv")
    df = df[df["graph_family"].isin(FAMILY_LABELS)]
    df["Graph family"] = df["graph_family"].map(FAMILY_LABELS)
    df["Prior strength"] = df["prior_strength"].str.capitalize()
    df = r(df, {"rmse_single_log10": 4, "rmse_graph_log10": 4, "delta_rmse_graph_minus_single": 5})
    df = pct(df, ["within_factor_2_single", "within_factor_2_graph"])
    df["Meets robust-improvement rule"] = df["improved_vs_single"].map({True: "Yes", False: "No"})
    out = df[["Graph family", "Prior strength", "n_nodes", "n_edges", "rmse_single_log10",
              "rmse_graph_log10", "delta_rmse_graph_minus_single", "within_factor_2_single",
              "within_factor_2_graph", "n_violations_before", "n_violations_after",
              "Meets robust-improvement rule"]]
    out.columns = ["Graph family", "Prior strength", "Nodes", "Edges", "log10 RMSE (single)",
                   "log10 RMSE (graph)", "Delta log10 RMSE", "Within factor 2, single (%)",
                   "Within factor 2, graph (%)", "Violations before", "Violations after",
                   "Meets robust-improvement rule"]
    out.to_csv(OUT / "Table4_display.csv", index=False)


def table5():
    df = pd.read_csv(DATA / "Manuscript_Table5_Mode_Comparison.csv")
    df = r(df, {"median_abs_log10_error": 4, "log10_rmse": 4, "log10_r2": 3})
    df = pct(df, ["within_factor_2"])
    out = df[["subset", "mode", "N", "median_abs_log10_error", "log10_rmse", "log10_r2",
              "within_factor_2"]]
    out.columns = ["Support", "Mode", "N", "Median |log10 error|", "log10 RMSE", "log10 R2",
                   "Within factor 2 (%)"]
    out.to_csv(OUT / "Table5_display.csv", index=False)


def table6():
    df = pd.read_csv(DATA / "Manuscript_Table6_Statistical_Significance.csv")

    def fmt_p(v):
        if pd.isna(v):
            return ""
        return "< 0.001" if v < 0.001 else f"{v:.3f}"

    df["p_value"] = df["p_value"].apply(fmt_p)
    for c in ("statistic", "ci_lower", "ci_upper"):
        df[c] = df[c].round(4)
    out = df[["test_or_metric", "comparison", "statistic", "p_value", "ci_lower", "ci_upper",
              "interpretation"]]
    out.columns = ["Test or metric", "Comparison", "Statistic", "p-value", "95% CI lower",
                   "95% CI upper", "Interpretation"]
    out.to_csv(OUT / "Table6_display.csv", index=False)


if __name__ == "__main__":
    table2()
    table3()
    table4()
    table5()
    table6()
    print("wrote display tables to", OUT)
