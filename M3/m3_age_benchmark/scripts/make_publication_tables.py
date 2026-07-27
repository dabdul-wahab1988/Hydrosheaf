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


def make_table4b_graph_replication() -> None:
    """Cross-aquifer replication of the graph benchmark.

    The primary graph benchmark runs on the national public-supply release
    (1,272 nodes).  Adding the 74 coordinate-bearing MRVA alluvial wells gives a
    structurally different aquifer -- contiguous alluvium, median well spacing
    15.5 km against 40-95 km nationally, with nested piezometers -- and so tests
    whether the graph conclusions hold across aquifer type rather than only on a
    larger sample.

    The Western Principal Aquifers release cannot appear here: it publishes no
    well coordinates, so no graph family can be constructed from it.
    """
    primary = _read_csv(RESULT_DIR / "m3_real_usgs_graph_benchmark.csv")
    repl = _read_csv(RESULT_DIR / "m3_graph_benchmark_national_plus_mrva.csv")
    if primary.empty or repl.empty:
        return
    keys = ["graph_family", "prior_strength"]
    cols = keys + ["n_nodes", "n_edges", "delta_rmse_graph_minus_single",
                   "n_violations_before", "n_violations_after", "improved_vs_single"]
    a = primary[[c for c in cols if c in primary.columns]].copy()
    b = repl[[c for c in cols if c in repl.columns]].copy()
    merged = a.merge(b, on=keys, suffixes=("_national", "_national_plus_mrva"))
    merged["support_national"] = "national public-supply only"
    merged["support_replication"] = "national public-supply + MRVA alluvial"
    merged["replication_note"] = (
        "Western Principal Aquifers is excluded: the release publishes no well "
        "coordinates, so no graph family can be built from it."
    )
    merged.to_csv(
        TABLE_DIR / "Manuscript_Table4b_Graph_Benchmark_Cross_Aquifer_Replication.csv",
        index=False,
    )


def _scenario_metrics_filtered(df: pd.DataFrame, scenario_id: str, site_ids: set | None = None) -> dict[str, float]:
    """Compute core parity metrics for a specific scenario, optionally filtered by site_ids."""
    subset = df[df["scenario_id"] == scenario_id].copy() if "scenario_id" in df.columns else df.copy()
    if site_ids is not None:
        subset = subset[subset["site_id"].isin(site_ids)].copy()
    
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


RELEASE_LABELS = {
    "national_public_supply": "National public-supply (Jurgens et al., 2022a)",
    "western_principal_aquifers": "Western Principal Aquifers (Faulkner & Jurgens, 2019)",
    "mrva_alluvial": "MRVA alluvial (Gratzer et al., 2025a)",
}

# What each release can actually support. Reporting a scenario on a release that
# cannot support it would be worse than omitting it.
RELEASE_CAPABILITY = {
    "national_public_supply": "coordinates + reported LPM: all scenarios",
    "western_principal_aquifers": (
        "no well coordinates published: parity and ablation scenarios only, "
        "excluded from all graph families"
    ),
    "mrva_alluvial": (
        "coordinates available; Table8_LPMs.csv unavailable from ScienceBase "
        "(HTTP 404), so reported-LPM parity cannot be reproduced"
    ),
}


def write_per_release_replication(design: pd.DataFrame, scenarios) -> None:
    """Report each scenario separately per source release.

    The three releases are treated as independent replications rather than as a
    single pooled sample. Pooling would be misleading here: the Western release
    contributes 290 rows that are all Central Valley, which would take that one
    aquifer from 6% to 22% of the benchmark and move aggregate metrics for
    reasons unrelated to method performance. Reporting side by side instead
    tests whether a conclusion holds across aquifer types, which is a stronger
    claim than a larger n.
    """
    if design.empty or "source_release" not in design.columns:
        return
    rows = []
    for release, grp in design.groupby("source_release"):
        for scenario_id, label in scenarios:
            if scenario_id not in set(grp.get("scenario_id", pd.Series(dtype=object))):
                continue
            metrics = _scenario_metrics_filtered(grp, scenario_id)
            if not metrics:
                continue
            rows.append({
                "source_release": RELEASE_LABELS.get(release, release),
                "release_capability": RELEASE_CAPABILITY.get(release, ""),
                "mode": label,
                **metrics,
            })
    if not rows:
        return
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "Manuscript_Table7_Per_Release_Replication.csv", index=False)


def _calibrated_metrics_filtered(df: pd.DataFrame, site_ids: set | None = None) -> dict[str, float]:
    """Compute core parity metrics for calibrated emulation, optionally filtered by site_ids."""
    subset = df.copy()
    if site_ids is not None:
        subset = subset[subset["site_id"].isin(site_ids)].copy()
        
    target = pd.to_numeric(subset.get("log10_reported_age"), errors="coerce")
    pred = pd.to_numeric(subset.get("log10_calibrated_age"), errors="coerce")
    valid = target.notna() & pred.notna()
    if not valid.any():
        return {}
    
    signed_residual = pred[valid] - target[valid]
    residual = signed_residual.abs()
    y = target[valid]
    ss_tot = float(((y - y.mean()) ** 2).sum()) if len(y) else float("nan")
    ss_res = float((signed_residual**2).sum()) if len(signed_residual) else float("nan")
    return {
        "N": int(valid.sum()),
        "finite_metrics": int(len(residual)),
        "median_abs_log10_error": float(residual.median()),
        "log10_rmse": float(np.sqrt(np.mean(residual**2))),
        "log10_r2": 1.0 - ss_res / ss_tot if ss_tot and np.isfinite(ss_tot) else float("nan"),
        "within_factor_2": float((residual <= np.log10(2.0)).mean()),
    }


def make_table5_mode_comparison() -> None:
    """Report strict parity, selection, and calibrated-emulation metrics across three subsets:
    1. Full Available
    2. High-N Common Support (size derived from the scenario intersection)
    3. Design Common Support (size derived from the scenario intersection)
    """
    design = _mode_results()
    calibrated = _read_csv(RESULT_DIR / "m3_usgs_calibrated_parity.csv")

    scenarios = [
        ("tracerlpm_strict_parity", "Strict TracerLPM parity"),
        ("tracerlpm_parity_hier_oldwater", "Strict parity + hierarchical old-water priors"),
        ("tracerlpm_parity_agefractions", "Strict parity + age-fraction constraints"),
        ("hydrosheaf_selection_corrected", "Hydrosheaf model selection"),
        ("screened_dgm_gases", "Screened young-gas correction"),
        ("parity_reported_corrected", "Reported-model parity"),
    ]

    # Independent-replication view: the three releases reported side by side.
    write_per_release_replication(design, scenarios)  # no-op unless design carries source_release

    # Step 1: Compute valid site IDs for each scenario
    valid_sites = {}
    for scenario_id, _ in scenarios:
        if not design.empty and "scenario_id" in design.columns:
            sub = design[design["scenario_id"] == scenario_id]
            valid_sub = sub[sub["ref_age"].notna() & sub["est_age_multi"].notna()]
            valid_sub = valid_sub[pd.to_numeric(valid_sub.get("log10_error"), errors="coerce").notna()]
            valid_sites[scenario_id] = set(valid_sub["site_id"].dropna().unique())
        else:
            valid_sites[scenario_id] = set()

    if not calibrated.empty:
        target = pd.to_numeric(calibrated.get("log10_reported_age"), errors="coerce")
        pred = pd.to_numeric(calibrated.get("log10_calibrated_age"), errors="coerce")
        valid_sub = calibrated[target.notna() & pred.notna()]
        valid_sites["calibrated_emulation"] = set(valid_sub["site_id"].dropna().unique())
    else:
        valid_sites["calibrated_emulation"] = set()

    # Step 2: Define intersections for the subsets
    # 2a. High-N intersection (N=655)
    high_n_scenarios = [
        "tracerlpm_strict_parity",
        "tracerlpm_parity_hier_oldwater",
        "tracerlpm_parity_agefractions",
        "parity_reported_corrected",
        "calibrated_emulation"
    ]
    high_n_sets = [valid_sites[sid] for sid in high_n_scenarios if sid in valid_sites and valid_sites[sid]]
    high_n_intersection = set.intersection(*high_n_sets) if high_n_sets else set()

    # 2b. Design intersection
    all_sets = [s for s in valid_sites.values() if s]
    design_intersection = set.intersection(*all_sets) if all_sets else set()

    # Label the common-support subsets from the intersections actually computed.
    # These were previously hard-coded (and had drifted: the label read N=43
    # while the table carried 272 rows). Deriving them keeps the label correct
    # whenever the underlying scenario coverage changes.
    high_n_label = f"High-N Common Support (N={len(high_n_intersection)})"
    design_label = f"Design Common Support (N={len(design_intersection)})"

    rows: list[dict] = []

    # --- Subset 1: Full Available ---
    for scenario_id, label in scenarios:
        if not design.empty and scenario_id in set(design["scenario_id"].dropna()):
            metrics = _scenario_metrics_filtered(design, scenario_id)
            if metrics:
                rows.append({"subset": "Full Available", "mode": label, **metrics})
    if not calibrated.empty:
        metrics = _calibrated_metrics_filtered(calibrated)
        if metrics:
            rows.append({"subset": "Full Available", "mode": "USGS-calibrated benchmark emulation", **metrics})

    # --- Subset 2: High-N Common Support ---
    for scenario_id, label in scenarios:
        if scenario_id in high_n_scenarios and not design.empty and scenario_id in set(design["scenario_id"].dropna()):
            metrics = _scenario_metrics_filtered(design, scenario_id, high_n_intersection)
            if metrics:
                rows.append({"subset": high_n_label, "mode": label, **metrics})
    if not calibrated.empty:
        metrics = _calibrated_metrics_filtered(calibrated, high_n_intersection)
        if metrics:
            rows.append({"subset": high_n_label, "mode": "USGS-calibrated benchmark emulation", **metrics})

    # --- Subset 3: Design Common Support ---
    for scenario_id, label in scenarios:
        if not design.empty and scenario_id in set(design["scenario_id"].dropna()):
            metrics = _scenario_metrics_filtered(design, scenario_id, design_intersection)
            if metrics:
                rows.append({"subset": design_label, "mode": label, **metrics})
    if not calibrated.empty:
        metrics = _calibrated_metrics_filtered(calibrated, design_intersection)
        if metrics:
            rows.append({"subset": design_label, "mode": "USGS-calibrated benchmark emulation", **metrics})

    if rows:
        pd.DataFrame(rows).to_csv(
            TABLE_DIR / "Manuscript_Table5_Mode_Comparison.csv", index=False
        )


def make_table6_statistical_significance() -> None:
    """Compute Wilcoxon signed-rank, paired t-test, and bootstrap CIs for strict vs age-fraction constraints."""
    df_strict = _read_csv(RESULT_DIR / "m3_tracerlpm_strict_parity_full.csv")
    df_frac = _read_csv(RESULT_DIR / "m3_tracerlpm_parity_agefractions_full.csv")
    
    if df_strict.empty or df_frac.empty:
        return

    df_strict = df_strict[df_strict["ref_age"].notna() & df_strict["est_age_multi"].notna()]
    df_strict["log10_err_strict"] = pd.to_numeric(df_strict["log10_error"], errors="coerce")
    df_strict = df_strict.dropna(subset=["log10_err_strict"])

    df_frac = df_frac[df_frac["ref_age"].notna() & df_frac["est_age_multi"].notna()]
    df_frac["log10_err_frac"] = pd.to_numeric(df_frac["log10_error"], errors="coerce")
    df_frac = df_frac.dropna(subset=["log10_err_frac"])

    merged = pd.merge(
        df_strict[["site_id", "log10_err_strict"]],
        df_frac[["site_id", "log10_err_frac"]],
        on="site_id",
        how="inner"
    )

    if merged.empty:
        return

    x_strict = merged["log10_err_strict"].values
    x_frac = merged["log10_err_frac"].values

    import scipy.stats as stats
    res_wilc = stats.wilcoxon(x_strict, x_frac)
    res_ttest = stats.ttest_rel(x_strict, x_frac)

    np.random.seed(42)
    n_boot = 5000
    boot_diff_mae = []
    boot_diff_rmse = []

    for _ in range(n_boot):
        idx = np.random.choice(len(merged), size=len(merged), replace=True)
        strict_b = x_strict[idx]
        frac_b = x_frac[idx]
        
        mae_strict = np.median(strict_b)
        mae_frac = np.median(frac_b)
        boot_diff_mae.append(mae_strict - mae_frac)
        
        rmse_strict = np.sqrt(np.mean(strict_b**2))
        rmse_frac = np.sqrt(np.mean(frac_b**2))
        boot_diff_rmse.append(rmse_strict - rmse_frac)

    ci_mae = np.percentile(boot_diff_mae, [2.5, 97.5])
    ci_rmse = np.percentile(boot_diff_rmse, [2.5, 97.5])

    rows = [
        {
            "test_or_metric": "Paired Wilcoxon signed-rank test",
            "comparison": "Strict Parity vs Age-Fraction Constraints",
            "statistic": float(res_wilc.statistic),
            "p_value": float(res_wilc.pvalue),
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "interpretation": "Significant difference in absolute error distributions (p < 0.05)." if res_wilc.pvalue < 0.05 else "No significant difference."
        },
        {
            "test_or_metric": "Paired t-test",
            "comparison": "Strict Parity vs Age-Fraction Constraints",
            "statistic": float(res_ttest.statistic),
            "p_value": float(res_ttest.pvalue),
            "ci_lower": np.nan,
            "ci_upper": np.nan,
            "interpretation": "Significant difference in mean absolute errors (p < 0.05)." if res_ttest.pvalue < 0.05 else "No significant difference."
        },
        {
            "test_or_metric": "Bootstrap Difference in MAE (strict - agefractions)",
            "comparison": "Strict Parity vs Age-Fraction Constraints",
            "statistic": float(np.mean(boot_diff_mae)),
            "p_value": np.nan,
            "ci_lower": float(ci_mae[0]),
            "ci_upper": float(ci_mae[1]),
            "interpretation": "Positive difference indicates age-fraction constraints reduce median absolute error (MAE)."
        },
        {
            "test_or_metric": "Bootstrap Difference in RMSE (strict - agefractions)",
            "comparison": "Strict Parity vs Age-Fraction Constraints",
            "statistic": float(np.mean(boot_diff_rmse)),
            "p_value": np.nan,
            "ci_lower": float(ci_rmse[0]),
            "ci_upper": float(ci_rmse[1]),
            "interpretation": "Positive difference indicates age-fraction constraints reduce RMSE."
        }
    ]

    pd.DataFrame(rows).to_csv(TABLE_DIR / "Manuscript_Table6_Statistical_Significance.csv", index=False)



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
        age_class.to_csv(TABLE_DIR / "Manuscript_Supp_TableS2_Age_Class_Performance.csv", index=False)

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
        old.to_csv(TABLE_DIR / "Manuscript_Supp_TableS3_Old_Groundwater_Diagnostics.csv", index=False)

    audit = _read_csv(RESULT_DIR / "m3_gas_correction_audit_summary.csv")
    if not audit.empty:
        audit.to_csv(TABLE_DIR / "Manuscript_Supp_TableS4_Gas_Correction_Audit.csv", index=False)


def main() -> int:
    make_table1()
    make_table2()
    make_table3()
    make_table4()
    make_table4b_graph_replication()
    make_table5_mode_comparison()
    make_table6_statistical_significance()
    make_supp_tables()
    print(f"Wrote M3 manuscript tables to {TABLE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
