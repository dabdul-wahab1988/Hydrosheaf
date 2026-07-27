"""Build M6 main and supplementary tables (CSV + Markdown) from results CSVs."""

from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as m6

RES = m6.RESULTS_DIR
TAB = m6.BENCH_DIR / "tables"
TAB.mkdir(parents=True, exist_ok=True)
FIELD_RESULTS = (
    m6.REPO_ROOT
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "results"
    / "supporting_validation"
)


def _to_md(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join(["---"] * len(cols)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append(
            "| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
        )
    return "\n".join(lines)


def _write(df: pd.DataFrame, name: str, title: str):
    df.to_csv(TAB / f"{name}.csv", index=False)
    (TAB / f"{name}.md").write_text(f"# {title}\n\n{_to_md(df)}\n", encoding="utf-8")
    print("wrote", name)


def rd(f):
    return pd.read_csv(RES / f)


# --- Table 1: dataset readiness & claim strength -----------------------------
def table1():
    r = rd("m6_dataset_readiness.csv")
    claim = {
        "northern_ghana": "M6 reaction-component workflow; class/equivalence inference",
        "manu": "External transfer; reaction inference, no Sr/SiO2 corroboration",
        "talensi": "External transfer; screening only (charge-balance limited)",
    }
    r["claim_strength"] = r["dataset"].map(claim)
    r["dataset"] = r["dataset"].map(
        {
            "northern_ghana": "Northern Ghana",
            "manu": "Lower Anayari",
            "talensi": "Talensi",
        }
    )
    cols = [
        "dataset",
        "n_samples",
        "n_sites",
        "seasonal",
        "native_tier",
        "n_quantitative",
        "n_screening",
        "n_exploratory",
        "cbe_median_abs",
        "claim_strength",
    ]
    _write(
        r[cols].round(2),
        "table1_dataset_readiness",
        "Table 1. Dataset readiness and claim strength",
    )


# --- Table 2: variable availability by dataset -------------------------------
def table2():
    a = rd("m6_variable_availability.csv")
    piv = a.pivot(index="variable", columns="dataset", values="status")
    order = [
        "Ca",
        "Mg",
        "Na",
        "K",
        "HCO3",
        "Cl",
        "SO4",
        "NO3",
        "F",
        "Fe",
        "d18O",
        "d2H",
        "Sr_mgL",
        "SiO2_mgL",
        "Calcite_SI",
        "Region",
    ]
    piv = piv.reindex(order)
    piv = piv[["northern_ghana", "manu", "talensi"]].rename(
        columns={
            "northern_ghana": "Northern Ghana",
            "manu": "Lower Anayari",
            "talensi": "Talensi",
        }
    )
    _write(
        piv.reset_index(),
        "table2_variable_availability",
        "Table 2. Variable availability by dataset",
    )


# --- Table 3: Hydrosheaf outputs under full vs reduced information ------------
def table3():
    tr = rd("m6_tier_ablation_transitions.csv")
    tlab = {
        "tier0_majors": "Tier 0 (majors)",
        "tier1_isotopes": "Tier 1 (+isotopes)",
        "tier2_fluoride": "Tier 2 (+F)",
        "tier3_sr_sio2": "Tier 3 (+Sr/SiO2)",
        "tier4_full_metadata": "Tier 4 (+SI, maximum M6 chemistry tier)",
    }
    tr["tier"] = tr["tier"].map(tlab)
    cols = [
        "tier",
        "n",
        "frac_non_identifiable",
        "frac_partial",
        "frac_class_changed_vs_tier4",
        "frac_family_changed_vs_tier4",
        "mean_mrs",
        "mean_stability",
    ]
    _write(
        tr[cols].round(3),
        "table3_tier_ablation",
        "Table 3. Northern Ghana Hydrosheaf outputs under full vs reduced information",
    )


# --- Table 4: external transfer performance ----------------------------------
def table4():
    s = rd("m6_external_summary.csv")
    _write(
        s.round(3),
        "table4_external_transfer",
        "Table 4. External transfer performance and uncertainty summary",
    )


def table5():
    """Truth-free Northern Ghana seasonal hold-forward performance."""

    summary = pd.read_csv(FIELD_RESULTS / "field_prequential_summary.csv")
    overall = summary[summary["ion"] == "ALL"].copy()
    overall["method"] = overall["method"].map(
        {
            "persistence": "Persistence",
            "expanding_mean_delta": "Expanding mean delta",
            "hydrosheaf_graph_ridge": "HydroSheaf graph ridge",
        }
    )
    columns = [
        "method",
        "n",
        "mae_log1p",
        "rmse_log1p",
        "bias_log1p",
        "coverage90",
    ]
    _write(
        overall[columns].round(4),
        "table5_field_prequential",
        "Table 5. Northern Ghana within-campaign seasonal hold-forward",
    )


# --- Supplementary tables ----------------------------------------------------
def supp_tables():
    # S2 missingness / QC
    cbe = rd("m6_cbe_distribution.csv")
    s2 = (
        cbe.groupby(["dataset", "cbe_class"])
        .size()
        .rename("n")
        .reset_index()
        .pivot(index="dataset", columns="cbe_class", values="n")
        .fillna(0)
        .astype(int)
    )
    _write(
        s2.reset_index(),
        "tableS2_qc_missingness",
        "Table S2. Charge-balance quality classes by dataset",
    )

    # S3 region/season summary
    _write(
        rd("m6_ng_aquifer_season_summary.csv").round(3),
        "tableS3_aquifer_summary",
        "Table S3. Northern Ghana region summary",
    )

    # S4 full tier ablation (per class counts)
    abl = rd("m6_tier_ablation.csv")
    s4 = (
        abl.groupby(["tier", "resolution_class"])
        .size()
        .rename("n")
        .reset_index()
        .pivot(index="tier", columns="resolution_class", values="n")
        .fillna(0)
        .astype(int)
    )
    _write(
        s4.reset_index(),
        "tableS4_tier_ablation_full",
        "Table S4. Full tier-ablation identifiability counts",
    )

    # S5 edge sensitivity
    _write(
        rd("m6_edge_network_summary.csv").round(3),
        "tableS5_edge_sensitivity",
        "Table S5. Edge-set sensitivity summary",
    )

    # S6 external per-edge outputs
    _write(
        rd("m6_external_transfer.csv").round(3),
        "tableS6_external_outputs",
        "Table S6. External-dataset transfer outputs (per edge)",
    )

    # S7 uncertainty metric definitions
    defs = pd.DataFrame(
        [
            (
                "Support stability",
                "Mean Jaccard of selected reaction support across 5% "
                "analytical-noise bootstrap resamples (0-1, higher = more stable)",
            ),
            ("MRS", "Mechanism Resolution Score (0-100) from frozen M5 calibration"),
            (
                "Held-out ion RMSE",
                "Leave-one-ion predictive error of the concentration change",
            ),
            (
                "Reactive RMSE",
                "Residual after Cl-conservative transport correction and inversion",
            ),
            (
                "Identifiability class",
                "non/partially/equivalence-class/identifiable, transferred "
                "classifier + evidence-lift gate",
            ),
            (
                "Evidence corroboration",
                "Whether the dominant process family is supported by an "
                "available tracer (isotopes/Sr/SiO2/F/SI)",
            ),
        ],
        columns=["Metric", "Definition"],
    )
    _write(
        defs,
        "tableS7_uncertainty_definitions",
        "Table S7. Uncertainty-metric definitions",
    )

    # S8 software environment
    import platform

    env = pd.DataFrame(
        [
            ("Python", platform.python_version()),
            ("numpy", np.__version__),
            ("pandas", pd.__version__),
            ("Inverse solver", "M5 fit_inverse (FISTA, column-normalised, SI-bounded)"),
            ("MRS calibration", "frozen M5 mrs_calibration_model.json (not re-fit)"),
            ("Seed", str(m6.SEED)),
        ],
        columns=["Component", "Version / setting"],
    )
    _write(
        env, "tableS8_environment", "Table S8. Software and computational environment"
    )

    # S9 competing no-flow explanations
    null_summary = rd("m6_null_sensitivity_summary.csv").copy()
    numeric = null_summary.select_dtypes(include=[np.number]).columns
    null_summary[numeric] = null_summary[numeric].round(3)
    _write(
        null_summary,
        "tableS9_null_sensitivity",
        "Table S9. Competing no-flow explanation sensitivity",
    )


def main():
    table1()
    table2()
    table3()
    table4()
    table5()
    supp_tables()
    print("M6 tables ->", TAB)


if __name__ == "__main__":
    main()
