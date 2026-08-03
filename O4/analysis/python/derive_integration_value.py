"""TAB-2/FIG-4 (M7 rows): entropy reduction (internal) vs predictive-skill
change (external) under native evidence combination and under adverse
controls, plus the M7.4/M7.5 sheaf-vs-graph calibration/ranking divergence.

Reads M7.3's own already-computed paired bootstrap contrasts
(`evidence_case_bootstrap_contrasts.csv`) and M7.4/M7.5's own already-computed
paired bootstrap contrasts (`paired_bootstrap_contrasts.csv`) and passes them
through with no recomputation of the underlying inference -- only relabelling
into the common internal-signal/external-signal schema.

Run:  .venv/Scripts/python.exe O4/analysis/python/derive_integration_value.py
"""

from __future__ import annotations

import pandas as pd

from _common import M7, write


def native_and_adverse_controls() -> pd.DataFrame:
    df = pd.read_csv(M7 / "results" / "m7_3_locked" / "evidence_case_bootstrap_contrasts.csv")
    wide = df.pivot_table(
        index=["contrast", "condition"], columns="metric", values=["mean_difference", "ci95_low", "ci95_high"]
    )
    wide.columns = [f"{a}__{b}" for a, b in wide.columns]
    wide = wide.reset_index()

    label = {
        "native_incremental_age": "native: age added to HC",
        "native_incremental_chemistry": "native: chemistry added to HA",
        "native_incremental_hydraulics": "native: hydraulics added to AC",
        "permuted_age_increment": "adverse control: age permuted",
    }
    keep = wide.loc[wide["contrast"].isin(label)].copy()
    keep["axis_x"] = keep["contrast"].map(label)
    keep["component"] = "M7 identifiability"
    keep["internal_signal_name"] = "mean_edge_entropy_change"
    keep["internal_signal_value"] = keep["mean_difference__mean_edge_entropy"]
    keep["external_signal_name"] = "pr_auc_change"
    keep["external_signal_value"] = keep["mean_difference__pr_auc"]
    keep["external_ci_low"] = keep["ci95_low__pr_auc"]
    keep["external_ci_high"] = keep["ci95_high__pr_auc"]
    keep["internal_ci_low"] = keep["ci95_low__mean_edge_entropy"]
    keep["internal_ci_high"] = keep["ci95_high__mean_edge_entropy"]
    keep["brier_change"] = keep["mean_difference__brier"]
    keep["log_loss_change"] = keep["mean_difference__log_loss"]
    keep["n_cases"] = 12
    keep["source"] = "results/m7_3_locked/evidence_case_bootstrap_contrasts.csv"
    out_cols = [
        "component", "axis_x", "internal_signal_name", "internal_signal_value",
        "internal_ci_low", "internal_ci_high", "external_signal_name",
        "external_signal_value", "external_ci_low", "external_ci_high",
        "brier_change", "log_loss_change", "n_cases", "source",
    ]
    return keep[out_cols]


def sheaf_vs_graph_calibration_divergence() -> pd.DataFrame:
    """M7.4, all scenarios pooled: ranking (PR-AUC, external) is
    statistically tied between the sheaf and the stronger edge-local
    weighted graph (95% CI crosses zero), while calibration (log-loss, an
    internal-consistency signal) is significantly worse for the sheaf (95%
    CI excludes zero, entirely positive) -- reported exactly as M7.4's own
    locked contrast file states it, with no recomputation."""
    df = pd.read_csv(
        M7 / "results" / "RUN-M7-SHEAF-VS-GRAPH-20260729-01" / "paired_bootstrap_contrasts.csv"
    )
    rows = df.loc[
        (df["scenario"] == "all")
        & (df["metric"].isin(["pr_auc", "log_loss", "brier"]))
        & (df["right"] == "weighted_graph")
    ].copy()
    rows["component"] = "M7 identifiability"
    rows["axis_x"] = "M7.4 all scenarios: sheaf vs weighted graph"
    rows["source"] = "results/RUN-M7-SHEAF-VS-GRAPH-20260729-01/paired_bootstrap_contrasts.csv"
    return rows


def main() -> None:
    write(native_and_adverse_controls(), "integration_value.csv")
    write(sheaf_vs_graph_calibration_divergence(), "integration_calibration_divergence.csv")


if __name__ == "__main__":
    main()
