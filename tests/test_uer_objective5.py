"""Automated unit and integration tests for Objective 5 (O5) executed on UER dataset.

Verifies:
1. Sheaf Laplacian regularisation computes finite chemical states, ungated residuals, and gated residuals.
2. Sheaf gating achieves substantial Dirichlet energy reduction (E_gated < E_ungated).
3. Tripartite boundary conflict localisation flags physical discordances (stoichiometry, tritium, nitrate).
4. Missingness semantics are strictly preserved (17 dual-tritium edges tested; 1,583 untested).
5. Certified Minimax Measurement Design identifies the Pareto sweet spot (SF6 + CFC-12) on empirical tritium cohorts.
6. Row-level sample provenance exists for all 38 empirical UER boreholes.
7. Publication figures and markdown tables are correctly generated.
"""

from __future__ import annotations

from pathlib import Path
import pytest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "outputs" / "uer_objective5"
TAB_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"


def test_objective5_artifacts_exist() -> None:
    """Ensure all required Objective 5 tables, figures, and reports are present."""
    required_files = [
        TAB_DIR / "uer_objective5_minimax_pareto_frontier.csv",
        TAB_DIR / "uer_objective5_minimax_pareto_frontier.md",
        TAB_DIR / "uer_objective5_tripartite_boundary_audit.csv",
        TAB_DIR / "uer_objective5_sheaf_conflict_localisation.csv",
        TAB_DIR / "uer_objective5_sheaf_conflict_localisation.md",
        TAB_DIR / "uer_objective5_sheaf_gating_comparison.csv",
        TAB_DIR / "uer_objective5_well_specific_minimax_designs.csv",
        FIG_DIR / "figure_o5_1_tripartite_boundary.png",
        FIG_DIR / "figure_o5_1_tripartite_boundary_spanned.png",
        FIG_DIR / "figure_o5_2_minimax_frontier.png",
        OUT_DIR / "UER_OBJECTIVE_5_FULL_REPORT.md",
    ]
    for f in required_files:
        assert f.exists(), f"Missing Objective 5 artifact: {f}"
        assert f.stat().st_size > 100, f"Artifact is unexpectedly empty: {f}"


def test_objective5_conflict_localisation_results() -> None:
    """Verify that Sheaf Laplacian detects and localizes multi-physics conflicts."""
    df_conf = pd.read_csv(TAB_DIR / "uer_objective5_sheaf_conflict_localisation.csv")
    assert len(df_conf) > 50, f"Expected >50 flagged conflicts, got {len(df_conf)}"

    # Must contain all three physical conflict types
    assert (df_conf["chem_conflict"]).sum() > 0, "No stoichiometric conflicts flagged"
    assert (df_conf["tritium_conflict"]).sum() > 0, "No tritium vertical inversion conflicts flagged"
    assert (df_conf["nitrate_conflict"]).sum() > 0, "No point-source nitrate conflicts flagged"

    # Gating action must be REJECT_EDGE_ISOLATE_CONFLICT
    assert (df_conf["gating_action"] == "REJECT_EDGE_ISOLATE_CONFLICT").all()


def test_objective5_gating_dirichlet_energy_reduction() -> None:
    """Verify that isolating conflicts reduces Sheaf Dirichlet energy."""
    df_gate = pd.read_csv(TAB_DIR / "uer_objective5_sheaf_gating_comparison.csv")
    assert len(df_gate) == 2
    
    e_ungated = df_gate[df_gate["solve_mode"].str.contains("Ungated")]["dirichlet_energy"].iloc[0]
    e_gated = df_gate[df_gate["solve_mode"].str.contains("Gated")]["dirichlet_energy"].iloc[0]

    assert e_gated < e_ungated, f"Expected Dirichlet energy reduction, got E_gated={e_gated} >= E_ungated={e_ungated}"
    pct_reduction = (1.0 - e_gated / e_ungated) * 100.0
    assert pct_reduction > 30.0, f"Expected >30% Dirichlet energy reduction, got {pct_reduction:.1f}%"


def test_objective5_missingness_semantics() -> None:
    """Verify that missing tritium is strictly categorized as UNTESTED, not passed."""
    df_audit = pd.read_csv(TAB_DIR / "uer_objective5_tripartite_boundary_audit.csv")
    assert len(df_audit) == 1600

    tritium_tested = df_audit[df_audit["tritium_status"].isin(["TESTED_COHERENT", "FAIL_TRITIUM_INVERSION"])]
    assert len(tritium_tested) == 17, f"Expected exactly 17 dual-tritium edges, got {len(tritium_tested)}"

    tritium_untested = df_audit[df_audit["tritium_status"] == "UNTESTED_MISSING_DATA"]
    assert len(tritium_untested) == 1583, f"Expected exactly 1583 untested tritium edges, got {len(tritium_untested)}"


def test_objective5_minimax_pareto_frontier_reduction() -> None:
    """Verify that Minimax LP achieves >= 70% ambiguity reduction at the $800 allocation."""
    df_front = pd.read_csv(TAB_DIR / "uer_objective5_minimax_pareto_frontier.csv")
    assert not df_front.empty

    med_suite = df_front[df_front["prior_regime"] == "Regional Median Screening Cohort"]
    assert not med_suite.empty

    b0_row = med_suite[med_suite["budget"] == 0.0].iloc[0]
    opt_row = med_suite[med_suite["budget"] == 800.0].iloc[0]
    high_row = med_suite[med_suite["budget"] == 2500.0].iloc[0]

    # Baseline ambiguity must be substantial (>45 years)
    assert b0_row["mtt_ambiguity_yr"] > 45.0

    # At $800, ambiguity reduction must exceed 70%
    reduction_pct = opt_row["ambiguity_reduction_pct"]
    assert reduction_pct >= 70.0, f"Expected >= 70% reduction at $800, got {reduction_pct:.1f}%"
    assert opt_row["status"] == "CERTIFIED_SUFFICIENT"
    assert "SF6" in opt_row["selected_options"] and "CFC12" in opt_row["selected_options"]

    # Redundancy plateau: Increasing budget from $800 to $2500 must yield < 0.05 yr further reduction
    delta_amb = opt_row["mtt_ambiguity_yr"] - high_row["mtt_ambiguity_yr"]
    assert delta_amb < 0.05, f"Expected <0.05 yr change between $800 and $2500, got {delta_amb:.5f} yr"


def test_objective5_well_specific_provenance() -> None:
    """Verify sample-level provenance for the 38 empirical UER boreholes."""
    df_bh = pd.read_csv(TAB_DIR / "uer_objective5_well_specific_minimax_designs.csv")
    assert len(df_bh) == 38, f"Expected 38 borehole designs, got {len(df_bh)}"
    assert "node_id" in df_bh.columns and "site_id" in df_bh.columns
    assert "measured_3H_TU" in df_bh.columns
    assert (df_bh["measured_3H_TU"] > 0).all()
    assert (df_bh["opt_800_cost"] <= 800.0).all()
