"""Unit tests verifying scientific claims, guardrails, and output ordering for public archives."""

from __future__ import annotations

import json
from pathlib import Path
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results"


def test_archive_summary_ordering():
    """Verify that external_modpath_archive_summary.csv contains three rows ordered: Savage, Great Miami, Long Island."""
    summary_path = RESULTS_DIR / "external_modpath_archive_summary.csv"
    assert summary_path.exists(), "Combined archive summary file is missing."

    df = pd.read_csv(summary_path)
    assert len(df) == 3, f"Expected 3 rows in combined summary, found {len(df)}"

    # Check validation tier ordering
    assert df.loc[0, "validation_tier"] == "tier_1_savage"
    assert df.loc[1, "validation_tier"] == "tier_2_great_miami"
    assert df.loc[2, "validation_tier"] == "tier_3_long_island"


def test_scientific_claims_and_guardrails():
    """Verify that scientific limits, allowed claims, and guardrails are written to output summaries."""
    summary_path = RESULTS_DIR / "external_modpath_archive_summary.csv"
    df = pd.read_csv(summary_path)

    # Check Savage guardrails and interpretations
    savage_row = df[df["validation_tier"] == "tier_1_savage"].iloc[0]
    
    assert "claim_guardrail" in savage_row
    assert "External archive validates MODPATH endpoint/pathline topology" in savage_row["claim_guardrail"]
    assert "does not validate field geochemistry" in savage_row["claim_guardrail"]

    assert "travel_time_rank_interpretation" in savage_row
    assert "supportive only when endpoint and pathline time definitions are harmonised" in savage_row["travel_time_rank_interpretation"]

    # Check Long Island stub claims
    li_row = df[df["validation_tier"] == "tier_3_long_island"].iloc[0]
    assert "Fallback stub" in li_row["claim_guardrail"]


def test_source_manifest_integrity():
    """Verify that the external_modpath_source_manifest.json contains metadata for all three validation tiers."""
    manifest_path = RESULTS_DIR / "external_modpath_source_manifest.json"
    assert manifest_path.exists(), "Source manifest JSON is missing."

    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    assert "archives" in manifest
    archives = manifest["archives"]
    assert "tier_1_savage" in archives
    assert "tier_2_great_miami" in archives
    assert "tier_3_long_island" in archives

    assert "m4_outputs" in manifest
    outputs = manifest["m4_outputs"]
    assert "archive_summary" in outputs
    assert "edge_agreement" in outputs


def test_negative_controls_are_not_independent_validation():
    """Diagnostic controls must not be labelled as independent graph inference."""
    metrics_path = RESULTS_DIR / "independent_graph_vs_modpath.csv"
    assert metrics_path.exists(), "Topology metrics file is missing."
    df = pd.read_csv(metrics_path)

    controls = df[df["scenario"].isin(["negative_random", "negative_wrong_direction", "negative_shortcut"])]
    assert not controls.empty
    assert set(controls["validation_mode"]) == {"diagnostic_negative_control"}
    assert not controls["independent_validation"].astype(bool).any()


def test_sparse_node_is_sensitivity_not_primary_confusion_metric():
    """Sparse-node rows are sensitivity evidence and should not invent TP/FP/FN counts."""
    metrics_path = RESULTS_DIR / "independent_graph_vs_modpath.csv"
    df = pd.read_csv(metrics_path)
    sparse = df[df["scenario"] == "sparse_node"]
    assert len(sparse) == 1
    row = sparse.iloc[0]
    assert row["validation_mode"] == "sensitivity_analysis"
    assert pd.isna(row["tp"])
    assert pd.isna(row["fp"])
    assert pd.isna(row["fn"])
