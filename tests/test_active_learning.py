"""Tests for active-learning measurement recommendation workflow.

Covers:
- Disagreement scoring (assumption-sensitive, calibration-removed, etc.)
- FP/FN edges rank above TP
- Missing-data measurement recommendations
- Selected-unlabeled reason mapping
- Evidence flags map to correct measurements
- top_k truncation
- Output files written (JSON/CSV/MD)
- Graceful handling of missing benchmark report
- VALIDATED edges have appropriate manuscript_safe_note
- Empty inputs produce valid empty structure
- Bootstrap instability modulation from benchmark CI width
- Priority scores in [0,1] range
"""

import json
import os
import tempfile
from pathlib import Path

import pytest

from hydrosheaf.calibration.active_learning import (
    _compute_disagreement_score,
    _compute_bootstrap_instability_modulation,
    _score_evidence_ambiguity,
    _score_validation_error,
    _score_missing_data,
    _recommended_measurements,
    _manuscript_safe_note,
    rank_next_measurements,
)
from hydrosheaf.config import Config as HConfig
from hydrosheaf.graph.types import Edge


# ── helpers ──────────────────────────────────────────────────────────

def _make_temp_dir():
    return tempfile.TemporaryDirectory()


def _make_edge(edge_id="E1", u="A", v="B", attrs=None):
    return Edge(edge_id=edge_id, u=u, v=v, attrs=attrs or {})


def _make_benchmark_report(variants=None, uncertainty=None):
    """Minimal benchmark report matching run_assumption_benchmark output."""
    return {
        "variants": variants or {},
        "improvement_summary": {},
        "independent_validation": True,
        "manuscript_claim_allowed": True,
        "uncertainty": uncertainty,
    }


def _make_validation_report(label_file=None):
    return {"validation_label_file": label_file}


# ── disagreement score tests ─────────────────────────────────────────

def test_disagreement_assumption_sensitive():
    """Calibrated selects, baseline does not → max disagreement score."""
    score, reason = _compute_disagreement_score(
        "E1",
        {"baseline": False, "null_model_defaults": False, "assumption_calibrated": True},
    )
    assert score == 1.0
    assert "assumption-sensitive" in reason


def test_disagreement_calibration_removed():
    """Baseline selects, calibrated removed → high score."""
    score, reason = _compute_disagreement_score(
        "E1",
        {"baseline": True, "null_model_defaults": False, "assumption_calibrated": False},
    )
    assert score == 0.8
    assert "calibration-removed" in reason


def test_disagreement_all_agree():
    """All variants agree → zero score."""
    score, reason = _compute_disagreement_score(
        "E1",
        {"baseline": True, "null_model_defaults": True, "assumption_calibrated": True},
    )
    assert score == 0.0
    assert "all-variants-agree" in reason


# ── bootstrap instability tests ──────────────────────────────────────

def test_bootstrap_instability_modulation_wide_ci():
    """Wide delta_f1 CI + low probability triggers instability boost."""
    report = _make_benchmark_report(uncertainty={
        "improvement_summary_ci": {
            "delta_f1": [0.1, 0.8],
            "probability_delta_f1_gt_0": 0.82,
        }
    })
    mod, reason = _compute_bootstrap_instability_modulation(report)
    assert mod > 0.0
    assert "wide delta_f1 CI" in reason
    assert "P(delta_f1>0)" in reason


def test_bootstrap_instability_no_uncertainty():
    """No uncertainty section → zero modulation."""
    report = _make_benchmark_report()
    mod, reason = _compute_bootstrap_instability_modulation(report)
    assert mod == 0.0
    assert reason == ""


# ── evidence ambiguity tests ─────────────────────────────────────────

def test_evidence_ambiguity_age_reversal():
    """age_reversal flag → highest evidence ambiguity score."""
    edge = _make_edge(attrs={"evidence_flags": "age_reversal"})
    score, reason = _score_evidence_ambiguity(edge)
    assert score == 0.7
    assert "age_reversal" in reason


def test_evidence_ambiguity_no_flags():
    """No flags → zero score."""
    edge = _make_edge()
    score, reason = _score_evidence_ambiguity(edge)
    assert score == 0.0


def test_evidence_flags_map_to_measurements():
    """age_reversal flag maps to correct measurement recommendations."""
    recs = _recommended_measurements(["age_reversal"])
    assert any("groundwater age tracer" in r for r in recs)
    assert any("tritium" in r for r in recs)


def test_evidence_flags_map_isotope_measurements():
    """Missing isotope flags map to d18O/d2H sampling."""
    recs = _recommended_measurements(["missing_isotopes"])
    assert any("d18O/d2H sampling" in r for r in recs)


# ── validation error tests ───────────────────────────────────────────

def test_fp_fn_rank_above_tp():
    """FP edges score 1.0, FN edges score 0.9, TP edges score 0.1."""
    # FP: selected but absent
    fp_score, fp_status, _ = _score_validation_error(
        "E_fp", None, {"E_fp": 0.0}, {"E_fp"},
    )
    assert fp_score == 1.0
    assert fp_status == "false_positive"

    # FN: present but not selected
    fn_score, fn_status, _ = _score_validation_error(
        "E_fn", None, {"E_fn": 1.0}, set(),
    )
    assert fn_score == 0.9
    assert fn_status == "false_negative"

    # TP: present and selected
    tp_score, tp_status, _ = _score_validation_error(
        "E_tp", None, {"E_tp": 1.0}, {"E_tp"},
    )
    assert tp_score == 0.1
    assert tp_status == "observed_present"

    assert fp_score > fn_score > tp_score


def test_selected_unlabeled_reason():
    """Unlabeled edge selected by calibrated variant → appropriate reason."""
    # val_obs_by_edge must be non-empty (has other edges) but not contain E_ul
    score, status, reason = _score_validation_error(
        "E_ul", None, {"E_other": 1.0}, {"E_ul"},
    )
    assert score == 0.5
    assert status == "selected_unlabeled"
    assert "validation label missing" in reason


def test_tn_edge_scores_zero():
    """Correctly excluded absent edge → zero score."""
    score, status, _ = _score_validation_error(
        "E_tn", None, {"E_tn": 0.0}, set(),
    )
    assert score == 0.0
    assert status == "true_negative"


# ── missing data tests ───────────────────────────────────────────────

def test_missing_data_recommendations():
    """Missing isotope data → recommendation includes d18O/d2H sampling."""
    edge = _make_edge(u="N1", v="N2")
    samples = [
        {"site_id": "N1", "hydraulic_head": 100.0, "d18O": None, "age": 50.0},
        {"site_id": "N2", "hydraulic_head": 95.0, "d18O": -5.0, "age": None},
    ]
    config = HConfig()
    # Align config field keys with sample dict keys
    config.edge_head_key = "hydraulic_head"
    config.edge_elevation_key = "elevation"
    config.isotope_d18o_key = "d18O"
    config.isotope_d2h_key = "d2H"
    score, reasons = _score_missing_data(edge, samples, config)
    assert score > 0.0
    # Should have missing_isotopes for N1 and missing_age for N2
    assert any("missing_isotopes" in r for r in reasons)
    assert any("missing_age" in r for r in reasons)


def test_missing_data_all_present():
    """All fields present → zero score."""
    edge = _make_edge(u="N1", v="N2")
    samples = [
        {"site_id": "N1", "hydraulic_head": 100.0, "elevation": 50.0,
         "d18O": -5.0, "d2H": -30.0, "Cl": 10.0,
         "age": 50.0, "lithology": "sand", "screen_interval": "10-20",
         "aquifer_layer": "A1"},
        {"site_id": "N2", "hydraulic_head": 95.0, "elevation": 45.0,
         "d18O": -4.5, "d2H": -28.0, "Cl": 12.0,
         "age": 45.0, "lithology": "clay", "screen_interval": "15-25",
         "aquifer_layer": "A1"},
    ]
    config = HConfig()
    # Align config field keys with sample dict keys
    config.edge_head_key = "hydraulic_head"
    config.edge_elevation_key = "elevation"
    config.isotope_d18o_key = "d18O"
    config.isotope_d2h_key = "d2H"
    score, reasons = _score_missing_data(edge, samples, config)
    assert score == 0.0
    assert reasons == []


# ── manuscript_safe_note tests ────────────────────────────────────────

def test_no_recommendation_claims_validated():
    """VALIDATED edge must have manuscript-safe note mentioning validated status."""
    note = _manuscript_safe_note("VALIDATED", "observed_present")
    assert "VALIDATED" in note
    assert "independent confirmation" in note.lower()


def test_falsified_manuscript_note():
    """FALSIFIED edge recommends testing alternative hypotheses."""
    note = _manuscript_safe_note("FALSIFIED", "observed_absent")
    assert "falsifies" in note.lower()
    assert "alternative hypotheses" in note.lower()


# ── rank_next_measurements integration tests ─────────────────────────

def test_top_k_truncates():
    """With many candidate edges, top_k limits results."""
    # Build many edges via variant selected sets
    variants = {
        "baseline": {"selected_edge_ids": [f"E{i}" for i in range(50)]},
        "null_model_defaults": {"selected_edge_ids": []},
        "assumption_calibrated": {"selected_edge_ids": [f"E{i}" for i in range(0, 50, 2)]},
    }
    report = _make_benchmark_report(variants=variants)

    result = rank_next_measurements(benchmark_report=report, top_k=5)
    assert result["summary"]["n_recommendations"] == 5
    assert len(result["recommendations"]) == 5
    assert result["recommendations"][0]["rank"] == 1
    assert result["recommendations"][4]["rank"] == 5


def test_output_files_written():
    """When output_dir provided, JSON/CSV/MD files are written."""
    with _make_temp_dir() as tmp:
        out_dir = str(Path(tmp) / "out")
        variants = {
            "baseline": {"selected_edge_ids": ["E1", "E2"]},
            "null_model_defaults": {"selected_edge_ids": ["E2"]},
            "assumption_calibrated": {"selected_edge_ids": ["E1"]},
        }
        report = _make_benchmark_report(variants=variants)

        result = rank_next_measurements(
            benchmark_report=report,
            top_k=10,
            output_dir=out_dir,
        )

        assert (Path(out_dir) / "active_learning_recommendations.json").exists()
        assert (Path(out_dir) / "active_learning_recommendations.csv").exists()
        assert (Path(out_dir) / "active_learning_report.md").exists()

        # JSON is valid and matches result
        with open(Path(out_dir) / "active_learning_recommendations.json") as f:
            loaded = json.load(f)
        assert loaded["summary"]["n_recommendations"] == result["summary"]["n_recommendations"]

        # MD report has expected sections
        md_text = (Path(out_dir) / "active_learning_report.md").read_text()
        assert "# Active Learning" in md_text
        assert "## Top Recommendations" in md_text
        assert "## Recommendation Details" in md_text


def test_active_learning_requires_benchmark():
    """Empty benchmark report → graceful handling with warning."""
    result = rank_next_measurements(benchmark_report={})
    assert result["recommendations"] == []
    assert result["summary"]["n_recommendations"] == 0
    assert "warning" in result["summary"]


def test_rank_next_measurements_empty_inputs():
    """All optional inputs None → valid empty structure."""
    result = rank_next_measurements(benchmark_report=_make_benchmark_report())
    assert "recommendations" in result
    assert "summary" in result
    assert "inputs_used" in result
    assert result["summary"]["n_recommendations"] == 0


def test_priority_score_range():
    """All priority scores must be in [0, 1]."""
    # Create edges with various evidence flags and selection patterns
    edges = [
        _make_edge("E_A", "A", "B", {"evidence_flags": "age_reversal"}),
        _make_edge("E_B", "B", "C", {"evidence_flags": "missing_isotopes"}),
        _make_edge("E_C", "C", "D"),
        _make_edge("E_D", "D", "E", {"evidence_class": "VALIDATED"}),
    ]
    variants = {
        "baseline": {"selected_edge_ids": ["E_A", "E_B"]},
        "null_model_defaults": {"selected_edge_ids": ["E_A", "E_C"]},
        "assumption_calibrated": {"selected_edge_ids": ["E_A", "E_C", "E_D"]},
    }
    report = _make_benchmark_report(
        variants=variants,
        uncertainty={
            "improvement_summary_ci": {
                "delta_f1": [0.1, 0.7],
                "probability_delta_f1_gt_0": 0.88,
            }
        },
    )

    result = rank_next_measurements(
        benchmark_report=report,
        candidate_edges=edges,
        top_k=20,
    )

    for rec in result["recommendations"]:
        ps = rec["priority_score"]
        assert 0.0 <= ps <= 1.0, f"Priority score {ps} out of range for {rec['edge_id']}"


def test_priority_score_with_validation_labels(tmp_path):
    """Integration: edges with validation errors rank higher."""
    import pandas as pd

    # Write a validation labels CSV
    val_csv = tmp_path / "val_labels.csv"
    rows = [
        {"edge_id": "E_A", "observed_present": 1, "weight": 1.0},
        {"edge_id": "E_B", "observed_present": 0, "weight": 1.0},
        {"edge_id": "E_C", "observed_present": 1, "weight": 1.0},
    ]
    pd.DataFrame(rows).to_csv(val_csv, index=False)

    validation_report = _make_validation_report(str(val_csv))

    variants = {
        "baseline": {"selected_edge_ids": ["E_A", "E_B"]},
        "null_model_defaults": {"selected_edge_ids": ["E_A"]},
        "assumption_calibrated": {"selected_edge_ids": ["E_A", "E_B", "E_D"]},
    }
    report = _make_benchmark_report(variants=variants)

    result = rank_next_measurements(
        benchmark_report=report,
        validation_report=validation_report,
        top_k=10,
    )

    scores = {rec["edge_id"]: rec["priority_score"] for rec in result["recommendations"]}

    # E_B is FP (selected but absent=0) → highest score
    # E_D is selected_unlabeled (not in validation) → middle
    # E_A is TP → lowest among selected
    assert scores.get("E_B", 0.0) > scores.get("E_A", 0.0), (
        f"FP edge E_B should rank above TP edge E_A: {scores}"
    )


def test_recommendations_include_manuscript_safe_note():
    """Every recommendation must have a manuscript_safe_note."""
    variants = {
        "baseline": {"selected_edge_ids": ["E1"]},
        "null_model_defaults": {"selected_edge_ids": ["E1"]},
        "assumption_calibrated": {"selected_edge_ids": ["E1", "E2"]},
    }
    edges = [
        _make_edge("E2", "A", "B", {"evidence_class": "VALIDATED"}),
    ]
    report = _make_benchmark_report(variants=variants)

    result = rank_next_measurements(
        benchmark_report=report,
        candidate_edges=edges,
        top_k=10,
    )

    for rec in result["recommendations"]:
        assert "manuscript_safe_note" in rec
        assert isinstance(rec["manuscript_safe_note"], str)
        assert len(rec["manuscript_safe_note"]) > 0

    # E2 has evidence_class VALIDATED → note must mention it
    e2 = next(r for r in result["recommendations"] if r["edge_id"] == "E2")
    assert "VALIDATED" in e2["manuscript_safe_note"]


def test_bootstrap_instability_modulation_in_ranking():
    """Bootstrap instability from benchmark propagates to edge rankings."""
    variants = {
        "baseline": {"selected_edge_ids": ["E1", "E2"]},
        "null_model_defaults": {"selected_edge_ids": ["E2"]},
        "assumption_calibrated": {"selected_edge_ids": ["E1"]},
    }

    # Without bootstrap instability
    report_stable = _make_benchmark_report(variants=variants)
    result_stable = rank_next_measurements(
        benchmark_report=report_stable, top_k=10,
    )

    # With bootstrap instability (wide delta CI)
    report_unstable = _make_benchmark_report(
        variants=variants,
        uncertainty={
            "improvement_summary_ci": {
                "delta_f1": [0.05, 0.75],
                "probability_delta_f1_gt_0": 0.80,
            }
        },
    )
    result_unstable = rank_next_measurements(
        benchmark_report=report_unstable, top_k=10,
    )

    # The unstable report should produce higher or equal priority scores
    # (bootstrap instability is a positive modulation, never negative)
    for rec_s, rec_u in zip(result_stable["recommendations"], result_unstable["recommendations"]):
        assert rec_u["priority_score"] >= rec_s["priority_score"], (
            f"Bootstrap instability should not decrease score for {rec_s['edge_id']}"
        )


def test_md_report_empty_recommendations(tmp_path):
    """MD report for empty recommendations should still be valid."""
    out_dir = str(tmp_path / "out")
    result = rank_next_measurements(
        benchmark_report=_make_benchmark_report(),
        output_dir=out_dir,
    )

    md_path = Path(out_dir) / "active_learning_report.md"
    md_text = md_path.read_text()
    assert "No recommendations generated" in md_text


def test_missing_data_reasons_map_to_measurements():
    """Reasons with node suffixes from _score_missing_data must map to
    actionable measurements via _recommended_measurements."""
    reasons = ["missing_age (node A)", "missing_isotopes (node B)"]
    recs = _recommended_measurements(reasons)
    assert any("groundwater age tracer" in r for r in recs), (
        f"missing_age should map to age tracer, got: {recs}"
    )
    assert any("d18O/d2H sampling" in r for r in recs), (
        f"missing_isotopes should map to d18O/d2H, got: {recs}"
    )


def test_csv_columns_match(tmp_path):
    """CSV output must contain expected columns."""
    variants = {
        "baseline": {"selected_edge_ids": ["E1"]},
        "null_model_defaults": {"selected_edge_ids": []},
        "assumption_calibrated": {"selected_edge_ids": ["E1"]},
    }
    report = _make_benchmark_report(variants=variants)
    out_dir = str(tmp_path / "out")

    rank_next_measurements(benchmark_report=report, output_dir=out_dir)

    csv_path = Path(out_dir) / "active_learning_recommendations.csv"
    with open(csv_path) as f:
        header = f.readline().strip().split(",")
    for col in ("rank", "edge_id", "u", "v", "priority_score",
                "uncertainty_reasons", "recommended_measurements",
                "expected_benefit", "validation_status"):
        assert col in header, f"Missing CSV column: {col}"
