"""Tests for the controlled synthetic age/directness benchmark."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from age_bayes_factor_benchmark import (  # noqa: E402
    CASE_STRATA,
    IDENTIFIABILITY_STRATA,
    METHODS,
    evaluate_predictions,
    generate_age_case,
    generate_cases,
    run_benchmark,
    score_case,
)


def test_inference_rows_keep_truth_in_a_separate_ledger() -> None:
    case = generate_age_case(
        2026090801,
        case_stratum="separable_transport",
        n_paths=2,
        nodes_per_path=6,
    )

    rows = case.inference_rows()
    assert rows
    assert len(rows) == len(case.truth.relation_by_candidate)
    assert {row["reference_type"] for row in rows} == {"independent_synthetic_truth"}
    assert not any(
        key in row
        for row in rows
        for key in ("is_direct", "truth_label", "relation", "direct_target")
    )
    assert {case.truth.relation(candidate.candidate_id) for candidate in case.candidates} == {
        "direct_adjacent",
        "indirect_reachable",
    }


def test_generation_is_deterministic_and_transport_is_edge_specific() -> None:
    first = generate_age_case(2026090802, case_stratum="mixed_transport")
    second = generate_age_case(2026090802, case_stratum="mixed_transport")

    assert first == second
    assert len({candidate.direct_travel_years for candidate in first.candidates}) > 1
    assert len({candidate.indirect_travel_years for candidate in first.candidates}) > 1
    assert set(candidate.identifiability_stratum for candidate in first.candidates) <= set(
        IDENTIFIABILITY_STRATA
    )


def test_controls_expose_ordering_failure_and_transport_pairing_signal() -> None:
    cases = generate_cases(n_cases=12)
    heldout = tuple(case for case in cases if case.split == "locked_test")
    predictions = [
        prediction
        for case in heldout
        for method in METHODS
        for prediction in score_case(case, method=method)
    ]
    metrics = evaluate_predictions(heldout, predictions)["methods"]

    assert metrics["no_age"]["n_scored"] == 0
    assert metrics["no_age"]["abstain_rate"] == pytest.approx(1.0)
    assert metrics["order_only"]["n_scored"] == metrics["order_only"]["n"]
    assert metrics["full_bf"]["n_scored"] > 0
    assert metrics["full_bf"]["directness_pr_auc"] > metrics["order_only"]["directness_pr_auc"]
    assert metrics["full_bf"]["directness_pr_auc"] > metrics["permuted_full_bf"]["directness_pr_auc"]
    assert metrics["full_bf"]["bf_sign_accuracy"] > metrics["permuted_full_bf"]["bf_sign_accuracy"]
    assert metrics["full_bf"]["identifiability_strata"]["overlapping"]["n_abstain"] > 0
    assert metrics["full_bf"]["case_strata"]["overlapping_transport"]["n_abstain"] > 0


def test_permuted_control_changes_transport_pairing_without_changing_ages() -> None:
    case = generate_age_case(
        2026090803,
        case_stratum="separable_transport",
        n_paths=2,
        nodes_per_path=6,
    )
    original = {row["candidate_id"]: row for row in case.inference_rows()}
    permuted = score_case(case, method="permuted_full_bf")

    assert all("truth_label" not in row for row in permuted)
    assert all(row["candidate_id"] in original for row in permuted)
    assert all(
        row["upstream_age_years"] == original[row["candidate_id"]]["upstream_age_years"]
        and row["downstream_age_years"] == original[row["candidate_id"]]["downstream_age_years"]
        for row in permuted
    )
    assert all("transport_source_candidate_id" in row for row in permuted)
    assert any(
        row["transport_source_candidate_id"] != row["candidate_id"] for row in permuted
    )


def test_run_writes_separate_truth_predictions_and_hashed_manifest(tmp_path: Path) -> None:
    manifest = run_benchmark(output=tmp_path / "age_bf_run", n_cases=6, n_development=3)
    output = tmp_path / "age_bf_run"

    assert set(manifest["methods"]) == set(METHODS)
    assert manifest["truth"]["inference_file_excludes_truth_labels"] is True
    assert set(manifest["outputs"]) == {
        "heldout_inference_candidates",
        "heldout_truth",
        "heldout_predictions",
        "complete_truth",
        "metrics",
    }
    assert set(manifest["output_sha256"]) == {
        "heldout_inference_candidates.csv",
        "heldout_truth.csv",
        "heldout_predictions.csv",
        "complete_truth.csv",
        "package.json",
        "truth_artifact.json",
        "qa.json",
        "metrics.json",
    }
    assert (output / "manifest.json").is_file()
    assert (output / "metrics.json").is_file()

    with (output / "heldout_inference_candidates.csv").open(newline="", encoding="utf-8") as handle:
        inference_fields = set(next(csv.reader(handle)))
    with (output / "heldout_truth.csv").open(newline="", encoding="utf-8") as handle:
        truth_fields = set(next(csv.reader(handle)))
    assert "truth_label" not in inference_fields
    assert "truth_label" in truth_fields

    saved_manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert saved_manifest["protocol"] == manifest["protocol"]
    assert saved_manifest["output_sha256"] == manifest["output_sha256"]


def test_case_schedule_covers_all_case_strata() -> None:
    cases = generate_cases(n_cases=24)
    assert {case.case_stratum for case in cases} == set(CASE_STRATA)
    assert {case.split for case in cases} == {"development", "locked_test"}
