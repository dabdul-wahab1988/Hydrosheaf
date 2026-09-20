"""Focused tests for the independent static graph-TTD virtual benchmark."""

from __future__ import annotations

import inspect
import json

import numpy as np

from hydrosheaf.benchmarks.ttd_graph_static import (
    TRUE_EDGES,
    generate_static_ttd_graph_case,
    run_all_pre_registered_controls,
    run_candidate_graph_estimator,
    run_local_baseline,
    score_static_ttd_graph_submission,
    write_static_ttd_graph_case,
)


def test_static_case_is_deterministic_sealed_and_manifested():
    first = generate_static_ttd_graph_case(seed=17)
    second = generate_static_ttd_graph_case(seed=17)

    assert first.observations.to_dict() == second.observations.to_dict()
    assert first.truth.to_dict() == second.truth.to_dict()
    assert first.manifest == second.manifest
    assert first.verify_manifest()
    assert (
        first.manifest["truth_release_policy"] == "sealed_until_submission_serialized"
    )
    assert first.manifest["truth_blind_reference_estimator"] is True

    public_payload = first.observations.to_dict()
    assert "true_edges" not in public_payload
    assert "candidate_controls" not in public_payload
    assert "edge_kernels" not in public_payload
    assert "node_ttd_masses" not in public_payload
    assert "held_out_values" not in public_payload
    assert first.truth.to_dict()["true_edges"] == [list(edge) for edge in TRUE_EDGES]


def test_truth_has_branch_merge_local_recharge_and_mass_conserving_serial_ttds():
    case = generate_static_ttd_graph_case(seed=31)
    truth = case.truth
    protocol = case.observations.protocol

    assert truth.true_edges == TRUE_EDGES
    assert set(("A", "B", "M", "O")).issubset(truth.node_mixing_weights)
    assert all(
        float(truth.node_mixing_weights[node]["local"]) > 0.0
        for node in ("A", "B", "M", "O")
    )
    assert all(
        np.isclose(sum(weights.values()), 1.0)
        for weights in truth.node_mixing_weights.values()
    )
    assert all(
        np.isclose(sum(kernel.masses), 1.0) and min(kernel.masses) >= 0.0
        for kernel in truth.edge_kernels
    )
    assert all(
        np.isclose(sum(masses), 1.0) and min(masses) >= 0.0
        for masses in truth.node_ttd_masses.values()
    )

    delta = np.zeros((len(protocol.age_grid_days),), dtype=float)
    delta[0] = 1.0
    kernel_ra = next(
        kernel
        for kernel in truth.edge_kernels
        if (kernel.source, kernel.target) == ("R", "A")
    )
    expected_a = (
        0.18 * delta
        + 0.82 * np.convolve(delta, np.asarray(kernel_ra.masses))[: delta.size]
    )
    assert np.allclose(np.asarray(truth.node_ttd_masses["A"]), expected_a)
    assert (
        truth.young_water_fractions["A"]
        > truth.young_water_fractions["M"]
        > truth.young_water_fractions["O"]
    )


def test_public_data_are_irregular_missing_and_have_truth_free_held_out_targets():
    case = generate_static_ttd_graph_case(seed=22)
    observations = case.observations

    assert observations.held_out_targets
    assert any(
        not row.observed and row.value is None for row in observations.input_rows
    )
    assert any(
        not row.observed and row.value is None for row in observations.calibration_rows
    )
    input_times = [
        row.time_day
        for row in observations.input_rows
        if row.node_id == "A" and row.tracer == "delta18O"
    ]
    assert any(right - left > 1 for left, right in zip(input_times, input_times[1:]))
    public_text = json.dumps(observations.to_dict(), sort_keys=True)
    assert "held_out_values" not in public_text
    assert "node_output_signals" not in public_text


def test_truth_blind_controls_produce_submissions_and_sealed_scores():
    case = generate_static_ttd_graph_case(seed=17)
    observations = case.observations
    signature = inspect.signature(run_candidate_graph_estimator)

    assert tuple(signature.parameters) == ("observations", "control")
    submissions = run_all_pre_registered_controls(observations, case.controls)
    assert set(submissions) == {
        "correct",
        "reversed",
        "random",
        "edge_removed",
        "local",
    }
    assert run_local_baseline(observations).candidate_edges == ()
    assert submissions["correct"].candidate_edges == TRUE_EDGES
    assert len(submissions["correct"].held_out_predictions) > 0
    assert any(
        interval.status == "ESTIMATED"
        for interval in submissions["correct"].node_intervals.values()
    )

    correct_score = score_static_ttd_graph_submission(
        case.truth, observations, submissions["correct"]
    )
    local_score = score_static_ttd_graph_submission(
        case.truth, observations, submissions["local"]
    )
    removed_score = score_static_ttd_graph_submission(
        case.truth, observations, submissions["edge_removed"]
    )
    reversed_score = score_static_ttd_graph_submission(
        case.truth, observations, submissions["reversed"]
    )

    assert correct_score.metrics["topology"]["candidate_graph"]["f1"] == 1.0
    assert removed_score.metrics["topology"]["candidate_graph"]["recall"] == 0.8
    assert reversed_score.metrics["topology"]["candidate_graph"]["true_positive"] == 0
    assert (
        correct_score.metrics["held_out_prediction"]["mae"]
        < local_score.metrics["held_out_prediction"]["mae"]
    )
    assert "conditional_coverage" in correct_score.metrics["interval_recovery"]
    assert "abstention_rate" in local_score.metrics["interval_recovery"]
    assert (
        correct_score.to_dict()["metrics"]["provenance"]["truth_used_by_submission"]
        is False
    )


def test_serialized_artifacts_keep_public_and_sealed_payloads_separate(tmp_path):
    case = generate_static_ttd_graph_case(seed=41)
    paths = write_static_ttd_graph_case(case, tmp_path)

    assert set(paths) == {"observations", "sealed_truth", "control_plan", "manifest"}
    public = json.loads(paths["observations"].read_text(encoding="utf-8"))
    sealed = json.loads(paths["sealed_truth"].read_text(encoding="utf-8"))
    controls = json.loads(paths["control_plan"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    assert "edge_kernels" not in public
    assert "edge_kernels" in sealed
    assert "candidate_controls" not in public
    assert controls["oracle_control_note"]
    assert manifest["observations_sha256"] != manifest["sealed_truth_sha256"]
    assert manifest["generator_family"] == "independent_analytic_source_mixing_network"
