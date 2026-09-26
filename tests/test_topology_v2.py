from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.validation.topology_v2 import (
    DECISION_ABSTAIN,
    DECISION_PRESENT,
    MODEL_FEATURE_SETS,
    TopologyV2Config,
    TopologyV2LogisticCalibrator,
    TopologyV2Scorer,
    TopologyV2ThresholdPolicy,
    build_topology_v2_feature_rows,
    bootstrap_case_metric_ci,
    generate_topology_v2_candidate_universe,
    topology_v2_metrics,
    tune_topology_thresholds,
)


def _observations() -> list[dict[str, object]]:
    return [
        {
            "site_id": "A",
            "x_m": 0.0,
            "y_m": 0.0,
            "head_meas": 100.0,
            "head_sigma_m": 0.1,
            "screen_depth": 10.0,
            "well_depth": 20.0,
            "aquifer_unit": "upper",
            "geology_symbol": "sand",
            "Ca": 1.0,
            "18O": -5.0,
        },
        {
            "site_id": "B",
            "x_m": 1000.0,
            "y_m": 0.0,
            "head_meas": 98.0,
            "head_sigma_m": 0.1,
            "screen_depth": 10.0,
            "well_depth": 20.0,
            "aquifer_unit": "upper",
            "geology_symbol": "sand",
            "Ca": 1.1,
            "18O": -5.1,
        },
        {
            "site_id": "C",
            "x_m": 2000.0,
            "y_m": 0.0,
            "head_meas": 96.0,
            "head_sigma_m": 0.1,
            "screen_depth": 10.0,
            "well_depth": 20.0,
            "aquifer_unit": "upper",
            "geology_symbol": "sand",
            "Ca": 1.2,
            "18O": -5.2,
        },
    ]


def _rows_and_labels():
    observations = _observations()
    universe = generate_topology_v2_candidate_universe(observations)
    rows = build_topology_v2_feature_rows(universe, observations, case_id="case-1")
    labels = np.asarray(
        [int(row.edge_id in {"A->B", "B->C"}) for row in rows],
        dtype=float,
    )
    return universe, rows, labels


def test_default_candidate_universe_is_all_directed_pairs_and_keeps_uphill_edges():
    observations = _observations()
    observations[1]["head_meas"] = 101.0
    universe = generate_topology_v2_candidate_universe(observations)

    assert len(universe.edges) == 6
    uphill = next(edge for edge in universe.edges if edge.edge_id == "A->B")
    assert uphill.attrs["direction_probability"] < 0.5
    assert universe.rejection_counts == {}


def test_distinct_nodes_with_coincident_coordinates_are_not_self_loop_rejections():
    observations = _observations()
    observations[1]["x_m"] = observations[0]["x_m"]
    observations[1]["y_m"] = observations[0]["y_m"]

    universe = generate_topology_v2_candidate_universe(observations)

    assert len(universe.edges) == 6
    assert "A->B" in universe.edge_ids
    assert "B->A" in universe.edge_ids


def test_hard_direction_is_explicit_and_audited():
    observations = _observations()
    observations[1]["head_meas"] = 101.0
    universe = generate_topology_v2_candidate_universe(
        observations,
        config=TopologyV2Config(hard_direction=True),
    )

    assert "A->B" not in universe.edge_ids
    assert universe.rejection_counts["hard_direction"] >= 1


def test_feature_signs_and_uncertainty_are_explicit():
    universe, rows, _labels = _rows_and_labels()
    row = next(row for row in rows if row.edge_id == "A->B")

    assert row.features["head_delta_m"] == pytest.approx(2.0)
    assert row.features["head_sigma_delta_m"] == pytest.approx(2**0.5 * 0.1)
    assert row.features["head_z_score"] > 0.0
    assert row.features["direction_probability"] > 0.5
    assert row.features["gradient_m_per_km"] == pytest.approx(2.0)
    assert universe.candidate_graph_recall([("A", "B"), ("B", "C")]) == 1.0


def test_potassium_chemistry_is_not_mistaken_for_hydraulic_conductivity():
    observations = _observations()
    observations[0]["K"] = 0.05
    observations[1]["K"] = 0.06
    universe = generate_topology_v2_candidate_universe(observations)
    rows = build_topology_v2_feature_rows(universe, observations)

    assert all(row.features["hydraulic_capacity_proxy"] is None for row in rows)


def test_truth_fields_are_rejected_before_candidate_generation():
    observations = _observations()
    observations[0]["heldout_truth"] = "A->B"

    with pytest.raises(ValueError, match="Truth/reference field"):
        generate_topology_v2_candidate_universe(observations)


def test_calibrator_is_deterministic_serializable_and_tri_state():
    _universe, rows, labels = _rows_and_labels()
    feature_names = MODEL_FEATURE_SETS["C_geometry_head_gradient"]
    first = TopologyV2LogisticCalibrator(feature_names=feature_names, l2=0.1).fit(
        rows,
        labels,
        scope="held_out_calibration",
        independent=True,
        generator_id="test-generator",
        split_id="development-case-1",
        dataset_hash="dataset-hash",
    )
    second = TopologyV2LogisticCalibrator(feature_names=feature_names, l2=0.1).fit(
        rows,
        labels,
        scope="held_out_calibration",
        independent=True,
        generator_id="test-generator",
        split_id="development-case-1",
        dataset_hash="dataset-hash",
    )

    assert first.deployment_eligible
    assert np.array_equal(first.predict_proba(rows), second.predict_proba(rows))
    restored = TopologyV2LogisticCalibrator.from_dict(first.to_dict())
    assert np.array_equal(first.predict_proba(rows), restored.predict_proba(rows))

    scorer = TopologyV2Scorer(
        restored,
        policy=TopologyV2ThresholdPolicy(present_threshold=0.55, absent_threshold=0.45),
    )
    records = scorer.score_rows(rows)
    assert all(record["probability"] is not None for record in records)
    assert any(record["decision"] in {DECISION_PRESENT, DECISION_ABSTAIN} for record in records)


def test_unfitted_or_non_deployable_calibration_fails_closed_to_abstain():
    _universe, rows, _labels = _rows_and_labels()
    calibrator = TopologyV2LogisticCalibrator(feature_names=MODEL_FEATURE_SETS["A_geometry"])
    records = TopologyV2Scorer(calibrator).score_rows(rows)

    assert records
    assert all(record["decision"] == DECISION_ABSTAIN for record in records)
    assert all(record["probability"] is None for record in records)


def test_zero_and_one_thresholds_are_invalid_for_selected_inference():
    with pytest.raises(ValueError):
        TopologyV2ThresholdPolicy(present_threshold=0.0, absent_threshold=0.0)
    with pytest.raises(ValueError):
        TopologyV2ThresholdPolicy(present_threshold=1.0, absent_threshold=0.5)


def test_metrics_report_pr_auc_calibration_and_abstention():
    policy = TopologyV2ThresholdPolicy(present_threshold=0.75, absent_threshold=0.25)
    metrics = topology_v2_metrics(
        [1, 1, 0, 0],
        [0.95, 0.55, 0.10, 0.45],
        policy=policy,
        candidate_recall=1.0,
    )

    assert metrics["pr_auc"] is not None
    assert 0.0 <= metrics["brier"] <= 1.0
    assert metrics["abstain_count"] == 2
    assert metrics["candidate_graph_recall"] == 1.0


def test_threshold_tuning_is_validation_only_and_uses_interior_grid():
    policy, audit = tune_topology_thresholds(
        [1, 1, 0, 0, 0, 1],
        [0.90, 0.72, 0.20, 0.30, 0.40, 0.60],
        target_fdr=0.50,
    )

    assert 0.0 < policy.absent_threshold < policy.present_threshold < 1.0
    assert audit
    assert all(
        0.0 < item["policy"]["absent_threshold"] < item["policy"]["present_threshold"] < 1.0
        for item in audit
    )


def test_case_bootstrap_declares_whole_case_resampling():
    result = bootstrap_case_metric_ci(
        [1, 0, 1, 0, 1, 0],
        [0.9, 0.1, 0.8, 0.2, 0.7, 0.3],
        ["a", "a", "b", "b", "c", "c"],
        n_bootstrap=50,
        seed=7,
    )

    assert result["resampling_unit"] == "whole_case"
    assert result["n_cases"] == 3
    assert result["ci95_low"] <= result["ci95_high"]
