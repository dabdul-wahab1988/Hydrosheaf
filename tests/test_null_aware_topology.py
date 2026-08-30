from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.inference.null_aware import (
    CALIBRATION_FITTED_MISSING_EVIDENCE,
    CALIBRATION_UNFITTED,
    DECISION_ABSTAIN,
    DECISION_ABSENT,
    DECISION_PRESENT,
    CalibrationNotFittedError,
    NullAwareFeatureRow,
    NullAwareLogisticCalibrator,
    NullAwareTopologyScorer,
    build_feature_rows,
    calibration_diagnostics,
)


def _training_rows() -> tuple[list[NullAwareFeatureRow], np.ndarray]:
    rows = [
        NullAwareFeatureRow(
            "flow-1",
            flow_features={"head_drop_m": 8.0, "age_direction_support": 0.9},
            null_features={"null_chemistry_similarity": 0.05, "common_lithology": 0.0},
        ),
        NullAwareFeatureRow(
            "flow-2",
            flow_features={"head_drop_m": 6.0, "age_direction_support": 0.8},
            null_features={"null_chemistry_similarity": 0.10, "common_lithology": 0.0},
        ),
        NullAwareFeatureRow(
            "flow-3",
            flow_features={"head_drop_m": 5.0, "age_direction_support": 0.7},
            null_features={"null_chemistry_similarity": 0.20, "common_lithology": 0.0},
        ),
        NullAwareFeatureRow(
            "null-1",
            flow_features={"head_drop_m": 8.0, "age_direction_support": 0.9},
            null_features={"null_chemistry_similarity": 0.95, "common_lithology": 1.0},
        ),
        NullAwareFeatureRow(
            "null-2",
            flow_features={"head_drop_m": 6.0, "age_direction_support": 0.8},
            null_features={"null_chemistry_similarity": 0.85, "common_lithology": 1.0},
        ),
        NullAwareFeatureRow(
            "null-3",
            flow_features={"head_drop_m": 1.0, "age_direction_support": 0.2},
            null_features={"null_chemistry_similarity": 0.90, "common_lithology": 1.0},
        ),
    ]
    return rows, np.asarray([1, 1, 1, 0, 0, 0], dtype=float)


def test_no_flow_confounding_is_separated_from_flow_support() -> None:
    rows, labels = _training_rows()
    calibrator = NullAwareLogisticCalibrator(l2=0.1).fit(rows, labels)
    candidates = [
        NullAwareFeatureRow(
            "supported",
            flow_features={"head_drop_m": 7.0, "age_direction_support": 0.85},
            null_features={"null_chemistry_similarity": 0.08, "common_lithology": 0.0},
        ),
        NullAwareFeatureRow(
            "confounded",
            flow_features={"head_drop_m": 7.0, "age_direction_support": 0.85},
            null_features={"null_chemistry_similarity": 0.92, "common_lithology": 1.0},
        ),
    ]

    records = NullAwareTopologyScorer(
        calibrator,
        present_threshold=0.65,
        absent_threshold=0.35,
    ).score_rows(candidates)

    assert records[0]["flow_probability"] > records[1]["flow_probability"]
    assert records[1]["null_score"] > records[0]["null_score"]
    assert records[0]["decision"] == DECISION_PRESENT
    assert records[1]["decision"] == DECISION_ABSENT


def test_truth_fields_do_not_change_features_or_predictions() -> None:
    edges = [
        {
            "edge_id": "e1",
            "u": "a",
            "v": "b",
            "attrs": {
                "flow_features": {"head_drop_m": 2.0},
                "null_features": {"null_chemistry_similarity": 0.2},
            },
        }
    ]
    observations = {
        "a": {"head": 100.0, "mean_age_years": 5.0, "truth": 999.0},
        "b": {"head": 98.0, "mean_age_years": 9.0, "ground_truth": -999.0},
    }
    rows_without_truth = build_feature_rows(edges, observations)
    edges_with_truth = [
        {
            **edges[0],
            "attrs": {
                **edges[0]["attrs"],
                "label": 0,
                "observed_present": 0,
                "is_true_flow": False,
                "reference_flow": "none",
                "arbitrary_numeric_truth": 1e9,
            },
        }
    ]
    observations_with_truth = {
        key: {**value, "target": 1, "actual": -1, "true_age_years": 1e6}
        for key, value in observations.items()
    }
    rows_with_truth = build_feature_rows(edges_with_truth, observations_with_truth)

    assert rows_without_truth[0].to_dict() == rows_with_truth[0].to_dict()

    rows, labels = _training_rows()
    calibrator = NullAwareLogisticCalibrator(l2=0.1).fit(rows, labels)
    assert np.array_equal(calibrator.predict(rows_without_truth), calibrator.predict(rows_with_truth))


def test_declared_observations_generate_null_features_without_truth() -> None:
    edges = [{"edge_id": "A_B", "u": "A", "v": "B"}]
    observations = {
        "A": {
            "Ca": 2.0,
            "Mg": 1.0,
            "18O": -5.0,
            "2H": -30.0,
            "aquifer_layer": "upper",
            "lat": 0.0,
            "lon": 0.0,
        },
        "B": {
            "Ca": 2.1,
            "Mg": 1.1,
            "18O": -5.1,
            "2H": -30.5,
            "aquifer_layer": "upper",
            "lat": 0.0,
            "lon": 0.005,
        },
    }

    row = build_feature_rows(edges, observations)[0]

    assert "null_chemistry_similarity" in row.null_features
    assert "null_isotope_similarity" in row.null_features
    assert "common_lithology" in row.null_features
    assert "spatial_proximity" in row.null_features
    assert "null_score" in row.null_features
    assert all("truth" not in field.lower() for field in row.source_fields)


def test_observation_head_drop_and_age_increment_use_explicit_sign_conventions() -> None:
    rows = build_feature_rows(
        [{"edge_id": "A_B", "u": "A", "v": "B"}],
        {
            "A": {"hydraulic_head": 100.0, "mean_age_years": 5.0},
            "B": {"hydraulic_head": 98.0, "mean_age_years": 9.0},
        },
    )

    assert rows[0].flow_features["head_drop_m"] == 2.0
    assert rows[0].flow_features["age_increment_years"] == 4.0


def test_missing_evidence_abstains_and_unfitted_calibrator_fails_closed() -> None:
    rows, labels = _training_rows()
    fitted = NullAwareLogisticCalibrator(l2=0.1).fit(rows, labels)
    incomplete = NullAwareFeatureRow(
        "incomplete",
        flow_features={"head_drop_m": 4.0},
        null_features={},
    )
    missing_record = NullAwareTopologyScorer(fitted).score_rows([incomplete])[0]
    assert missing_record["decision"] == DECISION_ABSTAIN
    assert missing_record["calibration_status"] == CALIBRATION_FITTED_MISSING_EVIDENCE
    assert any(channel.startswith("null:") for channel in missing_record["missing_channels"])
    assert missing_record["flow_probability"] is None

    unfitted = NullAwareLogisticCalibrator()
    with pytest.raises(CalibrationNotFittedError):
        unfitted.predict([rows[0]])
    unfitted_records = NullAwareTopologyScorer(unfitted).score_rows([rows[0], incomplete])
    assert len(unfitted_records) == 2
    assert all(record["decision"] == DECISION_ABSTAIN for record in unfitted_records)
    assert all(record["calibration_status"] == CALIBRATION_UNFITTED for record in unfitted_records)
    assert all(record["flow_probability"] is None for record in unfitted_records)


def test_fit_is_deterministic_and_serializable() -> None:
    rows, labels = _training_rows()
    first = NullAwareLogisticCalibrator(l2=0.1, max_iter=300).fit(rows, labels)
    second = NullAwareLogisticCalibrator(l2=0.1, max_iter=300).fit(rows, labels)

    assert first.to_dict() == second.to_dict()
    assert np.array_equal(first.predict(rows), second.predict(rows))
    assert first.to_dict()["calibration_status"] == "FITTED"
    assert "coefficients" in first.to_dict()
    restored = NullAwareLogisticCalibrator.from_dict(first.to_dict())
    assert np.array_equal(first.predict(rows), restored.predict(rows))


def test_probabilities_and_calibration_diagnostics_are_bounded() -> None:
    rows, labels = _training_rows()
    calibrator = NullAwareLogisticCalibrator(l2=0.1).fit(rows, labels)
    probabilities = calibrator.predict(rows)
    diagnostics = calibration_diagnostics(labels, probabilities, n_bins=5)

    assert np.all(np.isfinite(probabilities))
    assert np.all((probabilities >= 0.0) & (probabilities <= 1.0))
    assert 0.0 <= diagnostics["brier"] <= 1.0
    assert diagnostics["log_loss"] >= 0.0
    assert 0.0 <= diagnostics["ece"] <= 1.0
    assert diagnostics["n"] == len(rows)
    assert calibrator.diagnostics(rows, labels)["n"] == len(rows)
