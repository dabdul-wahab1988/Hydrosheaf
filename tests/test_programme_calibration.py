from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    AgeCalibrationObservation,
    apply_age_interval_calibrator,
    fit_age_interval_calibrator,
    score_calibrated_age_intervals,
    score_selective_risk,
)


def test_calibrator_is_development_only_and_expands_intervals() -> None:
    observations = [
        AgeCalibrationObservation("A", truth=10.0, estimate=10.0, lower=8.0, upper=12.0),
        AgeCalibrationObservation("B", truth=20.0, estimate=20.0, lower=19.0, upper=21.0),
        AgeCalibrationObservation("C", truth=30.0, estimate=30.0, lower=29.0, upper=31.0),
    ]
    calibrator = fit_age_interval_calibrator(observations, target_coverage=0.8)

    assert calibrator.fit_scope == "development_only"
    assert calibrator.additive_radius == 0.0
    applied = apply_age_interval_calibrator(
        calibrator,
        {"D": {"mean_age_years": 40.0, "age_95_low": 35.0, "age_95_high": 45.0}},
    )
    assert applied["D"]["calibration_status"] == "development_fitted"
    assert applied["D"]["calibrated_age_low"] == 35.0

    with pytest.raises(ValueError, match="development phase"):
        fit_age_interval_calibrator(observations, phase="locked_test")


def test_calibration_score_uses_only_post_fit_truth() -> None:
    calibrator = fit_age_interval_calibrator(
        [
            AgeCalibrationObservation("A", truth=10.0, estimate=10.0, lower=9.0, upper=11.0),
            AgeCalibrationObservation("B", truth=20.0, estimate=20.0, lower=19.0, upper=21.0),
        ],
        target_coverage=0.95,
    )
    predictions = apply_age_interval_calibrator(
        calibrator,
        {
            "T1": {"mean_age_years": 12.0, "age_95_low": 10.0, "age_95_high": 14.0},
            "T2": {"mean_age_years": 25.0, "age_95_low": 20.0, "age_95_high": 30.0},
        },
    )
    score = score_calibrated_age_intervals({"T1": 12.0, "T2": 40.0}, predictions)

    assert score["status"] == "scored"
    assert score["n"] == 2
    assert 0.0 <= float(score["coverage"]) <= 1.0


def test_selective_risk_accepts_least_uncertain_rows_first() -> None:
    curve = score_selective_risk(
        [
            {"target_id": "certain", "truth": 10.0, "estimate": 10.0, "lower": 8.0, "upper": 12.0, "uncertainty": 1.0},
            {"target_id": "uncertain", "truth": 30.0, "estimate": 10.0, "lower": 0.0, "upper": 20.0, "uncertainty": 5.0},
        ],
        acceptance_fractions=(0.5, 1.0),
    )

    assert curve[0].accepted == 1
    assert curve[0].mean_absolute_error == 0.0
    assert curve[0].abstention_rate == 0.5
    assert curve[1].mean_absolute_error == 10.0
