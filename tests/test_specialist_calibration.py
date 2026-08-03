from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    SpecialistAgeCalibrationObservation,
    SpecialistReactionCalibrationObservation,
    fit_specialist_age_calibrator,
    fit_specialist_reaction_calibrator,
    score_calibrated_specialist_reaction,
)


def test_age_calibration_is_case_blocked_and_truth_free_at_apply() -> None:
    observations = [
        SpecialistAgeCalibrationObservation("case-a", "a1", 20.0, 15.0, 10.0, 20.0),
        SpecialistAgeCalibrationObservation("case-a", "a2", 30.0, 25.0, 20.0, 30.0),
        SpecialistAgeCalibrationObservation("case-b", "b1", 40.0, 35.0, 30.0, 40.0),
        SpecialistAgeCalibrationObservation("case-b", "b2", 50.0, 45.0, 40.0, 50.0),
    ]
    calibrator = fit_specialist_age_calibrator(
        observations,
        target_coverage=0.95,
        phase="development",
    )

    assert calibrator.fit_scope == "development_only"
    assert calibrator.case_ids == ("case-a", "case-b")
    assert calibrator.bias_offset_years == pytest.approx(5.0)
    applied = calibrator.apply(
        {"locked-node": {"mean_age_years": 10.0, "age_95_low": 5.0, "age_95_high": 15.0}}
    )
    assert applied["locked-node"]["calibration_status"] == "development_fitted"
    assert applied["locked-node"]["calibrated_mean_age_years"] == pytest.approx(15.0)
    assert "truth" not in applied["locked-node"]

    with pytest.raises(ValueError, match="development data"):
        fit_specialist_age_calibrator(observations, phase="locked_test")


def test_age_calibration_selection_grid_is_fitted_only_on_development_truth() -> None:
    observations = [
        SpecialistAgeCalibrationObservation(
            "case-a", "a1", 20.0, 20.0, 19.0, 21.0, selection_score=1.0
        ),
        SpecialistAgeCalibrationObservation(
            "case-a", "a2", 22.0, 22.0, 21.0, 23.0, selection_score=2.0
        ),
        SpecialistAgeCalibrationObservation(
            "case-b", "b1", 40.0, 40.0, 39.0, 41.0, selection_score=10.0
        ),
        SpecialistAgeCalibrationObservation(
            "case-b", "b2", 42.0, 42.0, 41.0, 43.0, selection_score=11.0
        ),
    ]
    calibrator = fit_specialist_age_calibrator(
        observations,
        target_coverage=0.75,
        minimum_selection_rate=0.50,
        interval_scales=(0.0, 1.0),
        phase="development",
    )

    assert calibrator.fit_scope == "development_only"
    assert calibrator.calibration_kind.endswith("selection_grid")
    assert calibrator.selection_threshold == pytest.approx(11.0)
    assert calibrator.interval_scale == pytest.approx(0.0)
    applied = calibrator.apply(
        {
            "selected": {
                "mean_age_years": 30.0,
                "age_95_low": 20.0,
                "age_95_high": 40.0,
                "tracer_age_spread_years": 5.0,
            },
            "abstained": {
                "mean_age_years": 30.0,
                "age_95_low": 20.0,
                "age_95_high": 40.0,
                "tracer_age_spread_years": 20.0,
            },
        }
    )
    assert applied["selected"]["decision"] == "select"
    assert applied["abstained"]["decision"] == "abstain"
    assert "truth" not in applied["selected"]


def test_age_calibration_can_select_all_outputs_when_raw_abstentions_block_coverage() -> None:
    observations = [
        SpecialistAgeCalibrationObservation(
            "case-a", "a1", 20.0, 20.0, 19.0, 21.0, selection_score=1.0, selected=True
        ),
        SpecialistAgeCalibrationObservation(
            "case-b", "b1", 40.0, 40.0, 39.0, 41.0, selected=False
        ),
    ]
    calibrator = fit_specialist_age_calibrator(
        observations,
        target_coverage=0.95,
        minimum_selection_rate=0.95,
        interval_scales=(1.0,),
        phase="development",
    )
    assert calibrator.selection_all is True
    applied = calibrator.apply(
        {
            "locked": {
                "mean_age_years": 40.0,
                "age_95_low": 39.0,
                "age_95_high": 41.0,
                "decision": "abstain",
                "tracer_age_spread_years": 100.0,
            }
        }
    )
    assert applied["locked"]["decision"] == "select"


def test_reaction_temperature_scaling_calibrates_multiclass_probabilities() -> None:
    observations = [
        SpecialistReactionCalibrationObservation(
            "case-a", "a1", "denitrification", {"denitrification": 2.0, "none": 0.0}
        ),
        SpecialistReactionCalibrationObservation(
            "case-a", "a2", "none", {"denitrification": 0.0, "none": 2.0}
        ),
        SpecialistReactionCalibrationObservation(
            "case-b", "b1", "denitrification", {"denitrification": 1.5, "none": 0.0}
        ),
        SpecialistReactionCalibrationObservation(
            "case-b", "b2", "none", {"denitrification": 0.0, "none": 1.5}
        ),
    ]
    calibrator = fit_specialist_reaction_calibrator(observations, phase="development")
    assert calibrator.fit_scope == "development_only"
    assert calibrator.temperature > 0.0
    assert "grid_boundary_hit" in calibrator.provenance
    applied = calibrator.apply(
        {
            "locked-edge": {
                "family": "denitrification",
                "logits": {"denitrification": 1.0, "none": 0.0},
                "probabilities": {"denitrification": 0.73, "none": 0.27},
                "decision": "select",
            }
        }
    )
    record = applied["locked-edge"]
    assert sum(record["probabilities"].values()) == pytest.approx(1.0)
    assert record["calibration_status"] == "development_fitted"
    score = score_calibrated_specialist_reaction(
        {"locked-edge": "denitrification"}, applied
    )
    assert score["status"] == "scored"
    assert score["n"] == 1
    assert score["multiclass_log_loss"] >= 0.0
    missing = calibrator.apply({"missing-edge": {"decision": "abstain"}})
    assert missing["missing-edge"]["decision"] == "abstain"
    assert missing["missing-edge"]["calibration_status"] == "abstained_not_calibrated"


def test_specialist_calibration_rejects_empty_fit_data() -> None:
    with pytest.raises(ValueError, match="At least one"):
        fit_specialist_reaction_calibrator([], phase="development")
