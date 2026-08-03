from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    DiscrepancyCalibrationObservation,
    DecisionKind,
    IdentifiabilityStatus,
    ScenarioPrediction,
    apply_discrepancy_calibrator,
    apply_discrepancy_calibrator_to_record,
    build_discrepancy_report,
    fit_discrepancy_calibrator,
    score_locked_discrepancy,
)


def test_compatible_scenarios_return_a_set_report_and_conservative_envelope() -> None:
    report = build_discrepancy_report(
        "age_years",
        [
            ScenarioPrediction(
                "input_history_A",
                estimate=10.0,
                lower=7.0,
                upper=14.0,
                weight=2.0,
                scenario_family="age_forward_model",
            ),
            ScenarioPrediction(
                "input_history_B",
                estimate=11.0,
                lower=8.0,
                upper=15.0,
                weight=1.0,
                scenario_family="age_forward_model",
            ),
        ],
        compatible_hypotheses=("young_water", "mixed_young_water"),
    )

    assert report.identifiability is IdentifiabilityStatus.PARTIALLY_IDENTIFIED
    assert report.decision.decision_kind is DecisionKind.SET_REPORT
    assert report.estimate == pytest.approx(10.3333333333)
    assert report.lower == 7.0
    assert report.upper == 15.0
    assert report.interval_kind.endswith("not_calibrated")
    payload = report.to_dict()
    assert payload["calibrated"] is False
    assert payload["calibration_status"] == "NOT_CALIBRATED"
    assert payload["model_disagreement_status"] == "NO_MATERIAL_DISAGREEMENT"
    assert payload["scenario_families"] == ["age_forward_model"]
    assert sum(report.scenario_weights.values()) == pytest.approx(1.0)


def test_material_disagreement_forces_model_disagreement_abstention() -> None:
    report = build_discrepancy_report(
        "reaction_extent",
        [
            ScenarioPrediction(
                "nominal_topology",
                estimate=0.1,
                lower=0.0,
                upper=0.2,
                scenario_family="topology",
            ),
            ScenarioPrediction(
                "false_positive_edge",
                estimate=1.2,
                lower=1.0,
                upper=1.4,
                scenario_family="topology",
                discrepancy_tags=("false_positive_edge",),
            ),
        ],
        disagreement_threshold=0.25,
        scale_floor=1.0,
        compatible_hypotheses=("low_extent", "high_extent"),
    )

    assert report.identifiability is IdentifiabilityStatus.MODEL_DISAGREEMENT
    assert report.decision.decision_kind is DecisionKind.ABSTAIN
    assert report.decision.identifiability is IdentifiabilityStatus.MODEL_DISAGREEMENT
    assert report.decision.scenario_count == 2
    assert report.model_disagreement_status == "MODEL_DISAGREEMENT"
    assert report.identifiability_reason
    assert report.scenario_families == ("topology",)


def test_scenario_contract_rejects_invalid_or_ambiguous_inputs() -> None:
    with pytest.raises(ValueError, match="within its interval"):
        ScenarioPrediction(
            "invalid",
            estimate=3.0,
            lower=0.0,
            upper=2.0,
        )
    with pytest.raises(ValueError, match="unique"):
        build_discrepancy_report(
            "x",
            [
                ScenarioPrediction("same", 1.0, 0.0, 2.0),
                ScenarioPrediction("same", 1.0, 0.0, 2.0),
            ],
        )
    with pytest.raises(ValueError, match="positive weight"):
        build_discrepancy_report(
            "x",
            [ScenarioPrediction("zero", 1.0, 0.0, 2.0, weight=0.0)],
        )


def test_discrepancy_calibration_is_development_only_and_case_blocked() -> None:
    observations = [
        DiscrepancyCalibrationObservation(
            "case_a", "node_a", truth=12.0, estimate=10.0, lower=9.0, upper=11.0
        ),
        DiscrepancyCalibrationObservation(
            "case_b", "node_b", truth=10.5, estimate=10.0, lower=9.0, upper=11.0
        ),
    ]
    calibrator = fit_discrepancy_calibrator(observations, target_coverage=0.95)

    assert calibrator.fit_scope == "development_only"
    assert calibrator.calibration_kind.startswith("case_blocked")
    assert calibrator.scale_factor == pytest.approx(2.0)
    with pytest.raises(ValueError, match="development observations only"):
        fit_discrepancy_calibrator(
            [
                DiscrepancyCalibrationObservation(
                    "locked", "node", truth=10.0, estimate=10.0, lower=9.0, upper=11.0,
                    phase="locked_test",
                )
            ],
            phase="locked_test",
        )


def test_locked_discrepancy_scores_false_commitment_and_calibration_degradation() -> None:
    development = [
        DiscrepancyCalibrationObservation(
            "dev", "dev_node", truth=12.0, estimate=10.0, lower=9.0, upper=11.0
        )
    ]
    calibrator = fit_discrepancy_calibrator(development, target_coverage=0.95)
    raw_report = build_discrepancy_report(
        "age:node",
        [ScenarioPrediction("nominal", 10.0, 9.0, 11.0)],
        compatible_hypotheses=("young",),
    )
    calibrated_report = apply_discrepancy_calibrator(calibrator, raw_report)
    locked = [
        DiscrepancyCalibrationObservation(
            "locked", "locked|age:node", truth=12.0, estimate=10.0, lower=9.0, upper=11.0,
            phase="locked_test",
        )
    ]
    score = score_locked_discrepancy(
        locked,
        {"locked|age:node": raw_report},
        {"locked|age:node": calibrated_report},
        target_coverage=calibrator.target_coverage,
    )

    assert score["status"] == "scored"
    assert score["raw"]["false_commitment_rate"] == pytest.approx(1.0)
    assert score["calibrated"]["false_commitment_rate"] == pytest.approx(0.0)
    assert score["calibrated"]["coverage"] == pytest.approx(1.0)
    assert score["calibration_degradation"] < 0.0


def test_locked_discrepancy_missing_or_uncalibrated_records_abstain() -> None:
    locked = [
        DiscrepancyCalibrationObservation(
            "locked", "key", truth=1.0, estimate=1.0, lower=0.0, upper=2.0,
            phase="locked_test",
        )
    ]
    raw = build_discrepancy_report(
        "age:key", [ScenarioPrediction("nominal", 1.0, 0.0, 2.0)], compatible_hypotheses=("x",)
    )
    assert score_locked_discrepancy(
        locked, {"key": raw}, {"key": raw}, target_coverage=0.95
    )["status"] == "ABSTAIN_CALIBRATION_NOT_VERIFIED"
    assert score_locked_discrepancy(
        locked, {}, {}, target_coverage=0.95
    )["status"] == "ABSTAIN_MISSING_LOCKED_REPORTS"


def test_discrepancy_calibration_adapter_rejects_truth_bearing_records() -> None:
    calibrator = fit_discrepancy_calibrator(
        [
            DiscrepancyCalibrationObservation(
                "dev", "node", truth=12.0, estimate=10.0, lower=9.0, upper=11.0
            )
        ]
    )
    with pytest.raises(ValueError, match="Truth fields"):
        apply_discrepancy_calibrator_to_record(
            calibrator,
            {"estimate": 10.0, "lower": 9.0, "upper": 11.0, "truth": 12.0},
        )


def test_locked_discrepancy_scoring_abstains_on_truth_leakage() -> None:
    locked = [
        DiscrepancyCalibrationObservation(
            "locked", "key", truth=1.0, estimate=1.0, lower=0.0, upper=2.0,
            phase="locked_test",
        )
    ]
    raw = build_discrepancy_report(
        "age:key", [ScenarioPrediction("nominal", 1.0, 0.0, 2.0)], compatible_hypotheses=("x",)
    )
    calibrated = apply_discrepancy_calibrator(
        fit_discrepancy_calibrator(
            [
                DiscrepancyCalibrationObservation(
                    "dev", "key", truth=1.0, estimate=1.0, lower=0.0, upper=2.0
                )
            ]
        ),
        raw,
    )
    score = score_locked_discrepancy(
        locked,
        {"key": {**raw.to_dict(), "truth": 1.0}},
        {"key": calibrated},
        target_coverage=0.95,
    )
    assert score["status"] == "ABSTAIN_TRUTH_LEAKAGE"
