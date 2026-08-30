from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    DiscreteModelObservation,
    apply_discrete_model_average,
    fit_discrete_model_weights,
    score_locked_model_average,
)


def _development_records() -> tuple[DiscreteModelObservation, ...]:
    records = []
    for case_id, truths in (("case_a", ("present", "absent")), ("case_b", ("present",))):
        for index, truth in enumerate(truths):
            records.append(
                DiscreteModelObservation(
                    case_id=case_id,
                    target_id=f"edge_{index}",
                    truth=truth,
                    predictions={
                        "good": {
                            "present": 0.90 if truth == "present" else 0.10,
                            "absent": 0.10 if truth == "present" else 0.90,
                        },
                        "bad": {
                            "present": 0.60 if truth == "present" else 0.40,
                            "absent": 0.40 if truth == "present" else 0.60,
                        },
                    },
                )
            )
    return tuple(records)


def test_case_blocked_log_score_prefers_the_better_model_and_freezes_provenance() -> None:
    fit = fit_discrete_model_weights(_development_records())

    assert fit.fit_scope == "development_only"
    assert fit.score_rule == "case_blocked_mean_log_score"
    assert fit.weights["good"] > fit.weights["bad"]
    assert sum(fit.weights.values()) == pytest.approx(1.0)
    assert fit.case_ids == ("case_a", "case_b")
    assert len(fit.fit_data_hash) == 64
    assert fit.observation_count == 3
    assert fit.convergence_status == "CONVERGED"
    assert fit.converged is True
    assert fit.simplex_residual == pytest.approx(0.0)
    assert fit.kkt_residual <= fit.gradient_tolerance
    assert len(fit.objective_trace) >= 1


def test_locked_or_field_truth_cannot_change_the_fit() -> None:
    fit = fit_discrete_model_weights(_development_records())
    locked = (
        DiscreteModelObservation(
            case_id="locked",
            target_id="edge",
            truth="present",
            phase="locked_test",
            predictions={
                "good": {"present": 0.8, "absent": 0.2},
                "bad": {"present": 0.2, "absent": 0.8},
            },
        ),
    )
    first = score_locked_model_average(locked, fit)
    flipped = score_locked_model_average(
        (
            DiscreteModelObservation(
                case_id="locked",
                target_id="edge",
                truth="absent",
                phase="locked_test",
                predictions=locked[0].predictions,
            ),
        ),
        fit,
    )

    assert first["reports"] == flipped["reports"]
    with pytest.raises(ValueError, match="cannot be used for fitting"):
        fit_discrete_model_weights(locked)


def test_incomplete_model_matrix_fails_closed() -> None:
    records = list(_development_records())
    records.append(
        DiscreteModelObservation(
            case_id="case_c",
            target_id="edge",
            truth="present",
            predictions={"good": {"present": 1.0, "absent": 0.0}},
        )
    )
    with pytest.raises(ValueError, match="incomplete"):
        fit_discrete_model_weights(records)


def test_material_model_disagreement_suppresses_reportable_aggregate() -> None:
    fit = fit_discrete_model_weights(_development_records())
    report = apply_discrete_model_average(
        "locked_edge",
        {
            "good": {"present": 0.99, "absent": 0.01},
            "bad": {"present": 0.01, "absent": 0.99},
        },
        fit,
        disagreement_threshold=0.25,
    )

    assert report.reportable is False
    assert report.decision == "ABSTAIN"
    assert report.identifiability == "MODEL_DISAGREEMENT"
    assert report.compatible_outcomes == ()
    assert "not reportable" in report.reason
    assert report.disagreement_status == "MODEL_DISAGREEMENT"
    assert report.weighted_pairwise_total_variation > 0.0


def test_locked_score_reports_proper_score_and_disagreement_rate() -> None:
    fit = fit_discrete_model_weights(_development_records())
    locked = (
        DiscreteModelObservation(
            case_id="locked",
            target_id="agree",
            truth="present",
            phase="locked_test",
            predictions={
                "good": {"present": 0.8, "absent": 0.2},
                "bad": {"present": 0.7, "absent": 0.3},
            },
        ),
    )
    score = score_locked_model_average(locked, fit)

    assert score["status"] == "scored"
    assert score["n"] == 1
    assert score["mean_log_score"] < 0.0
    assert 0.0 <= score["disagreement_rate"] <= 1.0
    assert score["held_out_scoring"] is True
    assert set(score["per_model_mean_log_score"]) == {"good", "bad"}


def test_non_converged_fit_abstains_and_is_not_silently_scored() -> None:
    fit = fit_discrete_model_weights(
        _development_records(),
        max_iterations=1,
        tolerance=1.0e-14,
        gradient_tolerance=1.0e-14,
    )

    assert fit.converged is False
    assert fit.convergence_status == "ABSTAIN"
    report = apply_discrete_model_average(
        "locked_edge",
        {
            "good": {"present": 0.8, "absent": 0.2},
            "bad": {"present": 0.2, "absent": 0.8},
        },
        fit,
    )
    assert report.decision == "ABSTAIN"
    assert report.identifiability == "FIT_NOT_CONVERGED"
    assert report.fit_converged is False


def test_empty_application_or_mismatched_models_is_rejected() -> None:
    fit = fit_discrete_model_weights(_development_records())
    with pytest.raises(ValueError, match="At least one model"):
        apply_discrete_model_average("edge", {}, fit)
    with pytest.raises(ValueError, match="do not match"):
        apply_discrete_model_average(
            "edge",
            {"good": {"present": 1.0, "absent": 0.0}},
            fit,
        )
