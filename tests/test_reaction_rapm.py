"""Tests for the development-fitted RAPM-style reaction evidence model."""

from __future__ import annotations

import math

import pytest

from hydrosheaf.validation.reaction_rapm import (
    ABSTAIN,
    NULL_FAMILY,
    ReactionRAPM,
    ReactionRAPMCalibrationExample,
    ReactionRAPMTrainingExample,
    cross_fitted_reaction_rapm_calibration_examples,
    fit_reaction_rapm,
    fit_reaction_rapm_calibrator,
    score_reaction_rapm_outputs,
    training_examples_from_observations,
)


def _edge_examples() -> tuple[ReactionRAPMTrainingExample, ...]:
    examples: list[ReactionRAPMTrainingExample] = []
    for case_index in range(4):
        case_id = f"case-{case_index}"
        examples.extend(
            (
                ReactionRAPMTrainingExample(
                    case_id,
                    f"carbonate-{case_index}",
                    {"Ca": 1.0, "HCO3": 2.0},
                    {"Ca": 1.20, "HCO3": 2.20},
                    "carbonate_weathering",
                ),
                ReactionRAPMTrainingExample(
                    case_id,
                    f"denit-{case_index}",
                    {"NO3": 1.0, "HCO3": 2.0},
                    {"NO3": 0.80, "HCO3": 2.20},
                    "denitrification",
                ),
                ReactionRAPMTrainingExample(
                    case_id,
                    f"sulfate-{case_index}",
                    {"SO4": 1.0, "HCO3": 2.0},
                    {"SO4": 0.80, "HCO3": 2.20},
                    "sulfate_reduction",
                ),
                ReactionRAPMTrainingExample(
                    case_id,
                    f"mixing-{case_index}",
                    {"Ca": 1.0, "HCO3": 2.0},
                    {"Ca": 1.0, "HCO3": 2.0},
                    "mixing",
                ),
            )
        )
    return tuple(examples)


def test_fit_is_case_blocked_and_predicts_without_truth() -> None:
    model = fit_reaction_rapm(_edge_examples())

    assert model.fit_case_ids == ("case-0", "case-1", "case-2", "case-3")
    assert model.cv_case_count == 4
    assert model.to_audit_record()["fit_scope"] == "development_only"
    assert model.to_audit_record()["control"]["truth_blind_prediction"] is True
    assert "delta_SiO2" in model.feature_names
    assert model.to_audit_record()["identifiability"]["design_rank"] > 0

    observations = {
        "reaction_chemistry": {
            "edges": {
                "candidate": {
                    "upstream": {"Ca": 1.0, "HCO3": 2.0},
                    "downstream": {"Ca": 1.20, "HCO3": 2.20},
                }
            }
        }
    }
    first = model.predict(observations)
    second = model.predict(observations)

    assert first == second
    prediction = first["candidate"]
    assert prediction["candidate_family"] == "carbonate"
    assert prediction["family"] == "carbonate"
    assert prediction["decision"] == "select"
    assert sum(prediction["probabilities"].values()) == pytest.approx(1.0)
    assert prediction["uncertainty"]["calibrated"] is False
    assert prediction["diagnostics"]["on_off_scores"]


def test_null_mixing_control_is_explicit() -> None:
    model = fit_reaction_rapm(_edge_examples())
    observations = {
        "reaction_chemistry": {
            "edges": {
                "unchanged": {
                    "upstream": {"Ca": 1.0, "HCO3": 2.0},
                    "downstream": {"Ca": 1.0, "HCO3": 2.0},
                }
            }
        }
    }

    prediction = model.predict(observations)["unchanged"]

    assert NULL_FAMILY in prediction["probabilities"]
    assert prediction["candidate_family"] == "mixing"
    assert prediction["family"] == "mixing"
    assert prediction["reason"] == "regularized_adjusted_reaction_family_selected"


def test_unsupported_and_conflicting_inputs_abstain() -> None:
    model = fit_reaction_rapm(_edge_examples())
    observations = {
        "reaction_chemistry": {
            "edges": {
                "unsupported": {
                    "upstream": {"DOC": 1.0},
                    "downstream": {"DOC": 1.2},
                },
                "invalid": {
                    "upstream": {"Ca": 1.0, "HCO3": 2.0},
                    "downstream": {"Ca": -1.0, "HCO3": 2.0},
                },
            }
        }
    }

    predictions = model.predict(observations)
    assert predictions["unsupported"]["decision"] == ABSTAIN
    assert predictions["invalid"]["decision"] == ABSTAIN


def test_truth_fields_are_rejected_and_fit_is_development_only() -> None:
    model = fit_reaction_rapm(_edge_examples())
    with pytest.raises(ValueError, match="Truth/reference field"):
        model.predict(
            {
                "reaction_chemistry": {
                    "edges": {},
                    "reference_edges": {"E": "carbonate"},
                }
            }
        )
    with pytest.raises(ValueError, match="development"):
        ReactionRAPM.fit(_edge_examples(), phase="locked_test")


def test_training_adapter_keeps_truth_out_of_observations() -> None:
    observations = {
        "reaction_chemistry": {
            "edges": {
                "E": {
                    "upstream": {"Ca": 1.0, "HCO3": 2.0},
                    "downstream": {"Ca": 1.2, "HCO3": 2.2},
                }
            }
        }
    }
    examples = training_examples_from_observations(
        observations,
        {"E": "carbonate_weathering"},
        case_id="development-case",
    )
    assert len(examples) == 1
    assert examples[0].truth_family == "carbonate_weathering"
    assert "truth" not in observations["reaction_chemistry"]


def test_scoring_reports_accuracy_log_loss_and_abstention() -> None:
    predictions = {
        "E1": {
            "family": "denitrification",
            "decision": "select",
            "probabilities": {"denitrification": 0.8, "none": 0.2},
        },
        "E2": {
            "family": "none",
            "decision": ABSTAIN,
            "probabilities": {"none": 1.0},
        },
    }
    score = score_reaction_rapm_outputs(
        {"E1": "denitrification", "E2": "mixing", "missing": "carbonate"},
        predictions,
    )

    assert score["status"] == "scored"
    assert score["n"] == 3
    assert score["n_missing_outputs"] == 1
    assert score["n_abstain"] == 1
    assert score["accuracy"] == pytest.approx(1.0 / 3.0)
    assert score["selective_accuracy"] == pytest.approx(1.0)
    assert score["coverage"] == pytest.approx(1.0 / 3.0)
    assert score["outputs_complete"] is False
    assert math.isfinite(score["multiclass_log_loss"])
    assert math.isfinite(score["multiclass_brier"])


def test_scoring_can_include_explicit_decoy_edges() -> None:
    predictions = {
        "true": {
            "family": "carbonate",
            "decision": "select",
            "probabilities": {"carbonate": 0.9, "none": 0.1},
        },
        "decoy": {
            "family": "carbonate",
            "decision": "select",
            "probabilities": {"carbonate": 0.8, "none": 0.2},
        },
    }

    score = score_reaction_rapm_outputs(
        {"true": "carbonate_weathering"},
        predictions,
        candidate_edge_ids=("true", "decoy"),
    )

    assert score["n_reference_edges"] == 1
    assert score["n_decoy_edges"] == 1
    assert score["false_commitment_rate"] == pytest.approx(1.0)


def test_score_reports_calibration_partitions_missingness_and_equivalence() -> None:
    predictions = {
        "reference": {
            "family": "none",
            "decision": ABSTAIN,
            "probabilities": {"none": 0.55, "mixing": 0.45},
            "diagnostics": {"missing_fields": ["SO4"], "n_paired_fields": 2},
        },
        "mixing": {
            "family": "mixing",
            "decision": "select",
            "probabilities": {"none": 0.10, "mixing": 0.90},
            "diagnostics": {"missing_fields": [], "n_paired_fields": 5},
        },
        "decoy": {
            "family": "carbonate",
            "decision": "select",
            "probabilities": {"carbonate": 0.80, "none": 0.20},
            "diagnostics": {"unsupported_fields": ["DOC"], "n_paired_fields": 1},
        },
    }
    score = score_reaction_rapm_outputs(
        {"reference": "none", "mixing": "mixing"},
        predictions,
        candidate_edge_ids=("reference", "mixing", "decoy"),
        model_audit=fit_reaction_rapm(_edge_examples()).to_audit_record(),
    )

    assert score["classwise_reliability"]["none"]["bins"]
    assert score["reference_performance"]["n"] == 2
    assert score["decoy_performance"]["false_commitment_count"] == 1
    assert score["missingness"]["n_with_input_missing_or_invalid_fields"] == 2
    assert score["missingness"]["by_paired_field_count"] == {"1": 1, "2": 1, "5": 1}
    assert score["mixing_none_equivalence"]["selective_accuracy"] == pytest.approx(0.5)
    assert "design_rank" in score["identifiability_diagnostics_contract"]["required_from_model_audit"]
    assert score["identifiability"]["status"] == "supplied_from_model_audit"
    assert score["identifiability"]["design_rank"] > 0


def test_case_cross_fitted_calibrator_is_truth_blind_at_apply() -> None:
    examples = _edge_examples()
    calibration_examples = cross_fitted_reaction_rapm_calibration_examples(examples)
    calibrator = fit_reaction_rapm_calibrator(
        calibration_examples,
        classes=ReactionRAPM().config.classes,
        decision_threshold=0.60,
        decision_margin=0.10,
    )
    model = fit_reaction_rapm(examples)
    raw = model.predict(
        {
            "reaction_chemistry": {
                "edges": {
                    "candidate": {
                        "upstream": {"Ca": 1.0, "HCO3": 2.0},
                        "downstream": {"Ca": 1.20, "HCO3": 2.20},
                    }
                }
            }
        }
    )
    calibrated = calibrator.apply(raw)
    prediction = calibrated["candidate"]

    assert prediction["uncertainty"]["calibrated"] is True
    assert prediction["uncertainty"]["selection_probability_calibrated"] is True
    assert prediction["calibration_status"] == "development_fitted"
    assert prediction["diagnostics"]["calibration_threshold_fixed"] == pytest.approx(0.60)
    audit = calibrator.to_audit_record()
    assert audit["fit_scope"] == "development_only"
    assert audit["control"]["truth_blind_apply"] is True
    assert audit["control"]["threshold_tuned_on_locked_truth"] is False


def test_selection_rule_is_tuned_on_cross_fitted_development_rows() -> None:
    calibration_examples = cross_fitted_reaction_rapm_calibration_examples(_edge_examples())
    calibrator = fit_reaction_rapm_calibrator(
        calibration_examples,
        classes=ReactionRAPM().config.classes,
        decision_threshold=0.60,
        decision_margin=0.10,
        tune_selection_rule=True,
        target_coverage=0.25,
        max_selective_risk=0.40,
    )

    tuning = calibrator.provenance["threshold_tuning"]
    assert calibrator.calibration_kind.endswith("selective_rule_calibration")
    assert tuning["enabled"] is True
    assert tuning["target_met"] is True
    assert tuning["selected_development_coverage"] >= 0.25
    assert tuning["selected_development_selective_risk"] <= 0.40
    assert calibrator.decision_threshold == pytest.approx(tuning["selected_threshold"])
    assert calibrator.decision_margin == pytest.approx(tuning["selected_margin"])
    assert calibrator.to_audit_record()["control"]["threshold_tuned_on_locked_truth"] is False


def test_calibrator_rejects_in_sample_or_locked_fit() -> None:
    with pytest.raises(ValueError, match="case-held-out"):
        ReactionRAPMCalibrationExample(
            case_id="case-a",
            edge_id="edge-a",
            truth_family="carbonate",
            logits={"carbonate": 1.0, "none": 0.0},
            cross_fitted=False,
        )
    rows = tuple(
        ReactionRAPMCalibrationExample(
            case_id=case_id,
            edge_id=f"edge-{case_id}",
            truth_family="carbonate",
            logits={"carbonate": 1.0, "none": 0.0},
        )
        for case_id in ("case-a", "case-b")
    )
    with pytest.raises(ValueError, match="development"):
        fit_reaction_rapm_calibrator(rows, phase="locked_test")
