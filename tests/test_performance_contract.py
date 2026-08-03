from __future__ import annotations

from hydrosheaf.validation.performance_contract import (
    ABSTAIN,
    FAIL,
    PASS,
    aggregate_specialist_gates,
    assess_age_gate,
    assess_integrated_gate,
    assess_kinetics_gate,
    assess_reaction_gate,
)


def _age(**overrides: object) -> dict[str, object]:
    result: dict[str, object] = {
        "coverage_including_abstention": 0.96,
        "acceptance_rate": 0.80,
        "selective_risk": 0.20,
        "relative_width": 1.20,
        "specialist_mae": 4.0,
        "baseline_mae": 4.5,
        "held_out_generators": True,
        "competence_matched_baseline": True,
        "calibrated": True,
        "family_stratified": True,
        "baseline_noninferior": True,
    }
    result.update(overrides)
    return result


def _reaction(**overrides: object) -> dict[str, object]:
    result: dict[str, object] = {
        "coverage": 0.35,
        "selective_risk": 0.25,
        "false_commitment_rate": 0.05,
        "multiclass_log_loss": 1.5,
        "baseline_log_loss": 1.6,
        "max_classwise_ece": 0.20,
        "outputs_complete": True,
        "selection_rule_target_met": True,
        "held_out_generators": True,
        "competence_matched_baseline": True,
        "calibrated": True,
        "classwise_calibrated": True,
        "mechanism_stratified": True,
    }
    result.update(overrides)
    return result


def _kinetics(**overrides: object) -> dict[str, object]:
    result: dict[str, object] = {
        "interval_coverage": 0.93,
        "identifiability_rate": 0.85,
        "selective_risk": 0.20,
        "predictive_rmse": 0.1,
        "parameter_error": 0.2,
        "held_out_kinetic_regimes": True,
        "competence_matched_baseline": True,
        "calibrated": True,
        "ka_confounded_reported": True,
        "transport_stratified": True,
    }
    result.update(overrides)
    return result


def _integrated(**overrides: object) -> dict[str, object]:
    result: dict[str, object] = {
        "model_averaging_converged": True,
        "discrepancy_calibrated": True,
        "prospective_outcomes_complete": True,
        "prospective_evidence_sufficient": True,
        "prospective_case_count": 6,
        "prospective_outcomes_independent": True,
        "prospective_random_baseline_valid": True,
        "paired_uncertainty_available": True,
        "paired_random_delta_ci_low": 0.02,
        "paired_random_delta_ci_high": 0.40,
        "paired_specialist_delta_ci_low": -0.01,
        "source_snapshot_recorded": True,
        "improves_over_random": True,
        "noninferior_to_strongest_specialist": True,
        "no_material_calibration_degradation": True,
        "false_commitment_controlled": True,
        "raw_benchmark_verified": True,
        "hydrosheaf_mean_utility_per_cost": 0.30,
        "random_mean_utility_per_cost": 0.10,
        "strongest_specialist_mean_utility_per_cost": 0.28,
        "calibration_degradation": 0.01,
        "observed_false_commitment_rate": 0.04,
    }
    result.update(overrides)
    return result


def test_complete_specialist_contract_passes_only_with_all_evidence_layers():
    age = assess_age_gate(_age())
    reaction = assess_reaction_gate(_reaction())
    kinetics = assess_kinetics_gate(_kinetics())
    integrated = assess_integrated_gate(_integrated())

    assert age.status == PASS
    assert reaction.status == PASS
    assert kinetics.status == PASS
    assert integrated.status == PASS
    summary = aggregate_specialist_gates(
        age=age, reaction=reaction, kinetics=kinetics, integrated=integrated
    )
    assert summary["claim_status"] == PASS


def test_missing_competence_or_calibration_abstains_instead_of_passing():
    result = assess_age_gate(
        _age(competence_matched_baseline=None, calibrated=None)
    )

    assert result.status == ABSTAIN
    assert "missing_or_unverified:competence_matched_baseline" in result.reasons
    assert "missing_or_unverified:calibrated" in result.reasons


def test_reaction_gate_does_not_reward_low_coverage_abstention():
    result = assess_reaction_gate(_reaction(coverage=0.10, selective_risk=0.0))

    assert result.status == FAIL
    assert "failed:minimum_coverage" in result.reasons


def test_kinetics_gate_requires_explicit_k_times_area_boundary():
    result = assess_kinetics_gate(_kinetics(ka_confounded_reported=False))

    assert result.status == FAIL
    assert "failed:ka_confounded_reported" in result.reasons


def test_kinetics_gate_requires_predictive_and_parameter_error_limits():
    assert assess_kinetics_gate(_kinetics(predictive_rmse=None)).status == ABSTAIN
    assert assess_kinetics_gate(_kinetics(parameter_error=0.80)).status == FAIL


def test_integrated_gate_abstains_when_prospective_outcomes_are_missing():
    result = assess_integrated_gate(_integrated(prospective_outcomes_complete=None))

    assert result.status == ABSTAIN


def test_integrated_gate_recomputes_utility_comparisons_from_numeric_records():
    result = assess_integrated_gate(
        _integrated(
            improves_over_random=False,
            hydrosheaf_mean_utility_per_cost=0.05,
            random_mean_utility_per_cost=0.10,
        )
    )

    assert result.status == FAIL
    assert "failed:improves_over_random" in result.reasons
