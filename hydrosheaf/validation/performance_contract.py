"""Shared performance gates for the HydroSheaf specialist programme.

The synthetic programme produces many diagnostic records.  This module keeps
the scientific claim decision separate from those records and makes every
specialist gate fail closed when a required metric, comparator, or validation
layer is absent.  A ``PASS`` here is a bounded controlled-synthetic claim only;
it is never a field-validation or universal-superiority claim.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping


PASS = "PASS"
FAIL = "FAIL"
ABSTAIN = "ABSTAIN"


def _finite(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def _bool(value: object) -> bool | None:
    return value if isinstance(value, bool) else None


@dataclass(frozen=True)
class GateResult:
    """Auditable result of one predeclared specialist claim gate."""

    name: str
    status: str
    requirements: Mapping[str, object]
    observed: Mapping[str, object]
    reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.status not in {PASS, FAIL, ABSTAIN}:
            raise ValueError(f"Unsupported gate status: {self.status!r}")

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "status": self.status,
            "requirements": dict(self.requirements),
            "observed": dict(self.observed),
            "reasons": list(self.reasons),
        }


def _evaluate(
    name: str,
    conditions: Mapping[str, bool | None],
    *,
    requirements: Mapping[str, object],
    observed: Mapping[str, object],
) -> GateResult:
    missing = tuple(key for key, value in conditions.items() if value is None)
    failed = tuple(key for key, value in conditions.items() if value is False)
    if missing:
        status = ABSTAIN
        reasons = tuple(f"missing_or_unverified:{key}" for key in missing)
    elif failed:
        status = FAIL
        reasons = tuple(f"failed:{key}" for key in failed)
    else:
        status = PASS
        reasons = ()
    return GateResult(
        name=name,
        status=status,
        requirements=requirements,
        observed=observed,
        reasons=reasons,
    )


def assess_age_gate(
    metrics: Mapping[str, object],
    *,
    target_coverage: float = 0.95,
    minimum_acceptance_rate: float = 0.50,
    max_selective_risk: float = 0.30,
    max_selective_risk_years: float | None = None,
    max_relative_width: float = 1.50,
) -> GateResult:
    """Assess a calibrated, competence-matched age specialist result.

    Age selective risk is an age error in years, not a probability.  The
    legacy ``max_selective_risk`` argument remains accepted for compatibility;
    callers should use ``max_selective_risk_years`` for an unambiguous gate.
    """

    selective_risk_limit = (
        max_selective_risk
        if max_selective_risk_years is None
        else float(max_selective_risk_years)
    )

    coverage = _finite(
        metrics.get("coverage_including_abstention", metrics.get("coverage"))
    )
    acceptance_rate = _finite(
        metrics.get("acceptance_rate", metrics.get("coverage"))
    )
    selective_risk = _finite(metrics.get("selective_risk"))
    relative_width = _finite(metrics.get("relative_width"))
    specialist_mae = _finite(metrics.get("specialist_mae"))
    baseline_mae = _finite(metrics.get("baseline_mae"))
    observed = {
        "coverage_including_abstention": coverage,
        "acceptance_rate": acceptance_rate,
        "selective_risk": selective_risk,
        "mean_width": _finite(metrics.get("mean_width")),
        "relative_width": relative_width,
        "mae": _finite(metrics.get("mae", metrics.get("mean_absolute_error"))),
        "specialist_mae": specialist_mae,
        "baseline_mae": baseline_mae,
        "held_out_generators": _bool(metrics.get("held_out_generators")),
        "competence_matched_baseline": _bool(
            metrics.get("competence_matched_baseline")
        ),
        "calibrated": _bool(metrics.get("calibrated")),
        "family_stratified": _bool(metrics.get("family_stratified")),
        "baseline_noninferior": _bool(metrics.get("baseline_noninferior")),
    }
    conditions: dict[str, bool | None] = {
        "coverage": None if coverage is None else coverage >= target_coverage,
        "minimum_acceptance_rate": (
            None
            if acceptance_rate is None
            else acceptance_rate >= minimum_acceptance_rate
        ),
        "selective_risk": (
            None
            if selective_risk is None
            else selective_risk <= selective_risk_limit
        ),
        "relative_width": (
            None if relative_width is None else relative_width <= max_relative_width
        ),
        "held_out_generators": observed["held_out_generators"],
        "competence_matched_baseline": observed["competence_matched_baseline"],
        "calibrated": observed["calibrated"],
        "family_stratified": observed["family_stratified"],
        "baseline_noninferior": (
            None
            if specialist_mae is None or baseline_mae is None
            else specialist_mae <= baseline_mae
        ),
    }
    return _evaluate(
        "age_specialist_performance",
        conditions,
        requirements={
            "target_coverage": target_coverage,
            "minimum_acceptance_rate": minimum_acceptance_rate,
            "max_selective_risk_years": selective_risk_limit,
            "max_relative_width": max_relative_width,
            "claim_scope": "bounded_controlled_synthetic_age_inference",
        },
        observed=observed,
    )


def assess_reaction_gate(
    metrics: Mapping[str, object],
    *,
    minimum_coverage: float = 0.25,
    max_selective_risk: float = 0.40,
    max_false_commitment: float = 0.10,
) -> GateResult:
    """Assess a calibrated reaction-family specialist without rewarding abstention."""

    coverage = _finite(metrics.get("coverage"))
    selective_risk = _finite(metrics.get("selective_risk"))
    false_commitment = _finite(metrics.get("false_commitment_rate"))
    observed = {
        "coverage": coverage,
        "selective_risk": selective_risk,
        "false_commitment_rate": false_commitment,
        "selective_accuracy": _finite(metrics.get("selective_accuracy")),
        "multiclass_log_loss": _finite(metrics.get("multiclass_log_loss")),
        "baseline_log_loss": _finite(metrics.get("baseline_log_loss")),
        "max_classwise_ece": _finite(metrics.get("max_classwise_ece")),
        "outputs_complete": _bool(metrics.get("outputs_complete")),
        "selection_rule_target_met": _bool(
            metrics.get("selection_rule_target_met")
        ),
        "held_out_generators": _bool(metrics.get("held_out_generators")),
        "competence_matched_baseline": _bool(
            metrics.get("competence_matched_baseline")
        ),
        "calibrated": _bool(metrics.get("calibrated")),
        "classwise_calibrated": _bool(metrics.get("classwise_calibrated")),
        "mechanism_stratified": _bool(metrics.get("mechanism_stratified")),
    }
    baseline_log_loss = observed["baseline_log_loss"]
    model_log_loss = observed["multiclass_log_loss"]
    log_loss_improvement = (
        None
        if baseline_log_loss is None or model_log_loss is None
        else model_log_loss <= baseline_log_loss
    )
    observed["log_loss_noninferior"] = log_loss_improvement
    conditions: dict[str, bool | None] = {
        "minimum_coverage": None if coverage is None else coverage >= minimum_coverage,
        "selective_risk": (
            None if selective_risk is None else selective_risk <= max_selective_risk
        ),
        "false_commitment": (
            None
            if false_commitment is None
            else false_commitment <= max_false_commitment
        ),
        "log_loss_noninferior": log_loss_improvement,
        "outputs_complete": observed["outputs_complete"],
        "selection_rule_target_met": observed["selection_rule_target_met"],
        "held_out_generators": observed["held_out_generators"],
        "competence_matched_baseline": observed["competence_matched_baseline"],
        "calibrated": observed["calibrated"],
        "classwise_calibrated": observed["classwise_calibrated"],
        "classwise_ece": (
            None
            if observed["max_classwise_ece"] is None
            else observed["max_classwise_ece"] <= 0.35
        ),
        "mechanism_stratified": observed["mechanism_stratified"],
    }
    return _evaluate(
        "reaction_specialist_performance",
        conditions,
        requirements={
            "minimum_coverage": minimum_coverage,
            "max_selective_risk": max_selective_risk,
            "max_false_commitment": max_false_commitment,
            "max_classwise_ece": 0.35,
            "claim_scope": "bounded_controlled_synthetic_reaction_family_inference",
        },
        observed=observed,
    )


def assess_kinetics_gate(
    metrics: Mapping[str, object],
    *,
    minimum_interval_coverage: float = 0.90,
    minimum_identifiability_rate: float = 0.80,
    max_selective_risk: float = 0.40,
    max_predictive_rmse: float = 0.25,
    max_parameter_error: float = 0.50,
) -> GateResult:
    """Assess held-out M8 prediction and conditional k/A recovery.

    identifiability_rate is explicitly conditional on an independent
    surface-area measurement. The unconditional confounded fraction remains
    an observed diagnostic and must be reported alongside this gate.
    """

    interval_coverage = _finite(metrics.get("interval_coverage"))
    identifiability = _finite(metrics.get("identifiability_rate"))
    selective_risk = _finite(metrics.get("selective_risk"))
    observed = {
        "interval_coverage": interval_coverage,
        "identifiability_rate": identifiability,
        "identifiability_rate_overall": _finite(
            metrics.get("identifiability_rate_overall")
        ),
        "parameter_abstention_rate": _finite(
            metrics.get("parameter_abstention_rate")
        ),
        "selective_risk": selective_risk,
        "predictive_rmse": _finite(metrics.get("predictive_rmse")),
        "parameter_error": _finite(metrics.get("parameter_error")),
        "held_out_kinetic_regimes": _bool(metrics.get("held_out_kinetic_regimes")),
        "competence_matched_baseline": _bool(
            metrics.get("competence_matched_baseline")
        ),
        "calibrated": _bool(metrics.get("calibrated")),
        "ka_confounded_reported": _bool(metrics.get("ka_confounded_reported")),
        "transport_stratified": _bool(metrics.get("transport_stratified")),
    }
    conditions: dict[str, bool | None] = {
        "interval_coverage": (
            None
            if interval_coverage is None
            else interval_coverage >= minimum_interval_coverage
        ),
        "identifiability_rate": (
            None
            if identifiability is None
            else identifiability >= minimum_identifiability_rate
        ),
        "selective_risk": (
            None if selective_risk is None else selective_risk <= max_selective_risk
        ),
        "predictive_rmse": (
            None
            if observed["predictive_rmse"] is None
            else float(observed["predictive_rmse"]) <= max_predictive_rmse
        ),
        "parameter_error": (
            None
            if observed["parameter_error"] is None
            else float(observed["parameter_error"]) <= max_parameter_error
        ),
        "held_out_kinetic_regimes": observed["held_out_kinetic_regimes"],
        "competence_matched_baseline": observed["competence_matched_baseline"],
        "calibrated": observed["calibrated"],
        "ka_confounded_reported": observed["ka_confounded_reported"],
        "transport_stratified": observed["transport_stratified"],
    }
    return _evaluate(
        "m8_kinetics_performance",
        conditions,
        requirements={
            "minimum_interval_coverage": minimum_interval_coverage,
            "minimum_identifiability_rate": minimum_identifiability_rate,
            "max_selective_risk": max_selective_risk,
            "max_predictive_rmse": max_predictive_rmse,
            "max_parameter_error": max_parameter_error,
            "identifiability_rate_definition": (
                "identified k/A cases divided by locked cases with an independent "
                "surface-area measurement; all no-area cases remain structurally "
                "confounded and are reported as parameter abstentions"
            ),
            "claim_scope": "bounded_controlled_synthetic_kinetic_inference",
        },
        observed=observed,
    )


def assess_integrated_gate(
    metrics: Mapping[str, object],
    *,
    minimum_utility_improvement: float = 0.0,
    noninferiority_margin: float = 0.05,
    max_calibration_degradation: float = 0.05,
    max_false_commitment_rate: float = 0.10,
    minimum_prospective_case_count: int = 6,
) -> GateResult:
    """Assess the joint model/discrepancy/decision claim gate."""

    if int(minimum_prospective_case_count) < 1:
        raise ValueError("minimum_prospective_case_count must be positive")

    observed = {
        "model_averaging_converged": _bool(metrics.get("model_averaging_converged")),
        "discrepancy_calibrated": _bool(metrics.get("discrepancy_calibrated")),
        "prospective_outcomes_complete": _bool(
            metrics.get("prospective_outcomes_complete")
        ),
        "prospective_evidence_sufficient": _bool(
            metrics.get("prospective_evidence_sufficient")
        ),
        "prospective_case_count": _finite(metrics.get("prospective_case_count")),
        "prospective_outcomes_independent": _bool(
            metrics.get("prospective_outcomes_independent")
        ),
        "prospective_random_baseline_valid": _bool(
            metrics.get("prospective_random_baseline_valid")
        ),
        "paired_uncertainty_available": _bool(
            metrics.get("paired_uncertainty_available")
        ),
        "paired_random_delta_ci_low": _finite(
            metrics.get("paired_random_delta_ci_low")
        ),
        "paired_specialist_delta_ci_low": _finite(
            metrics.get("paired_specialist_delta_ci_low")
        ),
        "source_snapshot_recorded": _bool(metrics.get("source_snapshot_recorded")),
        "improves_over_random": _bool(metrics.get("improves_over_random")),
        "noninferior_to_strongest_specialist": _bool(
            metrics.get("noninferior_to_strongest_specialist")
        ),
        "no_material_calibration_degradation": _bool(
            metrics.get("no_material_calibration_degradation")
        ),
        "false_commitment_controlled": _bool(
            metrics.get("false_commitment_controlled")
        ),
        "raw_benchmark_verified": _bool(metrics.get("raw_benchmark_verified")),
        "hydrosheaf_mean_utility_per_cost": _finite(
            metrics.get("hydrosheaf_mean_utility_per_cost")
        ),
        "random_mean_utility_per_cost": _finite(
            metrics.get("random_mean_utility_per_cost")
        ),
        "strongest_specialist_mean_utility_per_cost": _finite(
            metrics.get("strongest_specialist_mean_utility_per_cost")
        ),
        "calibration_degradation": _finite(metrics.get("calibration_degradation")),
        "observed_false_commitment_rate": _finite(
            metrics.get("observed_false_commitment_rate")
        ),
    }
    hydro_utility = observed["hydrosheaf_mean_utility_per_cost"]
    random_utility = observed["random_mean_utility_per_cost"]
    specialist_utility = observed["strongest_specialist_mean_utility_per_cost"]
    calibration_degradation = observed["calibration_degradation"]
    false_commitment = observed["observed_false_commitment_rate"]
    improves_over_random = (
        None
        if hydro_utility is None or random_utility is None
        else float(hydro_utility) - float(random_utility)
        >= minimum_utility_improvement
    )
    noninferior_to_specialist = (
        None
        if hydro_utility is None or specialist_utility is None
        else float(hydro_utility) >= float(specialist_utility) - noninferiority_margin
    )
    no_calibration_degradation = (
        None
        if calibration_degradation is None
        else float(calibration_degradation) <= max_calibration_degradation
    )
    false_commitment_controlled = (
        None
        if false_commitment is None
        else float(false_commitment) <= max_false_commitment_rate
    )
    paired_random_ci_low = observed["paired_random_delta_ci_low"]
    paired_specialist_ci_low = observed["paired_specialist_delta_ci_low"]
    improves_over_random_with_uncertainty = (
        None
        if paired_random_ci_low is None
        else float(paired_random_ci_low) > minimum_utility_improvement
    )
    noninferior_to_specialist_with_uncertainty = (
        None
        if paired_specialist_ci_low is None
        else float(paired_specialist_ci_low) >= -noninferiority_margin
    )
    observed["recomputed_improves_over_random"] = improves_over_random
    observed["recomputed_noninferior_to_specialist"] = noninferior_to_specialist
    observed["recomputed_improves_over_random_with_uncertainty"] = (
        improves_over_random_with_uncertainty
    )
    observed["recomputed_noninferior_to_specialist_with_uncertainty"] = (
        noninferior_to_specialist_with_uncertainty
    )
    observed["recomputed_no_material_calibration_degradation"] = no_calibration_degradation
    observed["recomputed_false_commitment_controlled"] = false_commitment_controlled
    conditions = {
        "model_averaging_converged": observed["model_averaging_converged"],
        "discrepancy_calibrated": observed["discrepancy_calibrated"],
        "prospective_outcomes_complete": observed["prospective_outcomes_complete"],
        "prospective_evidence_sufficient": observed[
            "prospective_evidence_sufficient"
        ]
        and observed["prospective_case_count"] is not None
        and observed["prospective_case_count"] >= int(minimum_prospective_case_count),
        "prospective_outcomes_independent": observed[
            "prospective_outcomes_independent"
        ],
        "prospective_random_baseline_valid": observed[
            "prospective_random_baseline_valid"
        ],
        "paired_uncertainty_available": observed["paired_uncertainty_available"],
        "source_snapshot_recorded": observed["source_snapshot_recorded"],
        "raw_benchmark_verified": observed["raw_benchmark_verified"],
        "improves_over_random": improves_over_random,
        "improves_over_random_with_uncertainty": improves_over_random_with_uncertainty,
        "noninferior_to_strongest_specialist": noninferior_to_specialist,
        "noninferior_to_strongest_specialist_with_uncertainty": (
            noninferior_to_specialist_with_uncertainty
        ),
        "no_material_calibration_degradation": no_calibration_degradation,
        "false_commitment_controlled": false_commitment_controlled,
    }
    return _evaluate(
        "integrated_hydrosheaf_performance",
        conditions,
        requirements={
            "claim_scope": "bounded_controlled_synthetic_integrated_decision_performance",
            "field_validation": "not_established_by_this_gate",
            "minimum_utility_improvement": minimum_utility_improvement,
            "noninferiority_margin": noninferiority_margin,
            "max_calibration_degradation": max_calibration_degradation,
            "max_false_commitment_rate": max_false_commitment_rate,
            "prospective_minimum_locked_cases": int(minimum_prospective_case_count),
            "paired_ci_level": 0.95,
            "practical_improvement_requires_ci_lower_above": minimum_utility_improvement,
        },
        observed=observed,
    )


def aggregate_specialist_gates(
    *,
    age: GateResult,
    reaction: GateResult,
    kinetics: GateResult,
    integrated: GateResult,
) -> dict[str, object]:
    """Return one claim-safe summary for the locked performance package."""

    gates = {
        "age": age.to_dict(),
        "reaction": reaction.to_dict(),
        "kinetics": kinetics.to_dict(),
        "integrated": integrated.to_dict(),
    }
    statuses = [age.status, reaction.status, kinetics.status, integrated.status]
    if all(status == PASS for status in statuses):
        claim_status = PASS
    elif any(status == ABSTAIN for status in statuses):
        claim_status = ABSTAIN
    else:
        claim_status = FAIL
    return {
        "schema_version": "1.0",
        "claim_status": claim_status,
        "claim_scope": "controlled_synthetic_only; no field or universal-superiority claim",
        "gates": gates,
    }


__all__ = [
    "ABSTAIN",
    "FAIL",
    "PASS",
    "GateResult",
    "aggregate_specialist_gates",
    "assess_age_gate",
    "assess_integrated_gate",
    "assess_kinetics_gate",
    "assess_reaction_gate",
]
