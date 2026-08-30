"""Focused tests for the standalone synthetic claim and outcome scaffolding."""

from __future__ import annotations

import math

import pytest

from hydrosheaf.validation.synthetic_claims import (
    ClaimStatus,
    ComparatorAudit,
    GateThresholds,
    HeldOutEvidence,
    MetricThreshold,
    PreMeasurementPosteriorSummary,
    SyntheticObservationModel,
    TargetMetric,
    evaluate_synthetic_claims,
    simulate_post_measurement,
)


def _comparator(
    *,
    complete: bool = True,
    truth_blind: bool = True,
    independent_candidate_universe: bool = True,
    candidate_universe_hash: str | None = "sha256:independent-universe",
    metric: float = 0.78,
) -> ComparatorAudit:
    return ComparatorAudit(
        complete=complete,
        truth_blind=truth_blind,
        independent_candidate_universe=independent_candidate_universe,
        candidate_universe_scope="all_eligible_targets_v1",
        candidate_universe_hash=candidate_universe_hash,
        input_channels=("tracer_age", "reaction_chemistry"),
        metrics={"f1": metric},
    )


def _evidence(
    *,
    primary_metric: float = 0.82,
    comparator_metric: float = 0.78,
    brier_score: float | None = 0.12,
    log_loss: float | None = 0.30,
    generator_case_counts: dict[str, int] | None = None,
    calibration_coverage: float | None = 0.96,
    calibration_error: float | None = 0.04,
    mean_interval_width: float | None = 20.0,
    abstention_rate: float | None = 0.20,
    acceptance_rate: float | None = 0.80,
    selective_risk: float | None = 0.12,
    comparator: ComparatorAudit | None = None,
) -> HeldOutEvidence:
    return HeldOutEvidence(
        evaluation_split="locked_test",
        held_out_cases=4,
        generator_case_counts=generator_case_counts
        or {"analytic_lattice": 2, "independent_mixing": 2},
        held_out_metrics={
            "f1": primary_metric,
            "brier": brier_score,
            "log_loss": log_loss,
        },
        calibration_coverage=calibration_coverage,
        calibration_error=calibration_error,
        mean_interval_width=mean_interval_width,
        abstention_rate=abstention_rate,
        acceptance_rate=acceptance_rate,
        selective_risk=selective_risk,
        comparators={
            "local_specialist": comparator
            or _comparator(metric=comparator_metric)
        },
    )


def _thresholds(claim_name: str) -> GateThresholds:
    return GateThresholds(
        claim_name=claim_name,
        evaluation_split="locked_test",
        min_held_out_cases=4,
        min_cases_per_generator=2,
        required_generator_families=("analytic_lattice", "independent_mixing"),
        required_comparators=("local_specialist",),
        primary_metric="f1",
        held_out_metric_thresholds={
            "f1": MetricThreshold(value=0.75, direction="at_least"),
            "brier": MetricThreshold(value=0.20, direction="at_most"),
            "log_loss": MetricThreshold(value=0.45, direction="at_most"),
        },
        min_calibration_coverage=0.95,
        max_calibration_error=0.05,
        max_mean_interval_width=25.0,
        max_abstention_rate=0.25,
        min_acceptance_rate=0.75,
        max_selective_risk=0.15,
        comparator_metric="f1",
        max_comparator_degradation=0.05,
        require_truth_blind_comparators=True,
        require_independent_candidate_universes=True,
    )


def test_component_integrated_and_field_statuses_are_separate() -> None:
    evidence = _evidence()
    report = evaluate_synthetic_claims(
        evidence,
        evidence,
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    assert report.synthetic_component_claim.status is ClaimStatus.PASS
    assert report.synthetic_integrated_claim.status is ClaimStatus.PASS
    assert report.field_validation_claim.status is ClaimStatus.DEFERRED
    assert report.synthetic_status is ClaimStatus.PASS
    assert report.to_dict()["field_validation_claim"]["status"] == "DEFERRED"


def test_field_deferred_does_not_block_synthetic_status() -> None:
    evidence = _evidence()
    report = evaluate_synthetic_claims(
        evidence,
        evidence,
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    # There is intentionally no field evidence argument.  The two synthetic
    # gates can pass while the field claim remains hard-deferred.
    assert report.synthetic_status is ClaimStatus.PASS
    assert report.field_validation_claim.status is ClaimStatus.DEFERRED


def test_missing_calibration_or_selective_risk_abstains_instead_of_passing() -> None:
    report = evaluate_synthetic_claims(
        _evidence(calibration_coverage=None, selective_risk=None, brier_score=None),
        None,
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    assert report.synthetic_component_claim.status is ClaimStatus.ABSTAIN
    assert report.synthetic_integrated_claim.status is ClaimStatus.ABSTAIN
    assert any("missing" in reason for reason in report.synthetic_component_claim.reasons)
    assert report.field_validation_claim.status is ClaimStatus.DEFERRED


def test_explicit_primary_metric_failure_is_fail_not_abstain() -> None:
    report = evaluate_synthetic_claims(
        _evidence(primary_metric=0.70),
        _evidence(primary_metric=0.82),
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    assert report.synthetic_component_claim.status is ClaimStatus.FAIL
    assert report.synthetic_integrated_claim.status is ClaimStatus.PASS
    assert any("misses the declared threshold" in reason for reason in report.synthetic_component_claim.reasons)


@pytest.mark.parametrize(
    "audit",
    [
        _comparator(truth_blind=False),
        _comparator(independent_candidate_universe=False),
        _comparator(candidate_universe_hash=None),
        _comparator(complete=False),
    ],
)
def test_comparator_leakage_or_incompleteness_cannot_pass(audit: ComparatorAudit) -> None:
    report = evaluate_synthetic_claims(
        _evidence(comparator=audit),
        _evidence(comparator=audit),
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    assert report.synthetic_component_claim.status is ClaimStatus.ABSTAIN
    assert report.synthetic_integrated_claim.status is ClaimStatus.ABSTAIN
    assert any("not adjudicable" in reason for reason in report.synthetic_component_claim.reasons)


def test_all_abstention_is_not_low_selective_risk_success() -> None:
    report = evaluate_synthetic_claims(
        _evidence(abstention_rate=1.0, acceptance_rate=0.0, selective_risk=0.0),
        _evidence(abstention_rate=1.0, acceptance_rate=0.0, selective_risk=0.0),
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    assert report.synthetic_component_claim.status is ClaimStatus.FAIL
    assert report.synthetic_integrated_claim.status is ClaimStatus.FAIL
    assert any("below the declared minimum" in reason for reason in report.synthetic_component_claim.reasons)


def test_generator_coverage_is_explicit_and_denominator_checked() -> None:
    report = evaluate_synthetic_claims(
        _evidence(generator_case_counts={"analytic_lattice": 4, "independent_mixing": 0}),
        _evidence(generator_case_counts={"analytic_lattice": 4, "independent_mixing": 0}),
        _thresholds("synthetic_component_claim"),
        _thresholds("synthetic_integrated_claim"),
    )

    assert report.synthetic_component_claim.status is ClaimStatus.FAIL
    assert any("under-covered" in reason for reason in report.synthetic_component_claim.reasons)


def test_post_measurement_expected_improvement_uses_truth_conditioned_outcomes() -> None:
    posterior = PreMeasurementPosteriorSummary(
        posterior={"young": 0.5, "old": 0.5},
        decision="MEASURE",
        action="tracer_panel",
    )
    model = SyntheticObservationModel(
        likelihoods_by_action={
            "tracer_panel": {
                "young": {"low": 0.9, "high": 0.1},
                "old": {"low": 0.1, "high": 0.9},
            }
        },
        target_metric=TargetMetric.BAYES_ERROR.value,
    )

    result = simulate_post_measurement(
        posterior,
        model,
        "young",
        minimum_improvement=0.20,
    )

    assert result.evaluated is True
    assert result.improved is True
    assert result.pre_metric == pytest.approx(0.5)
    assert result.expected_post_metric == pytest.approx(0.1)
    assert result.expected_improvement == pytest.approx(0.4)
    assert [outcome.outcome for outcome in result.outcomes] == ["high", "low"]
    assert sum(outcome.probability_given_truth for outcome in result.outcomes) == pytest.approx(1.0)
    assert result.outcomes[1].posterior["young"] == pytest.approx(0.9)


def test_post_measurement_abstention_is_unevaluated_not_success() -> None:
    posterior = PreMeasurementPosteriorSummary(
        posterior={"young": 0.5, "old": 0.5},
        decision="ABSTAIN",
        action=None,
    )
    model = SyntheticObservationModel(
        likelihoods_by_action={
            "tracer_panel": {
                "young": {"low": 0.9, "high": 0.1},
                "old": {"low": 0.1, "high": 0.9},
            }
        },
        target_metric="bayes_error",
    )

    result = simulate_post_measurement(
        posterior,
        model,
        "young",
        minimum_improvement=0.0,
    )

    assert result.evaluated is False
    assert result.improved is False
    assert result.expected_post_metric is None
    assert result.expected_improvement is None
    assert result.outcomes == ()
    assert result.reason == "policy_abstained_no_post_measurement_score"


def test_noninformative_measurement_does_not_pass_strict_improvement_gate() -> None:
    posterior = PreMeasurementPosteriorSummary(
        posterior={"a": 0.5, "b": 0.5},
        decision="MEASURE",
        action="uninformative",
    )
    model = SyntheticObservationModel(
        likelihoods_by_action={
            "uninformative": {
                "a": {"same": 1.0},
                "b": {"same": 1.0},
            }
        },
        target_metric="bayes_error",
    )

    result = simulate_post_measurement(
        posterior,
        model,
        "a",
        minimum_improvement=0.0,
    )

    assert result.evaluated is True
    assert result.improved is False
    assert result.expected_improvement == pytest.approx(0.0)
    assert result.reason == "target_metric_did_not_meet_improvement_threshold"


def test_simulator_rejects_truth_or_model_states_that_cannot_be_scored() -> None:
    posterior = PreMeasurementPosteriorSummary(
        posterior={"a": 0.5, "b": 0.5},
        decision="MEASURE",
        action="measure",
    )
    model = SyntheticObservationModel(
        likelihoods_by_action={"measure": {"a": {"yes": 1.0}, "b": {"yes": 1.0}}},
        target_metric="posterior_entropy",
    )

    with pytest.raises(ValueError, match="true_hypothesis"):
        simulate_post_measurement(posterior, model, "missing", minimum_improvement=0.0)

    incomplete_model = SyntheticObservationModel(
        likelihoods_by_action={"measure": {"a": {"yes": 1.0}}},
        target_metric="posterior_entropy",
    )
    with pytest.raises(ValueError, match="missing posterior states"):
        simulate_post_measurement(posterior, incomplete_model, "a", minimum_improvement=0.0)


def test_outcome_serialisation_is_deterministic() -> None:
    posterior = PreMeasurementPosteriorSummary(
        posterior={"b": 0.5, "a": 0.5},
        decision="MEASURE",
        action="measure",
    )
    model = SyntheticObservationModel(
        likelihoods_by_action={
            "measure": {
                "b": {"no": 0.2, "yes": 0.8},
                "a": {"no": 0.8, "yes": 0.2},
            }
        },
        target_metric="posterior_entropy",
    )
    first = simulate_post_measurement(posterior, model, "a", minimum_improvement=0.01)
    second = simulate_post_measurement(posterior, model, "a", minimum_improvement=0.01)

    assert first.to_dict() == second.to_dict()
    assert math.isfinite(first.pre_metric)
