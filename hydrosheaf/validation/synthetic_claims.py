"""Synthetic-only claim gates and a bounded prospective-outcome simulator.

This module is deliberately independent of the programme runner.  It provides
an auditable claim boundary for controlled synthetic evidence and a small
truth-aware simulator for post-measurement outcomes.  It does not turn a
synthetic result into field validation.

The gate is fail-closed:

* thresholds are required explicitly; there are no performance defaults;
* missing held-out, calibration, selective-risk, generator, or comparator
  evidence yields ``ABSTAIN`` rather than a pass;
* comparator audits must declare truth-blind inputs and an independent,
  hashed candidate universe before comparative evidence is adjudicable; and
* field validation is always ``DEFERRED`` and is never part of the synthetic
  gate reduction.

The prospective simulator uses a declared likelihood model and a sealed truth
only to score possible measurement outcomes.  A policy that abstains is not
assigned the pre-measurement score, zero risk, or information gain: it is
reported as unevaluated with ``improved=False``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
import math
from typing import Any


_PROBABILITY_TOLERANCE = 1e-9
_MISSING = object()


class ClaimStatus(str, Enum):
    """Status values used by the synthetic and field claim records."""

    PASS = "PASS"
    FAIL = "FAIL"
    ABSTAIN = "ABSTAIN"
    DEFERRED = "DEFERRED"


class MetricDirection(str, Enum):
    """Direction for a preregistered scalar threshold."""

    AT_LEAST = "at_least"
    AT_MOST = "at_most"


class PolicyDecision(str, Enum):
    """Pre-measurement policy decisions supported by the simulator."""

    MEASURE = "MEASURE"
    ABSTAIN = "ABSTAIN"


class TargetMetric(str, Enum):
    """Lower-is-better posterior metrics supported by the simulator."""

    BAYES_ERROR = "bayes_error"
    POSTERIOR_ENTROPY = "posterior_entropy"


def _text(value: object, *, name: str) -> str:
    if isinstance(value, Enum):
        value = value.value
    result = str(value).strip()
    if not result:
        raise ValueError(f"{name} must be non-empty.")
    return result


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _optional_finite(value: object, *, name: str) -> float | None:
    if value is None:
        return None
    return _finite(value, name=name)


def _normalise_ids(values: object, *, name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be a sequence of IDs, not a string.")
    try:
        raw_values = tuple(values)  # type: ignore[arg-type]
    except TypeError as exc:
        raise ValueError(f"{name} must be a sequence of IDs.") from exc
    result = tuple(sorted({_text(value, name=name) for value in raw_values}))
    if len(result) != len(raw_values):
        raise ValueError(f"{name} must not contain duplicate IDs.")
    return result


def _normalise_optional_metrics(
    values: Mapping[str, object], *, name: str
) -> dict[str, float | None]:
    if not isinstance(values, Mapping):
        raise TypeError(f"{name} must be a mapping.")
    result: dict[str, float | None] = {}
    for key, value in values.items():
        metric_name = _text(key, name=f"{name} key")
        result[metric_name] = _optional_finite(value, name=f"{name}[{metric_name}]")
    return dict(sorted(result.items()))


def _jsonable(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))}
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    return value


@dataclass(frozen=True)
class MetricThreshold:
    """One explicit inclusive threshold for a held-out primary metric."""

    value: float
    direction: str

    def __post_init__(self) -> None:
        value = _finite(self.value, name="MetricThreshold.value")
        direction = _text(self.direction, name="MetricThreshold.direction")
        if direction not in {item.value for item in MetricDirection}:
            raise ValueError(
                "MetricThreshold.direction must be 'at_least' or 'at_most'."
            )
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "direction", direction)

    def meets(self, observed: float) -> bool:
        """Return whether an observed value meets this inclusive threshold."""

        if self.direction == MetricDirection.AT_LEAST.value:
            return float(observed) >= self.value
        return float(observed) <= self.value

    def to_dict(self) -> dict[str, object]:
        return {"value": self.value, "direction": self.direction}


@dataclass(frozen=True)
class GateThresholds:
    """All thresholds and completeness requirements for one claim.

    No field has a default.  A component and an integrated claim should be
    constructed with separate declarations, even when their thresholds happen
    to be equal.
    """

    claim_name: str
    evaluation_split: str
    min_held_out_cases: int
    min_cases_per_generator: int
    required_generator_families: tuple[str, ...]
    required_comparators: tuple[str, ...]
    primary_metric: str
    held_out_metric_thresholds: Mapping[str, MetricThreshold]
    min_calibration_coverage: float
    max_calibration_error: float
    max_mean_interval_width: float
    max_abstention_rate: float
    min_acceptance_rate: float
    max_selective_risk: float
    comparator_metric: str
    max_comparator_degradation: float
    require_truth_blind_comparators: bool
    require_independent_candidate_universes: bool

    def __post_init__(self) -> None:
        claim_name = _text(self.claim_name, name="GateThresholds.claim_name")
        evaluation_split = _text(
            self.evaluation_split, name="GateThresholds.evaluation_split"
        )
        if isinstance(self.min_held_out_cases, bool) or int(self.min_held_out_cases) < 1:
            raise ValueError("min_held_out_cases must be a positive integer.")
        if isinstance(self.min_cases_per_generator, bool) or int(self.min_cases_per_generator) < 1:
            raise ValueError("min_cases_per_generator must be a positive integer.")

        required_generators = _normalise_ids(
            self.required_generator_families,
            name="GateThresholds.required_generator_families",
        )
        required_comparators = _normalise_ids(
            self.required_comparators,
            name="GateThresholds.required_comparators",
        )
        if not required_generators:
            raise ValueError("At least one generator family is required.")
        if not required_comparators:
            raise ValueError("At least one comparator is required.")
        primary_metric = _text(self.primary_metric, name="GateThresholds.primary_metric")
        comparator_metric = _text(
            self.comparator_metric, name="GateThresholds.comparator_metric"
        )
        if not isinstance(self.held_out_metric_thresholds, Mapping):
            raise TypeError("held_out_metric_thresholds must be a mapping.")
        if not self.held_out_metric_thresholds:
            raise ValueError("At least one held-out metric threshold is required.")
        metric_thresholds: dict[str, MetricThreshold] = {}
        for key, value in self.held_out_metric_thresholds.items():
            metric_name = _text(key, name="held-out metric threshold key")
            if not isinstance(value, MetricThreshold):
                raise TypeError(
                    f"Held-out metric threshold {metric_name!r} must be a MetricThreshold."
                )
            metric_thresholds[metric_name] = value
        if primary_metric not in metric_thresholds:
            raise ValueError(
                "primary_metric must have a corresponding held-out metric threshold."
            )

        coverage = _finite(
            self.min_calibration_coverage,
            name="GateThresholds.min_calibration_coverage",
        )
        calibration_error = _finite(
            self.max_calibration_error,
            name="GateThresholds.max_calibration_error",
        )
        interval_width = _finite(
            self.max_mean_interval_width,
            name="GateThresholds.max_mean_interval_width",
        )
        abstention = _finite(
            self.max_abstention_rate,
            name="GateThresholds.max_abstention_rate",
        )
        acceptance = _finite(
            self.min_acceptance_rate,
            name="GateThresholds.min_acceptance_rate",
        )
        selective_risk = _finite(
            self.max_selective_risk,
            name="GateThresholds.max_selective_risk",
        )
        comparator_degradation = _finite(
            self.max_comparator_degradation,
            name="GateThresholds.max_comparator_degradation",
        )
        for value, field_name in (
            (coverage, "min_calibration_coverage"),
            (calibration_error, "max_calibration_error"),
            (abstention, "max_abstention_rate"),
            (acceptance, "min_acceptance_rate"),
            (selective_risk, "max_selective_risk"),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{field_name} must lie in [0, 1].")
        if interval_width < 0.0 or comparator_degradation < 0.0:
            raise ValueError(
                "max_mean_interval_width and max_comparator_degradation "
                "must be non-negative."
            )
        for field_name in (
            "require_truth_blind_comparators",
            "require_independent_candidate_universes",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ValueError(f"{field_name} must be boolean.")

        object.__setattr__(self, "claim_name", claim_name)
        object.__setattr__(self, "evaluation_split", evaluation_split)
        object.__setattr__(self, "min_held_out_cases", int(self.min_held_out_cases))
        object.__setattr__(self, "min_cases_per_generator", int(self.min_cases_per_generator))
        object.__setattr__(self, "required_generator_families", required_generators)
        object.__setattr__(self, "required_comparators", required_comparators)
        object.__setattr__(self, "primary_metric", primary_metric)
        object.__setattr__(
            self,
            "held_out_metric_thresholds",
            dict(sorted(metric_thresholds.items())),
        )
        object.__setattr__(self, "min_calibration_coverage", coverage)
        object.__setattr__(self, "max_calibration_error", calibration_error)
        object.__setattr__(self, "max_mean_interval_width", interval_width)
        object.__setattr__(self, "max_abstention_rate", abstention)
        object.__setattr__(self, "min_acceptance_rate", acceptance)
        object.__setattr__(self, "max_selective_risk", selective_risk)
        object.__setattr__(self, "comparator_metric", comparator_metric)
        object.__setattr__(self, "max_comparator_degradation", comparator_degradation)

    def to_dict(self) -> dict[str, object]:
        return {
            "claim_name": self.claim_name,
            "evaluation_split": self.evaluation_split,
            "min_held_out_cases": self.min_held_out_cases,
            "min_cases_per_generator": self.min_cases_per_generator,
            "required_generator_families": list(self.required_generator_families),
            "required_comparators": list(self.required_comparators),
            "primary_metric": self.primary_metric,
            "held_out_metric_thresholds": {
                name: threshold.to_dict()
                for name, threshold in self.held_out_metric_thresholds.items()
            },
            "min_calibration_coverage": self.min_calibration_coverage,
            "max_calibration_error": self.max_calibration_error,
            "max_mean_interval_width": self.max_mean_interval_width,
            "max_abstention_rate": self.max_abstention_rate,
            "min_acceptance_rate": self.min_acceptance_rate,
            "max_selective_risk": self.max_selective_risk,
            "comparator_metric": self.comparator_metric,
            "max_comparator_degradation": self.max_comparator_degradation,
            "require_truth_blind_comparators": self.require_truth_blind_comparators,
            "require_independent_candidate_universes": self.require_independent_candidate_universes,
        }


@dataclass(frozen=True)
class ComparatorAudit:
    """Fairness and scoring evidence for one synthetic comparator.

    ``independent_candidate_universe`` must be true only when the scored
    universe was declared before inference and is independent of HydroSheaf's
    proposed edges.  A node-level age comparator can satisfy this using an
    explicit all-eligible-node universe; a reaction comparator scored only on
    HydroSheaf's candidate edges must leave it false.
    """

    complete: bool
    truth_blind: bool
    independent_candidate_universe: bool
    candidate_universe_scope: str
    candidate_universe_hash: str | None
    input_channels: tuple[str, ...]
    metrics: Mapping[str, float | None]

    def __post_init__(self) -> None:
        for field_name in ("complete", "truth_blind", "independent_candidate_universe"):
            if not isinstance(getattr(self, field_name), bool):
                raise ValueError(f"ComparatorAudit.{field_name} must be boolean.")
        scope = _text(
            self.candidate_universe_scope,
            name="ComparatorAudit.candidate_universe_scope",
        )
        channels = _normalise_ids(
            self.input_channels,
            name="ComparatorAudit.input_channels",
        )
        universe_hash = None
        if self.candidate_universe_hash is not None:
            universe_hash = _text(
                self.candidate_universe_hash,
                name="ComparatorAudit.candidate_universe_hash",
            )
        metrics = _normalise_optional_metrics(
            self.metrics,
            name="ComparatorAudit.metrics",
        )
        object.__setattr__(self, "candidate_universe_scope", scope)
        object.__setattr__(self, "candidate_universe_hash", universe_hash)
        object.__setattr__(self, "input_channels", channels)
        object.__setattr__(self, "metrics", metrics)

    def to_dict(self) -> dict[str, object]:
        return {
            "complete": self.complete,
            "truth_blind": self.truth_blind,
            "independent_candidate_universe": self.independent_candidate_universe,
            "candidate_universe_scope": self.candidate_universe_scope,
            "candidate_universe_hash": self.candidate_universe_hash,
            "input_channels": list(self.input_channels),
            "metrics": dict(self.metrics),
        }


@dataclass(frozen=True)
class HeldOutEvidence:
    """Held-out evidence consumed by one synthetic claim gate.

    Generator counts are case counts and must sum to ``held_out_cases``.  This
    prevents a report from claiming family coverage with an unaccounted or
    duplicated denominator.  Calibration and selective-risk values are
    explicit scalar summaries; use ``None`` when a required result was not
    produced so the gate returns ``ABSTAIN``.
    """

    evaluation_split: str
    held_out_cases: int
    generator_case_counts: Mapping[str, int]
    held_out_metrics: Mapping[str, float | None]
    calibration_coverage: float | None
    calibration_error: float | None
    mean_interval_width: float | None
    abstention_rate: float | None
    acceptance_rate: float | None
    selective_risk: float | None
    comparators: Mapping[str, ComparatorAudit]

    def __post_init__(self) -> None:
        evaluation_split = _text(
            self.evaluation_split, name="HeldOutEvidence.evaluation_split"
        )
        if isinstance(self.held_out_cases, bool) or int(self.held_out_cases) < 0:
            raise ValueError("HeldOutEvidence.held_out_cases must be non-negative.")
        if not isinstance(self.generator_case_counts, Mapping):
            raise TypeError("HeldOutEvidence.generator_case_counts must be a mapping.")
        counts: dict[str, int] = {}
        for key, value in self.generator_case_counts.items():
            family = _text(key, name="generator family")
            if isinstance(value, bool) or int(value) < 0:
                raise ValueError(f"Generator count for {family!r} must be non-negative.")
            counts[family] = int(value)
        metrics = _normalise_optional_metrics(
            self.held_out_metrics,
            name="HeldOutEvidence.held_out_metrics",
        )
        scalar_names = (
            "calibration_coverage",
            "calibration_error",
            "mean_interval_width",
            "abstention_rate",
            "acceptance_rate",
            "selective_risk",
        )
        scalars = {
            name: _optional_finite(getattr(self, name), name=f"HeldOutEvidence.{name}")
            for name in scalar_names
        }
        if not isinstance(self.comparators, Mapping):
            raise TypeError("HeldOutEvidence.comparators must be a mapping.")
        comparators: dict[str, ComparatorAudit] = {}
        for key, value in self.comparators.items():
            comparator_name = _text(key, name="comparator name")
            if not isinstance(value, ComparatorAudit):
                raise TypeError(
                    f"Comparator {comparator_name!r} must be a ComparatorAudit."
                )
            comparators[comparator_name] = value
        object.__setattr__(self, "evaluation_split", evaluation_split)
        object.__setattr__(self, "held_out_cases", int(self.held_out_cases))
        object.__setattr__(self, "generator_case_counts", dict(sorted(counts.items())))
        object.__setattr__(self, "held_out_metrics", metrics)
        for name, value in scalars.items():
            object.__setattr__(self, name, value)
        object.__setattr__(self, "comparators", dict(sorted(comparators.items())))

    def to_dict(self) -> dict[str, object]:
        return {
            "evaluation_split": self.evaluation_split,
            "held_out_cases": self.held_out_cases,
            "generator_case_counts": dict(self.generator_case_counts),
            "held_out_metrics": dict(self.held_out_metrics),
            "calibration_coverage": self.calibration_coverage,
            "calibration_error": self.calibration_error,
            "mean_interval_width": self.mean_interval_width,
            "abstention_rate": self.abstention_rate,
            "acceptance_rate": self.acceptance_rate,
            "selective_risk": self.selective_risk,
            "comparators": {
                name: audit.to_dict() for name, audit in self.comparators.items()
            },
        }


@dataclass(frozen=True)
class GateCheck:
    """One auditable preregistered check."""

    name: str
    passed: bool | None
    observed: object
    requirement: object
    detail: str

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "passed": self.passed,
            "observed": _jsonable(self.observed),
            "requirement": _jsonable(self.requirement),
            "detail": self.detail,
        }


@dataclass(frozen=True)
class GateResult:
    """Result for one claim, including every check that determined it."""

    claim_name: str
    status: ClaimStatus
    checks: tuple[GateCheck, ...]
    reasons: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "claim_name": self.claim_name,
            "status": self.status.value,
            "checks": [check.to_dict() for check in self.checks],
            "reasons": list(self.reasons),
        }


@dataclass(frozen=True)
class SyntheticClaimReport:
    """Separate component, integrated, and field-validation claim statuses."""

    synthetic_component_claim: GateResult
    synthetic_integrated_claim: GateResult
    field_validation_claim: GateResult

    def __post_init__(self) -> None:
        if self.field_validation_claim.status is not ClaimStatus.DEFERRED:
            raise ValueError("field_validation_claim must always be DEFERRED.")

    @property
    def synthetic_status(self) -> ClaimStatus:
        """Reduce only the two synthetic statuses; field is intentionally ignored."""

        statuses = (
            self.synthetic_component_claim.status,
            self.synthetic_integrated_claim.status,
        )
        if ClaimStatus.FAIL in statuses:
            return ClaimStatus.FAIL
        if ClaimStatus.ABSTAIN in statuses:
            return ClaimStatus.ABSTAIN
        if all(status == ClaimStatus.PASS for status in statuses):
            return ClaimStatus.PASS
        return ClaimStatus.ABSTAIN

    def to_dict(self) -> dict[str, object]:
        return {
            "synthetic_component_claim": self.synthetic_component_claim.to_dict(),
            "synthetic_integrated_claim": self.synthetic_integrated_claim.to_dict(),
            "field_validation_claim": self.field_validation_claim.to_dict(),
            "synthetic_status": self.synthetic_status.value,
        }


def _check(
    name: str,
    passed: bool | None,
    observed: object,
    requirement: object,
    detail: str,
) -> GateCheck:
    return GateCheck(
        name=name,
        passed=passed,
        observed=observed,
        requirement=requirement,
        detail=detail,
    )


def _bounded_check(
    value: float | None,
    *,
    lower: float | None = None,
    upper: float | None = None,
) -> bool | None:
    if value is None or not math.isfinite(value):
        return None
    if lower is not None and value < lower:
        return False
    if upper is not None and value > upper:
        return False
    return True


def _evaluate_gate(
    claim_name: str,
    evidence: HeldOutEvidence | None,
    thresholds: GateThresholds,
) -> GateResult:
    if evidence is None:
        check = _check(
            "held_out_evidence_present",
            None,
            None,
            True,
            "held-out evidence is missing; the claim is not adjudicable",
        )
        return GateResult(claim_name, ClaimStatus.ABSTAIN, (check,), (check.detail,))

    checks: list[GateCheck] = []
    split_passed = evidence.evaluation_split == thresholds.evaluation_split
    checks.append(
        _check(
            "evaluation_split",
            split_passed,
            evidence.evaluation_split,
            thresholds.evaluation_split,
            (
                "held-out split matches the preregistration"
                if split_passed
                else "evaluation split does not match the preregistered split"
            ),
        )
    )

    case_count_passed = evidence.held_out_cases >= thresholds.min_held_out_cases
    checks.append(
        _check(
            "held_out_case_count",
            case_count_passed,
            evidence.held_out_cases,
            {"minimum": thresholds.min_held_out_cases},
            (
                "held-out denominator meets the declared minimum"
                if case_count_passed
                else "held-out denominator is below the declared minimum"
            ),
        )
    )
    accounted_cases = sum(evidence.generator_case_counts.values())
    accounting_passed = accounted_cases == evidence.held_out_cases
    checks.append(
        _check(
            "generator_case_accounting",
            accounting_passed,
            accounted_cases,
            {"held_out_cases": evidence.held_out_cases},
            (
                "generator counts account for the held-out denominator"
                if accounting_passed
                else "generator counts do not account for the held-out denominator"
            ),
        )
    )
    for family in thresholds.required_generator_families:
        observed_count = evidence.generator_case_counts.get(family, 0)
        family_passed = observed_count >= thresholds.min_cases_per_generator
        checks.append(
            _check(
                f"generator_coverage[{family}]",
                family_passed,
                observed_count,
                {"minimum_cases": thresholds.min_cases_per_generator},
                (
                    f"generator family {family!r} meets the declared coverage"
                    if family_passed
                    else f"generator family {family!r} is under-covered"
                ),
            )
        )

    primary_value = evidence.held_out_metrics.get(thresholds.primary_metric, _MISSING)
    for metric_name, metric_threshold in thresholds.held_out_metric_thresholds.items():
        metric_value = evidence.held_out_metrics.get(metric_name, _MISSING)
        metric_observed = None if metric_value is _MISSING else metric_value
        metric_passed: bool | None
        if metric_value is _MISSING or metric_value is None:
            metric_passed = None
        else:
            metric_passed = metric_threshold.meets(float(metric_value))
        metric_label = "primary" if metric_name == thresholds.primary_metric else "required"
        checks.append(
            _check(
                f"held_out_metric[{metric_name}]",
                metric_passed,
                metric_observed,
                metric_threshold.to_dict(),
                (
                    f"held-out {metric_label} metric meets the declared threshold"
                    if metric_passed is True
                    else (
                        f"held-out {metric_label} metric misses the declared threshold"
                        if metric_passed is False
                        else f"held-out {metric_label} metric is missing or non-finite"
                    )
                ),
            )
        )

    coverage_passed = _bounded_check(
        evidence.calibration_coverage,
        lower=thresholds.min_calibration_coverage,
        upper=1.0,
    )
    checks.append(
        _check(
            "calibration_coverage",
            coverage_passed,
            evidence.calibration_coverage,
            {"minimum": thresholds.min_calibration_coverage},
            (
                "calibration coverage meets the declared minimum"
                if coverage_passed is True
                else (
                    "calibration coverage misses the declared minimum"
                    if coverage_passed is False
                    else "calibration coverage is missing or invalid"
                )
            ),
        )
    )
    calibration_error_passed = _bounded_check(
        evidence.calibration_error,
        lower=0.0,
        upper=thresholds.max_calibration_error,
    )
    checks.append(
        _check(
            "calibration_error",
            calibration_error_passed,
            evidence.calibration_error,
            {"maximum": thresholds.max_calibration_error},
            (
                "calibration error is within the declared maximum"
                if calibration_error_passed is True
                else (
                    "calibration error exceeds the declared maximum"
                    if calibration_error_passed is False
                    else "calibration error is missing or invalid"
                )
            ),
        )
    )
    interval_width_passed = _bounded_check(
        evidence.mean_interval_width,
        lower=0.0,
        upper=thresholds.max_mean_interval_width,
    )
    checks.append(
        _check(
            "mean_interval_width",
            interval_width_passed,
            evidence.mean_interval_width,
            {"maximum": thresholds.max_mean_interval_width},
            (
                "mean interval width is within the declared maximum"
                if interval_width_passed is True
                else (
                    "mean interval width exceeds the declared maximum"
                    if interval_width_passed is False
                    else "mean interval width is missing or invalid"
                )
            ),
        )
    )

    rate_values = (evidence.abstention_rate, evidence.acceptance_rate)
    rate_domain_passed = all(
        value is not None and math.isfinite(value) and 0.0 <= value <= 1.0
        for value in rate_values
    )
    if not rate_domain_passed:
        rate_consistency: bool | None = None
    else:
        rate_consistency = math.isclose(
            float(evidence.abstention_rate) + float(evidence.acceptance_rate),
            1.0,
            rel_tol=0.0,
            abs_tol=_PROBABILITY_TOLERANCE,
        )
    checks.append(
        _check(
            "acceptance_abstention_accounting",
            rate_consistency,
            {
                "acceptance_rate": evidence.acceptance_rate,
                "abstention_rate": evidence.abstention_rate,
            },
            {"sum": 1.0},
            (
                "acceptance and abstention rates share one denominator"
                if rate_consistency is True
                else (
                    "acceptance and abstention rates do not sum to one"
                    if rate_consistency is False
                    else "acceptance or abstention rate is missing or invalid"
                )
            ),
        )
    )
    abstention_passed = _bounded_check(
        evidence.abstention_rate,
        lower=0.0,
        upper=thresholds.max_abstention_rate,
    )
    checks.append(
        _check(
            "abstention_rate",
            abstention_passed,
            evidence.abstention_rate,
            {"maximum": thresholds.max_abstention_rate},
            (
                "abstention rate is within the declared maximum"
                if abstention_passed is True
                else (
                    "abstention rate exceeds the declared maximum"
                    if abstention_passed is False
                    else "abstention rate is missing or invalid"
                )
            ),
        )
    )
    acceptance_passed = _bounded_check(
        evidence.acceptance_rate,
        lower=thresholds.min_acceptance_rate,
        upper=1.0,
    )
    checks.append(
        _check(
            "acceptance_rate",
            acceptance_passed,
            evidence.acceptance_rate,
            {"minimum": thresholds.min_acceptance_rate},
            (
                "acceptance rate meets the declared minimum"
                if acceptance_passed is True
                else (
                    "acceptance rate is below the declared minimum"
                    if acceptance_passed is False
                    else "acceptance rate is missing or invalid"
                )
            ),
        )
    )
    selective_risk_value = evidence.selective_risk
    if evidence.acceptance_rate is not None and evidence.acceptance_rate <= 0.0:
        selective_risk_value = None
    selective_risk_passed = _bounded_check(
        selective_risk_value,
        lower=0.0,
        upper=thresholds.max_selective_risk,
    )
    checks.append(
        _check(
            "selective_risk",
            selective_risk_passed,
            selective_risk_value,
            {"maximum": thresholds.max_selective_risk},
            (
                "selective risk is within the declared maximum"
                if selective_risk_passed is True
                else (
                    "selective risk is unavailable without accepted cases or exceeds the maximum"
                    if selective_risk_passed is False
                    else "selective risk is missing, invalid, or has no accepted cases"
                )
            ),
        )
    )

    for comparator_name in thresholds.required_comparators:
        audit = evidence.comparators.get(comparator_name)
        comparator_check_name = f"comparator_completeness[{comparator_name}]"
        if audit is None:
            checks.append(
                _check(
                    comparator_check_name,
                    None,
                    None,
                    True,
                    f"required comparator {comparator_name!r} is missing",
                )
            )
            continue
        audit_reasons: list[str] = []
        if not audit.complete:
            audit_reasons.append("output is incomplete")
        if thresholds.require_truth_blind_comparators and not audit.truth_blind:
            audit_reasons.append("truth-blind input boundary is not satisfied")
        if thresholds.require_independent_candidate_universes:
            if not audit.independent_candidate_universe:
                audit_reasons.append("candidate universe is not independent")
            if not audit.candidate_universe_hash:
                audit_reasons.append("candidate-universe hash is missing")
        if not audit.input_channels:
            audit_reasons.append("declared input channels are missing")
        comparator_metric_value = audit.metrics.get(thresholds.comparator_metric, _MISSING)
        if comparator_metric_value is _MISSING or comparator_metric_value is None:
            audit_reasons.append(
                f"comparator metric {thresholds.comparator_metric!r} is missing"
            )
        if primary_value is _MISSING or primary_value is None:
            audit_reasons.append(
                f"candidate metric {thresholds.primary_metric!r} is missing"
            )
        if audit_reasons:
            checks.append(
                _check(
                    comparator_check_name,
                    None,
                    {"reasons": tuple(audit_reasons)},
                    {
                        "complete": True,
                        "truth_blind": thresholds.require_truth_blind_comparators,
                        "independent_candidate_universe": thresholds.require_independent_candidate_universes,
                        "metric": thresholds.comparator_metric,
                    },
                    f"comparator {comparator_name!r} is not adjudicable: "
                    + "; ".join(audit_reasons),
                )
            )
            continue
        candidate_metric = float(primary_value)
        comparator_metric = float(comparator_metric_value)
        primary_threshold = thresholds.held_out_metric_thresholds[thresholds.primary_metric]
        if primary_threshold.direction == MetricDirection.AT_LEAST.value:
            comparator_passed = candidate_metric >= comparator_metric - thresholds.max_comparator_degradation
        else:
            comparator_passed = candidate_metric <= comparator_metric + thresholds.max_comparator_degradation
        checks.append(
            _check(
                comparator_check_name,
                comparator_passed,
                {
                    "candidate": candidate_metric,
                    "comparator": comparator_metric,
                },
                {"maximum_degradation": thresholds.max_comparator_degradation},
                (
                    f"comparator {comparator_name!r} is complete and the candidate is within the locked margin"
                    if comparator_passed
                    else f"candidate is outside the locked non-inferiority margin versus {comparator_name!r}"
                ),
            )
        )

    if any(check.passed is False for check in checks):
        status = ClaimStatus.FAIL
    elif any(check.passed is None for check in checks):
        status = ClaimStatus.ABSTAIN
    else:
        status = ClaimStatus.PASS
    reasons = tuple(check.detail for check in checks if check.passed is not True)
    return GateResult(claim_name, status, tuple(checks), reasons)


def evaluate_synthetic_claims(
    component_evidence: HeldOutEvidence | None,
    integrated_evidence: HeldOutEvidence | None,
    component_thresholds: GateThresholds,
    integrated_thresholds: GateThresholds,
) -> SyntheticClaimReport:
    """Evaluate separate component and integrated synthetic claim gates.

    ``field_validation_claim`` is created internally as ``DEFERRED``.  It has
    no evidence argument by design, so field readiness cannot accidentally
    become a synthetic pass/fail prerequisite.
    """

    component_name = "synthetic_component_claim"
    integrated_name = "synthetic_integrated_claim"
    if component_thresholds.claim_name != component_name:
        raise ValueError(
            f"component_thresholds.claim_name must be {component_name!r}."
        )
    if integrated_thresholds.claim_name != integrated_name:
        raise ValueError(
            f"integrated_thresholds.claim_name must be {integrated_name!r}."
        )
    component = _evaluate_gate(component_name, component_evidence, component_thresholds)
    integrated = _evaluate_gate(integrated_name, integrated_evidence, integrated_thresholds)
    field_check = GateCheck(
        name="field_validation_data",
        passed=None,
        observed=None,
        requirement="prospective field campaign with independent outcome measurements",
        detail="field validation is deferred until a preregistered field campaign exists",
    )
    field = GateResult(
        claim_name="field_validation_claim",
        status=ClaimStatus.DEFERRED,
        checks=(field_check,),
        reasons=(field_check.detail,),
    )
    return SyntheticClaimReport(
        synthetic_component_claim=component,
        synthetic_integrated_claim=integrated,
        field_validation_claim=field,
    )


def _normalise_distribution(
    values: Mapping[str, object], *, name: str
) -> dict[str, float]:
    if not isinstance(values, Mapping) or not values:
        raise ValueError(f"{name} must be a non-empty mapping.")
    result: dict[str, float] = {}
    for key, value in values.items():
        state = _text(key, name=f"{name} key")
        probability = _finite(value, name=f"{name}[{state}]")
        if probability < 0.0:
            raise ValueError(f"{name}[{state}] must be non-negative.")
        result[state] = probability
    total = sum(result.values())
    if total <= 0.0 or not math.isclose(
        total, 1.0, rel_tol=_PROBABILITY_TOLERANCE, abs_tol=_PROBABILITY_TOLERANCE
    ):
        raise ValueError(f"{name} probabilities must sum to one.")
    return {
        key: value / total for key, value in sorted(result.items())
    }


@dataclass(frozen=True)
class PreMeasurementPosteriorSummary:
    """Truth-free posterior and policy decision available before measurement."""

    posterior: Mapping[str, float]
    decision: str
    action: str | None

    def __post_init__(self) -> None:
        posterior = _normalise_distribution(
            self.posterior, name="PreMeasurementPosteriorSummary.posterior"
        )
        decision = _text(self.decision, name="PreMeasurementPosteriorSummary.decision")
        if decision not in {item.value for item in PolicyDecision}:
            raise ValueError("decision must be 'MEASURE' or 'ABSTAIN'.")
        action = None
        if self.action is not None:
            action = _text(self.action, name="PreMeasurementPosteriorSummary.action")
        if decision == PolicyDecision.MEASURE.value and action is None:
            raise ValueError("a MEASURE decision requires an action.")
        if decision == PolicyDecision.ABSTAIN.value and action is not None:
            raise ValueError("an ABSTAIN decision cannot carry an action.")
        object.__setattr__(self, "posterior", posterior)
        object.__setattr__(self, "decision", decision)
        object.__setattr__(self, "action", action)

    def to_dict(self) -> dict[str, object]:
        return {
            "posterior": dict(self.posterior),
            "decision": self.decision,
            "action": self.action,
        }


@dataclass(frozen=True)
class SyntheticObservationModel:
    """Declared action-specific likelihoods ``P(outcome | state, action)``."""

    likelihoods_by_action: Mapping[str, Mapping[str, Mapping[str, float]]]
    target_metric: str

    def __post_init__(self) -> None:
        if not isinstance(self.likelihoods_by_action, Mapping) or not self.likelihoods_by_action:
            raise ValueError("likelihoods_by_action must be non-empty.")
        actions: dict[str, dict[str, dict[str, float]]] = {}
        for raw_action, raw_by_state in self.likelihoods_by_action.items():
            action = _text(raw_action, name="observation action")
            if not isinstance(raw_by_state, Mapping) or not raw_by_state:
                raise ValueError(f"Action {action!r} needs state likelihoods.")
            by_state: dict[str, dict[str, float]] = {}
            for raw_state, raw_by_outcome in raw_by_state.items():
                state = _text(raw_state, name="observation state")
                by_state[state] = _normalise_distribution(
                    raw_by_outcome,
                    name=f"likelihoods_by_action[{action}][{state}]",
                )
            actions[action] = dict(sorted(by_state.items()))
        metric = _text(self.target_metric, name="SyntheticObservationModel.target_metric")
        if metric not in {item.value for item in TargetMetric}:
            raise ValueError(
                "target_metric must be 'bayes_error' or 'posterior_entropy'."
            )
        object.__setattr__(self, "likelihoods_by_action", dict(sorted(actions.items())))
        object.__setattr__(self, "target_metric", metric)

    def to_dict(self) -> dict[str, object]:
        return {
            "likelihoods_by_action": _jsonable(self.likelihoods_by_action),
            "target_metric": self.target_metric,
        }


def _posterior_metric(posterior: Mapping[str, float], metric: str) -> float:
    values = tuple(float(value) for value in posterior.values())
    if metric == TargetMetric.BAYES_ERROR.value:
        return 1.0 - max(values)
    return -sum(value * math.log(value) for value in values if value > 0.0)


@dataclass(frozen=True)
class PostMeasurementOutcome:
    """One possible truth-conditioned measurement outcome and updated posterior."""

    outcome: str
    probability_given_truth: float
    posterior: Mapping[str, float]
    post_metric: float
    improvement: float

    def to_dict(self) -> dict[str, object]:
        return {
            "outcome": self.outcome,
            "probability_given_truth": self.probability_given_truth,
            "posterior": dict(self.posterior),
            "post_metric": self.post_metric,
            "improvement": self.improvement,
        }


@dataclass(frozen=True)
class PostMeasurementEvaluation:
    """Truth-aware expected post-measurement result for one policy decision."""

    decision: str
    action: str | None
    evaluated: bool
    improved: bool
    reason: str
    true_hypothesis: str
    target_metric: str
    minimum_improvement: float
    pre_metric: float
    expected_post_metric: float | None
    expected_improvement: float | None
    outcomes: tuple[PostMeasurementOutcome, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "decision": self.decision,
            "action": self.action,
            "evaluated": self.evaluated,
            "improved": self.improved,
            "reason": self.reason,
            "true_hypothesis": self.true_hypothesis,
            "target_metric": self.target_metric,
            "minimum_improvement": self.minimum_improvement,
            "pre_metric": self.pre_metric,
            "expected_post_metric": self.expected_post_metric,
            "expected_improvement": self.expected_improvement,
            "outcomes": [outcome.to_dict() for outcome in self.outcomes],
        }


def simulate_post_measurement(
    posterior_summary: PreMeasurementPosteriorSummary,
    observation_model: SyntheticObservationModel,
    true_hypothesis: str,
    *,
    minimum_improvement: float,
) -> PostMeasurementEvaluation:
    """Evaluate truth-conditioned post-measurement improvement.

    The policy sees only ``posterior_summary``.  ``true_hypothesis`` is used
    after the decision to weight declared observation outcomes, never to choose
    the action.  The returned expected improvement is ``pre_metric -
    expected_post_metric`` for a lower-is-better metric.  Improvement must be
    strictly greater than the explicit ``minimum_improvement``; equality is
    not a pass.
    """

    minimum = _finite(minimum_improvement, name="minimum_improvement")
    if minimum < 0.0:
        raise ValueError("minimum_improvement must be non-negative.")
    truth = _text(true_hypothesis, name="true_hypothesis")
    if truth not in posterior_summary.posterior:
        raise ValueError("true_hypothesis must be present in the pre-measurement posterior.")
    if posterior_summary.posterior[truth] <= 0.0:
        raise ValueError("true_hypothesis must have positive prior posterior mass.")

    pre_metric = _posterior_metric(
        posterior_summary.posterior,
        observation_model.target_metric,
    )
    if posterior_summary.decision == PolicyDecision.ABSTAIN.value:
        return PostMeasurementEvaluation(
            decision=posterior_summary.decision,
            action=None,
            evaluated=False,
            improved=False,
            reason="policy_abstained_no_post_measurement_score",
            true_hypothesis=truth,
            target_metric=observation_model.target_metric,
            minimum_improvement=minimum,
            pre_metric=pre_metric,
            expected_post_metric=None,
            expected_improvement=None,
            outcomes=(),
        )

    assert posterior_summary.action is not None
    if posterior_summary.action not in observation_model.likelihoods_by_action:
        raise ValueError(
            f"No observation model is declared for action {posterior_summary.action!r}."
        )
    by_state = observation_model.likelihoods_by_action[posterior_summary.action]
    missing_states = set(posterior_summary.posterior) - set(by_state)
    if missing_states:
        raise ValueError(
            "Observation model is missing posterior states: "
            + ", ".join(sorted(missing_states))
        )
    if truth not in by_state:
        raise ValueError(f"Observation model is missing true hypothesis {truth!r}.")

    outcome_ids = sorted(
        {
            outcome
            for state in posterior_summary.posterior
            for outcome in by_state[state]
        }
    )
    outcomes: list[PostMeasurementOutcome] = []
    for outcome in outcome_ids:
        probability_given_truth = by_state[truth].get(outcome, 0.0)
        if probability_given_truth <= 0.0:
            continue
        evidence_mass = sum(
            posterior_summary.posterior[state]
            * by_state[state].get(outcome, 0.0)
            for state in posterior_summary.posterior
        )
        if evidence_mass <= 0.0 or not math.isfinite(evidence_mass):
            raise ValueError(
                f"Outcome {outcome!r} has positive truth likelihood but zero prior evidence mass."
            )
        updated = {
            state: (
                posterior_summary.posterior[state]
                * by_state[state].get(outcome, 0.0)
                / evidence_mass
            )
            for state in posterior_summary.posterior
        }
        post_metric = _posterior_metric(updated, observation_model.target_metric)
        outcomes.append(
            PostMeasurementOutcome(
                outcome=outcome,
                probability_given_truth=probability_given_truth,
                posterior=dict(sorted(updated.items())),
                post_metric=post_metric,
                improvement=pre_metric - post_metric,
            )
        )
    if not outcomes:
        raise ValueError("The true hypothesis assigns zero probability to every outcome.")

    expected_post_metric = sum(
        outcome.probability_given_truth * outcome.post_metric for outcome in outcomes
    )
    expected_improvement = pre_metric - expected_post_metric
    improved = expected_improvement > minimum
    reason = (
        "target_metric_improved_after_measurement"
        if improved
        else "target_metric_did_not_meet_improvement_threshold"
    )
    return PostMeasurementEvaluation(
        decision=posterior_summary.decision,
        action=posterior_summary.action,
        evaluated=True,
        improved=improved,
        reason=reason,
        true_hypothesis=truth,
        target_metric=observation_model.target_metric,
        minimum_improvement=minimum,
        pre_metric=pre_metric,
        expected_post_metric=expected_post_metric,
        expected_improvement=expected_improvement,
        outcomes=tuple(outcomes),
    )


__all__ = [
    "ClaimStatus",
    "ComparatorAudit",
    "GateCheck",
    "GateResult",
    "GateThresholds",
    "HeldOutEvidence",
    "MetricDirection",
    "MetricThreshold",
    "PolicyDecision",
    "PostMeasurementEvaluation",
    "PostMeasurementOutcome",
    "PreMeasurementPosteriorSummary",
    "SyntheticClaimReport",
    "SyntheticObservationModel",
    "TargetMetric",
    "evaluate_synthetic_claims",
    "simulate_post_measurement",
]
