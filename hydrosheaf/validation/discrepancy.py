"""Conservative model-discrepancy and disagreement reporting.

This module does not pretend that a scenario envelope is a calibrated
posterior.  It combines frozen scenario outputs for reporting, exposes the
between-scenario spread, and forces an explicit ``MODEL_DISAGREEMENT`` /
``ABSTAIN`` result when the scenarios are materially incompatible.  Scenario
weights must therefore be supplied by a predeclared development rule; locked
test truth is never accepted by this layer.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib
import json
import math
from typing import Iterable, Mapping, Sequence

from .programme_contract import (
    DecisionKind,
    IdentifiabilityStatus,
    ProgrammeDecision,
)


def _finite(value: object, *, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return result


def _truth_fields(record: Mapping[str, object]) -> tuple[str, ...]:
    return tuple(
        sorted(
            str(key).lower()
            for key in record
            if str(key).lower() in {"truth", "true_value", "ground_truth"}
            or str(key).lower().startswith(("true_", "truth_"))
        )
    )


def _json_hash(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ScenarioPrediction:
    """One prediction from one declared discrepancy scenario."""

    scenario_id: str
    estimate: float
    lower: float
    upper: float
    weight: float = 1.0
    scenario_family: str = "unspecified"
    discrepancy_tags: Sequence[str] = ()
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        scenario_id = str(self.scenario_id).strip()
        family = str(self.scenario_family).strip()
        if not scenario_id:
            raise ValueError("Scenario predictions require scenario_id.")
        if not family:
            raise ValueError("Scenario predictions require scenario_family.")
        estimate = _finite(self.estimate, name="estimate")
        lower = _finite(self.lower, name="lower")
        upper = _finite(self.upper, name="upper")
        weight = _finite(self.weight, name="weight")
        if lower > upper:
            raise ValueError("Scenario lower bound cannot exceed upper bound.")
        if not lower <= estimate <= upper:
            raise ValueError("Scenario estimate must lie within its interval.")
        if weight < 0.0:
            raise ValueError("Scenario weight must be non-negative.")
        tags = tuple(str(tag) for tag in self.discrepancy_tags)
        object.__setattr__(self, "scenario_id", scenario_id)
        object.__setattr__(self, "scenario_family", family)
        object.__setattr__(self, "estimate", estimate)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "discrepancy_tags", tags)

    def to_dict(self) -> dict[str, object]:
        return {
            "scenario_id": self.scenario_id,
            "scenario_family": self.scenario_family,
            "estimate": self.estimate,
            "lower": self.lower,
            "upper": self.upper,
            "weight": self.weight,
            "discrepancy_tags": list(self.discrepancy_tags),
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class DiscrepancyReport:
    """Scenario-combined report with an explicit safe decision."""

    target: str
    estimate: float
    lower: float
    upper: float
    scenario_range: float
    relative_disagreement: float
    intervals_overlap: bool
    identifiability: IdentifiabilityStatus
    decision: ProgrammeDecision
    scenario_weights: Mapping[str, float]
    scenarios: Sequence[ScenarioPrediction]
    interval_kind: str = "conservative_scenario_envelope_not_calibrated"
    calibration_status: str = "NOT_CALIBRATED"
    weight_source: str = "DECLARED_SCENARIO_WEIGHTS"
    model_disagreement_status: str = "UNKNOWN"
    identifiability_reason: str = ""
    scenario_families: tuple[str, ...] = ()
    calibration_scope: str = "none"
    calibration_fit_data_hash: str = ""
    calibration_factor: float | None = None
    calibration_target_coverage: float | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "target": self.target,
            "estimate": self.estimate,
            "lower": self.lower,
            "upper": self.upper,
            "scenario_range": self.scenario_range,
            "relative_disagreement": self.relative_disagreement,
            "intervals_overlap": self.intervals_overlap,
            "identifiability": self.identifiability.value,
            "decision": self.decision.to_dict(),
            "scenario_weights": dict(self.scenario_weights),
            "scenarios": [scenario.to_dict() for scenario in self.scenarios],
            "interval_kind": self.interval_kind,
            "calibration_status": self.calibration_status,
            "calibrated": self.calibration_status == "CALIBRATED",
            "weight_source": self.weight_source,
            "model_disagreement_status": self.model_disagreement_status,
            "identifiability_reason": self.identifiability_reason,
            "scenario_families": list(self.scenario_families),
            "scenario_family_count": len(self.scenario_families),
            "calibration_scope": self.calibration_scope,
            "calibration_fit_data_hash": self.calibration_fit_data_hash,
            "calibration_factor": self.calibration_factor,
            "calibration_target_coverage": self.calibration_target_coverage,
        }


@dataclass(frozen=True)
class DiscrepancyCalibrationObservation:
    """Truth-labelled interval used only for development calibration or scoring."""

    case_id: str
    target_id: str
    truth: float
    estimate: float
    lower: float
    upper: float
    phase: str = "development"

    def __post_init__(self) -> None:
        case_id = str(self.case_id).strip()
        target_id = str(self.target_id).strip()
        phase = str(self.phase).strip()
        if not case_id or not target_id or not phase:
            raise ValueError("Discrepancy calibration observations require IDs and phase.")
        truth = _finite(self.truth, name="truth")
        estimate = _finite(self.estimate, name="estimate")
        lower = _finite(self.lower, name="lower")
        upper = _finite(self.upper, name="upper")
        if lower > upper or not lower <= estimate <= upper:
            raise ValueError("Discrepancy calibration interval must contain its estimate.")
        object.__setattr__(self, "case_id", case_id)
        object.__setattr__(self, "target_id", target_id)
        object.__setattr__(self, "truth", truth)
        object.__setattr__(self, "estimate", estimate)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        object.__setattr__(self, "phase", phase)

    def to_dict(self, *, include_truth: bool = True) -> dict[str, object]:
        record: dict[str, object] = {
            "case_id": self.case_id,
            "target_id": self.target_id,
            "estimate": self.estimate,
            "lower": self.lower,
            "upper": self.upper,
            "phase": self.phase,
        }
        if include_truth:
            record["truth"] = self.truth
        return record


@dataclass(frozen=True)
class DiscrepancyCalibrator:
    """Development-only conservative dilation of discrepancy intervals."""

    target_coverage: float
    scale_factor: float
    fit_scope: str
    case_ids: tuple[str, ...]
    observation_count: int
    fit_data_hash: str
    calibration_kind: str = "case_blocked_conformal_interval_scale_v1"

    def __post_init__(self) -> None:
        coverage = _finite(self.target_coverage, name="target_coverage")
        factor = _finite(self.scale_factor, name="scale_factor")
        if not 0.0 < coverage < 1.0:
            raise ValueError("target_coverage must lie strictly between 0 and 1.")
        if factor < 1.0:
            raise ValueError("scale_factor must not contract a conservative interval.")
        if self.fit_scope != "development_only":
            raise ValueError("Discrepancy calibration may only use development data.")
        cases = tuple(sorted(str(case) for case in self.case_ids))
        if not cases or int(self.observation_count) <= 0:
            raise ValueError("A discrepancy calibrator requires development observations.")
        if len(str(self.fit_data_hash)) != 64:
            raise ValueError("fit_data_hash must be a SHA-256 digest.")
        object.__setattr__(self, "target_coverage", coverage)
        object.__setattr__(self, "scale_factor", factor)
        object.__setattr__(self, "case_ids", cases)
        object.__setattr__(self, "observation_count", int(self.observation_count))

    def to_dict(self) -> dict[str, object]:
        return {
            "target_coverage": self.target_coverage,
            "scale_factor": self.scale_factor,
            "fit_scope": self.fit_scope,
            "case_ids": list(self.case_ids),
            "observation_count": self.observation_count,
            "fit_data_hash": self.fit_data_hash,
            "calibration_kind": self.calibration_kind,
        }


def _required_scale(observation: DiscrepancyCalibrationObservation) -> float:
    if observation.truth >= observation.estimate:
        denominator = observation.upper - observation.estimate
        numerator = observation.truth - observation.estimate
    else:
        denominator = observation.estimate - observation.lower
        numerator = observation.estimate - observation.truth
    if numerator <= 0.0:
        return 1.0
    if denominator <= 0.0:
        raise ValueError(
            "A zero-width development discrepancy interval misses its truth; "
            "calibration cannot infer a finite dilation."
        )
    return max(1.0, numerator / denominator)


def fit_discrepancy_calibrator(
    observations: Iterable[DiscrepancyCalibrationObservation],
    *,
    target_coverage: float = 0.95,
    phase: str = "development",
) -> DiscrepancyCalibrator:
    """Fit one case-blocked dilation factor; locked truth is rejected."""

    records = tuple(observations)
    if not records:
        raise ValueError("At least one discrepancy development observation is required.")
    if any(record.phase != phase for record in records) or phase != "development":
        raise ValueError("Discrepancy calibration requires development observations only.")
    seen: set[tuple[str, str]] = set()
    by_case: dict[str, list[float]] = {}
    for record in records:
        key = (record.case_id, record.target_id)
        if key in seen:
            raise ValueError(f"Duplicate discrepancy calibration target: {key!r}")
        seen.add(key)
        by_case.setdefault(record.case_id, []).append(_required_scale(record))
    case_requirements = sorted(max(values) for values in by_case.values())
    coverage = _finite(target_coverage, name="target_coverage")
    if not 0.0 < coverage < 1.0:
        raise ValueError("target_coverage must lie strictly between 0 and 1.")
    # Conservative finite-sample order statistic: use the smallest order
    # statistic whose empirical coverage is at least the target.
    index = min(len(case_requirements) - 1, max(0, math.ceil(coverage * len(case_requirements)) - 1))
    factor = max(1.0, case_requirements[index])
    payload = [
        record.to_dict(include_truth=True)
        for record in sorted(records, key=lambda item: (item.case_id, item.target_id))
    ]
    return DiscrepancyCalibrator(
        target_coverage=coverage,
        scale_factor=factor,
        fit_scope="development_only",
        case_ids=tuple(by_case),
        observation_count=len(records),
        fit_data_hash=_json_hash(payload),
    )


def apply_discrepancy_calibrator(
    calibrator: DiscrepancyCalibrator,
    report: DiscrepancyReport,
) -> DiscrepancyReport:
    """Apply frozen development calibration without accessing truth."""

    if not isinstance(calibrator, DiscrepancyCalibrator):
        raise TypeError("calibrator must be a DiscrepancyCalibrator.")
    if not isinstance(report, DiscrepancyReport):
        raise TypeError("report must be a DiscrepancyReport.")
    lower = report.estimate - calibrator.scale_factor * (report.estimate - report.lower)
    upper = report.estimate + calibrator.scale_factor * (report.upper - report.estimate)
    return replace(
        report,
        lower=lower,
        upper=upper,
        interval_kind="development_fitted_discrepancy_interval",
        calibration_status="CALIBRATED",
        calibration_scope="development_only",
        calibration_fit_data_hash=calibrator.fit_data_hash,
        calibration_factor=calibrator.scale_factor,
        calibration_target_coverage=calibrator.target_coverage,
    )


def apply_discrepancy_calibrator_to_record(
    calibrator: DiscrepancyCalibrator,
    record: Mapping[str, object],
) -> dict[str, object]:
    """Apply frozen calibration to an archived JSON report record."""

    if not isinstance(record, Mapping):
        raise TypeError("record must be a mapping.")
    forbidden = _truth_fields(record)
    if forbidden:
        raise ValueError(
            "Truth fields cannot be present during discrepancy calibration application: "
            + ", ".join(sorted(forbidden))
        )
    estimate = _finite(record.get("estimate"), name="record.estimate")
    lower = _finite(record.get("lower"), name="record.lower")
    upper = _finite(record.get("upper"), name="record.upper")
    if lower > estimate or estimate > upper:
        raise ValueError("record interval must contain its estimate.")
    calibrated = dict(record)
    calibrated.update(
        {
            "lower": estimate - calibrator.scale_factor * (estimate - lower),
            "upper": estimate + calibrator.scale_factor * (upper - estimate),
            "interval_kind": "development_fitted_discrepancy_interval",
            "calibration_status": "CALIBRATED",
            "calibrated": True,
            "calibration_scope": "development_only",
            "calibration_fit_data_hash": calibrator.fit_data_hash,
            "calibration_factor": calibrator.scale_factor,
            "calibration_target_coverage": calibrator.target_coverage,
        }
    )
    return calibrated


def _score_report_set(
    observations: Sequence[DiscrepancyCalibrationObservation],
    reports: Mapping[str, DiscrepancyReport | Mapping[str, object]],
    *,
    target_coverage: float,
) -> dict[str, object]:
    covered = 0
    widths: list[float] = []
    committed = 0
    false_commitments = 0
    for observation in observations:
        report = reports.get(observation.target_id)
        if report is None:
            continue
        if isinstance(report, Mapping):
            lower = _finite(report.get("lower"), name="report.lower")
            upper = _finite(report.get("upper"), name="report.upper")
            decision = report.get("decision")
            if isinstance(decision, Mapping):
                decision_kind = str(decision.get("decision_kind", ""))
            else:
                decision_kind = str(decision)
        else:
            lower = report.lower
            upper = report.upper
            decision_kind = str(report.decision.decision_kind)
        inside = lower <= observation.truth <= upper
        covered += int(inside)
        widths.append(upper - lower)
        is_committed = decision_kind.endswith(DecisionKind.SET_REPORT.value)
        committed += int(is_committed)
        false_commitments += int(is_committed and not inside)
    n = len(widths)
    coverage = covered / n if n else None
    calibration_error = None if coverage is None else abs(coverage - target_coverage)
    return {
        "status": "scored" if n == len(observations) and n else "ABSTAIN_MISSING_RECORDS",
        "n": n,
        "coverage": coverage,
        "target_coverage": target_coverage,
        "calibration_error": calibration_error,
        "mean_interval_width": sum(widths) / n if n else None,
        "n_committed": committed,
        "false_commitment_count": false_commitments,
        "false_commitment_rate": (
            false_commitments / committed if committed else None
        ),
    }


def score_locked_discrepancy(
    observations: Iterable[DiscrepancyCalibrationObservation],
    raw_reports: Mapping[str, DiscrepancyReport],
    calibrated_reports: Mapping[str, DiscrepancyReport | Mapping[str, object]],
    *,
    target_coverage: float,
) -> dict[str, object]:
    """Score locked raw vs development-calibrated reports, fail closed."""

    records = tuple(observations)
    if not records:
        return {"status": "ABSTAIN_NO_LOCKED_OBSERVATIONS"}
    if any(record.phase != "locked_test" for record in records):
        raise ValueError("Locked discrepancy scoring accepts locked_test observations only.")
    for report_set in (raw_reports, calibrated_reports):
        for report in report_set.values():
            if isinstance(report, Mapping):
                forbidden = _truth_fields(report)
                if forbidden:
                    return {
                        "status": "ABSTAIN_TRUTH_LEAKAGE",
                        "truth_fields": list(forbidden),
                    }
    if not all(
        (
            (
                isinstance(report, DiscrepancyReport)
                and report.calibration_status == "CALIBRATED"
                and report.calibration_scope == "development_only"
            )
            or (
                isinstance(report, Mapping)
                and report.get("calibration_status") == "CALIBRATED"
                and report.get("calibration_scope") == "development_only"
            )
        )
        for report in calibrated_reports.values()
    ):
        return {
            "status": "ABSTAIN_CALIBRATION_NOT_VERIFIED",
            "n": 0,
        }
    raw = _score_report_set(records, raw_reports, target_coverage=target_coverage)
    calibrated = _score_report_set(
        records, calibrated_reports, target_coverage=target_coverage
    )
    if raw["status"] != "scored" or calibrated["status"] != "scored":
        return {
            "status": "ABSTAIN_MISSING_LOCKED_REPORTS",
            "raw": raw,
            "calibrated": calibrated,
        }
    raw_error = float(raw["calibration_error"])
    calibrated_error = float(calibrated["calibration_error"])
    return {
        "status": "scored",
        "n": len(records),
        "raw": raw,
        "calibrated": calibrated,
        "calibration_degradation": calibrated_error - raw_error,
        "calibration_degradation_definition": (
            "abs(calibrated.coverage-target)-abs(raw.coverage-target)"
        ),
        "calibration_improved": calibrated_error < raw_error,
        "false_commitment_rate_available": (
            calibrated["false_commitment_rate"] is not None
        ),
    }


def build_discrepancy_report(
    target: str,
    predictions: Iterable[ScenarioPrediction],
    *,
    disagreement_threshold: float = 0.25,
    scale_floor: float = 1.0,
    compatible_hypotheses: Sequence[str] = (),
) -> DiscrepancyReport:
    """Combine scenario outputs without hiding material disagreement.

    The returned interval is the min/max envelope of scenario intervals.  It
    is intentionally labelled non-calibrated; calibration and scenario-weight
    estimation belong to a development-only fitting stage that is not part of
    this reporting primitive.
    """

    target = str(target).strip()
    if not target:
        raise ValueError("Discrepancy reports require target.")
    threshold = _finite(disagreement_threshold, name="disagreement_threshold")
    floor = _finite(scale_floor, name="scale_floor")
    if threshold < 0.0:
        raise ValueError("disagreement_threshold must be non-negative.")
    if floor <= 0.0:
        raise ValueError("scale_floor must be positive.")

    scenarios = tuple(predictions)
    if not scenarios:
        raise ValueError("At least one scenario prediction is required.")
    ids = [scenario.scenario_id for scenario in scenarios]
    if len(ids) != len(set(ids)):
        raise ValueError("Scenario IDs must be unique within one report.")
    total_weight = sum(scenario.weight for scenario in scenarios)
    if total_weight <= 0.0:
        raise ValueError("At least one scenario must have positive weight.")

    weights = {
        scenario.scenario_id: scenario.weight / total_weight
        for scenario in scenarios
    }
    estimate = sum(
        weights[scenario.scenario_id] * scenario.estimate for scenario in scenarios
    )
    lower = min(scenario.lower for scenario in scenarios)
    upper = max(scenario.upper for scenario in scenarios)
    scenario_range = max(scenario.estimate for scenario in scenarios) - min(
        scenario.estimate for scenario in scenarios
    )
    relative_disagreement = scenario_range / max(floor, abs(estimate))
    intervals_overlap = max(scenario.lower for scenario in scenarios) <= min(
        scenario.upper for scenario in scenarios
    )
    model_disagreement = (
        not intervals_overlap or relative_disagreement > threshold
    )
    disagreement_status = (
        "MODEL_DISAGREEMENT" if model_disagreement else "NO_MATERIAL_DISAGREEMENT"
    )
    scenario_families = tuple(sorted({scenario.scenario_family for scenario in scenarios}))

    hypotheses = tuple(str(item) for item in compatible_hypotheses)
    if model_disagreement:
        identifiability = IdentifiabilityStatus.MODEL_DISAGREEMENT
        decision = ProgrammeDecision(
            decision_kind=DecisionKind.ABSTAIN,
            identifiability=identifiability,
            reason=(
                f"Scenario predictions for {target!r} disagree: relative "
                f"spread={relative_disagreement:.6g}, "
                f"intervals_overlap={intervals_overlap}."
            ),
            scenario_count=len(scenarios),
            provenance={"interval_kind": "scenario_envelope_not_calibrated"},
        )
    elif hypotheses:
        identifiability = IdentifiabilityStatus.PARTIALLY_IDENTIFIED
        decision = ProgrammeDecision(
            decision_kind=DecisionKind.SET_REPORT,
            identifiability=identifiability,
            reason=(
                "Declared discrepancy scenarios are compatible within the "
                "predeclared disagreement threshold."
            ),
            compatible_hypotheses=hypotheses,
            scenario_count=len(scenarios),
            provenance={"interval_kind": "scenario_envelope_not_calibrated"},
        )
    else:
        identifiability = IdentifiabilityStatus.UNKNOWN
        decision = ProgrammeDecision(
            decision_kind=DecisionKind.ABSTAIN,
            identifiability=identifiability,
            reason=(
                "Scenario predictions are compatible, but no reportable "
                "hypothesis set was supplied."
            ),
            scenario_count=len(scenarios),
            provenance={"interval_kind": "scenario_envelope_not_calibrated"},
        )

    return DiscrepancyReport(
        target=target,
        estimate=estimate,
        lower=lower,
        upper=upper,
        scenario_range=scenario_range,
        relative_disagreement=relative_disagreement,
        intervals_overlap=intervals_overlap,
        identifiability=identifiability,
        decision=decision,
        scenario_weights=weights,
        scenarios=scenarios,
        calibration_status="NOT_CALIBRATED",
        weight_source="DECLARED_SCENARIO_WEIGHTS",
        model_disagreement_status=disagreement_status,
        identifiability_reason=decision.reason,
        scenario_families=scenario_families,
    )


__all__ = [
    "DiscrepancyCalibrationObservation",
    "DiscrepancyCalibrator",
    "DiscrepancyReport",
    "ScenarioPrediction",
    "apply_discrepancy_calibrator",
    "apply_discrepancy_calibrator_to_record",
    "build_discrepancy_report",
    "fit_discrepancy_calibrator",
    "score_locked_discrepancy",
]
