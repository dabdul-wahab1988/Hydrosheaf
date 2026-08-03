"""Development-only calibration for specialist age and reaction outputs.

The specialist predictors are intentionally separated from this module.  A
predictor may be used truth-blind in a locked run; only this module receives
development truth, fits a frozen calibration object, and applies that object
to later predictions.  Case IDs are mandatory so wells from one synthetic
case cannot masquerade as independent calibration cases.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from collections.abc import Iterable, Mapping
from typing import Any

from .metrics import classification_metrics, interval_coverage, regression_metrics


ABSTAIN = "abstain"
SELECT = "select"
_AGE_MAX_YEARS = 200.0
_DEFAULT_REACTION_CLASSES = (
    "carbonate",
    "silicate_exchange",
    "sulfate_reduction",
    "iron_reduction",
    "denitrification",
    "sulfate_source",
    "other_redox",
    "other",
    "none",
)


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, got {value!r}.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return number


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(float(lower), min(float(upper), float(value)))


def _median(values: Iterable[float]) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("Median requires at least one value.")
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def _jsonable(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True, default=str))


def _fingerprint(record: Mapping[str, Any]) -> str:
    payload = json.dumps(record, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SpecialistAgeCalibrationObservation:
    """One age prediction/truth pair from one development case."""

    case_id: str
    target_id: str
    truth: float
    estimate: float
    lower: float
    upper: float
    selection_score: float | None = None
    selected: bool = True

    def __post_init__(self) -> None:
        case_id = str(self.case_id).strip()
        target_id = str(self.target_id).strip()
        if not case_id or not target_id:
            raise ValueError("Age calibration observations require case_id and target_id.")
        truth = _finite(self.truth, name="truth")
        estimate = _finite(self.estimate, name="estimate")
        lower = _finite(self.lower, name="lower")
        upper = _finite(self.upper, name="upper")
        if lower > upper or not lower <= estimate <= upper:
            raise ValueError("Age calibration interval is invalid.")
        selection_score = (
            None
            if self.selection_score is None
            else _finite(self.selection_score, name="selection_score")
        )
        if selection_score is not None and selection_score < 0.0:
            raise ValueError("selection_score must be non-negative.")
        object.__setattr__(self, "case_id", case_id)
        object.__setattr__(self, "target_id", target_id)
        object.__setattr__(self, "truth", truth)
        object.__setattr__(self, "estimate", estimate)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        object.__setattr__(self, "selection_score", selection_score)
        object.__setattr__(self, "selected", bool(self.selected))


@dataclass(frozen=True)
class SpecialistAgeCalibrator:
    """Frozen case-blocked bias and interval calibrator for age outputs."""

    target_coverage: float
    additive_radius: float
    bias_offset_years: float
    fit_count: int
    case_ids: tuple[str, ...]
    fit_scope: str = "development_only"
    calibration_kind: str = "case_blocked_median_bias_split_conformal"
    max_age_years: float = _AGE_MAX_YEARS
    provenance: Mapping[str, object] = field(default_factory=dict)
    interval_scale: float = 1.0
    selection_threshold: float | None = None
    selection_all: bool = False
    selection_score_field: str = "tracer_age_spread_years"
    minimum_selection_rate: float = 0.0

    def __post_init__(self) -> None:
        coverage = _finite(self.target_coverage, name="target_coverage")
        radius = _finite(self.additive_radius, name="additive_radius")
        bias = _finite(self.bias_offset_years, name="bias_offset_years")
        max_age = _finite(self.max_age_years, name="max_age_years")
        if not 0.0 < coverage < 1.0:
            raise ValueError("target_coverage must lie strictly between 0 and 1.")
        if radius < 0.0 or max_age <= 0.0:
            raise ValueError("Age calibrator bounds are invalid.")
        interval_scale = _finite(self.interval_scale, name="interval_scale")
        minimum_selection_rate = _finite(
            self.minimum_selection_rate,
            name="minimum_selection_rate",
        )
        selection_threshold = (
            None
            if self.selection_threshold is None
            else _finite(self.selection_threshold, name="selection_threshold")
        )
        if (
            interval_scale < 0.0
            or not 0.0 <= minimum_selection_rate <= 1.0
            or selection_threshold is not None and selection_threshold < 0.0
            or bool(self.selection_all) and selection_threshold is not None
        ):
            raise ValueError("Age interval/selection calibration parameters are invalid.")
        if int(self.fit_count) < 1:
            raise ValueError("Age calibration requires at least one observation.")
        if str(self.fit_scope) != "development_only":
            raise ValueError("Age calibration must be fitted on development data only.")
        case_ids = tuple(sorted({str(case_id) for case_id in self.case_ids}))
        if not case_ids:
            raise ValueError("Age calibration requires at least one case ID.")
        object.__setattr__(self, "target_coverage", coverage)
        object.__setattr__(self, "additive_radius", radius)
        object.__setattr__(self, "bias_offset_years", bias)
        object.__setattr__(self, "max_age_years", max_age)
        object.__setattr__(self, "fit_count", int(self.fit_count))
        object.__setattr__(self, "case_ids", case_ids)
        object.__setattr__(self, "interval_scale", interval_scale)
        object.__setattr__(self, "selection_threshold", selection_threshold)
        object.__setattr__(self, "selection_all", bool(self.selection_all))
        object.__setattr__(self, "selection_score_field", str(self.selection_score_field))
        object.__setattr__(self, "minimum_selection_rate", minimum_selection_rate)

    def to_dict(self) -> dict[str, object]:
        record: dict[str, object] = {
            "target_coverage": self.target_coverage,
            "additive_radius": self.additive_radius,
            "bias_offset_years": self.bias_offset_years,
            "fit_count": self.fit_count,
            "case_ids": list(self.case_ids),
            "fit_scope": self.fit_scope,
            "calibration_kind": self.calibration_kind,
            "max_age_years": self.max_age_years,
            "provenance": _jsonable(self.provenance),
            "interval_scale": self.interval_scale,
            "selection_threshold": self.selection_threshold,
            "selection_all": self.selection_all,
            "selection_score_field": self.selection_score_field,
            "minimum_selection_rate": self.minimum_selection_rate,
        }
        record["fingerprint"] = _fingerprint(record)
        return record

    def apply(
        self,
        predictions: Mapping[str, Mapping[str, object]],
    ) -> dict[str, dict[str, object]]:
        """Apply the frozen age calibration without truth."""

        calibrated: dict[str, dict[str, object]] = {}
        for target_id, raw in predictions.items():
            if not isinstance(raw, Mapping):
                raise TypeError(f"Age prediction {target_id!r} must be a mapping.")
            output = dict(raw)
            estimate = _finite(raw.get("mean_age_years"), name="mean_age_years")
            lower = _finite(raw.get("age_95_low"), name="age_95_low")
            upper = _finite(raw.get("age_95_high"), name="age_95_high")
            raw_decision = str(raw.get("decision", SELECT))
            raw_score = raw.get(self.selection_score_field)
            selection_score = (
                None
                if raw_score is None
                else _finite(raw_score, name=self.selection_score_field)
            )
            selected = raw_decision != ABSTAIN
            if self.selection_all:
                selected = True
            elif self.selection_threshold is not None and selection_score is not None:
                selected = selection_score <= self.selection_threshold
            if None in (estimate, lower, upper):
                output["calibration_status"] = (
                    "abstained_not_calibrated" if not selected else "invalid_not_calibrated"
                )
                calibrated[str(target_id)] = output
                continue
            corrected = _clamp(
                estimate + self.bias_offset_years,
                0.0,
                self.max_age_years,
            )
            shifted_lower = _clamp(lower + self.bias_offset_years, 0.0, self.max_age_years)
            shifted_upper = _clamp(upper + self.bias_offset_years, 0.0, self.max_age_years)
            calibrated_lower = _clamp(
                corrected
                - self.interval_scale * max(corrected - shifted_lower, 0.0)
                - self.additive_radius,
                0.0,
                self.max_age_years,
            )
            calibrated_upper = _clamp(
                corrected
                + self.interval_scale * max(shifted_upper - corrected, 0.0)
                + self.additive_radius,
                0.0,
                self.max_age_years,
            )
            calibrated_lower = min(calibrated_lower, corrected)
            calibrated_upper = max(calibrated_upper, corrected)
            output.update(
                {
                    "calibrated_mean_age_years": corrected,
                    "calibrated_age_low": calibrated_lower,
                    "calibrated_age_high": calibrated_upper,
                    "decision": SELECT if selected else ABSTAIN,
                    "calibration_status": (
                        "development_fitted" if selected else "development_fitted_abstain"
                    ),
                    "selection_calibration_status": "development_fitted",
                    "selection_score": selection_score,
                    "selection_threshold": self.selection_threshold,
                    "interval_scale": self.interval_scale,
                }
            )
            calibrated[str(target_id)] = output
        return calibrated


def fit_specialist_age_calibrator(
    observations: Iterable[SpecialistAgeCalibrationObservation],
    *,
    target_coverage: float = 0.95,
    phase: str = "development",
    max_age_years: float = _AGE_MAX_YEARS,
    interval_scales: Iterable[float] = (0.0, 0.25, 0.50, 0.75, 1.0),
    minimum_selection_rate: float = 0.80,
) -> SpecialistAgeCalibrator:
    """Fit a development-only interval-scale and selection calibrator.

    Candidate scales and selection thresholds are chosen using development
    truth only.  The objective is the narrowest selected interval satisfying
    the declared coverage target and minimum selection rate.  If no candidate
    satisfies both constraints, the legacy unscaled interval and raw
    selection decisions are retained.
    """

    if str(phase) != "development":
        raise ValueError("Specialist age calibration may only fit development data.")
    coverage = _finite(target_coverage, name="target_coverage")
    if not 0.0 < coverage < 1.0:
        raise ValueError("target_coverage must lie strictly between 0 and 1.")
    rows = tuple(observations)
    if not rows:
        raise ValueError("At least one specialist age observation is required.")
    minimum_rate = _finite(minimum_selection_rate, name="minimum_selection_rate")
    if not 0.0 < minimum_rate <= 1.0:
        raise ValueError("minimum_selection_rate must lie in (0, 1].")
    scales = tuple(
        sorted({_finite(value, name="interval_scale") for value in interval_scales})
    )
    if not scales or any(value < 0.0 for value in scales):
        raise ValueError("interval_scales must contain finite non-negative values.")
    by_case: dict[str, list[SpecialistAgeCalibrationObservation]] = {}
    for row in rows:
        by_case.setdefault(row.case_id, []).append(row)
    case_biases = {
        case_id: _median(row.truth - row.estimate for row in case_rows)
        for case_id, case_rows in by_case.items()
    }
    bias_offset = _median(case_biases.values())

    def bounds(
        row: SpecialistAgeCalibrationObservation,
        scale: float,
        radius: float,
    ) -> tuple[float, float]:
        center = _clamp(row.estimate + bias_offset, 0.0, max_age_years)
        shifted_lower = _clamp(row.lower + bias_offset, 0.0, max_age_years)
        shifted_upper = _clamp(row.upper + bias_offset, 0.0, max_age_years)
        lower = _clamp(
            center - scale * max(center - shifted_lower, 0.0) - radius,
            0.0,
            max_age_years,
        )
        upper = _clamp(
            center + scale * max(shifted_upper - center, 0.0) + radius,
            0.0,
            max_age_years,
        )
        return min(lower, center), max(upper, center)

    def conformal_radius(scale: float) -> float:
        case_scores = {
            case_id: max(
                max(
                    bounds(row, scale, 0.0)[0] - row.truth,
                    row.truth - bounds(row, scale, 0.0)[1],
                    0.0,
                )
                for row in case_rows
            )
            for case_id, case_rows in by_case.items()
        }
        ordered_scores = sorted(case_scores.values())
        rank = int(math.ceil((len(ordered_scores) + 1) * coverage)) - 1
        rank = min(max(rank, 0), len(ordered_scores) - 1)
        return float(ordered_scores[rank])

    score_values = sorted(
        {
            float(row.selection_score)
            for row in rows
            if row.selection_score is not None
        }
    )
    threshold_candidates: tuple[float | None, ...] = (
        tuple(score_values) + (float("inf"),) if score_values else (None,)
    )

    def selected(row: SpecialistAgeCalibrationObservation, threshold: float | None) -> bool:
        if threshold is not None and math.isinf(threshold):
            return True
        if threshold is None or row.selection_score is None:
            return bool(row.selected)
        return float(row.selection_score) <= threshold

    feasible: list[tuple[float, float, float | None, float, float, float, float]] = []
    for scale in scales:
        radius = conformal_radius(scale)
        for threshold in threshold_candidates:
            selected_rows = [row for row in rows if selected(row, threshold)]
            selection_rate = len(selected_rows) / len(rows)
            if not selected_rows or selection_rate < minimum_rate:
                continue
            covered_selected = [
                bounds(row, scale, radius)[0] <= row.truth <= bounds(row, scale, radius)[1]
                for row in selected_rows
            ]
            covered_all = sum(
                int(
                    selected(row, threshold)
                    and bounds(row, scale, radius)[0] <= row.truth <= bounds(row, scale, radius)[1]
                )
                for row in rows
            )
            conditional_coverage = sum(covered_selected) / len(covered_selected)
            inclusive_coverage = covered_all / len(rows)
            if conditional_coverage < coverage or inclusive_coverage < coverage:
                continue
            mean_width = sum(
                bounds(row, scale, radius)[1] - bounds(row, scale, radius)[0]
                for row in selected_rows
            ) / len(selected_rows)
            feasible.append(
                (
                    mean_width,
                    -selection_rate,
                    threshold,
                    scale,
                    radius,
                    conditional_coverage,
                    inclusive_coverage,
                )
            )

    if feasible:
        (
            mean_width,
            _negative_selection_rate,
            selection_threshold,
            interval_scale,
            radius,
            dev_conditional,
            dev_inclusive,
        ) = min(
            feasible,
            key=lambda item: (
                item[0],
                item[1],
                item[3],
                float("inf") if item[2] is None else item[2],
            ),
        )
        fit_status = "grid_candidate_selected"
        selected_count = sum(1 for row in rows if selected(row, selection_threshold))
    else:
        interval_scale = 1.0
        radius = conformal_radius(interval_scale)
        selection_threshold = None
        selected_count = sum(1 for row in rows if row.selected)
        dev_conditional = float("nan")
        dev_inclusive = float("nan")
        mean_width = float("nan")
        fit_status = "no_feasible_grid_candidate_legacy_fallback"
    selection_all = bool(selection_threshold is not None and math.isinf(selection_threshold))
    return SpecialistAgeCalibrator(
        target_coverage=coverage,
        additive_radius=float(radius),
        bias_offset_years=float(bias_offset),
        fit_count=len(rows),
        case_ids=tuple(sorted(by_case)),
        calibration_kind="case_blocked_bias_scaled_split_conformal_selection_grid",
        max_age_years=float(max_age_years),
        interval_scale=float(interval_scale),
        selection_threshold=(
            None
            if selection_threshold is None or selection_all
            else float(selection_threshold)
        ),
        selection_all=selection_all,
        minimum_selection_rate=float(minimum_rate),
        provenance={
            "fit_rule": (
                "case_equal_weight_bias_then_within_case_max_split_conformal_"
                "interval_scale_radius_selection_grid"
            ),
            "truth_used_for": "development_fit_only",
            "case_count": len(by_case),
            "interval_scales": list(scales),
            "selection_threshold_candidates": len(threshold_candidates),
            "selection_all": selection_all,
            "minimum_selection_rate": minimum_rate,
            "development_fit_status": fit_status,
            "development_selected_count": selected_count,
            "development_conditional_coverage": dev_conditional,
            "development_coverage_including_abstention": dev_inclusive,
            "development_mean_selected_interval_width": mean_width,
        },
    )


def score_calibrated_specialist_age(
    truth: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    """Score calibrated specialist age estimates after truth release."""

    observed: list[float] = []
    estimated: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    n_abstain = 0
    n_missing = 0
    for target_id, value in truth.items():
        row = predictions.get(str(target_id))
        if not isinstance(row, Mapping):
            n_missing += 1
            continue
        if str(row.get("decision", SELECT)) == ABSTAIN:
            n_abstain += 1
            continue
        required = (
            "calibrated_mean_age_years",
            "calibrated_age_low",
            "calibrated_age_high",
        )
        if not all(field_name in row for field_name in required):
            n_missing += 1
            continue
        observed.append(_finite(value, name="truth"))
        estimated.append(_finite(row[required[0]], name=required[0]))
        lower.append(_finite(row[required[1]], name=required[1]))
        upper.append(_finite(row[required[2]], name=required[2]))
    if not observed:
        return {
            "status": "no_comparable_outputs",
            "n": 0,
            "n_abstain": n_abstain,
            "n_missing": n_missing,
            "n_truth": len(truth),
        }
    interval = dict(interval_coverage(observed, lower, upper))
    covered = sum(
        int(low <= truth_value <= high)
        for truth_value, low, high in zip(observed, lower, upper)
    )
    return {
        "status": "scored",
        "n": len(observed),
        "n_abstain": n_abstain,
        "n_missing": n_missing,
        "n_truth": len(truth),
        "point": dict(regression_metrics(observed, estimated)),
        "interval": interval,
        "coverage_including_abstention": covered / len(truth) if truth else float("nan"),
    }


@dataclass(frozen=True)
class SpecialistReactionCalibrationObservation:
    """One reaction probability vector from one development case."""

    case_id: str
    edge_id: str
    truth_family: str
    logits: Mapping[str, float]

    def __post_init__(self) -> None:
        case_id = str(self.case_id).strip()
        edge_id = str(self.edge_id).strip()
        truth_family = str(self.truth_family).strip()
        if not case_id or not edge_id or not truth_family:
            raise ValueError("Reaction calibration identity fields are required.")
        values = {
            str(key): _finite(value, name=f"logit[{key}]")
            for key, value in self.logits.items()
        }
        if not values:
            raise ValueError("Reaction calibration requires non-empty logits.")
        object.__setattr__(self, "case_id", case_id)
        object.__setattr__(self, "edge_id", edge_id)
        object.__setattr__(self, "truth_family", truth_family)
        object.__setattr__(self, "logits", values)


@dataclass(frozen=True)
class SpecialistReactionCalibrator:
    """Frozen temperature-scaling calibrator for reaction probabilities."""

    temperature: float
    classes: tuple[str, ...]
    decision_threshold: float
    fit_count: int
    case_ids: tuple[str, ...]
    fit_scope: str = "development_only"
    calibration_kind: str = "case_blocked_temperature_scaling"
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        temperature = _finite(self.temperature, name="temperature")
        threshold = _finite(self.decision_threshold, name="decision_threshold")
        if temperature <= 0.0 or not 0.0 < threshold < 1.0:
            raise ValueError("Reaction calibration parameters are invalid.")
        classes = tuple(dict.fromkeys(str(item) for item in self.classes))
        if not classes:
            raise ValueError("Reaction calibration requires classes.")
        if int(self.fit_count) < 1 or not self.case_ids:
            raise ValueError("Reaction calibration requires observations and cases.")
        if str(self.fit_scope) != "development_only":
            raise ValueError("Reaction calibration must be fitted on development data only.")
        object.__setattr__(self, "temperature", temperature)
        object.__setattr__(self, "decision_threshold", threshold)
        object.__setattr__(self, "classes", classes)
        object.__setattr__(self, "fit_count", int(self.fit_count))
        object.__setattr__(self, "case_ids", tuple(sorted(set(self.case_ids))))

    def to_dict(self) -> dict[str, object]:
        record: dict[str, object] = {
            "temperature": self.temperature,
            "classes": list(self.classes),
            "decision_threshold": self.decision_threshold,
            "fit_count": self.fit_count,
            "case_ids": list(self.case_ids),
            "fit_scope": self.fit_scope,
            "calibration_kind": self.calibration_kind,
            "provenance": _jsonable(self.provenance),
        }
        record["fingerprint"] = _fingerprint(record)
        return record

    def apply(
        self,
        predictions: Mapping[str, Mapping[str, object]],
    ) -> dict[str, dict[str, object]]:
        """Apply temperature scaling to raw reaction logits without truth."""

        output: dict[str, dict[str, object]] = {}
        for edge_id, raw in predictions.items():
            if not isinstance(raw, Mapping):
                raise TypeError(f"Reaction prediction {edge_id!r} must be a mapping.")
            logits = raw.get("logits", raw.get("raw_scores"))
            if not isinstance(logits, Mapping):
                probabilities = raw.get("probabilities", {})
                logits = {
                    str(key): math.log(max(_finite(value, name="probability"), 1.0e-12))
                    for key, value in probabilities.items()
                } if isinstance(probabilities, Mapping) else {}
            if not logits:
                updated = dict(raw)
                updated.update(
                    {
                        "decision": ABSTAIN,
                        "calibration_status": "abstained_not_calibrated",
                        "calibration_reason": "missing_probability_or_logit_vector",
                    }
                )
                output[str(edge_id)] = updated
                continue
            scaled = {str(key): _finite(value, name=f"logit[{key}]") / self.temperature for key, value in logits.items()}
            probabilities = _softmax(scaled, self.classes)
            family = max(probabilities, key=probabilities.get)
            maximum = probabilities[family]
            updated = dict(raw)
            updated.update(
                {
                    "family": family,
                    "probabilities": probabilities,
                    "probability": maximum,
                    "decision": SELECT if maximum >= self.decision_threshold else ABSTAIN,
                    "calibration_status": "development_fitted",
                    "calibration_temperature": self.temperature,
                }
            )
            if bool(self.provenance.get("grid_boundary_hit")):
                updated["calibration_warning"] = "temperature_grid_boundary_hit"
            output[str(edge_id)] = updated
        return output


def fit_specialist_reaction_calibrator(
    observations: Iterable[SpecialistReactionCalibrationObservation],
    *,
    decision_threshold: float = 0.60,
    phase: str = "development",
) -> SpecialistReactionCalibrator:
    """Fit one temperature by case-blocked mean multiclass log score."""

    if str(phase) != "development":
        raise ValueError("Specialist reaction calibration may only fit development data.")
    threshold = _finite(decision_threshold, name="decision_threshold")
    if not 0.0 < threshold < 1.0:
        raise ValueError("decision_threshold must lie in (0, 1).")
    rows = tuple(observations)
    if not rows:
        raise ValueError("At least one specialist reaction observation is required.")
    classes = tuple(
        dict.fromkeys(
            [*_DEFAULT_REACTION_CLASSES]
            + [str(key) for row in rows for key in row.logits]
            + [row.truth_family for row in rows]
        )
    )
    by_case: dict[str, list[SpecialistReactionCalibrationObservation]] = {}
    for row in rows:
        by_case.setdefault(row.case_id, []).append(row)

    def loss(temperature: float) -> float:
        case_losses: list[float] = []
        for case_rows in by_case.values():
            values: list[float] = []
            for row in case_rows:
                probabilities = _softmax(row.logits, classes, temperature=temperature)
                values.append(-math.log(max(probabilities.get(row.truth_family, 0.0), 1.0e-12)))
            case_losses.append(sum(values) / len(values))
        return sum(case_losses) / len(case_losses)

    candidate_temperatures = tuple(
        math.exp(-2.5 + index * (5.0 / 160.0)) for index in range(161)
    )
    temperature = min(candidate_temperatures, key=lambda value: (loss(value), value))
    grid_boundary_hit = temperature in {
        candidate_temperatures[0],
        candidate_temperatures[-1],
    }
    return SpecialistReactionCalibrator(
        temperature=float(temperature),
        classes=classes,
        decision_threshold=threshold,
        fit_count=len(rows),
        case_ids=tuple(sorted(by_case)),
        provenance={
            "fit_rule": "case_equal_weight_multiclass_log_score_grid",
            "truth_used_for": "development_fit_only",
            "temperature_grid": [float(candidate_temperatures[0]), float(candidate_temperatures[-1]), len(candidate_temperatures)],
            "grid_boundary_hit": grid_boundary_hit,
        },
    )


def score_calibrated_specialist_reaction(
    truth_by_edge: Mapping[str, str],
    predictions: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    """Score calibrated multiclass reaction probabilities including ``none``."""

    expected: list[str] = []
    predicted: list[str] = []
    max_probabilities: list[float] = []
    correctness: list[float] = []
    brier_terms: list[float] = []
    log_losses: list[float] = []
    n_abstain = 0
    n_missing = 0
    for edge_id, truth_family in truth_by_edge.items():
        row = predictions.get(str(edge_id))
        if not isinstance(row, Mapping):
            n_missing += 1
            continue
        probabilities = row.get("probabilities", {})
        if not isinstance(probabilities, Mapping) or not probabilities:
            n_missing += 1
            continue
        normalised = {
            str(key): max(0.0, _finite(value, name="probability"))
            for key, value in probabilities.items()
        }
        total = sum(normalised.values())
        if total <= 0.0:
            continue
        normalised = {key: value / total for key, value in normalised.items()}
        expected.append(str(truth_family))
        family = max(normalised, key=normalised.get)
        predicted.append(family)
        maximum = normalised[family]
        max_probabilities.append(maximum)
        correctness.append(float(family == str(truth_family)))
        brier_terms.append(
            sum((probability - float(label == str(truth_family))) ** 2 for label, probability in normalised.items())
            + (0.0 if str(truth_family) in normalised else 1.0)
        )
        log_losses.append(-math.log(max(normalised.get(str(truth_family), 0.0), 1.0e-12)))
        if str(row.get("decision")) == ABSTAIN:
            n_abstain += 1
    if not expected:
        return {
            "status": "no_comparable_outputs",
            "n": 0,
            "n_abstain": n_abstain,
            "n_missing": n_missing,
            "n_truth": len(truth_by_edge),
        }
    return {
        "status": "scored",
        "n": len(expected),
        "n_abstain": n_abstain,
        "n_missing": n_missing,
        "n_truth": len(truth_by_edge),
        "coverage": (len(expected) - n_abstain) / len(expected),
        "coverage_including_abstention": (
            (len(expected) - n_abstain) / len(truth_by_edge)
            if truth_by_edge
            else float("nan")
        ),
        "metrics": dict(classification_metrics(expected, predicted)),
        "multiclass_log_loss": sum(log_losses) / len(log_losses),
        "multiclass_brier": sum(brier_terms) / len(brier_terms),
        "top_probability_ece": _ece(max_probabilities, correctness),
        "expected_families": sorted(set(expected)),
        "predicted_families": sorted(set(predicted)),
    }


def _softmax(
    logits: Mapping[str, float],
    classes: Iterable[str],
    *,
    temperature: float = 1.0,
) -> dict[str, float]:
    selected = {str(key): float(value) for key, value in logits.items()}
    for key in classes:
        selected.setdefault(str(key), 0.0)
    scaled = {key: value / float(temperature) for key, value in selected.items()}
    maximum = max(scaled.values())
    exponentials = {key: math.exp(value - maximum) for key, value in scaled.items()}
    total = sum(exponentials.values())
    return {key: value / total for key, value in exponentials.items()}


def _ece(probabilities: list[float], correctness: list[float], n_bins: int = 10) -> float:
    if not probabilities:
        return float("nan")
    error = 0.0
    n = len(probabilities)
    for index in range(n_bins):
        lower = index / n_bins
        upper = (index + 1) / n_bins
        members = [
            (probability, correct)
            for probability, correct in zip(probabilities, correctness)
            if lower <= probability < upper or (index == n_bins - 1 and probability <= upper)
        ]
        if members:
            mean_probability = sum(item[0] for item in members) / len(members)
            observed = sum(item[1] for item in members) / len(members)
            error += len(members) / n * abs(mean_probability - observed)
    return float(error)


__all__ = [
    "SpecialistAgeCalibrationObservation",
    "SpecialistAgeCalibrator",
    "fit_specialist_age_calibrator",
    "score_calibrated_specialist_age",
    "SpecialistReactionCalibrationObservation",
    "SpecialistReactionCalibrator",
    "fit_specialist_reaction_calibrator",
    "score_calibrated_specialist_reaction",
]
