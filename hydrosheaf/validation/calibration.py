"""Development-only uncertainty calibration and selective-risk scoring.

The calibration fit is deliberately separate from inference.  A development
truth table may fit an additive split-conformal expansion for age intervals;
the resulting frozen object can then be applied to blind predictions without
seeing truth.  Locked or field truth is rejected by the fit function.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable, Mapping, Sequence


def _finite(value: object, *, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return result


@dataclass(frozen=True)
class AgeCalibrationObservation:
    """One development-only truth/prediction pair used to fit calibration."""

    target_id: str
    truth: float
    estimate: float
    lower: float
    upper: float

    def __post_init__(self) -> None:
        target_id = str(self.target_id).strip()
        if not target_id:
            raise ValueError("Calibration observations require target_id.")
        truth = _finite(self.truth, name="truth")
        estimate = _finite(self.estimate, name="estimate")
        lower = _finite(self.lower, name="lower")
        upper = _finite(self.upper, name="upper")
        if lower > upper:
            raise ValueError("Calibration lower bound cannot exceed upper bound.")
        if not lower <= estimate <= upper:
            raise ValueError("Calibration estimate must lie within its interval.")
        object.__setattr__(self, "target_id", target_id)
        object.__setattr__(self, "truth", truth)
        object.__setattr__(self, "estimate", estimate)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True)
class AgeIntervalCalibrator:
    """Frozen development-fitted interval expansion rule."""

    target_coverage: float
    additive_radius: float
    fit_count: int
    fit_scope: str = "development_only"
    calibration_kind: str = "split_conformal_interval_expansion"
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        coverage = _finite(self.target_coverage, name="target_coverage")
        radius = _finite(self.additive_radius, name="additive_radius")
        if not 0.0 < coverage < 1.0:
            raise ValueError("target_coverage must lie strictly between 0 and 1.")
        if radius < 0.0:
            raise ValueError("additive_radius must be non-negative.")
        if int(self.fit_count) < 1:
            raise ValueError("fit_count must be positive.")
        if str(self.fit_scope) != "development_only":
            raise ValueError("Calibrators must be fitted on development data only.")
        object.__setattr__(self, "target_coverage", coverage)
        object.__setattr__(self, "additive_radius", radius)
        object.__setattr__(self, "fit_count", int(self.fit_count))

    def to_dict(self) -> dict[str, object]:
        return {
            "target_coverage": self.target_coverage,
            "additive_radius": self.additive_radius,
            "fit_count": self.fit_count,
            "fit_scope": self.fit_scope,
            "calibration_kind": self.calibration_kind,
            "provenance": dict(self.provenance),
        }


def fit_age_interval_calibrator(
    observations: Iterable[AgeCalibrationObservation],
    *,
    target_coverage: float = 0.95,
    phase: str = "development",
) -> AgeIntervalCalibrator:
    """Fit a split-conformal interval expansion from development truth.

    The nonconformity score is the distance from the supplied interval to the
    truth, or zero when the interval already contains the truth.  The finite
    sample order statistic is selected conservatively at the requested
    coverage.  ``phase`` is explicit so a locked-test caller cannot silently
    fit on evaluation truth.
    """

    if str(phase) != "development":
        raise ValueError("Age interval calibration may only be fitted in development phase.")
    coverage = _finite(target_coverage, name="target_coverage")
    if not 0.0 < coverage < 1.0:
        raise ValueError("target_coverage must lie strictly between 0 and 1.")
    rows = tuple(observations)
    if not rows:
        raise ValueError("At least one development observation is required.")
    scores = sorted(
        max(row.lower - row.truth, row.truth - row.upper, 0.0) for row in rows
    )
    rank = int(math.ceil((len(scores) + 1) * coverage)) - 1
    rank = min(max(rank, 0), len(scores) - 1)
    return AgeIntervalCalibrator(
        target_coverage=coverage,
        additive_radius=float(scores[rank]),
        fit_count=len(rows),
        provenance={
            "fit_rule": "finite_sample_conformal_order_statistic",
            "truth_used_for": "development_fit_only",
        },
    )


def apply_age_interval_calibrator(
    calibrator: AgeIntervalCalibrator,
    predictions: Mapping[str, Mapping[str, object]],
) -> dict[str, dict[str, float | str]]:
    """Apply a frozen calibrator to blind age predictions."""

    calibrated: dict[str, dict[str, float | str]] = {}
    for target_id, prediction in predictions.items():
        required = ("mean_age_years", "age_95_low", "age_95_high")
        if not all(field in prediction for field in required):
            raise KeyError(f"Prediction {target_id!r} lacks an age interval field.")
        estimate = _finite(prediction["mean_age_years"], name="mean_age_years")
        lower = _finite(prediction["age_95_low"], name="age_95_low")
        upper = _finite(prediction["age_95_high"], name="age_95_high")
        if lower > upper or not lower <= estimate <= upper:
            raise ValueError(f"Prediction {target_id!r} has an invalid age interval.")
        calibrated[str(target_id)] = {
            "mean_age_years": estimate,
            "age_95_low": lower,
            "age_95_high": upper,
            "calibrated_age_low": lower - calibrator.additive_radius,
            "calibrated_age_high": upper + calibrator.additive_radius,
            "calibration_status": "development_fitted",
        }
    return calibrated


def score_calibrated_age_intervals(
    truth: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    """Score calibrated intervals after the inference/calibration boundary."""

    comparable: list[tuple[float, float, float]] = []
    for target_id, true_value in truth.items():
        prediction = predictions.get(target_id)
        if not isinstance(prediction, Mapping):
            continue
        required = ("mean_age_years", "calibrated_age_low", "calibrated_age_high")
        if not all(field in prediction for field in required):
            continue
        comparable.append(
            (
                _finite(true_value, name="truth"),
                _finite(prediction["calibrated_age_low"], name="calibrated_age_low"),
                _finite(prediction["calibrated_age_high"], name="calibrated_age_high"),
            )
        )
    if not comparable:
        return {"status": "no_comparable_outputs", "n": 0}
    coverage = sum(low <= value <= high for value, low, high in comparable) / len(comparable)
    width = sum(high - low for _, low, high in comparable) / len(comparable)
    return {
        "status": "scored",
        "n": len(comparable),
        "coverage": float(coverage),
        "mean_width": float(width),
        "all_finite": all(math.isfinite(value) for row in comparable for value in row),
    }


@dataclass(frozen=True)
class SelectiveRiskPoint:
    """Risk at one retained fraction of the prediction set."""

    requested_acceptance: float
    accepted: int
    total: int
    acceptance_rate: float
    abstention_rate: float
    mean_absolute_error: float
    interval_coverage: float
    false_commitment_rate: float

    def to_dict(self) -> dict[str, object]:
        return {
            "requested_acceptance": self.requested_acceptance,
            "accepted": self.accepted,
            "total": self.total,
            "acceptance_rate": self.acceptance_rate,
            "abstention_rate": self.abstention_rate,
            "mean_absolute_error": self.mean_absolute_error,
            "interval_coverage": self.interval_coverage,
            "false_commitment_rate": self.false_commitment_rate,
        }


def score_selective_risk(
    rows: Iterable[Mapping[str, object]],
    *,
    acceptance_fractions: Sequence[float] = (0.25, 0.50, 0.75, 1.0),
) -> list[SelectiveRiskPoint]:
    """Score risk after accepting the least-uncertain predictions first.

    Each row is a scoring-only record with ``truth``, ``estimate``, ``lower``,
    ``upper``, and non-negative ``uncertainty`` fields.  The function does not
    alter the inference decision; it reports the trade-off between abstention
    and error for a locked acceptance rule.
    """

    prepared: list[tuple[str, float, float, float, float, float]] = []
    for index, row in enumerate(rows):
        required = ("truth", "estimate", "lower", "upper", "uncertainty")
        if not all(field in row for field in required):
            raise KeyError(f"Selective-risk row {index} lacks a required field.")
        truth = _finite(row["truth"], name="truth")
        estimate = _finite(row["estimate"], name="estimate")
        lower = _finite(row["lower"], name="lower")
        upper = _finite(row["upper"], name="upper")
        uncertainty = _finite(row["uncertainty"], name="uncertainty")
        if lower > upper or uncertainty < 0.0:
            raise ValueError(f"Selective-risk row {index} has invalid bounds or uncertainty.")
        target_id = str(row.get("target_id", index))
        prepared.append((target_id, truth, estimate, lower, upper, uncertainty))
    if not prepared:
        return []
    ordered = sorted(prepared, key=lambda item: (item[5], item[0]))
    points: list[SelectiveRiskPoint] = []
    for requested in acceptance_fractions:
        fraction = _finite(requested, name="acceptance_fraction")
        if not 0.0 < fraction <= 1.0:
            raise ValueError("acceptance fractions must lie in (0, 1].")
        accepted_count = max(1, int(math.ceil(fraction * len(ordered))))
        accepted = ordered[:accepted_count]
        errors = [abs(row[1] - row[2]) for row in accepted]
        covered = [row[3] <= row[1] <= row[4] for row in accepted]
        coverage = sum(covered) / len(covered)
        points.append(
            SelectiveRiskPoint(
                requested_acceptance=fraction,
                accepted=accepted_count,
                total=len(ordered),
                acceptance_rate=accepted_count / len(ordered),
                abstention_rate=1.0 - accepted_count / len(ordered),
                mean_absolute_error=sum(errors) / len(errors),
                interval_coverage=float(coverage),
                false_commitment_rate=1.0 - float(coverage),
            )
        )
    return points


__all__ = [
    "AgeCalibrationObservation",
    "AgeIntervalCalibrator",
    "SelectiveRiskPoint",
    "apply_age_interval_calibrator",
    "fit_age_interval_calibrator",
    "score_calibrated_age_intervals",
    "score_selective_risk",
]
