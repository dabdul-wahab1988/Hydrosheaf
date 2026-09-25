"""Conservative time-history-aware TTD inversion for stable-water isotopes.

This module deliberately adds a narrow workflow instead of changing the existing
scalar tracer kernels.  Stable-water-isotope rows are built by evaluating a dated
``InputHistory`` at ``sample_year - age``.  Radioactive and atmospheric tracers
continue to use the existing ``tracer_response_kernel`` implementation.

The returned estimate is conditional on the supplied source histories.  No
fractionation, evaporation, or source-history reconstruction is performed here,
and supplied histories are never extrapolated for stable-isotope rows.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
import math
from typing import Any

import numpy as np

from .input_history import InputHistory
from .joint_lpm import tracer_response_kernel
from .tracer_inputs import normalize_tracer_key
from .ttd_grid import TTDGrid, compute_trapezoidal_weights
from .ttd_kernel_builder import MultiTracerForwardSystem
from .ttd_network_solver import solve_single_node_ttd


# Public, stable reason-code strings for callers that need to branch on an
# abstention without parsing human-readable diagnostics.
REASON_MISSING_STABLE_ISOTOPE_SOURCE_HISTORY = "MISSING_STABLE_ISOTOPE_SOURCE_HISTORY"
REASON_SOURCE_HISTORY_OUT_OF_WINDOW = "SOURCE_HISTORY_OUT_OF_WINDOW"
REASON_MALFORMED_SOURCE_HISTORY = "MALFORMED_SOURCE_HISTORY"
REASON_NONFINITE_SOURCE_HISTORY = "NONFINITE_SOURCE_HISTORY"
REASON_INVALID_STABLE_ISOTOPE_TRANSFORM = "INVALID_STABLE_ISOTOPE_TRANSFORM"
REASON_INVALID_TRACER_WEIGHT = "INVALID_TRACER_WEIGHT"
REASON_INVALID_AGE_GRID = "INVALID_AGE_GRID"
REASON_INVALID_SAMPLE_YEAR = "INVALID_SAMPLE_YEAR"
REASON_INVALID_OBSERVATION = "INVALID_OBSERVATION"
REASON_NO_OBSERVATIONS = "NO_OBSERVATIONS"
REASON_UNSUPPORTED_TRACER = "UNSUPPORTED_TRACER"
REASON_NONFINITE_RESPONSE_MATRIX = "NONFINITE_RESPONSE_MATRIX"
REASON_ZERO_EFFECTIVE_WEIGHT = "ZERO_EFFECTIVE_WEIGHT"

STABLE_ISOTOPE_TRACERS = frozenset({"d18O", "d2H"})
RADIOACTIVE_TRACERS = frozenset({"3H", "3H/3He", "14C", "39Ar", "85Kr"})
HISTORY_DRIVEN_TRACERS = frozenset({"3H", "3H/3He", "SF6", "CFC11", "CFC12", "CFC113", "85Kr"})


class _HistoryInputError(ValueError):
    """Internal validation error carrying a machine-readable reason code."""

    def __init__(self, reason: str, message: str) -> None:
        self.reason = str(reason)
        self.message = str(message)
        super().__init__(f"{self.reason}: {self.message}")


def _canonical_tracer_key(tracer: Any) -> str:
    """Canonicalize stable-isotope aliases before using HydroSheaf aliases."""

    raw = str(tracer).strip()
    if not raw:
        return ""

    compact = (
        raw.replace("δ", "DELTA")
        .replace("Δ", "DELTA")
        .replace(" ", "")
        .replace("/", "")
        .replace("-", "")
        .replace("_", "")
        .upper()
    )
    if compact in {"D18O", "18O", "DELTA18O"}:
        return "d18O"
    if compact in {"D2H", "2H", "DELTA2H"}:
        return "d2H"
    if compact in {"3H3HE", "3H3HELIUM"}:
        return "3H/3He"
    return str(normalize_tracer_key(raw)).strip()


def _coerce_observations(
    observations: Sequence[HistoryTracerObservation] | HistoryTracerObservation,
) -> tuple[HistoryTracerObservation, ...]:
    if isinstance(observations, HistoryTracerObservation):
        items = (observations,)
    else:
        try:
            items = tuple(observations)
        except TypeError as exc:
            raise _HistoryInputError(REASON_INVALID_OBSERVATION, "observations must be a sequence.") from exc

    if not items:
        raise _HistoryInputError(REASON_NO_OBSERVATIONS, "At least one tracer observation is required.")

    coerced: list[HistoryTracerObservation] = []
    seen: set[str] = set()
    for item in items:
        try:
            observation = item if isinstance(item, HistoryTracerObservation) else HistoryTracerObservation(**item)
        except (TypeError, ValueError) as exc:
            raise _HistoryInputError(REASON_INVALID_OBSERVATION, str(exc)) from exc
        if observation.tracer in seen:
            raise _HistoryInputError(
                REASON_INVALID_OBSERVATION,
                f"Tracer {observation.tracer!r} appears more than once.",
            )
        seen.add(observation.tracer)
        coerced.append(observation)
    return tuple(coerced)


def _validate_age_grid(age_grid_years: Sequence[float]) -> np.ndarray:
    try:
        ages = np.asarray(age_grid_years, dtype=float)
    except (TypeError, ValueError) as exc:
        raise _HistoryInputError(REASON_INVALID_AGE_GRID, "age_grid_years must be numeric.") from exc
    if ages.ndim != 1 or ages.size < 3:
        raise _HistoryInputError(
            REASON_INVALID_AGE_GRID,
            "age_grid_years must be a one-dimensional array with at least three points.",
        )
    if not np.all(np.isfinite(ages)) or np.any(ages < 0.0):
        raise _HistoryInputError(
            REASON_INVALID_AGE_GRID,
            "age_grid_years must contain finite, non-negative ages.",
        )
    if np.any(np.diff(ages) <= 0.0):
        raise _HistoryInputError(
            REASON_INVALID_AGE_GRID,
            "age_grid_years must be strictly increasing.",
        )
    return ages


def _validate_sample_year(sample_year: float) -> float:
    try:
        year = float(sample_year)
    except (TypeError, ValueError) as exc:
        raise _HistoryInputError(REASON_INVALID_SAMPLE_YEAR, "sample_year must be numeric.") from exc
    if not math.isfinite(year):
        raise _HistoryInputError(REASON_INVALID_SAMPLE_YEAR, "sample_year must be finite.")
    return year


def _canonicalize_histories(
    source_histories: Mapping[str, InputHistory] | None,
) -> dict[str, InputHistory]:
    if source_histories is None:
        return {}
    if not isinstance(source_histories, Mapping):
        raise _HistoryInputError(REASON_MALFORMED_SOURCE_HISTORY, "source_histories must be a mapping.")

    histories: dict[str, InputHistory] = {}
    for raw_key, history in source_histories.items():
        key = _canonical_tracer_key(raw_key)
        if not key:
            raise _HistoryInputError(REASON_MALFORMED_SOURCE_HISTORY, "source_histories contains an empty tracer key.")
        if key in histories and histories[key] is not history:
            raise _HistoryInputError(
                REASON_MALFORMED_SOURCE_HISTORY,
                f"source_histories contains duplicate aliases for {key!r}.",
            )
        histories[key] = history

    # The existing scalar kernel looks up the tritium history under ``3H`` for
    # both 3H and 3H/3He.  Preserve that behavior for the new wrapper too.
    if "3H" not in histories and "3H/3He" in histories:
        histories["3H"] = histories["3H/3He"]
    return histories


def _history_arrays(history: InputHistory, tracer: str) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(history, InputHistory):
        raise _HistoryInputError(
            REASON_MALFORMED_SOURCE_HISTORY,
            f"Source history for {tracer!r} must be an InputHistory instance.",
        )
    try:
        years = np.asarray(history.years, dtype=float)
        values = np.asarray(history.values, dtype=float)
        sigma = np.asarray(history.sigma, dtype=float)
    except (AttributeError, TypeError, ValueError) as exc:
        raise _HistoryInputError(
            REASON_MALFORMED_SOURCE_HISTORY,
            f"Source history for {tracer!r} does not expose numeric years, values, and sigma arrays.",
        ) from exc

    if years.ndim != 1 or values.ndim != 1 or sigma.ndim != 1 or years.shape != values.shape:
        raise _HistoryInputError(
            REASON_MALFORMED_SOURCE_HISTORY,
            f"Source history for {tracer!r} has inconsistent array shapes.",
        )
    if years.size < 2:
        raise _HistoryInputError(
            REASON_MALFORMED_SOURCE_HISTORY,
            f"Source history for {tracer!r} must contain at least two dated values.",
        )
    if sigma.shape != values.shape:
        raise _HistoryInputError(
            REASON_MALFORMED_SOURCE_HISTORY,
            f"Source-history uncertainty for {tracer!r} has an inconsistent shape.",
        )
    if not np.all(np.isfinite(years)) or not np.all(np.isfinite(values)):
        raise _HistoryInputError(
            REASON_NONFINITE_SOURCE_HISTORY,
            f"Source history for {tracer!r} contains non-finite years or values.",
        )
    if not np.all(np.isfinite(sigma)) or np.any(sigma < 0.0):
        raise _HistoryInputError(
            REASON_NONFINITE_SOURCE_HISTORY,
            f"Source-history uncertainty for {tracer!r} is non-finite or negative.",
        )
    if np.any(np.diff(years) <= 0.0):
        raise _HistoryInputError(
            REASON_MALFORMED_SOURCE_HISTORY,
            f"Source history for {tracer!r} must have strictly increasing years.",
        )
    return years, values


def _validate_history_window(
    history: InputHistory,
    tracer: str,
    recharge_years: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    years, values = _history_arrays(history, tracer)
    target_min = float(np.min(recharge_years))
    target_max = float(np.max(recharge_years))
    if target_min < float(years[0]) or target_max > float(years[-1]):
        raise _HistoryInputError(
            REASON_SOURCE_HISTORY_OUT_OF_WINDOW,
            f"Source history for {tracer!r} covers {years[0]:g}-{years[-1]:g}, "
            f"but the requested recharge window is {target_min:g}-{target_max:g}; extrapolation is disabled.",
        )
    return years, values


def _numeric_mapping(
    values: Mapping[str, float] | None,
    *,
    reason: str,
    name: str,
) -> dict[str, float]:
    if values is None:
        return {}
    if not isinstance(values, Mapping):
        raise _HistoryInputError(reason, f"{name} must be a mapping from tracer to finite numbers.")
    result: dict[str, float] = {}
    for raw_key, raw_value in values.items():
        key = _canonical_tracer_key(raw_key)
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as exc:
            raise _HistoryInputError(reason, f"{name}[{raw_key!r}] must be numeric.") from exc
        if not math.isfinite(value):
            raise _HistoryInputError(reason, f"{name}[{raw_key!r}] must be finite.")
        if key in result and not math.isclose(result[key], value, rel_tol=0.0, abs_tol=0.0):
            raise _HistoryInputError(reason, f"{name} contains duplicate aliases for {key!r}.")
        result[key] = value
    return result


def _history_for_tracer(histories: Mapping[str, InputHistory], tracer: str) -> InputHistory | None:
    if tracer in histories:
        return histories[tracer]
    if tracer == "3H/3He":
        return histories.get("3H")
    return None


def _used_source_history_tracers(
    tracers: Sequence[str],
    histories: Mapping[str, InputHistory],
) -> list[str]:
    used: set[str] = set()
    for tracer in tracers:
        if tracer in STABLE_ISOTOPE_TRACERS or _history_for_tracer(histories, tracer) is not None:
            if tracer == "3H/3He" and "3H" in histories:
                used.add("3H")
            else:
                used.add(tracer)
    return sorted(used)


def _provenance(
    tracers: Sequence[str],
    histories: Mapping[str, InputHistory],
) -> dict[str, Any]:
    return {
        "inference_family": "time_history_ttd",
        "claim_scope": "conditional_inference",
        "field_validation_status": "not_performed",
        "source_history_tracers": _used_source_history_tracers(tracers, histories),
        "stable_isotopes_used": bool(set(tracers) & STABLE_ISOTOPE_TRACERS),
        "radioactive_tracers_used": bool(set(tracers) & RADIOACTIVE_TRACERS),
        "stable_isotope_treatment": "conservative_time_history_response",
        "fractionation_or_evaporation_model": "not_calibrated",
    }


def _result_age_array(age_grid_years: Sequence[float]) -> np.ndarray:
    try:
        ages = np.asarray(age_grid_years, dtype=float)
    except (TypeError, ValueError):
        return np.empty(0, dtype=float)
    return ages if ages.ndim == 1 else np.ravel(ages)


def _abstain_result(
    observations: Sequence[HistoryTracerObservation],
    age_grid_years: Sequence[float],
    reason: str,
    *,
    message: str,
    provenance: Mapping[str, Any] | None = None,
    diagnostics: Mapping[str, Any] | None = None,
) -> HistoryTTDResult:
    ages = _result_age_array(age_grid_years)
    tracers = tuple(observation.tracer for observation in observations)
    nan_predictions = np.full(len(tracers), np.nan, dtype=float)
    return HistoryTTDResult(
        status="ABSTAIN",
        g=np.zeros(ages.size, dtype=float),
        age_grid_years=ages,
        tracers=tracers,
        predicted=nan_predictions,
        residuals=nan_predictions.copy(),
        abstention_reasons=(reason,),
        diagnostics={"message": message, **dict(diagnostics or {})},
        provenance=dict(provenance or {}),
    )


@dataclass(frozen=True)
class HistoryTracerObservation:
    """One measured tracer value used by :func:`fit_history_ttd`.

    ``weight`` is an optional non-negative reliability multiplier.  A separate
    ``tracer_weights`` mapping supplied to :func:`fit_history_ttd` is multiplied
    by this observation-level weight.
    """

    tracer: str
    value: float
    sigma: float
    weight: float = 1.0

    def __post_init__(self) -> None:
        tracer = _canonical_tracer_key(self.tracer)
        if not tracer:
            raise ValueError("tracer must be non-empty.")
        try:
            value = float(self.value)
            sigma = float(self.sigma)
            weight = float(self.weight)
        except (TypeError, ValueError) as exc:
            raise ValueError("value, sigma, and weight must be numeric.") from exc
        if not math.isfinite(value):
            raise ValueError(f"Observed value for {tracer} must be finite.")
        if not math.isfinite(sigma) or sigma <= 0.0:
            raise ValueError(f"Uncertainty sigma for {tracer} must be strictly positive.")
        if not math.isfinite(weight) or weight < 0.0:
            raise ValueError(f"Weight for {tracer} must be finite and non-negative.")
        object.__setattr__(self, "tracer", tracer)
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "weight", weight)


@dataclass(frozen=True)
class HistoryTTDResult:
    """Result contract for the conditional time-history TTD workflow."""

    status: str
    g: np.ndarray
    age_grid_years: np.ndarray
    tracers: tuple[str, ...]
    predicted: np.ndarray
    residuals: np.ndarray
    abstention_reasons: tuple[str, ...] = ()
    diagnostics: Mapping[str, Any] = field(default_factory=dict)
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        status = str(self.status)
        if status not in {"ESTIMATED", "ABSTAIN"}:
            raise ValueError("status must be 'ESTIMATED' or 'ABSTAIN'.")
        g = np.asarray(self.g, dtype=float).copy()
        ages = np.asarray(self.age_grid_years, dtype=float).copy()
        predicted = np.asarray(self.predicted, dtype=float).copy()
        residuals = np.asarray(self.residuals, dtype=float).copy()
        tracers = tuple(_canonical_tracer_key(tracer) for tracer in self.tracers)
        if g.ndim != 1 or ages.ndim != 1 or g.shape != ages.shape:
            raise ValueError("g and age_grid_years must be one-dimensional arrays of equal shape.")
        if predicted.shape != (len(tracers),) or residuals.shape != (len(tracers),):
            raise ValueError("predicted and residuals must align with tracers.")
        g.setflags(write=False)
        ages.setflags(write=False)
        predicted.setflags(write=False)
        residuals.setflags(write=False)
        object.__setattr__(self, "status", status)
        object.__setattr__(self, "g", g)
        object.__setattr__(self, "age_grid_years", ages)
        object.__setattr__(self, "tracers", tracers)
        object.__setattr__(self, "predicted", predicted)
        object.__setattr__(self, "residuals", residuals)
        object.__setattr__(self, "abstention_reasons", tuple(str(reason) for reason in self.abstention_reasons))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))
        object.__setattr__(self, "provenance", dict(self.provenance))


def build_history_response_matrix(
    observations: Sequence[HistoryTracerObservation] | HistoryTracerObservation,
    sample_year: float,
    age_grid_years: Sequence[float],
    source_histories: Mapping[str, InputHistory] | None,
    *,
    stable_isotope_scales: Mapping[str, float] | None = None,
    stable_isotope_offsets: Mapping[str, float] | None = None,
) -> np.ndarray:
    """Build rows mapping a simplex TTD to tracer predictions.

    Stable-isotope rows are ``scale * C_in(sample_year - age) + offset`` and
    require the complete requested recharge-year window.  Other rows delegate
    to :func:`tracer_response_kernel`.  Validation errors are raised as
    ``ValueError`` with the public reason code in the message; the fitting
    wrapper converts them into an ``ABSTAIN`` result.
    """

    obs = _coerce_observations(observations)
    year = _validate_sample_year(sample_year)
    ages = _validate_age_grid(age_grid_years)
    histories = _canonicalize_histories(source_histories)
    scales = _numeric_mapping(
        stable_isotope_scales,
        reason=REASON_INVALID_STABLE_ISOTOPE_TRANSFORM,
        name="stable_isotope_scales",
    )
    offsets = _numeric_mapping(
        stable_isotope_offsets,
        reason=REASON_INVALID_STABLE_ISOTOPE_TRANSFORM,
        name="stable_isotope_offsets",
    )
    recharge_years = year - ages
    rows: list[np.ndarray] = []

    for observation in obs:
        tracer = observation.tracer
        if tracer in STABLE_ISOTOPE_TRACERS:
            history = histories.get(tracer)
            if history is None:
                raise _HistoryInputError(
                    REASON_MISSING_STABLE_ISOTOPE_SOURCE_HISTORY,
                    f"A dated source history is required for stable isotope {tracer!r}.",
                )
            years, values = _validate_history_window(history, tracer, recharge_years)
            # The window check above makes this existing InputHistory
            # interpolation non-extrapolating.
            source_values, _ = history.interpolate(recharge_years)
            scale = scales.get(tracer, 1.0)
            offset = offsets.get(tracer, 0.0)
            if scale <= 0.0:
                raise _HistoryInputError(
                    REASON_INVALID_STABLE_ISOTOPE_TRANSFORM,
                    f"stable_isotope_scales[{tracer!r}] must be positive.",
                )
            row = scale * source_values + offset
        else:
            history = _history_for_tracer(histories, tracer)
            if history is not None and tracer in HISTORY_DRIVEN_TRACERS:
                _validate_history_window(history, tracer, recharge_years)
            try:
                row = tracer_response_kernel(
                    tracer,
                    ages,
                    year,
                    histories=histories,
                )
            except ValueError as exc:
                raise _HistoryInputError(
                    REASON_UNSUPPORTED_TRACER,
                    f"Could not build a response row for tracer {tracer!r}: {exc}",
                ) from exc

        row = np.asarray(row, dtype=float)
        if row.shape != ages.shape or not np.all(np.isfinite(row)):
            raise _HistoryInputError(
                REASON_NONFINITE_RESPONSE_MATRIX,
                f"Response row for tracer {tracer!r} is non-finite or has the wrong shape.",
            )
        rows.append(row)

    return np.vstack(rows)


def fit_history_ttd(
    observations: Sequence[HistoryTracerObservation] | HistoryTracerObservation,
    sample_year: float,
    age_grid_years: Sequence[float],
    source_histories: Mapping[str, InputHistory] | None,
    *,
    stable_isotope_scales: Mapping[str, float] | None = None,
    stable_isotope_offsets: Mapping[str, float] | None = None,
    tracer_weights: Mapping[str, float] | None = None,
    lambda_smoothness: float = 0.05,
    max_condition_number: float = 1e8,
) -> HistoryTTDResult:
    """Fit a simplex TTD against dated stable-isotope and existing tracers.

    Stable isotopes are used only through their supplied dated source histories;
    no radioactive decay or calibrated fractionation model is applied.  The
    estimate is therefore conditional on those histories and should not be
    interpreted as field validation.
    """

    raw_observations: tuple[HistoryTracerObservation, ...]
    try:
        raw_observations = _coerce_observations(observations)
    except _HistoryInputError as exc:
        return _abstain_result(
            (),
            age_grid_years,
            exc.reason,
            message=exc.message,
            provenance={
                "inference_family": "time_history_ttd",
                "claim_scope": "conditional_inference",
                "field_validation_status": "not_performed",
                "source_history_tracers": [],
                "stable_isotopes_used": False,
                "radioactive_tracers_used": False,
            },
        )

    try:
        year = _validate_sample_year(sample_year)
        ages = _validate_age_grid(age_grid_years)
        histories = _canonicalize_histories(source_histories)
        provenance = _provenance(tuple(observation.tracer for observation in raw_observations), histories)
        multipliers = _numeric_mapping(
            tracer_weights,
            reason=REASON_INVALID_TRACER_WEIGHT,
            name="tracer_weights",
        )
        effective_weights = np.asarray(
            [observation.weight * multipliers.get(observation.tracer, 1.0) for observation in raw_observations],
            dtype=float,
        )
        if not np.all(np.isfinite(effective_weights)) or np.any(effective_weights < 0.0):
            raise _HistoryInputError(REASON_INVALID_TRACER_WEIGHT, "Effective tracer weights must be finite and non-negative.")
        if not np.any(effective_weights > 0.0):
            raise _HistoryInputError(REASON_ZERO_EFFECTIVE_WEIGHT, "At least one effective tracer weight must be positive.")
        if not math.isfinite(float(lambda_smoothness)) or float(lambda_smoothness) < 0.0:
            raise _HistoryInputError(REASON_INVALID_TRACER_WEIGHT, "lambda_smoothness must be finite and non-negative.")
        if not math.isfinite(float(max_condition_number)) or float(max_condition_number) <= 0.0:
            raise _HistoryInputError(REASON_INVALID_TRACER_WEIGHT, "max_condition_number must be finite and positive.")

        matrix = build_history_response_matrix(
            raw_observations,
            year,
            ages,
            histories,
            stable_isotope_scales=stable_isotope_scales,
            stable_isotope_offsets=stable_isotope_offsets,
        )
    except _HistoryInputError as exc:
        try:
            histories_for_provenance = _canonicalize_histories(source_histories)
            provenance = _provenance(tuple(observation.tracer for observation in raw_observations), histories_for_provenance)
        except _HistoryInputError:
            provenance = {
                "inference_family": "time_history_ttd",
                "claim_scope": "conditional_inference",
                "field_validation_status": "not_performed",
                "source_history_tracers": [],
                "stable_isotopes_used": bool(
                    set(observation.tracer for observation in raw_observations) & STABLE_ISOTOPE_TRACERS
                ),
                "radioactive_tracers_used": bool(
                    set(observation.tracer for observation in raw_observations) & RADIOACTIVE_TRACERS
                ),
            }
        return _abstain_result(
            raw_observations,
            age_grid_years,
            exc.reason,
            message=exc.message,
            provenance=provenance,
        )
    except (TypeError, ValueError) as exc:
        return _abstain_result(
            raw_observations,
            age_grid_years,
            REASON_INVALID_OBSERVATION,
            message=str(exc),
            provenance=provenance,
        )

    grid = TTDGrid(taus=ages, weights=compute_trapezoidal_weights(ages))
    observations_values = np.asarray([observation.value for observation in raw_observations], dtype=float)
    sigmas = np.asarray([observation.sigma for observation in raw_observations], dtype=float)
    system = MultiTracerForwardSystem(
        node_id="history_ttd",
        sample_year=year,
        grid=grid,
        tracers=tuple(observation.tracer for observation in raw_observations),
        matrix=matrix,
        observations=observations_values,
        sigmas=sigmas,
        weights=effective_weights,
    )
    solver_result = solve_single_node_ttd(
        system,
        grid,
        lambda_smoothness=float(lambda_smoothness),
        max_condition_number=float(max_condition_number),
    )
    g = np.asarray(solver_result.g, dtype=float)
    predicted = matrix @ g
    residuals = predicted - observations_values
    diagnostics = {
        **dict(solver_result.diagnostics),
        "condition_number": system.condition_number(),
        "effective_rank": system.effective_rank(),
        "n_observations": len(raw_observations),
        "lambda_smoothness": float(lambda_smoothness),
        "max_condition_number": float(max_condition_number),
        "weight_semantics": "observation_weight multiplied by tracer_weights multiplier",
        "observation_sigmas": sigmas.tolist(),
        "effective_tracer_weights": effective_weights.tolist(),
        "stable_isotope_response_semantics": "dated source history evaluated at sample_year_minus_age",
    }
    provenance = {
        **provenance,
        "stable_isotope_scales": {
            tracer: float(value)
            for tracer, value in _numeric_mapping(
                stable_isotope_scales,
                reason=REASON_INVALID_STABLE_ISOTOPE_TRANSFORM,
                name="stable_isotope_scales",
            ).items()
            if tracer in STABLE_ISOTOPE_TRACERS
        },
        "stable_isotope_offsets": {
            tracer: float(value)
            for tracer, value in _numeric_mapping(
                stable_isotope_offsets,
                reason=REASON_INVALID_STABLE_ISOTOPE_TRANSFORM,
                name="stable_isotope_offsets",
            ).items()
            if tracer in STABLE_ISOTOPE_TRACERS
        },
    }
    return HistoryTTDResult(
        status=solver_result.status,
        g=g,
        age_grid_years=ages,
        tracers=tuple(observation.tracer for observation in raw_observations),
        predicted=predicted,
        residuals=residuals,
        abstention_reasons=tuple(solver_result.abstention_reasons),
        diagnostics=diagnostics,
        provenance=provenance,
    )


__all__ = [
    "HistoryTracerObservation",
    "HistoryTTDResult",
    "build_history_response_matrix",
    "fit_history_ttd",
    "STABLE_ISOTOPE_TRACERS",
    "REASON_MISSING_STABLE_ISOTOPE_SOURCE_HISTORY",
    "REASON_SOURCE_HISTORY_OUT_OF_WINDOW",
    "REASON_MALFORMED_SOURCE_HISTORY",
    "REASON_NONFINITE_SOURCE_HISTORY",
]
