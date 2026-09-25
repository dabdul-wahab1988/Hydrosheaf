"""Joint multi-tracer dynamic graph inversion.

The M9 prototype fitted one anchor tracer and then reused its transport for all
other tracers.  This module fits one shared transport parameter vector against
all active tracer observations simultaneously.  Tracer response factors are
applied in the forward design (radioactive decay, SF6 exchange, optional stable
isotope fractionation, and the declared radiocarbon correction), and LOTO
sensitivity is obtained by actually refitting after each tracer is removed.

The implementation remains a controlled-synthetic inverse component.  It does
not turn a low residual into a claim of field identifiability; holdout gates and
explicit conflict abstention are retained.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import minimize

from ..models.ttd_losses import weighted_rmse
from .dynamic_edge_kernel import build_lag_curvature_matrix, build_phase_basis_matrix, build_temporal_smoothness_matrix
from .dynamic_kernel_inversion import DynamicTTDInversionConfig, DynamicTTDRecovery
from .graph_tracer_forward import edge_id_str
from .tracer_registry import (
    TracerSpec,
    build_default_tracer_registry,
    canonicalize_tracer_alias,
)


# Stable result-contract values for downstream consumers.  The legacy
# ``reason_codes`` entries remain unchanged for compatibility; these constants
# provide a non-natural-language machine-readable contract.
TRACER_CONFLICT_REASON_CODE = "TRACER_CONFLICT"
DUPLICATE_TRACER_ALIAS_REASON_CODE = "DUPLICATE_TRACER_ALIAS"
MISSING_TRACER_HISTORY_REASON_CODE = "MISSING_TRACER_HISTORY"
INVALID_TRACER_HISTORY_REASON_CODE = "INVALID_TRACER_HISTORY"
INVALID_STABLE_ISOTOPE_RESPONSE_REASON_CODE = "INVALID_STABLE_ISOTOPE_RESPONSE"
STABLE_ISOTOPE_WARNING_CODE = "STABLE_ISOTOPE_TIME_HISTORY_CONDITIONAL"
INFERENCE_FAMILY = "multi_tracer_dynamic_graph"
CLAIM_SCOPE = "controlled_synthetic"


@dataclass(frozen=True)
class MultiTracerGraphConfig:
    """Configuration for a shared-transport, multi-tracer inverse."""

    tracer_weights: Mapping[str, float] = field(
        default_factory=lambda: {
            "d18O": 1.0,
            "d2H": 0.5,
            "3H": 2.0,
            "SF6": 1.5,
            "14C": 1.0,
            "chemistry": 0.2,
        }
    )
    inversion_config: DynamicTTDInversionConfig = field(default_factory=DynamicTTDInversionConfig)
    q_carbon_correction: float = 0.85
    sf6_excess_air: float = 3.0e-3
    sf6_exchange_timescale_days: float = 365.25
    stable_isotope_fractionation: Mapping[str, float] = field(
        default_factory=lambda: {"d18O": 1.0, "d2H": 1.0}
    )
    conflict_z_score_threshold: float = 3.0
    min_required_tracers: int = 1
    enforce_holdout_gate: bool = False
    max_normalized_holdout_rmse: float = 3.0

    def __post_init__(self) -> None:
        if self.q_carbon_correction <= 0.0:
            raise ValueError("q_carbon_correction must be positive")
        if self.sf6_excess_air < 0.0 or self.sf6_exchange_timescale_days <= 0.0:
            raise ValueError("SF6 exchange parameters are invalid")
        if self.min_required_tracers < 1:
            raise ValueError("min_required_tracers must be at least one")
        if self.max_normalized_holdout_rmse <= 0.0:
            raise ValueError("max_normalized_holdout_rmse must be positive")


@dataclass(frozen=True)
class MultiTracerNodeRecovery:
    """Joint multi-tracer recovery result across a target node."""

    target_node: str
    status: str  # "RECOVERED" or "ABSTAIN"
    reason_codes: tuple[str, ...]
    edge_recoveries: Mapping[str, DynamicTTDRecovery]
    active_tracers: tuple[str, ...]
    per_tracer_rmse: Mapping[str, float]
    loto_sensitivity: Mapping[str, Optional[float]]
    effective_rank: float
    condition_number: float
    conflict_detected: bool
    held_out_per_tracer_rmse: Mapping[str, float] = field(default_factory=dict)
    held_out_per_tracer_r2: Mapping[str, Optional[float]] = field(default_factory=dict)
    held_out_per_tracer_normalized_rmse: Mapping[str, float] = field(default_factory=dict)
    diagnostics: Mapping[str, Any] = field(default_factory=dict)
    reason_code: Optional[str] = None
    provenance: Mapping[str, Any] = field(default_factory=dict)


@dataclass
class _JointFit:
    parent_ids: tuple[str, ...]
    lags: np.ndarray
    basis_mat: np.ndarray
    block_map: dict[str, list[int]]
    weights: np.ndarray
    condition_number: float
    effective_rank: float
    local_fraction: np.ndarray
    edge_recoveries: dict[str, DynamicTTDRecovery]
    calibration_predictions: dict[str, np.ndarray]
    calibration_values: dict[str, np.ndarray]
    full_predictions: dict[str, np.ndarray]
    calibration_times: np.ndarray
    t_total: int
    diagnostics: dict[str, Any]


def _fill_series(values: Sequence[float], length: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    finite = np.flatnonzero(np.isfinite(arr))
    if len(finite) == 0:
        raise ValueError("source series has no finite observations")
    if len(finite) == 1:
        return np.full(length, float(arr[finite[0]]))
    return np.interp(
        np.arange(length, dtype=float),
        finite.astype(float),
        arr[finite],
        left=float(arr[finite[0]]),
        right=float(arr[finite[-1]]),
    )


def _integer_times(values: Sequence[int | float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    rounded = np.rint(arr).astype(int)
    if not np.all(np.isfinite(arr)) or not np.allclose(arr, rounded) or np.any(rounded < 0):
        raise ValueError("time indices must be finite, non-negative integers")
    return rounded


def _series_at_times(values: Sequence[float], times: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if len(arr) == len(times):
        return arr.copy()
    if len(arr) <= int(np.max(times, initial=-1)):
        raise ValueError("series is shorter than requested time indices")
    return arr[times]


def _local_series_for_tracer(
    local_tracer_inputs: Optional[Mapping[str, Any]], tracer: str
) -> Optional[Sequence[float]]:
    if not local_tracer_inputs:
        return None
    # Preferred contract: {tracer_id: series}.
    if tracer in local_tracer_inputs and not isinstance(local_tracer_inputs[tracer], Mapping):
        return local_tracer_inputs[tracer]
    canonical = _canonical_tracer_id(tracer)
    for key, value in local_tracer_inputs.items():
        if _canonical_tracer_id(key) == canonical and not isinstance(value, Mapping):
            return value
    # Backwards-compatible node contract: {"R": {tracer_id: series}}.
    for value in local_tracer_inputs.values():
        if isinstance(value, Mapping):
            series = _history_value(value, tracer)
            if series is not None:
                return series
    return None


def _canonical_tracer_id(tracer: str) -> str:
    """Canonicalize aliases used in time-history mappings for this solver."""
    return canonicalize_tracer_alias(tracer)


def _canonicalize_tracer_mapping(
    observations: Optional[Mapping[str, Sequence[float]]],
) -> tuple[dict[str, Sequence[float]], tuple[str, ...]]:
    """Canonicalize tracer keys and report aliases that collide.

    A collision is not resolved by last-write-wins: two spellings for one
    tracer can carry different observations, and silently selecting one would
    make the inverse result depend on dictionary insertion order.
    """
    if observations is None:
        return {}, ()
    canonical: dict[str, Sequence[float]] = {}
    duplicates: list[str] = []
    for raw_tracer, values in observations.items():
        tracer = canonicalize_tracer_alias(raw_tracer)
        if tracer in canonical:
            duplicates.append(tracer)
            continue
        canonical[tracer] = values
    return canonical, tuple(sorted(set(duplicates)))


def _history_value(history: Mapping[str, Any], tracer: str) -> Any:
    """Return a tracer series from a mapping, accepting safe aliases."""
    if tracer in history:
        return history[tracer]
    canonical = _canonical_tracer_id(tracer)
    for key, value in history.items():
        if _canonical_tracer_id(key) == canonical:
            return value
    return None


def _history_status(value: Any) -> str:
    """Classify a supplied series without changing the forward numerical path."""
    if value is None:
        return "missing"
    try:
        arr = np.asarray(value, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return "invalid"
    if arr.size == 0 or not np.any(np.isfinite(arr)):
        return "invalid"
    return "available"


def _audit_tracer_histories(
    active_tracers: Sequence[str],
    registry: Mapping[str, TracerSpec],
    candidate_parents: Mapping[str, Mapping[str, Sequence[float]]],
    local_tracer_inputs: Optional[Mapping[str, Any]],
) -> dict[str, Any]:
    """Audit required source histories before the dynamic fit consumes them.

    The historical implementation converted a missing parent history into a
    zero series.  That is numerically convenient but scientifically ambiguous,
    especially for stable isotopes.  This audit fails closed for required
    tracer histories while leaving optional local recharge inputs optional.
    """
    required = [
        tracer
        for tracer in active_tracers
        if tracer in registry and registry[tracer].source_history_required
    ]
    status_by_tracer: dict[str, dict[str, Any]] = {}
    missing: list[str] = []
    invalid: list[str] = []

    if candidate_parents and required:
        for tracer in required:
            parent_status: dict[str, str] = {}
            for parent, histories in candidate_parents.items():
                if not isinstance(histories, Mapping):
                    status = "invalid"
                else:
                    status = _history_status(_history_value(histories, tracer))
                parent_status[str(parent)] = status
                if status == "missing":
                    missing.append(f"{tracer}:{parent}")
                elif status == "invalid":
                    invalid.append(f"{tracer}:{parent}")
            status_by_tracer[tracer] = {"parent_histories": parent_status}

    local_status: dict[str, str] = {}
    if local_tracer_inputs is not None and not isinstance(local_tracer_inputs, Mapping):
        invalid.append("local_tracer_inputs")
    elif isinstance(local_tracer_inputs, Mapping):
        for tracer in active_tracers:
            raw = _local_series_for_tracer(local_tracer_inputs, tracer)
            if raw is not None:
                status = _history_status(raw)
                local_status[tracer] = status
                if status == "invalid":
                    invalid.append(f"local:{tracer}")

    for tracer in active_tracers:
        status_by_tracer.setdefault(tracer, {})["local_history"] = local_status.get(tracer, "not_supplied")

    available: Optional[bool]
    if not required:
        available = None
    elif not candidate_parents:
        # The fit has its own stable ``no_candidate_parents`` failure reason.
        # Do not replace it with a source-history code when no history path was
        # supplied at all.
        available = None
    else:
        available = not missing and not invalid

    return {
        "required_tracers": list(required),
        "source_history_required": bool(required),
        "source_history_available": available,
        "missing": list(dict.fromkeys(missing)),
        "invalid": list(dict.fromkeys(invalid)),
        "by_tracer": status_by_tracer,
    }


def _build_provenance(
    registry: Mapping[str, TracerSpec],
    active_tracers: Sequence[str],
    *,
    holdout_evaluated: bool = False,
    history_audit: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Build the additive node-level provenance contract for this solver."""
    active = tuple(sorted(active_tracers))
    active_kinds = {tracer: registry[tracer].kind for tracer in active if tracer in registry}
    kinds = set(active_kinds.values())
    stable_used = "stable_isotope" in kinds
    required = bool(
        history_audit.get("source_history_required", False)
        if history_audit is not None
        else any(registry[tracer].source_history_required for tracer in active if tracer in registry)
    )
    available = (
        history_audit.get("source_history_available")
        if history_audit is not None
        else None
    )
    warnings: list[str] = []
    if stable_used:
        # This solver uses observed time-series histories and a scalar
        # fractionation factor.  It is not a field-validated isotope
        # fractionation model, and d-excess is deliberately not a conflict gate.
        warnings.append(STABLE_ISOTOPE_WARNING_CODE)
        if available is None:
            warnings.append("STABLE_ISOTOPE_SOURCE_HISTORY_UNVERIFIED")
        elif available is False:
            warnings.append("STABLE_ISOTOPE_SOURCE_HISTORY_UNAVAILABLE")
    return {
        "inference_family": INFERENCE_FAMILY,
        "claim_scope": CLAIM_SCOPE,
        "field_validation_status": "not_performed",
        "active_tracers": list(active),
        "active_tracer_count": len(active),
        "active_tracer_kinds": active_kinds,
        "stable_isotopes_used": stable_used,
        "radioactive_tracers_used": "radioactive" in kinds,
        "source_history_required": required,
        "source_history_available": available,
        "warnings": warnings,
        "history_audit": dict(history_audit or {}),
        "holdout_evaluated": bool(holdout_evaluated),
    }


def _response_factor(spec: TracerSpec, tracer: str, lags: np.ndarray, cfg: MultiTracerGraphConfig) -> np.ndarray:
    factor = np.ones(len(lags), dtype=float)
    # Stable water isotopes never receive a radioactive decay factor, even if
    # a caller supplies a malformed custom registry entry.
    if spec.kind != "stable_isotope" and spec.decay_constant_per_day is not None:
        factor *= np.exp(-float(spec.decay_constant_per_day) * lags)
    if tracer == "SF6" and spec.response_model == "gas_exchange":
        factor *= 1.0 + cfg.sf6_excess_air * np.exp(-lags / cfg.sf6_exchange_timescale_days)
    if spec.kind == "stable_isotope":
        factor *= float(cfg.stable_isotope_fractionation.get(tracer, 1.0))
    if tracer == "14C" and spec.response_model == "carbonate_corrected":
        factor *= cfg.q_carbon_correction
    return factor


def _check_tracer_conflict(
    tracer_observations: Mapping[str, np.ndarray],
    registry: Mapping[str, TracerSpec],
    z_threshold: float = 3.0,
) -> tuple[bool, str]:
    """Detect declared physically incompatible tracer evidence."""
    del registry, z_threshold  # reserved for future uncertainty-aware gates
    has_sf6 = "SF6" in tracer_observations
    has_3h = "3H" in tracer_observations
    has_14c = "14C" in tracer_observations

    if has_sf6 and has_3h:
        sf6_vals = np.asarray(tracer_observations["SF6"], dtype=float)
        h3_vals = np.asarray(tracer_observations["3H"], dtype=float)
        fs, fh = sf6_vals[np.isfinite(sf6_vals)], h3_vals[np.isfinite(h3_vals)]
        if len(fs) and len(fh) and float(np.mean(fs)) > 5.0 and float(np.mean(fh)) < 0.05:
            return True, "conflict: modern SF6 detected with zero/sub-detection tritium"

    if has_sf6 and has_14c:
        sf6_vals = np.asarray(tracer_observations["SF6"], dtype=float)
        c14_vals = np.asarray(tracer_observations["14C"], dtype=float)
        fs, fc = sf6_vals[np.isfinite(sf6_vals)], c14_vals[np.isfinite(c14_vals)]
        if len(fs) and len(fc) and float(np.mean(fs)) > 4.0 and float(np.mean(fc)) < 10.0:
            return True, "conflict: modern atmospheric SF6 co-occurs with dead radiocarbon"
    return False, ""


def _fit_joint_transport(
    target_node: str,
    target_times: np.ndarray,
    tracer_ids: Sequence[str],
    tracer_observations: Mapping[str, Sequence[float]],
    candidate_parents: Mapping[str, Mapping[str, Sequence[float]]],
    local_tracer_inputs: Optional[Mapping[str, Any]],
    cfg: MultiTracerGraphConfig,
    registry: Mapping[str, TracerSpec],
    holdout_times: Optional[np.ndarray] = None,
) -> tuple[Optional[_JointFit], Optional[str]]:
    """Fit shared dynamic weights for a specified active tracer subset."""
    inv = cfg.inversion_config
    parent_ids = tuple(sorted(candidate_parents))
    if not parent_ids:
        return None, "no_candidate_parents"
    if len(target_times) < inv.min_observed_samples:
        return None, "insufficient_calibration_samples"

    max_t = int(np.max(target_times))
    if holdout_times is not None and len(holdout_times):
        max_t = max(max_t, int(np.max(holdout_times)))
    for parent in parent_ids:
        for series in candidate_parents[parent].values():
            max_t = max(max_t, len(series) - 1)
    if local_tracer_inputs:
        for value in local_tracer_inputs.values():
            if isinstance(value, Mapping):
                max_t = max(max_t, max((len(np.asarray(v).reshape(-1)) for v in value.values()), default=1) - 1)
            elif hasattr(value, "__len__"):
                max_t = max(max_t, len(value) - 1)
    t_total = max_t + 1

    lags = np.asarray(inv.lag_grid_days, dtype=float)
    lag_steps = np.rint(lags / inv.step_days).astype(int)
    if inv.kernel_mode == "phase":
        n_bases = inv.n_phase_bins
        basis = build_phase_basis_matrix(np.arange(t_total), inv.season_period_steps, n_bases)
        phases = np.floor(np.mod(target_times, inv.season_period_steps) / inv.season_period_steps * n_bases).astype(int)
        counts = [int(np.sum(phases == b)) for b in range(n_bases)]
        if sum(c >= inv.min_pairs_per_phase for c in counts) < inv.min_identified_phases:
            return None, "insufficient_phase_coverage"
    elif inv.kernel_mode == "stationary":
        n_bases = 1
        basis = np.ones((t_total, 1), dtype=float)
        counts = [len(target_times)]
    else:
        n_bases = t_total
        basis = np.eye(t_total, dtype=float)
        counts = [1] * n_bases

    # Shared feature layout: local phase fractions followed by parent phase/lag
    # coefficients.  Every tracer receives the same columns, with its own
    # response factor and source history.
    local_available = any(_local_series_for_tracer(local_tracer_inputs, tr) is not None for tr in tracer_ids)
    block_map: dict[str, list[int]] = {}
    n_features = (n_bases if local_available else 0) + len(parent_ids) * n_bases * len(lags)
    next_index = 0
    if local_available:
        block_map["local"] = list(range(next_index, next_index + n_bases))
        next_index += n_bases
    for parent in parent_ids:
        block_map[parent] = list(range(next_index, next_index + n_bases * len(lags)))
        next_index += n_bases * len(lags)

    design_rows: list[np.ndarray] = []
    response_rows: list[float] = []
    tracer_row_slices: dict[str, slice] = {}
    unweighted_full_design: dict[str, np.ndarray] = {}
    calibration_values: dict[str, np.ndarray] = {}
    row_start = 0

    for tracer in tracer_ids:
        spec = registry[tracer]
        y = _series_at_times(tracer_observations[tracer], target_times)
        finite_y = np.isfinite(y)
        if int(np.sum(finite_y)) < inv.min_observed_samples:
            return None, f"insufficient_samples:{tracer}"
        local_raw = _local_series_for_tracer(local_tracer_inputs, tracer)
        local = _fill_series(local_raw, t_total) if local_raw is not None else None
        parent_full: dict[str, np.ndarray] = {}
        for parent in parent_ids:
            raw = (
                _history_value(candidate_parents[parent], tracer)
                if isinstance(candidate_parents[parent], Mapping)
                else None
            )
            if raw is None:
                parent_full[parent] = np.zeros(t_total, dtype=float)
            else:
                parent_full[parent] = _fill_series(raw, t_total)
        factors = _response_factor(spec, tracer, lags, cfg)

        tracer_design_full = np.zeros((t_total, n_features), dtype=float)
        if local is not None:
            local_factor = float(_response_factor(spec, tracer, np.array([0.0]), cfg)[0])
            for b in range(n_bases):
                tracer_design_full[:, block_map["local"][b]] = basis[:, b] * local * local_factor
        for parent in parent_ids:
            for b in range(n_bases):
                for j, lag_step in enumerate(lag_steps):
                    shifted = np.empty(t_total, dtype=float)
                    if lag_step == 0:
                        shifted[:] = parent_full[parent]
                    else:
                        shifted[:lag_step] = parent_full[parent][0]
                        shifted[lag_step:] = parent_full[parent][:-lag_step]
                    col = block_map[parent][b * len(lags) + j]
                    tracer_design_full[:, col] = basis[:, b] * shifted * factors[j]

        weight = max(float(cfg.tracer_weights.get(tracer, 1.0)), 0.0)
        weight /= max(float(spec.measurement_sd) ** 2, 1e-12)
        scale = math.sqrt(max(weight, 1e-12))
        valid_rows = np.flatnonzero(finite_y)
        design_rows.extend((tracer_design_full[valid_rows] * scale).tolist())
        response_rows.extend((y[valid_rows] * scale).tolist())
        tracer_row_slices[tracer] = slice(row_start, row_start + len(valid_rows))
        row_start += len(valid_rows)
        calibration_values[tracer] = y
        unweighted_full_design[tracer] = tracer_design_full

    design = np.asarray(design_rows, dtype=float)
    response = np.asarray(response_rows, dtype=float)
    if design.ndim != 2 or not np.all(np.isfinite(design)):
        return None, "nonfinite_joint_design"
    y_std = float(np.std(response))
    if y_std < 1e-10:
        return None, "no_joint_signal_variation"

    d2 = build_lag_curvature_matrix(lags)
    d2_pen = d2.T @ d2 if d2.size else np.zeros((len(lags), len(lags)))
    reg = np.zeros((n_features, n_features), dtype=float)
    for parent in parent_ids:
        idx = block_map[parent]
        for b in range(n_bases):
            chunk = idx[b * len(lags) : (b + 1) * len(lags)]
            reg[np.ix_(chunk, chunk)] += inv.lambda_lag * d2_pen
    if n_bases > 1 and inv.kernel_mode in {"phase", "time"}:
        dt = build_temporal_smoothness_matrix(np.arange(n_bases, dtype=float))
        dt_pen = dt.T @ dt if dt.size else np.zeros((n_bases, n_bases))
        for name in (["local", *parent_ids] if local_available else list(parent_ids)):
            idx = block_map[name]
            width = 1 if name == "local" else len(lags)
            for j in range(width):
                phase_idx = [idx[b * width + j] for b in range(n_bases)]
                reg[np.ix_(phase_idx, phase_idx)] += inv.lambda_time * dt_pen

    var_scale = max(len(response) * y_std**2, 1e-12)
    q = design.T @ design / var_scale + reg + 1e-5 * np.eye(n_features)
    b_vec = design.T @ response / var_scale
    for parent in parent_ids:
        b_vec[block_map[parent]] -= inv.lambda_edge_sparsity
    cond = float(np.linalg.cond(q)) if np.all(np.isfinite(q)) else float("inf")
    if not math.isfinite(cond) or cond > inv.max_condition_number:
        return None, "ill_conditioned_joint_design"

    svals = np.linalg.svd(design, compute_uv=False)
    rank = float(np.linalg.matrix_rank(design)) if len(svals) else 0.0
    a_eq = np.zeros((n_bases, n_features), dtype=float)
    for b in range(n_bases):
        if local_available:
            a_eq[b, block_map["local"][b]] = 1.0
        for parent in parent_ids:
            a_eq[b, block_map[parent][b * len(lags) : (b + 1) * len(lags)]] = 1.0

    w0 = np.zeros(n_features, dtype=float)
    rho0 = 0.15 if local_available else 0.0
    if local_available:
        for b in range(n_bases):
            w0[block_map["local"][b]] = rho0
    for parent in parent_ids:
        for b in range(n_bases):
            idx = block_map[parent][b * len(lags) : (b + 1) * len(lags)]
            w0[idx] = (1.0 - rho0) / max(len(parent_ids) * len(lags), 1)

    opt = minimize(
        lambda w: 0.5 * float(w @ q @ w) - float(b_vec @ w),
        w0,
        jac=lambda w: q @ w - b_vec,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * n_features,
        constraints=[
            {"type": "eq", "fun": lambda w, row=row: float(row @ w - 1.0), "jac": lambda w, row=row: row}
            for row in a_eq
        ],
        options={"maxiter": max(500, inv.max_iterations * 20), "ftol": 1e-10},
    )
    if not opt.success or not np.all(np.isfinite(opt.x)):
        return None, f"optimization_failed:{opt.message}"
    weights = np.clip(np.asarray(opt.x, dtype=float), 0.0, 1.0)
    for b in range(n_bases):
        idx = []
        if local_available:
            idx.append(block_map["local"][b])
        for parent in parent_ids:
            idx.extend(block_map[parent][b * len(lags) : (b + 1) * len(lags)])
        total = float(np.sum(weights[idx]))
        if total <= 1e-12:
            return None, "degenerate_zero_weights"
        weights[idx] /= total

    local_fraction = (
        basis @ np.array([weights[block_map["local"][b]] for b in range(n_bases)])
        if local_available
        else np.zeros(t_total, dtype=float)
    )
    full_predictions = {tracer: unweighted_full_design[tracer] @ weights for tracer in tracer_ids}
    calibration_predictions = {tracer: full_predictions[tracer][target_times] for tracer in tracer_ids}

    edge_recoveries: dict[str, DynamicTTDRecovery] = {}
    for parent in parent_ids:
        eid = edge_id_str(parent, target_node)
        idx = block_map[parent]
        basis_weights = weights[idx].reshape(n_bases, len(lags))
        phase_mass = np.sum(basis_weights, axis=1)
        pi_t = basis @ phase_mass
        if float(np.mean(pi_t)) < inv.edge_selection_threshold:
            edge_recoveries[eid] = DynamicTTDRecovery(
                method_id="hydrosheaf_joint_multitracer_v2",
                status="ABSTAIN",
                reason_codes=(f"edge_pruned_weight_below_{inv.edge_selection_threshold}",),
                edge_id=eid,
                target_node=target_node,
                condition_number=cond,
                effective_rank=rank,
            )
            continue
        raw = basis @ basis_weights
        row_sum = np.sum(raw, axis=1, keepdims=True)
        kernel = np.divide(raw, row_sum, out=np.full_like(raw, 1.0 / len(lags)), where=row_sum > 1e-12)
        kernel = np.maximum(kernel, 0.0)
        kernel /= np.sum(kernel, axis=1, keepdims=True)
        young = np.sum(kernel[:, lags <= inv.young_water_cutoff_days], axis=1)
        mean_age = kernel @ lags
        edge_recoveries[eid] = DynamicTTDRecovery(
            method_id="hydrosheaf_joint_multitracer_v2",
            status="RECOVERED",
            reason_codes=(),
            edge_id=eid,
            target_node=target_node,
            estimated_kernel=kernel,
            estimated_mixing_fractions=pi_t,
            local_recharge_fraction=float(np.mean(local_fraction)),
            local_recharge_fractions=local_fraction,
            young_water_fraction=float(np.mean(young)),
            young_water_interval=(float(np.min(young)), float(np.max(young))),
            mean_age=float(np.mean(mean_age)),
            mean_age_interval=(float(np.min(mean_age)), float(np.max(mean_age))),
            condition_number=cond,
            effective_rank=rank,
            residual_metrics={},
            diagnostics={"joint_fit": True, "phase_edge_mass": phase_mass.tolist()},
        )

    return _JointFit(
        parent_ids=parent_ids,
        lags=lags,
        basis_mat=basis,
        block_map=block_map,
        weights=weights,
        condition_number=cond,
        effective_rank=rank,
        local_fraction=local_fraction,
        edge_recoveries=edge_recoveries,
        calibration_predictions=calibration_predictions,
        calibration_values=calibration_values,
        full_predictions=full_predictions,
        calibration_times=target_times,
        t_total=t_total,
        diagnostics={
            "joint_fit": True,
            "n_features": n_features,
            "n_rows": len(response),
            "phase_sample_counts": counts,
            "response_models": {tracer: registry[tracer].response_model for tracer in tracer_ids},
            "tracer_weights": {tracer: float(cfg.tracer_weights.get(tracer, 1.0)) for tracer in tracer_ids},
        },
    ), None


def _holdout_metrics(
    fit: _JointFit,
    holdout_times: Optional[Sequence[int | float]],
    holdout_observations: Optional[Mapping[str, Sequence[float]]],
    tracer_ids: Sequence[str],
) -> tuple[dict[str, float], dict[str, Optional[float]]]:
    if holdout_times is None or not holdout_observations:
        return {}, {}
    times = _integer_times(holdout_times)
    rmse: dict[str, float] = {}
    r2: dict[str, Optional[float]] = {}
    for tracer in tracer_ids:
        if tracer not in holdout_observations:
            continue
        obs = _series_at_times(holdout_observations[tracer], times)
        pred = fit.full_predictions[tracer][times]
        valid = np.isfinite(obs) & np.isfinite(pred)
        if not np.any(valid):
            continue
        rmse[tracer] = weighted_rmse(pred[valid], obs[valid])
        tss = float(np.sum((obs[valid] - np.mean(obs[valid])) ** 2))
        r2[tracer] = float(1.0 - np.sum((pred[valid] - obs[valid]) ** 2) / tss) if tss > 1e-12 else None
    return rmse, r2


def solve_joint_multitracer_node_inversion(
    target_node: str,
    target_times: Sequence[int | float],
    tracer_observations: Mapping[str, Sequence[float]],
    candidate_parents: Mapping[str, Mapping[str, Sequence[float]]],
    local_tracer_inputs: Optional[Mapping[str, Any]] = None,
    config: Optional[MultiTracerGraphConfig] = None,
    tracer_registry: Optional[Mapping[str, TracerSpec]] = None,
    holdout_times: Optional[Sequence[int | float]] = None,
    holdout_tracer_observations: Optional[Mapping[str, Sequence[float]]] = None,
) -> MultiTracerNodeRecovery:
    """Fit shared transport parameters jointly across all active tracers."""
    cfg = config or MultiTracerGraphConfig()
    registry = tracer_registry or build_default_tracer_registry()
    canonical_observations, duplicate_aliases = _canonicalize_tracer_mapping(tracer_observations)
    canonical_holdout_observations, holdout_duplicate_aliases = _canonicalize_tracer_mapping(
        holdout_tracer_observations
    )
    duplicate_aliases = tuple(sorted(set(duplicate_aliases + holdout_duplicate_aliases)))
    active = tuple(sorted(t for t in canonical_observations if t in registry and registry[t].enabled))
    history_audit = _audit_tracer_histories(
        active,
        registry,
        candidate_parents,
        local_tracer_inputs,
    )

    if duplicate_aliases:
        provenance = _build_provenance(
            registry,
            active,
            history_audit=history_audit,
        )
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=(
                DUPLICATE_TRACER_ALIAS_REASON_CODE,
                f"canonical_tracers:{','.join(duplicate_aliases)}",
            ),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            reason_code=DUPLICATE_TRACER_ALIAS_REASON_CODE,
            provenance=provenance,
            diagnostics={
                "history_audit": history_audit,
                "duplicate_tracer_aliases": duplicate_aliases,
                "joint_fit": False,
            },
        )

    invalid_stable_response = tuple(
        tracer
        for tracer in active
        if registry[tracer].kind == "stable_isotope"
        and registry[tracer].decay_constant_per_day is not None
    )

    if invalid_stable_response:
        provenance = _build_provenance(
            registry,
            active,
            history_audit=history_audit,
        )
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=(
                INVALID_STABLE_ISOTOPE_RESPONSE_REASON_CODE,
                f"stable_isotope_decay_constant:{','.join(invalid_stable_response)}",
            ),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            reason_code=INVALID_STABLE_ISOTOPE_RESPONSE_REASON_CODE,
            provenance=provenance,
            diagnostics={"history_audit": history_audit, "joint_fit": False},
        )

    history_reasons: list[str] = []
    if history_audit["missing"]:
        history_reasons.append(MISSING_TRACER_HISTORY_REASON_CODE)
    if history_audit["invalid"]:
        history_reasons.append(INVALID_TRACER_HISTORY_REASON_CODE)
    if history_reasons:
        provenance = _build_provenance(
            registry,
            active,
            history_audit=history_audit,
        )
        detail_reasons = tuple(
            f"{kind.lower()}:{','.join(history_audit[kind.lower()])}"
            for kind in ("MISSING", "INVALID")
            if history_audit[kind.lower()]
        )
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=tuple(history_reasons) + detail_reasons,
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            reason_code=history_reasons[0],
            provenance=provenance,
            diagnostics={"history_audit": history_audit, "joint_fit": False},
        )

    try:
        times = _integer_times(target_times)
    except ValueError as exc:
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=(f"invalid_target_times:{exc}",),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            provenance=_build_provenance(registry, active, history_audit=history_audit),
            diagnostics={"history_audit": history_audit},
        )

    if len(active) < cfg.min_required_tracers:
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=("insufficient_tracers_available",),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            provenance=_build_provenance(registry, active, history_audit=history_audit),
            diagnostics={"history_audit": history_audit},
        )

    conflict, conflict_reason = _check_tracer_conflict(
        {tracer: np.asarray(canonical_observations[tracer], dtype=float) for tracer in active},
        registry,
        cfg.conflict_z_score_threshold,
    )
    if conflict:
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=(TRACER_CONFLICT_REASON_CODE, "tracer_conflict_detected", conflict_reason),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=True,
            reason_code=TRACER_CONFLICT_REASON_CODE,
            provenance=_build_provenance(registry, active, history_audit=history_audit),
            diagnostics={
                "conflict_reason": conflict_reason,
                "history_audit": history_audit,
                "joint_fit": False,
            },
        )

    try:
        holdout_time_indices = (
            _integer_times(holdout_times)
            if holdout_times is not None and len(holdout_times)
            else None
        )
        fit, failure = _fit_joint_transport(
            target_node,
            times,
            active,
            canonical_observations,
            candidate_parents,
            local_tracer_inputs,
            cfg,
            registry,
            holdout_times=holdout_time_indices,
        )
    except (IndexError, TypeError, ValueError) as exc:
        provenance = _build_provenance(
            registry,
            active,
            history_audit=history_audit,
        )
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=(INVALID_TRACER_HISTORY_REASON_CODE, f"history_validation:{exc}"),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            reason_code=INVALID_TRACER_HISTORY_REASON_CODE,
            provenance=provenance,
            diagnostics={"history_audit": history_audit, "joint_fit": False},
        )
    if fit is None:
        failure_code = (
            INVALID_TRACER_HISTORY_REASON_CODE
            if failure and (
                failure.startswith("insufficient_samples:")
                or failure.startswith("nonfinite_joint_design")
            )
            else None
        )
        failure_reasons = (failure_code, failure) if failure_code else (failure or "joint_fit_failed",)
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=tuple(reason for reason in failure_reasons if reason),
            edge_recoveries={},
            active_tracers=active,
            per_tracer_rmse={},
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=0.0,
            condition_number=float("inf"),
            conflict_detected=False,
            reason_code=failure_code,
            provenance=_build_provenance(registry, active, history_audit=history_audit),
            diagnostics={"history_audit": history_audit, "joint_fit": False},
        )

    per_rmse = {
        tracer: weighted_rmse(fit.calibration_predictions[tracer], fit.calibration_values[tracer])
        for tracer in active
    }
    hold_rmse, hold_r2 = _holdout_metrics(
        fit,
        holdout_times,
        canonical_holdout_observations,
        active,
    )
    hold_normalized: dict[str, float] = {}
    for tracer, value in hold_rmse.items():
        observed = _series_at_times(
            canonical_holdout_observations[tracer],
            _integer_times(holdout_times),
        )
        finite = observed[np.isfinite(observed)]
        scale = max(float(np.std(finite)) if len(finite) else 0.0, float(registry[tracer].measurement_sd), 1e-8)
        hold_normalized[tracer] = float(value / scale)
    failed_holdout = [
        tracer for tracer, value in hold_normalized.items() if value > cfg.max_normalized_holdout_rmse
    ]
    if failed_holdout and cfg.enforce_holdout_gate:
        return MultiTracerNodeRecovery(
            target_node=target_node,
            status="ABSTAIN",
            reason_codes=(f"negative_joint_holdout_skill:{','.join(failed_holdout)}",),
            edge_recoveries=fit.edge_recoveries,
            active_tracers=active,
            per_tracer_rmse=per_rmse,
            loto_sensitivity={tracer: None for tracer in active},
            effective_rank=fit.effective_rank,
            condition_number=fit.condition_number,
            conflict_detected=False,
            held_out_per_tracer_rmse=hold_rmse,
            held_out_per_tracer_r2=hold_r2,
            held_out_per_tracer_normalized_rmse=hold_normalized,
            provenance=_build_provenance(
                registry,
                active,
                holdout_evaluated=True,
                history_audit=history_audit,
            ),
            diagnostics={
                **fit.diagnostics,
                "history_audit": history_audit,
                "provenance_warnings": _build_provenance(
                    registry,
                    active,
                    history_audit=history_audit,
                )["warnings"],
                "holdout_gate_passed": False,
                "holdout_gate_enforced": True,
            },
        )

    # Actual LOTO refits.  Sensitivity is the absolute change in the mean
    # recovered young-water fraction across the retained edges.
    baseline_fy = [r.young_water_fraction for r in fit.edge_recoveries.values() if r.status == "RECOVERED" and r.young_water_fraction is not None]
    baseline = float(np.mean(baseline_fy)) if baseline_fy else float("nan")
    loto: dict[str, Optional[float]] = {}
    loto_status: dict[str, str] = {}
    for dropped in active:
        remaining = tuple(tracer for tracer in active if tracer != dropped)
        if not remaining:
            loto[dropped] = None
            loto_status[dropped] = "not_defined_single_tracer_drop"
            continue
        refit, refit_failure = _fit_joint_transport(
            target_node,
            times,
            remaining,
            canonical_observations,
            candidate_parents,
            local_tracer_inputs,
            cfg,
            registry,
            holdout_times=None,
        )
        if refit is None:
            loto[dropped] = None
            loto_status[dropped] = f"refit_failed:{refit_failure}"
            continue
        values = [r.young_water_fraction for r in refit.edge_recoveries.values() if r.status == "RECOVERED" and r.young_water_fraction is not None]
        dropped_baseline = float(np.mean(values)) if values else float("nan")
        loto[dropped] = abs(dropped_baseline - baseline) if math.isfinite(baseline) and math.isfinite(dropped_baseline) else None
        loto_status[dropped] = "refit_complete"

    # Attach the primary tracer's holdout metrics to each edge for compatibility
    # with the single-tracer evaluator, while the node record retains all tracer
    # metrics explicitly.
    primary = active[0]
    for eid, rec in list(fit.edge_recoveries.items()):
        if rec.status != "RECOVERED":
            continue
        fit.edge_recoveries[eid] = DynamicTTDRecovery(
            **{
                **rec.__dict__,
                "held_out_predictions": tuple(float(v) for v in fit.full_predictions[primary][holdout_time_indices])
                if holdout_times is not None and len(holdout_times)
                else (),
                "forecast_rmse": hold_rmse.get(primary),
                "forecast_r2": hold_r2.get(primary),
            }
        )

    return MultiTracerNodeRecovery(
        target_node=target_node,
        status="RECOVERED",
        reason_codes=(),
        edge_recoveries=fit.edge_recoveries,
        active_tracers=active,
        per_tracer_rmse=per_rmse,
        loto_sensitivity=loto,
        effective_rank=fit.effective_rank,
        condition_number=fit.condition_number,
        conflict_detected=False,
        held_out_per_tracer_rmse=hold_rmse,
        held_out_per_tracer_r2=hold_r2,
        held_out_per_tracer_normalized_rmse=hold_normalized,
        provenance=_build_provenance(
            registry,
            active,
            holdout_evaluated=bool(hold_rmse or hold_r2),
            history_audit=history_audit,
        ),
        diagnostics={
            **fit.diagnostics,
            "history_audit": history_audit,
            "provenance_warnings": _build_provenance(
                registry,
                active,
                history_audit=history_audit,
            )["warnings"],
            "holdout_gate_passed": True,
            "holdout_gate_enforced": cfg.enforce_holdout_gate,
            "holdout_gate_warning": failed_holdout,
            "loto_refit_status": loto_status,
            "baseline_young_water_fraction": baseline,
        },
    )


__all__ = [
    "CLAIM_SCOPE",
    "DUPLICATE_TRACER_ALIAS_REASON_CODE",
    "INVALID_STABLE_ISOTOPE_RESPONSE_REASON_CODE",
    "INVALID_TRACER_HISTORY_REASON_CODE",
    "INFERENCE_FAMILY",
    "MISSING_TRACER_HISTORY_REASON_CODE",
    "MultiTracerGraphConfig",
    "MultiTracerNodeRecovery",
    "STABLE_ISOTOPE_WARNING_CODE",
    "TRACER_CONFLICT_REASON_CODE",
    "solve_joint_multitracer_node_inversion",
]
