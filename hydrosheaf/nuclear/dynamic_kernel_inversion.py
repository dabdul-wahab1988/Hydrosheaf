"""Truth-blind dynamic graph TTD inversion.

The first implementation of the M9 extension exposed a useful dynamic-kernel
container, but the inverse result was not actually dynamic: parent weights were
summed over all phase bins, scenario labels were passed to the estimator, and
the serialized result omitted the fitted arrays.  This module keeps the public
API compatible while making the estimand explicit:

* every phase (or time basis) has a mass-conserving simplex;
* each edge gets a time-indexed mixing fraction ``pi_uv(t)`` and a normalized
  time-indexed lag kernel ``h_uv(tau, t)``;
* uncertainty intervals use a declared L1 residual cone and a real LP, rather
  than a hard-coded percentage interval;
* no scenario/truth hint is accepted by the estimator;
* fitted arrays and prediction/loss diagnostics are serializable.

The solver is still a regularized, finite-dimensional inverse model.  It does
not claim field validation or universal identifiability; abstention remains the
correct outcome when the declared design is not informative enough.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib
import json
import math
from typing import Any, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog, minimize

from ..models.ttd_losses import LossConfig, evaluate_composite_ttd_loss
from .dynamic_edge_kernel import (
    build_lag_curvature_matrix,
    build_phase_basis_matrix,
    build_temporal_smoothness_matrix,
)
from .graph_tracer_forward import edge_id_str


@dataclass(frozen=True)
class DynamicTTDInversionConfig:
    """Configuration for dynamic state-dependent graph TTD inversion.

    ``lag_grid_days`` is expressed in physical days.  ``step_days`` is the
    sampling interval used to convert a lag coordinate to an integer source
    index; non-grid-aligned lags are rejected instead of silently rounded.
    """

    lag_grid_days: tuple[float, ...] = (0.0, 7.0, 14.0, 28.0, 60.0, 90.0, 180.0, 365.0)
    young_water_cutoff_days: float = 90.0
    step_days: float = 1.0
    kernel_mode: str = "phase"  # "stationary", "phase", or "time"
    season_period_steps: int = 52
    n_phase_bins: int = 4
    lambda_lag: float = 0.05
    lambda_time: float = 0.05
    lambda_edge_sparsity: float = 0.01
    min_observed_samples: int = 8
    min_pairs_per_phase: int = 6
    min_identified_phases: int = 2
    min_r2: float = 0.15
    max_condition_number: float = 1e8
    residual_cone_tolerance: float = 0.05
    edge_selection_threshold: float = 0.05
    min_forecast_r2: float = 0.0
    max_iterations: int = 30
    tolerance: float = 1e-5
    require_local_input: bool = False

    def __post_init__(self) -> None:
        if self.kernel_mode not in {"stationary", "phase", "time"}:
            raise ValueError("kernel_mode must be 'stationary', 'phase', or 'time'")
        if self.step_days <= 0.0:
            raise ValueError("step_days must be positive")
        if len(self.lag_grid_days) < 2:
            raise ValueError("lag_grid_days must contain at least two values")
        lags = tuple(float(v) for v in self.lag_grid_days)
        if lags != tuple(sorted(set(lags))) or lags[0] < 0.0:
            raise ValueError("lag_grid_days must be strictly increasing and non-negative")
        if any(abs(lag / self.step_days - round(lag / self.step_days)) > 1e-7 for lag in lags):
            raise ValueError("lag_grid_days must lie on the declared step_days grid")
        if self.young_water_cutoff_days < 0.0:
            raise ValueError("young_water_cutoff_days cannot be negative")
        if self.n_phase_bins < 1 or self.season_period_steps < 1:
            raise ValueError("n_phase_bins and season_period_steps must be positive")
        if self.min_observed_samples < 1 or self.min_pairs_per_phase < 1:
            raise ValueError("sample and phase-count gates must be positive")
        if self.residual_cone_tolerance <= 0.0:
            raise ValueError("residual_cone_tolerance must be positive")


@dataclass(frozen=True)
class DynamicTTDRecovery:
    """Immutable result contract for one candidate incoming edge."""

    method_id: str
    status: str  # "RECOVERED" or "ABSTAIN"
    reason_codes: tuple[str, ...]
    edge_id: str
    target_node: str
    estimated_kernel: Optional[np.ndarray] = None  # (T, L)
    estimated_mixing_fractions: Optional[np.ndarray] = None  # pi_uv(t), shape (T,)
    local_recharge_fraction: Optional[float] = None  # scalar compatibility summary
    local_recharge_fractions: Optional[np.ndarray] = None  # rho_v(t), shape (T,)
    young_water_fraction: Optional[float] = None
    young_water_interval: Optional[tuple[float, float]] = None
    mean_age: Optional[float] = None
    mean_age_interval: Optional[tuple[float, float]] = None
    held_out_predictions: tuple[float, ...] = ()
    forecast_rmse: Optional[float] = None
    forecast_r2: Optional[float] = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)
    condition_number: Optional[float] = None
    effective_rank: Optional[float] = None
    residual_metrics: Mapping[str, float] = field(default_factory=dict)
    loss_decomposition: Mapping[str, Any] = field(default_factory=dict)
    configuration_hash: str = ""
    code_version: str = "1.1.0"
    truth_blindness_declared: bool = True

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-safe record containing all fitted arrays."""

        def serialise(value: Any) -> Any:
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, (np.floating, np.integer)):
                return value.item()
            if isinstance(value, Mapping):
                return {str(k): serialise(v) for k, v in value.items()}
            if isinstance(value, tuple):
                return [serialise(v) for v in value]
            return value

        return {
            "method_id": self.method_id,
            "status": self.status,
            "reason_codes": list(self.reason_codes),
            "edge_id": self.edge_id,
            "target_node": self.target_node,
            "estimated_kernel": serialise(self.estimated_kernel),
            "estimated_mixing_fractions": serialise(self.estimated_mixing_fractions),
            "local_recharge_fraction": self.local_recharge_fraction,
            "local_recharge_fractions": serialise(self.local_recharge_fractions),
            "young_water_fraction": self.young_water_fraction,
            "young_water_interval": serialise(self.young_water_interval),
            "mean_age": self.mean_age,
            "mean_age_interval": serialise(self.mean_age_interval),
            "held_out_predictions": serialise(self.held_out_predictions),
            "forecast_rmse": self.forecast_rmse,
            "forecast_r2": self.forecast_r2,
            "diagnostics": serialise(self.diagnostics),
            "condition_number": self.condition_number,
            "effective_rank": self.effective_rank,
            "residual_metrics": serialise(self.residual_metrics),
            "loss_decomposition": serialise(self.loss_decomposition),
            "configuration_hash": self.configuration_hash,
            "code_version": self.code_version,
            "truth_blindness_declared": self.truth_blindness_declared,
        }


def _configuration_hash(cfg: DynamicTTDInversionConfig) -> str:
    payload = {k: getattr(cfg, k) for k in cfg.__dataclass_fields__}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _fill_series(values: Sequence[float], length: int, name: str) -> np.ndarray:
    """Fill missing source measurements by linear interpolation."""

    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0 or not np.any(np.isfinite(arr)):
        raise ValueError(f"{name} has no finite observations")
    finite_idx = np.flatnonzero(np.isfinite(arr))
    finite_vals = arr[finite_idx]
    target = np.arange(length, dtype=float)
    if len(finite_idx) == 1:
        return np.full(length, float(finite_vals[0]), dtype=float)
    return np.interp(target, finite_idx.astype(float), finite_vals, left=finite_vals[0], right=finite_vals[-1])


def _as_integer_times(values: Sequence[int | float], name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain finite time indices")
    rounded = np.rint(arr).astype(int)
    if not np.allclose(arr, rounded):
        raise ValueError(f"{name} must contain integer time-step indices")
    if np.any(rounded < 0):
        raise ValueError(f"{name} cannot contain negative time indices")
    return rounded


def _build_uncertainty_bounds(
    design: np.ndarray,
    y_obs: np.ndarray,
    w_opt: np.ndarray,
    a_eq_weights: np.ndarray,
    b_eq: np.ndarray,
    objective: np.ndarray,
    residual_tolerance: float,
    bounds_weights: Sequence[tuple[float, float]],
) -> tuple[Optional[float], Optional[float], dict[str, Any]]:
    """Bound a linear functional over an exact L1 residual cone."""

    n_features = design.shape[1]
    n_obs = design.shape[0]
    residual = design @ w_opt - y_obs
    l1_star = float(np.sum(np.abs(residual)))
    max_l1 = l1_star * (1.0 + residual_tolerance) + 1e-10

    a_ub = np.vstack(
        [
            np.hstack([design, -np.eye(n_obs)]),
            np.hstack([-design, -np.eye(n_obs)]),
            np.hstack([np.zeros((1, n_features)), np.ones((1, n_obs))]),
        ]
    )
    b_ub = np.concatenate([y_obs, -y_obs, np.array([max_l1])])
    a_eq = np.hstack([a_eq_weights, np.zeros((a_eq_weights.shape[0], n_obs))])
    bounds = list(bounds_weights) + [(0.0, None)] * n_obs

    low = linprog(
        np.concatenate([objective, np.zeros(n_obs)]),
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=bounds,
        method="highs",
    )
    high = linprog(
        np.concatenate([-objective, np.zeros(n_obs)]),
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=bounds,
        method="highs",
    )
    diagnostics = {
        "uncertainty_set": "l1_residual_cone",
        "residual_l1_optimum": l1_star,
        "residual_l1_limit": max_l1,
        "lp_min_success": bool(low.success),
        "lp_max_success": bool(high.success),
    }
    if not low.success or not high.success:
        diagnostics["lp_message"] = f"min={low.message}; max={high.message}"
        return None, None, diagnostics
    return float(low.fun), float(-high.fun), diagnostics


def solve_dynamic_node_inversion(
    target_node: str,
    target_times: Sequence[int | float],
    target_values: Sequence[float],
    candidate_parents: Mapping[str, Sequence[float]],
    local_input: Optional[Sequence[float]] = None,
    config: Optional[DynamicTTDInversionConfig] = None,
    holdout_times: Optional[Sequence[int | float]] = None,
    holdout_values: Optional[Sequence[float]] = None,
    scenario_hint: Optional[str] = None,
    loss_config: Optional[LossConfig] = None,
) -> dict[str, DynamicTTDRecovery]:
    """Fit a truth-blind dynamic graph node and return one record per edge.

    ``scenario_hint`` is retained only to make accidental legacy calls fail
    loudly.  Scenario labels are benchmark metadata and must never enter the
    estimator; callers should omit the argument.
    """

    if scenario_hint is not None:
        raise ValueError(
            "scenario_hint is not accepted by the truth-blind estimator; "
            "perform scenario classification in the scorer"
        )

    cfg = config or DynamicTTDInversionConfig()
    cfg_hash = _configuration_hash(cfg)
    parent_ids = tuple(sorted(candidate_parents))
    edge_names = [edge_id_str(parent, target_node) for parent in parent_ids]

    def make_abstain(edge_id: str, reason: str, cond: Optional[float] = None) -> DynamicTTDRecovery:
        return DynamicTTDRecovery(
            method_id="hydrosheaf_dynamic_ttd_inversion_v2",
            status="ABSTAIN",
            reason_codes=(reason,),
            edge_id=edge_id,
            target_node=target_node,
            condition_number=cond,
            configuration_hash=cfg_hash,
        )

    if target_node == "R":
        return {eid: make_abstain(eid, "invalid_flowpath_recharge_boundary_cannot_have_parents") for eid in edge_names}
    if not edge_names:
        return {}
    if cfg.require_local_input and local_input is None:
        return {eid: make_abstain(eid, "local_input_unavailable_for_mass_conservation") for eid in edge_names}

    try:
        all_times = _as_integer_times(target_times, "target_times")
        y_all = np.asarray(target_values, dtype=float).reshape(-1)
    except (TypeError, ValueError) as exc:
        return {eid: make_abstain(eid, f"invalid_input:{exc}") for eid in edge_names}
    if len(all_times) != len(y_all):
        return {eid: make_abstain(eid, "target_times_values_length_mismatch") for eid in edge_names}

    valid = np.isfinite(y_all)
    t_obs = all_times[valid]
    y_obs = y_all[valid]
    if len(y_obs) < cfg.min_observed_samples:
        return {eid: make_abstain(eid, "insufficient_calibration_samples") for eid in edge_names}
    order = np.argsort(t_obs, kind="stable")
    t_obs = t_obs[order]
    y_obs = y_obs[order]
    y_std = float(np.std(y_obs))
    if y_std < 1e-8:
        return {eid: make_abstain(eid, "no_target_variation") for eid in edge_names}

    max_t = int(np.max(t_obs))
    if holdout_times is not None and len(holdout_times):
        ho_times_raw = _as_integer_times(holdout_times, "holdout_times")
        max_t = max(max_t, int(np.max(ho_times_raw)))
    else:
        ho_times_raw = np.array([], dtype=int)
    if local_input is not None:
        max_t = max(max_t, len(local_input) - 1)
    for series in candidate_parents.values():
        max_t = max(max_t, len(series) - 1)
    t_total = max_t + 1

    lags = np.asarray(cfg.lag_grid_days, dtype=float)
    n_lags = len(lags)
    lag_steps = np.rint(lags / cfg.step_days).astype(int)

    if cfg.kernel_mode == "phase":
        n_bases = cfg.n_phase_bins
        basis_mat = build_phase_basis_matrix(np.arange(t_total), cfg.season_period_steps, n_bases)
        obs_phases = np.floor(np.mod(t_obs, cfg.season_period_steps) / cfg.season_period_steps * n_bases).astype(int)
        counts = [int(np.sum(obs_phases == b)) for b in range(n_bases)]
        represented = sum(count >= cfg.min_pairs_per_phase for count in counts)
        if represented < cfg.min_identified_phases:
            return {eid: make_abstain(eid, "insufficient_phase_coverage") for eid in edge_names}
    elif cfg.kernel_mode == "stationary":
        n_bases = 1
        basis_mat = np.ones((t_total, 1), dtype=float)
        counts = [len(t_obs)]
    else:
        n_bases = t_total
        basis_mat = np.eye(t_total, dtype=float)
        counts = [1] * n_bases

    try:
        parent_series = {
            pid: _fill_series(candidate_parents[pid], t_total, f"parent {pid}") for pid in parent_ids
        }
        local_series = _fill_series(local_input, t_total, "local_input") if local_input is not None else None
    except ValueError as exc:
        return {eid: make_abstain(eid, f"invalid_source_series:{exc}") for eid in edge_names}

    columns: list[np.ndarray] = []
    full_columns: list[np.ndarray] = []
    block_map: dict[str, list[int]] = {}

    if local_series is not None:
        block_map["local"] = []
        for b in range(n_bases):
            full_col = basis_mat[:, b] * local_series
            block_map["local"].append(len(columns))
            columns.append(full_col[t_obs])
            full_columns.append(full_col)

    for pid in parent_ids:
        block_map[pid] = []
        source = parent_series[pid]
        for b in range(n_bases):
            for lag_step in lag_steps:
                shifted = np.empty(t_total, dtype=float)
                if lag_step == 0:
                    shifted[:] = source
                else:
                    shifted[:lag_step] = source[0]
                    shifted[lag_step:] = source[:-lag_step]
                full_col = basis_mat[:, b] * shifted
                block_map[pid].append(len(columns))
                columns.append(full_col[t_obs])
                full_columns.append(full_col)

    if not columns:
        return {eid: make_abstain(eid, "no_valid_design_columns") for eid in edge_names}
    design = np.column_stack(columns)
    full_design = np.column_stack(full_columns)
    n_features = design.shape[1]
    if not np.all(np.isfinite(design)) or not np.all(np.isfinite(full_design)):
        return {eid: make_abstain(eid, "nonfinite_design") for eid in edge_names}
    if sum(float(np.std(col)) for col in columns) < 1e-8:
        return {eid: make_abstain(eid, "forcing_has_no_variation") for eid in edge_names}

    d2 = build_lag_curvature_matrix(lags)
    d2_pen = d2.T @ d2 if d2.size else np.zeros((n_lags, n_lags), dtype=float)
    reg_mat = np.zeros((n_features, n_features), dtype=float)
    for pid in parent_ids:
        idx = block_map[pid]
        for b in range(n_bases):
            chunk = idx[b * n_lags : (b + 1) * n_lags]
            reg_mat[np.ix_(chunk, chunk)] += cfg.lambda_lag * d2_pen
    if n_bases > 1 and cfg.kernel_mode in {"phase", "time"}:
        d_time = build_temporal_smoothness_matrix(np.arange(n_bases, dtype=float))
        d_time_pen = d_time.T @ d_time if d_time.size else np.zeros((n_bases, n_bases))
        for pid in (["local", *parent_ids] if "local" in block_map else list(parent_ids)):
            idx = block_map[pid]
            width = 1 if pid == "local" else n_lags
            for j in range(width):
                phase_idx = [idx[b * width + j] for b in range(n_bases)]
                reg_mat[np.ix_(phase_idx, phase_idx)] += cfg.lambda_time * d_time_pen

    var_scale = max(len(y_obs) * y_std**2, 1e-12)
    q_mat = (design.T @ design) / var_scale + reg_mat + 1e-5 * np.eye(n_features)
    b_vec = (design.T @ y_obs) / var_scale
    for pid in parent_ids:
        b_vec[block_map[pid]] -= cfg.lambda_edge_sparsity

    try:
        cond = float(np.linalg.cond(q_mat))
    except Exception:
        cond = float("inf")
    if not math.isfinite(cond) or cond > cfg.max_condition_number:
        return {eid: make_abstain(eid, "ill_conditioned_design", cond) for eid in edge_names}

    s_vals = np.linalg.svd(design, compute_uv=False)
    eff_rank = float(np.sum(s_vals > 1e-5 * s_vals[0])) if len(s_vals) and s_vals[0] > 0 else 0.0

    def objective(w: np.ndarray) -> float:
        return 0.5 * float(w @ q_mat @ w) - float(b_vec @ w)

    def gradient(w: np.ndarray) -> np.ndarray:
        return q_mat @ w - b_vec

    a_eq = np.zeros((n_bases, n_features), dtype=float)
    for b in range(n_bases):
        if "local" in block_map:
            a_eq[b, block_map["local"][b]] = 1.0
        for pid in parent_ids:
            a_eq[b, block_map[pid][b * n_lags : (b + 1) * n_lags]] = 1.0
    b_eq = np.ones(n_bases, dtype=float)

    bounds = [(0.0, 1.0)] * n_features
    w0 = np.zeros(n_features, dtype=float)
    local_start = 0.15 if "local" in block_map else 0.0
    if "local" in block_map:
        for b in range(n_bases):
            w0[block_map["local"][b]] = local_start
    remaining = 1.0 - local_start
    for pid in parent_ids:
        for b in range(n_bases):
            chunk = block_map[pid][b * n_lags : (b + 1) * n_lags]
            w0[chunk] = remaining / max(len(parent_ids) * n_lags, 1)

    constraints = [
        {"type": "eq", "fun": lambda w, row=row: float(row @ w - 1.0), "jac": lambda w, row=row: row}
        for row in a_eq
    ]
    opt = minimize(
        objective,
        w0,
        jac=gradient,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": max(500, cfg.max_iterations * 20), "ftol": 1e-10},
    )
    if not opt.success or not np.all(np.isfinite(opt.x)):
        return {eid: make_abstain(eid, f"optimization_failed:{opt.message}", cond) for eid in edge_names}

    w_opt = np.clip(np.asarray(opt.x, dtype=float), 0.0, 1.0)
    for b in range(n_bases):
        idx = []
        if "local" in block_map:
            idx.append(block_map["local"][b])
        for pid in parent_ids:
            idx.extend(block_map[pid][b * n_lags : (b + 1) * n_lags])
        total = float(np.sum(w_opt[idx]))
        if total <= 1e-12:
            return {eid: make_abstain(eid, "degenerate_zero_weights", cond) for eid in edge_names}
        w_opt[idx] /= total

    y_pred_cal = design @ w_opt
    residual = y_obs - y_pred_cal
    rss = float(np.sum(residual**2))
    tss = float(np.sum((y_obs - np.mean(y_obs)) ** 2))
    r2 = float(1.0 - rss / tss) if tss > 1e-12 else 0.0
    rmse = float(np.sqrt(np.mean(residual**2)))
    mae = float(np.mean(np.abs(residual)))
    if r2 < cfg.min_r2:
        return {eid: make_abstain(eid, f"poor_fit_r2:{r2:.3f}", cond) for eid in edge_names}

    held_preds: list[float] = []
    fc_rmse: Optional[float] = None
    fc_r2: Optional[float] = None
    if holdout_times is not None and holdout_values is not None and len(holdout_times):
        ho_t = _as_integer_times(holdout_times, "holdout_times")
        ho_y = np.asarray(holdout_values, dtype=float).reshape(-1)
        if len(ho_t) != len(ho_y):
            return {eid: make_abstain(eid, "holdout_times_values_length_mismatch", cond) for eid in edge_names}
        hv = np.isfinite(ho_y) & (ho_t < t_total)
        if np.any(hv):
            pred_ho = full_design[ho_t[hv]] @ w_opt
            held_preds = [float(v) for v in pred_ho]
            err = pred_ho - ho_y[hv]
            fc_rmse = float(np.sqrt(np.mean(err**2)))
            ho_tss = float(np.sum((ho_y[hv] - np.mean(ho_y[hv])) ** 2))
            fc_r2 = float(1.0 - np.sum(err**2) / ho_tss) if ho_tss > 1e-12 else None
            if fc_r2 is not None and fc_r2 < cfg.min_forecast_r2:
                return {eid: make_abstain(eid, f"negative_forecast_skill_r2:{fc_r2:.3f}", cond) for eid in edge_names}

    loss_decomp: dict[str, Any] = {
        "time_rmse": rmse,
        "time_mae": mae,
        "rss": rss,
        "wasserstein_ttd": None,
        "ttd_component_status": "not_available_without_observed_ttd",
    }
    if loss_config is not None:
        eval_cfg = replace(loss_config, wasserstein_ttd_weight=0.0)
        comp = evaluate_composite_ttd_loss(
            predicted_signal=y_pred_cal,
            observed_signal=y_obs,
            time_grid=t_obs.astype(float),
            config=eval_cfg,
        )
        loss_decomp.update(comp.decomposition)

    recoveries: dict[str, DynamicTTDRecovery] = {}
    local_weights = (
        np.array([w_opt[block_map["local"][b]] for b in range(n_bases)], dtype=float)
        if "local" in block_map
        else np.zeros(n_bases, dtype=float)
    )
    rho_t = basis_mat @ local_weights

    for pid in parent_ids:
        eid = edge_id_str(pid, target_node)
        idx = block_map[pid]
        basis_weights = w_opt[idx].reshape(n_bases, n_lags)
        phase_mass = np.sum(basis_weights, axis=1)
        pi_t = basis_mat @ phase_mass
        if float(np.mean(pi_t)) < cfg.edge_selection_threshold:
            recoveries[eid] = make_abstain(eid, f"edge_pruned_weight_below_{cfg.edge_selection_threshold}", cond)
            continue

        raw_kernel = basis_mat @ basis_weights
        denom = np.sum(raw_kernel, axis=1, keepdims=True)
        h_mat = np.divide(raw_kernel, denom, out=np.full_like(raw_kernel, 1.0 / n_lags), where=denom > 1e-12)
        h_mat = np.maximum(h_mat, 0.0)
        h_mat /= np.sum(h_mat, axis=1, keepdims=True)

        young_mask = lags <= cfg.young_water_cutoff_days
        fy_series = np.sum(h_mat[:, young_mask], axis=1)
        mean_age_series = h_mat @ lags
        fy_mean = float(np.mean(fy_series))
        mean_age = float(np.mean(mean_age_series))

        point_mass = max(float(np.mean(phase_mass)), 1e-12)
        c_fy = np.zeros(n_features, dtype=float)
        c_age = np.zeros(n_features, dtype=float)
        for b in range(n_bases):
            for j, lag in enumerate(lags):
                c_fy[idx[b * n_lags + j]] = (1.0 / n_bases / point_mass) if lag <= cfg.young_water_cutoff_days else 0.0
                c_age[idx[b * n_lags + j]] = lag / n_bases / point_mass
        fy_low, fy_high, lp_fy_diag = _build_uncertainty_bounds(
            design, y_obs, w_opt, a_eq, b_eq, c_fy, cfg.residual_cone_tolerance, bounds
        )
        age_low, age_high, lp_age_diag = _build_uncertainty_bounds(
            design, y_obs, w_opt, a_eq, b_eq, c_age, cfg.residual_cone_tolerance, bounds
        )
        if fy_low is None or fy_high is None:
            fy_interval = (0.0, 1.0)
        else:
            fy_interval = (
                float(np.clip(min(fy_low, fy_mean), 0.0, 1.0)),
                float(np.clip(max(fy_high, fy_mean), 0.0, 1.0)),
            )
        if age_low is None or age_high is None:
            age_interval = (float(np.min(lags)), float(np.max(lags)))
        else:
            age_interval = (
                float(max(np.min(lags), min(age_low, mean_age))),
                float(min(np.max(lags), max(age_high, mean_age))),
            )

        recoveries[eid] = DynamicTTDRecovery(
            method_id="hydrosheaf_dynamic_ttd_inversion_v2",
            status="RECOVERED",
            reason_codes=(),
            edge_id=eid,
            target_node=target_node,
            estimated_kernel=h_mat,
            estimated_mixing_fractions=pi_t,
            local_recharge_fraction=float(np.mean(rho_t)),
            local_recharge_fractions=rho_t,
            young_water_fraction=fy_mean,
            young_water_interval=fy_interval,
            mean_age=mean_age,
            mean_age_interval=age_interval,
            held_out_predictions=tuple(held_preds),
            forecast_rmse=fc_rmse,
            forecast_r2=fc_r2,
            diagnostics={
                "kernel_mode": cfg.kernel_mode,
                "n_basis_functions": n_bases,
                "phase_sample_counts": counts,
                "phase_edge_mass": phase_mass.tolist(),
                "mixing_fraction_range": [float(np.min(pi_t)), float(np.max(pi_t))],
                "local_recharge_fraction_range": [float(np.min(rho_t)), float(np.max(rho_t))],
                "uncertainty_normalization": "point_edge_mass",
                "uncertainty_lp": {"young_water": lp_fy_diag, "mean_age": lp_age_diag},
                "optimization_iterations": int(getattr(opt, "nit", 0) or 0),
            },
            condition_number=cond,
            effective_rank=eff_rank,
            residual_metrics={"rmse_cal": rmse, "mae_cal": mae, "r2_cal": r2, "rss_cal": rss},
            loss_decomposition=loss_decomp,
            configuration_hash=cfg_hash,
        )

    return recoveries


__all__ = [
    "DynamicTTDInversionConfig",
    "DynamicTTDRecovery",
    "solve_dynamic_node_inversion",
]
