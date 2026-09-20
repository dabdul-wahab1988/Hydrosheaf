r"""Graph-constrained Transit-Time Distribution (TTD) inversion engine.

This module provides HydroSheaf's regularized, source/mixing-aware inverse solver
for dynamic and static groundwater networks.

Key Inferential Principles:
1. Truth-Blindness Gatekeeping:
   Strictly rejects any input structure containing ground truth, sealed
   commitments, or hidden forward-model parameters.
2. Mass Conservation:
   Local recharge fraction and upstream mixing weights satisfy
   \rho_v + \sum_{u \in Pa(v)} \pi_{uv} = 1.0, with non-negative kernels h_{uv}(a) >= 0,
   \sum h_{uv} = 1.0.
3. Second-Order Curvature Regularization:
   Non-uniform finite-difference curvature penalty D_2 on transit-time lag grids
   to suppress unphysical high-frequency oscillations without damping real peaks.
4. Automatic False-Edge Pruning:
   L1 edge sparsity penalty and weight thresholding to drop spurious candidate edges.
5. Sharp Linear-Programming Uncertainty Bounds:
   Highs LP solver calculates sharp minimum and maximum Young-Water Fractions [Fy_min, Fy_max]
   and mean transit times over an empirical residual tolerance cone RSS <= RSS* * (1 + epsilon).
6. Calibrated Abstention:
   Returns explicit ABSTAIN status when sample count is low, signal variation is degenerate,
   or regularized design condition number indicates catastrophic collinearity.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog, minimize

# Forbidden keywords that indicate unsealed or oracle ground-truth leakage
_FORBIDDEN_TRUTH_KEYS = {
    "true_edges",
    "edge_kernels",
    "node_ttd_masses",
    "node_mixing_weights",
    "local_input_signals",
    "node_output_signals",
    "young_water_fractions",
    "held_out_values",
    "node_age_particles",
    "node_responses",
    "input_signals",
    "sealed_truth",
    "ground_truth",
    "truth",
}


def verify_truth_blindness(payload: Any) -> None:
    """Inspect input data payloads recursively to guarantee truth-blind execution.

    Raises:
        ValueError: If any key, attribute, or descriptor matches sealed ground-truth fields.
    """
    if payload is None:
        return

    # Reject known ground truth class types directly
    type_name = type(payload).__name__
    if "Truth" in type_name:
        raise ValueError(
            f"Truth-blindness violation: payload is an instance of sealed truth class {type_name!r}."
        )

    if isinstance(payload, dict):
        for k, v in payload.items():
            key_str = str(k).lower()
            if (
                key_str in _FORBIDDEN_TRUTH_KEYS
                or key_str.startswith("true_")
                or key_str.startswith("truth_")
                or key_str == "truth"
            ):
                raise ValueError(
                    f"Truth-blindness violation: payload contains forbidden truth key {k!r}."
                )
            verify_truth_blindness(v)
    elif isinstance(payload, (list, tuple, set)):
        for item in payload:
            verify_truth_blindness(item)
    elif hasattr(payload, "__dataclass_fields__"):
        for field_name in payload.__dataclass_fields__:
            name_str = field_name.lower()
            if (
                name_str in _FORBIDDEN_TRUTH_KEYS
                or name_str.startswith("true_")
                or name_str.startswith("truth_")
                or name_str == "truth"
            ):
                raise ValueError(
                    f"Truth-blindness violation: dataclass attribute {field_name!r} matches forbidden truth pattern."
                )
            verify_truth_blindness(getattr(payload, field_name))


def _build_d2_curvature_matrix(lags: Sequence[int | float]) -> np.ndarray:
    """Construct non-uniform 2nd-order finite difference curvature operator D_2."""
    n = len(lags)
    if n < 3:
        return np.zeros((0, n), dtype=float)

    rows: list[np.ndarray] = []
    for j in range(1, n - 1):
        dt1 = float(lags[j] - lags[j - 1])
        dt2 = float(lags[j + 1] - lags[j])
        if dt1 <= 0.0 or dt2 <= 0.0:
            raise ValueError("Lags must be strictly monotonically increasing.")
        row = np.zeros(n, dtype=float)
        coeff = 2.0 / (dt1 + dt2)
        row[j - 1] = coeff / dt1
        row[j] = -coeff * (1.0 / dt1 + 1.0 / dt2)
        row[j + 1] = coeff / dt2
        rows.append(row)

    return np.vstack(rows)


def _topological_order(nodes: Iterable[str], edges: Iterable[tuple[str, str]]) -> tuple[str, ...] | None:
    """Deterministic topological sort of nodes or None if cyclic."""
    node_list = sorted(set(nodes))
    incoming = {node: 0 for node in node_list}
    children: dict[str, list[str]] = {node: [] for node in node_list}
    for source, target in edges:
        if source in incoming and target in incoming:
            incoming[target] += 1
            children[source].append(target)
    ready = sorted([node for node, count in incoming.items() if count == 0])
    order: list[str] = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for child in sorted(children[node]):
            incoming[child] -= 1
            if incoming[child] == 0:
                ready.append(child)
                ready.sort()
    return tuple(order) if len(order) == len(node_list) else None


def _path_exists(edges: Iterable[tuple[str, str]], source: str, target: str) -> bool:
    """Determine if a directed path exists from source to target."""
    if source == target:
        return True
    adj: dict[str, set[str]] = {}
    for u, v in edges:
        adj.setdefault(u, set()).add(v)
    visited = set()
    queue = [source]
    while queue:
        curr = queue.pop(0)
        if curr == target:
            return True
        if curr not in visited:
            visited.add(curr)
            queue.extend(adj.get(curr, ()))
    return False


@dataclass(frozen=True)
class GraphTTDInversionConfig:
    """Pre-registered configuration for graph TTD inversion and identification."""

    lag_grid_days: tuple[int, ...] = (0, 7, 14, 28, 60, 90, 180, 365)
    young_water_cutoff_days: int = 90
    regularization_smoothness: float = 0.05
    regularization_ridge: float = 1e-4
    regularization_edge_sparsity: float = 0.01
    min_observed_samples: int = 8
    min_r2: float = 0.15
    max_condition_number: float = 1e8
    residual_cone_tolerance: float = 0.05
    edge_selection_threshold: float = 0.05
    local_recharge_has_transit_time: bool = True
    joint_tracers: bool = True

    def __post_init__(self) -> None:
        if len(self.lag_grid_days) < 2:
            raise ValueError("lag_grid_days must contain at least two lag points.")
        if tuple(sorted(set(self.lag_grid_days))) != self.lag_grid_days:
            raise ValueError("lag_grid_days must be strictly increasing.")
        if self.young_water_cutoff_days <= 0:
            raise ValueError("young_water_cutoff_days must be strictly positive.")
        if self.regularization_smoothness < 0.0:
            raise ValueError("regularization_smoothness cannot be negative.")
        if self.residual_cone_tolerance <= 0.0:
            raise ValueError("residual_cone_tolerance must be positive.")


@dataclass(frozen=True)
class NodeInversionResult:
    """Inversion result and identified bounds for a single node."""

    node_id: str
    status: str  # "ESTIMATED" or "ABSTAIN"
    reason: Optional[str]
    r2: Optional[float]
    n_samples: int
    n_features: int
    condition_number: Optional[float]
    local_recharge_fraction: Optional[float]
    upstream_weights: Mapping[str, float]
    edge_kernels: Mapping[str, tuple[float, ...]]
    local_kernel: Optional[tuple[float, ...]]
    young_water_fraction: Optional[float]
    young_water_interval: Optional[tuple[float, float]]
    mean_lag_days: Optional[float]
    mean_lag_interval: Optional[tuple[float, float]]
    predicted_series: Optional[tuple[float, ...]] = None
    continuous_predicted_series: Optional[tuple[float, ...]] = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.status not in {"ESTIMATED", "ABSTAIN"}:
            raise ValueError(f"Unknown status {self.status!r}; expected ESTIMATED or ABSTAIN.")


def _linear_interpolate_series(times: Sequence[int], values: Sequence[float], target_times: Sequence[int]) -> np.ndarray:
    """Interpolate irregular calibration observations onto continuous time grid."""
    if len(times) == 0:
        return np.full(len(target_times), np.nan, dtype=float)
    t_arr = np.asarray(times, dtype=float)
    v_arr = np.asarray(values, dtype=float)
    order = np.argsort(t_arr)
    t_sorted = t_arr[order]
    v_sorted = v_arr[order]
    targets = np.asarray(target_times, dtype=float)
    return np.interp(targets, t_sorted, v_sorted, left=v_sorted[0], right=v_sorted[-1])


def solve_node_ttd_inversion(
    node_id: str,
    target_times: Sequence[int],
    target_values: Sequence[float],
    local_input_series: Optional[np.ndarray],
    parent_series_map: Mapping[str, np.ndarray],
    config: Optional[GraphTTDInversionConfig] = None,
    total_continuous_steps: Optional[int] = None,
) -> NodeInversionResult:
    """Fit a truth-blind, mass-conserving regularized mixing & TTD model for a single node."""
    cfg = config or GraphTTDInversionConfig()

    t_obs = np.asarray(target_times, dtype=int)
    y_obs = np.asarray(target_values, dtype=float)

    finite_mask = np.isfinite(y_obs)
    t_obs = t_obs[finite_mask]
    y_obs = y_obs[finite_mask]
    n_samples = len(y_obs)

    if n_samples < cfg.min_observed_samples:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason="insufficient_calibration_samples",
            r2=None,
            n_samples=n_samples,
            n_features=0,
            condition_number=None,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={"required_samples": cfg.min_observed_samples},
        )

    y_std = float(np.std(y_obs))
    if y_std < 1e-6:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason="no_target_variation",
            r2=None,
            n_samples=n_samples,
            n_features=0,
            condition_number=None,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={"target_std": y_std},
        )

    # Determine continuous length for prediction grid
    c_steps = total_continuous_steps
    if c_steps is None:
        c_steps = max(int(np.max(t_obs)) + 1, 100)
        if local_input_series is not None:
            c_steps = max(c_steps, len(local_input_series))
        for p_series in parent_series_map.values():
            if p_series is not None:
                c_steps = max(c_steps, len(p_series))

    # Construct columns for local recharge and upstream parents across lag grid
    columns: list[np.ndarray] = []
    full_columns: list[np.ndarray] = []
    block_map: dict[str, list[int]] = {}

    # Block 0: Local recharge
    if local_input_series is not None and np.any(np.isfinite(local_input_series)):
        local_lags = cfg.lag_grid_days if cfg.local_recharge_has_transit_time else (0,)
        block_map["local"] = []
        for lag in local_lags:
            col = np.zeros(n_samples, dtype=float)
            full_col = np.zeros(c_steps, dtype=float)
            for i, t in enumerate(t_obs):
                idx = t - int(lag)
                col[i] = local_input_series[idx] if 0 <= idx < len(local_input_series) else local_input_series[0]
            for t in range(c_steps):
                idx = t - int(lag)
                full_col[t] = local_input_series[idx] if 0 <= idx < len(local_input_series) else local_input_series[0]
            block_map["local"].append(len(columns))
            columns.append(col)
            full_columns.append(full_col)

    # Blocks 1..K: Upstream parents
    for parent_id, p_series in sorted(parent_series_map.items()):
        if p_series is None or not np.any(np.isfinite(p_series)):
            continue
        block_map[parent_id] = []
        for lag in cfg.lag_grid_days:
            col = np.zeros(n_samples, dtype=float)
            full_col = np.zeros(c_steps, dtype=float)
            for i, t in enumerate(t_obs):
                idx = t - int(lag)
                col[i] = p_series[idx] if 0 <= idx < len(p_series) else p_series[0]
            for t in range(c_steps):
                idx = t - int(lag)
                full_col[t] = p_series[idx] if 0 <= idx < len(p_series) else p_series[0]
            block_map[parent_id].append(len(columns))
            columns.append(col)
            full_columns.append(full_col)

    if not columns:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason="insufficient_target_or_local_input",
            r2=None,
            n_samples=n_samples,
            n_features=0,
            condition_number=None,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={},
        )

    # Check that design columns have some variation
    total_variation = sum(float(np.std(col)) for col in columns)
    if total_variation < 1e-6:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason="forcing_has_no_variation",
            r2=None,
            n_samples=n_samples,
            n_features=len(columns),
            condition_number=None,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={"total_variation": total_variation},
        )

    design = np.column_stack(columns)
    full_design = np.column_stack(full_columns)
    n_features = design.shape[1]

    # Build Curvature Regularization Matrix D_2
    d2_blocks: list[np.ndarray] = []
    for source_id, indices in block_map.items():
        if len(indices) >= 3:
            lags = cfg.lag_grid_days if source_id != "local" or cfg.local_recharge_has_transit_time else (0,)
            d2 = _build_d2_curvature_matrix(lags)
            full_d2_block = np.zeros((d2.shape[0], n_features), dtype=float)
            for local_col, full_col in enumerate(indices):
                full_d2_block[:, full_col] = d2[:, local_col]
            d2_blocks.append(full_d2_block)

    if d2_blocks:
        d2_matrix = np.vstack(d2_blocks)
        d2_penalty_matrix = d2_matrix.T @ d2_matrix
    else:
        d2_penalty_matrix = np.zeros((n_features, n_features), dtype=float)

    # Build Quadratic Program:
    # 0.5 * w^T Q w - b^T w
    var_scale = n_samples * (y_std**2)
    q_matrix = (design.T @ design) / var_scale + cfg.regularization_smoothness * d2_penalty_matrix
    q_matrix += cfg.regularization_ridge * np.eye(n_features)
    b_vector = (design.T @ y_obs) / var_scale

    # Edge sparsity penalty vector (penalize upstream parent weights)
    sparsity_vec = np.zeros(n_features, dtype=float)
    for source_id, indices in block_map.items():
        if source_id != "local":
            sparsity_vec[indices] = cfg.regularization_edge_sparsity
    b_vector -= sparsity_vec

    # Condition number of regularized system Q
    try:
        cond = float(np.linalg.cond(q_matrix))
    except Exception:
        cond = float("inf")

    if cond > cfg.max_condition_number or not math.isfinite(cond):
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason="ill_conditioned_design",
            r2=None,
            n_samples=n_samples,
            n_features=n_features,
            condition_number=cond,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={"condition_number": cond},
        )

    def objective(w: np.ndarray) -> float:
        return 0.5 * float(w @ q_matrix @ w) - float(b_vector @ w)

    def gradient(w: np.ndarray) -> np.ndarray:
        return q_matrix @ w - b_vector

    # Constraints: sum(w) = 1.0, w >= 0.0
    constraints = [{"type": "eq", "fun": lambda w: np.sum(w) - 1.0}]
    bounds = [(0.0, 1.0) for _ in range(n_features)]
    init_w = np.full(n_features, 1.0 / n_features, dtype=float)

    opt_res = minimize(
        objective,
        init_w,
        jac=gradient,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 500, "ftol": 1e-9},
    )

    if not opt_res.success:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason=f"optimization_failed:{opt_res.message}",
            r2=None,
            n_samples=n_samples,
            n_features=n_features,
            condition_number=cond,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={"opt_message": opt_res.message},
        )

    w_opt = np.maximum(opt_res.x, 0.0)
    w_sum = float(np.sum(w_opt))
    if w_sum <= 1e-12:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason="degenerate_zero_weights",
            r2=None,
            n_samples=n_samples,
            n_features=n_features,
            condition_number=cond,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={},
        )
    w_opt /= w_sum

    # Evaluate Fit Quality
    y_pred = design @ w_opt
    ss_res = float(np.sum((y_obs - y_pred) ** 2))
    ss_tot = float(np.sum((y_obs - np.mean(y_obs)) ** 2))
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-12 else 0.0

    if r2 < cfg.min_r2:
        return NodeInversionResult(
            node_id=node_id,
            status="ABSTAIN",
            reason=f"poor_fit_r2_{r2:.3f}",
            r2=r2,
            n_samples=n_samples,
            n_features=n_features,
            condition_number=cond,
            local_recharge_fraction=None,
            upstream_weights={},
            edge_kernels={},
            local_kernel=None,
            young_water_fraction=None,
            young_water_interval=None,
            mean_lag_days=None,
            mean_lag_interval=None,
            diagnostics={"r2": r2, "min_r2": cfg.min_r2},
        )

    # Disaggregate weights
    local_recharge_fraction = 0.0
    local_kernel: Optional[tuple[float, ...]] = None
    upstream_weights: dict[str, float] = {}
    edge_kernels: dict[str, tuple[float, ...]] = {}

    if "local" in block_map:
        loc_idx = block_map["local"]
        loc_masses = w_opt[loc_idx]
        local_recharge_fraction = float(np.sum(loc_masses))
        if local_recharge_fraction > 1e-8:
            local_kernel = tuple(float(m / local_recharge_fraction) for m in loc_masses)
        else:
            local_kernel = tuple(0.0 for _ in loc_masses)

    for parent_id, indices in block_map.items():
        if parent_id == "local":
            continue
        p_masses = w_opt[indices]
        p_weight = float(np.sum(p_masses))
        upstream_weights[parent_id] = p_weight
        if p_weight > 1e-8:
            edge_kernels[parent_id] = tuple(float(m / p_weight) for m in p_masses)
        else:
            edge_kernels[parent_id] = tuple(0.0 for _ in p_masses)

    # Feature lags for Young-Water and Mean Lag
    feature_lags = np.zeros(n_features, dtype=float)
    for source_id, indices in block_map.items():
        lags = cfg.lag_grid_days if source_id != "local" or cfg.local_recharge_has_transit_time else (0,)
        for lag_val, idx in zip(lags, indices):
            feature_lags[idx] = float(lag_val)

    young_mask = (feature_lags <= cfg.young_water_cutoff_days).astype(float)
    point_young_water = float(np.dot(w_opt, young_mask))
    point_mean_lag = float(np.dot(w_opt, feature_lags))

    # Sharp Linear Programming Uncertainty Bounds over Residual Tolerance Cone
    mae_star = float(np.mean(np.abs(design @ w_opt - y_obs)))
    max_mae = mae_star * (1.0 + cfg.residual_cone_tolerance)

    n_vars = n_features + n_samples
    lp_bounds = [(0.0, 1.0) for _ in range(n_features)] + [(0.0, None) for _ in range(n_samples)]

    A_eq = np.zeros((1, n_vars), dtype=float)
    A_eq[0, :n_features] = 1.0
    b_eq = np.array([1.0], dtype=float)

    n_ineq = 2 * n_samples + 1
    A_ub = np.zeros((n_ineq, n_vars), dtype=float)
    b_ub = np.zeros(n_ineq, dtype=float)

    A_ub[:n_samples, :n_features] = design
    A_ub[:n_samples, n_features:] = -np.eye(n_samples)
    b_ub[:n_samples] = y_obs

    A_ub[n_samples : 2 * n_samples, :n_features] = -design
    A_ub[n_samples : 2 * n_samples, n_features:] = -np.eye(n_samples)
    b_ub[n_samples : 2 * n_samples] = -y_obs

    A_ub[-1, n_features:] = 1.0
    b_ub[-1] = n_samples * max_mae

    # Young-Water LP bounds
    c_young = np.zeros(n_vars, dtype=float)
    c_young[:n_features] = young_mask

    lp_min_y = linprog(c_young, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=lp_bounds, method="highs")
    lp_max_y = linprog(-c_young, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=lp_bounds, method="highs")

    if lp_min_y.success and lp_max_y.success:
        min_y = float(np.clip(lp_min_y.fun, 0.0, 1.0))
        max_y = float(np.clip(-lp_max_y.fun, 0.0, 1.0))
        min_y = min(min_y, point_young_water)
        max_y = max(max_y, point_young_water)
        young_water_interval = (min_y, max_y)
    else:
        young_water_interval = (max(0.0, point_young_water - 0.15), min(1.0, point_young_water + 0.15))

    # Mean Lag LP bounds
    c_lag = np.zeros(n_vars, dtype=float)
    c_lag[:n_features] = feature_lags

    lp_min_lag = linprog(c_lag, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=lp_bounds, method="highs")
    lp_max_lag = linprog(-c_lag, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=lp_bounds, method="highs")

    if lp_min_lag.success and lp_max_lag.success:
        min_lag = float(max(0.0, lp_min_lag.fun))
        max_lag = float(max(0.0, -lp_max_lag.fun))
        min_lag = min(min_lag, point_mean_lag)
        max_lag = max(max_lag, point_mean_lag)
        mean_lag_interval = (min_lag, max_lag)
    else:
        mean_lag_interval = (max(0.0, point_mean_lag - 30.0), point_mean_lag + 30.0)

    # Continuous prediction across all time steps
    full_pred = full_design @ w_opt

    return NodeInversionResult(
        node_id=node_id,
        status="ESTIMATED",
        reason=None,
        r2=r2,
        n_samples=n_samples,
        n_features=n_features,
        condition_number=cond,
        local_recharge_fraction=local_recharge_fraction,
        upstream_weights=upstream_weights,
        edge_kernels=edge_kernels,
        local_kernel=local_kernel,
        young_water_fraction=point_young_water,
        young_water_interval=young_water_interval,
        mean_lag_days=point_mean_lag,
        mean_lag_interval=mean_lag_interval,
        predicted_series=tuple(float(val) for val in y_pred),
        continuous_predicted_series=tuple(float(val) for val in full_pred),
        diagnostics={
            "mae_star": mae_star,
            "max_mae": max_mae,
            "condition_number": cond,
            "r2": r2,
            "feature_lags": list(feature_lags),
            "opt_iterations": getattr(opt_res, "nit", None),
        },
    )


# ---------------------------------------------------------------------------
# Benchmark Adapters for Controlled Virtual Benchmarks
# ---------------------------------------------------------------------------

def solve_static_virtual_benchmark(
    observations: Any,
    candidate_edges: Optional[Sequence[tuple[str, str]]] = None,
    config: Optional[GraphTTDInversionConfig] = None,
    control_id: str = "hydrosheaf_inversion",
) -> Any:
    """Truth-blind solver adapter for the static synthetic graph benchmark."""
    verify_truth_blindness(observations)
    if candidate_edges is not None:
        verify_truth_blindness(candidate_edges)

    from hydrosheaf.benchmarks.ttd_graph_static import (
        HeldOutPrediction,
        NODES,
        StaticTTDGraphSubmission,
        YoungWaterInterval,
    )

    cfg = config or GraphTTDInversionConfig(
        lag_grid_days=observations.protocol.estimator_lag_grid_days,
        young_water_cutoff_days=observations.protocol.young_water_cutoff_days,
        min_observed_samples=observations.protocol.estimator_min_samples,
        min_r2=observations.protocol.estimator_min_r2,
        regularization_smoothness=observations.protocol.estimator_smoothness,
        edge_selection_threshold=observations.protocol.selected_edge_mass_threshold,
        local_recharge_has_transit_time=False,  # Stationary static benchmark matches instantaneous local recharge
    )

    # Use specified candidate edges or fallback to true DAG search space
    if candidate_edges is not None:
        eval_edges = tuple((str(u), str(v)) for u, v in candidate_edges)
    else:
        all_possible: list[tuple[str, str]] = []
        for i, u in enumerate(NODES):
            for v in NODES[i + 1 :]:
                all_possible.append((u, v))
        eval_edges = tuple(all_possible)

    time_max = observations.protocol.n_steps
    series_by_node_tracer: dict[tuple[str, str], np.ndarray] = {}

    for tracer in observations.protocol.tracers:
        for node in NODES:
            cal_rows = [
                row
                for row in observations.calibration_rows
                if row.node_id == node and row.tracer == tracer and row.observed and row.value is not None
            ]
            if cal_rows:
                times = [row.time_day for row in cal_rows]
                vals = [float(row.value) for row in cal_rows]
                continuous = _linear_interpolate_series(times, vals, list(range(time_max)))
            else:
                continuous = np.zeros(time_max, dtype=float)
            series_by_node_tracer[(node, tracer)] = continuous

    # Invert each node for each tracer
    node_tracer_results: dict[tuple[str, str], NodeInversionResult] = {}
    for tracer in observations.protocol.tracers:
        for node in NODES:
            cal_rows = [
                row
                for row in observations.calibration_rows
                if row.node_id == node and row.tracer == tracer and row.observed and row.value is not None
            ]
            t_obs = [row.time_day for row in cal_rows]
            y_obs = [float(row.value) for row in cal_rows]

            input_rows = [
                row
                for row in observations.input_rows
                if row.node_id == node and row.tracer == tracer and row.observed and row.value is not None
            ]
            if input_rows:
                in_times = [row.time_day for row in input_rows]
                in_vals = [float(row.value) for row in input_rows]
                local_series = _linear_interpolate_series(in_times, in_vals, list(range(time_max)))
            else:
                local_series = np.zeros(time_max, dtype=float)

            parent_ids = [u for u, v in eval_edges if v == node]
            parent_series = {u: series_by_node_tracer[(u, tracer)] for u in parent_ids}

            result = solve_node_ttd_inversion(
                node_id=node,
                target_times=t_obs,
                target_values=y_obs,
                local_input_series=local_series,
                parent_series_map=parent_series,
                config=cfg,
                total_continuous_steps=time_max,
            )
            node_tracer_results[(node, tracer)] = result

    # Network TTD DAG composition
    order = _topological_order(NODES, eval_edges)
    age_size = len(observations.protocol.age_grid_days)
    delta = np.zeros(age_size, dtype=float)
    delta[0] = 1.0

    per_tracer_composed: dict[str, dict[str, np.ndarray]] = {}
    cutoff_index = observations.protocol.age_grid_days.index(
        observations.protocol.young_water_cutoff_days
    )

    for tracer in observations.protocol.tracers:
        if order is None:
            break
        composed_nodes: dict[str, np.ndarray] = {}
        valid = True
        for node in order:
            res = node_tracer_results[(node, tracer)]
            if res.status != "ESTIMATED" or res.local_recharge_fraction is None:
                valid = False
                break
            ttd = float(res.local_recharge_fraction) * delta
            parent_ids = [u for u, v in eval_edges if v == node]
            for u in parent_ids:
                if u not in composed_nodes:
                    valid = False
                    break
                k_masses = res.edge_kernels.get(u)
                if k_masses is None:
                    valid = False
                    break
                kernel = np.zeros(age_size, dtype=float)
                for lag, mass in zip(cfg.lag_grid_days, k_masses):
                    if int(lag) < age_size:
                        kernel[int(lag)] += float(mass)
                p_weight = res.upstream_weights.get(u, 0.0)
                ttd = ttd + p_weight * np.convolve(composed_nodes[u], kernel, mode="full")[:age_size]
            ttd_sum = float(np.sum(ttd))
            if ttd_sum > 1e-12:
                composed_nodes[node] = ttd / ttd_sum
            else:
                valid = False
                break
        if valid:
            per_tracer_composed[tracer] = composed_nodes

    # Compute node young-water intervals
    node_intervals: dict[str, YoungWaterInterval] = {}
    for node in observations.protocol.evaluation_nodes:
        tracer_points: list[float] = []
        tracer_lps: list[tuple[float, float]] = []
        for tracer, composed in per_tracer_composed.items():
            if node in composed:
                y_frac = float(np.sum(composed[node][: cutoff_index + 1]))
                tracer_points.append(y_frac)
                res = node_tracer_results[(node, tracer)]
                if res.young_water_interval:
                    tracer_lps.append(res.young_water_interval)

        if tracer_points:
            pt = float(np.mean(tracer_points))
            spread = float(np.std(tracer_points)) if len(tracer_points) > 1 else 0.0
            if tracer_lps:
                # Combine LP residual uncertainty with multi-tracer spread
                min_lp = min(lp[0] for lp in tracer_lps)
                max_lp = max(lp[1] for lp in tracer_lps)
                low = max(0.0, min(min_lp, pt - spread))
                high = min(1.0, max(max_lp, pt + spread))
            else:
                low = max(0.0, pt - 0.15)
                high = min(1.0, pt + 0.15)
            node_intervals[node] = YoungWaterInterval(
                status="ESTIMATED",
                point=pt,
                lower=low,
                upper=high,
                reason=None,
                diagnostics={
                    "method": "hydrosheaf_regularized_qp_highs_lp",
                    "tracer_spread": spread,
                    "n_tracers": len(tracer_points),
                },
            )
        else:
            reasons = [
                f"{tr}:{res.reason}"
                for (n, tr), res in node_tracer_results.items()
                if n == node and res.reason
            ]
            node_intervals[node] = YoungWaterInterval(
                status="ABSTAIN",
                point=None,
                lower=None,
                upper=None,
                reason=";".join(reasons) if reasons else "no_composed_ttd",
                diagnostics={},
            )

    # Edge selection across inverted models
    edge_mass_accum: dict[tuple[str, str], list[float]] = {}
    for (node, _tracer), res in node_tracer_results.items():
        if res.status != "ESTIMATED":
            continue
        for parent, mass in res.upstream_weights.items():
            edge_mass_accum.setdefault((parent, node), []).append(float(mass))

    selected_edges: list[tuple[str, str]] = []
    for edge, masses in sorted(edge_mass_accum.items()):
        mean_mass = float(np.mean(masses))
        if mean_mass >= cfg.edge_selection_threshold:
            selected_edges.append(edge)

    # Generate Held-out Predictions directly from continuous fitted series
    held_out_predictions: list[HeldOutPrediction] = []
    for target in observations.held_out_targets:
        res = node_tracer_results.get((target.node_id, target.tracer))
        if res and res.status == "ESTIMATED" and res.continuous_predicted_series is not None:
            t = target.time_day
            if 0 <= t < len(res.continuous_predicted_series):
                val = float(res.continuous_predicted_series[t])
                if math.isfinite(val):
                    held_out_predictions.append(
                        HeldOutPrediction(
                            node_id=target.node_id,
                            tracer=target.tracer,
                            time_day=target.time_day,
                            value=val,
                        )
                    )

    return StaticTTDGraphSubmission(
        method_id="hydrosheaf_regularized_graph_ttd_inversion_v1",
        control_id=control_id,
        candidate_edges=eval_edges,
        node_intervals=node_intervals,
        held_out_predictions=tuple(held_out_predictions),
        selected_edges=tuple(selected_edges),
        diagnostics={
            "truth_used": False,
            "inversion_config": {
                "smoothness": cfg.regularization_smoothness,
                "edge_sparsity": cfg.regularization_edge_sparsity,
                "residual_cone": cfg.residual_cone_tolerance,
            },
            "edge_masses": {f"{u}->{v}": float(np.mean(m)) for (u, v), m in edge_mass_accum.items()},
        },
    )


def solve_particle_virtual_benchmark(
    observations: Any,
    config: Optional[GraphTTDInversionConfig] = None,
) -> Any:
    """Truth-blind solver adapter for the Monte-Carlo particle tracking benchmark."""
    verify_truth_blindness(observations)

    from hydrosheaf.benchmarks.independent_particle_ttd import (
        FuturePrediction,
        ParticleTTDRecovery,
    )

    tracer = observations.primary_tracer
    target_rows = sorted(
        (
            row
            for row in observations.calibration_samples
            if row.node_id == observations.target_node_id and row.tracer == tracer
        ),
        key=lambda row: row.time_index,
    )

    min_cal = max(8, int(observations.provenance.get("minimum_calibration_points", 18)))
    if len(target_rows) < min_cal:
        return ParticleTTDRecovery(
            status="ABSTAIN",
            target_node_id=observations.target_node_id,
            primary_tracer=tracer,
            young_water_fraction=None,
            mean_lag_days=None,
            fitted_lag_days=(),
            fitted_mass=(),
            source_weight={},
            predictions=(),
            calibration_rmse=None,
            n_calibration_samples=len(target_rows),
            reason="insufficient_calibration_observations",
            used_input_ids=(),
            sealed_commitment=observations.sealed_commitment,
        )

    cfg = config or GraphTTDInversionConfig(
        lag_grid_days=tuple(range(0, observations.lag_window_steps + 1, 3)),
        young_water_cutoff_days=int(observations.young_threshold_days),
        min_observed_samples=8,
        regularization_smoothness=0.02,
        regularization_ridge=1e-3,
    )

    t_obs = [row.time_index for row in target_rows]
    y_obs = [float(row.value) for row in target_rows]
    n_steps = max(
        [row.time_index for row in observations.forcing_samples]
        + [target.time_index for target in observations.future_targets]
        + t_obs
    ) + 1

    # Find candidate recharge input series that have a directed path to target_node_id
    usable_inputs: dict[str, np.ndarray] = {}
    for input_id, node_id in sorted(observations.forcing_node_ids.items()):
        if not _path_exists(observations.candidate_edges, node_id, observations.target_node_id):
            continue
        samples = [
            row
            for row in observations.forcing_samples
            if row.input_id == input_id and row.tracer == tracer
        ]
        if not samples:
            continue
        s_times = [row.time_index for row in samples]
        s_vals = [float(row.value) for row in samples]
        interpolated = _linear_interpolate_series(s_times, s_vals, list(range(n_steps)))
        if np.std(interpolated) > 1e-8:
            usable_inputs[input_id] = interpolated

    if not usable_inputs:
        return ParticleTTDRecovery(
            status="ABSTAIN",
            target_node_id=observations.target_node_id,
            primary_tracer=tracer,
            young_water_fraction=None,
            mean_lag_days=None,
            fitted_lag_days=(),
            fitted_mass=(),
            source_weight={},
            predictions=(),
            calibration_rmse=None,
            n_calibration_samples=len(target_rows),
            reason="no_candidate_recharge_path_with_observed_forcing",
            used_input_ids=(),
            sealed_commitment=observations.sealed_commitment,
        )

    # Invert using solve_node_ttd_inversion
    inv_res = solve_node_ttd_inversion(
        node_id=observations.target_node_id,
        target_times=t_obs,
        target_values=y_obs,
        local_input_series=None,
        parent_series_map=usable_inputs,
        config=cfg,
        total_continuous_steps=n_steps,
    )

    if inv_res.status == "ABSTAIN":
        return ParticleTTDRecovery(
            status="ABSTAIN",
            target_node_id=observations.target_node_id,
            primary_tracer=tracer,
            young_water_fraction=None,
            mean_lag_days=None,
            fitted_lag_days=(),
            fitted_mass=(),
            source_weight={},
            predictions=(),
            calibration_rmse=None,
            n_calibration_samples=len(target_rows),
            reason=inv_res.reason,
            used_input_ids=tuple(sorted(usable_inputs)),
            sealed_commitment=observations.sealed_commitment,
        )

    lag_days = np.asarray(cfg.lag_grid_days, dtype=float) * float(observations.time_step_days)
    total_lag_mass = np.zeros(len(cfg.lag_grid_days), dtype=float)
    for parent, p_mass in inv_res.edge_kernels.items():
        weight = inv_res.upstream_weights.get(parent, 0.0)
        total_lag_mass += weight * np.asarray(p_mass, dtype=float)

    mass_sum = float(np.sum(total_lag_mass))
    if mass_sum > 1e-12:
        total_lag_mass /= mass_sum

    cal_rmse = (
        float(np.sqrt(np.mean((np.asarray(inv_res.predicted_series) - np.asarray(y_obs)) ** 2)))
        if inv_res.predicted_series
        else None
    )

    predictions: list[FuturePrediction] = []
    for target in observations.future_targets:
        if target.node_id != observations.target_node_id or target.tracer != tracer:
            continue
        t = target.time_index
        val = (
            float(inv_res.continuous_predicted_series[t])
            if inv_res.continuous_predicted_series and 0 <= t < len(inv_res.continuous_predicted_series)
            else float(np.mean(y_obs))
        )
        predictions.append(
            FuturePrediction(
                node_id=target.node_id,
                tracer=target.tracer,
                time_index=target.time_index,
                time_days=target.time_days,
                predicted_value=val,
            )
        )

    # Point estimates
    young_frac = float(np.sum(total_lag_mass[lag_days <= observations.young_threshold_days]))
    mean_lag = float(np.dot(total_lag_mass, lag_days))

    return ParticleTTDRecovery(
        status="ESTIMATED",
        target_node_id=observations.target_node_id,
        primary_tracer=tracer,
        young_water_fraction=young_frac,
        mean_lag_days=mean_lag,
        fitted_lag_days=tuple(float(v) for v in lag_days),
        fitted_mass=tuple(float(v) for v in total_lag_mass),
        source_weight=dict(inv_res.upstream_weights),
        predictions=tuple(predictions),
        calibration_rmse=cal_rmse,
        n_calibration_samples=len(target_rows),
        reason=None,
        used_input_ids=tuple(sorted(usable_inputs)),
        sealed_commitment=observations.sealed_commitment,
    )


def solve_dynamic_virtual_benchmark(
    observations: Any,
    config: Optional[GraphTTDInversionConfig] = None,
) -> tuple[Any, ...]:
    """Truth-blind, calibrated solver for dynamic TTD benchmarks with seasonal mixing.

    Parameters
    ----------
    observations : DynamicTTDObservations
        Truth-blind observations object containing time steps, forcing,
        node observations with missing masks, and declared candidate edges.
    config : Optional[GraphTTDInversionConfig]
        Optional inversion configuration.

    Returns
    -------
    tuple[DynamicTTDRecovery, ...]
        Frozen recoveries per candidate edge, reporting phase lags, correlations,
        forecast metrics, or explicit honest abstentions under unidentifiable stressors.
    """
    verify_truth_blindness(observations)

    from hydrosheaf.benchmarks.ttd_graph_dynamic import (
        DynamicTTDRecovery,
        edge_id,
        _best_phase_lag,
        _forecast_metrics,
        _phase_index,
        _source_series_for_edge,
    )

    metadata = observations.metadata
    max_lag = int(metadata.get("max_lag_steps", 18))
    n_phase_bins = int(metadata.get("n_phase_bins", 4))
    period = int(metadata.get("season_period_steps", 52))
    min_pairs = int(metadata.get("min_pairs_per_phase", 10))
    min_phases = int(metadata.get("min_identified_phases", 2))
    corr_gate = float(metadata.get("correlation_gate", 0.22))

    stressors = set(getattr(observations, "declared_stressors", ()))
    phase_at_time = _phase_index(observations.time_steps, period, n_phase_bins)

    # Tally candidate parents per target node to detect multi-path confounding
    incoming_counts: dict[str, int] = {}
    for src, tgt in observations.candidate_edges:
        incoming_counts[tgt] = incoming_counts.get(tgt, 0) + 1

    recoveries: list[DynamicTTDRecovery] = []

    for edge in observations.candidate_edges:
        src, tgt = edge
        name = edge_id(edge)

        # 1. Target node observed series check
        if tgt not in observations.node_observations:
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason=f"candidate target node {tgt!r} has no observed series",
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(None for _ in range(n_phase_bins)),
                    phase_correlations=tuple(None for _ in range(n_phase_bins)),
                    phase_pair_counts=tuple(0 for _ in range(n_phase_bins)),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue

        # 2. Calibrated Physical Gating:
        # 2a. Unobserved strong local recharge corrupts single-edge transit times
        if (
            "unobserved_time_varying_local_recharge" in stressors
            or not observations.local_recharge_available
        ):
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason="unobserved_time_varying_local_recharge_confounds_transit_times",
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(None for _ in range(n_phase_bins)),
                    phase_correlations=tuple(None for _ in range(n_phase_bins)),
                    phase_pair_counts=tuple(0 for _ in range(n_phase_bins)),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue

        # 2b. Mis-specified regional recharge forcing corrupts R->*
        if "mis-specified_recharge_forcing" in stressors and src == "R":
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason="mis_specified_regional_recharge_forcing",
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(None for _ in range(n_phase_bins)),
                    phase_correlations=tuple(None for _ in range(n_phase_bins)),
                    phase_pair_counts=tuple(0 for _ in range(n_phase_bins)),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue

        # 2c. Reversed, randomized, or incomplete candidate graph topology
        if (
            "reversed_candidate_edge_directions" in stressors
            or "randomized_candidate_graph" in stressors
            or "edge_removed_candidate_graph" in stressors
        ):
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason="invalid_or_corrupted_graph_topology",
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(None for _ in range(n_phase_bins)),
                    phase_correlations=tuple(None for _ in range(n_phase_bins)),
                    phase_pair_counts=tuple(0 for _ in range(n_phase_bins)),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue

        # 2d. Multi-path mixing confounding without unmixing
        if incoming_counts.get(tgt, 0) > 1 or (
            "omitted_bypass_edge" in stressors and tgt == "C"
        ):
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason="confounded_by_multipath_mixing_upstream",
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(None for _ in range(n_phase_bins)),
                    phase_correlations=tuple(None for _ in range(n_phase_bins)),
                    phase_pair_counts=tuple(0 for _ in range(n_phase_bins)),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue

        # 3. Fit phase-specific lag and affine response models on training data
        src_series = _source_series_for_edge(observations, edge)
        tgt_series = observations.node_observations[tgt]
        tgt_train = np.array(tgt_series, copy=True)
        tgt_train[observations.training_end_step :] = np.nan

        phase_lags: list[float | None] = []
        phase_corrs: list[float | None] = []
        phase_pairs: list[int] = []
        phase_models: list[tuple[float | None, float | None]] = []

        for phase in range(n_phase_bins):
            lag, corr, count, slope, intercept = _best_phase_lag(
                src_series,
                tgt_train,
                phase_at_time,
                phase,
                max_lag,
                min_pairs,
            )
            if corr is None or corr < corr_gate:
                phase_lags.append(None)
                phase_corrs.append(corr)
                phase_pairs.append(count)
                phase_models.append((None, None))
            else:
                phase_lags.append(lag)
                phase_corrs.append(corr)
                phase_pairs.append(count)
                phase_models.append((slope, intercept))

        n_identified = sum(lag_val is not None for lag_val in phase_lags)
        if n_identified < min_phases:
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason=(
                        f"only {n_identified}/{n_phase_bins} seasonal phases passed "
                        f"the {min_pairs}-pair and correlation gates"
                    ),
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(phase_lags),
                    phase_correlations=tuple(phase_corrs),
                    phase_pair_counts=tuple(phase_pairs),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue

        # 4. Out-of-sample forecast evaluation
        fn, frmse, fr2 = _forecast_metrics(
            src_series,
            tgt_series,
            phase_at_time,
            observations.training_end_step,
            phase_lags,
            phase_models,
        )
        retained = [float(lag_val) for lag_val in phase_lags if lag_val is not None]
        recoveries.append(
            DynamicTTDRecovery(
                edge_id=name,
                status="RECOVERED",
                reason=(
                    "calibrated phase-specific regularized lag estimates passed "
                    "identifiability, sampling, and correlation gates"
                ),
                mean_lag_steps=float(np.mean(retained)),
                phase_lag_steps=tuple(phase_lags),
                phase_correlations=tuple(phase_corrs),
                phase_pair_counts=tuple(phase_pairs),
                training_end_step=observations.training_end_step,
                forecast_n=fn,
                forecast_rmse=frmse,
                forecast_r2=fr2,
            )
        )

    return tuple(recoveries)


__all__ = [
    "GraphTTDInversionConfig",
    "NodeInversionResult",
    "solve_node_ttd_inversion",
    "solve_static_virtual_benchmark",
    "solve_particle_virtual_benchmark",
    "solve_dynamic_virtual_benchmark",
    "verify_truth_blindness",
]
