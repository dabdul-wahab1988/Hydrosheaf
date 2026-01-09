"""
Temporal edge fitting - fit transport and reaction models across time series.
"""

from typing import List, Optional

import numpy as np

from . import TemporalEdgeResult, TemporalNode
from .interpolation import align_time_series
from .residence_time import estimate_residence_time


def fit_temporal_edge(
    node_u: TemporalNode,
    node_v: TemporalNode,
    config: "Config",  # type: ignore
    edge_id: str = "",
) -> TemporalEdgeResult:
    """
    Fit transport and reaction model across all time points.

    Parameters
    ----------
    node_u, node_v : TemporalNode
        Upstream and downstream nodes (must have same time grid)
    config : Config
        Hydrosheaf configuration
    edge_id : str
        Identifier for this edge

    Returns
    -------
    TemporalEdgeResult
        Contains time-averaged parameters and per-time extents

    Mathematical Implementation
    ---------------------------
    1. Estimate residence time τ via estimate_residence_time()

    2. Align time series: for each v sample at time t, find u sample at time t - τ
       (use interpolation if exact match not available)

    3. For each aligned time point k:
        a. Fit transport model (evap or mix) using existing fit_evaporation/fit_mixing
           γ_k = argmin_γ ||x_v(t_k) - γ * x_u(t_k - τ)||²_W

        b. Compute post-transport residual:
           r_k = x_v(t_k) - T(γ_k) * x_u(t_k - τ)

        c. Fit reactions using existing fit_reactions:
           ξ_k = argmin_ξ ||r_k - R*ξ||²_W + λ||ξ||_1

    4. Aggregate across time:
        γ_mean = (1/N) Σ_k γ_k
        γ_std = sqrt((1/N) Σ_k (γ_k - γ_mean)²)
        ξ_mean[j] = (1/N) Σ_k ξ_k[j]
        ξ_std[j] = sqrt((1/N) Σ_k (ξ_k[j] - ξ_mean[j])²)

    5. Compute total residual:
        L_total = Σ_k ||x_v(t_k) - T(γ_k)*x_u(t_k-τ) - R*ξ_k||²_W
    """
    from ..inference.edge_fit import fit_edge

    # Estimate residence time
    tau_days, tau_uncertainty, tau_method = estimate_residence_time(
        node_u,
        node_v,
        method=getattr(config, "residence_time_method", "cross_correlation"),
        tracer_ion=getattr(config, "residence_time_tracer", "Cl"),
        ion_order=config.ion_order,
    )

    # Align time series with lag
    aligned_u, aligned_v = align_time_series(node_u, node_v, tau_days)

    if not aligned_u or not aligned_v:
        # No aligned data, return empty result
        return TemporalEdgeResult(
            edge_id=edge_id,
            u=node_u.node_id,
            v=node_v.node_id,
            residence_time_days=tau_days,
            residence_time_method=tau_method,
            residence_time_uncertainty=tau_uncertainty,
            transport_model="none",
        )

    # Fit each time point
    gamma_series = []
    f_series = []
    extents_series = []
    residual_series = []
    timestamps = []

    for u_sample, v_sample in zip(aligned_u, aligned_v):
        # Fit edge for this time point
        try:
            result = fit_edge(
                u_sample.concentrations,
                v_sample.concentrations,
                config,
                edge_id=f"{edge_id}_t{len(timestamps)}",
            )

            # Store results
            if result.transport_model == "evap" and result.gamma is not None:
                gamma_series.append(result.gamma)
            elif result.transport_model == "mix" and result.f is not None:
                f_series.append(result.f)

            if result.z_extents:
                extents_series.append(result.z_extents)
            else:
                extents_series.append([])

            residual_series.append(result.transport_residual_norm)
            timestamps.append(v_sample.timestamp)

        except Exception:
            # Skip failed fits
            continue

    # Aggregate statistics
    result = TemporalEdgeResult(
        edge_id=edge_id,
        u=node_u.node_id,
        v=node_v.node_id,
        residence_time_days=tau_days,
        residence_time_method=tau_method,
        residence_time_uncertainty=tau_uncertainty,
        transport_model="evap" if gamma_series else ("mix" if f_series else "none"),
        timestamps=timestamps,
        reaction_extents_series=extents_series,
        per_time_residual=residual_series,
    )

    # Compute mean and std for gamma
    if gamma_series:
        gamma_arr = np.array(gamma_series)
        result.gamma_mean = float(np.mean(gamma_arr))
        result.gamma_std = float(np.std(gamma_arr, ddof=1) if len(gamma_arr) > 1 else 0.0)

    # Compute mean and std for f
    if f_series:
        f_arr = np.array(f_series)
        result.f_mean = float(np.mean(f_arr))
        result.f_std = float(np.std(f_arr, ddof=1) if len(f_arr) > 1 else 0.0)

    # Compute mean and std for extents
    if extents_series and extents_series[0]:
        n_reactions = len(extents_series[0])
        extents_arr = np.array([ext if ext else [0.0] * n_reactions for ext in extents_series])

        result.reaction_extents_mean = [float(np.mean(extents_arr[:, j])) for j in range(n_reactions)]
        result.reaction_extents_std = [
            float(np.std(extents_arr[:, j], ddof=1) if len(extents_arr) > 1 else 0.0) for j in range(n_reactions)
        ]

    # Total residual
    result.total_residual_norm = float(np.sum(residual_series)) if residual_series else 0.0

    return result


def compute_seasonal_decomposition(
    node: TemporalNode, ion_index: int, period_days: float = 365.0
) -> tuple[float, float, float, float]:
    """
    Decompose concentration into trend + seasonal + residual.

    C(t) = μ + β*t + α*cos(2πt/T) + β*sin(2πt/T) + ε(t)

    Parameters
    ----------
    node : TemporalNode
        Node with time series
    ion_index : int
        Index of ion to analyze
    period_days : float
        Period for seasonal component (default 365 days for annual)

    Returns
    -------
    tuple[float, float, float, float]
        (mean, trend, seasonal_amplitude, residual_std)
    """
    if len(node.samples) < 4:
        return 0.0, 0.0, 0.0, 0.0

    # Extract time series
    t0 = node.samples[0].timestamp
    times = [(s.timestamp - t0).total_seconds() / 86400.0 for s in node.samples]
    values = [s.concentrations[ion_index] for s in node.samples]

    times = np.array(times)
    values = np.array(values)

    # Build design matrix for least squares
    # C = a0 + a1*t + a2*cos(2πt/T) + a3*sin(2πt/T)
    omega = 2 * np.pi / period_days

    X = np.column_stack([np.ones_like(times), times, np.cos(omega * times), np.sin(omega * times)])

    # Solve normal equations
    try:
        coeffs = np.linalg.lstsq(X, values, rcond=None)[0]
    except np.linalg.LinAlgError:
        return float(np.mean(values)), 0.0, 0.0, float(np.std(values))

    mean = coeffs[0]
    trend = coeffs[1]
    seasonal_amplitude = np.sqrt(coeffs[2] ** 2 + coeffs[3] ** 2)

    # Compute residuals
    fitted = X @ coeffs
    residuals = values - fitted
    residual_std = np.std(residuals)

    return float(mean), float(trend), float(seasonal_amplitude), float(residual_std)
