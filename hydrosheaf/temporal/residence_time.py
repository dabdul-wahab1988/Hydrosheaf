"""
Residence time estimation methods.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np

from . import TemporalNode


def estimate_residence_time(
    node_u: TemporalNode,
    node_v: TemporalNode,
    method: str = "cross_correlation",
    tracer_ion: str = "Cl",
    ion_order: Optional[List[str]] = None,
    hydraulic_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float, str]:
    """
    Estimate groundwater travel time between two nodes.

    Parameters
    ----------
    node_u, node_v : TemporalNode
        Upstream and downstream nodes with time-series
    method : str
        "cross_correlation" - lag via signal correlation
        "gradient" - Darcy's law estimate
        "tracer_decay" - radioactive tracer
    tracer_ion : str
        Conservative tracer for correlation (typically "Cl")
    ion_order : Optional[List[str]]
        Ion order to find tracer index
    hydraulic_params : Optional[Dict]
        Required for "gradient" method:
        {"distance_m": float, "K_m_day": float, "gradient": float, "porosity": float}

    Returns
    -------
    Tuple[float, float, str]
        (estimated_tau_days, uncertainty_days, method_used)

    Mathematical Implementation
    ---------------------------
    Cross-correlation:
        1. Extract tracer time-series: u_t = C_u^{tracer}(t), v_t = C_v^{tracer}(t)
        2. Compute normalized cross-correlation:
           r(τ) = Σ_t (u_t - μ_u)(v_{t+τ} - μ_v) / (σ_u σ_v N)
        3. Find τ* = argmax_τ r(τ)
        4. Uncertainty from peak width at r(τ*) - 0.1

    Gradient:
        τ = (distance_m * porosity) / (K_m_day * gradient)

    Tracer decay (e.g., tritium, t_1/2 = 12.32 years):
        τ = -t_1/2 / ln(2) * ln(C_v / C_u)
    """
    if method == "cross_correlation":
        return _estimate_residence_time_cross_correlation(node_u, node_v, tracer_ion, ion_order)
    elif method == "gradient":
        return _estimate_residence_time_gradient(hydraulic_params)
    elif method == "tracer_decay":
        return _estimate_residence_time_tracer_decay(node_u, node_v, tracer_ion, ion_order)
    else:
        raise ValueError(f"Unknown method: {method}")


def _estimate_residence_time_cross_correlation(
    node_u: TemporalNode, node_v: TemporalNode, tracer_ion: str, ion_order: Optional[List[str]]
) -> Tuple[float, float, str]:
    """
    Estimate residence time via cross-correlation of tracer signal.
    """
    if not node_u.samples or not node_v.samples:
        return 0.0, 0.0, "cross_correlation_failed"

    # Find tracer index
    tracer_idx = 4  # Default to Cl (index 4 in standard order)
    if ion_order:
        try:
            tracer_idx = ion_order.index(tracer_ion)
        except ValueError:
            pass  # Use default

    # Extract tracer time series
    u_times = np.array([(s.timestamp - node_u.samples[0].timestamp).total_seconds() / 86400.0 for s in node_u.samples])
    v_times = np.array([(s.timestamp - node_v.samples[0].timestamp).total_seconds() / 86400.0 for s in node_v.samples])

    u_values = np.array([s.concentrations[tracer_idx] for s in node_u.samples])
    v_values = np.array([s.concentrations[tracer_idx] for s in node_v.samples])

    # Need overlapping time range
    if len(u_values) < 3 or len(v_values) < 3:
        return 0.0, 0.0, "cross_correlation_insufficient_data"

    # Normalize signals
    u_mean = np.mean(u_values)
    v_mean = np.mean(v_values)
    u_std = np.std(u_values)
    v_std = np.std(v_values)

    if u_std < 1e-6 or v_std < 1e-6:
        # No variation, can't compute correlation
        return 0.0, 0.0, "cross_correlation_no_variation"

    u_norm = (u_values - u_mean) / u_std
    v_norm = (v_values - v_mean) / v_std

    # Compute cross-correlation for different lags
    # We'll test lags from -max_lag to +max_lag days
    max_lag_days = min(365.0, (v_times[-1] - u_times[0]))  # Max 1 year
    lag_step_days = max(1.0, max_lag_days / 100)  # ~100 points

    lags = np.arange(0, max_lag_days, lag_step_days)
    correlations = []

    for lag in lags:
        # For each v sample, find u sample at time (t_v - lag)
        corr_sum = 0.0
        count = 0

        for i, v_t in enumerate(v_times):
            u_t_target = v_t - lag

            # Find nearest u sample
            if u_t_target < u_times[0] or u_t_target > u_times[-1]:
                continue

            # Linear interpolation
            u_interp = np.interp(u_t_target, u_times, u_norm)
            corr_sum += u_interp * v_norm[i]
            count += 1

        if count > 0:
            correlations.append(corr_sum / count)
        else:
            correlations.append(0.0)

    correlations = np.array(correlations)

    # Find peak
    # Center of Mass Calculation (First Moment)
    # Robust against dispersive smearing
    
    # Filter negative/low correlations to reduce noise
    positive_corr_mask = correlations > 0.1 * np.max(correlations)
    
    if np.sum(positive_corr_mask) < 2:
         # Fallback to peak if distribution is too sharp or noisy to integrate
        max_idx = np.argmax(correlations)
        best_lag = lags[max_idx]
        uncertainty = 10.0 # Default fallback
    else:
        # Weighted mean of lags
        valid_lags = lags[positive_corr_mask]
        valid_corrs = correlations[positive_corr_mask]
        
        sum_weights = np.sum(valid_corrs)
        center_of_mass = np.sum(valid_lags * valid_corrs) / sum_weights
        best_lag = center_of_mass
        
        # Weighted standard deviation for uncertainty
        variance = np.sum(((valid_lags - center_of_mass) ** 2) * valid_corrs) / sum_weights
        uncertainty = np.sqrt(variance)

    return float(best_lag), float(uncertainty), "cross_correlation"


def _estimate_residence_time_gradient(hydraulic_params: Optional[Dict[str, float]]) -> Tuple[float, float, str]:
    """
    Estimate residence time using Darcy's law.

    τ = (distance * porosity) / (K * gradient)
    """
    if not hydraulic_params:
        return 0.0, 0.0, "gradient_no_params"

    distance = hydraulic_params.get("distance_m", 0.0)
    K = hydraulic_params.get("K_m_day", 1.0)
    gradient = hydraulic_params.get("gradient", 0.001)
    porosity = hydraulic_params.get("porosity", 0.2)

    if K <= 0 or gradient <= 0 or distance <= 0:
        return 0.0, 0.0, "gradient_invalid_params"

    # Darcy velocity: v = K * i
    # Pore velocity: v_p = v / n_e
    # Travel time: τ = distance / v_p = distance * porosity / (K * gradient)

    tau_days = (distance * porosity) / (K * gradient)

    # Uncertainty: assume 50% uncertainty on parameters
    uncertainty = tau_days * 0.5

    return float(tau_days), float(uncertainty), "gradient"


def _estimate_residence_time_tracer_decay(
    node_u: TemporalNode, node_v: TemporalNode, tracer_ion: str, ion_order: Optional[List[str]]
) -> Tuple[float, float, str]:
    """
    Estimate residence time using radioactive tracer decay.

    τ = -(t_1/2 / ln(2)) * ln(C_v / C_u)

    For tritium: t_1/2 = 12.32 years = 4498 days
    """
    if not node_u.samples or not node_v.samples:
        return 0.0, 0.0, "tracer_decay_no_data"

    # Find tracer index (assume tritium is in isotopes dict)
    # This is a placeholder - would need actual tritium data
    # For now, return error

    return 0.0, 0.0, "tracer_decay_not_implemented"


def compute_residence_time_from_velocity(
    distance_m: float, velocity_m_day: float, porosity: float = 0.2
) -> float:
    """
    Simple residence time calculation from known velocity.

    Parameters
    ----------
    distance_m : float
        Distance along flow path
    velocity_m_day : float
        Darcy velocity in m/day
    porosity : float
        Effective porosity

    Returns
    -------
    float
        Residence time in days
    """
    if velocity_m_day <= 0:
        return 0.0

    pore_velocity = velocity_m_day / porosity
    tau_days = distance_m / pore_velocity

    return tau_days
