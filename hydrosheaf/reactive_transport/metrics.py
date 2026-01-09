"""
Forward-inverse consistency metrics for reactive transport validation.
"""

from typing import Dict, List, Optional

import numpy as np


def compute_consistency_metrics(
    x_forward: List[float],
    x_observed: List[float],
    weights: Optional[List[float]] = None,
) -> Dict[str, float]:
    """
    Compute inverse-forward consistency metrics.

    Parameters
    ----------
    x_forward : List[float]
        Forward model prediction
    x_observed : List[float]
        Observed downstream composition
    weights : Optional[List[float]]
        Ion weights for weighted metrics

    Returns
    -------
    Dict[str, float]
        {"rmse": float, "nse": float, "pbias": float, "r2": float}

    Mathematical Implementation
    ---------------------------
    n = len(x_observed)

    RMSE:
        rmse = sqrt(sum((x_f[j] - x_o[j])² for j in range(n)) / n)

    NSE (Nash-Sutcliffe Efficiency):
        ss_res = sum((x_f[j] - x_o[j])² for j in range(n))
        ss_tot = sum((x_o[j] - mean(x_o))² for j in range(n))
        nse = 1 - ss_res / ss_tot if ss_tot > 0 else -inf

    PBIAS (Percent Bias):
        pbias = 100 * sum(x_f[j] - x_o[j] for j in range(n)) / sum(x_o)

    R²:
        r2 = correlation(x_forward, x_observed)²
    """
    x_f = np.array(x_forward, dtype=float)
    x_o = np.array(x_observed, dtype=float)

    if len(x_f) != len(x_o):
        raise ValueError("Forward and observed vectors must have same length")

    n = len(x_f)

    if weights is not None:
        w = np.array(weights, dtype=float)
    else:
        w = np.ones(n)

    # RMSE (Root Mean Square Error)
    squared_errors = (x_f - x_o) ** 2
    rmse = float(np.sqrt(np.sum(w * squared_errors) / np.sum(w)))

    # NSE (Nash-Sutcliffe Efficiency)
    mean_observed = np.mean(x_o)
    ss_res = np.sum(w * (x_f - x_o) ** 2)
    ss_tot = np.sum(w * (x_o - mean_observed) ** 2)

    if ss_tot > 1e-10:
        nse = float(1.0 - ss_res / ss_tot)
    else:
        nse = float("-inf")

    # PBIAS (Percent Bias)
    sum_obs = np.sum(x_o)
    if sum_obs > 1e-10:
        pbias = float(100.0 * np.sum(x_f - x_o) / sum_obs)
    else:
        pbias = 0.0

    # R² (Coefficient of Determination)
    if np.std(x_f) > 1e-10 and np.std(x_o) > 1e-10:
        corr = np.corrcoef(x_f, x_o)[0, 1]
        r2 = float(corr**2)
    else:
        r2 = 0.0

    return {
        "rmse": rmse,
        "nse": nse,
        "pbias": pbias,
        "r2": r2,
    }


def compute_per_ion_metrics(
    x_forward: List[float],
    x_observed: List[float],
    ion_order: List[str],
) -> Dict[str, Dict[str, float]]:
    """
    Compute consistency metrics for each ion separately.

    Parameters
    ----------
    x_forward : List[float]
        Forward prediction
    x_observed : List[float]
        Observed concentrations
    ion_order : List[str]
        Ion names

    Returns
    -------
    Dict[str, Dict[str, float]]
        Mapping from ion name to metrics dict
    """
    metrics = {}

    for i, ion in enumerate(ion_order):
        if i < len(x_forward) and i < len(x_observed):
            error = x_forward[i] - x_observed[i]
            rel_error = error / x_observed[i] if x_observed[i] > 1e-10 else 0.0

            metrics[ion] = {
                "forward": x_forward[i],
                "observed": x_observed[i],
                "error": error,
                "rel_error": rel_error * 100.0,  # as percentage
            }

    return metrics


def check_thermodynamic_consistency(
    extents: List[float],
    si_initial: Dict[str, float],
    reaction_labels: List[str],
    si_threshold: float = 0.2,
) -> tuple[bool, List[str]]:
    """
    Check thermodynamic consistency of inverse results.

    Dissolution should only occur when undersaturated (SI < 0).
    Precipitation should only occur when supersaturated (SI > 0).

    Parameters
    ----------
    extents : List[float]
        Reaction extents from inverse model (mmol/L)
    si_initial : Dict[str, float]
        Initial saturation indices
    reaction_labels : List[str]
        Reaction names
    si_threshold : float
        Tolerance for SI near zero

    Returns
    -------
    tuple[bool, List[str]]
        (is_consistent, list_of_violations)

    Logic
    -----
    For each reaction k:
        If ξ_k > 0 (dissolution):
            Require SI_k(t=0) < τ (undersaturated initially)
        If ξ_k < 0 (precipitation):
            Require SI_k(t=0) > -τ (supersaturated initially)
    """
    violations = []
    is_consistent = True

    for extent, label in zip(extents, reaction_labels):
        # Skip if extent is negligible
        if abs(extent) < 1e-6:
            continue

        # Get initial SI (if available)
        si = si_initial.get(label, None)
        if si is None:
            continue

        # Check consistency
        if extent > 0:
            # Dissolution: should be undersaturated
            if si > si_threshold:
                violations.append(
                    f"{label}: dissolution (ξ={extent:.3f}) but supersaturated (SI={si:.2f})"
                )
                is_consistent = False

        elif extent < 0:
            # Precipitation: should be supersaturated
            if si < -si_threshold:
                violations.append(
                    f"{label}: precipitation (ξ={extent:.3f}) but undersaturated (SI={si:.2f})"
                )
                is_consistent = False

    return is_consistent, violations


def compute_mass_balance_error(
    x_forward: List[float],
    x_observed: List[float],
    ion_order: List[str],
) -> Dict[str, float]:
    """
    Compute mass balance error for key elements.

    Parameters
    ----------
    x_forward : List[float]
        Forward prediction
    x_observed : List[float]
        Observed concentrations
    ion_order : List[str]
        Ion names

    Returns
    -------
    Dict[str, float]
        Mass balance errors for key elements (%)
    """
    errors = {}

    # Create ion index map
    ion_idx = {ion: i for i, ion in enumerate(ion_order)}

    # Calcium balance
    if "Ca" in ion_idx:
        i = ion_idx["Ca"]
        if i < len(x_forward) and i < len(x_observed):
            if x_observed[i] > 1e-10:
                ca_error = 100.0 * (x_forward[i] - x_observed[i]) / x_observed[i]
                errors["Ca"] = ca_error

    # Sulfur balance (SO4)
    if "SO4" in ion_idx:
        i = ion_idx["SO4"]
        if i < len(x_forward) and i < len(x_observed):
            if x_observed[i] > 1e-10:
                s_error = 100.0 * (x_forward[i] - x_observed[i]) / x_observed[i]
                errors["S"] = s_error

    # Carbon balance (HCO3)
    if "HCO3" in ion_idx:
        i = ion_idx["HCO3"]
        if i < len(x_forward) and i < len(x_observed):
            if x_observed[i] > 1e-10:
                c_error = 100.0 * (x_forward[i] - x_observed[i]) / x_observed[i]
                errors["C"] = c_error

    # Nitrogen balance (NO3)
    if "NO3" in ion_idx:
        i = ion_idx["NO3"]
        if i < len(x_forward) and i < len(x_observed):
            if x_observed[i] > 1e-10:
                n_error = 100.0 * (x_forward[i] - x_observed[i]) / x_observed[i]
                errors["N"] = n_error

    return errors


def compute_reaction_rate_from_extent(
    extent_mmol_L: float, residence_time_days: float
) -> float:
    """
    Convert reaction extent to average rate.

    Parameters
    ----------
    extent_mmol_L : float
        Cumulative extent in mmol/L
    residence_time_days : float
        Residence time in days

    Returns
    -------
    float
        Average rate in mmol/L/day

    Formula
    -------
    r̄ = ξ / τ
    """
    if residence_time_days <= 0:
        return 0.0

    return extent_mmol_L / residence_time_days
