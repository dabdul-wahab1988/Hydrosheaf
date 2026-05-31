"""
Residual bootstrap methods for uncertainty quantification.

This module implements residual bootstrapping to estimate confidence intervals
for transport parameters and reaction extents.
"""

from typing import List, Mapping, Optional, Tuple

import numpy as np

from ..config import Config
from . import UncertaintyResult


def bootstrap_edge_fit(
    x_u: List[float],
    x_v: List[float],
    config: Config,
    obs_u: Optional[Mapping[str, float]] = None,
    obs_v: Optional[Mapping[str, float]] = None,
    bounds: Optional[dict] = None,
    n_resamples: int = 1000,
    random_state: Optional[int] = None,
) -> UncertaintyResult:
    """
    Perform residual bootstrap on edge fit.

    Parameters
    ----------
    x_u, x_v : List[float]
        Upstream and downstream concentration vectors
    config : Config
        Hydrosheaf configuration
    n_resamples : int
        Number of bootstrap samples
    random_state : Optional[int]
        Random seed for reproducibility

    Returns
    -------
    UncertaintyResult
        Bootstrap confidence intervals and standard errors

    Mathematical Implementation
    ---------------------------
    1. Fit original model:
       (γ̂, ξ̂) = fit_edge(x_u, x_v, config)
       r̂ = x_v - T(γ̂)·x_u - R·ξ̂

    2. For b = 1, ..., n_resamples:
       a. Sample residual indices: I_b = {i_1, ..., i_n} with replacement from {1,...,n}
       b. Construct pseudo-observation: x*_v = T(γ̂)·x_u + R·ξ̂ + [r̂_{i_1}, ..., r̂_{i_n}]
       c. Refit: (γ*_b, ξ*_b) = fit_edge(x_u, x*_v, config)

    3. Compute statistics:
       γ_mean = mean(γ*_b)
       γ_std = std(γ*_b)
       γ_ci = [quantile(γ*_b, 0.025), quantile(γ*_b, 0.975)]
       (same for each ξ_j)
    """
    # Import here to avoid circular imports
    from ..inference.edge_fit import fit_edge

    from copy import copy

    # Set random seed
    if random_state is not None:
        np.random.seed(random_state)

    # Convert to numpy arrays
    # x_u_vec = np.array(x_u, dtype=float)
    x_v_vec = np.array(x_v, dtype=float)

    # Fit original model
    base_config = copy(config)
    base_config.uncertainty_method = "none"
    original_result = fit_edge(
        x_u, x_v, base_config, obs_u=obs_u, obs_v=obs_v, bounds=bounds
    )
    
    # Fix the transport model for all bootstrap replicates to ensure consistency
    # in the parameter distributions (don't mix gamma and f)
    fixed_model = original_result.transport_model
    fixed_endmember = original_result.endmember_id
    
    # Create a temporary config with only the selected model enabled
    boot_config = copy(config)
    boot_config.uncertainty_method = "none"
    boot_config.transport_models_enabled = [fixed_model]
    if fixed_model == "mix" and fixed_endmember:
        # Filter to only the selected endmember
        if fixed_endmember in config.mixing_endmembers:
            boot_config.mixing_endmembers = {fixed_endmember: config.mixing_endmembers[fixed_endmember]}

    # Compute residuals

    residuals = np.array(original_result.residual_vector, dtype=float)
    n_ions = len(residuals)

    # Storage for bootstrap samples
    gamma_samples = []
    f_samples = []
    m_extents = (
        len(original_result.reaction_fit.extents) if original_result.reaction_fit else 0
    )
    extents_samples = [[] for _ in range(m_extents)]
    gamma_jackknife: List[float] = []
    f_jackknife: List[float] = []
    extents_jackknife: List[List[float]] = [[] for _ in range(m_extents)]

    if getattr(config, "bootstrap_ci_method", "percentile") == "bca":
        for leave_out in range(n_ions):
            jack_config = copy(boot_config)
            jack_config.weights = list(getattr(boot_config, "weights", []))
            jack_config.conservative_weights = list(
                getattr(boot_config, "conservative_weights", [])
            )
            if leave_out < len(jack_config.weights):
                jack_config.weights[leave_out] = 0.0
            if leave_out < len(jack_config.conservative_weights):
                jack_config.conservative_weights[leave_out] = 0.0
            try:
                jack_result = fit_edge(
                    x_u,
                    x_v,
                    jack_config,
                    obs_u=obs_u,
                    obs_v=obs_v,
                    bounds=bounds,
                )
                if jack_result.transport_model == "evap" and jack_result.gamma is not None:
                    gamma_jackknife.append(float(jack_result.gamma))
                elif jack_result.transport_model == "mix" and jack_result.f is not None:
                    f_jackknife.append(float(jack_result.f))
                if jack_result.reaction_fit:
                    for j, extent in enumerate(jack_result.reaction_fit.extents):
                        if j < len(extents_jackknife):
                            extents_jackknife[j].append(float(extent))
            except Exception:
                continue

    # Bootstrap iterations
    for _ in range(n_resamples):
        # Sample residuals with replacement
        resample_indices = np.random.choice(n_ions, size=n_ions, replace=True)
        resampled_residuals = residuals[resample_indices]

        # Reconstruct fitted values from original fit
        fitted_values = x_v_vec - residuals

        # Create pseudo-observation
        x_v_star = fitted_values + resampled_residuals

        # Refit model
        try:
            boot_result = fit_edge(
                x_u,
                x_v_star.tolist(),
                boot_config,
                obs_u=obs_u,
                obs_v=obs_v,
                bounds=bounds,
            )


            # Store parameters
            if boot_result.transport_model == "evap":
                gamma_samples.append(boot_result.gamma)
            elif boot_result.transport_model == "mix":
                f_samples.append(boot_result.f)

            if boot_result.reaction_fit:
                for j, extent in enumerate(boot_result.reaction_fit.extents):
                    extents_samples[j].append(extent)

        except Exception:
            # Skip failed bootstrap samples
            continue

    # Compute statistics
    result = UncertaintyResult(method="bootstrap", n_resamples=n_resamples)

    if gamma_samples:
        gamma_arr = np.array(gamma_samples)
        result.gamma_mean = float(np.mean(gamma_arr))
        result.gamma_std = float(np.std(gamma_arr, ddof=1))
        if getattr(config, "bootstrap_ci_method", "percentile") == "bca":
            result.gamma_ci_low, result.gamma_ci_high = compute_bca_ci(
                gamma_arr,
                float(original_result.gamma if original_result.gamma is not None else result.gamma_mean),
                jackknife_samples=np.array(gamma_jackknife) if gamma_jackknife else None,
            )
        else:
            result.gamma_ci_low = float(np.percentile(gamma_arr, 2.5))
            result.gamma_ci_high = float(np.percentile(gamma_arr, 97.5))
        result.gamma_samples = gamma_arr

    if f_samples:
        f_arr = np.array(f_samples)
        result.f_mean = float(np.mean(f_arr))
        result.f_std = float(np.std(f_arr, ddof=1))
        if getattr(config, "bootstrap_ci_method", "percentile") == "bca":
            result.f_ci_low, result.f_ci_high = compute_bca_ci(
                f_arr,
                float(original_result.f if original_result.f is not None else result.f_mean),
                jackknife_samples=np.array(f_jackknife) if f_jackknife else None,
            )
        else:
            result.f_ci_low = float(np.percentile(f_arr, 2.5))
            result.f_ci_high = float(np.percentile(f_arr, 97.5))

    if m_extents > 0:
        for j in range(m_extents):
            if extents_samples[j]:
                extent_arr = np.array(extents_samples[j])
                result.extents_mean.append(float(np.mean(extent_arr)))
                result.extents_std.append(float(np.std(extent_arr, ddof=1)))
                if getattr(config, "bootstrap_ci_method", "percentile") == "bca":
                    point = float(original_result.reaction_fit.extents[j]) if original_result.reaction_fit else float(np.mean(extent_arr))
                    ci_low, ci_high = compute_bca_ci(
                        extent_arr,
                        point,
                        jackknife_samples=np.array(extents_jackknife[j]) if extents_jackknife[j] else None,
                    )
                    result.extents_ci_low.append(ci_low)
                    result.extents_ci_high.append(ci_high)
                else:
                    result.extents_ci_low.append(float(np.percentile(extent_arr, 2.5)))
                    result.extents_ci_high.append(float(np.percentile(extent_arr, 97.5)))
            else:
                result.extents_mean.append(0.0)
                result.extents_std.append(0.0)
                result.extents_ci_low.append(0.0)
                result.extents_ci_high.append(0.0)

        # Store all extents samples as 2D array
        if extents_samples[0]:
            result.extents_samples = np.array(
                [
                    [
                        s[i] if i < len(s) else 0.0
                        for i in range(len(extents_samples[0]))
                    ]
                    for s in extents_samples
                ]
            ).T

    return result


def bootstrap_reaction_fit(
    residual: List[float],
    reaction_matrix: List[List[float]],
    weights: List[float],
    lambda_l1: float,
    config: Config,
    n_resamples: int = 1000,
    random_state: Optional[int] = None,
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """
    Bootstrap confidence intervals for reaction extents only.

    Parameters
    ----------
    residual : List[float]
        Post-transport residual vector
    reaction_matrix : List[List[float]]
        Stoichiometry matrix
    weights : List[float]
        Ion weights
    lambda_l1 : float
        L1 regularization parameter
    config : Config
        Configuration
    n_resamples : int
        Number of bootstrap samples
    random_state : Optional[int]
        Random seed

    Returns
    -------
    Tuple[List[float], List[float], List[float], List[float]]
        (extents_mean, extents_std, extents_ci_low, extents_ci_high)
    """
    # Import here to avoid circular imports
    from ..models.reactions import fit_reactions

    if random_state is not None:
        np.random.seed(random_state)

    residual_vec = np.array(residual, dtype=float)
    n_ions = len(residual_vec)

    # Fit original model
    original_fit = fit_reactions(residual, reaction_matrix, weights, lambda_l1, config)
    m_reactions = len(original_fit.extents)

    # Bootstrap samples storage
    extents_samples = [[] for _ in range(m_reactions)]

    for _ in range(n_resamples):
        # Sample residuals with replacement
        resample_indices = np.random.choice(n_ions, size=n_ions, replace=True)
        resampled_residual = residual_vec[resample_indices]

        # Refit
        try:
            boot_fit = fit_reactions(
                resampled_residual.tolist(), reaction_matrix, weights, lambda_l1, config
            )
            for j, extent in enumerate(boot_fit.extents):
                extents_samples[j].append(extent)
        except Exception:
            continue

    # Compute statistics
    extents_mean = []
    extents_std = []
    extents_ci_low = []
    extents_ci_high = []

    for j in range(m_reactions):
        if extents_samples[j]:
            extent_arr = np.array(extents_samples[j])
            extents_mean.append(float(np.mean(extent_arr)))
            extents_std.append(float(np.std(extent_arr, ddof=1)))
            extents_ci_low.append(float(np.percentile(extent_arr, 2.5)))
            extents_ci_high.append(float(np.percentile(extent_arr, 97.5)))
        else:
            extents_mean.append(0.0)
            extents_std.append(0.0)
            extents_ci_low.append(0.0)
            extents_ci_high.append(0.0)

    return extents_mean, extents_std, extents_ci_low, extents_ci_high


def compute_bca_ci(
    samples: np.ndarray,
    point_estimate: float,
    alpha: float = 0.05,
    jackknife_samples: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    """
    Compute bias-corrected accelerated (BCa) bootstrap confidence interval.

    Parameters
    ----------
    samples : np.ndarray
        Bootstrap samples
    point_estimate : float
        Original point estimate
    alpha : float
        Significance level (default 0.05 for 95% CI)

    Returns
    -------
    Tuple[float, float]
        (lower_bound, upper_bound)

    Mathematical Implementation
    ---------------------------
    1. Compute bias correction:
       z0 = Φ^{-1}(#{θ* < θ̂} / B)

    2. Compute acceleration from delete-one jackknife estimates:
       a = Σ(θ̄_jack - θ_jack_i)^3 / [6 * {Σ(θ̄_jack - θ_jack_i)^2}^{3/2}]

    3. Adjust percentiles:
       α_low = Φ(z0 + (z0 + z_{α/2}) / (1 - a(z0 + z_{α/2})))
       α_high = Φ(z0 + (z0 + z_{1-α/2}) / (1 - a(z0 + z_{1-α/2})))
    """
    from scipy.stats import norm

    clean_samples = np.asarray(samples, dtype=float)
    clean_samples = clean_samples[np.isfinite(clean_samples)]
    n_boot = len(clean_samples)
    if n_boot == 0:
        return float("nan"), float("nan")

    # Bias correction
    n_less = np.sum(clean_samples < point_estimate)
    p_less = (n_less + 0.5) / (n_boot + 1.0)
    p_less = float(np.clip(p_less, 1e-12, 1.0 - 1e-12))
    z0 = norm.ppf(p_less)

    a = 0.0
    if jackknife_samples is not None:
        jack = np.asarray(jackknife_samples, dtype=float)
        jack = jack[np.isfinite(jack)]
        if len(jack) >= 2:
            jack_mean = float(np.mean(jack))
            diffs = jack_mean - jack
            numerator = float(np.sum(diffs**3))
            denominator = float(6.0 * (np.sum(diffs**2) ** 1.5))
            if denominator > 0.0:
                a = numerator / denominator

    # Adjusted percentiles
    z_alpha = norm.ppf(alpha / 2)
    z_1alpha = norm.ppf(1 - alpha / 2)

    denom_low = 1 - a * (z0 + z_alpha)
    denom_high = 1 - a * (z0 + z_1alpha)
    if abs(denom_low) < 1e-12:
        alpha_low = 0.0 if denom_low < 0 else 1.0
    else:
        alpha_low = norm.cdf(z0 + (z0 + z_alpha) / denom_low)
    if abs(denom_high) < 1e-12:
        alpha_high = 0.0 if denom_high < 0 else 1.0
    else:
        alpha_high = norm.cdf(z0 + (z0 + z_1alpha) / denom_high)

    # Clip to valid range
    alpha_low = np.clip(alpha_low, 0, 1)
    alpha_high = np.clip(alpha_high, 0, 1)

    ci_low = float(np.percentile(clean_samples, 100 * alpha_low))
    ci_high = float(np.percentile(clean_samples, 100 * alpha_high))

    return ci_low, ci_high
