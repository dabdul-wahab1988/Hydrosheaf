"""
Monte Carlo error propagation for uncertainty quantification.

This module implements Monte Carlo sampling to propagate input measurement
uncertainty through the transport and reaction inference pipeline.
"""

from typing import List, Optional

import numpy as np

from . import UncertaintyResult


def monte_carlo_propagate(
    x_u: List[float],
    x_v: List[float],
    config: "Config",  # type: ignore
    input_uncertainty_pct: float = 5.0,
    n_samples: int = 1000,
    random_state: Optional[int] = None,
) -> UncertaintyResult:
    """
    Propagate input measurement uncertainty through the model.

    Parameters
    ----------
    x_u, x_v : List[float]
        Observed concentration vectors
    config : Config
        Hydrosheaf configuration
    input_uncertainty_pct : float
        Relative uncertainty as percentage (e.g., 5 means 5%)
    n_samples : int
        Number of Monte Carlo samples
    random_state : Optional[int]
        Random seed

    Returns
    -------
    UncertaintyResult
        Parameter distributions reflecting input uncertainty

    Mathematical Implementation
    ---------------------------
    1. Compute input standard deviations:
       σ_u[j] = (input_uncertainty_pct / 100) * x_u[j]
       σ_v[j] = (input_uncertainty_pct / 100) * x_v[j]

    2. For k = 1, ..., n_samples:
       a. Sample perturbed inputs:
          x_u^(k)[j] ~ Normal(x_u[j], σ_u[j])
          x_v^(k)[j] ~ Normal(x_v[j], σ_v[j])
       b. Ensure non-negativity: x^(k)[j] = max(0, x^(k)[j])
       c. Fit model: (γ^(k), ξ^(k)) = fit_edge(x_u^(k), x_v^(k), config)

    3. Compute statistics from {γ^(k), ξ^(k)}
    """
    # Import here to avoid circular imports
    from ..inference.edge_fit import fit_edge

    if random_state is not None:
        np.random.seed(random_state)

    # Convert to numpy
    x_u_vec = np.array(x_u, dtype=float)
    x_v_vec = np.array(x_v, dtype=float)
    n_ions = len(x_u_vec)

    # Compute standard deviations from relative uncertainty
    sigma_u = (input_uncertainty_pct / 100.0) * x_u_vec
    sigma_v = (input_uncertainty_pct / 100.0) * x_v_vec

    # Avoid zero uncertainty
    sigma_u = np.maximum(sigma_u, 1e-6)
    sigma_v = np.maximum(sigma_v, 1e-6)

    # Storage for samples
    gamma_samples = []
    f_samples = []
    extents_samples = []
    n_reactions = 0

    # Monte Carlo sampling
    for _ in range(n_samples):
        # Sample perturbed inputs
        x_u_pert = np.random.normal(x_u_vec, sigma_u)
        x_v_pert = np.random.normal(x_v_vec, sigma_v)

        # Enforce non-negativity
        x_u_pert = np.maximum(x_u_pert, 0.0)
        x_v_pert = np.maximum(x_v_pert, 0.0)

        # Fit model
        try:
            result = fit_edge(x_u_pert.tolist(), x_v_pert.tolist(), {}, "", config)

            # Store parameters
            if result.transport_model == "evap":
                gamma_samples.append(result.gamma)
            elif result.transport_model == "mix":
                f_samples.append(result.f)

            if result.reaction_fit:
                extents_samples.append(result.reaction_fit.extents)
                n_reactions = len(result.reaction_fit.extents)

        except Exception:
            # Skip failed samples
            continue

    # Compute statistics
    result = UncertaintyResult(method="monte_carlo", n_resamples=n_samples)

    if gamma_samples:
        gamma_arr = np.array(gamma_samples)
        result.gamma_mean = float(np.mean(gamma_arr))
        result.gamma_std = float(np.std(gamma_arr, ddof=1))
        result.gamma_ci_low = float(np.percentile(gamma_arr, 2.5))
        result.gamma_ci_high = float(np.percentile(gamma_arr, 97.5))
        result.gamma_samples = gamma_arr

    if f_samples:
        f_arr = np.array(f_samples)
        result.f_mean = float(np.mean(f_arr))
        result.f_std = float(np.std(f_arr, ddof=1))
        result.f_ci_low = float(np.percentile(f_arr, 2.5))
        result.f_ci_high = float(np.percentile(f_arr, 97.5))

    if extents_samples:
        extents_arr = np.array(extents_samples)  # shape: (n_samples, n_reactions)

        for j in range(n_reactions):
            extent_j = extents_arr[:, j]
            result.extents_mean.append(float(np.mean(extent_j)))
            result.extents_std.append(float(np.std(extent_j, ddof=1)))
            result.extents_ci_low.append(float(np.percentile(extent_j, 2.5)))
            result.extents_ci_high.append(float(np.percentile(extent_j, 97.5)))

        result.extents_samples = extents_arr

    return result


def propagate_variance_decomposition(
    x_u: List[float],
    x_v: List[float],
    config: "Config",  # type: ignore
    input_uncertainty_pct: float = 5.0,
    n_monte_carlo: int = 1000,
    n_bootstrap: int = 1000,
    random_state: Optional[int] = None,
) -> dict:
    """
    Decompose total variance into aleatory (input) and epistemic (model) components.

    Parameters
    ----------
    x_u, x_v : List[float]
        Concentration vectors
    config : Config
        Configuration
    input_uncertainty_pct : float
        Input uncertainty percentage
    n_monte_carlo : int
        Monte Carlo samples for aleatory uncertainty
    n_bootstrap : int
        Bootstrap samples for epistemic uncertainty
    random_state : Optional[int]
        Random seed

    Returns
    -------
    dict
        {
            "aleatory_variance": float,
            "epistemic_variance": float,
            "total_variance": float,
            "aleatory_fraction": float
        }

    Mathematical Implementation
    ---------------------------
    Total variance = Aleatory (input) + Epistemic (model):
    Var(γ) = Var_aleatory(γ) + Var_epistemic(γ)

    Estimate aleatory via Monte Carlo, epistemic via bootstrap.
    """
    from .bootstrap import bootstrap_edge_fit

    if random_state is not None:
        np.random.seed(random_state)

    # Monte Carlo for aleatory uncertainty
    mc_result = monte_carlo_propagate(
        x_u, x_v, config, input_uncertainty_pct, n_monte_carlo, random_state
    )

    # Bootstrap for epistemic uncertainty (with mean inputs)
    boot_result = bootstrap_edge_fit(x_u, x_v, config, n_bootstrap, random_state)

    decomposition = {}

    # Gamma variance decomposition
    if mc_result.gamma_std is not None and boot_result.gamma_std is not None:
        var_aleatory = mc_result.gamma_std**2
        var_epistemic = boot_result.gamma_std**2
        var_total = var_aleatory + var_epistemic

        decomposition["gamma"] = {
            "aleatory_variance": var_aleatory,
            "epistemic_variance": var_epistemic,
            "total_variance": var_total,
            "aleatory_fraction": var_aleatory / var_total if var_total > 0 else 0.0,
            "total_std": np.sqrt(var_total),
        }

    # Extents variance decomposition
    if mc_result.extents_std and boot_result.extents_std:
        n_reactions = len(mc_result.extents_std)
        decomposition["extents"] = []

        for j in range(n_reactions):
            var_aleatory = mc_result.extents_std[j] ** 2
            var_epistemic = boot_result.extents_std[j] ** 2
            var_total = var_aleatory + var_epistemic

            decomposition["extents"].append(
                {
                    "reaction_index": j,
                    "aleatory_variance": var_aleatory,
                    "epistemic_variance": var_epistemic,
                    "total_variance": var_total,
                    "aleatory_fraction": var_aleatory / var_total if var_total > 0 else 0.0,
                    "total_std": np.sqrt(var_total),
                }
            )

    return decomposition


def compute_sensitivity_indices(
    x_u: List[float],
    x_v: List[float],
    config: "Config",  # type: ignore
    perturbation_pct: float = 1.0,
) -> dict:
    """
    Compute local sensitivity indices for each input ion.

    Parameters
    ----------
    x_u, x_v : List[float]
        Concentration vectors
    config : Config
        Configuration
    perturbation_pct : float
        Perturbation percentage for finite differences

    Returns
    -------
    dict
        Sensitivity indices for gamma and reaction extents w.r.t. each input ion

    Mathematical Implementation
    ---------------------------
    For parameter θ and input x_i:
    S_i = ∂θ/∂x_i ≈ (θ(x_i + δ) - θ(x_i)) / δ

    Normalized: S_i^* = (x_i / θ) * (∂θ/∂x_i)
    """
    from ..inference.edge_fit import fit_edge

    # Baseline fit
    baseline = fit_edge(x_u, x_v, {}, "", config)

    n_ions = len(x_u)
    delta = perturbation_pct / 100.0

    sensitivity = {
        "gamma_sensitivity": [],
        "extents_sensitivity": [],
    }

    # Perturb each upstream ion
    for i in range(n_ions):
        x_u_pert = list(x_u)
        x_u_pert[i] *= 1 + delta

        try:
            pert_result = fit_edge(x_u_pert, x_v, {}, "", config)

            # Gamma sensitivity
            if baseline.transport_model == "evap":
                d_gamma = (pert_result.gamma - baseline.gamma) / (x_u[i] * delta)
                sensitivity["gamma_sensitivity"].append(
                    {"ion_index": i, "sensitivity": d_gamma}
                )

            # Extents sensitivity
            if baseline.reaction_fit and pert_result.reaction_fit:
                for j, (base_ext, pert_ext) in enumerate(
                    zip(baseline.reaction_fit.extents, pert_result.reaction_fit.extents)
                ):
                    d_extent = (pert_ext - base_ext) / (x_u[i] * delta)
                    sensitivity["extents_sensitivity"].append(
                        {
                            "ion_index": i,
                            "reaction_index": j,
                            "sensitivity": d_extent,
                        }
                    )

        except Exception:
            pass

    return sensitivity
