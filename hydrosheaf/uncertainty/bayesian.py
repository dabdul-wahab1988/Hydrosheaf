"""
Bayesian MCMC methods for uncertainty quantification.

This module implements Bayesian inference using MCMC (Hamiltonian Monte Carlo / NUTS)
to estimate posterior distributions of transport parameters and reaction extents.

Requirements:
    pip install pymc>=5.0 arviz>=0.15

Alternative:
    pip install numpyro>=0.12 jax>=0.4
"""

from typing import Dict, List, Optional, Tuple

import numpy as np

from . import UncertaintyResult


def bayesian_edge_fit(
    x_u: List[float],
    x_v: List[float],
    reaction_matrix: List[List[float]],
    reaction_labels: List[str],
    config: "Config",  # type: ignore
    n_samples: int = 5000,
    n_chains: int = 4,
    target_accept: float = 0.95,
    bounds: Optional[List[Tuple[Optional[float], Optional[float]]]] = None,
) -> UncertaintyResult:
    """
    Bayesian inference for edge parameters using MCMC.

    Parameters
    ----------
    x_u, x_v : List[float]
        Concentration vectors
    reaction_matrix : List[List[float]]
        Stoichiometry matrix R (m reactions x n ions)
    reaction_labels : List[str]
        Names of reactions
    config : Config
        Hydrosheaf configuration
    n_samples : int
        Number of MCMC draws per chain
    n_chains : int
        Number of independent chains
    target_accept : float
        Target acceptance rate for NUTS
    bounds : Optional[List[Tuple]]
        (lower, upper) bounds for each reaction extent

    Returns
    -------
    UncertaintyResult
        Posterior summaries and samples

    Mathematical Implementation
    ---------------------------
    1. Define priors:
       γ ~ TruncatedNormal(1.0, 0.5, lower=1.0)
       ξ_j ~ Laplace(0, b=1/λ) for j = 1,...,m

    2. Define likelihood:
       x_v | γ, ξ ~ Normal(γ·x_u + R·ξ, σ²I)
       σ ~ HalfNormal(0.1)

    3. Apply thermodynamic bounds as hard constraints via Potential

    4. Run NUTS sampler with n_chains chains, n_samples draws each

    5. Discard first n_samples/2 as warmup

    6. Compute posterior summaries:
       - Mean, std for each parameter
       - 2.5%, 97.5% quantiles for CI
       - R̂ and ESS for convergence
    """
    try:
        import pymc as pm
        import arviz as az
    except ImportError:
        raise ImportError(
            "PyMC is required for Bayesian inference. "
            "Install with: pip install pymc>=5.0 arviz>=0.15"
        )

    # Convert to numpy
    x_u_vec = np.array(x_u, dtype=float)
    x_v_vec = np.array(x_v, dtype=float)
    R = np.array(reaction_matrix, dtype=float).T  # n_ions x n_reactions
    n_ions, m_reactions = R.shape

    # Get prior hyperparameters from config
    prior_gamma_mu = getattr(config, "prior_gamma_mu", 1.0)
    prior_gamma_sigma = getattr(config, "prior_gamma_sigma", 0.5)
    prior_xi_scale = getattr(config, "prior_xi_scale", 1.0)
    prior_sigma_scale = getattr(config, "prior_sigma_scale", 0.1)

    # Build model
    with pm.Model() as model:
        # Priors for transport parameter (evaporation model)
        gamma = pm.TruncatedNormal(
            "gamma", mu=prior_gamma_mu, sigma=prior_gamma_sigma, lower=1.0
        )

        # Priors for reaction extents (Laplace for sparsity)
        xi = pm.Laplace("xi", mu=0, b=prior_xi_scale, shape=m_reactions)

        # Apply thermodynamic bounds if provided
        if bounds is not None:
            for j, (lb, ub) in enumerate(bounds):
                if lb is not None:
                    pm.Potential(f"lb_{j}", pm.math.switch(xi[j] < lb, -1e10, 0))
                if ub is not None:
                    pm.Potential(f"ub_{j}", pm.math.switch(xi[j] > ub, -1e10, 0))

        # Predicted downstream concentration
        x_pred = gamma * x_u_vec + pm.math.dot(R, xi)

        # Likelihood with heteroscedastic noise
        sigma = pm.HalfNormal("sigma", sigma=prior_sigma_scale, shape=n_ions)
        pm.Normal("x_v_obs", mu=x_pred, sigma=sigma, observed=x_v_vec)

        # Sample
        trace = pm.sample(
            draws=n_samples,
            chains=n_chains,
            cores=min(n_chains, 4),
            target_accept=target_accept,
            return_inferencedata=True,
            progressbar=False,
        )

    # Extract posterior samples
    gamma_samples = trace.posterior["gamma"].values.flatten()
    xi_samples = trace.posterior["xi"].values  # shape: (chains, draws, m_reactions)
    xi_samples = xi_samples.reshape(-1, m_reactions)

    # Compute summaries
    result = UncertaintyResult(method="bayesian")

    # Gamma statistics
    result.gamma_mean = float(np.mean(gamma_samples))
    result.gamma_std = float(np.std(gamma_samples, ddof=1))
    result.gamma_ci_low = float(np.percentile(gamma_samples, 2.5))
    result.gamma_ci_high = float(np.percentile(gamma_samples, 97.5))
    result.gamma_samples = gamma_samples

    # Reaction extent statistics
    for j in range(m_reactions):
        xi_j = xi_samples[:, j]
        result.extents_mean.append(float(np.mean(xi_j)))
        result.extents_std.append(float(np.std(xi_j, ddof=1)))
        result.extents_ci_low.append(float(np.percentile(xi_j, 2.5)))
        result.extents_ci_high.append(float(np.percentile(xi_j, 97.5)))

    result.extents_samples = xi_samples

    # Convergence diagnostics
    result.r_hat = {}
    result.ess = {}

    try:
        # R-hat (Gelman-Rubin statistic)
        r_hat_gamma = float(az.rhat(trace.posterior["gamma"]).values)
        result.r_hat["gamma"] = r_hat_gamma

        for j, label in enumerate(reaction_labels):
            r_hat_xi = float(az.rhat(trace.posterior["xi"].sel(xi_dim_0=j)).values)
            result.r_hat[f"xi_{label}"] = r_hat_xi

        # ESS (Effective Sample Size)
        ess_gamma = float(az.ess(trace.posterior["gamma"]).values)
        result.ess["gamma"] = ess_gamma

        for j, label in enumerate(reaction_labels):
            ess_xi = float(az.ess(trace.posterior["xi"].sel(xi_dim_0=j)).values)
            result.ess[f"xi_{label}"] = ess_xi

    except Exception as e:
        # Convergence diagnostics may fail in some edge cases
        print(f"Warning: Could not compute convergence diagnostics: {e}")

    return result


def bayesian_reaction_fit(
    residual: List[float],
    reaction_matrix: List[List[float]],
    reaction_labels: List[str],
    weights: List[float],
    lambda_l1: float,
    config: "Config",  # type: ignore
    n_samples: int = 5000,
    n_chains: int = 4,
    bounds: Optional[List[Tuple[Optional[float], Optional[float]]]] = None,
) -> Tuple[List[float], List[float], List[float], List[float], Dict[str, float], Dict[str, float]]:
    """
    Bayesian inference for reaction extents only (no transport parameter).

    Parameters
    ----------
    residual : List[float]
        Post-transport residual vector
    reaction_matrix : List[List[float]]
        Stoichiometry matrix
    reaction_labels : List[str]
        Reaction names
    weights : List[float]
        Ion weights
    lambda_l1 : float
        L1 regularization parameter (used to set Laplace prior scale)
    config : Config
        Configuration
    n_samples : int
        MCMC samples
    n_chains : int
        Number of chains
    bounds : Optional[List[Tuple]]
        Reaction extent bounds

    Returns
    -------
    Tuple
        (extents_mean, extents_std, extents_ci_low, extents_ci_high, r_hat, ess)
    """
    try:
        import pymc as pm
        import arviz as az
    except ImportError:
        raise ImportError(
            "PyMC is required for Bayesian inference. "
            "Install with: pip install pymc>=5.0 arviz>=0.15"
        )

    residual_vec = np.array(residual, dtype=float)
    R = np.array(reaction_matrix, dtype=float).T  # n_ions x n_reactions
    W = np.array(weights, dtype=float)
    n_ions, m_reactions = R.shape

    # Laplace prior scale from L1 lambda
    laplace_scale = 1.0 / lambda_l1 if lambda_l1 > 0 else 1.0

    with pm.Model() as model:
        # Priors
        xi = pm.Laplace("xi", mu=0, b=laplace_scale, shape=m_reactions)

        # Apply bounds
        if bounds is not None:
            for j, (lb, ub) in enumerate(bounds):
                if lb is not None:
                    pm.Potential(f"lb_{j}", pm.math.switch(xi[j] < lb, -1e10, 0))
                if ub is not None:
                    pm.Potential(f"ub_{j}", pm.math.switch(xi[j] > ub, -1e10, 0))

        # Predicted residual (should be near zero after reaction fit)
        residual_pred = pm.math.dot(R, xi)

        # Weighted likelihood
        sigma = pm.math.sqrt(1.0 / W)  # weights are 1/sigma^2
        pm.Normal("residual_obs", mu=residual_pred, sigma=sigma, observed=residual_vec)

        # Sample
        trace = pm.sample(
            draws=n_samples,
            chains=n_chains,
            cores=min(n_chains, 4),
            target_accept=0.95,
            return_inferencedata=True,
            progressbar=False,
        )

    # Extract samples
    xi_samples = trace.posterior["xi"].values.reshape(-1, m_reactions)

    # Compute statistics
    extents_mean = []
    extents_std = []
    extents_ci_low = []
    extents_ci_high = []

    for j in range(m_reactions):
        xi_j = xi_samples[:, j]
        extents_mean.append(float(np.mean(xi_j)))
        extents_std.append(float(np.std(xi_j, ddof=1)))
        extents_ci_low.append(float(np.percentile(xi_j, 2.5)))
        extents_ci_high.append(float(np.percentile(xi_j, 97.5)))

    # Convergence diagnostics
    r_hat = {}
    ess = {}

    try:
        for j, label in enumerate(reaction_labels):
            r_hat[label] = float(az.rhat(trace.posterior["xi"].sel(xi_dim_0=j)).values)
            ess[label] = float(az.ess(trace.posterior["xi"].sel(xi_dim_0=j)).values)
    except Exception:
        pass

    return extents_mean, extents_std, extents_ci_low, extents_ci_high, r_hat, ess


def compute_r_hat(chains: np.ndarray) -> float:
    """
    Compute Gelman-Rubin R-hat statistic for convergence.

    Parameters
    ----------
    chains : np.ndarray
        Shape (n_chains, n_samples)

    Returns
    -------
    float
        R-hat value. Should be < 1.01 for convergence.

    Mathematical Implementation
    ---------------------------
    1. Compute within-chain variance:
       W = (1/m) Σ_j s_j²
       where s_j² is variance of chain j

    2. Compute between-chain variance:
       B = (n/(m-1)) Σ_j (θ̄_j - θ̄_total)²

    3. Estimate marginal variance:
       Var(θ|y) = ((n-1)/n)W + (1/n)B

    4. R̂ = sqrt(Var(θ|y) / W)
    """
    n_chains, n_samples = chains.shape

    # Chain means
    chain_means = np.mean(chains, axis=1)
    overall_mean = np.mean(chain_means)

    # Within-chain variance
    chain_variances = np.var(chains, axis=1, ddof=1)
    W = np.mean(chain_variances)

    # Between-chain variance
    B = n_samples * np.var(chain_means, ddof=1)

    # Marginal posterior variance
    var_plus = ((n_samples - 1) / n_samples) * W + (1 / n_samples) * B

    # R-hat
    r_hat = np.sqrt(var_plus / W) if W > 0 else 1.0

    return float(r_hat)


def compute_ess(samples: np.ndarray, max_lag: int = 100) -> float:
    """
    Compute effective sample size accounting for autocorrelation.

    Parameters
    ----------
    samples : np.ndarray
        MCMC samples (1D array)
    max_lag : int
        Maximum lag for autocorrelation computation

    Returns
    -------
    float
        Effective sample size

    Mathematical Implementation
    ---------------------------
    ESS = N / (1 + 2Σ_{k=1}^∞ ρ_k)

    Where ρ_k is the lag-k autocorrelation.
    """
    n = len(samples)

    # Demean
    samples_demeaned = samples - np.mean(samples)

    # Autocorrelation
    autocorr = np.correlate(samples_demeaned, samples_demeaned, mode="full")
    autocorr = autocorr[n - 1 :]  # positive lags only
    autocorr = autocorr / autocorr[0]  # normalize

    # Sum autocorrelations until they become negligible
    rho_sum = 0.0
    for k in range(1, min(max_lag, len(autocorr))):
        if autocorr[k] < 0.05:  # stop when autocorrelation is small
            break
        rho_sum += autocorr[k]

    ess = n / (1 + 2 * rho_sum)

    return float(ess)
