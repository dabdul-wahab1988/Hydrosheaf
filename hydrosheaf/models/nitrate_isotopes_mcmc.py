"""
MCMC Bayesian Isotope Mixing Model for Nitrate Source Apportionment.

This module provides full MCMC posterior sampling using PyMC for quantifying
uncertainty in nitrate source fractions based on dual isotope (d15N, d18O) data.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any

import numpy as np

from .nitrate_isotopes import IsotopeSample, SourceIsotopes

# Try to import PyMC and ArviZ
try:
    import pymc as pm
    import arviz as az

    PYMC_AVAILABLE = True
except ImportError:
    pm = None
    az = None
    PYMC_AVAILABLE = False


@dataclass
class MCMCMixingResult:
    """Result of MCMC Bayesian isotope mixing analysis."""

    # Point estimates (posterior means)
    source_fractions: Dict[str, float]

    # Posterior samples (n_samples x n_sources)
    posterior_samples: Optional[np.ndarray] = None
    source_names: List[str] = field(default_factory=list)

    # Credible intervals (95% by default)
    ci_lower: Dict[str, float] = field(default_factory=dict)
    ci_upper: Dict[str, float] = field(default_factory=dict)

    # Diagnostics
    r_hat: Dict[str, float] = field(default_factory=dict)
    ess_bulk: Dict[str, float] = field(default_factory=dict)
    ess_tail: Dict[str, float] = field(default_factory=dict)

    # Model fit
    log_likelihood: Optional[float] = None
    waic: Optional[float] = None

    # Convergence flag
    converged: bool = True
    warnings: List[str] = field(default_factory=list)


def check_pymc_available() -> bool:
    """Check if PyMC is available for MCMC sampling."""
    return PYMC_AVAILABLE


def run_mcmc_mixing(
    sample: IsotopeSample,
    sources: List[SourceIsotopes],
    n_samples: int = 2000,
    n_chains: int = 4,
    target_accept: float = 0.9,
    warmup: int = 1000,
    random_seed: Optional[int] = None,
    prior_alpha: Optional[List[float]] = None,
) -> MCMCMixingResult:
    """
    Run MCMC Bayesian mixing model for nitrate source apportionment.

    The model assumes the observed isotope signature is a linear mixture
    of source signatures:

        d15N_obs ~ N(sum(f_i * d15N_i), sigma_N)
        d18O_obs ~ N(sum(f_i * d18O_i), sigma_O)

    where f_i are source fractions (Dirichlet prior) and sigma are
    observation errors.

    Parameters
    ----------
    sample : IsotopeSample
        Observed isotope values (d15N, d18O)
    sources : List[SourceIsotopes]
        List of source endmember definitions
    n_samples : int
        Number of posterior samples per chain
    n_chains : int
        Number of MCMC chains
    target_accept : float
        Target acceptance rate for NUTS sampler
    warmup : int
        Number of warmup/burn-in samples
    random_seed : int, optional
        Random seed for reproducibility
    prior_alpha : List[float], optional
        Dirichlet prior concentrations (default: uniform)

    Returns
    -------
    MCMCMixingResult
        MCMC results with posterior samples and diagnostics
    """
    if not PYMC_AVAILABLE:
        raise ImportError(
            "PyMC is required for MCMC isotope mixing. "
            "Install with: pip install pymc>=5.0 arviz>=0.15"
        )

    if len(sources) < 2:
        raise ValueError("At least 2 sources required for mixing model.")

    n_sources = len(sources)
    source_names = [s.name for s in sources]

    # Default uniform Dirichlet prior
    if prior_alpha is None:
        prior_alpha = [1.0] * n_sources

    # Extract source parameters
    d15N_means = np.array([s.d15N_mean for s in sources])
    d15N_stds = np.array([s.d15N_std for s in sources])
    d18O_means = np.array([s.d18O_mean for s in sources])
    d18O_stds = np.array([s.d18O_std for s in sources])

    # Observed values
    obs_d15N = sample.d15N
    obs_d18O = sample.d18O

    warnings_list = []

    with pm.Model():
        # Prior on source fractions (Dirichlet ensures sum-to-one constraint)
        fractions = pm.Dirichlet("fractions", a=np.array(prior_alpha))

        # Predicted mixture signatures
        pred_d15N = pm.math.dot(fractions, d15N_means)
        pred_d18O = pm.math.dot(fractions, d18O_means)

        # Observation error (combining source uncertainty and measurement error)
        # Use pooled variance from sources as prior
        sigma_N = pm.HalfNormal("sigma_N", sigma=np.mean(d15N_stds))
        sigma_O = pm.HalfNormal("sigma_O", sigma=np.mean(d18O_stds))

        # Likelihood
        pm.Normal("obs_d15N", mu=pred_d15N, sigma=sigma_N, observed=obs_d15N)
        pm.Normal("obs_d18O", mu=pred_d18O, sigma=sigma_O, observed=obs_d18O)

        # Sample
        try:
            trace = pm.sample(
                draws=n_samples,
                tune=warmup,
                chains=n_chains,
                target_accept=target_accept,
                random_seed=random_seed,
                return_inferencedata=True,
                progressbar=False,
            )
        except Exception as e:
            warnings_list.append(f"Sampling error: {str(e)}")
            # Return fallback uniform result
            return MCMCMixingResult(
                source_fractions={name: 1.0 / n_sources for name in source_names},
                source_names=source_names,
                converged=False,
                warnings=warnings_list,
            )

    # Extract posterior samples
    posterior_samples = trace.posterior["fractions"].values
    # Shape: (n_chains, n_samples, n_sources) -> flatten chains
    posterior_flat = posterior_samples.reshape(-1, n_sources)

    # Compute point estimates (posterior means)
    means = posterior_flat.mean(axis=0)
    source_fractions = {name: float(means[i]) for i, name in enumerate(source_names)}

    # Compute 95% credible intervals
    ci_lower = {}
    ci_upper = {}
    for i, name in enumerate(source_names):
        ci_lower[name] = float(np.percentile(posterior_flat[:, i], 2.5))
        ci_upper[name] = float(np.percentile(posterior_flat[:, i], 97.5))

    # Diagnostics
    summary = az.summary(trace, var_names=["fractions"])

    r_hat = {}
    ess_bulk = {}
    ess_tail = {}

    for i, name in enumerate(source_names):
        idx = f"fractions[{i}]"
        if idx in summary.index:
            r_hat[name] = float(summary.loc[idx, "r_hat"])
            ess_bulk[name] = float(summary.loc[idx, "ess_bulk"])
            ess_tail[name] = float(summary.loc[idx, "ess_tail"])

    # Check convergence
    converged = True
    for name, rh in r_hat.items():
        if rh > 1.05:
            converged = False
            warnings_list.append(f"R-hat for {name} = {rh:.3f} > 1.05")

    for name, ess in ess_bulk.items():
        if ess < 400:
            warnings_list.append(f"Low ESS for {name}: {ess:.0f}")

    # Compute WAIC if available
    waic_value = None
    try:
        waic_result = az.waic(trace)
        waic_value = float(waic_result.elpd_waic)
    except Exception:
        pass

    return MCMCMixingResult(
        source_fractions=source_fractions,
        posterior_samples=posterior_flat,
        source_names=source_names,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        r_hat=r_hat,
        ess_bulk=ess_bulk,
        ess_tail=ess_tail,
        waic=waic_value,
        converged=converged,
        warnings=warnings_list,
    )


def run_mcmc_mixing_batch(
    samples: List[IsotopeSample],
    sources: List[SourceIsotopes],
    n_samples: int = 2000,
    n_chains: int = 4,
    target_accept: float = 0.9,
    warmup: int = 1000,
    random_seed: Optional[int] = None,
) -> List[MCMCMixingResult]:
    """
    Run MCMC mixing model for multiple samples.

    Parameters
    ----------
    samples : List[IsotopeSample]
        List of observed isotope samples
    sources : List[SourceIsotopes]
        List of source endmember definitions
    n_samples : int
        Number of posterior samples per chain
    n_chains : int
        Number of MCMC chains
    target_accept : float
        Target acceptance rate
    warmup : int
        Number of warmup samples
    random_seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    List[MCMCMixingResult]
        List of results for each sample
    """
    results = []
    for i, sample in enumerate(samples):
        seed = random_seed + i if random_seed is not None else None
        result = run_mcmc_mixing(
            sample=sample,
            sources=sources,
            n_samples=n_samples,
            n_chains=n_chains,
            target_accept=target_accept,
            warmup=warmup,
            random_seed=seed,
        )
        results.append(result)
    return results


def summarize_mcmc_results(results: List[MCMCMixingResult]) -> Dict[str, Any]:
    """
    Summarize MCMC results across multiple samples.

    Returns aggregated statistics including:
    - Mean source fractions across samples
    - Standard deviation of fractions
    - Fraction of samples with convergence issues
    """
    if not results:
        return {}

    source_names = results[0].source_names
    n_samples = len(results)

    # Collect fractions
    all_fractions = {name: [] for name in source_names}
    n_converged = 0
    all_warnings = []

    for res in results:
        for name in source_names:
            if name in res.source_fractions:
                all_fractions[name].append(res.source_fractions[name])
        if res.converged:
            n_converged += 1
        all_warnings.extend(res.warnings)

    summary = {
        "n_samples": n_samples,
        "n_converged": n_converged,
        "convergence_rate": n_converged / n_samples if n_samples > 0 else 0.0,
        "mean_fractions": {},
        "std_fractions": {},
        "warnings": list(set(all_warnings)),
    }

    for name in source_names:
        if all_fractions[name]:
            arr = np.array(all_fractions[name])
            summary["mean_fractions"][name] = float(np.mean(arr))
            summary["std_fractions"][name] = float(np.std(arr))

    return summary
