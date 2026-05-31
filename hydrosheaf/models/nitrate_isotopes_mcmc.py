"""
MCMC Bayesian Isotope Mixing Model for Nitrate Source Apportionment.

This module provides full MCMC posterior sampling using PyMC for quantifying
uncertainty in nitrate source fractions based on dual isotope (d15N, d18O) data.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from .nitrate_isotopes import IsotopeSample, SourceIsotopes, compute_isotope_prob



def _load_mcmc_dependencies():
    """Load optional MCMC dependencies lazily."""
    import numba  # noqa: F401 - required backend dependency for nutpie
    import pymc as pm
    import arviz as az
    import nutpie

    return pm, az, nutpie


def check_pymc_available() -> bool:
    """Check if PyMC and Nutpie are available and working."""
    try:
        _load_mcmc_dependencies()
        return True
    except Exception:
        return False

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


@dataclass
class HierarchicalMCMCMixingResult:
    """Result of hierarchical multi-sample isotope mixing analysis."""

    sample_results: Dict[str, MCMCMixingResult]
    global_prior_mean: Dict[str, float]
    source_names: List[str] = field(default_factory=list)
    converged: bool = True
    warnings: List[str] = field(default_factory=list)


def _extract_source_parameters(sources: List[SourceIsotopes]) -> Dict[str, np.ndarray]:
    d15N_means = np.array([s.d15N_mean for s in sources], dtype=float)
    d15N_stds = np.array([max(s.d15N_std, 1e-12) for s in sources], dtype=float)
    d18O_means = np.array([s.d18O_mean for s in sources], dtype=float)
    d18O_stds = np.array([max(s.d18O_std, 1e-12) for s in sources], dtype=float)
    d15N_vars = np.square(d15N_stds)
    d18O_vars = np.square(d18O_stds)
    source_cov_no = []
    for src, var_n, var_o in zip(sources, d15N_vars, d18O_vars):
        cov = float(src.covariance_d15N_d18O())
        max_abs_cov = float(np.sqrt(var_n * var_o) * 0.999999)
        source_cov_no.append(max(-max_abs_cov, min(max_abs_cov, cov)))
    source_cov_no = np.array(source_cov_no, dtype=float)
    source_means = np.column_stack([d15N_means, d18O_means])
    return {
        "d15N_means": d15N_means,
        "d15N_stds": d15N_stds,
        "d18O_means": d18O_means,
        "d18O_stds": d18O_stds,
        "d15N_vars": d15N_vars,
        "d18O_vars": d18O_vars,
        "source_cov_no": source_cov_no,
        "source_means": source_means,
    }


def _normalize_prior_vector(prior_values: np.ndarray) -> np.ndarray:
    clipped = np.maximum(np.asarray(prior_values, dtype=float), 1e-12)
    total = float(np.sum(clipped))
    if total <= 0.0:
        return np.ones_like(clipped) / float(len(clipped))
    return clipped / total


def _empirical_bayes_hierarchical(
    samples: List[IsotopeSample],
    sources: List[SourceIsotopes],
    sample_ids: List[str],
    prior_alpha: Optional[List[float]] = None,
    n_steps: int = 2,
    warnings_list: Optional[List[str]] = None,
) -> HierarchicalMCMCMixingResult:
    source_names = [s.name for s in sources]
    n_sources = len(sources)
    if warnings_list is None:
        warnings_list = []

    if prior_alpha is None:
        prior_probs = np.ones(n_sources, dtype=float) / float(n_sources)
    else:
        prior_probs = _normalize_prior_vector(np.array(prior_alpha, dtype=float))

    posterior_by_sample: List[Dict[str, float]] = []
    n_steps = max(int(n_steps), 1)
    for _ in range(n_steps):
        posterior_by_sample = [
            compute_isotope_prob(
                sample,
                sources,
                prior_probs={name: float(prior_probs[i]) for i, name in enumerate(source_names)},
            )
            for sample in samples
        ]
        pooled = np.zeros(n_sources, dtype=float)
        for probs in posterior_by_sample:
            for idx, name in enumerate(source_names):
                pooled[idx] += float(probs.get(name, 0.0))
        prior_probs = _normalize_prior_vector(pooled)

    sample_results: Dict[str, MCMCMixingResult] = {}
    for sample_id, probs in zip(sample_ids, posterior_by_sample):
        sample_results[sample_id] = MCMCMixingResult(
            source_fractions=probs,
            source_names=source_names,
            converged=False,
            warnings=["Hierarchical MCMC fallback used empirical Bayes approximation."],
        )

    warnings_out = list(warnings_list)
    warnings_out.append("Hierarchical MCMC fallback used empirical Bayes approximation.")
    global_prior = {name: float(prior_probs[i]) for i, name in enumerate(source_names)}
    return HierarchicalMCMCMixingResult(
        sample_results=sample_results,
        global_prior_mean=global_prior,
        source_names=source_names,
        converged=False,
        warnings=warnings_out,
    )


def run_mcmc_mixing_hierarchical(
    samples: List[IsotopeSample],
    sources: List[SourceIsotopes],
    sample_ids: Optional[List[str]] = None,
    n_samples: int = 2000,
    n_chains: int = 4,
    target_accept: float = 0.9,
    warmup: int = 1000,
    random_seed: Optional[int] = None,
    prior_alpha: Optional[List[float]] = None,
) -> HierarchicalMCMCMixingResult:
    """Run hierarchical multi-sample MCMC with shared source priors."""
    if not samples:
        raise ValueError("At least one sample is required for hierarchical mixing.")
    if len(sources) < 2:
        raise ValueError("At least 2 sources required for hierarchical mixing model.")

    n_obs = len(samples)
    n_sources = len(sources)
    source_names = [s.name for s in sources]
    if sample_ids is None:
        sample_ids = [f"sample_{idx + 1}" for idx in range(n_obs)]
    if len(sample_ids) != n_obs:
        raise ValueError("sample_ids length must match number of samples.")

    obs = np.array([[float(s.d15N), float(s.d18O)] for s in samples], dtype=float)
    if prior_alpha is None:
        prior_alpha_arr = np.ones(n_sources, dtype=float)
    else:
        prior_alpha_arr = np.array(prior_alpha, dtype=float)
        if prior_alpha_arr.size != n_sources:
            raise ValueError("prior_alpha length must match number of sources.")
    prior_alpha_arr = np.maximum(prior_alpha_arr, 1e-3)

    warnings_list: List[str] = []
    try:
        pm, az, nutpie = _load_mcmc_dependencies()
    except Exception as exc:
        warnings_list.append(f"MCMC dependencies unavailable: {exc}")
        return _empirical_bayes_hierarchical(
            samples=samples,
            sources=sources,
            sample_ids=sample_ids,
            prior_alpha=prior_alpha,
            warnings_list=warnings_list,
        )

    source_params = _extract_source_parameters(sources)
    d15N_stds = source_params["d15N_stds"]
    d18O_stds = source_params["d18O_stds"]
    d15N_vars = source_params["d15N_vars"]
    d18O_vars = source_params["d18O_vars"]
    source_cov_no = source_params["source_cov_no"]
    source_means = source_params["source_means"]

    log_two_pi = float(np.log(2.0 * np.pi))

    with pm.Model():
        alpha_shared = pm.Gamma("alpha_shared", alpha=prior_alpha_arr + 1.0, beta=1.0)
        fractions = pm.Dirichlet("fractions", a=alpha_shared, shape=(n_obs, n_sources))
        pred_isotopes = pm.math.dot(fractions, source_means)

        sigma_N = pm.HalfNormal("sigma_N", sigma=max(float(np.mean(d15N_stds)), 1e-6))
        sigma_O = pm.HalfNormal("sigma_O", sigma=max(float(np.mean(d18O_stds)), 1e-6))

        fractions_sq = fractions**2
        mix_var_n = pm.math.dot(fractions_sq, d15N_vars) + sigma_N**2
        mix_var_o = pm.math.dot(fractions_sq, d18O_vars) + sigma_O**2
        mix_cov_no = pm.math.dot(fractions_sq, source_cov_no)

        det = pm.math.maximum(mix_var_n * mix_var_o - mix_cov_no**2, 1e-12)
        dx_n = obs[:, 0] - pred_isotopes[:, 0]
        dx_o = obs[:, 1] - pred_isotopes[:, 1]
        quad = (
            mix_var_o * dx_n**2 - 2.0 * mix_cov_no * dx_n * dx_o + mix_var_n * dx_o**2
        ) / det
        logp = -0.5 * (2.0 * log_two_pi + pm.math.log(det) + quad)
        pm.Potential("obs_logp", pm.math.sum(logp))

        try:
            compiled_model = nutpie.compile_pymc_model(pm.Model.get_context())
            trace = nutpie.sample(
                compiled_model,
                draws=n_samples,
                tune=warmup,
                chains=n_chains,
                seed=random_seed,
                progress_bar=False,
            )
        except Exception as exc:
            warnings_list.append(f"Sampling error: {exc}")
            return _empirical_bayes_hierarchical(
                samples=samples,
                sources=sources,
                sample_ids=sample_ids,
                prior_alpha=prior_alpha,
                warnings_list=warnings_list,
            )

    posterior = trace.posterior["fractions"].values.reshape(-1, n_obs, n_sources)
    alpha_samples = trace.posterior["alpha_shared"].values.reshape(-1, n_sources)
    alpha_mean = np.maximum(alpha_samples.mean(axis=0), 1e-12)
    alpha_prob = _normalize_prior_vector(alpha_mean)
    global_prior = {name: float(alpha_prob[idx]) for idx, name in enumerate(source_names)}

    summary = az.summary(trace, var_names=["fractions", "alpha_shared"])
    sample_results: Dict[str, MCMCMixingResult] = {}
    global_converged = True

    for obs_idx, sample_id in enumerate(sample_ids):
        draws = posterior[:, obs_idx, :]
        means = draws.mean(axis=0)
        fractions_mean = {
            source_names[idx]: float(means[idx]) for idx in range(n_sources)
        }
        ci_lower = {
            source_names[idx]: float(np.percentile(draws[:, idx], 2.5))
            for idx in range(n_sources)
        }
        ci_upper = {
            source_names[idx]: float(np.percentile(draws[:, idx], 97.5))
            for idx in range(n_sources)
        }
        r_hat: Dict[str, float] = {}
        ess_bulk: Dict[str, float] = {}
        ess_tail: Dict[str, float] = {}
        row_warnings: List[str] = []
        row_converged = True
        for src_idx, src_name in enumerate(source_names):
            summary_idx = f"fractions[{obs_idx}, {src_idx}]"
            if summary_idx not in summary.index:
                summary_idx = f"fractions[{obs_idx},{src_idx}]"
            if summary_idx in summary.index:
                r_value = float(summary.loc[summary_idx, "r_hat"])
                b_value = float(summary.loc[summary_idx, "ess_bulk"])
                t_value = float(summary.loc[summary_idx, "ess_tail"])
                r_hat[src_name] = r_value
                ess_bulk[src_name] = b_value
                ess_tail[src_name] = t_value
                if r_value > 1.05:
                    row_converged = False
                    row_warnings.append(f"R-hat for {sample_id}:{src_name} = {r_value:.3f}")
                if b_value < 400:
                    row_warnings.append(
                        f"Low ESS for {sample_id}:{src_name}: {b_value:.0f}"
                    )

        if not row_converged:
            global_converged = False
            warnings_list.extend(row_warnings)

        sample_results[sample_id] = MCMCMixingResult(
            source_fractions=fractions_mean,
            posterior_samples=draws,
            source_names=source_names,
            ci_lower=ci_lower,
            ci_upper=ci_upper,
            r_hat=r_hat,
            ess_bulk=ess_bulk,
            ess_tail=ess_tail,
            converged=row_converged,
            warnings=row_warnings,
        )

    return HierarchicalMCMCMixingResult(
        sample_results=sample_results,
        global_prior_mean=global_prior,
        source_names=source_names,
        converged=global_converged,
        warnings=list(dict.fromkeys(warnings_list)),
    )


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
    of source signatures with covariance-aware uncertainties:

        [d15N_obs, d18O_obs] ~ MVN(mu_mix, Sigma_mix)

    where mu_mix is the fraction-weighted source mean vector and Sigma_mix
    includes fraction-weighted source covariance plus residual noise terms.

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
    if len(sources) < 2:
        raise ValueError("At least 2 sources required for mixing model.")


    n_sources = len(sources)
    source_names = [s.name for s in sources]
    warnings_list = []

    try:
        pm, az, nutpie = _load_mcmc_dependencies()
    except Exception as exc:
        warnings_list.append(f"MCMC dependencies unavailable: {exc}")
        return MCMCMixingResult(
            source_fractions={name: 1.0 / n_sources for name in source_names},
            source_names=source_names,
            converged=False,
            warnings=warnings_list,
        )

    # Default uniform Dirichlet prior
    if prior_alpha is None:
        prior_alpha = [1.0] * n_sources

    # Extract source parameters
    source_params = _extract_source_parameters(sources)
    d15N_stds = source_params["d15N_stds"]
    d18O_stds = source_params["d18O_stds"]
    d15N_vars = source_params["d15N_vars"]
    d18O_vars = source_params["d18O_vars"]
    source_cov_no = source_params["source_cov_no"]
    source_means = source_params["source_means"]

    # Observed values
    obs_d15N = sample.d15N
    obs_d18O = sample.d18O

    with pm.Model():
        # Prior on source fractions (Dirichlet ensures sum-to-one constraint)
        fractions = pm.Dirichlet("fractions", a=np.array(prior_alpha))

        # Predicted mixture signatures
        pred_isotopes = pm.math.dot(fractions, source_means)

        # Residual observation noise on top of source uncertainty.
        sigma_N = pm.HalfNormal("sigma_N", sigma=max(float(np.mean(d15N_stds)), 1e-6))
        sigma_O = pm.HalfNormal("sigma_O", sigma=max(float(np.mean(d18O_stds)), 1e-6))

        fractions_sq = fractions**2
        mix_var_n = pm.math.dot(fractions_sq, d15N_vars) + sigma_N**2
        mix_var_o = pm.math.dot(fractions_sq, d18O_vars) + sigma_O**2
        mix_cov_no = pm.math.dot(fractions_sq, source_cov_no)

        cov_row_1 = pm.math.stack([mix_var_n, mix_cov_no])
        cov_row_2 = pm.math.stack([mix_cov_no, mix_var_o])
        mix_covariance = pm.math.stack([cov_row_1, cov_row_2])

        pm.MvNormal(
            "obs_isotopes",
            mu=pred_isotopes,
            cov=mix_covariance,
            observed=np.array([obs_d15N, obs_d18O], dtype=float),
        )

        # Sample
        try:
            # Compile and sample using Nutpie (Numba-based NUTS sampler)
            compiled_model = nutpie.compile_pymc_model(pm.Model.get_context())
            
            trace = nutpie.sample(
                compiled_model,
                draws=n_samples,
                tune=warmup,
                chains=n_chains,
                seed=random_seed,
                progress_bar=False,
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
