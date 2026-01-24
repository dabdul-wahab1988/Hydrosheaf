"""Bayesian and heuristic head estimation utilities.

This module supports "hierarchical head estimation" by combining multiple
head-related data sources:
  A) measured hydraulic head (head_key)
  B) elevation - depth-to-water (dtw_key + elevation_key)
  C) topographic proxy (elevation_key) with large uncertainty (sigma_topo)

The Bayesian implementation uses a linear-Gaussian hierarchical model:
  - latent head h_i per node
  - shared mean depth-to-water parameter mu_dtw
  - topographic prior: h_i + mu_dtw ≈ elevation_i
  - observations: head_meas and/or (elevation - dtw)

Posterior inference is closed-form (multivariate normal), computed via a
precision-matrix solve (no PyMC required).
"""

from __future__ import annotations

from dataclasses import dataclass

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pymc as pm
import sys
import os

if sys.platform == "win32":
    # Fix for PyMC/MKL crash on Windows (Heap Corruption)
    # Must be set before pymc/numpy BLAS initialization
    os.environ["MKL_THREADING_LAYER"] = "GNU"

@dataclass(frozen=True)
class HeadPosterior:
    node_ids: List[str]
    head_mean: np.ndarray  # shape: (n,)
    head_cov: np.ndarray  # shape: (n, n)
    tiers: List[str]
    mu_dtw_mean: float
    mu_dtw_sigma: float

    def index(self) -> Dict[str, int]:
        return {node_id: idx for idx, node_id in enumerate(self.node_ids)}


def _present(sample: Mapping[str, object], key: str) -> bool:
    return key in sample and sample.get(key) not in (None, "")


def _parse_float(sample: Mapping[str, object], key: str) -> Optional[float]:
    value = sample.get(key)
    if value in (None, ""):
        return None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _tier_label(
    sample: Mapping[str, object],
    head_key: str,
    dtw_key: str,
    elevation_key: str,
) -> str:
    tiers: List[str] = []
    if _present(sample, head_key):
        tiers.append("A")
    if _present(sample, dtw_key) and _present(sample, elevation_key):
        tiers.append("B")
    elif _present(sample, elevation_key):
        tiers.append("C")
    return "+".join(tiers) if tiers else "missing"


def infer_heads_bayesian_linear(
    samples: Sequence[Mapping[str, object]],
    *,
    node_id_key: str,
    head_key: str,
    dtw_key: str,
    elevation_key: str,
    sigma_meas: float,
    sigma_dtw: float,
    sigma_elev: float,
    sigma_topo: float,
    dtw_prior_mu: float,
    dtw_prior_sigma: float,
    head_prior_mu: float,
    head_prior_sigma: float,
    jitter: float = 1e-9,
) -> HeadPosterior:
    """Infer a posterior over heads using a linear-Gaussian hierarchical model.

    Nodes with no usable node_id are ignored. Nodes with no elevation and no
    head observations are ignored (uninformative for flow direction inference).
    """
    if sigma_meas <= 0 or sigma_dtw <= 0 or sigma_elev <= 0 or sigma_topo <= 0:
        raise ValueError("All sigma values must be positive.")
    if dtw_prior_sigma <= 0 or head_prior_sigma <= 0:
        raise ValueError("Prior sigmas must be positive.")

    node_ids: List[str] = []
    elevations: List[Optional[float]] = []
    observations: List[List[Tuple[float, float]]] = []
    tiers: List[str] = []

    for sample in samples:
        node_id = sample.get(node_id_key)
        if node_id in (None, ""):
            continue
        node_id_str = str(node_id)

        elev = _parse_float(sample, elevation_key)
        head_obs = _parse_float(sample, head_key)
        dtw = _parse_float(sample, dtw_key)

        obs: List[Tuple[float, float]] = []
        if head_obs is not None:
            obs.append((head_obs, sigma_meas))
        if elev is not None and dtw is not None:
            obs.append((elev - dtw, float(np.sqrt(sigma_elev**2 + sigma_dtw**2))))

        tier = _tier_label(sample, head_key, dtw_key, elevation_key)

        # If there is no elevation prior and no observation, the node is not
        # identifiable and would yield an arbitrary direction probability.
        if elev is None and not obs:
            continue

        node_ids.append(node_id_str)
        elevations.append(elev)
        observations.append(obs)
        tiers.append(tier)

    n = len(node_ids)
    if n == 0:
        return HeadPosterior(
            node_ids=[],
            head_mean=np.array([]),
            head_cov=np.zeros((0, 0)),
            tiers=[],
            mu_dtw_mean=float(dtw_prior_mu),
            mu_dtw_sigma=float(dtw_prior_sigma),
        )

    # x = [h_1..h_n, mu_dtw]
    q = np.zeros((n + 1, n + 1), dtype=float)
    b = np.zeros((n + 1,), dtype=float)

    # Prior: mu_dtw ~ Normal(dtw_prior_mu, dtw_prior_sigma^2)
    mu_idx = n
    mu_prec = 1.0 / (dtw_prior_sigma**2)
    q[mu_idx, mu_idx] += mu_prec
    b[mu_idx] += dtw_prior_mu * mu_prec

    topo_prec = 1.0 / (sigma_topo**2)
    head_prior_prec = 1.0 / (head_prior_sigma**2)

    for idx, (elev, obs_list) in enumerate(zip(elevations, observations)):
        # Topographic prior (if elevation is available):
        #   h_i + mu_dtw ≈ elevation_i  with sigma_topo
        if elev is not None:
            q[idx, idx] += topo_prec
            q[mu_idx, mu_idx] += topo_prec
            q[idx, mu_idx] += topo_prec
            q[mu_idx, idx] += topo_prec
            b[idx] += elev * topo_prec
            b[mu_idx] += elev * topo_prec
        else:
            # Weak independent prior on h_i if there is no elevation.
            q[idx, idx] += head_prior_prec
            b[idx] += head_prior_mu * head_prior_prec

        # Observations of h_i
        for value, sigma in obs_list:
            obs_prec = 1.0 / (sigma**2)
            q[idx, idx] += obs_prec
            b[idx] += value * obs_prec

    if jitter > 0:
        q += np.eye(n + 1) * float(jitter)

    mean = np.linalg.solve(q, b)
    cov = np.linalg.inv(q)

    head_mean = mean[:n]
    head_cov = cov[:n, :n]
    mu_mean = float(mean[mu_idx])
    mu_sigma = float(np.sqrt(cov[mu_idx, mu_idx]))

    return HeadPosterior(
        node_ids=node_ids,
        head_mean=head_mean,
        head_cov=head_cov,
        tiers=tiers,
        mu_dtw_mean=mu_mean,
        mu_dtw_sigma=mu_sigma,
    )


def infer_heads_bayesian_mcmc(
    samples: Sequence[Mapping[str, object]],
    *,
    node_id_key: str,
    head_key: str,
    dtw_key: str,
    elevation_key: str,
    sigma_meas: float,
    sigma_dtw: float,
    sigma_elev: float,
    sigma_topo: float,
    dtw_prior_mu: float,
    dtw_prior_sigma: float,
    head_prior_mu: float,
    head_prior_sigma: float,
    mcmc_draws: int = 1000,
    mcmc_chains: int = 2,
    mcmc_target_accept: float = 0.9,
    mcmc_warmup_fraction: float = 0.5,
    random_seed: Optional[int] = None,
    cores: int = 1,
) -> HeadPosterior:
    """Run MCMC with PyMC if available; otherwise fall back to closed form.

    Returns a HeadPosterior whose mean/cov are computed from posterior samples
    when PyMC is available. If PyMC import fails, this uses
    infer_heads_bayesian_linear instead.
    """
    # Force single core on Windows to avoid multiprocessing crashes (0xc0000374)
    if sys.platform == "win32":
        cores = 1



    # Build a consistent node list and observation arrays
    node_ids: List[str] = []
    elevations: List[Optional[float]] = []
    head_obs_vals: List[float] = []
    head_obs_idx: List[int] = []
    dtw_obs_vals: List[float] = []
    dtw_obs_idx: List[int] = []
    tiers: List[str] = []

    rows: List[Tuple[str, Mapping[str, object]]] = []
    for sample in samples:
        node_id = sample.get(node_id_key)
        if node_id in (None, ""):
            continue
        rows.append((str(node_id), sample))

    for node_id, sample in rows:
        elev = _parse_float(sample, elevation_key)
        head_obs = _parse_float(sample, head_key)
        dtw = _parse_float(sample, dtw_key)

        # Ignore nodes with no usable signal at all
        if (
            elev is None
            and head_obs is None
            and not (elev is not None and dtw is not None)
        ):
            continue

        idx = len(node_ids)
        node_ids.append(node_id)
        elevations.append(elev)
        tiers.append(_tier_label(sample, head_key, dtw_key, elevation_key))

        if head_obs is not None:
            head_obs_idx.append(idx)
            head_obs_vals.append(head_obs)
        if elev is not None and dtw is not None:
            dtw_obs_idx.append(idx)
            dtw_obs_vals.append(elev - dtw)

    n = len(node_ids)
    if n == 0:
        return HeadPosterior(
            node_ids=[],
            head_mean=np.array([]),
            head_cov=np.zeros((0, 0)),
            tiers=[],
            mu_dtw_mean=float(dtw_prior_mu),
            mu_dtw_sigma=float(dtw_prior_sigma),
        )

    elev_idx = [i for i, e in enumerate(elevations) if e is not None]
    elev_vals = [float(elevations[i]) for i in elev_idx]

    head_obs_vals_arr = np.asarray(head_obs_vals, dtype=float)
    head_obs_idx_arr = np.asarray(head_obs_idx, dtype=int)
    dtw_obs_vals_arr = np.asarray(dtw_obs_vals, dtype=float)
    dtw_obs_idx_arr = np.asarray(dtw_obs_idx, dtype=int)
    elev_vals_arr = np.asarray(elev_vals, dtype=float)
    elev_idx_arr = np.asarray(elev_idx, dtype=int)

    tune = int(max(0, round(float(mcmc_draws) * float(mcmc_warmup_fraction))))
    draws = int(max(10, mcmc_draws))
    chains = int(max(1, mcmc_chains))

    with pm.Model():
        mu_dtw = pm.Normal(
            "mu_dtw", mu=float(dtw_prior_mu), sigma=float(dtw_prior_sigma)
        )
        h = pm.Normal(
            "h", mu=float(head_prior_mu), sigma=float(head_prior_sigma), shape=n
        )

        # Topographic prior as an observed relationship when elevation exists:
        # elevation_i ~ Normal(h_i + mu_dtw, sigma_topo)
        if elev_idx_arr.size:
            pm.Normal(
                "elevation_obs",
                mu=h[elev_idx_arr] + mu_dtw,
                sigma=float(sigma_topo),
                observed=elev_vals_arr,
            )

        # Direct head measurements
        if head_obs_idx_arr.size:
            pm.Normal(
                "head_obs",
                mu=h[head_obs_idx_arr],
                sigma=float(sigma_meas),
                observed=head_obs_vals_arr,
            )

        # Derived head observations from elevation - dtw
        if dtw_obs_idx_arr.size:
            pm.Normal(
                "dtw_head_obs",
                mu=h[dtw_obs_idx_arr],
                sigma=float(np.sqrt(sigma_elev**2 + sigma_dtw**2)),
                observed=dtw_obs_vals_arr,
            )

        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=float(mcmc_target_accept),
            random_seed=random_seed,
            cores=int(max(1, cores)),
            progressbar=False,
            compute_convergence_checks=False,
        )

    # Extract posterior draws and compute mean/cov
    h_draws = idata.posterior["h"].values  # (chain, draw, n)
    h_samples = h_draws.reshape(-1, n)
    head_mean = np.mean(h_samples, axis=0)
    head_cov = np.cov(h_samples, rowvar=False)

    mu_draws = idata.posterior["mu_dtw"].values.reshape(-1)
    mu_mean = float(np.mean(mu_draws))
    mu_sigma = float(np.std(mu_draws, ddof=1)) if mu_draws.size > 1 else 0.0

    return HeadPosterior(
        node_ids=node_ids,
        head_mean=head_mean,
        head_cov=head_cov,
        tiers=tiers,
        mu_dtw_mean=mu_mean,
        mu_dtw_sigma=mu_sigma,
    )


def infer_head_heuristic(
    sample: Mapping[str, object],
    *,
    head_key: str,
    dtw_key: str,
    elevation_key: str,
    sigma_meas: float,
    sigma_dtw: float,
    sigma_elev: float,
    sigma_topo: float,
) -> Tuple[Optional[float], Optional[float], str]:
    """Return (head_estimate, sigma, tier) using the existing heuristic tiers."""
    if _present(sample, head_key):
        return float(sample[head_key]), sigma_meas, "A"  # type: ignore[arg-type]
    if _present(sample, dtw_key) and _present(sample, elevation_key):
        head = float(sample[elevation_key]) - float(sample[dtw_key])  # type: ignore[arg-type]
        sigma = float(np.sqrt(sigma_elev**2 + sigma_dtw**2))
        return head, sigma, "B"
    if _present(sample, elevation_key):
        return float(sample[elevation_key]), sigma_topo, "C"  # type: ignore[arg-type]
    return None, None, "missing"
