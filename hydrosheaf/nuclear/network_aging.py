"""
Network-Enhanced Bayesian Groundwater Dating.

This module implements a hierarchical Bayesian model to infer groundwater age distributions
across a flow network. It resolves the "post-bomb tritium ambiguity" by constraining
individual node ages with network topology (upstream -> downstream aging) and mixing logic.

Novelty:
- Treats age as a network field, not isolated point estimates.
- Uses Hamiltonian Monte Carlo (Nutpie/PyMC) to sample the joint posterior of all ages.
- Explicitly handles "mixing" vs "piston flow" via network structure.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import networkx as nx
import numpy as np

from hydrosheaf.nuclear.nuclides import Nuclide, TRITIUM
from hydrosheaf.nuclear.input_history import InputHistory, get_input_history
from hydrosheaf.nuclear.lpm import convolve_input
from ..log import get_logger

logger = get_logger(__name__)


def _load_network_aging_dependencies() -> Tuple[Any, Any, Any]:
    """Load optional Bayesian dependencies lazily."""
    try:
        import nutpie
        import pymc as pm
        import pytensor.tensor as pt
    except Exception as exc:
        raise ImportError(
            "Network aging requires pymc, pytensor, and nutpie with writable runtime directories."
        ) from exc
    return pm, pt, nutpie


def _linear_interp(x, x_ref, y_ref, pt_module: Any):
    """
    Differentiable linear interpolation in PyTensor.
    Assumes x_ref is sorted and effectively constant.
    """
    x_min = x_ref[0]
    x_max = x_ref[-1]
    x_clipped = pt_module.clip(x, x_min, x_max)
    
    dx = x_ref[1] - x_ref[0]
    idx_float = (x_clipped - x_min) / dx
    idx_floor = pt_module.floor(idx_float).astype("int64")
    idx_ceil = idx_floor + 1
    
    n_points = x_ref.shape[0]
    idx_floor = pt_module.clip(idx_floor, 0, n_points - 2)
    idx_ceil = pt_module.clip(idx_ceil, 1, n_points - 1)
    
    y_floor = y_ref[idx_floor]
    y_ceil = y_ref[idx_ceil]
    
    w_ceil = idx_float - idx_floor
    w_floor = 1.0 - w_ceil
    
    return y_floor * w_floor + y_ceil * w_ceil

def infer_network_ages_bayesian(
    graph: nx.DiGraph,
    node_observations: Dict[str, float],
    node_sigmas: Dict[str, float],
    sample_date: float,
    nuclide: Nuclide = TRITIUM,
    input_hist: Optional[InputHistory] = None,
    n_samples: int = 1000,
    n_chains: int = 4,
    random_seed: int = 42,
    model_type: str = "PFM",  # "PFM" (Piston) or "EM" (Exponential/Mean age)
    velocity_prior_mu: float = 10.0, # m/year
    velocity_prior_sigma: float = 50.0,
    input_scale_prior_mu: float = 1.0,
    input_scale_prior_sigma: float = 0.5,
    c14_initial_prior_mu: float = 85.0, # Standard approx for carbonates
    c14_initial_prior_sigma: float = 10.0,
) -> Dict[str, Dict[str, float]]:
    """
    Infer groundwater mean residence times (ages) for an entire network simultaneously.
    
    The network model supports multiple transit-time assumptions, input-history
    scaling, velocity-distance priors, and optional joint inference of the
    radiocarbon dead-carbon fraction.
    """
    pm, pt, nutpie = _load_network_aging_dependencies()
    
    if input_hist is None:
        input_hist = get_input_history("global")
        
    lambda_y = np.log(2) / nuclide.half_life_years
    nodes = list(graph.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    n_nodes = len(nodes)
    
    # Identify observed nodes
    obs_indices = []
    obs_values = []
    obs_errors = []
    
    for n in nodes:
        if n in node_observations:
            obs_indices.append(node_idx[n])
            val = max(node_observations[n], 0.0) 
            obs_values.append(val)
            obs_errors.append(node_sigmas.get(n, max(0.5, val*0.1)))
            
    obs_indices = np.array(obs_indices, dtype=int)
    obs_values = np.array(obs_values, dtype=float)
    obs_errors = np.array(obs_errors, dtype=float)
    
    # Build Look-up Table for Forward Model
    max_age = 500.0
    age_grid = np.linspace(0, max_age, 250)
    conc_grid = []
    for t in age_grid:
        c = convolve_input(sample_date, t, input_hist.years, input_hist.values, lambda_y, model_type=model_type)
        conc_grid.append(c)
    conc_grid = np.array(conc_grid)
    
    # Extract edges with lengths
    us, vs, lens = [], [], []
    for u, v, data in graph.edges(data=True):
        if u in node_idx and v in node_idx:
            us.append(node_idx[u])
            vs.append(node_idx[v])
            lens.append(data.get("length_m", 1.0))
    
    edge_us = np.array(us, dtype=int)
    edge_vs = np.array(vs, dtype=int)
    edge_lens = np.array(lens, dtype=float)

    logger.info(f"Starting Bayesian Network Dating for {n_nodes} nodes, {len(edge_us)} edges.")
    logger.info(f"Settings: n_samples={n_samples}, nuclide={nuclide.symbol}, model_type={model_type}")

    with pm.Model() as model:
        # 1. Latent Variables
        # Input Scaling (Accounts for regional variability / Tritium scaling)
        input_scale = pm.LogNormal("input_scale", mu=np.log(input_scale_prior_mu), sigma=input_scale_prior_sigma)
        
        # C14 Geochemical Correction (Hierarchical Dead Carbon Fraction)
        # Allows for local geochemical heterogeneity (Partial Pooling)
        if nuclide.symbol == "14C":
            # Global mean A0 for the aquifer system
            c14_a0_mu = pm.Normal("c14_a0_mu", mu=c14_initial_prior_mu, sigma=c14_initial_prior_sigma)
            
            # Local A0 per node (constrained by global mean)
            # sigma=5.0 represents local geochemical variability (approx +/- 10 pmc range)
            c14_a0_local = pm.Normal("c14_a0_local", mu=c14_a0_mu, sigma=5.0, shape=n_nodes)
        else:
            # Placeholder to avoid "variable not found" errors in non-C14 runs
            c14_a0_local = None
        
        # Velocity (Accounts for physical travel time)
        velocity_mu = max(float(velocity_prior_mu), 1e-6)
        velocity_sigma = max(float(velocity_prior_sigma), 1e-6)
        velocity_sigma_log = float(np.log1p(velocity_sigma / velocity_mu))
        velocity = pm.LogNormal(
            "velocity", mu=np.log(velocity_mu), sigma=velocity_sigma_log
        )

        # Node Ages.
        # NOTE (identifiability caveat): for tracers that are uninformative at old
        # age (dead tritium beyond ~60-70 yr), node age is only weakly identified —
        # the joint velocity/input-scale/age posterior is degenerate and the model
        # cannot recover accurate old ages from tritium alone. Bayesian network
        # dating is therefore reliable in the tracer-INFORMATIVE / ambiguity-
        # resolution regime; old-water accuracy requires an independent velocity
        # (hydraulic) constraint or a longer-lived tracer (14C). Ordering-based
        # alias resolution (see M3 run_m3_network_dating_demo) is the robust
        # network-dating method for the tritium bomb-peak regime.
        log_ages = pm.Normal("log_ages", mu=np.log(20), sigma=2.0, shape=n_nodes)
        ages = pm.Deterministic("ages", pm.math.exp(log_ages))

        # 2. Physics/Topology Constraint
        if edge_us.size:
            age_u = ages[edge_us]
            age_v = ages[edge_vs]
            expected_dt = edge_lens / velocity
            diffs = age_v - (age_u + expected_dt)
            pm.Potential("network_flow_constraint", -20.0 * pm.math.sum(pm.math.switch(diffs < 0, (diffs)**2, 0)))

        # 3. Observation Likelihood
        grid_x_pt = pt.as_tensor_variable(age_grid)
        grid_y_pt = pt.as_tensor_variable(conc_grid)
        c_base = _linear_interp(ages[obs_indices], grid_x_pt, grid_y_pt, pt)
        
        # Combine base concentration with latent scaling/corrections
        if nuclide.symbol == "14C":
            # For C14: Expected = Local_A0 * exp(-lambda * tau)
            # Use the specific A0 for each observed node
            a0_obs = c14_a0_local[obs_indices]
            c_expected = c_base * (a0_obs / 100.0)
        else:
            # For 3H: Expected = Lookup(tau) * regional_scale
            c_expected = c_base * input_scale
        
        pm.Normal("obs_likelihood", mu=c_expected, sigma=obs_errors, observed=obs_values)
        
        # 4. Sampling
        compiled = nutpie.compile_pymc_model(model)
        trace = nutpie.sample(compiled, draws=n_samples, chains=n_chains, seed=random_seed, progress_bar=False)
        
    # Process
    results = {}
    flat_ages = trace.posterior["ages"].values.reshape(-1, n_nodes)
    
    # Diagnostics
    diag = {
        "inferred_input_scale": float(np.mean(trace.posterior["input_scale"])) if "input_scale" in trace.posterior else 1.0,
        "inferred_velocity_m_y": float(np.mean(trace.posterior["velocity"])),
        "model_type": model_type
    }
    try:
        import arviz as az

        age_summary = az.summary(
            trace,
            var_names=["ages"],
            round_to=None,
        )
        r_hat_values = np.asarray(age_summary["r_hat"], dtype=float)
        ess_bulk_values = np.asarray(age_summary["ess_bulk"], dtype=float)
        ess_tail_values = np.asarray(age_summary["ess_tail"], dtype=float)
        mcse_mean_values = np.asarray(age_summary["mcse_mean"], dtype=float)
        diag.update(
            {
                "age_r_hat_max": float(np.nanmax(r_hat_values)),
                "age_ess_bulk_min": float(np.nanmin(ess_bulk_values)),
                "age_ess_tail_min": float(np.nanmin(ess_tail_values)),
                "age_mcse_mean_max": float(np.nanmax(mcse_mean_values)),
            }
        )
    except Exception as exc:
        logger.warning("Could not compute Bayesian-age diagnostics: %s", exc)
        diag.update(
            {
                "age_r_hat_max": None,
                "age_ess_bulk_min": None,
                "age_ess_tail_min": None,
                "age_mcse_mean_max": None,
            }
        )
    sample_stats = getattr(trace, "sample_stats", None)
    divergences = None
    if sample_stats is not None:
        for key in ("diverging", "divergence"):
            if key in sample_stats:
                divergences = int(np.asarray(sample_stats[key]).sum())
                break
    diag["divergences"] = divergences
    if "c14_a0_mu" in trace.posterior:
        diag["inferred_c14_a0_global"] = float(np.mean(trace.posterior["c14_a0_mu"]))
        
    logger.info(f"MCMC Finished. Diagnostics: {diag}")
    results["_diagnostics"] = diag
    
    for i, n in enumerate(nodes):
        node_ages = flat_ages[:, i]
        results[n] = {
            "mean_age_years": float(np.mean(node_ages)),
            "std_age_years": float(np.std(node_ages)),
            "p_modern": float(np.mean(node_ages < 60.0)),
            "age_95_low": float(np.percentile(node_ages, 2.5)),
            "age_95_high": float(np.percentile(node_ages, 97.5))
        }
        
    return results
