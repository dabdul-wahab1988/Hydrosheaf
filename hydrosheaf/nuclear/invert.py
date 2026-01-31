"""
Inverse modeling for radioactive tracers.
"""
from __future__ import annotations

import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from .nuclides import Nuclide, TRITIUM
from .input_history import InputHistory, get_input_history
from .lpm import convolve_input
from ..log import get_logger

logger = get_logger(__name__)

def infer_age_from_tracer(
    obs_value: float,
    obs_sigma: float,
    sample_date_year: float,
    nuclide: Nuclide,
    input_hist: Optional[InputHistory] = None,
    model: str = "PFM",
    tau_range: Tuple[float, float] = (0.0, 200.0),
    search_steps: int = 200,
    param2: float = 0.05,  # e.g., dispersion parameter for DM
    censoring_limit: Optional[float] = None,
    initial_activity: float = 100.0, # For C14: A0 (in pmc). Default 100 pmc.
) -> Dict[str, Any]:
    """
    Infer mean residence time (tau) from a single tracer observation.
    
    Parameters
    ----------
    obs_value : float
        Observed concentration (e.g., TU or pmc).
    obs_sigma : float
        Measurement uncertainty (1 sigma).
    sample_date_year : float
        Date of sampling (decimal year, e.g. 2023.4).
    nuclide : Nuclide
        Nuclide definition (constants).
    input_hist : Optional[InputHistory]
        Input function history (For 3H). If None, uses default/global.
        For C14, this is usually ignored in favor of 'initial_activity'.
    model : str
        "PFM" (Piston Flow) or "EM" (Exponential Mixing).
    tau_range : Tuple[float, float]
        Min and max age to search (years).
    initial_activity: float
        For C14: The initial activity A0 (q * A_atm). 
        Default 100.0 (no dead carbon). Adjust to ~85.0 or lower for carbonates.
    
    Returns
    -------
    Dict with keys:
        - tau_map_years: Most probable age
        - tau_ci_low_years: 2.5 percentile
        - tau_ci_high_years: 97.5 percentile
        - multimodal: bool (true if multiple probability peaks)
        - posterior_grid: dict of {tau, p}
    """
    
    # Pre-calculate lambda in 1/years
    lambda_y = np.log(2) / nuclide.half_life_years

    # Special handling for Radiocarbon (C14)
    # C14 usually uses a constant input A0, not a time-series input_hist
    if nuclide.symbol == "14C":
        # Create a constant input history at A0
        # A0 = initial_activity (user provided, accounting for q)
        # Note: This overrides input_hist for C14
        years = np.array([sample_date_year - tau_range[1] - 100, sample_date_year + 10])
        values = np.array([initial_activity, initial_activity])
        input_hist = InputHistory(years, values)
    elif input_hist is None:
        input_hist = get_input_history("global")
        
    # Grid search
    taus = np.linspace(tau_range[0], tau_range[1], search_steps)
    preds = []
    
    for t in taus:

        c_pred = convolve_input(
            sample_date_year,
            t,
            input_hist.years,
            input_hist.values,
            lambda_y,
            model_type=model,
            param2=param2
        )
        preds.append(c_pred)
        
    preds = np.array(preds)
    
    # Compute Likelihood
    # L(tau) = P(obs | tau) ~ N(pred(tau), obs_sigma)
    # Handle censoring if provided (e.g. < 0.1 TU)
    
    log_like = np.zeros_like(taus)
    sigma_eff = max(obs_sigma, 0.1) # avoid div by zero
    
    if censoring_limit is not None and obs_value <= censoring_limit:
        # Censored: L(tau) = P(pred(tau) <= limit) 
        # Actually, obs is "less than limit".
        # We can approximate with a step function or error function.
        # Let's say Prob(measured < limit | true = pred)
        # = CDF_normal((limit - pred) / sigma_inst)
        # But usually we just model measurement as N(pred, sigma).
        # If censored, we compute prob that measurement would fall below limit given pred.
        from scipy.special import erf
        # prob = 0.5 * (1 + erf((limit - pred) / (sigma * sqrt(2))))
        for i, pred in enumerate(preds):
            z = (censoring_limit - pred) / (sigma_eff * np.sqrt(2))
            prob = 0.5 * (1.0 + erf(z))
            log_like[i] = np.log(max(prob, 1e-9))
    else:
        # Standard Gaussian likelihood
        resid = obs_value - preds
        log_like = -0.5 * (resid / sigma_eff)**2
        
    # Prior? Uniform over search range usually sufficient for this specific tool.
    # Posterior
    # Normalize
    max_ll = np.max(log_like)
    post = np.exp(log_like - max_ll)
    norm = np.sum(post)
    if norm == 0:
        post = np.ones_like(post) / len(post)
    else:
        post = post / norm
        
    # Analyze Posterior
    # MAP
    idx_map = np.argmax(post)
    tau_map = float(taus[idx_map])
    
    # Credible Interval
    cdf = np.cumsum(post)
    p025 = np.searchsorted(cdf, 0.025)
    p975 = np.searchsorted(cdf, 0.975)
    tau_lo = float(taus[p025])
    tau_hi = float(taus[p975])
    
    # Multimodality check
    # Find local maxima
    peaks = []
    for i in range(1, len(post)-1):
        if post[i] > post[i-1] and post[i] > post[i+1]:
            if post[i] > 0.05 * np.max(post): # threshold
                peaks.append(i)
                
    is_multimodal = len(peaks) > 1
    
    if is_multimodal:
        logger.warning(f"Ambiguous age inference for {nuclide.symbol}={obs_value}: found {len(peaks)} peaks.")
    else:
        logger.debug(f"Inferred age for {nuclide.symbol}={obs_value}: {tau_map:.1f} years (95% CI: {tau_lo:.1f}-{tau_hi:.1f})")

    return {
        "tau_map_years": tau_map,
        "tau_ci_low_years": tau_lo,
        "tau_ci_high_years": tau_hi,
        "multimodal": is_multimodal,
        "posterior_grid": {"tau": taus.tolist(), "p": post.tolist()},
        "predicted_conc_at_map": float(preds[idx_map]),
        "model_used": model
    }
