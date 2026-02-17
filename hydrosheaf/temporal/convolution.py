"""Temporal Convolution Utilities."""

import math
from typing import List, Tuple

def convolve_reactive_ttd(
    c_in: float,
    ttd_weights: List[float],
    ttd_lags: List[float],
    k: float
) -> float:
    """
    Convolve input concentration with a reactive TTD.
    
    C_out = Sum( w_i * C_in * exp(-k * tau_i) )
    
    Args:
        c_in: Input concentration (assumed constant or time-averaged for this step)
        ttd_weights: Discrete weights of the TTD (sum to 1)
        ttd_lags: Corresponding lag times (days)
        k: First-order decay constant (1/day)
        
    Returns:
        Reacted output concentration
    """
    c_out = 0.0
    for w, tau in zip(ttd_weights, ttd_lags):
        # First order decay along flowpath
        decay = math.exp(-k * tau)
        c_out += w * c_in * decay
    return c_out

def compute_effective_reaction_factor(
    ttd_weights: List[float],
    ttd_lags: List[float],
    k: float
) -> float:
    """
    Compute the effective attenuation factor R = C_out / C_in.
    
    This replaces the simple exp(-k * mean_age) used in piston flow models.
    """
    return convolve_reactive_ttd(1.0, ttd_weights, ttd_lags, k)
