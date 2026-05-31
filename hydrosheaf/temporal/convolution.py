"""Temporal Convolution Utilities."""

import math
from typing import List, Optional, Sequence

import numpy as np

def convolve_reactive_ttd(
    c_in: float,
    ttd_weights: Sequence[float],
    ttd_lags: Sequence[float],
    k: float,
    *,
    input_times: Optional[Sequence[float]] = None,
    input_values: Optional[Sequence[float]] = None,
    t_sample: Optional[float] = None,
) -> float:
    """
    Convolve input concentration with a reactive TTD.
    
    C_out = Sum( w_i * C_in(t_sample - tau_i) * exp(-k * tau_i) )
    
    Args:
        c_in: Input concentration used when no input time series is supplied
        ttd_weights: Discrete weights of the TTD (sum to 1)
        ttd_lags: Corresponding lag times (days)
        k: First-order decay constant (1/day)
        input_times: Optional input concentration timestamps (days)
        input_values: Optional input concentrations at input_times
        t_sample: Output sample time (days), required with input_times/input_values
        
    Returns:
        Reacted output concentration
    """
    weights = np.asarray(ttd_weights, dtype=float)
    lags = np.asarray(ttd_lags, dtype=float)
    if weights.shape != lags.shape:
        raise ValueError("ttd_weights and ttd_lags must have the same length")

    use_series = input_times is not None or input_values is not None or t_sample is not None
    if use_series:
        if input_times is None or input_values is None or t_sample is None:
            raise ValueError("input_times, input_values, and t_sample must be supplied together")
        times = np.asarray(input_times, dtype=float)
        values = np.asarray(input_values, dtype=float)
        if times.ndim != 1 or values.ndim != 1 or times.size != values.size:
            raise ValueError("input_times and input_values must be one-dimensional arrays of equal length")
        if times.size == 0:
            raise ValueError("input time series is empty")
        order = np.argsort(times)
        times = times[order]
        values = values[order]
        recharge_times = float(t_sample) - lags
        c_inputs = np.interp(
            recharge_times,
            times,
            values,
            left=float(values[0]),
            right=float(values[-1]),
        )
    else:
        c_inputs = np.full_like(lags, float(c_in), dtype=float)

    decay = np.exp(-float(k) * lags)
    return float(np.sum(weights * c_inputs * decay))

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
