"""
Analytical solutions for 1D solute transport.
"""

import numpy as np
from scipy import signal

def ade_1d_green(t: np.ndarray, x: float, v: float, D: float, k: float = 0.0) -> np.ndarray:
    """
    1D ADE Green's function for flux concentration (Jury & Roth, 1990).
    Calculates the response at distance x and time t to a Dirac delta injection at x=0, t=0.
    
    Args:
        t: Time array.
        x: Distance from source.
        v: Pore water velocity.
        D: Dispersion coefficient.
        k: First-order decay rate (1/t).
        
    Returns:
        Green's function g(t).
    """
    # Avoid division by zero
    t_safe = np.maximum(t, 1e-9)
    
    # Term 1: Advection-Dispersion
    # g(t) = (x / (2 * sqrt(pi * D * t^3))) * exp( - (x - v*t)^2 / (4*D*t) )
    term1 = x / (2 * np.sqrt(np.pi * D * t_safe**3))
    term2 = np.exp( -((x - v*t_safe)**2) / (4*D*t_safe) )
    
    # Term 3: Decay
    term3 = np.exp(-k * t_safe)
    
    g = term1 * term2 * term3
    
    # Zero out negative or zero times (causality)
    g[t <= 0] = 0.0
    return g

def simulate_1d_ade(
    t_grid: np.ndarray, 
    c_in: np.ndarray, 
    x: float, 
    v: float, 
    D: float, 
    k: float = 0.0,
    dilution: float = 1.0,
    baseflow: float = 0.0
) -> np.ndarray:
    """
    Simulate 1D transport using convolution of input signal with ADE Green's function.
    
    Args:
        t_grid: Regular time grid (must be uniform dt).
        c_in: Input concentration series on t_grid.
        x: Distance.
        v: Velocity.
        D: Dispersion.
        k: Decay.
        dilution: Scaling factor for output (C_out = C_conv * dilution).
        baseflow: Additive background concentration.
        
    Returns:
        Output concentration series on t_grid.
    """
    g = ade_1d_green(t_grid, x, v, D, k)
    
    # Convolve (dt is assumed 1.0 if not scaled, here we assume t_grid is index or dt is handled in parameters)
    # Ideally, we should multiply by dt. Assuming t_grid has spacing 1.
    dt = t_grid[1] - t_grid[0] if len(t_grid) > 1 else 1.0
    
    # scipy.signal.fftconvolve is efficient
    # We multiply by dt because convolution integral is Sum(a*b)*dt
    c_out_full = signal.fftconvolve(c_in, g, mode='full') * dt
    
    # Trim to original length
    c_out = c_out_full[:len(t_grid)]
    
    # Apply mixing
    return (c_out * dilution) + baseflow
