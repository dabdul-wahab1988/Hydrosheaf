"""
Lumped Parameter Models (LPMs) for groundwater dating.
"""
from __future__ import annotations

import numpy as np

def piston_flow_model(tau_m: float, t: np.ndarray) -> np.ndarray:
    """
    Green's function for Piston Flow Model (PFM).
    g(t) = delta(t - tau_m)
    
    Since we do discrete convolution, this is handled logically rather than functionally usually.
    But for consistency, returns PDF values if t matches tau_m (discrete approx).
    """
    # Not used directly in convolution usually, but for completeness.
    return np.zeros_like(t) # Delta function handling is special

def exponential_model(tau_m: float, t: np.ndarray) -> np.ndarray:
    """
    Green's function for Exponential Model (EM).
    g(t) = (1/tau_m) * exp(-t/tau_m)
    """
    if tau_m <= 1e-6:
        return np.zeros_like(t) # Should approach delta at 0
    return (1.0 / tau_m) * np.exp(-t / tau_m)

def dispersion_model(tau_m: float, pd: float, t: np.ndarray) -> np.ndarray:
    """
    Dispersion Model (DM) approx (1D advection-dispersion).
    pd = Dispersion Parameter (1/Pe)
    """
    if tau_m <= 1e-6 or pd <= 1e-6:
         return np.zeros_like(t)
    
    # Simple approx for flux concentration
    pi = np.pi
    term1 = 1.0 / (tau_m * np.sqrt(4.0 * pi * pd * (t / tau_m)))
    term2 = np.exp(-((1.0 - t / tau_m) ** 2) / (4.0 * pd * (t / tau_m)))
    return term1 * term2


def convolve_input(
    t_sample_year: float,
    tau_mean_years: float,
    input_years: np.ndarray,
    input_values: np.ndarray,
    decay_lambda_inv_year: float,
    model_type: str = "EM",
    param2: float = 0.0
) -> float:
    """
    Compute output concentration at t_sample_year given an input history and LPM.
    
    C_out(t) = integral_0_inf C_in(t - tau) * g(tau) * exp(-lambda * tau) dtau
    """
    
    if model_type == "PFM":
        # Piston Flow: C_out(t) = C_in(t - tau) * exp(-lambda * tau)
        t_recharge = t_sample_year - tau_mean_years
        
        # Warning for physically questionable PFM results
        if tau_mean_years < 0:
            import warnings
            warnings.warn(f"Negative residence time ({tau_mean_years:.2f} years) encountered in PFM. This is physically impossible.")
        if t_recharge < input_years[0]:
             import warnings
             warnings.warn(f"Recharge date ({t_recharge:.1f}) for PFM exceeds input history start ({input_years[0]:.1f}). Using pre-bomb baseline (0.5 TU).")
             
        # Use a low baseline for pre-bomb water instead of the first entry in input_history
        c_in = np.interp(t_recharge, input_years, input_values, left=0.5)
        return float(c_in * np.exp(-decay_lambda_inv_year * tau_mean_years))
    
    # For distributed models (EM, DM, etc.), we need to integrate.
    # We'll discretize the transit time tau.
    
    # 1. Define integration grid for tau
    # Go out to 5 * tau_mean or until weight is negligible
    max_tau = max(100.0, tau_mean_years * 5.0) 
    # Resolution
    d_tau = 0.25 # quarter-year steps
    taus = np.arange(0.0, max_tau, d_tau)
    
    # 2. Compute g(tau)
    if model_type == "EM":
        # g(tau) = 1/Tm * exp(-tau/Tm)
        if tau_mean_years < 1e-3:
             # Treat as piston flow at 0
             t_recharge = t_sample_year
             c_in = np.interp(t_recharge, input_years, input_values)
             return float(c_in)
             
        g = (1.0 / tau_mean_years) * np.exp(-taus / tau_mean_years)
        
    elif model_type == "DM":
        # Dispersion model
        # param2 is Pd
        # Avoid division by zero
        ts = taus.copy()
        ts[ts==0] = 1e-9
        term1 = 1.0 / (tau_mean_years * np.sqrt(4.0 * np.pi * param2 * (ts / tau_mean_years)))
        term2 = np.exp(-((1.0 - ts / tau_mean_years) ** 2) / (4.0 * param2 * (ts / tau_mean_years)))
        g = term1 * term2
        g[taus==0] = 0 # boundary condition
        
    else:
        raise ValueError(f"Unknown model type: {model_type}")
        
    # Normalize g to ensure mass conservation (sum * d_tau approx 1)
    area = np.sum(g) * d_tau
    if area > 0:
        g = g / area
        
    # 3. Get C_in(t_sample - tau)
    t_recharge_grid = t_sample_year - taus
    c_in_grid = np.interp(t_recharge_grid, input_years, input_values)
    
    # 4. Radioactive decay
    decay = np.exp(-decay_lambda_inv_year * taus)
    
    # 5. Integrate
    # integrand = C_in * g * decay
    integrand = c_in_grid * g * decay
    c_out = np.sum(integrand) * d_tau
    
    return float(c_out)
