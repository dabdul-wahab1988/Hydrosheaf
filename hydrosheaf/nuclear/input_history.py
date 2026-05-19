"""
Input history handling for transient tracers (e.g., Tritium).
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple, Union
import numpy as np
import pandas as pd
from pathlib import Path
from ..log import get_logger

logger = get_logger(__name__)

class InputHistory:
    """
    Represents the time-variant input concentration of a tracer (C_in(t)).
    """
    def __init__(self, dates: np.ndarray, values: np.ndarray, sigma: Optional[np.ndarray] = None):
        """
        Parameters
        ----------
        dates : np.ndarray
            Dates as fractional years (e.g., 1963.5) or timestamps.
            For internal calculation, we store as fractional years or days?
            Let's stick to fractional years for input but convert to absolute days for convolution if needed.
            Actually, the rest of Hydrosheaf uses timestamps/days.
            Let's store as days relative to an epoch or just fractional years if that's standard for nuclear.
            Hydrosheaf uses 'days since epoch' internally often.
            Let's store as fractional years (standard for 3H input curves) and provide conversion.
        values : np.ndarray
            Concentration values (e.g., TU).
        sigma : Optional[np.ndarray]
            Uncertainty (1 std dev).
        """
        self.years = np.array(dates, dtype=float)
        self.values = np.array(values, dtype=float)
        if sigma is not None:
            self.sigma = np.array(sigma, dtype=float)
        else:
            self.sigma = np.zeros_like(self.values)
        
        # Sort by date
        idx = np.argsort(self.years)
        self.years = self.years[idx]
        self.values = self.values[idx]
        self.sigma = self.sigma[idx]

    def interpolate(self, target_years: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Interpolate input history to target dates.
        Returns (values, sigmas).
        """
        vals = np.interp(target_years, self.years, self.values)
        sigs = np.interp(target_years, self.years, self.sigma)
        return vals, sigs


def load_input_series_csv(path: Union[str, Path]) -> InputHistory:
    """
    Load input history from a CSV file.
    Expected columns: 'year', 'value', 'sigma' (optional).
    """
    df = pd.read_csv(path)
    # Simple normalization
    cols = [c.lower() for c in df.columns]
    df.columns = cols
    
    if 'year' not in df.columns or 'value' not in df.columns:
        raise ValueError("CSV must contain 'year' and 'value' columns.")
    
    sigma = df['sigma'].values if 'sigma' in df.columns else None
    return InputHistory(df['year'].values, df['value'].values, sigma)


def build_default_tritium_input(region: str = "global") -> InputHistory:
    """
    Generate a default tritium input function.
    
    Parameters
    ----------
    region : str
        "global" / "northern_continental": Standard Ottawa-style peak (~3000 TU).
        "tropical" / "coastal": Dampened peak (~500 TU) for equatorial/coastal regions (e.g., Ghana).
    """
    # Base Curve: Stylized Ottawa/Vienna-like curve
    years = np.arange(1950.0, 2025.0, 0.5)
    values = []
    
    # Scale factor based on region
    # Northern Continental is the reference (1.0)
    # Tropics/Coastal is much lower due to ocean dilution (approx 0.15 - 0.2)
    scale = 1.0
    r_lower = region.lower()
    if "tropic" in r_lower or "ghana" in r_lower or "coast" in r_lower:
        scale = 0.2  # Approx scaling for Accra/tropical stations vs Ottawa
        logger.info(f"Using Tropical/Coastal Tritium scaling (factor={scale}) for region '{region}'")
    else:
        logger.info(f"Using Standard Northern Continental Tritium scaling (factor={scale}) for region '{region}'")
    
    for y in years:
        if y < 1952:
            v = 5.0  # Pre-bomb background
        elif 1952 <= y < 1963:
            # Rise
            v = 5.0 + 100.0 * np.exp((y - 1952) / 2.0)
        elif 1963 <= y < 1964:
            # Peak
            v = 3000.0  # approximate peak TU
        else:
            # Decay + washout
            dt = y - 1964
            v = 3000.0 * np.exp(-dt / 3.0) + 10.0 # decaying to modern background
        
        # Apply regional scaling
        # Note: Background (5-10 TU) is less affected by scaling than the peak, 
        # but for a crude approx, scaling everything is safer than overestimating the peak.
        # A more refined model would scale the 'excess' only.
        v_excess = max(0, v - 3.0) # Assume cosmogenic base ~3 TU
        v_final = 3.0 + v_excess * scale
        
        values.append(v_final)
        
    return InputHistory(years, np.array(values), sigma=np.array(values)*0.1)

# Global cache for input functions to avoid reloading
_INPUT_CACHE: Dict[str, InputHistory] = {}

def get_input_history(region_or_path: str) -> InputHistory:
    if region_or_path in _INPUT_CACHE:
        return _INPUT_CACHE[region_or_path]
    
    path = Path(region_or_path)
    if path.exists():
        hist = load_input_series_csv(path)
    else:
        # Fallback to generator
        hist = build_default_tritium_input(region_or_path)
    
    _INPUT_CACHE[region_or_path] = hist
    return hist
