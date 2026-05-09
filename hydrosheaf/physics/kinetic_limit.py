"""Kinetic feasibility filtering for inverse geochemical models."""

from typing import Dict, List, Mapping, Optional, Tuple
import math

from ..config import Config

# Standard Mineral Kinetic Data (Mol/m2/s at 25C)
# Values adapted from Palandri and Kharaka (2004) or similar.
# log10(k)
MINERAL_KINETICS = {
    "calcite": -6.0,
    "dolomite": -7.5,
    "gypsum": -4.0,
    "anhydrite": -4.0,
    "halite": -2.0,
    "albite": -12.0,
    "anorthite": -11.0,
    "k_feldspar": -12.5,
    "quartz": -13.0,
    "pyrite_net": -8.5,
}

def calculate_kinetic_limit(
    mineral: str,
    residence_time_days: float,
    config: Config,
) -> float:
    """
    Calculate the maximum theoretical dissolution/precipitation extent (mmol/L)
    based on residence time and mineral kinetics.
    
    Beta_max = k * A * tau
    k: rate constant (mol/m2/s)
    A: reactive surface area (m2/L)
    tau: residence time (s)
    """
    if residence_time_days <= 0:
        return 0.0
    
    # Get rate constant from config or default
    k_log = config.mineral_rate_constants.get(mineral.lower())
    if k_log is None:
        k_log = MINERAL_KINETICS.get(mineral.lower(), -10.0)
    
    k = 10**k_log
    area = config.kinetic_surface_area_m2_L
    tau_sec = residence_time_days * 86400.0
    
    # Limit in mol/L
    limit_mol = k * area * tau_sec
    
    # Convert to mmol/L
    return limit_mol * 1000.0

def apply_kinetic_penalties(
    extents: List[float],
    labels: List[str],
    residence_time_days: Optional[float],
    config: Config,
) -> float:
    """
    Returns a penalty score if reaction extents exceed kinetic limits.
    """
    if not config.kinetic_filter_enabled or residence_time_days is None:
        return 0.0
    
    penalty = 0.0
    for label, extent in zip(labels, extents):
        limit = calculate_kinetic_limit(label, residence_time_days, config)
        
        # If extent is significantly higher than physical limit
        if abs(extent) > limit * 1.5:
            # Quadratic penalty for violation
            violation = abs(extent) - limit
            penalty += 10.0 * (violation ** 2)
            
    return penalty
