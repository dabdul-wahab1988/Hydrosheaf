"""CoDA SBP-ilr coordinate system for 7-ion composition.

Implements the 7-ion CoDA composition (Ca, Mg, Na, K, HCO3, Cl, SO4)
defined in Hydrosheaf v2 Nitrate Source Discrimination Layer specs.
NO3 is intentionally excluded from the composition.
"""

import math
import numpy as np
from typing import Dict, List, Optional, Sequence, Tuple

# 7-ion order for CoDA
IONS_7_CODA = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4"]

# Sequential Binary Partition (SBP) Tree Definition
# Tuple structure: (Name, Group 1 Ions, Group 2 Ions)
# R_j = Group 1 (numerator), S_j = Group 2 (denominator)
# Groups defined by indices in IONS_7_CODA
# 0:Ca, 1:Mg, 2:Na, 3:K, 4:HCO3, 5:Cl, 6:SO4

SBP_DEFINITIONS = [
    # 1. Cations vs Anions
    # R1 = {Ca,Mg,Na,K} (0,1,2,3), S1 = {HCO3,Cl,SO4} (4,5,6)
    ("ilr1_cat_an", [0, 1, 2, 3], [4, 5, 6]),
    
    # 2. Alkaline earth vs Alkali (inside cations)
    # R2 = {Ca,Mg} (0,1), S2 = {Na,K} (2,3)
    ("ilr2_alkearth_alkali", [0, 1], [2, 3]),
    
    # 3. Ca vs Mg
    # R3 = {Ca} (0), S3 = {Mg} (1)
    ("ilr3_ca_mg", [0], [1]),
    
    # 4. Na vs K
    # R4 = {Na} (2), S4 = {K} (3)
    ("ilr4_na_k", [2], [3]),
    
    # 5. Bicarbonate vs (Cl+SO4) (inside anions)
    # R5 = {HCO3} (4), S5 = {Cl,SO4} (5,6)
    ("ilr5_hco3_salinity", [4], [5, 6]),
    
    # 6. Cl vs SO4
    # R6 = {Cl} (5), S6 = {SO4} (6)
    ("ilr6_cl_so4", [5], [6]),
]


def geometric_mean(values: Sequence[float]) -> float:
    """Compute geometric mean of a sequence of positive numbers."""
    if not values:
        return 0.0
    log_sum = sum(math.log(v) for v in values)
    return math.exp(log_sum / len(values))


def multiplicative_replace(
    values: Sequence[float], epsilon: float = 1e-12
) -> List[float]:
    """Replace zero or negative values with small epsilon."""
    return [max(v, epsilon) for v in values]


def compute_ilr_coordinate(
    values: Sequence[float], r_indices: List[int], s_indices: List[int]
) -> float:
    """Compute single ilr coordinate for a balance."""
    r_vals = [values[i] for i in r_indices]
    s_vals = [values[i] for i in s_indices]
    
    r = len(r_vals)
    s = len(s_vals)
    
    if r == 0 or s == 0:
        return 0.0
        
    g_r = geometric_mean(r_vals)
    g_s = geometric_mean(s_vals)
    
    if g_r <= 0 or g_s <= 0:
        # Should be handled by replacement, but safety check
        return 0.0
        
    coeff = math.sqrt((r * s) / (r + s))
    return coeff * math.log(g_r / g_s)


def ilr_from_sbp(
    sample: Dict[str, float], epsilon: float = 1e-12
) -> Tuple[Optional[List[float]], bool]:
    """Compute 6 ilr coordinates from sample dictionary.
    
    Args:
        sample: Dictionary containing ion concentrations
        epsilon: Small value for zero replacement
        
    Returns:
        (ilr_vector, is_valid)
        ilr_vector: List of 6 floats or None if insufficient data
        is_valid: boolean indicating if calculation was successful
    """
    # 1. Extract 7 ions
    raw_values = []
    for ion in IONS_7_CODA:
        val = sample.get(ion)
        if val is None or math.isnan(val):
            return None, False
        raw_values.append(float(val))
        
    # 2. Replaced values
    values = multiplicative_replace(raw_values, epsilon)
    
    # 3. Compute coords
    ilr_vec = []
    for _, r_idx, s_idx in SBP_DEFINITIONS:
        coord = compute_ilr_coordinate(values, r_idx, s_idx)
        ilr_vec.append(coord)
        
    return ilr_vec, True


def robust_zscore(
    value: float, median: float, mad: float, epsilon: float = 1e-12
) -> float:
    """Compute robust z-score using Median and MAD.
    
    z = (x - median) / (1.4826 * MAD)
    """
    scale = 1.4826 * mad
    if scale < epsilon:
        return 0.0
    return (value - median) / scale
