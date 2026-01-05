"""Ion exchange reaction calculations and Chloro-Alkaline Indices."""

from typing import Dict, Mapping, Optional, Tuple


def compute_cai_indices(
    sample: Mapping[str, object],
    na_key: str = "Na",
    k_key: str = "K",
    cl_key: str = "Cl",
    so4_key: str = "SO4",
    hco3_key: str = "HCO3",
    no3_key: str = "NO3",
) -> Tuple[Optional[float], Optional[float]]:
    """
    Compute Chloro-Alkaline Indices (CAI-1 and CAI-2).
    
    CAI-1 = (Cl - (Na + K)) / Cl
    CAI-2 = (Cl - (Na + K)) / (SO4 + HCO3 + NO3)
    
    Interpretation:
    - CAI > 0: Reverse ion exchange (salinization) - Na released, Ca absorbed
    - CAI < 0: Direct ion exchange (freshening) - Ca released, Na absorbed
    - CAI ≈ 0: No significant ion exchange
    
    Args:
        sample: Sample dictionary with ion concentrations (meq/L preferred)
        
    Returns:
        Tuple of (CAI-1, CAI-2), None values if required ions missing
    """
    try:
        na = float(sample.get(na_key) or 0)
        k = float(sample.get(k_key) or 0)
        cl = float(sample.get(cl_key) or 0)
        so4 = float(sample.get(so4_key) or 0)
        hco3 = float(sample.get(hco3_key) or 0)
        no3 = float(sample.get(no3_key) or 0)
    except (TypeError, ValueError):
        return None, None
    
    # CAI-1: Cl-based
    if cl > 0:
        cai1 = (cl - (na + k)) / cl
    else:
        cai1 = None
    
    # CAI-2: Anion-based
    anion_sum = so4 + hco3 + no3
    if anion_sum > 0:
        cai2 = (cl - (na + k)) / anion_sum
    else:
        cai2 = None
    
    return cai1, cai2


def classify_exchange_direction(
    sample: Mapping[str, object],
    threshold: float = 0.1,
    **kwargs,
) -> str:
    """
    Classify ion exchange direction based on CAI indices.
    
    Args:
        sample: Sample dictionary with ion concentrations
        threshold: CAI magnitude threshold for classification
        
    Returns:
        "freshening" (CAI < -threshold), "salinization" (CAI > threshold), or "neutral"
    """
    cai1, cai2 = compute_cai_indices(sample, **kwargs)
    
    if cai1 is None and cai2 is None:
        return "unknown"
    
    # Use average of available indices
    indices = [i for i in [cai1, cai2] if i is not None]
    avg_cai = sum(indices) / len(indices)
    
    if avg_cai < -threshold:
        return "freshening"
    elif avg_cai > threshold:
        return "salinization"
    else:
        return "neutral"


def exchange_reaction_vectors(ion_order: list) -> Dict[str, list]:
    """
    Get ion exchange reaction stoichiometry vectors.
    
    Exchange reactions follow 1:2 ratio (Ca/Mg : Na):
    - CaNa_exch: Ca²⁺ + 2Na-X → 2Na⁺ + Ca-X₂  
    - MgNa_exch: Mg²⁺ + 2Na-X → 2Na⁺ + Mg-X₂
    
    Returns:
        Dict mapping reaction name to coefficient vector
    """
    def _make_vector(coeffs: dict) -> list:
        return [float(coeffs.get(ion, 0.0)) for ion in ion_order]
    
    return {
        # Direct exchange (freshening): Ca/Mg enters solution, Na leaves
        "CaNa_exch": _make_vector({"Ca": 1, "Na": -2}),
        "MgNa_exch": _make_vector({"Mg": 1, "Na": -2}),
    }


def suggest_exchange_bounds(
    sample_u: Mapping[str, object],
    sample_v: Mapping[str, object],
    threshold: float = 0.1,
) -> Dict[str, Tuple[float, float]]:
    """
    Suggest reaction bounds based on CAI evolution along flow path.
    
    If downstream sample shows freshening trend, allow positive CaNa/MgNa exchange.
    If downstream sample shows salinization, allow negative exchange (reverse).
    
    Returns:
        Dict mapping reaction name to (lb, ub) bounds
    """
    dir_u = classify_exchange_direction(sample_u, threshold)
    dir_v = classify_exchange_direction(sample_v, threshold)
    
    inf = float("inf")
    bounds = {}
    
    # Default: allow both directions (signed reaction)
    if dir_v == "freshening":
        # Freshening at v: Ca/Mg released, Na absorbed -> positive extent
        bounds["CaNa_exch"] = (0.0, inf)
        bounds["MgNa_exch"] = (0.0, inf)
    elif dir_v == "salinization":
        # Salinization at v: Na released, Ca/Mg absorbed -> negative extent
        bounds["CaNa_exch"] = (-inf, 0.0)
        bounds["MgNa_exch"] = (-inf, 0.0)
    else:
        # Neutral or unknown: allow both directions
        bounds["CaNa_exch"] = (-inf, inf)
        bounds["MgNa_exch"] = (-inf, inf)
    
    return bounds
