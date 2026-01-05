"""Gibbs diagram classification for hydrogeochemical process identification."""

from typing import Dict, Mapping, Optional, Tuple


def compute_gibbs_ratios(
    sample: Mapping[str, object],
    na_key: str = "Na",
    ca_key: str = "Ca",
    cl_key: str = "Cl",
    hco3_key: str = "HCO3",
    tds_key: str = "TDS",
) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """
    Compute Gibbs diagram ratios.
    
    Returns:
        Tuple of (TDS, Na/(Na+Ca), Cl/(Cl+HCO3))
        None values if required ions missing
    """
    try:
        na = float(sample.get(na_key) or 0)
        ca = float(sample.get(ca_key) or 0)
        cl = float(sample.get(cl_key) or 0)
        hco3 = float(sample.get(hco3_key) or 0)
        tds = sample.get(tds_key)
        if tds is not None:
            tds = float(tds)
    except (TypeError, ValueError):
        return None, None, None
    
    # Cation ratio: Na/(Na+Ca)
    cation_sum = na + ca
    ratio_cation = na / cation_sum if cation_sum > 0 else None
    
    # Anion ratio: Cl/(Cl+HCO3)
    anion_sum = cl + hco3
    ratio_anion = cl / anion_sum if anion_sum > 0 else None
    
    return tds, ratio_cation, ratio_anion


def classify_gibbs_dominance(
    sample: Mapping[str, object],
    tds_precipitation: float = 100.0,
    tds_evaporation: float = 1000.0,
    ratio_threshold: float = 0.5,
    **kwargs,
) -> str:
    """
    Classify dominant hydrogeochemical process using Gibbs diagram.
    
    Args:
        sample: Sample dictionary with ion concentrations
        tds_precipitation: Upper TDS bound for precipitation dominance (mg/L)
        tds_evaporation: Lower TDS bound for evaporation dominance (mg/L)
        ratio_threshold: Ratio threshold for classification
        
    Returns:
        "precipitation", "rock", or "evaporation"
    """
    tds, ratio_cation, ratio_anion = compute_gibbs_ratios(sample, **kwargs)
    
    if tds is None:
        # Fall back to ratio-based classification
        if ratio_cation is None or ratio_anion is None:
            return "unknown"
        avg_ratio = (ratio_cation + ratio_anion) / 2
        if avg_ratio < 0.3:
            return "precipitation"
        elif avg_ratio > 0.7:
            return "evaporation"
        else:
            return "rock"
    
    # TDS-based classification with ratio confirmation
    if tds < tds_precipitation:
        return "precipitation"
    elif tds > tds_evaporation:
        if ratio_cation is not None and ratio_cation > ratio_threshold:
            return "evaporation"
        if ratio_anion is not None and ratio_anion > ratio_threshold:
            return "evaporation"
        return "evaporation"
    else:
        return "rock"


def gibbs_transport_weights(
    sample: Mapping[str, object],
    tds_precipitation: float = 100.0,
    tds_evaporation: float = 1000.0,
    **kwargs,
) -> Dict[str, float]:
    """
    Compute transport model weights based on Gibbs classification.
    
    Returns probability-like weights for transport model selection.
    These weights are used when isotope data is unavailable.
    
    Returns:
        Dict with "evap" and "mix" weights (sum to 1.0)
    """
    dominance = classify_gibbs_dominance(
        sample,
        tds_precipitation=tds_precipitation,
        tds_evaporation=tds_evaporation,
        **kwargs,
    )
    
    if dominance == "evaporation":
        # High evaporation signal -> favor evaporation model
        return {"evap": 0.85, "mix": 0.15}
    elif dominance == "precipitation":
        # Precipitation dominance -> could be recent recharge, mixing more likely
        return {"evap": 0.3, "mix": 0.7}
    else:
        # Rock dominance -> neutral, slight favor to evaporation
        return {"evap": 0.6, "mix": 0.4}


def gibbs_evaporation_penalty(
    sample: Mapping[str, object],
    tds_precipitation: float = 100.0,
    tds_evaporation: float = 1000.0,
    **kwargs,
) -> float:
    """
    Compute penalty for evaporation model based on Gibbs classification.
    
    Lower penalty means evaporation is more consistent with Gibbs classification.
    This is used to supplement isotope penalties.
    
    Returns:
        Penalty value (0 = consistent with evaporation, higher = inconsistent)
    """
    dominance = classify_gibbs_dominance(
        sample,
        tds_precipitation=tds_precipitation,
        tds_evaporation=tds_evaporation,
        **kwargs,
    )
    
    if dominance == "evaporation":
        return 0.0
    elif dominance == "precipitation":
        return 1.0
    else:  # rock dominance
        return 0.3


def compute_gibbs_metrics(
    sample: Mapping[str, object],
    **kwargs,
) -> Dict[str, object]:
    """
    Compute all Gibbs-related metrics for a sample.
    
    Returns:
        Dict with TDS, ratios, classification, and weights
    """
    tds, ratio_cation, ratio_anion = compute_gibbs_ratios(sample, **kwargs)
    dominance = classify_gibbs_dominance(sample, **kwargs)
    weights = gibbs_transport_weights(sample, **kwargs)
    penalty = gibbs_evaporation_penalty(sample, **kwargs)
    
    return {
        "tds": tds,
        "ratio_na_naca": ratio_cation,
        "ratio_cl_clhco3": ratio_anion,
        "gibbs_dominance": dominance,
        "gibbs_evap_weight": weights.get("evap"),
        "gibbs_mix_weight": weights.get("mix"),
        "gibbs_evap_penalty": penalty,
    }
