"""Gibbs diagram classification for hydrogeochemical process identification."""

import math
from typing import Dict, Mapping, Optional, Tuple


def _logistic(x: float, x0: float, k: float) -> float:
    """Logistic function: 1 / (1 + exp(-k*(x-x0)))"""
    try:
        return 1.0 / (1.0 + math.exp(-k * (x - x0)))
    except OverflowError:
        return 1.0 if x > x0 else 0.0


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
    Updated with probabilistic logic to prevent high-TDS misclassification.
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

    # 1. Probability of Precipitation Dominance
    # Characterized by LOW TDS and HIGH Ratios
    p_precip_tds = 1.0 - _logistic(tds, tds_precipitation, 0.05)
    p_precip_ratio = _logistic(max(ratio_cation or 0, ratio_anion or 0), 0.7, 10.0)
    p_precip = p_precip_tds * p_precip_ratio

    # 2. Probability of Evaporation Dominance
    # Characterized by HIGH TDS and HIGH Ratios
    p_evap_tds = _logistic(tds, tds_evaporation, 0.005)
    p_evap_ratio = _logistic(max(ratio_cation or 0, ratio_anion or 0), 0.6, 10.0)
    p_evap = p_evap_tds * p_evap_ratio

    # 3. Probability of Rock Dominance
    # Characterized by Intermediate TDS or High TDS with LOW Ratios
    # (High TDS with low Na ratio = Evaporite dissolution, which is 'Rock')
    p_rock = max(0.0, 1.0 - p_precip - p_evap)

    # Hard classification based on max probability
    probs = {"precipitation": p_precip, "evaporation": p_evap, "rock": p_rock}
    return max(probs, key=lambda k: probs[k])


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
    Compute smooth penalty for evaporation model based on Gibbs probabilities.
    Penalty = 1.0 - P(Evaporation).
    """
    tds, ratio_cation, ratio_anion = compute_gibbs_ratios(sample, **kwargs)
    if tds is None or ratio_cation is None:
        return 0.3 # Default neutral penalty

    p_evap_tds = _logistic(tds, tds_evaporation, 0.005)
    p_evap_ratio = _logistic(max(ratio_cation or 0, ratio_anion or 0), 0.6, 10.0)
    p_evap = p_evap_tds * p_evap_ratio
    
    return 1.0 - p_evap


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
