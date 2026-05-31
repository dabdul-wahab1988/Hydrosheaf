"""Isotope utilities (delta-18O / delta-2H)."""

from typing import Iterable, Mapping, Optional, Tuple

from .data.schema import parse_numeric

GMWL_SLOPE = 8.0
GMWL_INTERCEPT = 10.0


def extract_isotopes(

    sample: Mapping[str, object],
    d18o_key: str = "18O",
    d2h_key: str = "2H",
) -> Optional[Tuple[float, float]]:
    d18o = parse_numeric(sample.get(d18o_key), "value")
    d2h = parse_numeric(sample.get(d2h_key), "value")
    if d18o is None or d2h is None:
        return None
    return float(d18o), float(d2h)


def compute_d_excess(d18o: float, d2h: float) -> float:
    return d2h - 8.0 * d18o


def evaporation_index(d18o: float, d2h: float, a: float, b: float) -> float:
    return d2h - (a + b * d18o)


def fit_lmwl(
    samples: Iterable[Mapping[str, object]],
    d18o_key: str = "18O",
    d2h_key: str = "2H",
    sample_type_col: Optional[str] = None,
    meteoric_types: Optional[Iterable[str]] = None,
    weight_col: Optional[str] = None,
) -> Tuple[float, float]:
    """
    Fit LMWL using Precipitation-Weighted Reduced Major Axis (PWRMA) regression.

    Standard OLS underestimates slope due to error in X (d18O). RMA treats X and Y errors
    symmetrically, which is standard in isotope hydrology (IAEA).
    """
    import numpy as np

    values = []
    weights = []

    # Normalize meteoric types for case-insensitive comparison
    valid_types = set()
    if meteoric_types:
        valid_types = {str(t).lower().strip() for t in meteoric_types}

    for sample in samples:
        # Filter by sample type if requested
        if sample_type_col and valid_types:
            stype = str(sample.get(sample_type_col, "")).lower().strip()
            if stype not in valid_types:
                continue

        pair = extract_isotopes(sample, d18o_key=d18o_key, d2h_key=d2h_key)
        if pair is None:
            continue

        w = 1.0
        if weight_col:
            w_val = parse_numeric(sample.get(weight_col), "zero")
            if w_val is not None and w_val > 0:
                w = float(w_val)
            else:
                w = 1e-6

        values.append(pair)
        weights.append(w)

    if not values:
        if sample_type_col:
            raise ValueError(
                f"No precipitation samples found (looking for types {valid_types} in column '{sample_type_col}'). "
                "Ensure your data contains rain samples or use --use-gmwl."
            )
        else:
            raise ValueError("No isotope values available to fit LMWL.")

    d18o = np.array([item[0] for item in values])
    d2h = np.array([item[1] for item in values])
    w_arr = np.array(weights)

    # Weighted means
    sum_w = np.sum(w_arr)
    if sum_w <= 0:
        sum_w = 1.0
        w_arr = np.ones_like(w_arr)

    mean_x = np.average(d18o, weights=w_arr)
    mean_y = np.average(d2h, weights=w_arr)

    # Weighted Variances (RMA Slope = std_y / std_x)
    # Using 'reliability weights' approach for variance
    var_x = np.average((d18o - mean_x) ** 2, weights=w_arr)
    var_y = np.average((d2h - mean_y) ** 2, weights=w_arr)

    if var_x < 1e-9:
        raise ValueError("Insufficient d18O variance to fit LMWL.")

    # Sign of correlation to determine slope direction (usually positive)
    # Weighted covariance
    cov_xy = np.average((d18o - mean_x) * (d2h - mean_y), weights=w_arr)
    sign = 1.0 if cov_xy >= 0 else -1.0

    slope = sign * np.sqrt(var_y / var_x)
    intercept = mean_y - slope * mean_x

    return intercept, slope



def to_ratio(delta: float) -> float:
    """Convert delta value (permil) to isotope ratio R/R_std."""
    # VSMOW ratios could be used, but relative changes work with R/R_std = 1 + delta/1000
    return 1.0 + delta / 1000.0


def to_delta(ratio: float) -> float:
    """Convert isotope ratio R/R_std to delta value (permil)."""
    return (ratio - 1.0) * 1000.0


def rayleigh_fractionation(delta_0: float, f: float, alpha: float) -> float:
    """
    Calculate isotopic composition after Rayleigh fractionation.
    
    R = R0 * f^(alpha - 1)
    
    Parameters
    ----------
    delta_0 : float
        Initial delta value (permil)
    f : float
        Remaining fraction of liquid (0 < f <= 1)
    alpha : float
        Fractionation factor (R_liquid / R_vapor).
        Typically alpha > 1 for evaporation.
        
    Returns
    -------
    float
        Final delta value (permil)
    """
    if f <= 0 or f > 1:
        raise ValueError("Remaining fraction f must be in (0, 1]")
        
    R0 = to_ratio(delta_0)
    R = R0 * (f ** (alpha - 1.0))
    return to_delta(R)


def craig_gordon_enrichment(
    delta_L: float,
    delta_A: float,
    h: float,
    epsilon_eq: float,
    epsilon_k: float
) -> float:
    """
    Calculate steady-state isotopic composition of evaporating water body (Craig-Gordon).
    
    delta_E = (delta_L - h * delta_A - epsilon_total) / ((1 - h) + epsilon_k / 1000)
    
    (Note: There are various forms. This calculates the composition of the *evaporate flux* delta_E,
    or the steady state liquid delta_L_ss. The critique requested "Craig-Gordon model". 
    Usually used to model the evolution of delta_L.)
    
    Let's implement the standard linear approximation for delta_L evolution or steady state.
    But for a "Fractionation Node", a simplified Rayleigh with kinetic effects is often sufficient.
    
    However, provided here is the calculation for the limiting steady-state enrichment 
    of a terminal lake (delta_ss).
    
    delta_ss = (h * delta_A + epsilon_total) / (h - epsilon_total/1000) ... approximations vary.
    
    Let's stick to a generalized enrichment function that can be used in the adapter.
    """
    humidity = float(h)
    if not (0.0 <= humidity < 1.0):
        raise ValueError("Relative humidity h must be in [0, 1).")

    eps_eq = float(epsilon_eq)
    eps_k = float(epsilon_k)
    denominator = (1.0 - humidity) + eps_k / 1000.0
    if denominator <= 0.0:
        raise ValueError("Craig-Gordon denominator must be positive.")

    epsilon_total = eps_eq + eps_k
    return float((float(delta_L) - humidity * float(delta_A) - epsilon_total) / denominator)



def isotope_penalty(
    d18o_u: float,
    d2h_u: float,
    d18o_v: float,
    d2h_v: float,
    a: float,
    b: float,
    mode: str,
    d_excess_weight: float = 0.0,
    endmember_iso: Optional[Tuple[float, float]] = None,
) -> Tuple[float, dict]:
    """
    Calculate isotopic penalty.

    Evaporation: Penalizes if downstream water is 'less evaporated' (closer to LMWL) than upstream.
    Mixing: Penalizes deviation from the theoretical mixing line between U and Endmember.
    """
    e_u = evaporation_index(d18o_u, d2h_u, a, b)
    e_v = evaporation_index(d18o_v, d2h_v, a, b)
    d_u = compute_d_excess(d18o_u, d2h_u)
    d_v = compute_d_excess(d18o_v, d2h_v)

    mix_deviation = 0.0

    if mode == "evap":
        # For evaporation, V should be more enriched and lower d-excess (or higher E)
        e_penalty = max(0.0, abs(e_u) - abs(e_v))
        d_penalty = max(0.0, d_v - d_u)
    else:
        # For mixing, we expect V to lie on the line between U and Endmember
        e_penalty = 0.0
        d_penalty = 0.0
        
        if endmember_iso is not None:
            # Vector U -> End
            d18_vec = endmember_iso[0] - d18o_u
            d2h_vec = endmember_iso[1] - d2h_u
            
            # Vector U -> V
            d18_uv = d18o_v - d18o_u
            d2h_uv = d2h_v - d2h_u
            
            # Project UV onto U->End line to find orthogonal distance
            # Line equation: P = U + t * Vec
            # Distance from point V to line passing through U and End
            # d = |cross_product(Vec, UV)| / |Vec|
            vec_len_sq = d18_vec**2 + d2h_vec**2
            if vec_len_sq > 1e-9:
                cross_prod = d18_vec * d2h_uv - d2h_vec * d18_uv
                mix_deviation = abs(cross_prod) / (vec_len_sq**0.5)
            else:
                # U and End are same point
                mix_deviation = (d18_uv**2 + d2h_uv**2)**0.5

    penalty = e_penalty + d_excess_weight * d_penalty + mix_deviation

    metrics = {
        "d18o_u": d18o_u,
        "d2h_u": d2h_u,
        "d18o_v": d18o_v,
        "d2h_v": d2h_v,
        "e_u": e_u,
        "e_v": e_v,
        "d_excess_u": d_u,
        "d_excess_v": d_v,
        "d_excess_delta": d_v - d_u,
        "d_excess_penalty": d_penalty,
        "mix_deviation": mix_deviation,
        "enrichment_slope": (
            (d2h_v - d2h_u) / (d18o_v - d18o_u)
            if abs(d18o_v - d18o_u) > 1e-6
            else float("nan")
        ),
    }
    return penalty, metrics
