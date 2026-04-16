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
    Fit LMWL using provided samples.
    
    If sample_type_col is provided, only samples where sample[sample_type_col]
    is in meteoric_types (case-insensitive) are used.
    
    If weight_col is provided and present in samples, performs precipitation-weighted 
    regression (PWLMWL), which is scientifically preferred for meteoric lines.
    """
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
                # If weights are requested but missing/zero, skip or treat as minimal?
                # Standard practice: skip or negligible weight. Let's use negligible.
                w = 1e-6

        values.append(pair)
        weights.append(w)

    if not values:
        if sample_type_col:
            raise ValueError(
                f"No precipitation samples found (looking for types {valid_types} in column '{sample_type_col}') to fit LMWL. "
                "Ensure your data contains rain samples or use --use-gmwl."
            )
        else:
            raise ValueError("No isotope values available to fit LMWL.")

    xs = [item[0] for item in values]
    ys = [item[1] for item in values]
    
    # Weighted means
    sum_w = sum(weights)
    if sum_w <= 0:
        sum_w = 1.0
        weights = [1.0] * len(weights)
        
    x_mean = sum(x * w for x, w in zip(xs, weights)) / sum_w
    y_mean = sum(y * w for y, w in zip(ys, weights)) / sum_w
    
    # Weighted covariance and variance
    # Slope b = Sum(w * (x - x_bar) * (y - y_bar)) / Sum(w * (x - x_bar)^2)
    
    numer = sum(w * (x - x_mean) * (y - y_mean) for x, y, w in zip(xs, ys, weights))
    denom = sum(w * (x - x_mean) ** 2 for x, w in zip(xs, weights))
    
    if abs(denom) < 1e-9:
        # Insufficient x-variance to fit a line (vertical line or single point)
        # Fallback: return a flat line through the mean, or raise?
        # Raising is safer for scientific integrity.
        raise ValueError(
            "Insufficient d18O variance in precipitation samples to fit a valid LMWL. "
            "Samples may be identical or too clustered. Use --use-gmwl or fixed parameters."
        )

    slope = numer / denom
    intercept = y_mean - slope * x_mean
    return intercept, slope



def isotope_penalty(
    d18o_u: float,
    d2h_u: float,
    d18o_v: float,
    d2h_v: float,
    a: float,
    b: float,
    mode: str,
    d_excess_weight: float = 0.0,
) -> Tuple[float, dict]:
    e_u = evaporation_index(d18o_u, d2h_u, a, b)
    e_v = evaporation_index(d18o_v, d2h_v, a, b)
    d_u = compute_d_excess(d18o_u, d2h_u)
    d_v = compute_d_excess(d18o_v, d2h_v)

    if mode == "evap":
        e_penalty = max(0.0, abs(e_u) - abs(e_v))
        d_penalty = max(0.0, d_v - d_u)
    else:
        e_penalty = max(0.0, abs(e_v) - abs(e_u))
        d_penalty = 0.0

    penalty = e_penalty + d_excess_weight * d_penalty

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
        "enrichment_slope": (
            (d2h_v - d2h_u) / (d18o_v - d18o_u)
            if abs(d18o_v - d18o_u) > 1e-6
            else float("nan")
        ),
    }
    return penalty, metrics
