"""Isotope utilities (delta-18O / delta-2H)."""

from typing import Iterable, Mapping, Optional, Tuple

from .data.schema import parse_numeric


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
) -> Tuple[float, float]:
    values = []
    for sample in samples:
        pair = extract_isotopes(sample, d18o_key=d18o_key, d2h_key=d2h_key)
        if pair is None:
            continue
        values.append(pair)

    if not values:
        raise ValueError("No isotope values available to fit LMWL.")

    xs = [item[0] for item in values]
    ys = [item[1] for item in values]
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    denom = sum((x - x_mean) ** 2 for x in xs)
    if denom == 0:
        return 0.0, y_mean
    slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denom
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
        "enrichment_slope": (d2h_v - d2h_u) / (d18o_v - d18o_u) if abs(d18o_v - d18o_u) > 1e-6 else float("nan"),
    }
    return penalty, metrics
