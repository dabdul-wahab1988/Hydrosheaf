"""Isotope-derived metrics for sheaf scoring."""

from dataclasses import dataclass
from typing import Iterable, Mapping, Optional, Tuple
import math

import numpy as np

from ..config import Config
from ..data.schema import parse_numeric
from ..isotopes import compute_d_excess, evaporation_index, extract_isotopes


@dataclass
class IsotopeStats:
    d_excess_p25: float = 10.0
    evap_index_p25: float = 0.0
    evap_index_p75: float = 2.0


def compute_isotope_stats(
    samples: Iterable[Mapping[str, object]],
    config: Config,
) -> IsotopeStats:
    d_excess_vals = []
    evap_vals = []
    for sample in samples:
        isotopes = extract_isotopes(
            sample,
            d18o_key=config.isotope_d18o_key,
            d2h_key=config.isotope_d2h_key,
        )
        if isotopes is None:
            continue
        d18o, d2h = isotopes
        d_excess_vals.append(compute_d_excess(d18o, d2h))
        evap_vals.append(evaporation_index(d18o, d2h, config.lmwl_a, config.lmwl_b))

    stats = IsotopeStats()
    if len(d_excess_vals) >= 3:
        stats.d_excess_p25 = float(np.quantile(d_excess_vals, 0.25))
    if len(evap_vals) >= 3:
        stats.evap_index_p25 = float(np.quantile(evap_vals, 0.25))
        stats.evap_index_p75 = float(np.quantile(evap_vals, 0.75))
    if stats.evap_index_p75 < stats.evap_index_p25 + 1e-6:
        stats.evap_index_p75 = stats.evap_index_p25 + 1.0
    return stats


def sample_depth_m(sample: Mapping[str, object], config: Config) -> Optional[float]:
    detection_policy = config.detection_limit_policy
    top_key = getattr(config, "screen_top_key", "screen_top")
    bottom_key = getattr(config, "screen_bottom_key", "screen_bottom")
    top = parse_numeric(sample.get(top_key), detection_policy)
    bottom = parse_numeric(sample.get(bottom_key), detection_policy)
    if top is not None and bottom is not None:
        return 0.5 * (float(top) + float(bottom))
    if top is not None:
        return float(top)
    if bottom is not None:
        return float(bottom)

    depth_key = getattr(config, "edge_screen_depth_key", "screen_depth")
    depth = parse_numeric(sample.get(depth_key), detection_policy)
    if depth is not None:
        return float(depth)

    well_depth_key = getattr(config, "edge_well_depth_key", "well_depth")
    depth = parse_numeric(sample.get(well_depth_key), detection_policy)
    if depth is not None:
        return float(depth)

    z_key = getattr(config, "z_coordinate_key", None)
    if z_key:
        z_val = parse_numeric(sample.get(z_key), detection_policy)
        if z_val is not None and getattr(config, "z_coordinate_positive_down", True):
            return float(z_val)
    return None


def compute_evaporation_probability(
    evap_index: Optional[float],
    d_excess: Optional[float],
    depth_m: Optional[float],
    stats: IsotopeStats,
    config: Config,
) -> float:
    if evap_index is None or d_excess is None:
        return 0.0
    evap_mag = abs(evap_index)
    evap_scale = max(stats.evap_index_p75 - stats.evap_index_p25, 1e-6)
    evap_score = (evap_mag - stats.evap_index_p25) / evap_scale

    d_excess_scale = max(abs(stats.d_excess_p25), 1.0)
    d_excess_score = (stats.d_excess_p25 - d_excess) / d_excess_scale

    depth_score = 0.0
    shallow_depth = float(getattr(config, "sheaf_shallow_depth_m", 30.0))
    if depth_m is not None and shallow_depth > 0:
        depth_score = (shallow_depth - depth_m) / shallow_depth

    combined = 0.6 * evap_score + 0.3 * d_excess_score + 0.2 * depth_score
    strength = float(getattr(config, "sheaf_evap_gate_strength", 1.0))
    base = 1.0 / (1.0 + math.exp(-strength * combined))
    if depth_m is not None and shallow_depth > 0 and depth_m > shallow_depth:
        depth_gate = shallow_depth / depth_m
        base *= max(0.0, min(1.0, depth_gate))
    return base
