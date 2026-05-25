"""Endmember and recharge-based null model.

Similar isotope signatures can arise from shared recharge, common
endmember mixing, regional evaporation/isotope trends, spatial
autocorrelation, or common anthropogenic sources — without direct
groundwater flow between two specific wells.
"""

from __future__ import annotations

from typing import List, Mapping, Optional, Tuple

from ..config import Config
from ..log import get_logger

logger = get_logger("null_models.endmembers")


def _maybe_float(value: object) -> Optional[float]:
    """Safe float conversion, returns None on failure."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def endmember_null_score(
    sample_a: Mapping[str, object],
    sample_b: Mapping[str, object],
    config: Config,
) -> Tuple[float, List[str]]:
    """Compute null-model score from shared recharge/endmember evidence.

    Checks:
    1. Isotope proximity to LMWL (shared precipitation source)
    2. Spatial proximity (nearby wells may share recharge)
    3. Common anthropogenic source flags (nitrate, chloride)

    Returns (null_score, flags).
    """
    flags: List[str] = []
    null_score = 0.0

    d18o_key = getattr(config, "isotope_d18o_key", "18O")
    d2h_key = getattr(config, "isotope_d2h_key", "2H")

    d18o_a = _maybe_float(sample_a.get(d18o_key))
    d2h_a = _maybe_float(sample_a.get(d2h_key))
    d18o_b = _maybe_float(sample_b.get(d18o_key))
    d2h_b = _maybe_float(sample_b.get(d2h_key))

    # --- Shared recharge via isotope similarity ---
    if d18o_a is not None and d18o_b is not None and d2h_a is not None and d2h_b is not None:
        lmwl_a = float(getattr(config, "lmwl_a", 8.66))
        lmwl_b = float(getattr(config, "lmwl_b", 7.22))

        # Compute departures from LMWL (deuterium excess proxy)
        dex_a = d2h_a - (lmwl_a + lmwl_b * d18o_a)
        dex_b = d2h_b - (lmwl_a + lmwl_b * d18o_b)

        # If both wells have similar LMWL departure (within 3 per mil)
        # they may share the same recharge source
        if abs(dex_a - dex_b) < 3.0:
            null_score += 0.3
            flags.append("null_shared_recharge")

        # Isotope proximity: if d18O values are within 0.5 per mil
        if abs(d18o_a - d18o_b) < 0.5:
            null_score += 0.2
            flags.append("null_isotope_proximity")

    # --- Spatial autocorrelation ---
    spatial_weight = float(getattr(config, "null_spatial_weight", 0.2))
    lat_a = _maybe_float(sample_a.get("lat"))
    lon_a = _maybe_float(sample_a.get("lon"))
    lat_b = _maybe_float(sample_b.get("lat"))
    lon_b = _maybe_float(sample_b.get("lon"))

    if lat_a is not None and lon_a is not None and lat_b is not None and lon_b is not None:
        # Approximate distance in degrees (crude but sufficient for flagging)
        dist_deg = ((lat_a - lat_b) ** 2 + (lon_a - lon_b) ** 2) ** 0.5
        # If within ~0.01 degrees (~1 km), spatial autocorrelation is plausible
        if dist_deg < 0.01:
            null_score += spatial_weight * 0.5
            flags.append("null_spatial_autocorr")
        elif dist_deg < 0.05:  # ~5 km
            null_score += spatial_weight * 0.2

    # --- Common anthropogenic source ---
    anthro_weight = float(getattr(config, "null_anthropogenic_weight", 0.2))
    no3_a = _maybe_float(sample_a.get("NO3"))
    no3_b = _maybe_float(sample_b.get("NO3"))
    cl_a = _maybe_float(sample_a.get("Cl"))
    cl_b = _maybe_float(sample_b.get("Cl"))

    anthro_hits = 0
    # Both have elevated nitrate (>10 mg/L proxy): could be common fertilizer
    if no3_a is not None and no3_b is not None:
        if no3_a > 10.0 and no3_b > 10.0 and abs(no3_a - no3_b) / max(no3_a + no3_b, 0.1) < 0.5:
            anthro_hits += 1

    # Both have elevated chloride: could be common salinity source
    if cl_a is not None and cl_b is not None:
        if cl_a > 50.0 and cl_b > 50.0 and abs(cl_a - cl_b) / max(cl_a + cl_b, 0.1) < 0.5:
            anthro_hits += 1

    if anthro_hits >= 1:
        null_score += anthro_weight
        flags.append("null_common_anthropogenic")

    return min(null_score, 1.0), flags
