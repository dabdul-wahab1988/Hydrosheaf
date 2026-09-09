"""Lithology-based null model.

If two wells are in the same aquifer layer or share common lithology,
their chemical similarity could arise from water-rock interaction with
the same mineral assemblage — not necessarily direct flow between them.
"""

from __future__ import annotations

from typing import List, Mapping, Optional, Tuple

from ..config import Config


def _text(sample: Mapping[str, object], key: str) -> Optional[str]:
    value = sample.get(key)
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return None
    return text.lower()


def mapped_geology_similarity(
    sample_a: Mapping[str, object],
    sample_b: Mapping[str, object],
) -> Tuple[Optional[float], List[str]]:
    """Compare mapped surface-geology metadata as a *null* covariate.

    The score is deliberately graded: an exact code/symbol match is 1, a
    shared stratigraphic unit is 0.75, and a shared tectonic domain is 0.5.
    ``None`` is returned when either record is unmatched/unknown.  This is a
    descriptive alternative explanation for chemistry similarity, never a
    hard candidate-edge filter and never a substitute for screened-interval
    lithology or hydraulic truth.
    """

    status_a = (_text(sample_a, "geology_join_status") or "").upper()
    status_b = (_text(sample_b, "geology_join_status") or "").upper()
    invalid = {"", "NO_MATCH", "UNMATCHED", "BOUNDARY_REVIEW", "OVERLAP_REVIEW"}
    if status_a in invalid or status_b in invalid:
        return None, ["mapped_geology_unknown"]

    flags: List[str] = []
    code_a = _text(sample_a, "geology_code_1000")
    code_b = _text(sample_b, "geology_code_1000")
    symbol_a = _text(sample_a, "geology_symbol")
    symbol_b = _text(sample_b, "geology_symbol")
    if code_a and code_b and code_a == code_b:
        flags.append("null_mapped_geology_exact")
        return 1.0, flags
    if symbol_a and symbol_b and symbol_a == symbol_b:
        flags.append("null_mapped_geology_symbol")
        return 0.9, flags

    unit_a = _text(sample_a, "geology_stratigraphic_unit")
    unit_b = _text(sample_b, "geology_stratigraphic_unit")
    if unit_a and unit_b and unit_a == unit_b:
        flags.append("null_mapped_geology_unit")
        return 0.75, flags

    domain_a = _text(sample_a, "geology_tectonic_domain")
    domain_b = _text(sample_b, "geology_tectonic_domain")
    if domain_a and domain_b and domain_a == domain_b:
        flags.append("null_mapped_geology_domain")
        return 0.5, flags
    return 0.0, ["mapped_geology_different_or_unresolved"]


def lithology_null_score(
    sample_a: Mapping[str, object],
    sample_b: Mapping[str, object],
    config: Config,
) -> Tuple[float, List[str]]:
    """Compute null-model score from common lithology / aquifer layer.

    Returns (null_score, flags).
    """
    flags: List[str] = []

    layer_key = getattr(config, "layer_key", "aquifer_layer")
    aquifer_key = getattr(config, "edge_aquifer_key", "aquifer_unit")
    lithology_key = "lithology"

    layer_a = sample_a.get(layer_key) or sample_a.get(aquifer_key)
    layer_b = sample_b.get(layer_key) or sample_b.get(aquifer_key)

    # Check for common lithology tag
    lith_a = sample_a.get(lithology_key)
    lith_b = sample_b.get(lithology_key)

    null_score = 0.0

    if layer_a is not None and layer_b is not None:
        if str(layer_a) == str(layer_b):
            null_score = max(null_score, 0.5)
            flags.append("null_common_lithology")

    if lith_a is not None and lith_b is not None:
        if str(lith_a).strip().lower() == str(lith_b).strip().lower():
            null_score = max(null_score, 0.7)
            flags.append("null_common_lithology_explicit")

    return null_score, flags


__all__ = ["lithology_null_score", "mapped_geology_similarity"]
