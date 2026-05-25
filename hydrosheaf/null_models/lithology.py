"""Lithology-based null model.

If two wells are in the same aquifer layer or share common lithology,
their chemical similarity could arise from water-rock interaction with
the same mineral assemblage — not necessarily direct flow between them.
"""

from __future__ import annotations

from typing import List, Mapping, Tuple

from ..config import Config


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
