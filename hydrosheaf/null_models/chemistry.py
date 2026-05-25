"""Chemistry-similarity null model.

Computes the probability/cost that two wells share similar major-ion
chemistry without a direct groundwater flow connection.

Rationale: similar chemistry can arise from shared lithology,
regional groundwater evolution, or common recharge — not necessarily
direct flow between two specific wells.
"""

from __future__ import annotations

import math
from typing import List, Mapping, Tuple

from ..config import Config
from ..log import get_logger

logger = get_logger("null_models.chemistry")


def _extract_ion_vector(
    sample: Mapping[str, object],
    config: Config,
) -> list[float]:
    """Extract major-ion concentrations as a numeric vector."""
    ions = getattr(config, "ion_order", ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"])
    vec = []
    for ion in ions:
        val = sample.get(ion)
        if val is None:
            vec.append(0.0)
        else:
            try:
                vec.append(float(val))
            except (TypeError, ValueError):
                vec.append(0.0)
    return vec


def _normalize_vector(vec: list[float]) -> list[float]:
    """Normalize to unit sum for compositional comparison."""
    total = sum(vec)
    if total <= 0.0:
        return vec
    return [v / total for v in vec]


def chemistry_null_score(
    sample_a: Mapping[str, object],
    sample_b: Mapping[str, object],
    config: Config,
) -> Tuple[float, List[str]]:
    """Compute the null-model score from chemistry similarity.

    Returns (null_score, flags) where:
    - null_score in [0, 1]: higher means the chemistry similarity
      is more plausibly explained by common source/lithology
      rather than direct flow.
    - flags: descriptive labels.
    """
    flags: List[str] = []

    vec_a = _extract_ion_vector(sample_a, config)
    vec_b = _extract_ion_vector(sample_b, config)

    # Check if we have enough chemistry data
    nonzero_a = sum(1 for v in vec_a if v > 0)
    nonzero_b = sum(1 for v in vec_b if v > 0)
    if nonzero_a < 2 or nonzero_b < 2:
        # Too few ions for meaningful comparison
        return 0.0, flags

    norm_a = _normalize_vector(vec_a)
    norm_b = _normalize_vector(vec_b)

    # Euclidean distance on normalized composition
    sq_diff = sum((a - b) ** 2 for a, b in zip(norm_a, norm_b))
    euclidean_dist = math.sqrt(sq_diff)

    # Also compute cosine similarity (1 = identical direction)
    dot = sum(a * b for a, b in zip(norm_a, norm_b))
    norm_a_mag = math.sqrt(sum(a * a for a in norm_a))
    norm_b_mag = math.sqrt(sum(b * b for b in norm_b))
    if norm_a_mag > 0 and norm_b_mag > 0:
        cos_sim = dot / (norm_a_mag * norm_b_mag)
    else:
        cos_sim = 0.0

    threshold = float(getattr(config, "null_chemistry_similarity_threshold", 0.3))

    # If chemistry is very similar (low Euclidean distance, high cosine similarity)
    # the null score rises — similar chemistry could be regional, not flow-driven
    if euclidean_dist < threshold and cos_sim > (1.0 - threshold):
        # Linear ramp: closer = higher null score
        null_score = 1.0 - (euclidean_dist / threshold)
        null_score = max(0.0, min(1.0, null_score))
        flags.append("null_chemistry_similar")
        return null_score, flags

    return 0.0, flags
