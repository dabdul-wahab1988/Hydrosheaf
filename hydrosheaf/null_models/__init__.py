"""Null models for non-connectivity explanations.

Each null model computes a score in [0, 1] representing the plausibility
that observed chemical or isotopic similarity between two wells could arise
without a direct groundwater flow connection.

Higher null_score -> more plausible no-flow explanation -> chemistry/isotope
evidence for a flow edge should be downgraded.
"""

from __future__ import annotations

from typing import List, Mapping, Tuple

from ..config import Config
from ..log import get_logger

logger = get_logger("null_models")


def compute_null_penalty(
    sample_a: Mapping[str, object],
    sample_b: Mapping[str, object],
    config: Config,
) -> Tuple[float, List[str]]:
    """Compute combined null-model score and flags for an edge pair.

    Returns
    -------
    null_score : float in [0, 1]
        Plausibility that similarity arises without direct flow.
    flags : List[str]
        Descriptive flags for each active null explanation.
    """
    flags: List[str] = []
    total_score = 0.0

    # --- Chemistry null ---
    try:
        from .chemistry import chemistry_null_score
        chem_score, chem_flags = chemistry_null_score(sample_a, sample_b, config)
        if chem_score > 0.0:
            total_score += chem_score
            flags.extend(chem_flags)
    except Exception:
        logger.warning("Chemistry null model failed for edge pair; flagging.")
        flags.append("null_chemistry_error")

    # --- Lithology null ---
    try:
        from .lithology import lithology_null_score
        lith_score, lith_flags = lithology_null_score(sample_a, sample_b, config)
        if lith_score > 0.0:
            total_score += lith_score * float(getattr(config, "null_lithology_weight", 0.3))
            flags.extend(lith_flags)
    except Exception:
        logger.warning("Lithology null model failed for edge pair; flagging.")
        flags.append("null_lithology_error")

    # --- Endmember null ---
    try:
        from .endmembers import endmember_null_score
        endm_score, endm_flags = endmember_null_score(sample_a, sample_b, config)
        if endm_score > 0.0:
            total_score += endm_score * float(getattr(config, "null_endmember_weight", 0.4))
            flags.extend(endm_flags)
    except Exception:
        logger.warning("Endmember null model failed for edge pair; flagging.")
        flags.append("null_endmember_error")

    return min(total_score, 1.0), flags
