"""Edge-level evidence classification for assumption-audited topology.

Distinct from `claims.py` EvidenceLevel, which governs manuscript wording.
This module classifies individual candidate edges by the strength and
independence of their supporting observations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import List, Optional, Tuple

from ..config import Config


class EdgeEvidenceClass(IntEnum):
    """Ordered edge evidence classes.

    Higher values mean stronger, more independent support.
    These are stamped onto every final edge's attributes.
    """

    FALSIFIED = 0
    """Rejected by null model, benchmark, tracer contradiction, or impossible
    head/age relation."""

    AMBIGUOUS = 1
    """Insufficient independent evidence to confirm or reject flow."""

    PRIOR_ASSISTED = 2
    """Supported mainly by an imported model, MODPATH prior, or head prior
    without corroborating chemistry/isotope agreement."""

    PROBABLE = 3
    """Multiple independent evidence streams support the edge and no strong
    null explanation dominates."""

    VALIDATED = 4
    """Supported by an independent benchmark, tracer test, or field-confirmed
    topology."""


@dataclass
class EdgeEvidenceAnnotation:
    """Evidence classification for a single candidate edge."""

    evidence_class: EdgeEvidenceClass
    evidence_score: float
    """Overall confidence score in [0, 1]. Higher = more confident."""

    null_score: float
    """Plausibility that similarity arose without direct flow [0, 1]."""

    evidence_flags: List[str] = field(default_factory=list)
    """Detailed diagnostic flags (e.g. iso_missing, null_chemistry)."""

    reason: str = ""
    """Short human-readable reason for the classification."""


def _evidence_score(local_score: float) -> float:
    """Map a local cost to an evidence confidence score in [0, 1].

    Lower local_score (better fit) → higher evidence_score.
    """
    return 1.0 / (1.0 + max(0.0, local_score))


def classify_edge_evidence(
    local_score: float,
    null_score: float,
    flags: Optional[List[str]] = None,
    config: Optional[Config] = None,
    validated_by: str = "",
) -> Tuple[str, str]:
    """Classify an edge into an evidence class based on score, null score, and flags.

    Parameters
    ----------
    local_score : float
        The combined edge cost (prior + iso + cl + age). Lower is better.
    null_score : float in [0, 1]
        Null-model plausibility. Higher means similarity could be no-flow.
    flags : list of str, optional
        Diagnostic flags from scoring (e.g. 'iso_missing', 'age_reversal').
    config : Config, optional
        For thresholds. Defaults apply when None.
    validated_by : str
        If non-empty, names the independent validation source (e.g. 'modpath',
        'tracer_test', 'field'). Required for VALIDATED classification.

    Returns
    -------
    (class_name, reason) : (str, str)
        Name of the EdgeEvidenceClass and a short human-readable reason.
    """
    evidence_score = _evidence_score(local_score)
    flag_set = set(flags or [])

    threshold_probable = (
        float(getattr(config, "evidence_threshold_probable", 0.6))
        if config
        else 0.6
    )
    threshold_validated = (
        float(getattr(config, "evidence_threshold_validated", 0.8))
        if config
        else 0.8
    )

    # --- FALSIFIED gates ---
    if "age_reversal" in flag_set:
        return (
            EdgeEvidenceClass.FALSIFIED.name,
            "age_reversal: downstream water is younger than upstream",
        )
    if null_score > 0.8:
        return (
            EdgeEvidenceClass.FALSIFIED.name,
            f"null_model_score={null_score:.2f} exceeds 0.8: similarity plausibly no-flow",
        )
    if evidence_score < 0.1:
        return (
            EdgeEvidenceClass.FALSIFIED.name,
            f"evidence_score={evidence_score:.3f} below minimum",
        )

    # --- AMBIGUOUS gates ---
    if "iso_missing_u" in flag_set or "iso_missing_v" in flag_set:
        return (
            EdgeEvidenceClass.AMBIGUOUS.name,
            "isotope data missing on one or both nodes",
        )
    if "missing_evidence" in flag_set:
        return (
            EdgeEvidenceClass.AMBIGUOUS.name,
            "key chemical/isotope observations missing",
        )
    if null_score > 0.5:
        return (
            EdgeEvidenceClass.AMBIGUOUS.name,
            f"null_model_score={null_score:.2f} between 0.5-0.8: competing no-flow explanation",
        )
    if evidence_score < threshold_probable:
        return (
            EdgeEvidenceClass.AMBIGUOUS.name,
            f"evidence_score={evidence_score:.3f} below PROBABLE threshold {threshold_probable}",
        )

    # --- PRIOR_ASSISTED ---
    corroborating = any(
        f in flag_set
        for f in ("evap_candidate", "age_consistent")
    )
    if not corroborating and evidence_score < threshold_validated:
        return (
            EdgeEvidenceClass.PRIOR_ASSISTED.name,
            "head/distance prior dominates; no corroborating chemistry or isotope signal",
        )

    # --- VALIDATED: requires external confirmation ---
    if validated_by:
        if evidence_score >= threshold_validated:
            return (
                EdgeEvidenceClass.VALIDATED.name,
                f"validated_by={validated_by} with evidence_score={evidence_score:.3f}",
            )
        return (
            EdgeEvidenceClass.PROBABLE.name,
            f"validation source '{validated_by}' present but evidence_score={evidence_score:.3f} below threshold; classified PROBABLE",
        )

    return (
        EdgeEvidenceClass.PROBABLE.name,
        f"evidence_score={evidence_score:.3f} supports connectivity; no external validation provided",
    )
