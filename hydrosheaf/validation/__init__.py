"""Validation helpers for manuscript-grade Hydrosheaf evidence.

The validation package keeps benchmark metrics and claim guardrails close to
the code that produces Hydrosheaf results. It is intentionally lightweight so
benchmark scripts and manuscript tables can share the same definitions.
"""

from .claims import (
    ClaimRecord,
    EvidenceLevel,
    assess_claim_records,
    evidence_level_allows,
)
from .metrics import (
    classification_metrics,
    interval_coverage,
    regression_metrics,
    topology_metrics,
)

__all__ = [
    "ClaimRecord",
    "EvidenceLevel",
    "assess_claim_records",
    "classification_metrics",
    "evidence_level_allows",
    "interval_coverage",
    "regression_metrics",
    "topology_metrics",
]
