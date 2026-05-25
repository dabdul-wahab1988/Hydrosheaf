"""Causal discovery layer for groundwater flow direction support."""

from .discovery import (
    CausalResult,
    assess_edge_causality,
    compute_causal_support,
)

__all__ = ["CausalResult", "assess_edge_causality", "compute_causal_support"]
