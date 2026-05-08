"""Topology validation helpers for Hydrosheaf graph manuscripts.

M4 uses two scientifically distinct modes:

1. independent graph inference compared against a MODPATH reference graph;
2. MODPATH-informed graph priors used inside Hydrosheaf.

These helpers keep those modes separate so validation reports do not overstate
evidence from graph priors as independent predictive performance.
"""
from __future__ import annotations

import math
from statistics import median
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from ..graph.types import Edge
from ..physics.priors import PhysicsPrior, apply_physics_priors


DirectedEdge = Tuple[str, str]


def _safe_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _edge_tuple(edge: Any) -> Optional[DirectedEdge]:
    if isinstance(edge, Edge):
        return str(edge.u), str(edge.v)
    if isinstance(edge, PhysicsPrior):
        return str(edge.u), str(edge.v)
    if isinstance(edge, Mapping):
        u = edge.get("u") or edge.get("upstream") or edge.get("source") or edge.get("from")
        v = edge.get("v") or edge.get("downstream") or edge.get("target") or edge.get("to")
        if u is None or v is None:
            edge_id = str(edge.get("edge_id") or "")
            if "->" in edge_id:
                u, v = edge_id.split("->", 1)
        if u is None or v is None:
            return None
        return str(u), str(v)
    if isinstance(edge, Sequence) and not isinstance(edge, (str, bytes)) and len(edge) >= 2:
        return str(edge[0]), str(edge[1])
    return None


def normalize_directed_edges(edges: Iterable[Any]) -> Set[DirectedEdge]:
    """Return a set of ``(u, v)`` directed edges from common Hydrosheaf formats."""

    out: Set[DirectedEdge] = set()
    for edge in edges:
        parsed = _edge_tuple(edge)
        if parsed is None:
            continue
        u, v = parsed
        if u and v and u != v:
            out.add((u, v))
    return out


def edge_confusion(
    reference_edges: Iterable[Any],
    inferred_edges: Iterable[Any],
    *,
    candidate_edges: Optional[Iterable[Any]] = None,
) -> Dict[str, Any]:
    """Compute directed-edge confusion diagnostics.

    ``candidate_edges`` should contain the possible edge universe when true
    negatives or false-positive rate are needed. Without it, the report still
    gives TP/FP/FN, precision, recall, F1, and false-negative rate.
    """

    reference = normalize_directed_edges(reference_edges)
    inferred = normalize_directed_edges(inferred_edges)
    universe = normalize_directed_edges(candidate_edges or []) if candidate_edges is not None else set()
    if candidate_edges is None:
        universe = reference | inferred
    else:
        universe |= reference | inferred

    true_positive = reference & inferred
    false_positive = inferred - reference
    false_negative = reference - inferred
    true_negative = universe - reference - inferred

    tp = len(true_positive)
    fp = len(false_positive)
    fn = len(false_negative)
    tn = len(true_negative) if candidate_edges is not None else None
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    fpr = fp / (fp + tn) if tn is not None and fp + tn else float("nan")
    fnr = fn / (fn + tp) if fn + tp else 0.0

    return {
        "reference_edges": sorted(reference),
        "inferred_edges": sorted(inferred),
        "candidate_edges": sorted(universe),
        "true_positives": sorted(true_positive),
        "false_positives": sorted(false_positive),
        "false_negatives": sorted(false_negative),
        "true_negatives": sorted(true_negative) if candidate_edges is not None else [],
        "n_reference_edges": float(len(reference)),
        "n_inferred_edges": float(len(inferred)),
        "tp": float(tp),
        "fp": float(fp),
        "fn": float(fn),
        "tn": float(tn) if tn is not None else float("nan"),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "false_positive_rate": fpr,
        "false_negative_rate": fnr,
    }


def _edge_length_lookup(edge_lengths: Optional[Mapping[Any, Any]]) -> Dict[DirectedEdge, float]:
    out: Dict[DirectedEdge, float] = {}
    for key, value in dict(edge_lengths or {}).items():
        parsed = _edge_tuple(key)
        if parsed is None and isinstance(key, str) and "->" in key:
            parsed = tuple(key.split("->", 1))  # type: ignore[assignment]
        length = _safe_float(value)
        if parsed is not None and length is not None and length > 0:
            out[(str(parsed[0]), str(parsed[1]))] = length
    return out


def scale_mismatch_diagnostics(
    reference_edges: Iterable[Any],
    inferred_edges: Iterable[Any],
    *,
    edge_lengths: Optional[Mapping[Any, Any]] = None,
    mismatch_ratio_threshold: float = 2.0,
) -> Dict[str, Any]:
    """Check whether inferred and reference edges operate at different scales."""

    lengths = _edge_length_lookup(edge_lengths)
    reference = normalize_directed_edges(reference_edges)
    inferred = normalize_directed_edges(inferred_edges)
    ref_lengths = [lengths[edge] for edge in reference if edge in lengths]
    inf_lengths = [lengths[edge] for edge in inferred if edge in lengths]
    if not ref_lengths or not inf_lengths:
        return {
            "has_length_data": False,
            "median_reference_length": float("nan"),
            "median_inferred_length": float("nan"),
            "max_inferred_length": float("nan"),
            "length_ratio": float("nan"),
            "max_inferred_to_reference_median_ratio": float("nan"),
            "scale_mismatch": False,
        }
    ref_median = float(median(ref_lengths))
    inf_median = float(median(inf_lengths))
    max_inferred = float(max(inf_lengths))
    ratio = max(ref_median, inf_median) / max(min(ref_median, inf_median), 1.0e-12)
    shortcut_ratio = max_inferred / max(ref_median, 1.0e-12)
    return {
        "has_length_data": True,
        "median_reference_length": ref_median,
        "median_inferred_length": inf_median,
        "max_inferred_length": max_inferred,
        "length_ratio": ratio,
        "max_inferred_to_reference_median_ratio": shortcut_ratio,
        "scale_mismatch": ratio >= mismatch_ratio_threshold or shortcut_ratio >= mismatch_ratio_threshold,
    }


def validate_independent_graph_against_modpath(
    inferred_edges: Iterable[Any],
    modpath_reference_edges: Iterable[Any],
    *,
    candidate_edges: Optional[Iterable[Any]] = None,
    edge_lengths: Optional[Mapping[Any, Any]] = None,
) -> Dict[str, Any]:
    """Validate independently inferred Hydrosheaf topology against MODPATH."""

    confusion = edge_confusion(
        modpath_reference_edges,
        inferred_edges,
        candidate_edges=candidate_edges,
    )
    scale = scale_mismatch_diagnostics(
        modpath_reference_edges,
        inferred_edges,
        edge_lengths=edge_lengths,
    )
    return {
        "validation_mode": "independent_graph_inference",
        "reference": "MODPATH advective connectivity",
        "metrics": confusion,
        "scale_mismatch": scale,
        "failure_modes": {
            "false_positive_edges": confusion["false_positives"],
            "false_negative_edges": confusion["false_negatives"],
            "scale_mismatch": scale["scale_mismatch"],
        },
        "claim_guardrail": (
            "These metrics support only reduced-order graph reproduction of the "
            "specified MODPATH benchmark, not general multi-archive MODFLOW/MODPATH performance."
        ),
    }


def build_modpath_informed_graph_priors(
    modpath_edges: Iterable[Any],
    *,
    default_p_uv: float = 1.0,
    travel_time_days: Optional[Mapping[Any, Any]] = None,
    source: str = "modpath_informed_prior",
) -> List[PhysicsPrior]:
    """Convert MODPATH connectivity into Hydrosheaf graph priors.

    This is a prior-construction mode, not independent validation.
    """

    times = _edge_length_lookup(travel_time_days)
    priors: List[PhysicsPrior] = []
    for u, v in sorted(normalize_directed_edges(modpath_edges)):
        tau = times.get((u, v))
        priors.append(
            PhysicsPrior(
                u=u,
                v=v,
                p_uv=float(default_p_uv),
                tt_mean_days=tau,
                source=source,
            )
        )
    return priors


def apply_modpath_informed_graph_priors(
    hydrosheaf_edges: Iterable[Any],
    modpath_edges: Iterable[Any],
    *,
    mode: str = "merge",
    default_p_uv: float = 1.0,
    travel_time_days: Optional[Mapping[Any, Any]] = None,
) -> Dict[str, Any]:
    """Apply MODPATH-informed priors while reporting that this is not validation."""

    base_edges: List[Edge] = []
    for u, v in sorted(normalize_directed_edges(hydrosheaf_edges)):
        base_edges.append(Edge(edge_id=f"{u}->{v}", u=u, v=v))
    priors = build_modpath_informed_graph_priors(
        modpath_edges,
        default_p_uv=default_p_uv,
        travel_time_days=travel_time_days,
    )
    updated_edges = apply_physics_priors(base_edges, priors, mode=mode)
    return {
        "validation_mode": "modpath_informed_graph_prior",
        "mode": mode,
        "n_input_hydrosheaf_edges": len(base_edges),
        "n_modpath_prior_edges": len(priors),
        "n_output_edges": len(updated_edges),
        "prior_edges": [prior.edge_id() for prior in priors],
        "output_edges": [edge.edge_id for edge in updated_edges],
        "not_independent_validation": True,
        "claim_guardrail": (
            "Do not use these results as independent topology accuracy. MODPATH has "
            "entered the Hydrosheaf graph as prior information."
        ),
    }
