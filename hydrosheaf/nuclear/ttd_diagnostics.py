"""Quality gates, physical consistency audits, and calibrated abstention diagnostics for network TTDs.

This module enforces pre-registered gating rules from TTD_GRAPH_EXTENSION_PROTOCOL:
1. Physical tracer conflict gates (e.g. modern gas with dead 14C and zero 3H);
2. Forward matrix condition number and effective rank;
3. Graph acyclicity and hydraulic head gradient consistency;
4. Machine-readable abstention reason codes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping, Optional, Tuple

import networkx as nx

from .ttd_kernel_builder import MultiTracerForwardSystem


# Pre-registered standard abstention reason codes
CODE_NO_TRACERS = "NO_TRACERS_OBSERVED"
CODE_INSUFFICIENT_TRACERS = "INSUFFICIENT_TRACERS"
CODE_CONDITION_NUMBER = "EXCESSIVE_CONDITION_NUMBER"
CODE_TRACER_CONFLICT = "PHYSICAL_TRACER_CONFLICT"
CODE_GRAPH_CYCLES = "GRAPH_CONTAINS_CYCLES"
CODE_HEAD_CONFLICT = "HYDRAULIC_HEAD_CONFLICT"
CODE_POOR_FIT = "POOR_CALIBRATION_FIT"
CODE_UNIDENTIFIABLE = "AGE_UNIDENTIFIABLE"
CODE_MISSING_FORWARD_SYSTEM = "MISSING_FORWARD_SYSTEM"
CODE_MISSING_MIXING_SPEC = "MISSING_MIXING_SPEC"
CODE_MISSING_TRANSPORT_OPERATOR = "MISSING_TRANSPORT_OPERATOR"
CODE_OPTIMIZATION_FAILURE = "OPTIMIZATION_FAILED"
CODE_NODE_ID_COLLISION = "NODE_ID_COLLISION"


@dataclass(frozen=True)
class DiagnosticGateReport:
    """Report on pre-inversion and post-inversion physical consistency gates."""

    can_proceed: bool
    status: str  # "PROCEED" or "ABSTAIN"
    reason_codes: Tuple[str, ...]
    messages: Tuple[str, ...]
    metrics: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "reason_codes", tuple(self.reason_codes))
        object.__setattr__(self, "messages", tuple(self.messages))
        object.__setattr__(self, "metrics", dict(self.metrics))


def audit_node_physical_evidence(
    system: MultiTracerForwardSystem,
    *,
    max_condition_number: float = 1e8,
    min_required_tracers: int = 2,
    min_effective_rank: int = 2,
) -> DiagnosticGateReport:
    """Audit pre-inversion evidence at a single node for physical conflicts and conditioning."""
    reasons: list[str] = []
    messages: list[str] = []
    metrics: dict[str, float] = {}

    n_tracers = system.n_tracers
    metrics["n_tracers"] = float(n_tracers)

    if n_tracers == 0:
        return DiagnosticGateReport(
            can_proceed=False,
            status="ABSTAIN",
            reason_codes=(CODE_NO_TRACERS,),
            messages=(f"Node {system.node_id} has zero observed tracers.",),
            metrics=metrics,
        )

    if n_tracers < min_required_tracers:
        reasons.append(CODE_INSUFFICIENT_TRACERS)
        messages.append(f"Node {system.node_id} has {n_tracers} tracers, requires {min_required_tracers}.")

    # Matrix conditioning
    cond = system.condition_number()
    metrics["condition_number"] = cond
    if math.isinf(cond) or cond > max_condition_number:
        reasons.append(CODE_CONDITION_NUMBER)
        messages.append(f"Condition number ({cond:.2e}) exceeds threshold ({max_condition_number:.2e}).")

    effective_rank = system.effective_rank()
    metrics["effective_rank"] = effective_rank
    if effective_rank < float(min_effective_rank):
        reasons.append(CODE_UNIDENTIFIABLE)
        messages.append(
            f"Forward system effective rank ({effective_rank:.0f}) is below the "
            f"minimum required rank ({min_effective_rank})."
        )

    # Physical tracer contradiction screening
    # Map tracer values for convenient lookup
    obs_map = {tracer: val for tracer, val in zip(system.tracers, system.observations)}

    has_tritium = "3H" in obs_map
    has_sf6 = "SF6" in obs_map
    has_c14 = "14C" in obs_map

    # Contradiction rule: Modern SF6 (> 3.0 pptv) co-occurring with zero tritium (< 0.2 TU)
    # and depleted radiocarbon (< 10 pmc) indicates serious contamination or instrument failure
    if has_sf6 and has_tritium and has_c14:
        sf6_val = obs_map["SF6"]
        h3_val = obs_map["3H"]
        c14_val = obs_map["14C"]
        if sf6_val > 3.0 and h3_val < 0.2 and c14_val < 10.0:
            reasons.append(CODE_TRACER_CONFLICT)
            messages.append(
                f"Contradictory evidence: high modern SF6 ({sf6_val:.2f} pptv) "
                f"with sub-detection 3H ({h3_val:.2f} TU) and dead 14C ({c14_val:.1f} pmc)."
            )

    can_proceed = len(reasons) == 0
    status = "PROCEED" if can_proceed else "ABSTAIN"

    return DiagnosticGateReport(
        can_proceed=can_proceed,
        status=status,
        reason_codes=tuple(reasons),
        messages=tuple(messages),
        metrics=metrics,
    )


def audit_graph_topology(
    graph: nx.DiGraph,
    node_heads: Optional[Mapping[str, float]] = None,
    *,
    head_tolerance_m: float = 1.0,
) -> DiagnosticGateReport:
    """Audit graph DAG acyclicity and hydraulic head gradient consistency."""
    reasons: list[str] = []
    messages: list[str] = []
    metrics: dict[str, float] = {}

    metrics["n_nodes"] = float(graph.number_of_nodes())
    metrics["n_edges"] = float(graph.number_of_edges())

    if not nx.is_directed_acyclic_graph(graph):
        reasons.append(CODE_GRAPH_CYCLES)
        messages.append("Candidate flow graph contains directed cycles.")

    # Head consistency check: along edge (u, v), head_u should be >= head_v - tolerance
    if node_heads:
        adverse_edges = 0
        max_adverse_head = 0.0
        for u, v in graph.edges():
            h_u = node_heads.get(str(u))
            h_v = node_heads.get(str(v))
            if h_u is not None and h_v is not None:
                head_drop = h_u - h_v
                if head_drop < -head_tolerance_m:
                    adverse_edges += 1
                    max_adverse_head = max(max_adverse_head, abs(head_drop))

        metrics["adverse_head_edges"] = float(adverse_edges)
        metrics["max_adverse_head_m"] = float(max_adverse_head)

        if adverse_edges > 0:
            reasons.append(CODE_HEAD_CONFLICT)
            messages.append(
                f"{adverse_edges} edges flow against hydraulic head gradient "
                f"(max adverse drop = {max_adverse_head:.2f} m)."
            )

    can_proceed = len(reasons) == 0
    status = "PROCEED" if can_proceed else "ABSTAIN"

    return DiagnosticGateReport(
        can_proceed=can_proceed,
        status=status,
        reason_codes=tuple(reasons),
        messages=tuple(messages),
        metrics=metrics,
    )
