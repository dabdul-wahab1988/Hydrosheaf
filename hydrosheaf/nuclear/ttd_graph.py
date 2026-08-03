"""Graph compatibility audits for local TTD identified sets.

Version 1 deliberately does not condition or tighten local evidence.  It asks
whether two independently constructed local feasible sets can satisfy a
declared non-negative, mass-conserving transport map, and returns the smallest
L1 obstruction when they cannot.  This supports graph falsification without
turning an unvalidated edge into additional age information.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.optimize import linprog

from .ttd_identified import (
    TtdEvidenceReport,
    _audit_witness,
    _compile_constraints,
)


@dataclass(frozen=True)
class MassTransportMap:
    """Declared source-to-target probability-mass transport operator."""

    edge_id: str
    matrix: Sequence[Sequence[float]]
    provenance_tier: str
    source: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        matrix = np.asarray(self.matrix, dtype=float).copy()
        if not str(self.edge_id).strip() or not str(self.source).strip():
            raise ValueError("edge_id and source must be non-empty.")
        if matrix.ndim != 2 or min(matrix.shape) == 0:
            raise ValueError("transport matrix must be two-dimensional and non-empty.")
        if not np.all(np.isfinite(matrix)) or np.any(matrix < -1.0e-12):
            raise ValueError("transport matrix must be finite and non-negative.")
        if not np.allclose(matrix.sum(axis=0), 1.0, rtol=0.0, atol=1.0e-8):
            raise ValueError(
                "transport matrix columns must sum to one to conserve source mass."
            )
        matrix.setflags(write=False)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class GraphCompatibilityAudit:
    """Minimum graph obstruction between two frozen local evidence sets."""

    status: str
    edge_id: str
    minimum_l1_gap: float
    maximum_bin_gap: float
    source_witness: np.ndarray | None
    target_witness: np.ndarray | None
    mapped_source_mass: np.ndarray | None
    provenance_tier: str
    transport_source: str
    tightening_authorized: bool = False
    message: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self, *, include_witnesses: bool = True) -> dict[str, Any]:
        payload = {
            "status": self.status,
            "edge_id": self.edge_id,
            "minimum_l1_gap": float(self.minimum_l1_gap),
            "maximum_bin_gap": float(self.maximum_bin_gap),
            "provenance_tier": self.provenance_tier,
            "transport_source": self.transport_source,
            "tightening_authorized": False,
            "message": self.message,
            "metadata": dict(self.metadata),
        }
        if include_witnesses:
            payload.update(
                {
                    "source_witness": (
                        self.source_witness.tolist()
                        if self.source_witness is not None
                        else None
                    ),
                    "target_witness": (
                        self.target_witness.tolist()
                        if self.target_witness is not None
                        else None
                    ),
                    "mapped_source_mass": (
                        self.mapped_source_mass.tolist()
                        if self.mapped_source_mass is not None
                        else None
                    ),
                }
            )
        return payload


def _padded_constraint_rows(
    report: TtdEvidenceReport,
    total_variables: int,
    offset: int,
) -> tuple[list[np.ndarray], list[float]]:
    matrix, limits = _compile_constraints(
        report.constraints,
        report.age_grid_years.size,
        report.sigma_multiplier,
    )
    if matrix is None or limits is None:
        return [], []
    rows: list[np.ndarray] = []
    for row in matrix:
        padded = np.zeros(total_variables, dtype=float)
        padded[offset : offset + row.size] = row
        rows.append(padded)
    return rows, limits.tolist()


def audit_ttd_graph_compatibility(
    source_report: TtdEvidenceReport,
    target_report: TtdEvidenceReport,
    transport: MassTransportMap,
    *,
    compatibility_tolerance: float = 1.0e-7,
) -> GraphCompatibilityAudit:
    """Minimize ``||M w_source - w_target||_1`` over both local sets."""
    if not math.isfinite(float(compatibility_tolerance)) or compatibility_tolerance < 0:
        raise ValueError("compatibility_tolerance must be finite and non-negative.")
    if source_report.feasible_witness is None or target_report.feasible_witness is None:
        return GraphCompatibilityAudit(
            status="ABSTAIN",
            edge_id=transport.edge_id,
            minimum_l1_gap=float("nan"),
            maximum_bin_gap=float("nan"),
            source_witness=None,
            target_witness=None,
            mapped_source_mass=None,
            provenance_tier=transport.provenance_tier,
            transport_source=transport.source,
            message="At least one local identified set is infeasible or unsupported.",
            metadata=transport.metadata,
        )

    n_source = source_report.age_grid_years.size
    n_target = target_report.age_grid_years.size
    matrix = np.asarray(transport.matrix, dtype=float)
    if matrix.shape != (n_target, n_source):
        raise ValueError(
            "transport matrix shape must be (target age bins, source age bins)."
        )
    total = n_source + n_target + n_target
    source_slice = slice(0, n_source)
    target_slice = slice(n_source, n_source + n_target)
    gap_slice = slice(n_source + n_target, total)

    rows, limits = _padded_constraint_rows(source_report, total, 0)
    target_rows, target_limits = _padded_constraint_rows(target_report, total, n_source)
    rows.extend(target_rows)
    limits.extend(target_limits)
    identity = np.eye(n_target, dtype=float)
    positive = np.zeros((n_target, total), dtype=float)
    positive[:, source_slice] = matrix
    positive[:, target_slice] = -identity
    positive[:, gap_slice] = -identity
    negative = np.zeros((n_target, total), dtype=float)
    negative[:, source_slice] = -matrix
    negative[:, target_slice] = identity
    negative[:, gap_slice] = -identity
    rows.extend(positive)
    rows.extend(negative)
    limits.extend([0.0] * (2 * n_target))

    equalities = np.zeros((2, total), dtype=float)
    equalities[0, source_slice] = 1.0
    equalities[1, target_slice] = 1.0
    objective = np.zeros(total, dtype=float)
    objective[gap_slice] = 1.0
    result = linprog(
        objective,
        A_ub=np.vstack(rows) if rows else None,
        b_ub=np.asarray(limits, dtype=float) if limits else None,
        A_eq=equalities,
        b_eq=np.ones(2, dtype=float),
        bounds=(0.0, None),
        method="highs",
    )
    if not result.success:
        return GraphCompatibilityAudit(
            status="ABSTAIN",
            edge_id=transport.edge_id,
            minimum_l1_gap=float("nan"),
            maximum_bin_gap=float("nan"),
            source_witness=None,
            target_witness=None,
            mapped_source_mass=None,
            provenance_tier=transport.provenance_tier,
            transport_source=transport.source,
            message=f"Joint compatibility optimization failed: {result.message}",
            metadata=transport.metadata,
        )

    source_witness = np.asarray(result.x[source_slice], dtype=float)
    target_witness = np.asarray(result.x[target_slice], dtype=float)
    source_a, source_b = _compile_constraints(
        source_report.constraints,
        n_source,
        source_report.sigma_multiplier,
    )
    target_a, target_b = _compile_constraints(
        target_report.constraints,
        n_target,
        target_report.sigma_multiplier,
    )
    _audit_witness(
        source_witness,
        source_a,
        source_b,
        source_report.feasibility_tolerance,
    )
    _audit_witness(
        target_witness,
        target_a,
        target_b,
        target_report.feasibility_tolerance,
    )
    mapped = matrix @ source_witness
    differences = np.abs(mapped - target_witness)
    minimum_gap = float(differences.sum())
    maximum_gap = float(differences.max(initial=0.0))
    compatible = minimum_gap <= float(compatibility_tolerance)
    return GraphCompatibilityAudit(
        status="COMPATIBLE" if compatible else "INCOMPATIBLE",
        edge_id=transport.edge_id,
        minimum_l1_gap=minimum_gap,
        maximum_bin_gap=maximum_gap,
        source_witness=source_witness,
        target_witness=target_witness,
        mapped_source_mass=mapped,
        provenance_tier=transport.provenance_tier,
        transport_source=transport.source,
        tightening_authorized=False,
        message=(
            "The declared map intersects both frozen local identified sets."
            if compatible
            else "No exact mapped pair exists; the reported gap is a graph obstruction."
        ),
        metadata={
            **transport.metadata,
            "compatibility_tolerance": float(compatibility_tolerance),
            "conditioning_performed": False,
        },
    )


__all__ = [
    "GraphCompatibilityAudit",
    "MassTransportMap",
    "audit_ttd_graph_compatibility",
]
