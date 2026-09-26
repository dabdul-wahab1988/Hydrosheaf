"""Bounded polyhedral state spaces for the robust ACMD loop.

The TTD design solver uses a probability simplex.  Chemical reaction extents
may be signed and have independent physical bounds, so they require a different
feasible set while retaining ACMD's paired-state ambiguity criterion.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import linprog

from .nuclear.ttd_certified_design import AmbiguityEvaluation, CertifiedCandidateTracer
from .nuclear.ttd_identified import AgeFunctional, TracerConstraint, _compile_constraints


def _matrix(values: Optional[Sequence[Sequence[float]]], name: str, width: int) -> Optional[np.ndarray]:
    if values is None:
        return None
    matrix = np.asarray(values, dtype=float).copy()
    if matrix.ndim != 2 or matrix.shape[1] != width:
        raise ValueError(f"{name} must have {width} columns.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must contain only finite values.")
    if matrix.shape[0] == 0:
        return None
    matrix.setflags(write=False)
    return matrix


def _vector(values: Optional[Sequence[float]], name: str, length: int) -> Optional[np.ndarray]:
    if values is None:
        return None
    vector = np.asarray(values, dtype=float).copy()
    if vector.shape != (length,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must be a finite vector of length {length}.")
    vector.setflags(write=False)
    return vector


@dataclass(frozen=True)
class PolyhedralStateSpace:
    """A finite, auditable LP domain for a non-TTD ACMD state vector."""

    state_names: Sequence[str]
    bounds: Sequence[tuple[float, float]]
    a_ub: Optional[Sequence[Sequence[float]]] = None
    b_ub: Optional[Sequence[float]] = None
    a_eq: Optional[Sequence[Sequence[float]]] = None
    b_eq: Optional[Sequence[float]] = None

    def __post_init__(self) -> None:
        names = tuple(str(name) for name in self.state_names)
        if not names or any(not name.strip() for name in names) or len(set(names)) != len(names):
            raise ValueError("state_names must be non-empty and unique.")
        raw_bounds = np.asarray(self.bounds, dtype=float)
        if raw_bounds.shape != (len(names), 2) or not np.all(np.isfinite(raw_bounds)):
            raise ValueError("bounds must contain one finite (lower, upper) pair per state.")
        if np.any(raw_bounds[:, 0] > raw_bounds[:, 1]):
            raise ValueError("Each state lower bound must not exceed its upper bound.")
        a_ub = _matrix(self.a_ub, "a_ub", len(names))
        a_eq = _matrix(self.a_eq, "a_eq", len(names))
        b_ub = _vector(self.b_ub, "b_ub", 0 if a_ub is None else a_ub.shape[0])
        b_eq = _vector(self.b_eq, "b_eq", 0 if a_eq is None else a_eq.shape[0])
        if (a_ub is None) != (b_ub is None) or (a_eq is None) != (b_eq is None):
            raise ValueError("Each constraint matrix must have a matching right-hand side.")
        raw_bounds.setflags(write=False)
        object.__setattr__(self, "state_names", names)
        object.__setattr__(self, "bounds", raw_bounds)
        object.__setattr__(self, "a_ub", a_ub)
        object.__setattr__(self, "b_ub", b_ub)
        object.__setattr__(self, "a_eq", a_eq)
        object.__setattr__(self, "b_eq", b_eq)

    @property
    def dimension(self) -> int:
        return len(self.state_names)

    def to_dict(self) -> dict[str, object]:
        """Canonical numerical specification for problem provenance."""
        return {
            "state_names": list(self.state_names),
            "bounds": self.bounds.tolist(),
            "a_ub": None if self.a_ub is None else self.a_ub.tolist(),
            "b_ub": None if self.b_ub is None else self.b_ub.tolist(),
            "a_eq": None if self.a_eq is None else self.a_eq.tolist(),
            "b_eq": None if self.b_eq is None else self.b_eq.tolist(),
        }


def _audit_state(
    state: np.ndarray,
    domain: PolyhedralStateSpace,
    a_ub: Optional[np.ndarray],
    b_ub: Optional[np.ndarray],
    tolerance: float,
) -> None:
    if not np.all(np.isfinite(state)):
        raise RuntimeError("LP witness contains non-finite values.")
    if np.any(state < domain.bounds[:, 0] - tolerance) or np.any(state > domain.bounds[:, 1] + tolerance):
        raise RuntimeError("LP witness violates state bounds.")
    if a_ub is not None and b_ub is not None and np.any(a_ub @ state > b_ub + tolerance):
        raise RuntimeError("LP witness violates inequality constraints.")
    if domain.a_eq is not None and domain.b_eq is not None:
        if np.any(np.abs(domain.a_eq @ state - domain.b_eq) > tolerance):
            raise RuntimeError("LP witness violates equality constraints.")


def evaluate_polyhedral_ambiguity(
    state_space: PolyhedralStateSpace,
    prior_constraints: Sequence[TracerConstraint],
    candidate_subset: Sequence[CertifiedCandidateTracer],
    target_functional: AgeFunctional,
    *,
    sigma_multiplier: float = 1.96,
    feasibility_tolerance: float = 1.0e-7,
) -> AmbiguityEvaluation:
    """Maximize q·(z-z') with z,z' feasible and candidates indistinguishable.

    Candidate error bounds are symmetric absolute errors, hence two states that
    could produce the same future observation may differ by at most 2*epsilon.
    """
    n = state_space.dimension
    q = np.asarray(target_functional.coefficients, dtype=float)
    if q.shape != (n,):
        raise ValueError("target_functional length does not match state space.")
    if not np.isfinite(sigma_multiplier) or sigma_multiplier <= 0:
        raise ValueError("sigma_multiplier must be finite and positive.")
    if not np.isfinite(feasibility_tolerance) or feasibility_tolerance <= 0:
        raise ValueError("feasibility_tolerance must be finite and positive.")

    obs_a, obs_b = _compile_constraints(tuple(prior_constraints), n, sigma_multiplier)
    pieces_a = [a for a in (state_space.a_ub, obs_a) if a is not None]
    pieces_b = [b for b in (state_space.b_ub, obs_b) if b is not None]
    single_a = np.vstack(pieces_a) if pieces_a else None
    single_b = np.concatenate(pieces_b) if pieces_b else None
    lp_bounds = [tuple(pair) for pair in state_space.bounds]
    subset_ids = tuple(candidate.option_id for candidate in candidate_subset)

    feasible = linprog(
        np.zeros(n), A_ub=single_a, b_ub=single_b,
        A_eq=state_space.a_eq, b_eq=state_space.b_eq,
        bounds=lp_bounds, method="highs",
    )
    if not feasible.success:
        status = "INFEASIBLE_PRIOR" if feasible.status == 2 else f"NUMERICAL_ERROR_{feasible.status}"
        return AmbiguityEvaluation(subset_ids, float("nan"), status, tolerance=feasibility_tolerance)

    joint_a: list[np.ndarray] = []
    joint_b: list[np.ndarray] = []
    if single_a is not None and single_b is not None:
        zero = np.zeros_like(single_a)
        joint_a.extend((np.hstack((single_a, zero)), np.hstack((zero, single_a))))
        joint_b.extend((single_b, single_b))

    # The coordinates are only a shape carrier here; every chemical action must
    # provide its explicit linear response and never calls an age-tracer kernel.
    coordinates = np.arange(n, dtype=float)
    resolved: list[tuple[np.ndarray, float]] = []
    for candidate in candidate_subset:
        if candidate.response is None:
            raise ValueError("Polyhedral ACMD candidates require an explicit response.")
        h = candidate.resolve_response(coordinates)
        limit = 2.0 * float(candidate.error_bound)
        resolved.append((h, limit))
        joint_a.extend((np.r_[h, -h][None, :], np.r_[-h, h][None, :]))
        joint_b.extend((np.asarray([limit]), np.asarray([limit])))

    eq_a = None
    eq_b = None
    if state_space.a_eq is not None and state_space.b_eq is not None:
        zero = np.zeros_like(state_space.a_eq)
        eq_a = np.vstack((np.hstack((state_space.a_eq, zero)), np.hstack((zero, state_space.a_eq))))
        eq_b = np.r_[state_space.b_eq, state_space.b_eq]

    result = linprog(
        np.r_[-q, q],
        A_ub=np.vstack(joint_a) if joint_a else None,
        b_ub=np.concatenate(joint_b) if joint_b else None,
        A_eq=eq_a, b_eq=eq_b,
        bounds=lp_bounds * 2,
        method="highs",
    )
    if not result.success:
        status = "INFEASIBLE_JOINT" if result.status == 2 else f"NUMERICAL_ERROR_{result.status}"
        return AmbiguityEvaluation(subset_ids, float("nan"), status, tolerance=feasibility_tolerance)

    upper = np.asarray(result.x[:n], dtype=float)
    lower = np.asarray(result.x[n:], dtype=float)
    _audit_state(upper, state_space, single_a, single_b, feasibility_tolerance)
    _audit_state(lower, state_space, single_a, single_b, feasibility_tolerance)
    for h, limit in resolved:
        if abs(float(h @ (upper - lower))) > limit + feasibility_tolerance:
            raise RuntimeError("LP witness violates candidate indistinguishability.")
    width = float(max(0.0, q @ (upper - lower)))
    return AmbiguityEvaluation(
        subset_ids, width, "FEASIBLE", lower_witness=lower,
        upper_witness=upper, tolerance=feasibility_tolerance,
    )
