"""Certified measurement design for partially identified groundwater systems.

Provides exact minimax guarantees for transit-time-distribution (TTD) experimental
design over the convex polyhedral compatible set P defined by prior tracer data.
Unlike Bayesian optimal experimental design, this module requires NO prior probability
distribution over hypotheses or age bins, refusing to invent probabilities where
none are defensible.

Key concepts:
  - Feasible set P: Polytope of age-bin mass vectors compatible with prior observations.
  - Decision functional q^T x: Linear scalar target of scientific or policy interest
    (e.g., Anthropocene young-water fraction, Holocene fraction).
  - Worst-case ambiguity W(S): The maximum possible difference |q^T (x - x')| between
    any two states x, x' in P that future candidate measurements S cannot distinguish.
  - Sufficiency Certificate: Guarantees that measuring subset S will resolve target q^T x
    within specified tolerance delta, regardless of the true state in P.
  - Impossibility Witness: Explicit pair of indistinguishable states (x, x') proving that
    even measuring every candidate tracer cannot achieve tolerance delta.
  - Inconsistent Evidence: Detected when P is empty (data contradiction or contamination).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import itertools
import json
import math
from typing import Any, Callable, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog

from .joint_lpm import tracer_response_kernel
from .ttd_identified import (
    AgeFunctional,
    TracerConstraint,
    _compile_constraints,
    _readonly_vector,
)


@dataclass(frozen=True)
class CertifiedCandidateTracer:
    """One candidate future tracer measurement available for selection."""

    option_id: str
    tracer: str
    sample_year: float
    error_bound: float
    cost: float = 1.0
    response: Optional[Sequence[float]] = None
    kernel_kwargs: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not str(self.option_id).strip():
            raise ValueError("option_id must be non-empty.")
        if not str(self.tracer).strip():
            raise ValueError("tracer must be non-empty.")
        if not math.isfinite(float(self.sample_year)):
            raise ValueError("sample_year must be finite.")
        if not math.isfinite(float(self.error_bound)) or float(self.error_bound) <= 0.0:
            raise ValueError("error_bound must be finite and strictly positive.")
        if not math.isfinite(float(self.cost)) or float(self.cost) <= 0.0:
            raise ValueError("cost must be finite and strictly positive.")
        if self.response is not None:
            object.__setattr__(
                self, "response", _readonly_vector(self.response, name="response")
            )
        object.__setattr__(self, "kernel_kwargs", dict(self.kernel_kwargs))
        object.__setattr__(self, "metadata", dict(self.metadata))

    def resolve_response(self, age_grid_years: np.ndarray) -> np.ndarray:
        """Return the response vector h_a on the given age grid."""
        if self.response is not None:
            resp = np.asarray(self.response, dtype=float)
            if resp.shape != age_grid_years.shape:
                raise ValueError(
                    f"Candidate {self.option_id!r} response length {resp.size} "
                    f"does not match age grid size {age_grid_years.size}."
                )
            return resp
        kernel = tracer_response_kernel(
            self.tracer,
            age_grid_years,
            float(self.sample_year),
            **self.kernel_kwargs,
        )
        return np.asarray(kernel, dtype=float)


@dataclass(frozen=True)
class AmbiguityEvaluation:
    """Result of evaluating worst-case target ambiguity for a measurement subset."""

    subset_ids: tuple[str, ...]
    worst_case_ambiguity: float
    status: str
    lower_witness: Optional[np.ndarray] = None
    upper_witness: Optional[np.ndarray] = None
    tolerance: float = 1.0e-7

    def __post_init__(self) -> None:
        object.__setattr__(self, "subset_ids", tuple(self.subset_ids))
        if self.lower_witness is not None:
            object.__setattr__(
                self,
                "lower_witness",
                _readonly_vector(self.lower_witness, name="lower_witness"),
            )
        if self.upper_witness is not None:
            object.__setattr__(
                self,
                "upper_witness",
                _readonly_vector(self.upper_witness, name="upper_witness"),
            )

    def to_dict(self, *, include_witnesses: bool = True) -> dict[str, Any]:
        data: dict[str, Any] = {
            "subset_ids": list(self.subset_ids),
            "worst_case_ambiguity": float(self.worst_case_ambiguity),
            "status": self.status,
            "tolerance": float(self.tolerance),
        }
        if include_witnesses and self.lower_witness is not None and self.upper_witness is not None:
            data["lower_witness"] = self.lower_witness.tolist()
            data["upper_witness"] = self.upper_witness.tolist()
        return data


@dataclass(frozen=True)
class CertifiedDesignCertificate:
    """Auditable certificate of measurement sufficiency, impossibility, or inconsistency."""

    status: str
    target_functional_name: str
    target_tolerance: float
    initial_ambiguity: float
    achieved_ambiguity: float
    all_candidates_ambiguity: float
    selected_option_ids: tuple[str, ...]
    total_cost: float
    candidate_count: int
    feasibility_status: str
    lower_witness: Optional[np.ndarray] = None
    upper_witness: Optional[np.ndarray] = None
    evaluations_count: int = 0
    certificate_hash: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "selected_option_ids", tuple(self.selected_option_ids))
        if self.lower_witness is not None:
            object.__setattr__(
                self,
                "lower_witness",
                _readonly_vector(self.lower_witness, name="lower_witness"),
            )
        if self.upper_witness is not None:
            object.__setattr__(
                self,
                "upper_witness",
                _readonly_vector(self.upper_witness, name="upper_witness"),
            )
        object.__setattr__(self, "metadata", dict(self.metadata))
        if not self.certificate_hash:
            payload = json.dumps(
                {
                    "status": self.status,
                    "target": self.target_functional_name,
                    "target_tolerance": float(self.target_tolerance),
                    "initial_ambiguity": float(self.initial_ambiguity),
                    "achieved_ambiguity": float(self.achieved_ambiguity),
                    "selected_option_ids": list(self.selected_option_ids),
                    "total_cost": float(self.total_cost),
                },
                sort_keys=True,
            )
            h = hashlib.sha256(payload.encode("utf-8")).hexdigest()
            object.__setattr__(self, "certificate_hash", h)

    def to_dict(self, *, include_witnesses: bool = True) -> dict[str, Any]:
        data: dict[str, Any] = {
            "status": self.status,
            "target_functional_name": self.target_functional_name,
            "target_tolerance": float(self.target_tolerance),
            "initial_ambiguity": float(self.initial_ambiguity),
            "achieved_ambiguity": float(self.achieved_ambiguity),
            "all_candidates_ambiguity": float(self.all_candidates_ambiguity),
            "selected_option_ids": list(self.selected_option_ids),
            "total_cost": float(self.total_cost),
            "candidate_count": int(self.candidate_count),
            "feasibility_status": self.feasibility_status,
            "evaluations_count": int(self.evaluations_count),
            "certificate_hash": self.certificate_hash,
            "metadata": dict(self.metadata),
        }
        if include_witnesses and self.lower_witness is not None and self.upper_witness is not None:
            data["lower_witness"] = self.lower_witness.tolist()
            data["upper_witness"] = self.upper_witness.tolist()
        return data


def _audit_witness_pair(
    x: np.ndarray,
    x_prime: np.ndarray,
    a_ub: Optional[np.ndarray],
    b_ub: Optional[np.ndarray],
    candidate_responses: Mapping[str, tuple[np.ndarray, float]],
    tolerance: float,
) -> None:
    """Verify that (x, x') are both in P and indistinguishable under candidate measurements."""
    for name, w in (("x", x), ("x_prime", x_prime)):
        if not np.all(np.isfinite(w)):
            raise RuntimeError(f"Linear-program witness {name} contains non-finite values.")
        if float(np.min(w)) < -tolerance:
            raise RuntimeError(f"Linear-program witness {name} violates non-negativity: {float(np.min(w))}.")
        if abs(float(np.sum(w)) - 1.0) > tolerance:
            raise RuntimeError(f"Linear-program witness {name} violates unit mass: {float(np.sum(w))}.")
        if a_ub is not None and b_ub is not None:
            violation = float(np.max(a_ub @ w - b_ub))
            if violation > tolerance:
                raise RuntimeError(
                    f"Linear-program witness {name} violates prior constraints by {violation:.3e}."
                )

    for opt_id, (resp, bound) in candidate_responses.items():
        diff = abs(float(resp @ (x - x_prime)))
        if diff > bound + tolerance:
            raise RuntimeError(
                f"Witness pair violates candidate {opt_id!r} error bound: {diff:.3e} > {bound:.3e}."
            )


def evaluate_worst_case_ambiguity(
    age_grid_years: Sequence[float],
    prior_constraints: Sequence[TracerConstraint],
    candidate_subset: Sequence[CertifiedCandidateTracer],
    target_functional: AgeFunctional,
    *,
    sigma_multiplier: float = 1.96,
    feasibility_tolerance: float = 1.0e-7,
) -> AmbiguityEvaluation:
    """Compute the worst-case target ambiguity W(S) over indistinguishable pairs in P.

    Formulates and solves the joint linear program over (x, x') in P x P:
        max  q^T (x - x')
        s.t. x, x' in P
             |h_a^T (x - x')| <= 2 * epsilon_a   for all a in candidate_subset

    Returns:
        AmbiguityEvaluation containing W(S) and the maximizing witness pair (x, x').
    """
    ages = _readonly_vector(age_grid_years, name="age_grid_years")
    n_bins = ages.size
    q = np.asarray(target_functional.coefficients, dtype=float)
    if q.shape != (n_bins,):
        raise ValueError(
            f"Target functional {target_functional.name!r} coefficients size {q.size} "
            f"does not match age grid size {n_bins}."
        )

    # 1. Compile prior polytope P constraints for a single x
    a_ub_single, b_ub_single = _compile_constraints(
        tuple(prior_constraints), n_bins, float(sigma_multiplier)
    )

    # Check feasibility of single x in P first
    test_res = linprog(
        np.zeros(n_bins),
        A_ub=a_ub_single,
        b_ub=b_ub_single,
        A_eq=np.ones((1, n_bins), dtype=float),
        b_eq=np.ones(1, dtype=float),
        bounds=(0.0, None),
        method="highs",
    )
    if not test_res.success:
        status_label = "INFEASIBLE_PRIOR" if test_res.status == 2 else f"NUMERICAL_ERROR_{test_res.status}"
        return AmbiguityEvaluation(
            subset_ids=tuple(c.option_id for c in candidate_subset),
            worst_case_ambiguity=float("nan"),
            status=status_label,
            tolerance=feasibility_tolerance,
        )

    # 2. Joint variable z = [x, x'] of size 2 * n_bins
    # Objective: maximize q^T x - q^T x'  ==> minimize [-q, q]^T z
    c_obj = np.concatenate([-q, q])

    # Equality constraints: sum(x) = 1, sum(x') = 1
    # [1 ... 1  0 ... 0]
    # [0 ... 0  1 ... 1]
    a_eq = np.zeros((2, 2 * n_bins), dtype=float)
    a_eq[0, :n_bins] = 1.0
    a_eq[1, n_bins:] = 1.0
    b_eq = np.ones(2, dtype=float)

    # Inequality constraints from P:
    # A_ub x <= b_ub  ==>  [A_ub, 0] z <= b_ub
    # A_ub x' <= b_ub ==>  [0, A_ub] z <= b_ub
    ub_rows: list[np.ndarray] = []
    ub_limits: list[float] = []

    if a_ub_single is not None and b_ub_single is not None:
        n_p_rows = a_ub_single.shape[0]
        # x constraints
        row_x = np.zeros((n_p_rows, 2 * n_bins), dtype=float)
        row_x[:, :n_bins] = a_ub_single
        ub_rows.append(row_x)
        ub_limits.extend(b_ub_single.tolist())

        # x' constraints
        row_xp = np.zeros((n_p_rows, 2 * n_bins), dtype=float)
        row_xp[:, n_bins:] = a_ub_single
        ub_rows.append(row_xp)
        ub_limits.extend(b_ub_single.tolist())

    # Inequality constraints from candidate measurements S:
    # |h_a^T (x - x')| <= 2 * epsilon_a
    # ==>  h_a^T x - h_a^T x' <= 2 * epsilon_a
    # ==> -h_a^T x + h_a^T x' <= 2 * epsilon_a
    resolved_candidates: dict[str, tuple[np.ndarray, float]] = {}
    for cand in candidate_subset:
        h_a = cand.resolve_response(ages)
        bound = 2.0 * float(cand.error_bound)
        resolved_candidates[cand.option_id] = (h_a, bound)

        r_pos = np.concatenate([h_a, -h_a])[None, :]
        r_neg = np.concatenate([-h_a, h_a])[None, :]

        ub_rows.append(r_pos)
        ub_limits.append(bound)

        ub_rows.append(r_neg)
        ub_limits.append(bound)

    a_ub = np.vstack(ub_rows) if ub_rows else None
    b_ub = np.asarray(ub_limits, dtype=float) if ub_limits else None

    res = linprog(
        c_obj,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=a_eq,
        b_eq=b_eq,
        bounds=(0.0, None),
        method="highs",
    )

    if not res.success:
        status_label = "INFEASIBLE_JOINT" if res.status == 2 else f"NUMERICAL_ERROR_{res.status}"
        return AmbiguityEvaluation(
            subset_ids=tuple(c.option_id for c in candidate_subset),
            worst_case_ambiguity=float("nan"),
            status=status_label,
            tolerance=feasibility_tolerance,
        )

    z_star = np.asarray(res.x, dtype=float)
    x_star = z_star[:n_bins]
    xp_star = z_star[n_bins:]

    # Due to numerical tolerance, clip tiny negative values
    x_star = np.clip(x_star, 0.0, None)
    if x_star.sum() > 0.0:
        x_star /= x_star.sum()

    xp_star = np.clip(xp_star, 0.0, None)
    if xp_star.sum() > 0.0:
        xp_star /= xp_star.sum()

    # The maximum ambiguity is q^T x - q^T x' == -res.fun
    ambiguity = float(max(0.0, q @ x_star - q @ xp_star))

    _audit_witness_pair(
        x_star,
        xp_star,
        a_ub_single,
        b_ub_single,
        resolved_candidates,
        feasibility_tolerance,
    )

    return AmbiguityEvaluation(
        subset_ids=tuple(c.option_id for c in candidate_subset),
        worst_case_ambiguity=ambiguity,
        status="FEASIBLE",
        lower_witness=xp_star,
        upper_witness=x_star,
        tolerance=feasibility_tolerance,
    )


def solve_certified_measurement_design(
    age_grid_years: Sequence[float],
    prior_constraints: Sequence[TracerConstraint],
    candidates: Sequence[CertifiedCandidateTracer],
    target_functional: AgeFunctional,
    target_tolerance: float,
    *,
    sigma_multiplier: float = 1.96,
    feasibility_tolerance: float = 1.0e-7,
    max_exhaustive_subset_size: int = 12,
    metadata: Optional[Mapping[str, Any]] = None,
) -> CertifiedDesignCertificate:
    """Find the minimum-cost candidate tracer subset S* achieving W(S*) <= target_tolerance.

    Provides three definitive mathematical outcomes:
      1. INCONSISTENT_EVIDENCE: The prior constraints P are mutually contradictory.
      2. ALREADY_RESOLVED: Prior constraints already resolve the target within target_tolerance.
      3. IMPOSSIBILITY_WITNESS: Even measuring all candidates cannot achieve target_tolerance.
      4. CERTIFIED_SUFFICIENT: Found minimal-cost subset S* certified to achieve W(S*) <= delta.
    """
    if not math.isfinite(float(target_tolerance)) or float(target_tolerance) <= 0.0:
        raise ValueError("target_tolerance must be finite and strictly positive.")

    cand_list = list(candidates)
    cand_ids = [c.option_id for c in cand_list]
    if len(set(cand_ids)) != len(cand_ids):
        raise ValueError("Candidate option_ids must be unique.")

    # 1. Evaluate W(empty): initial ambiguity without any candidate measurements
    eval_empty = evaluate_worst_case_ambiguity(
        age_grid_years,
        prior_constraints,
        (),
        target_functional,
        sigma_multiplier=sigma_multiplier,
        feasibility_tolerance=feasibility_tolerance,
    )

    if eval_empty.status != "FEASIBLE":
        status_label = "INCONSISTENT_EVIDENCE" if eval_empty.status == "INFEASIBLE_PRIOR" else "NUMERICAL_ERROR"
        return CertifiedDesignCertificate(
            status=status_label,
            target_functional_name=target_functional.name,
            target_tolerance=float(target_tolerance),
            initial_ambiguity=float("nan"),
            achieved_ambiguity=float("nan"),
            all_candidates_ambiguity=float("nan"),
            selected_option_ids=(),
            total_cost=0.0,
            candidate_count=len(cand_list),
            feasibility_status=eval_empty.status,
            evaluations_count=1,
            metadata=dict(metadata or {}),
        )

    initial_w = eval_empty.worst_case_ambiguity

    # If prior ambiguity already satisfies tolerance, no measurements are needed
    if initial_w <= float(target_tolerance) + feasibility_tolerance:
        return CertifiedDesignCertificate(
            status="ALREADY_RESOLVED",
            target_functional_name=target_functional.name,
            target_tolerance=float(target_tolerance),
            initial_ambiguity=initial_w,
            achieved_ambiguity=initial_w,
            all_candidates_ambiguity=initial_w,
            selected_option_ids=(),
            total_cost=0.0,
            candidate_count=len(cand_list),
            feasibility_status="FEASIBLE",
            lower_witness=eval_empty.lower_witness,
            upper_witness=eval_empty.upper_witness,
            evaluations_count=1,
            metadata=dict(metadata or {}),
        )

    # 2. Evaluate W(All): ambiguity if ALL candidate measurements are taken
    eval_all = evaluate_worst_case_ambiguity(
        age_grid_years,
        prior_constraints,
        cand_list,
        target_functional,
        sigma_multiplier=sigma_multiplier,
        feasibility_tolerance=feasibility_tolerance,
    )
    # A feasible prior does not guarantee that the all-candidate optimisation
    # completed.  In particular, an iteration limit or other numerical
    # failure must not be treated as an impossibility witness and must not
    # seed the subset search with NaN/non-finite values.
    if eval_all.status != "FEASIBLE":
        status_label = (
            "INCONSISTENT_EVIDENCE"
            if eval_all.status.startswith("INFEASIBLE")
            else "NUMERICAL_ERROR"
        )
        error_metadata = dict(metadata or {})
        error_metadata["all_candidates_evaluation_status"] = eval_all.status
        return CertifiedDesignCertificate(
            status=status_label,
            target_functional_name=target_functional.name,
            target_tolerance=float(target_tolerance),
            initial_ambiguity=initial_w,
            achieved_ambiguity=float("nan"),
            all_candidates_ambiguity=float("nan"),
            selected_option_ids=(),
            total_cost=0.0,
            candidate_count=len(cand_list),
            feasibility_status=eval_all.status,
            evaluations_count=2,
            metadata=error_metadata,
        )
    all_w = eval_all.worst_case_ambiguity

    # If even measuring everything fails to reach delta, emit IMPOSSIBILITY_WITNESS
    if all_w > float(target_tolerance) + feasibility_tolerance:
        return CertifiedDesignCertificate(
            status="IMPOSSIBILITY_WITNESS",
            target_functional_name=target_functional.name,
            target_tolerance=float(target_tolerance),
            initial_ambiguity=initial_w,
            achieved_ambiguity=all_w,
            all_candidates_ambiguity=all_w,
            selected_option_ids=tuple(cand_ids),
            total_cost=sum(c.cost for c in cand_list),
            candidate_count=len(cand_list),
            feasibility_status="FEASIBLE",
            lower_witness=eval_all.lower_witness,
            upper_witness=eval_all.upper_witness,
            evaluations_count=2,
            metadata=dict(metadata or {}),
        )

    # 3. Find minimum-cost subset S* such that W(S*) <= target_tolerance
    # For small sets (<= max_exhaustive_subset_size), exact branch-by-cost enumeration guarantees global optimality
    evaluations = 2
    best_subset: Optional[tuple[CertifiedCandidateTracer, ...]] = tuple(cand_list)
    best_cost = sum(c.cost for c in cand_list)
    best_ambiguity = all_w
    best_witnesses = (eval_all.lower_witness, eval_all.upper_witness)

    n_cands = len(cand_list)
    if n_cands <= max_exhaustive_subset_size:
        # Build all non-empty subsets sorted by cost ascending
        all_subsets = []
        for r in range(1, n_cands):
            for combo in itertools.combinations(cand_list, r):
                cost = sum(c.cost for c in combo)
                all_subsets.append((cost, combo))
        all_subsets.sort(key=lambda item: (item[0], len(item[1])))

        for cost, subset in all_subsets:
            # Prune: if this subset already costs as much or more than our best known feasible solution, skip
            if cost >= best_cost:
                continue

            ev = evaluate_worst_case_ambiguity(
                age_grid_years,
                prior_constraints,
                subset,
                target_functional,
                sigma_multiplier=sigma_multiplier,
                feasibility_tolerance=feasibility_tolerance,
            )
            evaluations += 1

            if ev.status == "FEASIBLE" and ev.worst_case_ambiguity <= float(target_tolerance) + feasibility_tolerance:
                best_cost = cost
                best_subset = subset
                best_ambiguity = ev.worst_case_ambiguity
                best_witnesses = (ev.lower_witness, ev.upper_witness)
    else:
        # For larger candidate suites: greedy forward selection baseline with cost-awareness
        current_subset: list[CertifiedCandidateTracer] = []
        remaining = list(cand_list)
        while remaining:
            best_cand = None
            best_reduction_per_cost = -1.0
            best_ev = None

            for cand in remaining:
                trial = [*current_subset, cand]
                ev = evaluate_worst_case_ambiguity(
                    age_grid_years,
                    prior_constraints,
                    trial,
                    target_functional,
                    sigma_multiplier=sigma_multiplier,
                    feasibility_tolerance=feasibility_tolerance,
                )
                evaluations += 1
                curr_w = (
                    best_ambiguity
                    if current_subset
                    else initial_w
                )
                reduction = max(0.0, curr_w - ev.worst_case_ambiguity)
                ratio = reduction / cand.cost
                if ratio > best_reduction_per_cost:
                    best_reduction_per_cost = ratio
                    best_cand = cand
                    best_ev = ev

            if best_cand is None or best_ev is None:
                break

            current_subset.append(best_cand)
            remaining.remove(best_cand)
            best_ambiguity = best_ev.worst_case_ambiguity

            if best_ambiguity <= float(target_tolerance) + feasibility_tolerance:
                best_subset = tuple(current_subset)
                best_cost = sum(c.cost for c in best_subset)
                best_witnesses = (best_ev.lower_witness, best_ev.upper_witness)
                break

    return CertifiedDesignCertificate(
        status="CERTIFIED_SUFFICIENT",
        target_functional_name=target_functional.name,
        target_tolerance=float(target_tolerance),
        initial_ambiguity=initial_w,
        achieved_ambiguity=best_ambiguity,
        all_candidates_ambiguity=all_w,
        selected_option_ids=tuple(c.option_id for c in (best_subset or ())),
        total_cost=best_cost,
        candidate_count=len(cand_list),
        feasibility_status="FEASIBLE",
        lower_witness=best_witnesses[0],
        upper_witness=best_witnesses[1],
        evaluations_count=evaluations,
        metadata=dict(metadata or {}),
    )


def solve_budgeted_minimax_design(
    age_grid_years: Sequence[float],
    prior_constraints: Sequence[TracerConstraint],
    candidates: Sequence[CertifiedCandidateTracer],
    target_functional: AgeFunctional,
    budget: float,
    *,
    target_tolerance: Optional[float] = None,
    sigma_multiplier: float = 1.96,
    feasibility_tolerance: float = 1.0e-7,
    metadata: Optional[Mapping[str, Any]] = None,
    cost_function: Optional[Callable[[tuple[CertifiedCandidateTracer, ...]], float]] = None,
) -> CertifiedDesignCertificate:
    """Find the affordable subset that minimizes worst-case ambiguity ``W(S)``.

    By default, a subset costs the sum of its candidate ``cost`` values. A
    ``cost_function`` may be supplied when a field campaign has shared or
    fixed charges (for example, one purging visit, one courier shipment or a
    batch customs fee). It receives the selected candidates as a tuple and
    must return a finite, non-negative amount in the same currency as
    ``budget``. This keeps the minimax calculation exact while preventing
    shared campaign costs from being double-counted across tracers.
    """
    if not math.isfinite(float(budget)) or float(budget) < 0.0:
        raise ValueError("budget must be finite and non-negative.")

    cand_list = list(candidates)
    cand_ids = [c.option_id for c in cand_list]
    if len(set(cand_ids)) != len(cand_ids):
        raise ValueError("Candidate option_ids must be unique.")

    def subset_cost(subset: Sequence[CertifiedCandidateTracer]) -> float:
        subset_tuple = tuple(subset)
        raw_cost = (
            sum(c.cost for c in subset_tuple)
            if cost_function is None
            else cost_function(subset_tuple)
        )
        try:
            cost_value = float(raw_cost)
        except (TypeError, ValueError) as exc:
            raise ValueError("cost_function must return a numeric cost") from exc
        if not math.isfinite(cost_value) or cost_value < 0.0:
            raise ValueError("cost_function must return a finite, non-negative cost")
        return cost_value

    target_tol = (
        float(target_tolerance)
        if target_tolerance is not None
        else (
            float(target_functional.maximum_reportable_width)
            if target_functional.maximum_reportable_width is not None
            else None
        )
    )

    eval_empty = evaluate_worst_case_ambiguity(
        age_grid_years,
        prior_constraints,
        (),
        target_functional,
        sigma_multiplier=sigma_multiplier,
        feasibility_tolerance=feasibility_tolerance,
    )

    if eval_empty.status != "FEASIBLE":
        status_label = "INCONSISTENT_EVIDENCE" if eval_empty.status == "INFEASIBLE_PRIOR" else "NUMERICAL_ERROR"
        return CertifiedDesignCertificate(
            status=status_label,
            target_functional_name=target_functional.name,
            target_tolerance=target_tol if target_tol is not None else float("nan"),
            initial_ambiguity=float("nan"),
            achieved_ambiguity=float("nan"),
            all_candidates_ambiguity=float("nan"),
            selected_option_ids=(),
            total_cost=0.0,
            candidate_count=len(cand_list),
            feasibility_status=eval_empty.status,
            evaluations_count=1,
            metadata=dict(metadata or {}),
        )

    initial_w = eval_empty.worst_case_ambiguity
    eval_all = evaluate_worst_case_ambiguity(
        age_grid_years,
        prior_constraints,
        cand_list,
        target_functional,
        sigma_multiplier=sigma_multiplier,
        feasibility_tolerance=feasibility_tolerance,
    )
    # Keep solver failure distinct from a mathematical result.  ``all_w`` is
    # exported in the certificate and therefore must never be populated from
    # an unsuccessful LP.  The budgeted search cannot establish an optimum
    # when the all-candidate reference solve is unverified.
    if eval_all.status != "FEASIBLE":
        status_label = (
            "INCONSISTENT_EVIDENCE"
            if eval_all.status.startswith("INFEASIBLE")
            else "NUMERICAL_ERROR"
        )
        error_metadata = dict(metadata or {})
        error_metadata["all_candidates_evaluation_status"] = eval_all.status
        return CertifiedDesignCertificate(
            status=status_label,
            target_functional_name=target_functional.name,
            target_tolerance=target_tol if target_tol is not None else float("nan"),
            initial_ambiguity=initial_w,
            achieved_ambiguity=float("nan"),
            all_candidates_ambiguity=float("nan"),
            selected_option_ids=(),
            total_cost=0.0,
            candidate_count=len(cand_list),
            feasibility_status=eval_all.status,
            evaluations_count=2,
            metadata=error_metadata,
        )
    all_w = eval_all.worst_case_ambiguity

    best_subset: tuple[CertifiedCandidateTracer, ...] = ()
    best_cost = 0.0
    best_w = initial_w
    best_witnesses = (eval_empty.lower_witness, eval_empty.upper_witness)
    evaluations = 2

    # Check all affordable subsets
    for r in range(1, len(cand_list) + 1):
        for combo in itertools.combinations(cand_list, r):
            cost = subset_cost(combo)
            if cost > float(budget) + 1.0e-12:
                continue

            ev = evaluate_worst_case_ambiguity(
                age_grid_years,
                prior_constraints,
                combo,
                target_functional,
                sigma_multiplier=sigma_multiplier,
                feasibility_tolerance=feasibility_tolerance,
            )
            evaluations += 1

            if ev.status == "FEASIBLE" and (
                ev.worst_case_ambiguity < best_w - 1.0e-9
                or (
                    abs(ev.worst_case_ambiguity - best_w) <= 1.0e-9
                    and cost < best_cost
                )
            ):
                best_w = ev.worst_case_ambiguity
                best_cost = cost
                best_subset = combo
                best_witnesses = (ev.lower_witness, ev.upper_witness)

    if target_tol is not None:
        if best_w <= target_tol + feasibility_tolerance:
            status = "CERTIFIED_SUFFICIENT" if len(best_subset) > 0 else "ALREADY_RESOLVED"
        else:
            if best_w < initial_w - 1.0e-9:
                status = "BUDGET_EXHAUSTED_INSUFFICIENT"
            else:
                status = "NO_ACTION_AFFORDABLE"
    else:
        status = "CERTIFIED_SUFFICIENT" if best_w < initial_w - 1.0e-9 else "ALREADY_RESOLVED"

    return CertifiedDesignCertificate(
        status=status,
        target_functional_name=target_functional.name,
        target_tolerance=target_tol if target_tol is not None else float("nan"),
        initial_ambiguity=initial_w,
        achieved_ambiguity=best_w,
        all_candidates_ambiguity=all_w,
        selected_option_ids=tuple(c.option_id for c in best_subset),
        total_cost=best_cost,
        candidate_count=len(cand_list),
        feasibility_status="FEASIBLE",
        lower_witness=best_witnesses[0],
        upper_witness=best_witnesses[1],
        evaluations_count=evaluations,
        metadata=dict(metadata or {}),
    )


__all__ = [
    "AmbiguityEvaluation",
    "CertifiedCandidateTracer",
    "CertifiedDesignCertificate",
    "evaluate_worst_case_ambiguity",
    "solve_budgeted_minimax_design",
    "solve_certified_measurement_design",
]
