"""Active and Certified Measurement Design (ACMD) for partially identified groundwater systems.

This module unifies robust minimax ambiguity certification and probabilistic
Bayesian optimal experimental design into a single, closed-loop active learning
architecture.

Architecture:
    [Active Certified Measurement Design (ACMD)]
    Two complementary decision criteria:
        1. Robust criterion: minimise worst-case ambiguity W(S) over compatible set P_t
        2. Probabilistic criterion: maximise expected information gain (EIG) / utility
    Subject to:
        Cost (base + travel mobilization), feasibility, measurement error,
        accessibility, and identifiability.

Active closed-loop cycle:
    P_t  --> candidate measurements --> [ Minimax Certification / Bayesian EIG ]
         --> select x_{t+1} --> observe y_{t+1} --> update to P_{t+1} --> ...
    until either:
        W_t <= delta  --> declare CERTIFIED (with cryptographic certificate), OR
        candidate set cannot attain delta --> emit constructive IMPOSSIBILITY_WITNESS.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import hashlib
import hmac
import json
import math
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from .calibration.bayesian_active_learning import (
    AcquisitionConfig,
    MeasurementOption,
    PredictiveScenario,
    expected_information_gain,
)
from .nuclear.joint_lpm import tracer_response_kernel
from .nuclear.ttd_certified_design import (
    AmbiguityEvaluation,
    CertifiedCandidateTracer,
    CertifiedDesignCertificate,
    evaluate_worst_case_ambiguity,
    solve_budgeted_minimax_design,
    solve_certified_measurement_design,
)
from .nuclear.ttd_design import (
    TtdHypothesisEnsemble,
    _probability_gate,
)
from .nuclear.ttd_identified import (
    AgeFunctional,
    TracerConstraint,
    _readonly_vector,
)
from .acmd_polytope import PolyhedralStateSpace, evaluate_polyhedral_ambiguity


def _certificate_json_value(value: Any) -> Any:
    """Represent non-finite diagnostic values without nonstandard JSON NaN tokens."""
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, np.ndarray):
        return [_certificate_json_value(item) for item in value.tolist()]
    if isinstance(value, Mapping):
        return {str(key): _certificate_json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_certificate_json_value(item) for item in value]
    return value


class ACMDMode(str, Enum):
    """Decision criteria supported by the ACMD architecture."""

    ROBUST_MINIMAX = "robust_minimax"
    PROBABILISTIC_EIG = "probabilistic_eig"
    HYBRID = "hybrid"


class ACMDStatus(str, Enum):
    """Auditable state of the ACMD active learning cycle."""

    INITIALIZED = "INITIALIZED"
    ACTIVE = "ACTIVE"
    CERTIFIED = "CERTIFIED"
    ALREADY_RESOLVED = "ALREADY_RESOLVED"
    IMPOSSIBILITY_WITNESS = "IMPOSSIBILITY_WITNESS"
    BUDGET_EXHAUSTED = "BUDGET_EXHAUSTED"
    INCONSISTENT_EVIDENCE = "INCONSISTENT_EVIDENCE"
    ABSTAIN = "ABSTAIN"


@dataclass(frozen=True)
class ACMDAction:
    """A concrete candidate measurement action evaluated by ACMD."""

    action_id: str
    measurement_type: str
    target_id: str
    cost: float = 1.0
    travel_cost: float = 0.0
    accessibility: float = 1.0
    feasible: bool = True
    error_bound: float = 0.1
    standard_deviation: float = 0.05
    sample_year: float = 2026.0
    response: Optional[Sequence[float]] = None
    scenarios: Sequence[PredictiveScenario] = ()
    kernel_kwargs: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not str(self.action_id).strip():
            raise ValueError("action_id must be non-empty.")
        if not str(self.measurement_type).strip():
            raise ValueError("measurement_type must be non-empty.")
        if not str(self.target_id).strip():
            raise ValueError("target_id must be non-empty.")
        if not math.isfinite(float(self.cost)) or float(self.cost) <= 0.0:
            raise ValueError("cost must be finite and strictly positive.")
        if not math.isfinite(float(self.travel_cost)) or float(self.travel_cost) < 0.0:
            raise ValueError("travel_cost must be finite and non-negative.")
        if not math.isfinite(float(self.accessibility)) or float(self.accessibility) <= 0.0:
            raise ValueError("accessibility must be finite and strictly positive.")
        if not math.isfinite(float(self.error_bound)) or float(self.error_bound) <= 0.0:
            raise ValueError("error_bound must be finite and strictly positive.")
        if not math.isfinite(float(self.standard_deviation)) or float(self.standard_deviation) <= 0.0:
            raise ValueError("standard_deviation must be finite and strictly positive.")
        if self.response is not None:
            object.__setattr__(
                self, "response", _readonly_vector(self.response, name="response")
            )
        object.__setattr__(self, "scenarios", tuple(self.scenarios))
        object.__setattr__(self, "kernel_kwargs", dict(self.kernel_kwargs))
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def standalone_cost(self) -> float:
        """Total standalone cost including mobilization without shared logistics."""
        return float((self.cost + self.travel_cost) / self.accessibility)

    def incremental_cost(self, visited_targets: Set[str]) -> float:
        """Incremental cost accounting for already-paid travel mobilization."""
        travel = 0.0 if self.target_id in visited_targets else self.travel_cost
        return float((self.cost + travel) / self.accessibility)

    def resolve_response(self, age_grid_years: np.ndarray) -> np.ndarray:
        """Resolve linear kernel response vector on the given age grid."""
        if self.response is not None:
            resp = np.asarray(self.response, dtype=float)
            if resp.shape != age_grid_years.shape:
                raise ValueError(
                    f"Action {self.action_id!r} response length {resp.size} "
                    f"does not match age grid size {age_grid_years.size}."
                )
            return resp
        kernel = tracer_response_kernel(
            self.measurement_type,
            age_grid_years,
            float(self.sample_year),
            **self.kernel_kwargs,
        )
        return np.asarray(kernel, dtype=float)

    def to_certified_candidate(self, age_grid_years: np.ndarray) -> CertifiedCandidateTracer:
        """Convert into a CertifiedCandidateTracer for minimax linear programming."""
        resp = self.resolve_response(age_grid_years)
        return CertifiedCandidateTracer(
            option_id=self.action_id,
            tracer=self.measurement_type,
            sample_year=float(self.sample_year),
            error_bound=float(self.error_bound),
            cost=float(self.standalone_cost),
            response=resp,
            kernel_kwargs=self.kernel_kwargs,
            metadata={
                **self.metadata,
                "target_id": self.target_id,
                "accessibility": float(self.accessibility),
            },
        )

    def to_measurement_option(self, default_scenarios: Sequence[PredictiveScenario]) -> MeasurementOption:
        """Convert into a MeasurementOption for Bayesian EIG acquisition."""
        scenarios = self.scenarios if self.scenarios else default_scenarios
        return MeasurementOption(
            option_id=self.action_id,
            measurement_type=self.measurement_type,
            target_id=self.target_id,
            cost=float(self.standalone_cost),
            scenarios=scenarios,
            feasible=bool(self.feasible),
            metadata={
                **self.metadata,
                "sample_year": float(self.sample_year),
                "accessibility": float(self.accessibility),
            },
        )


@dataclass(frozen=True)
class ACMDStepRecord:
    """Telemetry and outcome of one discrete step in the ACMD active cycle."""

    step_index: int
    selected_action_id: str
    target_id: str
    observed_value: float
    error_bound: float
    incremental_cost: float
    cumulative_cost: float
    prior_ambiguity: float
    posterior_ambiguity: float
    ambiguity_reduction: float
    status_after_step: ACMDStatus
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "step_index": int(self.step_index),
            "selected_action_id": self.selected_action_id,
            "target_id": self.target_id,
            "observed_value": float(self.observed_value),
            "error_bound": float(self.error_bound),
            "incremental_cost": float(self.incremental_cost),
            "cumulative_cost": float(self.cumulative_cost),
            "prior_ambiguity": float(self.prior_ambiguity),
            "posterior_ambiguity": float(self.posterior_ambiguity),
            "ambiguity_reduction": float(self.ambiguity_reduction),
            "status_after_step": str(self.status_after_step.value),
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class ACMDCertificate:
    """Cryptographically verifiable certificate of sufficiency, impossibility, or state."""

    status: ACMDStatus
    target_functional_name: str
    target_tolerance: float
    initial_ambiguity: float
    final_ambiguity: float
    candidate_attainable_ambiguity: float
    total_cost: float
    budget: Optional[float]
    steps: Tuple[ACMDStepRecord, ...]
    selected_action_ids: Tuple[str, ...]
    visited_targets: Tuple[str, ...]
    lower_witness: Optional[np.ndarray] = None
    upper_witness: Optional[np.ndarray] = None
    certificate_hash: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "steps", tuple(self.steps))
        object.__setattr__(self, "selected_action_ids", tuple(self.selected_action_ids))
        object.__setattr__(self, "visited_targets", tuple(self.visited_targets))
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
            object.__setattr__(self, "certificate_hash", self.compute_hash())

    def compute_hash(self) -> str:
        """Recompute the digest from the result and compiled problem fingerprint."""
        payload = json.dumps(
            _certificate_json_value({
                "status": str(self.status.value),
                "target_functional": self.target_functional_name,
                "target_tolerance": float(self.target_tolerance),
                "initial_ambiguity": self.initial_ambiguity,
                "final_ambiguity": self.final_ambiguity,
                "candidate_attainable_ambiguity": self.candidate_attainable_ambiguity,
                "selected_action_ids": list(self.selected_action_ids),
                "visited_targets": list(self.visited_targets),
                "total_cost": float(self.total_cost),
                "budget": self.budget,
                "steps": [step.to_dict() for step in self.steps],
                "lower_witness": None if self.lower_witness is None else self.lower_witness.tolist(),
                "upper_witness": None if self.upper_witness is None else self.upper_witness.tolist(),
                "problem_hash": self.metadata.get("problem_hash"),
            }),
            sort_keys=True,
            allow_nan=False,
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def verify_hash(self) -> bool:
        """Check whether result fields still match the recorded SHA-256 digest."""
        return hmac.compare_digest(self.certificate_hash, self.compute_hash())

    def to_dict(self, *, include_witnesses: bool = True) -> dict[str, Any]:
        data: dict[str, Any] = {
            "status": str(self.status.value),
            "target_functional_name": self.target_functional_name,
            "target_tolerance": float(self.target_tolerance),
            "initial_ambiguity": float(self.initial_ambiguity),
            "final_ambiguity": float(self.final_ambiguity),
            "candidate_attainable_ambiguity": float(self.candidate_attainable_ambiguity),
            "total_cost": float(self.total_cost),
            "budget": float(self.budget) if self.budget is not None else None,
            "selected_action_ids": list(self.selected_action_ids),
            "visited_targets": list(self.visited_targets),
            "steps": [s.to_dict() for s in self.steps],
            "certificate_hash": self.certificate_hash,
            "metadata": dict(self.metadata),
        }
        if include_witnesses and self.lower_witness is not None and self.upper_witness is not None:
            data["lower_witness"] = self.lower_witness.tolist()
            data["upper_witness"] = self.upper_witness.tolist()
        return data


class ACMDLoop:
    """Stateful dynamic active cycle runner for Active Certified Measurement Design.

    Implements the closed loop:
        P_t --> evaluate W_t --> recommend x_{t+1} --> observe y_{t+1} --> update to P_{t+1}
    until:
        W_t <= delta (CERTIFIED),
        or complete remaining candidate set cannot attain delta (IMPOSSIBILITY_WITNESS),
        or budget exhausted, or inconsistent evidence detected.
    """

    def __init__(
        self,
        age_grid_years: Sequence[float],
        initial_constraints: Sequence[TracerConstraint],
        candidates: Sequence[ACMDAction],
        target_functional: AgeFunctional,
        target_tolerance: float,
        *,
        budget: Optional[float] = None,
        cost_exponent: float = 1.0,
        mode: ACMDMode = ACMDMode.ROBUST_MINIMAX,
        sigma_multiplier: float = 1.96,
        feasibility_tolerance: float = 1.0e-7,
        hypothesis_ensemble: Optional[TtdHypothesisEnsemble] = None,
        acquisition_config: Optional[AcquisitionConfig] = None,
        metadata: Optional[Mapping[str, Any]] = None,
        state_space: Optional[PolyhedralStateSpace] = None,
    ) -> None:
        self.age_grid = _readonly_vector(age_grid_years, name="age_grid_years")
        if state_space is not None:
            if state_space.dimension != self.age_grid.size:
                raise ValueError("state_space dimension must match the response vector length.")
            if mode != ACMDMode.ROBUST_MINIMAX:
                raise ValueError("Polyhedral state spaces currently support robust_minimax mode only.")
            if any(candidate.response is None for candidate in candidates):
                raise ValueError("Polyhedral state spaces require explicit candidate responses.")
        self.state_space = state_space
        self.constraints: List[TracerConstraint] = list(initial_constraints)
        self.candidates: Dict[str, ACMDAction] = {c.action_id: c for c in candidates}
        self.target_functional = target_functional
        self.target_tolerance = float(target_tolerance)
        if not math.isfinite(self.target_tolerance) or self.target_tolerance <= 0.0:
            raise ValueError("target_tolerance must be finite and strictly positive.")

        self.budget = float(budget) if budget is not None else None
        if self.budget is not None and (not math.isfinite(self.budget) or self.budget < 0.0):
            raise ValueError("budget must be finite and non-negative.")

        self.cost_exponent = float(cost_exponent)
        self.mode = mode
        self.sigma_multiplier = float(sigma_multiplier)
        if not math.isfinite(self.sigma_multiplier) or self.sigma_multiplier <= 0.0:
            raise ValueError("sigma_multiplier must be finite and positive.")
        self.feasibility_tolerance = float(feasibility_tolerance)
        if not math.isfinite(self.feasibility_tolerance) or self.feasibility_tolerance <= 0.0:
            raise ValueError("feasibility_tolerance must be finite and positive.")
        self.hypothesis_ensemble = hypothesis_ensemble
        self.acquisition_config = acquisition_config or AcquisitionConfig()
        self.metadata = dict(metadata or {})

        # Runtime state
        self.steps: List[ACMDStepRecord] = []
        self.visited_targets: Set[str] = set()
        self.cumulative_cost: float = 0.0
        self.status = ACMDStatus.INITIALIZED

        # Evaluate initial ambiguity W_0
        initial_eval = self.evaluate_polytope_ambiguity(candidate_subset=())
        if initial_eval.status != "FEASIBLE":
            self.status = (
                ACMDStatus.INCONSISTENT_EVIDENCE
                if initial_eval.status.startswith("INFEASIBLE")
                else ACMDStatus.ABSTAIN
            )
            self.initial_ambiguity = float("nan")
            self.current_ambiguity = float("nan")
        else:
            self.initial_ambiguity = initial_eval.worst_case_ambiguity
            self.current_ambiguity = initial_eval.worst_case_ambiguity
            if self.current_ambiguity <= self.target_tolerance + self.feasibility_tolerance:
                self.status = ACMDStatus.ALREADY_RESOLVED
            else:
                self.status = ACMDStatus.ACTIVE

    @property
    def is_terminated(self) -> bool:
        """Return True if the active cycle has reached a definitive stopping condition."""
        return self.status not in (ACMDStatus.INITIALIZED, ACMDStatus.ACTIVE)

    def evaluate_polytope_ambiguity(
        self, candidate_subset: Sequence[CertifiedCandidateTracer] = ()
    ) -> AmbiguityEvaluation:
        """Solve joint linear program to find worst-case ambiguity W(S) over P_t."""
        if self.state_space is not None:
            return evaluate_polyhedral_ambiguity(
                self.state_space,
                self.constraints,
                candidate_subset,
                self.target_functional,
                sigma_multiplier=self.sigma_multiplier,
                feasibility_tolerance=self.feasibility_tolerance,
            )
        return evaluate_worst_case_ambiguity(
            self.age_grid,
            self.constraints,
            candidate_subset,
            self.target_functional,
            sigma_multiplier=self.sigma_multiplier,
            feasibility_tolerance=self.feasibility_tolerance,
        )

    def check_all_candidates_attainability(self) -> Tuple[float, Optional[AmbiguityEvaluation]]:
        """Evaluate whether measuring ALL remaining candidates can achieve target_tolerance."""
        feasible_remaining = [
            c.to_certified_candidate(self.age_grid)
            for c in self.candidates.values()
            if c.feasible
        ]
        eval_all = self.evaluate_polytope_ambiguity(feasible_remaining)
        if eval_all.status != "FEASIBLE":
            return float("nan"), eval_all
        return eval_all.worst_case_ambiguity, eval_all

    def recommend_next_action(self) -> Tuple[Optional[ACMDAction], Dict[str, Any]]:
        """Recommend optimal next measurement action x_{t+1} using declared criterion."""
        if self.is_terminated:
            return None, {"status": self.status.value, "reason": "loop_already_terminated"}

        # 1. Check mode-specific preconditions (e.g. valid probability gate for Bayesian EIG)
        if self.mode == ACMDMode.PROBABILISTIC_EIG:
            if self.hypothesis_ensemble is None:
                self.status = ACMDStatus.ABSTAIN
                return None, {
                    "status": ACMDStatus.ABSTAIN.value,
                    "reason": "no_hypothesis_ensemble_provided_for_eig_mode",
                }
            gate_check = _probability_gate(self.hypothesis_ensemble)
            if gate_check is not None:
                self.status = ACMDStatus.ABSTAIN
                return None, gate_check

        # 2. Check if remaining candidates can attain delta (for minimax and hybrid certification)
        if self.mode in (ACMDMode.ROBUST_MINIMAX, ACMDMode.HYBRID):
            attainable_w, eval_all = self.check_all_candidates_attainability()
            if (
                math.isfinite(attainable_w)
                and attainable_w > self.target_tolerance + self.feasibility_tolerance
            ):
                self.status = ACMDStatus.IMPOSSIBILITY_WITNESS
                return None, {
                    "status": ACMDStatus.IMPOSSIBILITY_WITNESS.value,
                    "reason": "candidates_cannot_achieve_target_tolerance",
                    "attainable_ambiguity": attainable_w,
                    "target_tolerance": self.target_tolerance,
                    "witness_pair": (
                        eval_all.upper_witness.tolist()
                        if eval_all and eval_all.upper_witness is not None
                        else None,
                        eval_all.lower_witness.tolist()
                        if eval_all and eval_all.lower_witness is not None
                        else None,
                    ),
                }

        # Filter feasible candidates
        feasible_candidates = [c for c in self.candidates.values() if c.feasible]
        if not feasible_candidates:
            self.status = ACMDStatus.IMPOSSIBILITY_WITNESS
            return None, {
                "status": ACMDStatus.IMPOSSIBILITY_WITNESS.value,
                "reason": "no_feasible_candidates_remain",
            }

        # Check budget availability
        affordable_candidates = []
        for c in feasible_candidates:
            inc_cost = c.incremental_cost(self.visited_targets)
            if self.budget is None or (self.cumulative_cost + inc_cost <= self.budget + 1.0e-12):
                affordable_candidates.append(c)

        if not affordable_candidates:
            self.status = ACMDStatus.BUDGET_EXHAUSTED
            return None, {
                "status": ACMDStatus.BUDGET_EXHAUSTED.value,
                "reason": "insufficient_remaining_budget",
            }

        # Score candidates based on selected mode
        scores: Dict[str, float] = {}
        details: Dict[str, Any] = {}

        if self.mode == ACMDMode.ROBUST_MINIMAX:
            # Score by worst-case ambiguity reduction per unit cost
            for cand in affordable_candidates:
                cert_cand = cand.to_certified_candidate(self.age_grid)
                eval_cand = self.evaluate_polytope_ambiguity([cert_cand])
                if eval_cand.status != "FEASIBLE":
                    continue
                post_w = eval_cand.worst_case_ambiguity
                reduction = max(0.0, self.current_ambiguity - post_w)
                inc_cost = cand.incremental_cost(self.visited_targets)
                score = reduction / (inc_cost ** self.cost_exponent)
                scores[cand.action_id] = score
                details[cand.action_id] = {
                    "reduction": reduction,
                    "post_ambiguity": post_w,
                    "cost": inc_cost,
                    "score": score,
                }

        elif self.mode == ACMDMode.PROBABILISTIC_EIG:
            # Check probability gate
            if self.hypothesis_ensemble is None:
                self.status = ACMDStatus.ABSTAIN
                return None, {
                    "status": ACMDStatus.ABSTAIN.value,
                    "reason": "no_hypothesis_ensemble_provided_for_eig_mode",
                }
            gate_check = _probability_gate(self.hypothesis_ensemble)
            if gate_check is not None:
                self.status = ACMDStatus.ABSTAIN
                return None, gate_check

            options = [
                c.to_measurement_option(default_scenarios=())
                for c in affordable_candidates
            ]
            options_by_id = {opt.option_id: opt for opt in options}
            prior_probs = np.asarray(self.hypothesis_ensemble.probabilities, dtype=float)

            for cand in affordable_candidates:
                opt = options_by_id[cand.action_id]
                # If explicit scenarios exist, compute nominal EIG
                if opt.scenarios:
                    scenario = opt.scenarios[0]
                    means = np.asarray(scenario.means, dtype=float)
                    sds = np.asarray(scenario.standard_deviations, dtype=float)
                    eig_val = expected_information_gain(
                        prior_probs,
                        means,
                        sds,
                        quadrature_order=self.acquisition_config.quadrature_order,
                    )
                else:
                    # Construct forward projection from hypothesis masses
                    resp = cand.resolve_response(self.age_grid)
                    means = np.asarray(self.hypothesis_ensemble.masses) @ resp
                    sds = float(cand.standard_deviation)
                    eig_val = expected_information_gain(
                        prior_probs,
                        means,
                        sds,
                        quadrature_order=self.acquisition_config.quadrature_order,
                    )
                inc_cost = cand.incremental_cost(self.visited_targets)
                score = eig_val / (inc_cost ** self.cost_exponent)
                scores[cand.action_id] = score
                details[cand.action_id] = {
                    "eig": eig_val,
                    "cost": inc_cost,
                    "score": score,
                }

        elif self.mode == ACMDMode.HYBRID:
            # Hybrid: filter candidates that strictly reduce worst-case ambiguity,
            # then prioritize by EIG among ambiguity reducers.
            for cand in affordable_candidates:
                cert_cand = cand.to_certified_candidate(self.age_grid)
                eval_cand = self.evaluate_polytope_ambiguity([cert_cand])
                if eval_cand.status != "FEASIBLE":
                    continue
                reduction = max(0.0, self.current_ambiguity - eval_cand.worst_case_ambiguity)
                inc_cost = cand.incremental_cost(self.visited_targets)

                eig_val = 0.0
                if self.hypothesis_ensemble is not None and self.hypothesis_ensemble.probabilities is not None:
                    prior_probs = np.asarray(self.hypothesis_ensemble.probabilities, dtype=float)
                    resp = cand.resolve_response(self.age_grid)
                    means = np.asarray(self.hypothesis_ensemble.masses) @ resp
                    eig_val = expected_information_gain(
                        prior_probs,
                        means,
                        float(cand.standard_deviation),
                        quadrature_order=self.acquisition_config.quadrature_order,
                    )

                # Lexicographic score: primary weight on ambiguity reduction, secondary on EIG
                score = (reduction + 0.1 * eig_val) / (inc_cost ** self.cost_exponent)
                scores[cand.action_id] = score
                details[cand.action_id] = {
                    "reduction": reduction,
                    "eig": eig_val,
                    "cost": inc_cost,
                    "score": score,
                }

        if not scores:
            self.status = ACMDStatus.IMPOSSIBILITY_WITNESS
            return None, {
                "status": ACMDStatus.IMPOSSIBILITY_WITNESS.value,
                "reason": "no_candidates_produced_valid_score",
            }

        # Pick candidate with highest score
        best_id = max(scores, key=lambda k: scores[k])
        best_action = self.candidates[best_id]
        return best_action, {
            "selected_action_id": best_id,
            "score": scores[best_id],
            "mode": self.mode.value,
            "all_scores": scores,
            "details": details,
        }

    def step(
        self,
        action_id: str,
        observed_value: float,
        error_bound: Optional[float] = None,
        sigma: Optional[float] = None,
    ) -> ACMDStepRecord:
        """Ingest new observation y_{t+1}, update polytope P_{t+1}, and update loop state."""
        if action_id not in self.candidates:
            raise KeyError(f"Action {action_id!r} is not an available candidate in ACMDLoop.")
        if self.is_terminated:
            raise RuntimeError("Cannot add an observation to a terminated ACMD loop.")
        if not math.isfinite(float(observed_value)):
            raise ValueError("observed_value must be finite.")
        action = self.candidates[action_id]
        if not action.feasible:
            raise ValueError(f"Action {action_id!r} is marked infeasible.")
        bound = float(error_bound) if error_bound is not None else float(action.error_bound)
        if not math.isfinite(bound) or bound <= 0.0:
            raise ValueError("error_bound must be finite and positive.")
        if sigma is not None and (not math.isfinite(float(sigma)) or float(sigma) <= 0.0):
            raise ValueError("sigma must be finite and positive.")
        if sigma is not None:
            sigma_bound = float(sigma) * self.sigma_multiplier
            if error_bound is not None and not math.isclose(
                bound, sigma_bound, rel_tol=1.0e-9, abs_tol=1.0e-12
            ):
                raise ValueError("error_bound and sigma imply different LP intervals.")
            bound = sigma_bound
        polytope_sig = float(sigma) if sigma is not None else (bound / self.sigma_multiplier)
        bayes_sig = float(sigma) if sigma is not None else float(action.standard_deviation)
        inc_cost = action.incremental_cost(self.visited_targets)
        if self.budget is not None and self.cumulative_cost + inc_cost > self.budget + 1.0e-12:
            raise ValueError("Action exceeds remaining ACMD budget.")

        self.candidates.pop(action_id)

        # 1. Update costs and logistics
        self.cumulative_cost += inc_cost
        self.visited_targets.add(action.target_id)

        # 2. Update polyhedral constraints P_{t+1} = P_t \cap {x : |h^T x - y| <= epsilon}
        resp = action.resolve_response(self.age_grid)
        new_constraint = TracerConstraint(
            tracer=f"{action.measurement_type}_{action.action_id}",
            response=resp,
            observed=float(observed_value),
            sigma=polytope_sig,
            units=action.metadata.get("units", ""),
            metadata={"action_id": action.action_id, "target_id": action.target_id},
        )
        self.constraints.append(new_constraint)

        # 3. In Bayesian mode, update hypothesis posterior probabilities if available
        if (
            self.hypothesis_ensemble is not None
            and self.hypothesis_ensemble.probabilities is not None
        ):
            prior_p = np.asarray(self.hypothesis_ensemble.probabilities, dtype=float)
            means = np.asarray(self.hypothesis_ensemble.masses) @ resp
            # Gaussian likelihood
            log_lik = -0.5 * ((float(observed_value) - means) / bayes_sig) ** 2 - np.log(bayes_sig * np.sqrt(2 * np.pi))
            log_post = np.log(np.maximum(prior_p, 1e-300)) + log_lik
            log_post -= np.max(log_post)
            post_p = np.exp(log_post)
            post_p /= post_p.sum()
            object.__setattr__(self.hypothesis_ensemble, "probabilities", post_p)

        # 4. Re-evaluate worst-case ambiguity W_{t+1}
        prior_w = self.current_ambiguity
        eval_post = self.evaluate_polytope_ambiguity(candidate_subset=())

        if eval_post.status != "FEASIBLE":
            self.status = (
                ACMDStatus.INCONSISTENT_EVIDENCE
                if eval_post.status.startswith("INFEASIBLE")
                else ACMDStatus.ABSTAIN
            )
            post_w = float("nan")
            reduction = float("nan")
            self.current_ambiguity = float("nan")
        else:
            post_w = eval_post.worst_case_ambiguity
            reduction = max(0.0, prior_w - post_w)
            self.current_ambiguity = post_w

            # Check convergence to target tolerance delta
            if post_w <= self.target_tolerance + self.feasibility_tolerance:
                self.status = ACMDStatus.CERTIFIED
            elif self.budget is not None and self.cumulative_cost >= self.budget:
                self.status = ACMDStatus.BUDGET_EXHAUSTED
            else:
                self.status = ACMDStatus.ACTIVE

        step_record = ACMDStepRecord(
            step_index=len(self.steps) + 1,
            selected_action_id=action.action_id,
            target_id=action.target_id,
            observed_value=float(observed_value),
            error_bound=bound,
            incremental_cost=inc_cost,
            cumulative_cost=self.cumulative_cost,
            prior_ambiguity=prior_w,
            posterior_ambiguity=post_w,
            ambiguity_reduction=reduction,
            status_after_step=self.status,
            metadata={"measurement_type": action.measurement_type},
        )
        self.steps.append(step_record)
        return step_record

    def run_simulation(
        self,
        oracle_fn: Callable[[ACMDAction], float],
        *,
        max_steps: Optional[int] = None,
    ) -> ACMDCertificate:
        """Run the dynamic active cycle automatically against an observation oracle."""
        steps_taken = 0
        while not self.is_terminated:
            if max_steps is not None and steps_taken >= max_steps:
                break
            action, _ = self.recommend_next_action()
            if action is None:
                break
            obs_y = float(oracle_fn(action))
            self.step(action.action_id, obs_y)
            steps_taken += 1

        return self.get_certificate()

    def get_certificate(self) -> ACMDCertificate:
        """Compile final auditable ACMDCertificate with witnesses and SHA-256 hash."""
        # Find all-candidates attainability
        attainable_w, eval_all = self.check_all_candidates_attainability()
        lower_w = eval_all.lower_witness if eval_all is not None else None
        upper_w = eval_all.upper_witness if eval_all is not None else None

        return ACMDCertificate(
            status=self.status,
            target_functional_name=self.target_functional.name,
            target_tolerance=self.target_tolerance,
            initial_ambiguity=self.initial_ambiguity,
            final_ambiguity=self.current_ambiguity,
            candidate_attainable_ambiguity=attainable_w,
            total_cost=self.cumulative_cost,
            budget=self.budget,
            steps=tuple(self.steps),
            selected_action_ids=tuple(s.selected_action_id for s in self.steps),
            visited_targets=tuple(sorted(self.visited_targets)),
            lower_witness=lower_w,
            upper_witness=upper_w,
            metadata=self.metadata,
        )


class ACMD:
    """Unified Facade for Active and Certified Measurement Design (ACMD)."""

    @staticmethod
    def create_loop(
        age_grid_years: Sequence[float],
        initial_constraints: Sequence[TracerConstraint],
        candidates: Sequence[ACMDAction],
        target_functional: AgeFunctional,
        target_tolerance: float,
        *,
        budget: Optional[float] = None,
        cost_exponent: float = 1.0,
        mode: ACMDMode = ACMDMode.ROBUST_MINIMAX,
        sigma_multiplier: float = 1.96,
        feasibility_tolerance: float = 1.0e-7,
        hypothesis_ensemble: Optional[TtdHypothesisEnsemble] = None,
        acquisition_config: Optional[AcquisitionConfig] = None,
        metadata: Optional[Mapping[str, Any]] = None,
        state_space: Optional[PolyhedralStateSpace] = None,
    ) -> ACMDLoop:
        """Create a stateful, iterative active cycle runner."""
        return ACMDLoop(
            age_grid_years,
            initial_constraints,
            candidates,
            target_functional,
            target_tolerance,
            budget=budget,
            cost_exponent=cost_exponent,
            mode=mode,
            sigma_multiplier=sigma_multiplier,
            feasibility_tolerance=feasibility_tolerance,
            hypothesis_ensemble=hypothesis_ensemble,
            acquisition_config=acquisition_config,
            metadata=metadata,
            state_space=state_space,
        )

    @staticmethod
    def solve_static_certification(
        age_grid_years: Sequence[float],
        initial_constraints: Sequence[TracerConstraint],
        candidates: Sequence[CertifiedCandidateTracer],
        target_functional: AgeFunctional,
        target_tolerance: float,
        *,
        sigma_multiplier: float = 1.96,
        feasibility_tolerance: float = 1.0e-7,
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> CertifiedDesignCertificate:
        """Solve one-shot offline minimax certified design."""
        return solve_certified_measurement_design(
            age_grid_years,
            initial_constraints,
            candidates,
            target_functional,
            target_tolerance,
            sigma_multiplier=sigma_multiplier,
            feasibility_tolerance=feasibility_tolerance,
            metadata=metadata,
        )

    @staticmethod
    def solve_budgeted_minimax(
        age_grid_years: Sequence[float],
        initial_constraints: Sequence[TracerConstraint],
        candidates: Sequence[CertifiedCandidateTracer],
        budget: float,
        *,
        target_tolerance: Optional[float] = None,
        sigma_multiplier: float = 1.96,
        feasibility_tolerance: float = 1.0e-7,
        cost_function: Optional[Callable[[Sequence[CertifiedCandidateTracer]], float]] = None,
        metadata: Optional[Mapping[str, Any]] = None,
    ) -> AmbiguityEvaluation:
        """Solve one-shot offline budgeted minimax design."""
        return solve_budgeted_minimax_design(
            age_grid_years,
            initial_constraints,
            candidates,
            budget,
            target_tolerance=target_tolerance,
            sigma_multiplier=sigma_multiplier,
            feasibility_tolerance=feasibility_tolerance,
            cost_function=cost_function,
            metadata=metadata,
        )
