"""Truth-blind one-step decision utility for candidate measurements.

This module ranks next-measurement actions from declared uncertainty models,
not from realised labels or hidden truth.  Each action carries a cost,
feasibility flag, and either outcome likelihoods, outcome-specific posterior
distributions, or a scenario-specific posterior distribution.  The selector
computes expected information gain (entropy reduction), divides by cost, and
uses the worst scenario as the robust decision utility.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
import random
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence, Tuple


Distribution = Mapping[str, float]


def _stable_hash(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _paired_bootstrap_mean_interval(
    values: Sequence[float],
    *,
    seed_material: object,
    replicates: int = 2000,
) -> tuple[float | None, float | None, int | None]:
    """Return a deterministic percentile interval for paired mean deltas."""

    observed = [float(value) for value in values if math.isfinite(float(value))]
    if len(observed) < 2 or int(replicates) < 1:
        return None, None, None
    seed_digest = _stable_hash(seed_material)
    seed = int(seed_digest[:16], 16)
    rng = random.Random(seed)
    n = len(observed)
    means = [
        sum(observed[rng.randrange(n)] for _ in range(n)) / n
        for _ in range(int(replicates))
    ]
    means.sort()

    def percentile(probability: float) -> float:
        index = int(round(probability * (len(means) - 1)))
        return float(means[min(max(index, 0), len(means) - 1)])

    return percentile(0.025), percentile(0.975), seed


def _contains_forbidden_truth_key(value: object) -> bool:
    forbidden = {"truth", "true_hypothesis", "ground_truth", "reference_label"}
    if isinstance(value, Mapping):
        return any(
            str(key).strip().lower() in forbidden
            or _contains_forbidden_truth_key(nested)
            for key, nested in value.items()
        )
    if isinstance(value, (tuple, list)):
        return any(_contains_forbidden_truth_key(item) for item in value)
    return False


def _numeric_or_negative_infinity(value: float | None) -> float:
    """Keep a valid zero utility distinct from a missing utility."""

    return -float("inf") if value is None else float(value)


def _as_distribution(values: Mapping[str, object], *, name: str) -> dict[str, float]:
    """Validate and normalise a finite non-negative discrete distribution."""

    if not values:
        raise ValueError(f"{name} must contain at least one entry")
    out: dict[str, float] = {}
    total = 0.0
    for key, raw in values.items():
        if key is None or str(key) == "":
            raise ValueError(f"{name} contains an empty state/outcome id")
        try:
            value = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name}[{key!r}] is not numeric") from exc
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{name}[{key!r}] must be finite and non-negative")
        out[str(key)] = value
        total += value
    if total <= 0.0:
        raise ValueError(f"{name} must have positive total probability")
    return {key: value / total for key, value in out.items()}


def entropy_bits(distribution: Mapping[str, object]) -> float:
    """Return Shannon entropy in bits for a discrete distribution."""

    probs = _as_distribution(distribution, name="distribution")
    entropy = 0.0
    for probability in probs.values():
        if probability > 0.0:
            entropy -= probability * math.log2(probability)
    return entropy


@dataclass(frozen=True)
class ScenarioBelief:
    """Declared belief model for an action under one structural scenario.

    Exactly one of these evidence specifications is needed:

    * ``outcome_likelihoods``: mapping ``outcome -> hypothesis -> P(outcome |
      hypothesis)``.  Outcome probabilities and posteriors are derived by
      Bayes' rule from ``prior``.
    * ``posterior_by_outcome`` with ``outcome_probabilities``: declared
      posterior distribution for each possible measurement outcome.
    * ``posterior``: a single scenario-specific posterior distribution.  This
      is treated as one deterministic declared outcome.
    """

    prior: Distribution
    outcome_likelihoods: Optional[Mapping[str, Distribution]] = None
    outcome_probabilities: Optional[Distribution] = None
    posterior_by_outcome: Optional[Mapping[str, Distribution]] = None
    posterior: Optional[Distribution] = None

    def resolved_outcomes(self) -> tuple[dict[str, float], dict[str, dict[str, float]]]:
        """Return declared outcome probabilities and posteriors.

        No realised outcome or ground truth is accepted here; all quantities are
        the model declaration available before the measurement is chosen.
        """

        prior = _as_distribution(self.prior, name="prior")
        specs = [
            self.outcome_likelihoods is not None,
            self.posterior_by_outcome is not None,
            self.posterior is not None,
        ]
        if sum(specs) != 1:
            raise ValueError(
                "provide exactly one of outcome_likelihoods, "
                "posterior_by_outcome, or posterior"
            )

        if self.outcome_likelihoods is not None:
            return _resolve_likelihood_outcomes(prior, self.outcome_likelihoods)

        if self.posterior_by_outcome is not None:
            if self.outcome_probabilities is None:
                raise ValueError(
                    "outcome_probabilities are required with posterior_by_outcome"
                )
            outcome_probs = _as_distribution(
                self.outcome_probabilities,
                name="outcome_probabilities",
            )
            missing = set(self.posterior_by_outcome) ^ set(outcome_probs)
            if missing:
                raise ValueError(
                    "posterior_by_outcome and outcome_probabilities must use "
                    "the same outcome ids"
                )
            posteriors = {
                str(outcome): _as_distribution(posterior, name=f"posterior[{outcome}]")
                for outcome, posterior in self.posterior_by_outcome.items()
            }
            return outcome_probs, posteriors

        return {"posterior": 1.0}, {
            "posterior": _as_distribution(self.posterior or {}, name="posterior")
        }


def _resolve_likelihood_outcomes(
    prior: Mapping[str, float],
    outcome_likelihoods: Mapping[str, Distribution],
) -> tuple[dict[str, float], dict[str, dict[str, float]]]:
    if not outcome_likelihoods:
        raise ValueError("outcome_likelihoods must contain at least one outcome")
    likelihoods: dict[str, dict[str, float]] = {}
    for outcome, by_state in outcome_likelihoods.items():
        if outcome is None or str(outcome) == "":
            raise ValueError("outcome_likelihoods contains an empty outcome id")
        likelihoods[str(outcome)] = {}
        for state in prior:
            try:
                likelihood = float(by_state[state])
            except KeyError as exc:
                raise ValueError(
                    f"outcome_likelihoods[{outcome!r}] is missing state {state!r}"
                ) from exc
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"outcome_likelihoods[{outcome!r}][{state!r}] is not numeric"
                ) from exc
            if not math.isfinite(likelihood) or likelihood < 0.0:
                raise ValueError(
                    f"outcome_likelihoods[{outcome!r}][{state!r}] must be "
                    "finite and non-negative"
                )
            likelihoods[str(outcome)][state] = likelihood

    for state in prior:
        total = sum(by_state[state] for by_state in likelihoods.values())
        if not math.isclose(total, 1.0, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"likelihoods across outcomes for state {state!r} must sum to 1"
            )

    outcome_probs: dict[str, float] = {}
    posteriors: dict[str, dict[str, float]] = {}
    for outcome, by_state in likelihoods.items():
        outcome_probability = sum(prior[state] * by_state[state] for state in prior)
        outcome_probs[outcome] = outcome_probability
        if outcome_probability > 0.0:
            posteriors[outcome] = {
                state: prior[state] * by_state[state] / outcome_probability
                for state in prior
            }
        else:
            posteriors[outcome] = dict(prior)
    return outcome_probs, posteriors


@dataclass(frozen=True)
class CandidateMeasurementAction:
    """One candidate measurement action available before field/lab execution."""

    action_id: str
    cost: float
    feasible: bool = True
    scenarios: Mapping[str, ScenarioBelief] = field(default_factory=dict)
    prior: Optional[Distribution] = None
    outcome_likelihoods: Optional[Mapping[str, Distribution]] = None
    outcome_probabilities: Optional[Distribution] = None
    posterior_by_outcome: Optional[Mapping[str, Distribution]] = None
    posterior: Optional[Distribution] = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def resolved_scenarios(self) -> dict[str, ScenarioBelief]:
        if not self.action_id:
            raise ValueError("action_id must be non-empty")
        if not math.isfinite(float(self.cost)) or float(self.cost) <= 0.0:
            raise ValueError(f"action {self.action_id!r} cost must be positive")
        if self.scenarios:
            return {str(key): value for key, value in self.scenarios.items()}
        if self.prior is None:
            raise ValueError(
                f"action {self.action_id!r} needs scenarios or a top-level prior"
            )
        return {
            "declared": ScenarioBelief(
                prior=self.prior,
                outcome_likelihoods=self.outcome_likelihoods,
                outcome_probabilities=self.outcome_probabilities,
                posterior_by_outcome=self.posterior_by_outcome,
                posterior=self.posterior,
            )
        }


@dataclass(frozen=True)
class ScenarioUtility:
    scenario_id: str
    prior_entropy_bits: float
    expected_posterior_entropy_bits: float
    information_gain_bits: float
    utility_per_cost: float
    outcome_probabilities: Mapping[str, float]
    posterior_entropy_by_outcome: Mapping[str, float]
    identifiability_status: str = "SCENARIO_ENVELOPE_NOT_CALIBRATED"

    def to_dict(self) -> dict[str, Any]:
        return {
            "scenario_id": self.scenario_id,
            "prior_entropy_bits": self.prior_entropy_bits,
            "expected_posterior_entropy_bits": self.expected_posterior_entropy_bits,
            "information_gain_bits": self.information_gain_bits,
            "utility_per_cost": self.utility_per_cost,
            "outcome_probabilities": dict(self.outcome_probabilities),
            "posterior_entropy_by_outcome": dict(self.posterior_entropy_by_outcome),
            "identifiability_status": self.identifiability_status,
        }


@dataclass(frozen=True)
class ActionAuditRecord:
    action_id: str
    feasible: bool
    cost: float
    status: str
    rank: Optional[int]
    selected: bool
    robust_utility_per_cost: Optional[float]
    mean_utility_per_cost: Optional[float]
    information_gain_bits: Optional[float]
    scenario_disagreement: Optional[float]
    scenario_utilities: Tuple[ScenarioUtility, ...]
    metadata: Mapping[str, Any] = field(default_factory=dict)
    identifiability: str = "UNKNOWN"

    def to_dict(self) -> dict[str, Any]:
        return {
            "action_id": self.action_id,
            "feasible": self.feasible,
            "cost": self.cost,
            "status": self.status,
            "rank": self.rank,
            "selected": self.selected,
            "robust_utility_per_cost": self.robust_utility_per_cost,
            "mean_utility_per_cost": self.mean_utility_per_cost,
            "information_gain_bits": self.information_gain_bits,
            "scenario_disagreement": self.scenario_disagreement,
            "scenario_utilities": [
                utility.to_dict() for utility in self.scenario_utilities
            ],
            "metadata": dict(self.metadata),
            "identifiability": self.identifiability,
        }


@dataclass(frozen=True)
class DecisionUtilityReport:
    decision: str
    selected_action_id: Optional[str]
    threshold_utility_per_cost: float
    audit_records: Tuple[ActionAuditRecord, ...]
    truth_blind: bool = True
    selection_contract: str = "PRE_MEASUREMENT_TRUTH_BLIND"
    input_hash: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "selected_action_id": self.selected_action_id,
            "threshold_utility_per_cost": self.threshold_utility_per_cost,
            "audit_records": [record.to_dict() for record in self.audit_records],
            "truth_blind": self.truth_blind,
            "selection_contract": self.selection_contract,
            "input_hash": self.input_hash,
        }

    def to_json(self) -> str:
        """Return deterministic JSON for hashing or archival audit logs."""

        return json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))


def expected_information_gain(
    scenario: ScenarioBelief,
) -> tuple[float, float, dict[str, float], dict[str, float]]:
    """Compute entropy reduction for a declared scenario.

    Returns ``prior_entropy``, ``expected_posterior_entropy``,
    ``outcome_probabilities``, and ``posterior_entropy_by_outcome``.
    """

    prior_entropy = entropy_bits(scenario.prior)
    outcome_probs, posteriors = scenario.resolved_outcomes()
    posterior_entropies = {
        outcome: entropy_bits(posterior) for outcome, posterior in posteriors.items()
    }
    expected_posterior_entropy = sum(
        outcome_probs[outcome] * posterior_entropies[outcome]
        for outcome in outcome_probs
    )
    return prior_entropy, expected_posterior_entropy, outcome_probs, posterior_entropies


def evaluate_action(action: CandidateMeasurementAction) -> ActionAuditRecord:
    """Evaluate a single action without selecting it."""

    scenarios = action.resolved_scenarios()
    if not action.feasible:
        return ActionAuditRecord(
            action_id=action.action_id,
            feasible=False,
            cost=float(action.cost),
            status="INFEASIBLE",
            rank=None,
            selected=False,
            robust_utility_per_cost=None,
            mean_utility_per_cost=None,
            information_gain_bits=None,
            scenario_disagreement=None,
            scenario_utilities=(),
            metadata=action.metadata,
            identifiability="INFEASIBLE",
        )

    scenario_utilities = []
    for scenario_id in sorted(scenarios):
        scenario = scenarios[scenario_id]
        (
            prior_entropy,
            expected_posterior_entropy,
            outcome_probs,
            posterior_entropies,
        ) = expected_information_gain(scenario)
        information_gain = prior_entropy - expected_posterior_entropy
        scenario_utilities.append(
            ScenarioUtility(
                scenario_id=str(scenario_id),
                prior_entropy_bits=prior_entropy,
                expected_posterior_entropy_bits=expected_posterior_entropy,
                information_gain_bits=information_gain,
                utility_per_cost=information_gain / float(action.cost),
                outcome_probabilities=outcome_probs,
                posterior_entropy_by_outcome=posterior_entropies,
                identifiability_status="SCENARIO_ENVELOPE_NOT_CALIBRATED",
            )
        )

    utilities = [item.utility_per_cost for item in scenario_utilities]
    information_gains = [item.information_gain_bits for item in scenario_utilities]
    return ActionAuditRecord(
        action_id=action.action_id,
        feasible=True,
        cost=float(action.cost),
        status="EVALUATED",
        rank=None,
        selected=False,
        robust_utility_per_cost=min(utilities),
        mean_utility_per_cost=sum(utilities) / len(utilities),
        information_gain_bits=min(information_gains),
        scenario_disagreement=max(utilities) - min(utilities),
        scenario_utilities=tuple(scenario_utilities),
        metadata=action.metadata,
        identifiability=(
            "SCENARIO_DISAGREEMENT"
            if max(utilities) - min(utilities) > 0.0
            else "SCENARIOS_AGREE"
        ),
    )


def _validate_truth_blind_actions(actions: Sequence[CandidateMeasurementAction]) -> None:
    for action in actions:
        if _contains_forbidden_truth_key(action.metadata):
            raise ValueError(
                f"Action {action.action_id!r} contains a truth field in metadata; "
                "policy selection must be truth-blind."
            )


def _ranked_records(
    records: Sequence[ActionAuditRecord],
    *,
    selected_id: str | None,
    threshold: float,
    input_hash: str,
) -> DecisionUtilityReport:
    evaluated = [
        record for record in records
        if record.status == "EVALUATED" and record.robust_utility_per_cost is not None
    ]
    ranked_ids = {
        record.action_id: rank
        for rank, record in enumerate(
            sorted(
                evaluated,
                key=lambda item: (
                    -_numeric_or_negative_infinity(item.robust_utility_per_cost),
                    -_numeric_or_negative_infinity(item.information_gain_bits),
                    item.cost,
                    item.action_id,
                ),
            ),
            start=1,
        )
    }
    if selected_id is not None and selected_id not in ranked_ids:
        selected_id = None
    ranked_records = tuple(
        ActionAuditRecord(
            action_id=record.action_id,
            feasible=record.feasible,
            cost=record.cost,
            status=record.status,
            rank=ranked_ids.get(record.action_id),
            selected=record.action_id == selected_id,
            robust_utility_per_cost=record.robust_utility_per_cost,
            mean_utility_per_cost=record.mean_utility_per_cost,
            information_gain_bits=record.information_gain_bits,
            scenario_disagreement=record.scenario_disagreement,
            scenario_utilities=record.scenario_utilities,
            metadata=record.metadata,
            identifiability=record.identifiability,
        )
        for record in records
    )
    return DecisionUtilityReport(
        decision="MEASURE" if selected_id is not None else "ABSTAIN",
        selected_action_id=selected_id,
        threshold_utility_per_cost=threshold,
        audit_records=ranked_records,
        truth_blind=True,
        selection_contract="PRE_MEASUREMENT_TRUTH_BLIND",
        input_hash=input_hash,
    )


def select_next_measurement(
    actions: Iterable[CandidateMeasurementAction],
    *,
    min_utility_per_cost: float,
) -> DecisionUtilityReport:
    """Rank actions and select the best robust utility above a fixed threshold.

    The threshold is a predeclared stopping/abstention rule.  If every feasible
    action falls below it, the report returns ``decision="ABSTAIN"`` and keeps
    all serialisable audit records.
    """

    if not math.isfinite(float(min_utility_per_cost)):
        raise ValueError("min_utility_per_cost must be finite")

    action_list = list(actions)
    _validate_truth_blind_actions(action_list)
    seen: set[str] = set()
    for action in action_list:
        if action.action_id in seen:
            raise ValueError(f"duplicate action_id {action.action_id!r}")
        seen.add(action.action_id)
    action_list.sort(key=lambda action: action.action_id)
    records = [evaluate_action(action) for action in action_list]
    evaluated = [
        record for record in records
        if record.status == "EVALUATED" and record.robust_utility_per_cost is not None
    ]
    best = min(
        evaluated,
        key=lambda item: (
            -_numeric_or_negative_infinity(item.robust_utility_per_cost),
            -_numeric_or_negative_infinity(item.information_gain_bits),
            item.cost,
            item.action_id,
        ),
        default=None,
    )
    selected_id = (
        best.action_id
        if best is not None
        and best.robust_utility_per_cost is not None
        and best.robust_utility_per_cost >= float(min_utility_per_cost)
        else None
    )
    input_hash = _stable_hash([record.to_dict() for record in records])
    return _ranked_records(
        records,
        selected_id=selected_id,
        threshold=float(min_utility_per_cost),
        input_hash=input_hash,
    )


def select_random_measurement(
    actions: Iterable[CandidateMeasurementAction],
    *,
    seed: int = 0,
) -> DecisionUtilityReport:
    """Select one feasible action with a stable, truth-blind random baseline."""

    action_list = sorted(list(actions), key=lambda action: action.action_id)
    _validate_truth_blind_actions(action_list)
    ids = [action.action_id for action in action_list]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate action_id in random policy input")
    records = [evaluate_action(action) for action in action_list]
    eligible = [record for record in records if record.status == "EVALUATED"]
    input_hash = _stable_hash([record.to_dict() for record in records])
    if not eligible:
        return _ranked_records(
            records, selected_id=None, threshold=0.0, input_hash=input_hash
        )
    digest = _stable_hash({"seed": int(seed), "action_ids": ids})
    selected = eligible[int(digest[:16], 16) % len(eligible)].action_id
    return _ranked_records(
        records, selected_id=selected, threshold=-1.0e308, input_hash=input_hash
    )


def select_specialist_measurement(
    actions: Iterable[CandidateMeasurementAction],
    scores: Mapping[str, float],
) -> DecisionUtilityReport:
    """Select using a declared specialist score without seeing synthetic truth."""

    action_list = sorted(list(actions), key=lambda action: action.action_id)
    _validate_truth_blind_actions(action_list)
    records = [evaluate_action(action) for action in action_list]
    eligible = [record for record in records if record.status == "EVALUATED"]
    if not isinstance(scores, Mapping):
        raise TypeError("specialist scores must be a mapping")
    finite_scores = {str(key): float(value) for key, value in scores.items()}
    if any(not math.isfinite(value) for value in finite_scores.values()):
        raise ValueError("specialist scores must be finite")
    missing = [record.action_id for record in eligible if record.action_id not in finite_scores]
    if missing:
        raise ValueError(
            "specialist scores are incomplete for eligible actions: "
            + ", ".join(sorted(missing))
        )
    selected = max(
        eligible,
        key=lambda record: (
            finite_scores[record.action_id],
            -record.cost,
            record.action_id,
        ),
        default=None,
    )
    input_hash = _stable_hash(
        {"records": [record.to_dict() for record in records], "specialist_scores": finite_scores}
    )
    return _ranked_records(
        records,
        selected_id=selected.action_id if selected is not None else None,
        threshold=-1.0e308,
        input_hash=input_hash,
    )


def select_declared_utility_measurement(
    actions: Iterable[CandidateMeasurementAction],
    scores: Mapping[str, float],
    *,
    min_utility_per_cost: float = -float("inf"),
) -> DecisionUtilityReport:
    """Select using a declared truth-blind expected-value score.

    This is distinct from entropy-based information gain and from the
    specialist comparator.  The score is supplied by a pre-measurement model
    and is recorded in the same action audit, so it cannot consume realised
    outcomes or hidden labels.
    """

    if not math.isfinite(float(min_utility_per_cost)) and not math.isinf(
        float(min_utility_per_cost)
    ):
        raise ValueError("min_utility_per_cost must be finite or infinite")
    action_list = sorted(list(actions), key=lambda action: action.action_id)
    _validate_truth_blind_actions(action_list)
    ids = [action.action_id for action in action_list]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate action_id in declared-utility input")
    if not isinstance(scores, Mapping):
        raise TypeError("declared utility scores must be a mapping")
    finite_scores = {str(key): float(value) for key, value in scores.items()}
    if any(not math.isfinite(value) for value in finite_scores.values()):
        raise ValueError("declared utility scores must be finite")
    records = [evaluate_action(action) for action in action_list]
    evaluated = [record for record in records if record.status == "EVALUATED"]
    missing = [record.action_id for record in evaluated if record.action_id not in finite_scores]
    if missing:
        raise ValueError(
            "declared utility scores are incomplete for eligible actions: "
            + ", ".join(sorted(missing))
        )
    scored_records = tuple(
        ActionAuditRecord(
            action_id=record.action_id,
            feasible=record.feasible,
            cost=record.cost,
            status=record.status,
            rank=None,
            selected=False,
            robust_utility_per_cost=(
                finite_scores[record.action_id]
                if record.status == "EVALUATED"
                else None
            ),
            mean_utility_per_cost=(
                finite_scores[record.action_id]
                if record.status == "EVALUATED"
                else None
            ),
            information_gain_bits=record.information_gain_bits,
            scenario_disagreement=record.scenario_disagreement,
            scenario_utilities=record.scenario_utilities,
            metadata=record.metadata,
            identifiability="DECLARED_EXPECTED_UTILITY",
        )
        for record in records
    )
    best = max(
        (
            record
            for record in scored_records
            if record.status == "EVALUATED"
            and record.robust_utility_per_cost is not None
        ),
        key=lambda record: (
            float(record.robust_utility_per_cost),
            -record.cost,
            record.action_id,
        ),
        default=None,
    )
    selected_id = (
        best.action_id
        if best is not None
        and best.robust_utility_per_cost is not None
        and best.robust_utility_per_cost >= float(min_utility_per_cost)
        else None
    )
    input_hash = _stable_hash(
        {"records": [record.to_dict() for record in scored_records], "scores": finite_scores}
    )
    return _ranked_records(
        scored_records,
        selected_id=selected_id,
        threshold=float(min_utility_per_cost),
        input_hash=input_hash,
    )


@dataclass(frozen=True)
class ProspectiveMeasurementCase:
    """Synthetic case; truth is consumed only during post-policy scoring.

    ``true_state`` supports the original one-hypothesis contract.  For a
    candidate universe containing several possible targets,
    ``true_state_by_action`` records the hidden state of each action instead.
    This keeps target selection honest: a selector cannot see the state map,
    while the scorer can evaluate the state of the action that was actually
    chosen.
    """

    case_id: str
    actions: Tuple[CandidateMeasurementAction, ...]
    benefit_by_action_and_state: Mapping[str, Mapping[str, float]]
    true_state: str = ""
    true_state_by_action: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        case_id = str(self.case_id).strip()
        true_state = str(self.true_state).strip()
        state_by_action = {
            str(action_id): str(state).strip()
            for action_id, state in self.true_state_by_action.items()
        }
        if not case_id or (not true_state and not state_by_action):
            raise ValueError(
                "Prospective cases require case_id and either true_state or "
                "true_state_by_action."
            )
        actions = tuple(sorted(self.actions, key=lambda action: action.action_id))
        ids = [action.action_id for action in actions]
        if len(ids) != len(set(ids)):
            raise ValueError("Prospective case action IDs must be unique.")
        _validate_truth_blind_actions(actions)
        benefits: dict[str, dict[str, float]] = {}
        for action_id, by_state in self.benefit_by_action_and_state.items():
            action_key = str(action_id)
            if not isinstance(by_state, Mapping):
                raise TypeError(f"Benefits for {action_key!r} must be a mapping.")
            benefits[action_key] = {}
            for state, raw_value in by_state.items():
                value = float(raw_value)
                if not math.isfinite(value):
                    raise ValueError("Prospective benefits must be finite.")
                benefits[action_key][str(state)] = value
        available_states = {
            state for by_state in benefits.values() for state in by_state
        }
        if true_state and true_state not in available_states:
            raise ValueError("true_state must have a truth-released benefit record.")
        if state_by_action:
            missing_states = sorted(
                set(state_by_action.values()) - available_states
            )
            if missing_states:
                raise ValueError(
                    "true_state_by_action contains states absent from the benefit "
                    f"table: {missing_states}"
                )
            missing_actions = sorted(set(ids) - set(state_by_action))
            if missing_actions:
                raise ValueError(
                    "true_state_by_action must cover every candidate action: "
                    + ", ".join(missing_actions)
                )
        object.__setattr__(self, "case_id", case_id)
        object.__setattr__(self, "true_state", true_state)
        object.__setattr__(self, "actions", actions)
        object.__setattr__(self, "benefit_by_action_and_state", benefits)
        object.__setattr__(self, "true_state_by_action", state_by_action)

    def state_for_action(self, action_id: str) -> str:
        """Return the hidden state for one selected action after selection."""

        return self.true_state_by_action.get(str(action_id), self.true_state)


PolicySelector = Callable[[Tuple[CandidateMeasurementAction, ...]], DecisionUtilityReport]


@dataclass(frozen=True)
class ProspectivePolicy:
    """Named policy whose selector receives pre-measurement actions only."""

    policy_id: str
    selector: PolicySelector
    truth_blind: bool = True
    scoring_mode: str = "selected_action"

    def __post_init__(self) -> None:
        policy_id = str(self.policy_id).strip()
        if not policy_id:
            raise ValueError("Prospective policies require policy_id.")
        if not callable(self.selector):
            raise TypeError("Prospective policy selector must be callable.")
        scoring_mode = str(self.scoring_mode).strip()
        if scoring_mode not in {"selected_action", "uniform_action_expectation"}:
            raise ValueError(
                "Prospective policy scoring_mode must be selected_action or "
                "uniform_action_expectation."
            )
        object.__setattr__(self, "policy_id", policy_id)
        object.__setattr__(self, "scoring_mode", scoring_mode)


@dataclass(frozen=True)
class ProspectivePolicyScore:
    policy_id: str
    status: str
    case_count: int
    scored_case_count: int
    selection_coverage: float
    abstention_rate: float
    mean_cost_adjusted_utility: float | None
    mean_regret: float | None
    mean_downside_risk: float | None
    negative_utility_rate: float | None
    worst_case_cost_adjusted_utility: float | None
    evaluation_mode: str = "selected_action"

    def to_dict(self) -> dict[str, Any]:
        return {
            "policy_id": self.policy_id,
            "status": self.status,
            "case_count": self.case_count,
            "scored_case_count": self.scored_case_count,
            "selection_coverage": self.selection_coverage,
            "abstention_rate": self.abstention_rate,
            "mean_cost_adjusted_utility": self.mean_cost_adjusted_utility,
            "mean_regret": self.mean_regret,
            "mean_downside_risk": self.mean_downside_risk,
            "negative_utility_rate": self.negative_utility_rate,
            "worst_case_cost_adjusted_utility": self.worst_case_cost_adjusted_utility,
            "evaluation_mode": self.evaluation_mode,
        }


@dataclass(frozen=True)
class ProspectivePolicyBenchmark:
    status: str
    claim_status: str
    reason: str
    case_count: int
    cost_penalty: float
    policies: Mapping[str, ProspectivePolicyScore]
    pairwise: Mapping[str, Mapping[str, Any]]
    calibration_status: str = "NOT_PROVIDED"
    calibration_sufficient: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "claim_status": self.claim_status,
            "reason": self.reason,
            "case_count": self.case_count,
            "cost_penalty": self.cost_penalty,
            "policies": {key: value.to_dict() for key, value in self.policies.items()},
            "pairwise": {key: dict(value) for key, value in self.pairwise.items()},
            "calibration_status": self.calibration_status,
            "calibration_sufficient": self.calibration_sufficient,
        }


def evaluate_prospective_policies(
    cases: Iterable[ProspectiveMeasurementCase],
    policies: Iterable[ProspectivePolicy],
    *,
    cost_penalty: float = 1.0,
    required_policy_ids: Sequence[str] = ("random", "specialist"),
    calibration_sufficient: bool = False,
) -> ProspectivePolicyBenchmark:
    """Score truth-blind choices after synthetic truth is released.

    Selectors receive only candidate actions.  The truth and independent
    benefit table are accessed only after every selector has returned.  Scores
    are therefore available for held-out comparison, but this helper never
    turns them into an integrated superiority claim.

    A policy with ``scoring_mode="uniform_action_expectation"`` is scored as
    the exact expected utility of a uniformly sampled feasible action.  This
    avoids treating one random seed as a stable comparator while preserving a
    truth-blind common-action baseline.
    """

    penalty = float(cost_penalty)
    if not math.isfinite(penalty) or penalty < 0.0:
        raise ValueError("cost_penalty must be finite and non-negative")
    case_list = sorted(tuple(cases), key=lambda case: case.case_id)
    policy_list = sorted(tuple(policies), key=lambda policy: policy.policy_id)
    policy_ids = [policy.policy_id for policy in policy_list]
    if len(policy_ids) != len(set(policy_ids)):
        raise ValueError("Prospective policy IDs must be unique.")
    if not case_list:
        return ProspectivePolicyBenchmark(
            status="ABSTAIN_NO_CASES",
            claim_status="ABSTAIN",
            reason="No prospective synthetic cases were supplied.",
            case_count=0,
            cost_penalty=penalty,
            policies={},
            pairwise={},
            calibration_status="NOT_PROVIDED",
            calibration_sufficient=False,
        )

    required = tuple(str(policy_id) for policy_id in required_policy_ids)
    missing_required = sorted(set(required) - set(policy_ids))
    if missing_required:
        selection_status = "ABSTAIN_REQUIRED_POLICY_MISSING"
    elif any(not policy.truth_blind for policy in policy_list):
        selection_status = "ABSTAIN_POLICY_NOT_TRUTH_BLIND"
    else:
        selection_status = "SCORED"

    selected_by_policy: dict[str, dict[str, str | None]] = {
        policy.policy_id: {} for policy in policy_list
    }
    uniform_expectation_marker = "__uniform_action_expectation__"
    for case in case_list:
        for policy in policy_list:
            if policy.scoring_mode == "uniform_action_expectation":
                selected_by_policy[policy.policy_id][case.case_id] = (
                    uniform_expectation_marker
                )
                continue
            if not policy.truth_blind:
                selected_by_policy[policy.policy_id][case.case_id] = None
                continue
            report = policy.selector(case.actions)
            if not isinstance(report, DecisionUtilityReport) or not report.truth_blind:
                selected_by_policy[policy.policy_id][case.case_id] = None
                continue
            selected = report.selected_action_id if report.decision == "MEASURE" else None
            feasible_ids = {action.action_id for action in case.actions if action.feasible}
            selected_by_policy[policy.policy_id][case.case_id] = (
                selected if selected in feasible_ids else None
            )

    case_metrics: dict[str, dict[str, tuple[float, float, float]]] = {
        policy.policy_id: {} for policy in policy_list
    }
    case_negative_rates: dict[str, dict[str, float]] = {
        policy.policy_id: {} for policy in policy_list
    }
    incomplete = False
    for case in case_list:
        feasible_actions = [action for action in case.actions if action.feasible]
        oracle_values = []
        for action in feasible_actions:
            truth = case.state_for_action(action.action_id)
            raw_benefit = case.benefit_by_action_and_state.get(action.action_id, {}).get(truth)
            if raw_benefit is not None:
                oracle_values.append(float(raw_benefit) - penalty * float(action.cost))
        for policy in policy_list:
            selected_id = selected_by_policy[policy.policy_id][case.case_id]
            if policy.scoring_mode == "uniform_action_expectation":
                action_utilities: list[float] = []
                for action in feasible_actions:
                    truth = case.state_for_action(action.action_id)
                    raw_benefit = case.benefit_by_action_and_state.get(
                        action.action_id, {}
                    ).get(truth)
                    if raw_benefit is None:
                        action_utilities = []
                        break
                    action_utilities.append(
                        float(raw_benefit) - penalty * float(action.cost)
                    )
                if not action_utilities:
                    incomplete = True
                    continue
                expected_utility = sum(action_utilities) / len(action_utilities)
                oracle_utility = max(action_utilities)
                case_metrics[policy.policy_id][case.case_id] = (
                    expected_utility,
                    max(0.0, oracle_utility - expected_utility),
                    sum(max(0.0, -value) for value in action_utilities)
                    / len(action_utilities),
                )
                case_negative_rates[policy.policy_id][case.case_id] = sum(
                    value < 0.0 for value in action_utilities
                ) / len(action_utilities)
                continue
            if selected_id is None:
                incomplete = True
                continue
            action = next(action for action in case.actions if action.action_id == selected_id)
            truth = case.state_for_action(action.action_id)
            raw_benefit = case.benefit_by_action_and_state.get(action.action_id, {}).get(truth)
            if raw_benefit is None or not oracle_values:
                incomplete = True
                continue
            utility = float(raw_benefit) - penalty * float(action.cost)
            regret = max(0.0, max(oracle_values) - utility)
            downside = max(0.0, -utility)
            case_metrics[policy.policy_id][case.case_id] = (utility, regret, downside)
            case_negative_rates[policy.policy_id][case.case_id] = float(utility < 0.0)

    scores: dict[str, ProspectivePolicyScore] = {}
    for policy in policy_list:
        metrics = list(case_metrics[policy.policy_id].values())
        scored_count = len(metrics)
        scores[policy.policy_id] = ProspectivePolicyScore(
            policy_id=policy.policy_id,
            status=(
                "SCORED" if scored_count == len(case_list)
                else "ABSTAIN_INCOMPLETE_OUTCOMES" if scored_count
                else "ABSTAIN_NO_SELECTED_ACTION"
            ),
            case_count=len(case_list),
            scored_case_count=scored_count,
            selection_coverage=(
                1.0
                if policy.scoring_mode == "uniform_action_expectation"
                and scored_count == len(case_list)
                else scored_count / len(case_list)
            ),
            abstention_rate=(
                0.0
                if policy.scoring_mode == "uniform_action_expectation"
                and scored_count == len(case_list)
                else 1.0 - scored_count / len(case_list)
            ),
            mean_cost_adjusted_utility=(sum(item[0] for item in metrics) / scored_count if metrics else None),
            mean_regret=(sum(item[1] for item in metrics) / scored_count if metrics else None),
            mean_downside_risk=(sum(item[2] for item in metrics) / scored_count if metrics else None),
            negative_utility_rate=(
                sum(case_negative_rates[policy.policy_id].values()) / scored_count
                if metrics
                and len(case_negative_rates[policy.policy_id]) == scored_count
                else None
            ),
            worst_case_cost_adjusted_utility=(min(item[0] for item in metrics) if metrics else None),
            evaluation_mode=policy.scoring_mode,
        )

    pairwise: dict[str, dict[str, Any]] = {}
    for left in policy_ids:
        for right in policy_ids:
            if left >= right:
                continue
            common = sorted(set(case_metrics[left]) & set(case_metrics[right]))
            if not common:
                continue
            utility_deltas = [case_metrics[left][case_id][0] - case_metrics[right][case_id][0] for case_id in common]
            regret_deltas = [case_metrics[left][case_id][1] - case_metrics[right][case_id][1] for case_id in common]
            ci_low, ci_high, bootstrap_seed = _paired_bootstrap_mean_interval(
                utility_deltas,
                seed_material={
                    "left": left,
                    "right": right,
                    "case_ids": common,
                    "cost_penalty": penalty,
                },
            )
            pairwise[f"{left}_vs_{right}"] = {
                "left_policy": left,
                "right_policy": right,
                "paired_case_count": len(common),
                "mean_cost_adjusted_utility_delta": sum(utility_deltas) / len(common),
                "mean_regret_delta": sum(regret_deltas) / len(common),
                "left_utility_win_rate": sum(delta > 0.0 for delta in utility_deltas) / len(common),
                "paired_delta_ci_level": 0.95,
                "paired_delta_ci_low": ci_low,
                "paired_delta_ci_high": ci_high,
                "paired_bootstrap_replicates": 2000 if bootstrap_seed is not None else 0,
                "paired_bootstrap_seed": bootstrap_seed,
                "paired_uncertainty_available": ci_low is not None and ci_high is not None,
            }

    if selection_status == "SCORED" and not incomplete:
        status = "SCORED"
        if calibration_sufficient:
            reason = "Held-out policy scores are complete; no superiority gate is applied here."
            claim_status = "ABSTAIN_NO_SUPERIORITY_GATE"
        else:
            reason = "Held-out outcomes are complete, but calibration sufficiency was not declared."
            claim_status = "ABSTAIN_CALIBRATION_INSUFFICIENT"
    elif selection_status == "ABSTAIN_REQUIRED_POLICY_MISSING":
        status = selection_status
        reason = "Required random and specialist policy comparators are missing."
        claim_status = "ABSTAIN"
    elif selection_status == "ABSTAIN_POLICY_NOT_TRUTH_BLIND":
        status = selection_status
        reason = "At least one policy was not declared truth-blind."
        claim_status = "ABSTAIN"
    else:
        status = "ABSTAIN_INCOMPLETE_OUTCOMES"
        reason = "At least one policy/case lacks a post-measurement outcome score."
        claim_status = "ABSTAIN"
    return ProspectivePolicyBenchmark(
        status=status,
        claim_status=claim_status,
        reason=reason,
        case_count=len(case_list),
        cost_penalty=penalty,
        policies=scores,
        pairwise=pairwise,
        calibration_status=("SUFFICIENT_DECLARED" if calibration_sufficient else "NOT_PROVIDED"),
        calibration_sufficient=bool(calibration_sufficient),
    )
