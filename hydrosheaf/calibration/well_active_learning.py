"""Bayesian active learning and experimental design at the well/campaign level.

Connects the discrete topology posterior ensemble to the Bayesian acquisition
engine, making the decision unit a well measurement action:
    ``well_id + measurement_type + sampling_time``
with predicted observations under each plausible topology, measurement error,
cost, accessibility, and feasibility.

Key capabilities:
  1. Consumes a probability-bearing ensemble of posterior graphs.
  2. Evaluates robust expected topology-information gain and validation-decision
     risk reduction targeting False Positives, False Negatives, and
     Selected-Unlabeled candidate edges.
  3. Selects non-redundant batches at the well level, enforcing well-level
     deduplication and discounting shared travel/mobilization costs.
  4. Enforces an explicit abstention gate when MCMC diagnostics (R-hat, ESS),
     validation labels, or predictive measurement models are inadequate.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple, Union

import numpy as np

from ..graph.types import Edge
from ..inference.topology_posterior import extract_topology_hypotheses
from ..log import get_logger
from .bayesian_active_learning import (
    AcquisitionConfig,
    MeasurementOption,
    PredictiveScenario,
    rank_measurement_options,
    shannon_entropy,
    _prior_decision_risk,
    _robust_joint_eig,
    _scenario_eigs,
)

logger = get_logger("calibration.well_active_learning")

DEFAULT_MEASUREMENT_COSTS: Dict[str, float] = {
    "head_monitoring": 1.0,
    "major_ions": 2.5,
    "isotopes_d18O_d2H": 3.0,
    "age_tracer": 6.0,
    "connectivity_tracer": 10.0,
}

DEFAULT_TRAVEL_COST_PER_WELL: float = 3.0

SUPPORTED_MEASUREMENT_TYPES = tuple(DEFAULT_MEASUREMENT_COSTS.keys())

# A prediction provider receives one explicit action/hypothesis context and
# returns ``(mean, standard_deviation)`` for that observation.  Keeping this
# contract callable-based lets a field/transport model supply predictions
# without making the acquisition engine depend on a particular simulator.
PredictiveModel = Callable[[Mapping[str, Any]], Tuple[float, float]]


def _action_id(well_id: str, measurement_type: str, sampling_time: float) -> str:
    return f"{well_id}@{measurement_type}@t{float(sampling_time):.1f}"


def _validate_prediction(value: Any) -> Tuple[float, float]:
    """Validate one predictive-model response before it enters EIG scoring."""
    if isinstance(value, Mapping):
        mean = value.get("mean")
        sd = value.get("sd", value.get("standard_deviation"))
    else:
        try:
            mean, sd = value
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Predictive model must return (mean, sd) or a mapping with "
                "'mean' and 'sd'."
            ) from exc
    if mean is None or sd is None:
        raise ValueError("Predictive model response is missing mean or sd.")
    mean_f = float(mean)
    sd_f = float(sd)
    if not math.isfinite(mean_f) or not math.isfinite(sd_f) or sd_f <= 0.0:
        raise ValueError("Predictive model mean must be finite and sd must be positive.")
    return mean_f, sd_f


def topology_contrast_surrogate(context: Mapping[str, Any]) -> Tuple[float, float]:
    """Return a clearly labelled topology-contrast surrogate prediction.

    This is retained for synthetic/unit-test use and explicit exploratory
    campaigns. Production campaigns should pass a calibrated predictive model
    or a predictive-scenario file instead of silently relying on these values.
    """
    well_id = str(context.get("well_id", ""))
    measurement_type = str(context.get("measurement_type", ""))
    sampling_time = float(context.get("sampling_time", 0.0))
    active_edges = set(str(e) for e in context.get("active_edge_ids", ()))
    candidate_edges = context.get("candidate_edges", ())
    incoming = [
        edge for edge in candidate_edges
        if str(getattr(edge, "v", "")).strip() == well_id
    ]
    outgoing = [
        edge for edge in candidate_edges
        if str(getattr(edge, "u", "")).strip() == well_id
    ]
    active_in = [edge for edge in incoming if str(edge.edge_id) in active_edges]
    active_out = [edge for edge in outgoing if str(edge.edge_id) in active_edges]

    if measurement_type == "connectivity_tracer":
        if active_in:
            mean = math.exp(-0.5 * ((sampling_time - 30.0) / 10.0) ** 2)
        else:
            mean = 0.0
        return float(mean), 0.08
    if measurement_type == "head_monitoring":
        return (97.5 if active_in else (101.5 if active_out else 100.0)), 0.50
    if measurement_type == "isotopes_d18O_d2H":
        return (-8.5 if active_in else -7.0), 0.25
    if measurement_type == "age_tracer":
        return (35.0 if active_in else 50.0), 2.0
    if measurement_type == "major_ions":
        return (14.0 if active_in else 10.0), 0.8
    return (1.0 if (active_in or active_out) else 0.0), 0.20


def load_predictive_scenarios_file(path: str) -> PredictiveModel:
    """Load an explicit action/hypothesis predictive model from JSON.

    The accepted shape is either ``{"actions": {...}}`` or the action mapping
    directly. Each action maps hypothesis IDs to ``{"mean": ..., "sd": ...}``.
    """
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if isinstance(payload, Mapping) and isinstance(payload.get("actions"), Mapping):
        payload = payload["actions"]
    if not isinstance(payload, Mapping):
        raise ValueError("Predictive scenario file must contain an action mapping.")

    def _file_model(context: Mapping[str, Any]) -> Tuple[float, float]:
        action_id = str(
            context.get("action_id")
            or _action_id(
                str(context.get("well_id", "")),
                str(context.get("measurement_type", "")),
                float(context.get("sampling_time", 0.0)),
            )
        )
        entry = payload.get(action_id)
        if entry is None:
            # Permit a time-independent key for campaign planners that use one
            # predictive model for all sampling times.
            entry = payload.get(
                f"{context.get('well_id')}@{context.get('measurement_type')}"
            )
        if not isinstance(entry, Mapping):
            raise KeyError(f"No predictive scenarios found for action {action_id!r}.")
        hypothesis_id = str(context.get("hypothesis_id", ""))
        prediction = entry.get(hypothesis_id, entry.get("default"))
        if prediction is None and "mean" in entry:
            prediction = entry
        if prediction is None:
            raise KeyError(
                f"No predictive scenario found for action {action_id!r}, "
                f"hypothesis {hypothesis_id!r}."
            )
        return _validate_prediction(prediction)

    return _file_model


@dataclass(frozen=True)
class HypothesisScenario:
    """Predicted observation under one topology hypothesis."""

    hypothesis_id: str
    mean: float
    sd: float
    weight: float = 1.0

    @property
    def standard_deviation(self) -> float:
        return self.sd

    @property
    def name(self) -> str:
        return self.hypothesis_id


@dataclass(frozen=True)
class WellAction:
    """A physical measurement action at a designated well and sampling time."""

    well_id: str
    measurement_type: str
    sampling_time: float = 0.0
    base_cost: float = 2.0
    travel_cost: float = 3.0
    accessibility: float = 1.0
    feasible: bool = True
    connected_edge_ids: Tuple[str, ...] = ()
    target_edges: Tuple[str, ...] = ()
    scenarios: Tuple[HypothesisScenario, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.target_edges and not self.connected_edge_ids:
            object.__setattr__(self, "connected_edge_ids", tuple(self.target_edges))
        elif self.connected_edge_ids and not self.target_edges:
            object.__setattr__(self, "target_edges", tuple(self.connected_edge_ids))

    @property
    def standalone_cost(self) -> float:
        """Effective standalone cost adjusted for accessibility."""
        acc = max(0.05, min(1.0, float(self.accessibility)))
        return float((self.base_cost + self.travel_cost) / acc)

    @property
    def incremental_test_cost(self) -> float:
        """Incremental cost when travel mobilization is already paid."""
        acc = max(0.05, min(1.0, float(self.accessibility)))
        return float(self.base_cost / acc)

    @property
    def action_id(self) -> str:
        return f"{self.well_id}@{self.measurement_type}@t{self.sampling_time:.1f}"

    def __str__(self) -> str:
        return f"WellAction({self.action_id}, standalone_cost={self.standalone_cost:.2f})"

    def __repr__(self) -> str:
        return (
            f"WellAction(well_id={self.well_id!r}, type={self.measurement_type!r}, "
            f"time={self.sampling_time}, cost={self.standalone_cost:.2f})"
        )


@dataclass(frozen=True)
class CampaignConfig:
    """Configuration for well-level Bayesian experimental design."""

    acquisition_config: AcquisitionConfig = field(
        default_factory=lambda: AcquisitionConfig(decision_weight=0.5)
    )
    batch_size: int = 3
    budget: Optional[float] = None
    max_actions_per_well: int = 1
    minimum_robust_eig: float = 1.0e-4
    max_r_hat: float = 1.10
    max_mcmc_r_hat: float = 1.10
    min_ess: float = 20.0
    min_mcmc_ess: float = 20.0
    min_entropy: float = 1.0e-4
    require_mcmc_diagnostics: bool = True
    require_independent_validation: bool = True
    allow_surrogate_predictive_model: bool = True
    predictive_model: Optional[PredictiveModel] = None
    predictive_model_name: str = "topology_contrast_surrogate"
    travel_cost_per_well: float = 3.0
    measurement_costs: Mapping[str, float] = field(
        default_factory=lambda: dict(DEFAULT_MEASUREMENT_COSTS)
    )
    measurement_types: Sequence[str] = SUPPORTED_MEASUREMENT_TYPES
    sampling_times: Sequence[float] = (0.0,)
    cost_exponent: float = 1.0

    def __post_init__(self):
        if self.max_mcmc_r_hat is not None and (
            self.max_r_hat == 1.10 or self.max_mcmc_r_hat != 1.10
        ):
            object.__setattr__(self, "max_r_hat", float(self.max_mcmc_r_hat))
        if self.min_mcmc_ess is not None and (
            self.min_ess == 20.0 or self.min_mcmc_ess != 20.0
        ):
            object.__setattr__(self, "min_ess", float(self.min_mcmc_ess))
        if self.cost_exponent != 1.0 and self.acquisition_config.cost_exponent == 1.0:
            object.__setattr__(
                self,
                "acquisition_config",
                AcquisitionConfig(
                    quadrature_order=self.acquisition_config.quadrature_order,
                    robustness_weight=self.acquisition_config.robustness_weight,
                    decision_weight=self.acquisition_config.decision_weight,
                    cost_exponent=float(self.cost_exponent),
                    minimum_expected_information_gain=self.acquisition_config.minimum_expected_information_gain,
                    batch_qmc_samples=self.acquisition_config.batch_qmc_samples,
                    random_seed=self.acquisition_config.random_seed,
                ),
            )
        if self.predictive_model is not None and self.predictive_model_name == "topology_contrast_surrogate":
            object.__setattr__(self, "predictive_model_name", "configured_callback")

    def to_dict(self) -> Dict[str, Any]:
        return {
            "batch_size": self.batch_size,
            "budget": self.budget,
            "max_actions_per_well": self.max_actions_per_well,
            "minimum_robust_eig": self.minimum_robust_eig,
            "max_r_hat": self.max_r_hat,
            "min_ess": self.min_ess,
            "min_entropy": self.min_entropy,
            "require_mcmc_diagnostics": self.require_mcmc_diagnostics,
            "require_independent_validation": self.require_independent_validation,
            "travel_cost_per_well": self.travel_cost_per_well,
            "measurement_costs": dict(self.measurement_costs),
            "measurement_types": list(self.measurement_types),
            "sampling_times": list(self.sampling_times),
            "cost_exponent": self.cost_exponent,
            "acquisition_config": {
                "quadrature_order": self.acquisition_config.quadrature_order,
                "robustness_weight": self.acquisition_config.robustness_weight,
                "decision_weight": self.acquisition_config.decision_weight,
                "cost_exponent": self.acquisition_config.cost_exponent,
                "minimum_expected_information_gain": self.acquisition_config.minimum_expected_information_gain,
                "batch_qmc_samples": self.acquisition_config.batch_qmc_samples,
                "random_seed": self.acquisition_config.random_seed,
            },
            "allow_surrogate_predictive_model": self.allow_surrogate_predictive_model,
            "predictive_model_configured": self.predictive_model is not None,
            "predictive_model_name": self.predictive_model_name,
        }

    @classmethod
    def from_mapping(
        cls,
        settings: Optional[Mapping[str, Any]] = None,
        *,
        predictive_model: Optional[PredictiveModel] = None,
        allow_surrogate_default: bool = False,
    ) -> "CampaignConfig":
        """Build a campaign configuration from calibration ``model`` settings."""
        values = dict(settings or {})

        def _first(*keys: str, default: Any = None) -> Any:
            for key in keys:
                if key in values:
                    return values[key]
            return default

        measurement_types = _first(
            "active_learning_measurement_types",
            "measurement_types",
            default=SUPPORTED_MEASUREMENT_TYPES,
        )
        if isinstance(measurement_types, str):
            measurement_types = (measurement_types,)
        sampling_times = _first(
            "active_learning_sampling_times",
            "sampling_times",
            default=(0.0,),
        )
        if isinstance(sampling_times, (int, float)):
            sampling_times = (float(sampling_times),)

        costs = dict(DEFAULT_MEASUREMENT_COSTS)
        configured_costs = _first(
            "active_learning_measurement_costs",
            "measurement_costs",
            default={},
        )
        if isinstance(configured_costs, Mapping):
            costs.update({str(k): float(v) for k, v in configured_costs.items()})

        acquisition = AcquisitionConfig(
            quadrature_order=int(_first("active_learning_quadrature_order", default=21)),
            robustness_weight=float(_first("active_learning_robustness_weight", default=0.5)),
            decision_weight=float(_first("active_learning_decision_weight", default=0.5)),
            cost_exponent=float(_first("active_learning_cost_exponent", "cost_exponent", default=1.0)),
            minimum_expected_information_gain=float(
                _first("active_learning_minimum_eig", default=1.0e-4)
            ),
            batch_qmc_samples=int(_first("active_learning_batch_qmc_samples", default=2048)),
            random_seed=int(_first("active_learning_random_seed", default=20260728)),
        )

        budget = _first("active_learning_budget", "budget", default=None)
        return cls(
            acquisition_config=acquisition,
            batch_size=int(_first("active_learning_batch_size", "batch_size", default=3)),
            budget=None if budget is None else float(budget),
            max_actions_per_well=int(
                _first("active_learning_max_actions_per_well", "max_actions_per_well", default=1)
            ),
            minimum_robust_eig=float(
                _first("active_learning_minimum_robust_eig", default=1.0e-4)
            ),
            max_r_hat=float(_first("active_learning_max_r_hat", default=1.10)),
            max_mcmc_r_hat=float(_first("active_learning_max_r_hat", default=1.10)),
            min_ess=float(_first("active_learning_min_ess", default=20.0)),
            min_mcmc_ess=float(_first("active_learning_min_ess", default=20.0)),
            min_entropy=float(_first("active_learning_min_entropy", default=1.0e-4)),
            require_mcmc_diagnostics=bool(
                _first("active_learning_require_mcmc_diagnostics", default=True)
            ),
            require_independent_validation=bool(
                _first("active_learning_require_independent_validation", default=True)
            ),
            travel_cost_per_well=float(
                _first("active_learning_travel_cost_per_well", "travel_cost_per_well", default=3.0)
            ),
            measurement_costs=costs,
            measurement_types=tuple(str(v) for v in measurement_types),
            sampling_times=tuple(float(v) for v in sampling_times),
            cost_exponent=float(acquisition.cost_exponent),
            allow_surrogate_predictive_model=bool(
                _first(
                    "active_learning_allow_surrogate",
                    default=allow_surrogate_default,
                )
            ),
            predictive_model=predictive_model,
            predictive_model_name=(
                "configured_file" if predictive_model is not None
                else "topology_contrast_surrogate"
            ),
        )


# ── action generator & predictive scenarios ──────────────────────────


def build_well_actions(
    sample_nodes: Optional[Sequence[Any]] = None,
    candidate_edges: Sequence[Edge] = (),
    *,
    hypothesis_ids: Optional[Sequence[str]] = None,
    prior_probs: Optional[Sequence[float]] = None,
    edge_tuples: Optional[Sequence[Tuple[str, ...]]] = None,
    validation_decision_matrix: Optional[Any] = None,
    config: Optional[CampaignConfig] = None,
    accessibility_map: Optional[Mapping[str, float]] = None,
    feasibility_map: Optional[Mapping[str, bool]] = None,
    measurement_types: Optional[Sequence[str]] = None,
    sampling_times: Optional[Sequence[float]] = None,
    predictive_model: Optional[PredictiveModel] = None,
) -> List[WellAction]:
    """Generate candidate well measurement actions across nodes and measurement types.

    Parameters
    ----------
    sample_nodes : sequence of dict or Edge, optional
        Well/sample node dictionaries with ``site_id`` / ``well_id``.
        Can also be omitted if candidate_edges is passed.
    candidate_edges : sequence of Edge
        Candidate network edges.
    config : CampaignConfig, optional
        Campaign configuration.
    accessibility_map : mapping of str -> float, optional
        Well ID to accessibility factor in (0, 1].
    feasibility_map : mapping of str -> bool, optional
        Well ID to feasibility bool.
    measurement_types : sequence of str, optional
        Types of measurements to consider.
    """
    # Handle argument flexibility if candidate_edges passed first
    if sample_nodes is not None and len(sample_nodes) > 0 and isinstance(sample_nodes[0], Edge):
        candidate_edges = sample_nodes
        sample_nodes = None

    cfg = config or CampaignConfig()
    types = tuple(measurement_types or cfg.measurement_types or SUPPORTED_MEASUREMENT_TYPES)
    times = tuple(sampling_times or cfg.sampling_times or (0.0,))
    costs = cfg.measurement_costs
    travel_cost = cfg.travel_cost_per_well
    predictor = predictive_model or cfg.predictive_model
    predictor_name = cfg.predictive_model_name
    if predictor is None and cfg.allow_surrogate_predictive_model:
        predictor = topology_contrast_surrogate
        predictor_name = "topology_contrast_surrogate"

    # Map well IDs to incident edge IDs
    incident_edges: Dict[str, List[str]] = {}
    for edge in candidate_edges:
        u = str(edge.u or "").strip()
        v = str(edge.v or "").strip()
        if u:
            incident_edges.setdefault(u, []).append(edge.edge_id)
        if v:
            incident_edges.setdefault(v, []).append(edge.edge_id)

    # Collect well IDs from samples
    well_ids: List[str] = []
    seen_wells: Set[str] = set()
    sample_dict_lookup: Dict[str, Dict[str, Any]] = {}
    if sample_nodes:
        sample_iterable = (
            list(sample_nodes.values())
            if isinstance(sample_nodes, Mapping)
            else sample_nodes
        )
        for s in sample_iterable:
            if isinstance(s, dict):
                for key in ("site_id", "well_id", "node_id", "sample_id"):
                    val = s.get(key)
                    if val is not None and str(val).strip():
                        wid = str(val).strip()
                        sample_dict_lookup[wid] = s
                        if wid not in seen_wells:
                            seen_wells.add(wid)
                            well_ids.append(wid)
                        break

    # If samples didn't yield wells, extract from candidate edges
    if not well_ids:
        for edge in candidate_edges:
            for w in (str(edge.u or "").strip(), str(edge.v or "").strip()):
                if w and w not in seen_wells:
                    seen_wells.add(w)
                    well_ids.append(w)

    well_ids.sort()

    # Pre-parse hypotheses if provided
    hyp_ids = list(hypothesis_ids) if hypothesis_ids is not None else []
    hyp_edge_sets = [set(t) for t in edge_tuples] if edge_tuples is not None else []

    actions: List[WellAction] = []
    for wid in well_ids:
        # Feasibility check: if explicitly infeasible, skip
        if feasibility_map and wid in feasibility_map and not feasibility_map[wid]:
            continue

        connected = tuple(sorted(incident_edges.get(wid, [])))
        accessibility = 1.0
        if accessibility_map and wid in accessibility_map:
            accessibility = float(accessibility_map[wid])

        for m_type in types:
            base_c = float(costs.get(m_type, 2.0))
            for t_samp in times:
                t_val = float(t_samp)

                # Generate hypothesis scenarios if hypotheses provided
                scenarios_list: List[HypothesisScenario] = []
                prediction_error: Optional[str] = None
                if hyp_ids and hyp_edge_sets:
                    for k, h_id in enumerate(hyp_ids):
                        active_set = hyp_edge_sets[k] if k < len(hyp_edge_sets) else set()
                        if predictor is None:
                            prediction_error = (
                                "No predictive measurement model was supplied; "
                                "configure predictive_model or explicitly allow the surrogate."
                            )
                            scenarios_list = []
                            break
                        try:
                            mean_val, sd_val = _validate_prediction(
                                predictor(
                                    {
                                        "action_id": _action_id(wid, m_type, t_val),
                                        "well_id": wid,
                                        "measurement_type": m_type,
                                        "sampling_time": t_val,
                                        "hypothesis_id": str(h_id),
                                        "active_edge_ids": tuple(sorted(active_set)),
                                        "candidate_edges": tuple(candidate_edges),
                                        "sample": sample_dict_lookup.get(wid),
                                        "connected_edge_ids": connected,
                                    }
                                )
                            )
                        except Exception as exc:
                            prediction_error = str(exc)
                            scenarios_list = []
                            break
                        scenarios_list.append(
                            HypothesisScenario(
                                hypothesis_id=h_id,
                                mean=mean_val,
                                sd=sd_val,
                            )
                        )

                actions.append(
                    WellAction(
                        well_id=wid,
                        measurement_type=m_type,
                        sampling_time=t_val,
                        base_cost=base_c,
                        travel_cost=travel_cost,
                        accessibility=accessibility,
                        feasible=True,
                        connected_edge_ids=connected,
                        target_edges=connected,
                        scenarios=tuple(scenarios_list),
                        metadata={
                            "n_incident_edges": len(connected),
                            "incident_edges": list(connected),
                            "predictive_model": predictor_name if predictor is not None else "missing",
                            "predictive_model_error": prediction_error,
                        },
                    )
                )

    return actions


def build_well_measurement_options(
    actions: Sequence[WellAction],
    hypothesis_ids: Sequence[str],
    topology_edge_tuples: Sequence[Tuple[str, ...]],
    candidate_edges: Sequence[Edge],
    samples: Optional[Sequence[Dict[str, Any]]] = None,
) -> List[MeasurementOption]:
    """Convert WellActions into MeasurementOptions with predictive Gaussian scenarios."""
    n_hypotheses = len(hypothesis_ids)
    options: List[MeasurementOption] = []

    for action in actions:
        wid = action.well_id
        m_type = action.measurement_type

        if action.scenarios and len(action.scenarios) == n_hypotheses:
            means = np.array([s.mean for s in action.scenarios], dtype=float)
            sds = np.array([s.sd for s in action.scenarios], dtype=float)
        else:
            # Do not invent a predictive observation model here.  The generic
            # Bayesian engine must abstain when the action has no complete,
            # probability-bearing scenario set.
            options.append(
                MeasurementOption(
                    option_id=action.action_id,
                    measurement_type=action.measurement_type,
                    target_id=f"{action.well_id}@{action.measurement_type}",
                    cost=action.standalone_cost,
                    scenarios=(),
                    feasible=False,
                    metadata={
                        "well_id": action.well_id,
                        "sampling_time": action.sampling_time,
                        "base_cost": action.base_cost,
                        "travel_cost": action.travel_cost,
                        "accessibility": action.accessibility,
                        "connected_edges": list(action.connected_edge_ids),
                        "predictive_model": action.metadata.get("predictive_model", "missing"),
                        "predictive_model_error": action.metadata.get("predictive_model_error"),
                    },
                )
            )
            continue

        # Build robust scenarios: nominal, separation_stress, noise_stress
        grand_mean = float(np.mean(means))
        separation_means = grand_mean + 0.65 * (means - grand_mean)
        separation_sds = sds * 1.25
        noise_sds = sds * 1.60

        robust_scenarios = (
            PredictiveScenario("nominal", means.tolist(), sds.tolist(), weight=0.50),
            PredictiveScenario("separation_stress", separation_means.tolist(), separation_sds.tolist(), weight=0.25),
            PredictiveScenario("noise_stress", means.tolist(), noise_sds.tolist(), weight=0.25),
        )

        options.append(
            MeasurementOption(
                option_id=action.action_id,
                measurement_type=action.measurement_type,
                target_id=f"{action.well_id}@{action.measurement_type}",
                cost=action.standalone_cost,
                scenarios=robust_scenarios,
                feasible=action.feasible,
                metadata={
                    "well_id": action.well_id,
                    "sampling_time": action.sampling_time,
                    "base_cost": action.base_cost,
                    "travel_cost": action.travel_cost,
                    "accessibility": action.accessibility,
                    "connected_edges": list(action.connected_edge_ids),
                    "predictive_model": action.metadata.get("predictive_model", "unknown"),
                    "predictive_model_error": action.metadata.get("predictive_model_error"),
                },
            )
        )

    return options


# ── decision target matrix (targeting FP / FN / Unlabeled edges) ──────


class ValidationDecisionResult:
    """Composite result supporting both 3-tuple unpacking and dictionary access."""

    def __init__(
        self,
        matrix: np.ndarray,
        target_ids: List[str],
        weights: Dict[str, float],
        details: Dict[str, Dict[str, Any]],
    ):
        self.matrix = matrix
        self.target_ids = target_ids
        self.weights = weights
        self.details = details

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, int):
            return (self.matrix, self.target_ids, self.weights)[key]
        return self.details[key]

    def __iter__(self):
        return iter((self.matrix, self.target_ids, self.weights))

    def __len__(self) -> int:
        return 3

    def get(self, key: str, default: Any = None) -> Any:
        return self.details.get(key, default)

    def __contains__(self, key: Any) -> bool:
        return key in self.details


def build_validation_decision_matrix(
    topology_edge_tuples_or_candidates: Any = None,
    candidate_edges_or_report: Any = None,
    validation_status_map: Optional[Mapping[str, Any]] = None,
    *,
    candidate_edges: Optional[Sequence[Edge]] = None,
    topology_edge_tuples: Optional[Sequence[Tuple[str, ...]]] = None,
) -> ValidationDecisionResult:
    """Construct a decision target matrix targeting FP, FN, and unlabeled edges.

    Weights are explicitly scaled by validation status:
      - False Positive (FP): W = 2.0 (high priority to refute false link)
      - False Negative (FN): W = 1.8 (high priority to recover missed flow path)
      - Selected-Unlabeled / Ambiguous: W = 1.4 (needs independent verification)
      - True Positive / True Negative / Validated: W = 0.2 (already confirmed)
      - Other / Default: W = 1.0

    Supports both (topology_edge_tuples, candidate_edges, status_map) and
    (candidate_edges, validation_report) calling styles.
    """
    edges: Sequence[Edge] = ()
    tuples: Sequence[Tuple[str, ...]] = ()
    status_input: Any = validation_status_map

    # Disambiguate arguments
    if isinstance(topology_edge_tuples_or_candidates, Sequence) and len(topology_edge_tuples_or_candidates) > 0 and isinstance(topology_edge_tuples_or_candidates[0], Edge):
        edges = topology_edge_tuples_or_candidates
        status_input = candidate_edges_or_report
        tuples = topology_edge_tuples or ()
    elif isinstance(candidate_edges_or_report, Sequence) and len(candidate_edges_or_report) > 0 and isinstance(candidate_edges_or_report[0], Edge):
        edges = candidate_edges_or_report
        tuples = topology_edge_tuples_or_candidates or ()
    elif candidate_edges is not None:
        edges = candidate_edges
        tuples = topology_edge_tuples or ()

    target_edge_ids = [e.edge_id for e in edges]
    if not target_edge_ids:
        target_edge_ids = ["mock_target"]

    # Extract status mapping
    raw_status_map: Dict[str, str] = {}
    if isinstance(status_input, dict):
        if "labels" in status_input:
            for k, v in status_input["labels"].items():
                if isinstance(v, dict):
                    raw_status_map[k] = str(v.get("status", "unlabeled"))
                else:
                    raw_status_map[k] = str(v)
        elif "validation_status_map" in status_input:
            raw_status_map = dict(status_input["validation_status_map"])
        else:
            for k, v in status_input.items():
                if isinstance(v, dict):
                    raw_status_map[k] = str(v.get("validation_status", v.get("status", "unlabeled")))
                else:
                    raw_status_map[k] = str(v)

    weights_dict: Dict[str, float] = {}
    details_dict: Dict[str, Dict[str, Any]] = {}

    for eid in target_edge_ids:
        raw = raw_status_map.get(eid)
        if raw is not None:
            r_str = str(raw).strip()
            if r_str in {"false_positive", "FP"}:
                w = 2.0
                v_stat = "FP"
            elif r_str in {"false_negative", "FN"}:
                w = 1.8
                v_stat = "FN"
            elif r_str in {"selected_unlabeled", "ambiguous", "AMBIGUOUS"}:
                w = 1.4
                v_stat = "selected_unlabeled"
            elif r_str in {"observed_present", "VALIDATED", "TP"}:
                w = 0.2
                v_stat = "TP"
            elif r_str in {"observed_absent", "TN"}:
                w = 0.2
                v_stat = "TN"
            else:
                w = 1.0
                v_stat = r_str
        else:
            # Unlabeled edge
            w = 1.4
            v_stat = "selected_unlabeled"

        weights_dict[eid] = w
        details_dict[eid] = {
            "edge_id": eid,
            "validation_status": v_stat,
            "risk_weight": w,
        }

    n_hypotheses = max(1, len(tuples))
    n_targets = len(target_edge_ids)
    decision_matrix = np.zeros((n_hypotheses, n_targets), dtype=float)

    if tuples:
        edge_sets = [set(t) for t in tuples]
        for k, active_set in enumerate(edge_sets):
            for j, eid in enumerate(target_edge_ids):
                is_pres = float(eid in active_set)
                decision_matrix[k, j] = weights_dict[eid] * is_pres
    else:
        for j, eid in enumerate(target_edge_ids):
            decision_matrix[0, j] = weights_dict[eid]

    return ValidationDecisionResult(
        matrix=decision_matrix,
        target_ids=target_edge_ids,
        weights=weights_dict,
        details=details_dict,
    )


# ── abstention gate ──────────────────────────────────────────────────


class AbstentionReasons(list):
    """List of reasons supporting substring containment checks."""

    def __contains__(self, item: Any) -> bool:
        if super().__contains__(item):
            return True
        s = str(item).lower()
        return any(s in str(r).lower() for r in self)

    def __str__(self) -> str:
        return "; ".join(self)

    def lower(self) -> str:
        return str(self).lower()


def evaluate_abstention_gate(
    posterior_result: Dict[str, Any],
    benchmark_report: Optional[Dict[str, Any]] = None,
    candidate_options: Optional[Sequence[Any]] = None,
    config: Optional[CampaignConfig] = None,
    prior_entropy: Optional[float] = None,
    max_robust_eig: Optional[float] = None,
    prior_probabilities: Optional[Sequence[float]] = None,
) -> Tuple[bool, AbstentionReasons]:
    """Evaluate diagnostic, validation, and predictive model gates.

    Returns (should_abstain: bool, abstention_reasons: AbstentionReasons).
    """
    cfg = config or CampaignConfig()
    reasons = AbstentionReasons()

    # 1. Posterior completion & convergence gate
    status = str(posterior_result.get("status", "")).lower()
    if status and status not in {"completed", "converged", "success"}:
        reasons.append(f"Posterior sampling did not complete successfully (status='{posterior_result.get('status')}').")

    r_hat = posterior_result.get("n_edges_r_hat")
    if r_hat is None and "attrs" in posterior_result:
        r_hat = posterior_result["attrs"].get("mcmc_r_hat")
    if r_hat is None:
        r_hat = posterior_result.get("r_hat")

    if r_hat is not None and not math.isnan(float(r_hat)):
        if float(r_hat) > cfg.max_r_hat:
            reasons.append(f"MCMC chains did not converge: R-hat={float(r_hat):.3f} > {cfg.max_r_hat}.")

    edge_rhats = posterior_result.get("edge_r_hat", {})
    if edge_rhats:
        high_rhat_edges = [
            eid for eid, rh in edge_rhats.items()
            if rh is not None and math.isfinite(float(rh)) and float(rh) > cfg.max_r_hat + 0.05
        ]
        if len(high_rhat_edges) > max(3, int(0.25 * len(edge_rhats))):
            reasons.append(
                f"{len(high_rhat_edges)} candidate edges exceed R-hat convergence threshold {cfg.max_r_hat}."
            )

    ess = posterior_result.get("n_edges_ess")
    if ess is None and "attrs" in posterior_result:
        ess = posterior_result["attrs"].get("mcmc_ess")
    if ess is None:
        ess = posterior_result.get("ess")

    if ess is not None and float(ess) < cfg.min_ess:
        reasons.append(f"MCMC effective sample size is inadequate (ESS={float(ess):.1f} < {cfg.min_ess}).")

    if cfg.require_mcmc_diagnostics and (r_hat is None or ess is None):
        reasons.append(
            "MCMC convergence diagnostics are unavailable; run at least two "
            "posterior chains or explicitly disable this gate."
        )

    if prior_entropy is not None:
        entropy = float(prior_entropy)
    else:
        ent = posterior_result.get("joint_graph_entropy")
        if ent is None and "attrs" in posterior_result:
            ent = posterior_result["attrs"].get("posterior_joint_graph_entropy")
        if ent is None:
            ent = posterior_result.get("entropy", 0.0)
        entropy = float(ent) if ent is not None else 0.0

    if entropy < cfg.min_entropy:
        reasons.append(f"Posterior topology entropy is near zero ({entropy:.6f} nats); no uncertainty remains to resolve.")

    accept_rate = float(posterior_result.get("acceptance_rate", 0.25))
    if accept_rate < 0.005 or accept_rate > 0.99:
        reasons.append(f"MCMC acceptance rate is degenerate ({accept_rate:.4f}).")

    if max_robust_eig is not None and float(max_robust_eig) < cfg.minimum_robust_eig:
        reasons.append(f"Maximum robust information gain ({float(max_robust_eig):.6f}) is below campaign threshold ({cfg.minimum_robust_eig}).")

    # 2. Benchmark & independent validation gate
    if benchmark_report is not None:
        if not benchmark_report.get("variants"):
            reasons.append("Benchmark report is missing or contains no variants.")
        if cfg.require_independent_validation and not benchmark_report.get("independent_validation", False):
            reasons.append("Independent validation labels required but not satisfied in benchmark report.")

    # 3. Measurement options & predictive model gate
    if candidate_options is not None:
        missing_predictions = [
            str(getattr(opt, "option_id", "unknown"))
            for opt in candidate_options
            if getattr(opt, "feasible", True) and not getattr(opt, "scenarios", ())
        ]
        if missing_predictions:
            reasons.append(
                "Predictive measurement scenarios are missing for "
                f"{len(missing_predictions)} feasible actions."
            )
        feasible_options = [opt for opt in candidate_options if getattr(opt, "feasible", True)]
        if not feasible_options:
            reasons.append("No feasible measurement actions are available in the campaign domain.")
        if max_robust_eig is None and prior_probabilities is not None and feasible_options:
            try:
                prior = np.asarray(prior_probabilities, dtype=float)
                eig_values = [
                    _robust_joint_eig(prior, [opt], cfg.acquisition_config)[2]
                    for opt in feasible_options
                    if getattr(opt, "scenarios", ())
                ]
                if eig_values and max(eig_values) < cfg.minimum_robust_eig:
                    reasons.append(
                        "Maximum robust information gain "
                        f"({max(eig_values):.6f}) is below campaign threshold "
                        f"({cfg.minimum_robust_eig:.6f})."
                    )
            except Exception as exc:
                reasons.append(f"Predictive-model adequacy could not be evaluated: {exc}")

    should_abstain = (len(reasons) > 0)
    return should_abstain, reasons


# ── well-level batch selection with shared travel costs ───────────────


def select_well_campaign_batch(
    candidate_options: Sequence[Union[MeasurementOption, WellAction]],
    prior_probabilities: Sequence[float],
    config: Optional[CampaignConfig] = None,
    batch_size: Optional[int] = None,
    budget: Optional[float] = None,
    max_actions_per_well: Optional[int] = None,
    decision_values: Optional[Sequence[Sequence[float]]] = None,
) -> Dict[str, Any]:
    """Greedy forward batch selection of well measurement actions."""
    cfg = config or CampaignConfig()
    batch_size = batch_size if batch_size is not None else cfg.batch_size
    budget = budget if budget is not None else cfg.budget
    max_actions_per_well = (
        max_actions_per_well if max_actions_per_well is not None else cfg.max_actions_per_well
    )
    prior = list(prior_probabilities)
    ac_cfg = cfg.acquisition_config
    decision_array = None
    prior_decision_risk = 0.0
    if decision_values is not None:
        decision_array = np.asarray(decision_values, dtype=float)
        if decision_array.ndim == 1:
            decision_array = decision_array[:, None]
        if decision_array.ndim != 2 or decision_array.shape[0] != len(prior):
            raise ValueError("decision_values must have one row per hypothesis.")
        prior_array = np.asarray(prior, dtype=float)
        prior_array /= prior_array.sum()
        prior_decision_risk = _prior_decision_risk(prior_array, decision_array)

    # Convert WellAction to MeasurementOption if needed
    converted_options: List[MeasurementOption] = []
    for item in candidate_options:
        if isinstance(item, WellAction):
            if item.scenarios and len(item.scenarios) == len(prior):
                m_vec = [s.mean for s in item.scenarios]
                s_vec = [s.sd for s in item.scenarios]
                scenarios = (
                    PredictiveScenario("nominal", m_vec, s_vec, weight=0.50),
                    PredictiveScenario("separation_stress", m_vec, [s * 1.25 for s in s_vec], weight=0.25),
                    PredictiveScenario("noise_stress", m_vec, [s * 1.60 for s in s_vec], weight=0.25),
                )
            else:
                scenarios = ()

            converted_options.append(
                MeasurementOption(
                    option_id=item.action_id,
                    measurement_type=item.measurement_type,
                    target_id=f"{item.well_id}@{item.measurement_type}",
                    cost=item.standalone_cost,
                    scenarios=scenarios,
                    feasible=bool(item.feasible and scenarios),
                    metadata={
                        "well_id": item.well_id,
                        "sampling_time": item.sampling_time,
                        "base_cost": item.base_cost,
                        "travel_cost": item.travel_cost,
                        "accessibility": item.accessibility,
                        "connected_edges": list(item.connected_edge_ids),
                    },
                )
            )
        else:
            converted_options.append(item)

    available = [opt for opt in converted_options if opt.feasible and opt.scenarios]
    if not available:
        return {
            "status": "ABSTAIN",
            "selected_options": [],
            "selected_option_details": [],
            "visited_wells": [],
            "total_cost": 0.0,
            "travel_cost_savings": 0.0,
            "joint_robust_information_gain": 0.0,
            "prior_entropy": shannon_entropy(prior),
            "reason": "No feasible measurement actions available.",
        }

    selected: List[MeasurementOption] = []
    selection_rows: List[Dict[str, Any]] = []
    visited_wells: Set[str] = set()
    well_action_counts: Dict[str, int] = {}
    total_cost = 0.0
    nominal_unshared_cost = 0.0
    current_joint = 0.0

    for _ in range(min(batch_size, len(available))):
        candidates = []
        for opt in available:
            if opt in selected:
                continue

            wid = str(opt.metadata.get("well_id", opt.target_id.split("@")[0]))
            count_at_well = well_action_counts.get(wid, 0)
            if count_at_well >= max_actions_per_well:
                continue

            accessibility = max(0.05, float(opt.metadata.get("accessibility", 1.0)))
            base_c = float(opt.metadata.get("base_cost", opt.cost))
            travel_c = float(opt.metadata.get("travel_cost", cfg.travel_cost_per_well))

            if wid in visited_wells:
                incremental_cost = base_c / accessibility
            else:
                incremental_cost = (base_c + travel_c) / accessibility

            proposed_cost = total_cost + incremental_cost
            if budget is not None and proposed_cost > float(budget) + 1.0e-12:
                continue

            mean_eig, worst_eig, robust_eig, scenarios = _robust_joint_eig(
                prior, [*selected, opt], ac_cfg
            )
            marginal_eig = max(0.0, robust_eig - current_joint)
            marginal_decision_risk = 0.0
            if decision_array is not None:
                _, _, _, _, _, marginal_decision_risk, _ = _scenario_eigs(
                    np.asarray(prior, dtype=float),
                    opt,
                    ac_cfg,
                    decision_array,
                    prior_decision_risk,
                )
            information_fraction = marginal_eig / max(
                shannon_entropy(prior), 1.0e-15
            )
            decision_fraction = marginal_decision_risk / max(
                prior_decision_risk, 1.0e-15
            )
            selection_utility = (
                (1.0 - ac_cfg.decision_weight) * information_fraction
                + ac_cfg.decision_weight * decision_fraction
            )
            score = selection_utility / (incremental_cost ** ac_cfg.cost_exponent)

            candidates.append((
                -score,
                opt.option_id,
                opt,
                incremental_cost,
                wid,
                mean_eig,
                worst_eig,
                robust_eig,
                marginal_eig,
                marginal_decision_risk,
                selection_utility,
                scenarios,
            ))

        if not candidates:
            break

        (
            _,
            _,
            best_opt,
            inc_cost,
            wid,
            mean_eig,
            worst_eig,
            robust_eig,
            marginal_eig,
            marginal_decision_risk,
            selection_utility,
            scenarios,
        ) = min(candidates, key=lambda x: (x[0], x[1]))

        if marginal_eig < cfg.minimum_robust_eig and len(selected) > 0:
            break

        selected.append(best_opt)
        visited_wells.add(wid)
        well_action_counts[wid] = well_action_counts.get(wid, 0) + 1
        total_cost += inc_cost
        nominal_unshared_cost += float(best_opt.cost)
        current_joint = robust_eig

        selection_rows.append({
            "batch_rank": len(selected),
            "option_id": best_opt.option_id,
            "well_id": wid,
            "measurement_type": best_opt.measurement_type,
            "marginal_cost": float(inc_cost),
            "standalone_cost": float(best_opt.cost),
            "cumulative_cost": float(total_cost),
            "marginal_robust_information_gain": float(marginal_eig),
            "marginal_decision_risk_reduction": float(marginal_decision_risk),
            "selection_utility": float(selection_utility),
            "joint_robust_information_gain": float(robust_eig),
            "scenario_scores": scenarios,
        })

    travel_savings = max(0.0, nominal_unshared_cost - total_cost)
    actionable = (len(selected) > 0)

    return {
        "status": "ACTIONABLE" if actionable else "ABSTAIN",
        "selected_options": [opt.option_id for opt in selected],
        "selected_option_details": selection_rows,
        "visited_wells": sorted(visited_wells),
        "total_cost": float(total_cost),
        "travel_cost_savings": float(travel_savings),
        "joint_robust_information_gain": float(current_joint),
        "prior_entropy": shannon_entropy(prior),
        "claim_guardrail": (
            "Well-level recommendations are experimental-design decisions conditional on the "
            "posterior topology ensemble, observation error model, and specified access constraints."
        ),
    }


# ── main entry point ────────────────────────────────────────────────


def rank_campaign_measurements(
    posterior_result: Dict[str, Any],
    benchmark_report: Dict[str, Any],
    candidate_edges: Sequence[Edge],
    samples: Optional[Sequence[Dict[str, Any]]] = None,
    validation_report: Optional[Dict[str, Any]] = None,
    *,
    config: Optional[CampaignConfig] = None,
    accessibility_map: Optional[Mapping[str, float]] = None,
    feasibility_map: Optional[Mapping[str, bool]] = None,
    predictive_model: Optional[PredictiveModel] = None,
    output_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """Rank well-level measurement actions and select non-redundant campaign batches."""
    cfg = config or CampaignConfig()

    # 1. Extract discrete topology hypotheses from posterior result
    hypothesis_ids, prior_probs, edge_tuples = extract_topology_hypotheses(
        posterior_result,
        min_probability=0.0,
    )
    if len(hypothesis_ids) < 2:
        result = {
            "status": "ABSTAIN",
            "abstention_reasons": [
                "A probability-bearing topology ensemble with at least two "
                "distinct hypotheses is required for active learning."
            ],
            "prior_topology_entropy": 0.0,
            "n_hypotheses": len(hypothesis_ids),
            "rankings": [],
            "campaign_batch": {
                "status": "ABSTAIN",
                "selected_options": [],
                "visited_wells": [],
                "total_cost": 0.0,
            },
            "summary": {
                "status": "ABSTAIN",
                "n_recommendations": 0,
                "top_priority_score": 0.0,
                "n_actions_scored": 0,
                "n_hypotheses": len(hypothesis_ids),
                "reasons": [
                    "A probability-bearing topology ensemble with at least two "
                    "distinct hypotheses is required for active learning."
                ],
            },
            "claim_guardrail": (
                "Abstention triggered: a single MAP topology cannot support "
                "information-gain-based measurement selection."
            ),
        }
        if output_dir:
            _write_campaign_outputs(result, output_dir)
        return result

    # 2. Build well measurement actions
    if isinstance(samples, Mapping):
        sample_list = list(samples.values())
    else:
        sample_list = list(samples) if samples is not None else []
    actions = build_well_actions(
        sample_nodes=sample_list,
        candidate_edges=candidate_edges,
        hypothesis_ids=hypothesis_ids,
        prior_probs=prior_probs,
        edge_tuples=edge_tuples,
        config=cfg,
        accessibility_map=accessibility_map,
        feasibility_map=feasibility_map,
        predictive_model=predictive_model,
    )

    # 3. Convert actions to MeasurementOptions with predictive scenarios
    options = build_well_measurement_options(
        actions=actions,
        hypothesis_ids=hypothesis_ids,
        topology_edge_tuples=edge_tuples,
        candidate_edges=candidate_edges,
        samples=sample_list,
    )

    # 4. Build validation status map and decision matrix
    val_status_map: Dict[str, str] = {}
    cal_selected = set(
        benchmark_report.get("variants", {})
        .get("assumption_calibrated", {})
        .get("selected_edge_ids", [])
    )
    if validation_report:
        obs_by_edge: Dict[str, float] = {}
        val_label_file = validation_report.get("validation_label_file")
        candidate_label_paths = []
        if val_label_file:
            candidate_label_paths.append(Path(str(val_label_file)))
            if output_dir and not Path(str(val_label_file)).is_absolute():
                candidate_label_paths.append(Path(output_dir) / str(val_label_file))
        existing_label_path = next(
            (path for path in candidate_label_paths if path.is_file()), None
        )
        if existing_label_path is not None:
            try:
                from .validation_workflow import _load_topology_observations
                val_obs = _load_topology_observations(str(existing_label_path))
                for o in val_obs:
                    obs_by_edge[o.edge_id] = o.observed_present
            except Exception:
                pass

        # Also accept an in-memory validation report. This keeps the adapter
        # usable outside the CLI and prevents inline labels from being silently
        # treated as unlabeled edges.
        labels = validation_report.get("labels", validation_report.get("validation_labels", {}))
        if isinstance(labels, Mapping):
            for edge_id, label in labels.items():
                if isinstance(label, Mapping):
                    if "status" in label:
                        val_status_map[str(edge_id)] = str(label["status"])
                    elif "expected" in label and "inferred" in label:
                        expected = bool(label["expected"])
                        inferred = bool(label["inferred"])
                        if expected and inferred:
                            val_status_map[str(edge_id)] = "observed_present"
                        elif expected and not inferred:
                            val_status_map[str(edge_id)] = "false_negative"
                        elif not expected and inferred:
                            val_status_map[str(edge_id)] = "false_positive"
                        else:
                            val_status_map[str(edge_id)] = "observed_absent"
                else:
                    val_status_map[str(edge_id)] = str(label)

    for e in candidate_edges:
        eid = e.edge_id
        if eid in val_status_map:
            continue
        if validation_report and eid in obs_by_edge:
            obs = obs_by_edge[eid]
            if obs >= 0.5:
                val_status_map[eid] = "observed_present" if eid in cal_selected else "false_negative"
            else:
                val_status_map[eid] = "false_positive" if eid in cal_selected else "observed_absent"
        elif eid in cal_selected:
            val_status_map[eid] = "selected_unlabeled"
        else:
            val_status_map[eid] = "unlabeled"

    decision_res = build_validation_decision_matrix(
        topology_edge_tuples=edge_tuples,
        candidate_edges=candidate_edges,
        validation_status_map=val_status_map,
    )
    decision_matrix, target_ids, weights = decision_res.matrix, decision_res.target_ids, decision_res.weights

    # 5. Abstention Gate
    should_abstain, abstention_reasons = evaluate_abstention_gate(
        posterior_result=posterior_result,
        benchmark_report=benchmark_report,
        candidate_options=options,
        config=cfg,
        prior_probabilities=prior_probs,
    )

    prior_entropy = shannon_entropy(prior_probs)

    if should_abstain:
        result = {
            "status": "ABSTAIN",
            "abstention_reasons": list(abstention_reasons),
            "prior_topology_entropy": prior_entropy,
            "n_hypotheses": len(hypothesis_ids),
            "rankings": [],
            "campaign_batch": {
                "status": "ABSTAIN",
                "selected_options": [],
                "visited_wells": [],
                "total_cost": 0.0,
            },
            "summary": {
                "status": "ABSTAIN",
                "n_recommendations": 0,
                "top_priority_score": 0.0,
                "n_actions_scored": len(options),
                "n_hypotheses": len(hypothesis_ids),
                "reasons": list(abstention_reasons),
            },
            "claim_guardrail": (
                "Abstention triggered: active learning recommendations are decision-support "
                "guidance only and cannot proceed under degenerate or unvalidated inference conditions."
            ),
        }
        if output_dir:
            _write_campaign_outputs(result, output_dir)
        return result

    # 6. Rank individual measurement options
    ranking_res = rank_measurement_options(
        hypothesis_ids=hypothesis_ids,
        prior_probabilities=prior_probs,
        options=options,
        decision_values=decision_matrix,
        config=cfg.acquisition_config,
    )
    rankings = ranking_res.get("rankings", [])
    for r in rankings:
        r["acquisition_score"] = r.get("cost_adjusted_score", 0.0)
        r["expected_decision_risk_reduction"] = r.get("expected_brier_risk_reduction", 0.0)

    if ranking_res.get("status") != "ACTIONABLE":
        result = {
            "status": "ABSTAIN",
            "abstention_reasons": [
                "No measurement action exceeded the robust acquisition threshold."
            ],
            "prior_topology_entropy": prior_entropy,
            "n_hypotheses": len(hypothesis_ids),
            "rankings": rankings,
            "campaign_batch": {
                "status": "ABSTAIN",
                "selected_options": [],
                "visited_wells": [],
                "total_cost": 0.0,
            },
            "summary": {
                "status": "ABSTAIN",
                "n_actions_scored": len(options),
                "n_recommendations": 0,
                "top_priority_score": 0.0,
                "n_hypotheses": len(hypothesis_ids),
                "reasons": [
                    "No measurement action exceeded the robust acquisition threshold."
                ],
            },
            "claim_guardrail": ranking_res.get("claim_guardrail"),
        }
        if output_dir:
            _write_campaign_outputs(result, output_dir)
        return result

    # Attach well-level metadata to rankings
    opt_by_id = {opt.option_id: opt for opt in options}
    augmented_rankings = []
    for r in rankings:
        opt = opt_by_id.get(r["option_id"])
        augmented_r = dict(r)
        if opt:
            augmented_r["well_id"] = opt.metadata.get("well_id")
            augmented_r["sampling_time"] = opt.metadata.get("sampling_time")
            augmented_r["base_cost"] = opt.metadata.get("base_cost")
            augmented_r["travel_cost"] = opt.metadata.get("travel_cost")
            augmented_r["accessibility"] = opt.metadata.get("accessibility")
            augmented_r["connected_edges"] = opt.metadata.get("connected_edges", [])
        augmented_rankings.append(augmented_r)

    # 7. Select non-redundant well campaign batch
    campaign_batch = select_well_campaign_batch(
        candidate_options=options,
        prior_probabilities=prior_probs,
        config=cfg,
        decision_values=decision_matrix,
    )

    summary = {
        "status": campaign_batch.get("status", "ABSTAIN"),
        "n_recommendations": len(augmented_rankings),
        "n_hypotheses": len(hypothesis_ids),
        "n_actions_scored": len(options),
        "prior_topology_entropy": float(prior_entropy),
        "top_action_id": augmented_rankings[0]["option_id"] if augmented_rankings else None,
        "top_acquisition_score": augmented_rankings[0]["acquisition_score"] if augmented_rankings else 0.0,
        "top_priority_score": augmented_rankings[0]["acquisition_score"] if augmented_rankings else 0.0,
        "batch_size_selected": len(campaign_batch["selected_options"]),
        "visited_wells_count": len(campaign_batch["visited_wells"]),
        "total_campaign_cost": campaign_batch["total_cost"],
        "travel_cost_savings": campaign_batch["travel_cost_savings"],
        "predictive_model": cfg.predictive_model_name if predictive_model is None else "configured_callback",
        "decision_weight": cfg.acquisition_config.decision_weight,
    }

    result = {
        "status": summary["status"],
        "summary": summary,
        "rankings": augmented_rankings,
        "campaign_batch": campaign_batch,
        "validation_targeting": {
            "target_edges": target_ids,
            "weights": weights,
            "statuses": {eid: decision_res[eid]["validation_status"] for eid in target_ids},
        },
        "posterior_diagnostics": {
            "n_hypotheses": len(hypothesis_ids),
            "prior_entropy": float(prior_entropy),
            "r_hat": posterior_result.get("n_edges_r_hat", posterior_result.get("attrs", {}).get("mcmc_r_hat")),
            "ess": posterior_result.get("n_edges_ess", posterior_result.get("attrs", {}).get("mcmc_ess")),
            "joint_graph_entropy": posterior_result.get(
                "joint_graph_entropy",
                posterior_result.get("attrs", {}).get("posterior_joint_graph_entropy"),
            ),
            "n_unique_graphs": posterior_result.get("n_unique_graphs"),
        },
        "claim_guardrail": (
            "Well-level recommendations are decision-support tools for targeted field sampling, "
            "conditional on the declared predictive model. "
            "They do not validate or falsify any groundwater topology edge until field data are acquired and tested."
        ),
    }

    if output_dir:
        _write_campaign_outputs(result, output_dir)

    return result


def _write_campaign_outputs(report: Dict[str, Any], output_dir: str) -> None:
    """Write JSON, CSV, and Markdown campaign reports."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. JSON report
    json_path = out / "well_campaign_recommendations.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    # 2. CSV report (rankings)
    csv_path = out / "well_campaign_recommendations.csv"
    rankings = report.get("rankings", [])
    if rankings:
        fieldnames = [
            "rank",
            "action_id",
            "well_id",
            "measurement_type",
            "sampling_time",
            "cost",
            "expected_information_gain",
            "expected_decision_risk_reduction",
            "acquisition_score",
            "worst_case_information_gain",
        ]
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for r in rankings:
                writer.writerow(r)

    # 3. Markdown summary
    md_path = out / "well_campaign_recommendations.md"
    _write_campaign_markdown(report, md_path)


def _write_campaign_markdown(report: Dict[str, Any], path: Path) -> None:
    """Generate manuscript-ready Markdown report for campaign design."""
    lines = [
        "# Well-Level Active Learning Campaign Design",
        "",
        "> **Decision-Support Notice**: Measurement actions recommended below are experimental-design",
        "> decisions targeting maximum topology uncertainty reduction. They do not validate any",
        "> candidate edge without independent field observation.",
        "",
    ]

    status = report.get("status", "UNKNOWN")
    lines += [
        f"**Status**: `{status}`",
        "",
    ]

    if status == "ABSTAIN":
        lines += [
            "## Campaign Abstention",
            "",
            "The active learning engine abstained from recommending measurements for the following reasons:",
            "",
        ]
        for reason in report.get("abstention_reasons", []):
            lines.append(f"- {reason}")
        lines.append("")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    summary = report.get("summary", {})
    batch = report.get("campaign_batch", {})
    diag = report.get("posterior_diagnostics", {})

    lines += [
        "## Executive Summary",
        "",
        f"- **Posterior Hypotheses Visited**: {summary.get('n_hypotheses')}",
        f"- **Prior Topology Entropy**: {summary.get('prior_topology_entropy', 0.0):.4f} nats",
        f"- **MCMC Diagnostics**: R-hat={diag.get('r_hat', 'N/A')}, ESS={diag.get('ess', 'N/A')}",
        f"- **Predictive Model**: {summary.get('predictive_model', 'not reported')}",
        f"- **Decision-Risk Weight**: {summary.get('decision_weight', 'not reported')}",
        f"- **Total Campaign Cost**: {summary.get('total_campaign_cost', 0.0):.2f}",
        f"- **Mobilization / Travel Savings**: {summary.get('travel_cost_savings', 0.0):.2f}",
        f"- **Distinct Wells Visited**: {summary.get('visited_wells_count')} ({', '.join(batch.get('visited_wells', []))})",
        "",
        "## Selected Campaign Batch",
        "",
        "| Batch Rank | Action ID | Well ID | Measurement Type | Marginal Cost | Standalone Cost | Cumulative Cost | Marginal Robust EIG |",
        "|---|---|---|---|---|---|---|---|",
    ]

    for item in batch.get("selected_option_details", []):
        lines.append(
            f"| {item['batch_rank']} | `{item['option_id']}` | **{item['well_id']}** "
            f"| {item['measurement_type']} | {item['marginal_cost']:.2f} | {item['standalone_cost']:.2f} "
            f"| {item['cumulative_cost']:.2f} | {item['marginal_robust_information_gain']:.4f} |"
        )

    lines += [
        "",
        "## Top Individual Actions Ranked",
        "",
        "| Rank | Action ID | Well | Type | Standalone Cost | Robust EIG | Decision Risk Red. | Score |",
        "|---|---|---|---|---|---|---|---|",
    ]

    for r in report.get("rankings", [])[:10]:
        lines.append(
            f"| {r.get('rank')} | `{r.get('option_id')}` | **{r.get('well_id')}** "
            f"| {r.get('measurement_type')} | {r.get('cost', 0.0):.2f} "
            f"| {r.get('expected_information_gain', 0.0):.4f} "
            f"| {r.get('expected_decision_risk_reduction', 0.0):.4f} "
            f"| {r.get('acquisition_score', 0.0):.4f} |"
        )

    lines += [
        "",
        "---",
        f"*{report.get('claim_guardrail')}*",
        "",
    ]

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
