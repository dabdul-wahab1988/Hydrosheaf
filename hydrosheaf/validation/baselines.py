"""Auditable specialist baselines for controlled-synthetic validation.

This module is intentionally disjoint from the HydroSheaf inference pipeline.
It defines truth-blind baseline scorers that consume only declared observation
channels and a caller-supplied candidate universe.  The goal is fair comparison:
every registered baseline exposes the same audit metadata for inputs, tuning,
uncertainty, abstention, and decision cost.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from collections.abc import Callable, Iterable, Mapping
from typing import Any


Edge = tuple[str, str]
Decision = str
ObservationMap = Mapping[str, Any]
Scorer = Callable[[Iterable[Any], ObservationMap], tuple["BaselinePrediction", ...]]

SELECT = "select"
REJECT = "reject"
ABSTAIN = "abstain"

FORBIDDEN_TRUTH_KEYS = frozenset(
    {
        "truth",
        "true",
        "label",
        "labels",
        "reference",
        "reference_edges",
        "target",
        "y",
        "is_true",
        "ground_truth",
    }
)


@dataclass(frozen=True)
class BaselinePrediction:
    """Per-edge baseline output used by the validation benchmark."""

    source: str
    target: str
    probability: float
    decision: Decision
    evidence_channel: str
    reason: str
    score: float | None = None
    cost: Mapping[str, float] = field(default_factory=dict)

    @property
    def edge(self) -> Edge:
        return (self.source, self.target)

    def to_audit_record(self) -> dict[str, Any]:
        """Return a JSON-serialisable record for locked benchmark outputs."""

        return {
            "edge": [self.source, self.target],
            "probability": self.probability,
            "decision": self.decision,
            "evidence_channel": self.evidence_channel,
            "reason": self.reason,
            "score": self.score,
            "cost": dict(self.cost),
        }


@dataclass(frozen=True)
class BaselineSpec:
    """Complete, auditable description of a validation baseline."""

    name: str
    version: str
    family: str
    input_channels: tuple[str, ...]
    tuning: Mapping[str, Any]
    uncertainty: Mapping[str, Any]
    abstention: Mapping[str, Any]
    cost: Mapping[str, Any]
    scorer: Scorer
    description: str = ""
    control: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("BaselineSpec.name must be non-empty.")
        if not self.input_channels:
            raise ValueError("BaselineSpec.input_channels must be non-empty.")
        for section_name in ("tuning", "uncertainty", "abstention", "cost"):
            section = getattr(self, section_name)
            if not isinstance(section, Mapping):
                raise TypeError(f"BaselineSpec.{section_name} must be a mapping.")

    def score(
        self,
        candidate_universe: Iterable[Any],
        observations: ObservationMap,
    ) -> tuple[BaselinePrediction, ...]:
        """Score exactly the supplied candidate universe without truth labels."""

        assert_truth_blind_observations(observations)
        expected_edges = tuple(_normalise_edge(edge) for edge in candidate_universe)
        if len(expected_edges) != len(set(expected_edges)):
            raise ValueError("candidate_universe contains duplicate edges.")
        predictions = tuple(self.scorer(expected_edges, observations))
        seen: set[Edge] = set()
        for prediction in predictions:
            if prediction.edge in seen:
                raise ValueError(f"Duplicate prediction for edge {prediction.edge!r}.")
            seen.add(prediction.edge)
            if prediction.evidence_channel not in self.input_channels:
                raise ValueError(
                    f"{self.name} emitted channel {prediction.evidence_channel!r}, "
                    f"which is not declared in {self.input_channels!r}."
                )
            if prediction.decision not in {SELECT, REJECT, ABSTAIN}:
                raise ValueError(
                    f"{self.name} emitted unsupported decision {prediction.decision!r}."
                )
            _finite_probability(prediction.probability)
        expected_edge_set = set(expected_edges)
        if seen != expected_edge_set:
            missing = sorted(expected_edge_set - seen)
            extra = sorted(seen - expected_edge_set)
            raise ValueError(
                f"{self.name} must emit one prediction per candidate edge; "
                f"missing={missing!r}, extra={extra!r}."
            )
        return predictions

    def to_audit_record(self) -> dict[str, Any]:
        """Return stable metadata excluding the callable scorer."""

        record: dict[str, Any] = {
            "name": self.name,
            "version": self.version,
            "family": self.family,
            "description": self.description,
            "input_channels": list(self.input_channels),
            "tuning": _jsonable(self.tuning),
            "uncertainty": _jsonable(self.uncertainty),
            "abstention": _jsonable(self.abstention),
            "cost": _jsonable(self.cost),
            "control": _jsonable(self.control),
        }
        record["fingerprint"] = _fingerprint(record)
        return record


class BaselineRegistry:
    """Name-addressed registry with a shared audit contract."""

    required_metadata_sections = (
        "input_channels",
        "tuning",
        "uncertainty",
        "abstention",
        "cost",
    )

    def __init__(self, specs: Iterable[BaselineSpec] = ()) -> None:
        self._specs: dict[str, BaselineSpec] = {}
        for spec in specs:
            self.register(spec)

    def register(self, spec: BaselineSpec) -> BaselineSpec:
        if spec.name in self._specs:
            raise ValueError(f"Baseline already registered: {spec.name}")
        # Force metadata validation at registration time.
        audit = spec.to_audit_record()
        missing = [
            section
            for section in self.required_metadata_sections
            if section not in audit
        ]
        if missing:
            raise ValueError(f"{spec.name} is missing metadata sections: {missing}")
        self._specs[spec.name] = spec
        return spec

    def get(self, name: str) -> BaselineSpec:
        try:
            return self._specs[name]
        except KeyError as exc:
            raise KeyError(f"Unknown baseline: {name}") from exc

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._specs))

    def specs(self) -> tuple[BaselineSpec, ...]:
        return tuple(self._specs[name] for name in self.names())

    def audit_table(self) -> tuple[dict[str, Any], ...]:
        return tuple(self._specs[name].to_audit_record() for name in self.names())


def hydraulic_only_baseline_spec() -> BaselineSpec:
    """Deterministic hydraulic-only topology baseline."""

    abstention = {
        "rule": "abstain when probability is inside the uncertainty band",
        "select_threshold": 0.65,
        "reject_threshold": 0.35,
        "missing_evidence_decision": ABSTAIN,
    }
    return BaselineSpec(
        name="hydraulic_only_directional_gradient",
        version="1.0",
        family="hydraulic_only",
        input_channels=("hydraulic",),
        tuning={
            "gradient_scale": 0.05,
            "conductance_weight": 0.10,
            "threshold_policy": "fixed_predeclared",
        },
        uncertainty={
            "type": "heuristic_probability",
            "calibrated": False,
            "source": "directional hydraulic gradient only",
        },
        abstention=abstention,
        cost={
            "false_positive": 1.0,
            "false_negative": 1.0,
            "abstain": 0.10,
            "measurement": 0.0,
        },
        scorer=lambda candidates, observations: _score_hydraulic_only(
            candidates,
            observations,
            channel="hydraulic",
            abstention=abstention,
            cost={
                "false_positive": 1.0,
                "false_negative": 1.0,
                "abstain": 0.10,
                "measurement": 0.0,
            },
        ),
        description=(
            "Truth-blind specialist baseline using only hydraulic head/drop "
            "evidence for each candidate edge."
        ),
    )


def edge_local_topology_baseline_spec() -> BaselineSpec:
    """Deterministic topology-only edge-local baseline."""

    abstention = {
        "rule": "abstain when local support is weak or unavailable",
        "select_threshold": 0.70,
        "reject_threshold": 0.30,
        "missing_evidence_decision": ABSTAIN,
    }
    return BaselineSpec(
        name="edge_local_topology_support",
        version="1.0",
        family="topology_only",
        input_channels=("topology",),
        tuning={
            "support_scale": 1.0,
            "threshold_policy": "fixed_predeclared",
            "global_graph_features": False,
        },
        uncertainty={
            "type": "local_support_probability",
            "calibrated": False,
            "source": "edge-local candidate support only",
        },
        abstention=abstention,
        cost={
            "false_positive": 1.0,
            "false_negative": 1.0,
            "abstain": 0.10,
            "measurement": 0.0,
        },
        scorer=lambda candidates, observations: _score_edge_local_topology(
            candidates,
            observations,
            channel="topology",
            abstention=abstention,
            cost={
                "false_positive": 1.0,
                "false_negative": 1.0,
                "abstain": 0.10,
                "measurement": 0.0,
            },
        ),
        description=(
            "Truth-blind specialist baseline using only edge-local topology "
            "support from the declared topology channel."
        ),
    )


def default_baseline_registry() -> BaselineRegistry:
    """Return the default fair-comparison baseline registry."""

    return BaselineRegistry(
        (
            hydraulic_only_baseline_spec(),
            edge_local_topology_baseline_spec(),
        )
    )


def permutation_control_baseline(
    spec: BaselineSpec,
    *,
    evidence_channel: str,
    name: str | None = None,
) -> BaselineSpec:
    """Wrap a baseline so only its declared evidence channel is permuted.

    The wrapper presents ``evidence_channel`` to the scorer as the original
    baseline channel.  Candidate edges, tuning, uncertainty, abstention, and
    cost are left unchanged.  This makes channel-permutation controls auditable
    without silently changing the decision rule.
    """

    if len(spec.input_channels) != 1:
        raise ValueError("Permutation control currently supports one-channel baselines.")
    if not evidence_channel:
        raise ValueError("evidence_channel must be non-empty.")
    original_channel = spec.input_channels[0]
    control_name = name or f"{spec.name}__permute_{evidence_channel}_as_{original_channel}"

    def scorer(
        candidate_universe: Iterable[Any],
        observations: ObservationMap,
    ) -> tuple[BaselinePrediction, ...]:
        assert_truth_blind_observations(observations)
        if evidence_channel not in observations:
            remapped = dict(observations)
            remapped.pop(original_channel, None)
        else:
            remapped = dict(observations)
            remapped[original_channel] = observations[evidence_channel]
        predictions = spec.score(candidate_universe, remapped)
        return tuple(
            BaselinePrediction(
                source=prediction.source,
                target=prediction.target,
                probability=prediction.probability,
                decision=prediction.decision,
                evidence_channel=evidence_channel,
                reason=f"permutation_control:{prediction.reason}",
                score=prediction.score,
                cost=prediction.cost,
            )
            for prediction in predictions
        )

    return BaselineSpec(
        name=control_name,
        version=spec.version,
        family=spec.family,
        input_channels=(evidence_channel,),
        tuning=spec.tuning,
        uncertainty=spec.uncertainty,
        abstention=spec.abstention,
        cost=spec.cost,
        scorer=scorer,
        description=(
            f"Permutation control for {spec.name}: use channel "
            f"{evidence_channel!r} as {original_channel!r} and change no other "
            "scoring metadata."
        ),
        control={
            "type": "evidence_channel_permutation",
            "baseline": spec.name,
            "original_channel": original_channel,
            "permuted_channel": evidence_channel,
        },
    )


def assert_truth_blind_observations(observations: ObservationMap) -> None:
    """Reject observations containing obvious truth/reference labels."""

    def visit(value: Any, path: tuple[str, ...]) -> None:
        if isinstance(value, Mapping):
            for key, item in value.items():
                key_text = str(key)
                lowered = key_text.lower()
                if (
                    lowered in FORBIDDEN_TRUTH_KEYS
                    or lowered.startswith("true_")
                    or lowered.startswith("truth_")
                ):
                    dotted = ".".join((*path, key_text))
                    raise ValueError(f"Truth/reference field is forbidden: {dotted}")
                visit(item, (*path, key_text))
        elif isinstance(value, (list, tuple)):
            for index, item in enumerate(value):
                visit(item, (*path, str(index)))

    visit(observations, ())


def _score_hydraulic_only(
    candidate_universe: Iterable[Any],
    observations: ObservationMap,
    *,
    channel: str,
    abstention: Mapping[str, Any],
    cost: Mapping[str, float],
) -> tuple[BaselinePrediction, ...]:
    channel_data = observations.get(channel, {})
    predictions: list[BaselinePrediction] = []
    for candidate in candidate_universe:
        edge = _normalise_edge(candidate)
        features = _edge_features(channel_data, edge)
        head_drop = _first_finite(
            features,
            ("head_drop", "delta_head", "dh", "source_minus_target_head"),
        )
        if head_drop is None:
            head_drop = _node_head_drop(channel_data, edge)
        gradient = _first_finite(
            features,
            ("gradient", "hydraulic_gradient", "head_gradient"),
        )
        distance = _positive_or_none(
            _first_finite(features, ("distance", "length", "edge_length"))
        )
        if gradient is None and head_drop is not None:
            gradient = head_drop / distance if distance else head_drop
        conductance = _first_finite(
            features,
            ("conductance", "transmissivity", "hydraulic_conductivity"),
        )
        if gradient is None:
            probability = 0.5
            score = None
            reason = "missing_hydraulic_gradient"
        else:
            gradient_scale = 0.05
            conductance_weight = 0.10
            conductance_term = 0.0
            if conductance is not None and conductance > 0:
                conductance_term = conductance_weight * math.tanh(
                    math.log1p(conductance)
                )
            score = (gradient / gradient_scale) + conductance_term
            probability = _logistic(score)
            reason = "down_gradient_hydraulic_support"
        predictions.append(
            _prediction(
                edge,
                probability,
                channel=channel,
                reason=reason,
                score=score,
                abstention=abstention,
                cost=cost,
            )
        )
    return tuple(predictions)


def _score_edge_local_topology(
    candidate_universe: Iterable[Any],
    observations: ObservationMap,
    *,
    channel: str,
    abstention: Mapping[str, Any],
    cost: Mapping[str, float],
) -> tuple[BaselinePrediction, ...]:
    channel_data = observations.get(channel, {})
    predictions: list[BaselinePrediction] = []
    for candidate in candidate_universe:
        edge = _normalise_edge(candidate)
        features = _edge_features(channel_data, edge)
        probability = _first_finite(
            features,
            ("probability", "p", "local_probability", "posterior"),
        )
        score = _first_finite(
            features,
            (
                "local_support",
                "support",
                "similarity",
                "connectivity_prior",
                "candidate_weight",
                "weight",
            ),
        )
        if probability is None and score is not None:
            if 0.0 <= score <= 1.0:
                probability = score
            else:
                probability = _logistic(float(score))
        if probability is None:
            probability = 0.5
            reason = "missing_edge_local_topology_support"
        else:
            probability = _clamp_probability(probability)
            reason = "edge_local_topology_support"
        predictions.append(
            _prediction(
                edge,
                probability,
                channel=channel,
                reason=reason,
                score=score,
                abstention=abstention,
                cost=cost,
            )
        )
    return tuple(predictions)


def _prediction(
    edge: Edge,
    probability: float,
    *,
    channel: str,
    reason: str,
    score: float | None,
    abstention: Mapping[str, Any],
    cost: Mapping[str, float],
) -> BaselinePrediction:
    probability = _clamp_probability(probability)
    decision = _decision_from_probability(
        probability,
        select_threshold=float(abstention["select_threshold"]),
        reject_threshold=float(abstention["reject_threshold"]),
    )
    return BaselinePrediction(
        source=edge[0],
        target=edge[1],
        probability=probability,
        decision=decision,
        evidence_channel=channel,
        reason=reason,
        score=score,
        cost=dict(cost),
    )


def _decision_from_probability(
    probability: float,
    *,
    select_threshold: float,
    reject_threshold: float,
) -> Decision:
    if probability >= select_threshold:
        return SELECT
    if probability <= reject_threshold:
        return REJECT
    return ABSTAIN


def _normalise_edge(edge: Any) -> Edge:
    if hasattr(edge, "u") and hasattr(edge, "v"):
        return str(getattr(edge, "u")), str(getattr(edge, "v"))
    if hasattr(edge, "source") and hasattr(edge, "target"):
        return str(getattr(edge, "source")), str(getattr(edge, "target"))
    if isinstance(edge, Mapping):
        for left_key, right_key in (
            ("source", "target"),
            ("u", "v"),
            ("from", "to"),
            ("src", "dst"),
        ):
            if left_key in edge and right_key in edge:
                return str(edge[left_key]), str(edge[right_key])
        raise ValueError(f"Candidate edge mapping lacks source/target keys: {edge!r}")
    values = tuple(edge)
    if len(values) < 2:
        raise ValueError(f"Candidate edge must contain at least two values: {edge!r}")
    return str(values[0]), str(values[1])


def _edge_features(channel_data: Any, edge: Edge) -> Mapping[str, Any]:
    if isinstance(channel_data, Mapping):
        edge_maps = []
        for key in ("edges", "edge_features", "candidate_edges"):
            value = channel_data.get(key)
            if isinstance(value, Mapping):
                edge_maps.append(value)
        edge_maps.append(channel_data)
        for edge_map in edge_maps:
            feature = _lookup_edge_value(edge_map, edge)
            if isinstance(feature, Mapping):
                return feature
            if feature is not None:
                return {"support": feature}
    return {}


def _lookup_edge_value(edge_map: Mapping[Any, Any], edge: Edge) -> Any:
    keys = (
        edge,
        list(edge),
        f"{edge[0]}->{edge[1]}",
        f"{edge[0]},{edge[1]}",
        f"{edge[0]}|{edge[1]}",
    )
    for key in keys:
        try:
            if key in edge_map:
                return edge_map[key]  # type: ignore[index]
        except TypeError:
            continue
    return None


def _node_head_drop(channel_data: Any, edge: Edge) -> float | None:
    if not isinstance(channel_data, Mapping):
        return None
    node_heads = channel_data.get("node_heads") or channel_data.get("heads")
    if not isinstance(node_heads, Mapping):
        return None
    source_head = _coerce_float(node_heads.get(edge[0]))
    target_head = _coerce_float(node_heads.get(edge[1]))
    if source_head is None or target_head is None:
        return None
    return source_head - target_head


def _first_finite(features: Mapping[str, Any], keys: Iterable[str]) -> float | None:
    for key in keys:
        value = _coerce_float(features.get(key))
        if value is not None:
            return value
    return None


def _positive_or_none(value: float | None) -> float | None:
    return value if value is not None and value > 0 else None


def _coerce_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _finite_probability(value: float) -> None:
    if not math.isfinite(float(value)) or not 0.0 <= float(value) <= 1.0:
        raise ValueError(f"Probability must be finite and in [0, 1]: {value!r}")


def _clamp_probability(value: float) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Probability must be finite: {value!r}")
    return max(0.0, min(1.0, number))


def _logistic(value: float) -> float:
    if value >= 0:
        exp_neg = math.exp(-value)
        return 1.0 / (1.0 + exp_neg)
    exp_pos = math.exp(value)
    return exp_pos / (1.0 + exp_pos)


def _jsonable(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True, default=str))


def _fingerprint(record: Mapping[str, Any]) -> str:
    without_fingerprint = {
        key: value for key, value in record.items() if key != "fingerprint"
    }
    payload = json.dumps(without_fingerprint, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


__all__ = [
    "ABSTAIN",
    "REJECT",
    "SELECT",
    "BaselinePrediction",
    "BaselineRegistry",
    "BaselineSpec",
    "assert_truth_blind_observations",
    "default_baseline_registry",
    "edge_local_topology_baseline_spec",
    "hydraulic_only_baseline_spec",
    "permutation_control_baseline",
]
