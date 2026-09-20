"""Independent, static virtual benchmark for source/mixing-aware graph TTDs.

The module is deliberately a *benchmark fixture*, not a change to HydroSheaf's
production graph-conditioning semantics.  Its forward model is implemented
here from elementary convolution and mixing operations and imports no
HydroSheaf inference module.  The simulator holds the true topology, stationary
edge kernels, local recharge fractions, exact forcing series, and noiseless
node signals in :class:`StaticTTDGraphTruth`.  The public
:class:`StaticTTDGraphObservations` object contains only noisy, irregular,
partly missing inputs and calibration observations plus timestamps (but not
values) for held-out targets.

That separation is intentional.  A recovery method receives observations only;
the scorer receives sealed truth only after the method has produced a serialised
submission.  The included reference estimator is truth-blind: it uses a
non-negative source/mixing regression for each declared candidate graph and
can return an explicit ``ABSTAIN`` result when the public data are inadequate.
It is a benchmark control, not a claim of a generally identified field inverse.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

SCHEMA_VERSION = "hydrosheaf-static-ttd-graph-benchmark-v1"
GENERATOR_FAMILY = "independent_analytic_source_mixing_network"
NODES: tuple[str, ...] = ("R", "A", "B", "M", "O")
TRUE_EDGES: tuple[tuple[str, str], ...] = (
    ("R", "A"),
    ("R", "B"),
    ("A", "M"),
    ("B", "M"),
    ("M", "O"),
)
TRACERS: tuple[str, ...] = ("delta18O", "delta2H")


def _edge_key(source: str, target: str) -> str:
    return f"{source}->{target}"


def _held_out_key(node_id: str, tracer: str, time_day: int) -> str:
    return f"{node_id}|{tracer}|{int(time_day)}"


def _canonical_sha256(payload: Mapping[str, Any] | Sequence[Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _source_sha256() -> str:
    source_path = Path(__file__)
    return hashlib.sha256(source_path.read_bytes()).hexdigest()


def _finite_float(value: float | np.floating[Any]) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Benchmark payloads must contain finite floats.")
    return result


def _topological_order(
    nodes: Iterable[str], edges: Iterable[tuple[str, str]]
) -> tuple[str, ...] | None:
    """Return a deterministic topological order, or ``None`` for a cycle."""

    node_list = tuple(dict.fromkeys(str(node) for node in nodes))
    incoming = {node: 0 for node in node_list}
    children: dict[str, list[str]] = {node: [] for node in node_list}
    for source, target in edges:
        if source not in incoming or target not in incoming or source == target:
            return None
        children[source].append(target)
        incoming[target] += 1
    ready = sorted(node for node, degree in incoming.items() if degree == 0)
    order: list[str] = []
    while ready:
        node = ready.pop(0)
        order.append(node)
        for target in sorted(children[node]):
            incoming[target] -= 1
            if incoming[target] == 0:
                ready.append(target)
                ready.sort()
    return tuple(order) if len(order) == len(node_list) else None


def _safe_ratio(numerator: int | float, denominator: int | float) -> float | None:
    if float(denominator) <= 0.0:
        return None
    return float(numerator) / float(denominator)


@dataclass(frozen=True)
class StaticTTDGraphProtocol:
    """Public, pre-registered benchmark and scoring protocol."""

    schema_version: str
    n_steps: int
    age_grid_days: tuple[int, ...]
    estimator_lag_grid_days: tuple[int, ...]
    young_water_cutoff_days: int
    tracers: tuple[str, ...]
    evaluation_nodes: tuple[str, ...]
    input_sampling_fraction: float
    output_sampling_fraction: float
    input_missing_fraction: float
    output_missing_fraction: float
    held_out_fraction: float
    estimator_min_samples: int
    estimator_min_r2: float
    estimator_smoothness: float
    selected_edge_mass_threshold: float
    nominal_interval_coverage: float

    def __post_init__(self) -> None:
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError("Unexpected static graph TTD benchmark schema version.")
        if self.n_steps < 100:
            raise ValueError("n_steps must be at least 100 for the declared benchmark.")
        if len(self.age_grid_days) < 2 or self.age_grid_days[0] != 0:
            raise ValueError(
                "age_grid_days must start at zero and have at least two bins."
            )
        if tuple(sorted(set(self.age_grid_days))) != self.age_grid_days:
            raise ValueError("age_grid_days must be strictly increasing.")
        if not self.estimator_lag_grid_days or self.estimator_lag_grid_days[0] != 0:
            raise ValueError("estimator_lag_grid_days must start at zero.")
        if any(lag not in self.age_grid_days for lag in self.estimator_lag_grid_days):
            raise ValueError("Estimator lags must lie on the public age grid.")
        if self.young_water_cutoff_days not in self.age_grid_days:
            raise ValueError("young_water_cutoff_days must lie on the age grid.")
        fractions = (
            self.input_sampling_fraction,
            self.output_sampling_fraction,
            self.input_missing_fraction,
            self.output_missing_fraction,
            self.held_out_fraction,
            self.nominal_interval_coverage,
        )
        if any(not 0.0 < float(value) < 1.0 for value in fractions):
            raise ValueError(
                "Declared sampling, missingness, and coverage fractions need be in (0, 1)."
            )
        if self.estimator_min_samples < 8:
            raise ValueError("estimator_min_samples must be at least eight.")
        if not 0.0 <= self.estimator_min_r2 < 1.0:
            raise ValueError("estimator_min_r2 must be in [0, 1).")
        if self.estimator_smoothness < 0.0:
            raise ValueError("estimator_smoothness must be non-negative.")
        if not 0.0 < self.selected_edge_mass_threshold < 1.0:
            raise ValueError("selected_edge_mass_threshold must be in (0, 1).")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "n_steps": int(self.n_steps),
            "age_grid_days": [int(value) for value in self.age_grid_days],
            "estimator_lag_grid_days": [
                int(value) for value in self.estimator_lag_grid_days
            ],
            "young_water_cutoff_days": int(self.young_water_cutoff_days),
            "tracers": list(self.tracers),
            "evaluation_nodes": list(self.evaluation_nodes),
            "input_sampling_fraction": float(self.input_sampling_fraction),
            "output_sampling_fraction": float(self.output_sampling_fraction),
            "input_missing_fraction": float(self.input_missing_fraction),
            "output_missing_fraction": float(self.output_missing_fraction),
            "held_out_fraction": float(self.held_out_fraction),
            "estimator_min_samples": int(self.estimator_min_samples),
            "estimator_min_r2": float(self.estimator_min_r2),
            "estimator_smoothness": float(self.estimator_smoothness),
            "selected_edge_mass_threshold": float(self.selected_edge_mass_threshold),
            "nominal_interval_coverage": float(self.nominal_interval_coverage),
            "primary_estimand": "node_young_water_fraction",
            "primary_estimand_definition": (
                "Probability mass of the node-level stationary TTD at ages no greater "
                "than young_water_cutoff_days."
            ),
            "pre_registered_metrics": [
                "conditional_interval_coverage",
                "mean_interval_width",
                "explicit_abstention_rate",
                "held_out_prediction_mae",
                "held_out_prediction_rmse",
                "candidate_and_selected_topology_precision_recall_f1",
            ],
        }


def default_static_ttd_graph_protocol(*, n_steps: int = 300) -> StaticTTDGraphProtocol:
    """Return the fixed public protocol used by the default virtual case."""

    return StaticTTDGraphProtocol(
        schema_version=SCHEMA_VERSION,
        n_steps=int(n_steps),
        age_grid_days=tuple(range(0, 241)),
        estimator_lag_grid_days=tuple(range(0, 81, 4)),
        # Thirty days deliberately separates the branch/merge age structure;
        # at sixty days this particular small network would be nearly all
        # "young" and would make the declared estimand uninformative.
        young_water_cutoff_days=30,
        tracers=TRACERS,
        evaluation_nodes=("A", "B", "M", "O"),
        input_sampling_fraction=0.78,
        output_sampling_fraction=0.76,
        input_missing_fraction=0.11,
        output_missing_fraction=0.12,
        held_out_fraction=0.25,
        estimator_min_samples=60,
        estimator_min_r2=0.05,
        estimator_smoothness=0.035,
        selected_edge_mass_threshold=0.055,
        nominal_interval_coverage=0.90,
    )


@dataclass(frozen=True)
class CandidateGraphControl:
    """A public candidate topology used as a pre-registered control arm."""

    control_id: str
    description: str
    edges: tuple[tuple[str, str], ...]

    def __post_init__(self) -> None:
        if not self.control_id.strip() or not self.description.strip():
            raise ValueError(
                "Candidate graph controls need non-empty identifiers and descriptions."
            )
        if len(set(self.edges)) != len(self.edges):
            raise ValueError("Candidate graph control edges must be unique.")
        for source, target in self.edges:
            if source not in NODES or target not in NODES or source == target:
                raise ValueError(
                    "Candidate graph control has an invalid directed edge."
                )

    def to_dict(self) -> dict[str, Any]:
        return {
            "control_id": self.control_id,
            "description": self.description,
            "edges": [[source, target] for source, target in self.edges],
        }


@dataclass(frozen=True)
class InputObservationRow:
    """One noisy or missing public local-recharge tracer input observation."""

    node_id: str
    tracer: str
    time_day: int
    value: float | None
    sigma: float
    observed: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "tracer": self.tracer,
            "time_day": int(self.time_day),
            "value": None if self.value is None else _finite_float(self.value),
            "sigma": _finite_float(self.sigma),
            "observed": bool(self.observed),
        }


@dataclass(frozen=True)
class CalibrationObservationRow:
    """One public output observation available to a recovery method."""

    node_id: str
    tracer: str
    time_day: int
    value: float | None
    sigma: float
    observed: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "tracer": self.tracer,
            "time_day": int(self.time_day),
            "value": None if self.value is None else _finite_float(self.value),
            "sigma": _finite_float(self.sigma),
            "observed": bool(self.observed),
            "split": "calibration",
        }


@dataclass(frozen=True)
class HeldOutTarget:
    """A public target timestamp whose value remains in the sealed truth payload."""

    node_id: str
    tracer: str
    time_day: int
    sigma: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "tracer": self.tracer,
            "time_day": int(self.time_day),
            "sigma": _finite_float(self.sigma),
            "split": "held_out",
            "value": None,
        }


@dataclass(frozen=True)
class EdgeKernelTruth:
    """A sealed, normalized stationary edge kernel on the declared age grid."""

    source: str
    target: str
    masses: tuple[float, ...]

    def __post_init__(self) -> None:
        values = np.asarray(self.masses, dtype=float)
        if (
            self.source not in NODES
            or self.target not in NODES
            or self.source == self.target
        ):
            raise ValueError("Sealed edge kernel has an invalid directed edge.")
        if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
            raise ValueError("Sealed edge kernel masses must be a finite vector.")
        if np.any(values < -1.0e-12) or not np.isclose(
            float(values.sum()), 1.0, atol=1.0e-9
        ):
            raise ValueError(
                "Sealed edge kernel masses must be non-negative and normalized."
            )

    @property
    def edge_key(self) -> str:
        return _edge_key(self.source, self.target)

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "target": self.target,
            "masses": [_finite_float(value) for value in self.masses],
        }


@dataclass(frozen=True)
class StaticTTDGraphObservations:
    """Serializable public payload with no latent topology or TTD truth.

    Candidate graph controls are intentionally *not* embedded here.  The
    benchmark harness supplies a declared candidate graph as a separate input
    to the estimator, exactly as a field workflow would receive an externally
    proposed topology.  This keeps an oracle ``correct`` control from leaking
    the simulator's latent topology into the public observation artifact.
    """

    benchmark_id: str
    seed: int
    protocol: StaticTTDGraphProtocol
    input_rows: tuple[InputObservationRow, ...]
    calibration_rows: tuple[CalibrationObservationRow, ...]
    held_out_targets: tuple[HeldOutTarget, ...]
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "provenance", dict(self.provenance))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "benchmark_id": self.benchmark_id,
            "seed": int(self.seed),
            "protocol": self.protocol.to_dict(),
            "input_rows": [row.to_dict() for row in self.input_rows],
            "calibration_rows": [row.to_dict() for row in self.calibration_rows],
            "held_out_targets": [target.to_dict() for target in self.held_out_targets],
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class StaticTTDGraphTruth:
    """Serializable sealed truth.  Do not pass this to a recovery method."""

    benchmark_id: str
    seed: int
    true_edges: tuple[tuple[str, str], ...]
    edge_kernels: tuple[EdgeKernelTruth, ...]
    node_mixing_weights: Mapping[str, Mapping[str, float]]
    local_input_signals: Mapping[str, Mapping[str, tuple[float, ...]]]
    node_output_signals: Mapping[str, Mapping[str, tuple[float, ...]]]
    node_ttd_masses: Mapping[str, tuple[float, ...]]
    young_water_fractions: Mapping[str, float]
    held_out_values: Mapping[str, float]
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "node_mixing_weights",
            {node: dict(weights) for node, weights in self.node_mixing_weights.items()},
        )
        object.__setattr__(
            self,
            "local_input_signals",
            {
                node: {
                    tracer: tuple(float(value) for value in values)
                    for tracer, values in by_tracer.items()
                }
                for node, by_tracer in self.local_input_signals.items()
            },
        )
        object.__setattr__(
            self,
            "node_output_signals",
            {
                node: {
                    tracer: tuple(float(value) for value in values)
                    for tracer, values in by_tracer.items()
                }
                for node, by_tracer in self.node_output_signals.items()
            },
        )
        object.__setattr__(
            self,
            "node_ttd_masses",
            {
                node: tuple(float(value) for value in masses)
                for node, masses in self.node_ttd_masses.items()
            },
        )
        object.__setattr__(
            self, "young_water_fractions", dict(self.young_water_fractions)
        )
        object.__setattr__(self, "held_out_values", dict(self.held_out_values))
        object.__setattr__(self, "provenance", dict(self.provenance))

        if tuple(self.true_edges) != TRUE_EDGES:
            raise ValueError(
                "This fixed virtual benchmark expects the declared branch/merge truth edges."
            )
        if set(self.node_mixing_weights) != set(NODES):
            raise ValueError(
                "Sealed truth needs mixing weights for every benchmark node."
            )
        for node, weights in self.node_mixing_weights.items():
            values = np.asarray(list(weights.values()), dtype=float)
            if not np.all(np.isfinite(values)) or np.any(values < -1.0e-12):
                raise ValueError(
                    f"Sealed mixing weights at {node} must be finite and non-negative."
                )
            if not np.isclose(float(values.sum()), 1.0, atol=1.0e-9):
                raise ValueError(f"Sealed mixing weights at {node} must sum to one.")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "benchmark_id": self.benchmark_id,
            "seed": int(self.seed),
            "true_edges": [[source, target] for source, target in self.true_edges],
            "edge_kernels": [kernel.to_dict() for kernel in self.edge_kernels],
            "node_mixing_weights": {
                node: {key: _finite_float(value) for key, value in weights.items()}
                for node, weights in self.node_mixing_weights.items()
            },
            "local_input_signals": {
                node: {
                    tracer: [_finite_float(value) for value in values]
                    for tracer, values in by_tracer.items()
                }
                for node, by_tracer in self.local_input_signals.items()
            },
            "node_output_signals": {
                node: {
                    tracer: [_finite_float(value) for value in values]
                    for tracer, values in by_tracer.items()
                }
                for node, by_tracer in self.node_output_signals.items()
            },
            "node_ttd_masses": {
                node: [_finite_float(value) for value in masses]
                for node, masses in self.node_ttd_masses.items()
            },
            "young_water_fractions": {
                node: _finite_float(value)
                for node, value in self.young_water_fractions.items()
            },
            "held_out_values": {
                key: _finite_float(value) for key, value in self.held_out_values.items()
            },
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class StaticTTDGraphCase:
    """One case with public observations, sealed truth, and a harness control plan."""

    observations: StaticTTDGraphObservations
    truth: StaticTTDGraphTruth
    controls: tuple[CandidateGraphControl, ...]
    manifest: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "manifest", dict(self.manifest))
        if self.observations.benchmark_id != self.truth.benchmark_id:
            raise ValueError("Observation and truth benchmark identifiers must match.")
        if self.observations.seed != self.truth.seed:
            raise ValueError("Observation and truth seeds must match.")
        _validate_controls(self.controls)

    def control(self, control_id: str) -> CandidateGraphControl:
        for control in self.controls:
            if control.control_id == control_id:
                return control
        raise KeyError(f"Unknown benchmark-harness control {control_id!r}.")

    def verify_manifest(self) -> bool:
        expected = _case_manifest(self.observations, self.truth, self.controls)
        return dict(self.manifest) == expected


@dataclass(frozen=True)
class YoungWaterInterval:
    """One declared interval or an explicit abstention for a node estimand."""

    status: str
    point: float | None
    lower: float | None
    upper: float | None
    reason: str | None = None
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.status not in {"ESTIMATED", "ABSTAIN"}:
            raise ValueError(
                "Young-water interval status must be ESTIMATED or ABSTAIN."
            )
        if self.status == "ESTIMATED":
            values = (self.point, self.lower, self.upper)
            if any(
                value is None or not math.isfinite(float(value)) for value in values
            ):
                raise ValueError(
                    "An estimated interval needs finite point, lower, and upper values."
                )
            if (
                not 0.0
                <= float(self.lower)
                <= float(self.point)
                <= float(self.upper)
                <= 1.0
            ):
                raise ValueError(
                    "An estimated young-water interval must lie in [0, 1]."
                )
        elif self.reason is None:
            raise ValueError("An abstention needs a machine-readable reason.")
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "point": None if self.point is None else _finite_float(self.point),
            "lower": None if self.lower is None else _finite_float(self.lower),
            "upper": None if self.upper is None else _finite_float(self.upper),
            "reason": self.reason,
            "diagnostics": dict(self.diagnostics),
        }


@dataclass(frozen=True)
class HeldOutPrediction:
    """One truth-blind prediction submitted for a public held-out timestamp."""

    node_id: str
    tracer: str
    time_day: int
    value: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "tracer": self.tracer,
            "time_day": int(self.time_day),
            "value": _finite_float(self.value),
        }


@dataclass(frozen=True)
class StaticTTDGraphSubmission:
    """Serializable, truth-blind recovery submission for one control arm."""

    method_id: str
    control_id: str
    candidate_edges: tuple[tuple[str, str], ...]
    node_intervals: Mapping[str, YoungWaterInterval]
    held_out_predictions: tuple[HeldOutPrediction, ...]
    selected_edges: tuple[tuple[str, str], ...]
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "node_intervals", dict(self.node_intervals))
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))
        if not self.method_id.strip() or not self.control_id.strip():
            raise ValueError(
                "A benchmark submission needs method and control identifiers."
            )
        if len(set(self.candidate_edges)) != len(self.candidate_edges):
            raise ValueError("Submitted candidate edges must be unique.")
        if len(set(self.selected_edges)) != len(self.selected_edges):
            raise ValueError("Submitted selected edges must be unique.")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "method_id": self.method_id,
            "control_id": self.control_id,
            "candidate_edges": [
                [source, target] for source, target in self.candidate_edges
            ],
            "node_intervals": {
                node: interval.to_dict()
                for node, interval in self.node_intervals.items()
            },
            "held_out_predictions": [
                prediction.to_dict() for prediction in self.held_out_predictions
            ],
            "selected_edges": [
                [source, target] for source, target in self.selected_edges
            ],
            "diagnostics": dict(self.diagnostics),
            "truth_used": False,
        }


@dataclass(frozen=True)
class StaticTTDGraphScore:
    """Sealed-truth score under the public pre-registered protocol."""

    benchmark_id: str
    method_id: str
    control_id: str
    metrics: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "metrics", dict(self.metrics))

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "benchmark_id": self.benchmark_id,
            "method_id": self.method_id,
            "control_id": self.control_id,
            "metrics": dict(self.metrics),
        }


@dataclass(frozen=True)
class _TracerFit:
    """Internal data-only fit used by the benchmark control estimator."""

    status: str
    reason: str | None
    node_id: str
    tracer: str
    local_mass: float | None
    edge_lag_masses: Mapping[str, tuple[float, ...]]
    r2: float | None
    n_observed: int
    n_features: int
    condition_number: float | None
    predicted_grid: tuple[float, ...] | None
    diagnostics: Mapping[str, Any]


def _default_controls() -> tuple[CandidateGraphControl, ...]:
    """Return benchmark-harness controls, including a labelled oracle arm.

    These controls are a separate experiment-plan artifact rather than public
    observations.  In particular, ``correct`` is an oracle sensitivity control,
    not a graph that a blind recovery method may discover from a case payload.
    """

    return (
        CandidateGraphControl(
            "correct",
            "Declared true branch/merge topology; a positive-control topology arm.",
            TRUE_EDGES,
        ),
        CandidateGraphControl(
            "reversed",
            "All true edge directions reversed; tests direction-sensitive recovery.",
            tuple((target, source) for source, target in TRUE_EDGES),
        ),
        CandidateGraphControl(
            "random",
            "Deterministic wrong DAG with no true directed edge.",
            (("R", "M"), ("R", "O"), ("B", "A"), ("A", "O"), ("B", "O")),
        ),
        CandidateGraphControl(
            "edge_removed",
            "Correct topology with the B-to-M branch deliberately removed.",
            tuple(edge for edge in TRUE_EDGES if edge != ("B", "M")),
        ),
        CandidateGraphControl(
            "local",
            "Local-only baseline: no declared inter-node transport edges.",
            (),
        ),
    )


def _validate_controls(controls: Sequence[CandidateGraphControl]) -> None:
    ids = [control.control_id for control in controls]
    if len(ids) != len(set(ids)):
        raise ValueError(
            "Benchmark-harness candidate-control identifiers must be unique."
        )
    required = {"correct", "reversed", "random", "edge_removed", "local"}
    missing = required.difference(ids)
    if missing:
        raise ValueError(f"Missing required pre-registered controls: {sorted(missing)}")


def _gamma_like_kernel(
    age_grid_days: Sequence[int], shape: float, scale: float
) -> np.ndarray:
    """Construct a smooth positive kernel without using an inverse-model family."""

    age = np.asarray(age_grid_days, dtype=float)
    support = np.minimum(age, 80.0)
    raw = np.power(support + 0.35, shape - 1.0) * np.exp(-support / scale)
    raw[age > 80.0] = 0.0
    raw[0] *= 0.20
    if float(raw.sum()) <= 0.0:
        raise RuntimeError("Independent virtual kernel construction failed.")
    return raw / float(raw.sum())


def _true_edge_kernels(protocol: StaticTTDGraphProtocol) -> tuple[EdgeKernelTruth, ...]:
    parameters = {
        ("R", "A"): (2.25, 4.5),
        ("R", "B"): (3.10, 6.4),
        ("A", "M"): (2.05, 5.2),
        ("B", "M"): (3.35, 6.0),
        ("M", "O"): (2.65, 5.0),
    }
    return tuple(
        EdgeKernelTruth(
            source,
            target,
            tuple(
                _gamma_like_kernel(
                    protocol.age_grid_days, *parameters[(source, target)]
                )
            ),
        )
        for source, target in TRUE_EDGES
    )


def _mixing_weights() -> dict[str, dict[str, float]]:
    """True source contributions; every non-root node has local recharge."""

    return {
        "R": {"local": 1.0},
        "A": {"local": 0.18, "R": 0.82},
        "B": {"local": 0.29, "R": 0.71},
        "M": {"local": 0.13, "A": 0.47, "B": 0.40},
        "O": {"local": 0.08, "M": 0.92},
    }


def _make_local_input_signals(
    protocol: StaticTTDGraphProtocol,
    seed: int,
) -> dict[str, dict[str, tuple[float, ...]]]:
    """Generate exact source/recharge forcing over a warm-up plus observation window."""

    rng = np.random.default_rng(int(seed) * 16127 + 97)
    warmup = int(protocol.age_grid_days[-1]) + 4
    time = np.arange(-warmup, protocol.n_steps, dtype=float)
    values: dict[str, dict[str, tuple[float, ...]]] = {}
    for node_index, node in enumerate(NODES):
        phase = 0.77 * node_index + rng.uniform(-0.14, 0.14)
        annual = np.sin(2.0 * np.pi * time / 47.0 + phase)
        short = np.cos(2.0 * np.pi * time / 19.0 + 0.43 * node_index)
        pulse_a = np.exp(-(((time - (57.0 + 13.0 * node_index)) / 7.0) ** 2))
        pulse_b = np.exp(-(((time - (170.0 - 9.0 * node_index)) / 11.0) ** 2))
        d18o = (
            -8.3
            + 0.43 * node_index
            + 1.28 * annual
            + 0.42 * short
            + (0.72 - 0.06 * node_index) * pulse_a
            - 0.38 * pulse_b
        )
        d2h = 8.0 * d18o + 10.4 + 0.72 * np.cos(2.0 * np.pi * time / 31.0 + phase)
        values[node] = {
            "delta18O": tuple(float(value) for value in d18o),
            "delta2H": tuple(float(value) for value in d2h),
        }
    return values


def _forward_node_signals(
    local_inputs: Mapping[str, Mapping[str, Sequence[float]]],
    mixing: Mapping[str, Mapping[str, float]],
    kernels: Mapping[str, Sequence[float]],
) -> dict[str, dict[str, tuple[float, ...]]]:
    """Forward the exact stationary convolution/mixing model on a DAG."""

    order = _topological_order(NODES, TRUE_EDGES)
    if order is None:
        raise RuntimeError("The fixed virtual truth graph must be acyclic.")
    outputs: dict[str, dict[str, tuple[float, ...]]] = {}
    parents: dict[str, list[str]] = {node: [] for node in NODES}
    for source, target in TRUE_EDGES:
        parents[target].append(source)
    for node in order:
        outputs[node] = {}
        for tracer in TRACERS:
            local = np.asarray(local_inputs[node][tracer], dtype=float)
            signal = float(mixing[node]["local"]) * local
            for source in parents[node]:
                edge_mass = float(mixing[node][source])
                kernel = np.asarray(kernels[_edge_key(source, node)], dtype=float)
                source_signal = np.asarray(outputs[source][tracer], dtype=float)
                signal = (
                    signal
                    + edge_mass
                    * np.convolve(source_signal, kernel, mode="full")[: local.size]
                )
            outputs[node][tracer] = tuple(float(value) for value in signal)
    return outputs


def _forward_node_ttds(
    mixing: Mapping[str, Mapping[str, float]],
    kernels: Mapping[str, Sequence[float]],
    age_size: int,
) -> dict[str, tuple[float, ...]]:
    """Compute the exact node-level TTD implied by serial convolution and mixing."""

    order = _topological_order(NODES, TRUE_EDGES)
    if order is None:
        raise RuntimeError("The fixed virtual truth graph must be acyclic.")
    parents: dict[str, list[str]] = {node: [] for node in NODES}
    for source, target in TRUE_EDGES:
        parents[target].append(source)
    node_ttds: dict[str, tuple[float, ...]] = {}
    delta = np.zeros((age_size,), dtype=float)
    delta[0] = 1.0
    for node in order:
        ttd = float(mixing[node]["local"]) * delta
        for source in parents[node]:
            transferred = np.convolve(
                np.asarray(node_ttds[source], dtype=float),
                np.asarray(kernels[_edge_key(source, node)], dtype=float),
                mode="full",
            )[:age_size]
            ttd = ttd + float(mixing[node][source]) * transferred
        total = float(ttd.sum())
        if total <= 0.0:
            raise RuntimeError(f"Non-positive truth TTD mass at node {node}.")
        node_ttds[node] = tuple(float(value) for value in (ttd / total))
    return node_ttds


def _sample_irregular_times(
    rng: np.random.Generator, n_steps: int, fraction: float
) -> np.ndarray:
    n = max(2, min(n_steps, int(round(float(n_steps) * float(fraction)))))
    return np.sort(rng.choice(np.arange(n_steps, dtype=int), size=n, replace=False))


def _input_sigma(tracer: str) -> float:
    return 0.075 if tracer == "delta18O" else 0.58


def _output_sigma(tracer: str) -> float:
    return 0.095 if tracer == "delta18O" else 0.72


def _build_public_observations(
    protocol: StaticTTDGraphProtocol,
    seed: int,
    local_inputs: Mapping[str, Mapping[str, Sequence[float]]],
    node_outputs: Mapping[str, Mapping[str, Sequence[float]]],
) -> tuple[
    tuple[InputObservationRow, ...],
    tuple[CalibrationObservationRow, ...],
    tuple[HeldOutTarget, ...],
    dict[str, float],
]:
    """Sample noisy public data and retain exact held-out values separately."""

    rng = np.random.default_rng(int(seed) * 31337 + 101)
    warmup = (
        len(next(iter(next(iter(local_inputs.values())).values()))) - protocol.n_steps
    )
    input_rows: list[InputObservationRow] = []
    calibration_rows: list[CalibrationObservationRow] = []
    held_out_targets: list[HeldOutTarget] = []
    held_out_values: dict[str, float] = {}

    for node in NODES:
        for tracer in protocol.tracers:
            sigma_input = _input_sigma(tracer)
            input_times = _sample_irregular_times(
                rng, protocol.n_steps, protocol.input_sampling_fraction
            )
            input_signal = np.asarray(local_inputs[node][tracer], dtype=float)[warmup:]
            for time_day in input_times:
                missing = bool(rng.random() < protocol.input_missing_fraction)
                value = (
                    None
                    if missing
                    else float(input_signal[time_day] + rng.normal(0.0, sigma_input))
                )
                input_rows.append(
                    InputObservationRow(
                        node, tracer, int(time_day), value, sigma_input, not missing
                    )
                )

            sigma_output = _output_sigma(tracer)
            output_times = _sample_irregular_times(
                rng, protocol.n_steps, protocol.output_sampling_fraction
            )
            output_signal = np.asarray(node_outputs[node][tracer], dtype=float)[warmup:]
            for time_day in output_times:
                missing = bool(rng.random() < protocol.output_missing_fraction)
                if missing:
                    calibration_rows.append(
                        CalibrationObservationRow(
                            node, tracer, int(time_day), None, sigma_output, False
                        )
                    )
                elif rng.random() < protocol.held_out_fraction:
                    held_out_targets.append(
                        HeldOutTarget(node, tracer, int(time_day), sigma_output)
                    )
                    held_out_values[_held_out_key(node, tracer, int(time_day))] = float(
                        output_signal[time_day]
                    )
                else:
                    value = float(
                        output_signal[time_day] + rng.normal(0.0, sigma_output)
                    )
                    calibration_rows.append(
                        CalibrationObservationRow(
                            node, tracer, int(time_day), value, sigma_output, True
                        )
                    )
    return (
        tuple(input_rows),
        tuple(calibration_rows),
        tuple(held_out_targets),
        held_out_values,
    )


def _case_manifest(
    observations: StaticTTDGraphObservations,
    truth: StaticTTDGraphTruth,
    controls: Sequence[CandidateGraphControl],
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "benchmark_id": observations.benchmark_id,
        "seed": int(observations.seed),
        "generator_family": GENERATOR_FAMILY,
        "generator_module": __name__,
        "generator_source_sha256": _source_sha256(),
        "observations_sha256": _canonical_sha256(observations.to_dict()),
        "sealed_truth_sha256": _canonical_sha256(truth.to_dict()),
        "control_plan_sha256": _canonical_sha256(
            {"controls": [control.to_dict() for control in controls]}
        ),
        "protocol_sha256": _canonical_sha256(observations.protocol.to_dict()),
        "truth_release_policy": "sealed_until_submission_serialized",
        "truth_blind_reference_estimator": True,
        "provenance": {
            "simulator_imports_hydrosheaf_inference": False,
            "stationary_edge_kernels": True,
            "local_recharge_at_each_node": True,
            "branch_and_merge_topology": True,
            "irregular_missing_observations": True,
        },
    }


def generate_static_ttd_graph_case(
    seed: int = 20260917,
    *,
    n_steps: int | None = None,
    protocol: StaticTTDGraphProtocol | None = None,
) -> StaticTTDGraphCase:
    """Generate one deterministic branch/merge virtual graph-TTD case.

    The returned object keeps public observations and sealed truth separate. A
    recovery call receives ``case.observations`` plus an externally declared
    candidate graph, such as ``case.control("correct")`` in the benchmark
    harness.  Use :func:`score_static_ttd_graph_submission` only afterward with
    ``case.truth``.
    """

    if protocol is None:
        protocol = default_static_ttd_graph_protocol(
            n_steps=300 if n_steps is None else int(n_steps)
        )
    elif n_steps is not None and int(n_steps) != protocol.n_steps:
        raise ValueError(
            "n_steps cannot disagree with an explicitly supplied static benchmark protocol."
        )
    kernels = _true_edge_kernels(protocol)
    kernel_map = {kernel.edge_key: kernel.masses for kernel in kernels}
    mixing = _mixing_weights()
    local_inputs = _make_local_input_signals(protocol, int(seed))
    outputs = _forward_node_signals(local_inputs, mixing, kernel_map)
    node_ttds = _forward_node_ttds(mixing, kernel_map, len(protocol.age_grid_days))
    cutoff_index = protocol.age_grid_days.index(protocol.young_water_cutoff_days)
    young_fractions = {
        node: float(np.sum(np.asarray(masses, dtype=float)[: cutoff_index + 1]))
        for node, masses in node_ttds.items()
    }
    input_rows, calibration_rows, held_out_targets, held_out_values = (
        _build_public_observations(protocol, int(seed), local_inputs, outputs)
    )
    benchmark_id = f"static-ttd-graph-{int(seed)}"
    controls = _default_controls()
    observations = StaticTTDGraphObservations(
        benchmark_id=benchmark_id,
        seed=int(seed),
        protocol=protocol,
        input_rows=input_rows,
        calibration_rows=calibration_rows,
        held_out_targets=held_out_targets,
        provenance={
            "generator_family": GENERATOR_FAMILY,
            "public_truth_policy": "No topology, edge kernels, mixing fractions, exact signals, or held-out values are in this payload.",
            "forcing": "Noisy local recharge time series are public; exact inputs remain sealed.",
            "time_units": "days",
        },
    )
    truth = StaticTTDGraphTruth(
        benchmark_id=benchmark_id,
        seed=int(seed),
        true_edges=TRUE_EDGES,
        edge_kernels=kernels,
        node_mixing_weights=mixing,
        local_input_signals=local_inputs,
        node_output_signals=outputs,
        node_ttd_masses=node_ttds,
        young_water_fractions=young_fractions,
        held_out_values=held_out_values,
        provenance={
            "generator_family": GENERATOR_FAMILY,
            "truth_is_sealed": True,
            "forward_model": "local recharge plus non-negative upstream mixing and stationary serial convolution",
            "time_units": "days",
        },
    )
    return StaticTTDGraphCase(
        observations,
        truth,
        controls,
        _case_manifest(observations, truth, controls),
    )


def write_static_ttd_graph_case(
    case: StaticTTDGraphCase, directory: str | Path
) -> dict[str, Path]:
    """Write public, sealed-truth, control-plan, and manifest JSON artifacts."""

    if not case.verify_manifest():
        raise ValueError("Refusing to write a case with a failed provenance manifest.")
    output_dir = Path(directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "observations": output_dir / "ttd_graph_static_observations.json",
        "sealed_truth": output_dir / "ttd_graph_static_sealed_truth.json",
        "control_plan": output_dir / "ttd_graph_static_control_plan.json",
        "manifest": output_dir / "ttd_graph_static_manifest.json",
    }
    payloads = {
        "observations": case.observations.to_dict(),
        "sealed_truth": case.truth.to_dict(),
        "control_plan": {
            "schema_version": SCHEMA_VERSION,
            "benchmark_id": case.observations.benchmark_id,
            "controls": [control.to_dict() for control in case.controls],
            "oracle_control_note": (
                "The correct control is an oracle topology sensitivity arm and is "
                "not embedded in the public observation payload."
            ),
        },
        "manifest": dict(case.manifest),
    }
    for name, path in paths.items():
        path.write_text(
            json.dumps(payloads[name], ensure_ascii=True, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
    return paths


def _series_from_input_rows(
    observations: StaticTTDGraphObservations,
    node_id: str,
    tracer: str,
) -> np.ndarray:
    series = np.full((observations.protocol.n_steps,), np.nan, dtype=float)
    for row in observations.input_rows:
        if (
            row.node_id == node_id
            and row.tracer == tracer
            and row.observed
            and row.value is not None
        ):
            series[row.time_day] = float(row.value)
    return series


def _series_from_calibration_rows(
    observations: StaticTTDGraphObservations,
    node_id: str,
    tracer: str,
) -> np.ndarray:
    series = np.full((observations.protocol.n_steps,), np.nan, dtype=float)
    for row in observations.calibration_rows:
        if (
            row.node_id == node_id
            and row.tracer == tracer
            and row.observed
            and row.value is not None
        ):
            series[row.time_day] = float(row.value)
    return series


def _interpolate_public_series(series: np.ndarray) -> np.ndarray | None:
    observed = np.flatnonzero(np.isfinite(series))
    if observed.size < 2:
        return None
    return np.interp(
        np.arange(series.size, dtype=float), observed.astype(float), series[observed]
    ).astype(float)


def _lagged_feature(series: np.ndarray, lag: int) -> np.ndarray:
    output = np.empty_like(series, dtype=float)
    if lag <= 0:
        output[:] = series
    else:
        output[:lag] = series[0]
        output[lag:] = series[:-lag]
    return output


def _fit_nonnegative_mixing_model(
    observations: StaticTTDGraphObservations,
    control: CandidateGraphControl,
    node_id: str,
    tracer: str,
) -> _TracerFit:
    """Fit one public-only local-recharge plus incoming-edge mixing model.

    The response uses *only calibration output rows*.  Parent output and local
    recharge inputs are reconstructed from their own public calibration/input
    rows.  No truth object is accepted or read by this function.
    """

    protocol = observations.protocol
    raw_target = _series_from_calibration_rows(observations, node_id, tracer)
    target = _interpolate_public_series(raw_target)
    local = _interpolate_public_series(
        _series_from_input_rows(observations, node_id, tracer)
    )
    parents = tuple(
        source for source, target_node in control.edges if target_node == node_id
    )
    if target is None or local is None:
        return _TracerFit(
            "ABSTAIN",
            "insufficient_target_or_local_input",
            node_id,
            tracer,
            None,
            {},
            None,
            int(np.isfinite(raw_target).sum()),
            0,
            None,
            None,
            {},
        )

    parent_series: dict[str, np.ndarray] = {}
    for source in parents:
        source_series = _interpolate_public_series(
            _series_from_calibration_rows(observations, source, tracer)
        )
        if source_series is None:
            return _TracerFit(
                "ABSTAIN",
                f"insufficient_parent_calibration:{source}",
                node_id,
                tracer,
                None,
                {},
                None,
                int(np.isfinite(raw_target).sum()),
                0,
                None,
                None,
                {},
            )
        parent_series[source] = source_series

    columns = [local]
    edge_column_indices: dict[str, list[int]] = {}
    for source in parents:
        edge_key = _edge_key(source, node_id)
        edge_column_indices[edge_key] = []
        for lag in protocol.estimator_lag_grid_days:
            edge_column_indices[edge_key].append(len(columns))
            columns.append(_lagged_feature(parent_series[source], int(lag)))
    design = np.column_stack(columns)
    fit_mask = np.isfinite(raw_target) & np.all(np.isfinite(design), axis=1)
    n_observed = int(fit_mask.sum())
    n_features = int(design.shape[1])
    minimum = max(protocol.estimator_min_samples, n_features + 12)
    if n_observed < minimum:
        return _TracerFit(
            "ABSTAIN",
            "insufficient_calibration_samples",
            node_id,
            tracer,
            None,
            {},
            None,
            n_observed,
            n_features,
            None,
            None,
            {"minimum_required": int(minimum), "parents": list(parents)},
        )
    y = raw_target[fit_mask]
    x = design[fit_mask]
    if float(np.std(y)) < 1.0e-8:
        return _TracerFit(
            "ABSTAIN",
            "no_target_variation",
            node_id,
            tracer,
            None,
            {},
            None,
            n_observed,
            n_features,
            None,
            None,
            {},
        )

    # A non-negative constrained fit gives interpretable local and edge masses.
    # Sum-to-one and edge-kernel smoothness are penalties, not truth constraints.
    scale = max(float(np.std(y)) * math.sqrt(float(n_observed)), 1.0)
    augmented_x = [x, (2.25 * scale) * np.ones((1, n_features), dtype=float)]
    augmented_y = [y, np.asarray([2.25 * scale], dtype=float)]
    if protocol.estimator_smoothness > 0.0 and parents:
        rows: list[np.ndarray] = []
        for indices in edge_column_indices.values():
            for left, right in zip(indices, indices[1:]):
                row = np.zeros((n_features,), dtype=float)
                row[left] = -1.0
                row[right] = 1.0
                rows.append(row)
        if rows:
            penalty = math.sqrt(protocol.estimator_smoothness) * scale
            augmented_x.append(penalty * np.vstack(rows))
            augmented_y.append(np.zeros((len(rows),), dtype=float))
    try:
        from scipy.optimize import nnls

        coefficients, _ = nnls(np.vstack(augmented_x), np.concatenate(augmented_y))
    except (
        Exception
    ) as error:  # pragma: no cover - hard to force scipy failure portably
        return _TracerFit(
            "ABSTAIN",
            f"nonnegative_fit_failed:{type(error).__name__}",
            node_id,
            tracer,
            None,
            {},
            None,
            n_observed,
            n_features,
            None,
            None,
            {},
        )
    coefficient_sum = float(np.sum(coefficients))
    if coefficient_sum <= 1.0e-12:
        return _TracerFit(
            "ABSTAIN",
            "zero_nonnegative_mass",
            node_id,
            tracer,
            None,
            {},
            None,
            n_observed,
            n_features,
            None,
            None,
            {},
        )
    coefficients = np.asarray(coefficients, dtype=float) / coefficient_sum
    fitted = design @ coefficients
    fitted_at_observations = fitted[fit_mask]
    sse = float(np.sum((y - fitted_at_observations) ** 2))
    sst = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 - sse / max(sst, 1.0e-12)
    condition_number = float(np.linalg.cond(x))
    if not math.isfinite(condition_number):
        condition_number = float("inf")
    if r2 < protocol.estimator_min_r2:
        return _TracerFit(
            "ABSTAIN",
            "below_declared_r2_threshold",
            node_id,
            tracer,
            None,
            {},
            float(r2),
            n_observed,
            n_features,
            condition_number,
            None,
            {
                "r2_threshold": float(protocol.estimator_min_r2),
                "parents": list(parents),
            },
        )
    edge_masses = {
        edge_key: tuple(float(coefficients[index]) for index in indices)
        for edge_key, indices in edge_column_indices.items()
    }
    return _TracerFit(
        "ESTIMATED",
        None,
        node_id,
        tracer,
        float(coefficients[0]),
        edge_masses,
        float(r2),
        n_observed,
        n_features,
        condition_number,
        tuple(float(value) for value in fitted),
        {
            "parents": list(parents),
            "coefficient_sum_before_normalization": coefficient_sum,
            "interval_semantics": "diagnostic_not_a_confidence_guarantee",
        },
    )


def _compose_estimated_ttds(
    observations: StaticTTDGraphObservations,
    control: CandidateGraphControl,
    fits: Mapping[tuple[str, str], _TracerFit],
    tracer: str,
) -> tuple[dict[str, np.ndarray], str | None]:
    """Compose fitted edge masses on the candidate DAG into node-level TTDs."""

    order = _topological_order(NODES, control.edges)
    if order is None:
        return {}, "cyclic_candidate_graph"
    parents: dict[str, list[str]] = {node: [] for node in NODES}
    for source, target in control.edges:
        parents[target].append(source)
    age_size = len(observations.protocol.age_grid_days)
    lag_grid = observations.protocol.estimator_lag_grid_days
    delta = np.zeros((age_size,), dtype=float)
    delta[0] = 1.0
    results: dict[str, np.ndarray] = {}
    for node in order:
        fit = fits[(node, tracer)]
        if fit.status != "ESTIMATED" or fit.local_mass is None:
            return {}, f"upstream_fit_abstained:{node}"
        result = float(fit.local_mass) * delta
        for source in parents[node]:
            if source not in results:
                return {}, f"missing_parent_ttd:{source}"
            edge_key = _edge_key(source, node)
            lag_masses = fit.edge_lag_masses.get(edge_key)
            if lag_masses is None:
                return {}, f"missing_edge_fit:{edge_key}"
            kernel = np.zeros((age_size,), dtype=float)
            for lag, mass in zip(lag_grid, lag_masses):
                kernel[int(lag)] += float(mass)
            result = (
                result + np.convolve(results[source], kernel, mode="full")[:age_size]
            )
        total = float(result.sum())
        if total <= 1.0e-12:
            return {}, f"nonpositive_composed_mass:{node}"
        results[node] = result / total
    return results, None


def _estimate_interval(
    values: Sequence[float], fits: Sequence[_TracerFit]
) -> YoungWaterInterval:
    if not values or not fits:
        return YoungWaterInterval(
            "ABSTAIN", None, None, None, "no_valid_tracer_ttd", {}
        )
    point = float(np.mean(np.asarray(values, dtype=float)))
    r2s = np.asarray([float(fit.r2 or 0.0) for fit in fits], dtype=float)
    sample_ratios = np.asarray(
        [
            min(1.0, float(fit.n_observed) / max(1.0, 2.0 * float(fit.n_features)))
            for fit in fits
        ],
        dtype=float,
    )
    tracer_spread = float(np.std(np.asarray(values, dtype=float), ddof=0))
    # This is deliberately marked as a diagnostic interval.  The sealed scorer
    # measures its realised coverage rather than asserting a nominal guarantee.
    half_width = float(
        np.clip(
            0.09
            + 0.22 * (1.0 - float(np.clip(np.mean(r2s), 0.0, 1.0)))
            + 0.18 * (1.0 - float(np.mean(sample_ratios)))
            + 1.15 * tracer_spread,
            0.06,
            0.48,
        )
    )
    lower = max(0.0, point - half_width)
    upper = min(1.0, point + half_width)
    return YoungWaterInterval(
        "ESTIMATED",
        point,
        lower,
        upper,
        None,
        {
            "accepted_tracers": [fit.tracer for fit in fits],
            "mean_r2": float(np.mean(r2s)),
            "tracer_point_spread": tracer_spread,
            "interval_method": "pre_registered_residual_diagnostic_heuristic",
            "nominal_coverage_is_not_assumed": True,
        },
    )


def run_candidate_graph_estimator(
    observations: StaticTTDGraphObservations,
    control: CandidateGraphControl,
) -> StaticTTDGraphSubmission:
    """Run the truth-blind source/mixing-aware recovery control for one graph.

    The only data payload accepted is ``observations``. The candidate graph is
    separately supplied, as it would be from a hydrogeologic graph hypothesis.
    The function does not accept a :class:`StaticTTDGraphTruth`, a case object,
    or any hidden signal. It should therefore be suitable as a conservative
    common baseline when new inverse methods are evaluated with the same sealed
    scorer.
    """

    fits: dict[tuple[str, str], _TracerFit] = {}
    for tracer in observations.protocol.tracers:
        for node in NODES:
            fits[(node, tracer)] = _fit_nonnegative_mixing_model(
                observations, control, node, tracer
            )

    node_estimates: dict[str, YoungWaterInterval] = {}
    per_tracer_ttds: dict[str, dict[str, np.ndarray]] = {}
    composition_reasons: dict[str, str] = {}
    for tracer in observations.protocol.tracers:
        composed, reason = _compose_estimated_ttds(observations, control, fits, tracer)
        if reason is None:
            per_tracer_ttds[tracer] = composed
        else:
            composition_reasons[tracer] = reason
    cutoff_index = observations.protocol.age_grid_days.index(
        observations.protocol.young_water_cutoff_days
    )
    for node in observations.protocol.evaluation_nodes:
        values: list[float] = []
        supporting_fits: list[_TracerFit] = []
        for tracer, composed in per_tracer_ttds.items():
            if node in composed:
                values.append(float(np.sum(composed[node][: cutoff_index + 1])))
                supporting_fits.append(fits[(node, tracer)])
        if not values:
            reasons = sorted(
                {
                    fit.reason
                    for tracer in observations.protocol.tracers
                    for fit in [fits[(node, tracer)]]
                    if fit.reason
                }
            )
            if composition_reasons:
                reasons.extend(
                    f"composition:{value}"
                    for value in sorted(set(composition_reasons.values()))
                )
            node_estimates[node] = YoungWaterInterval(
                "ABSTAIN",
                None,
                None,
                None,
                ";".join(reasons) if reasons else "no_composed_ttd",
                {},
            )
        else:
            node_estimates[node] = _estimate_interval(values, supporting_fits)

    predictions: list[HeldOutPrediction] = []
    for target in observations.held_out_targets:
        fit = fits[(target.node_id, target.tracer)]
        if fit.status == "ESTIMATED" and fit.predicted_grid is not None:
            value = float(fit.predicted_grid[target.time_day])
            if math.isfinite(value):
                predictions.append(
                    HeldOutPrediction(
                        target.node_id, target.tracer, target.time_day, value
                    )
                )

    edge_masses: dict[str, list[float]] = {}
    for fit in fits.values():
        if fit.status != "ESTIMATED":
            continue
        for edge_key, lag_masses in fit.edge_lag_masses.items():
            edge_masses.setdefault(edge_key, []).append(float(sum(lag_masses)))
    selected_edges = tuple(
        tuple(edge_key.split("->", 1))
        for edge_key, masses in sorted(edge_masses.items())
        if float(np.mean(masses)) >= observations.protocol.selected_edge_mass_threshold
    )
    fit_summary = {
        _edge_key(node, tracer): {
            "status": fit.status,
            "reason": fit.reason,
            "r2": fit.r2,
            "n_observed": fit.n_observed,
            "n_features": fit.n_features,
            "condition_number": fit.condition_number,
        }
        for (node, tracer), fit in fits.items()
    }
    return StaticTTDGraphSubmission(
        method_id="truth_blind_nonnegative_source_mixing_graph_control_v1",
        control_id=control.control_id,
        candidate_edges=control.edges,
        node_intervals=node_estimates,
        held_out_predictions=tuple(predictions),
        selected_edges=selected_edges,
        diagnostics={
            "truth_used": False,
            "control_description": control.description,
            "fit_summary": fit_summary,
            "composition_reasons": composition_reasons,
            "selected_edge_mean_masses": {
                edge_key: float(np.mean(masses))
                for edge_key, masses in edge_masses.items()
            },
        },
    )


def run_local_baseline(
    observations: StaticTTDGraphObservations,
) -> StaticTTDGraphSubmission:
    """Run the public local-only benchmark baseline without any graph edges."""

    return run_candidate_graph_estimator(
        observations,
        CandidateGraphControl(
            "local", "Local-only baseline: no declared inter-node transport edges.", ()
        ),
    )


def run_all_pre_registered_controls(
    observations: StaticTTDGraphObservations,
    controls: Sequence[CandidateGraphControl],
) -> dict[str, StaticTTDGraphSubmission]:
    """Run a supplied control plan without accessing sealed truth."""

    _validate_controls(controls)

    return {
        control.control_id: run_candidate_graph_estimator(observations, control)
        for control in controls
    }


def _topology_metrics(
    predicted_edges: Iterable[tuple[str, str]],
    true_edges: Iterable[tuple[str, str]],
) -> dict[str, Any]:
    predicted = set(tuple(edge) for edge in predicted_edges)
    truth = set(tuple(edge) for edge in true_edges)
    true_positive = len(predicted.intersection(truth))
    false_positive = len(predicted.difference(truth))
    false_negative = len(truth.difference(predicted))
    precision = _safe_ratio(true_positive, true_positive + false_positive)
    recall = _safe_ratio(true_positive, true_positive + false_negative)
    f1 = (
        None
        if precision is None or recall is None or precision + recall <= 0.0
        else 2.0 * precision * recall / (precision + recall)
    )
    unoriented_overlap = sum(
        1
        for source, target in predicted
        if frozenset((source, target)) in {frozenset(edge) for edge in truth}
    )
    orientation_accuracy = _safe_ratio(true_positive, unoriented_overlap)
    return {
        "n_predicted": len(predicted),
        "n_true": len(truth),
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "orientation_accuracy_among_unoriented_matches": orientation_accuracy,
    }


def score_static_ttd_graph_submission(
    truth: StaticTTDGraphTruth,
    observations: StaticTTDGraphObservations,
    submission: StaticTTDGraphSubmission,
) -> StaticTTDGraphScore:
    """Score a serialised truth-blind submission against the sealed virtual truth."""

    if truth.benchmark_id != observations.benchmark_id:
        raise ValueError("Truth and observations belong to different benchmark cases.")
    intervals = [
        submission.node_intervals.get(node)
        for node in observations.protocol.evaluation_nodes
    ]
    abstentions = [
        interval
        for interval in intervals
        if interval is None or interval.status == "ABSTAIN"
    ]
    estimates = [
        interval
        for interval in intervals
        if interval is not None and interval.status == "ESTIMATED"
    ]
    covered = 0
    widths: list[float] = []
    for node, interval in zip(observations.protocol.evaluation_nodes, intervals):
        if interval is None or interval.status != "ESTIMATED":
            continue
        truth_value = float(truth.young_water_fractions[node])
        covered += int(float(interval.lower) <= truth_value <= float(interval.upper))
        widths.append(float(interval.upper) - float(interval.lower))

    prediction_map = {
        _held_out_key(
            prediction.node_id, prediction.tracer, prediction.time_day
        ): float(prediction.value)
        for prediction in submission.held_out_predictions
    }
    errors: list[float] = []
    for target in observations.held_out_targets:
        key = _held_out_key(target.node_id, target.tracer, target.time_day)
        predicted = prediction_map.get(key)
        if predicted is not None:
            errors.append(float(predicted - truth.held_out_values[key]))
    metrics = {
        "protocol": {
            "primary_estimand": "node_young_water_fraction",
            "nominal_interval_coverage": observations.protocol.nominal_interval_coverage,
            "coverage_semantics": "conditional_on_non_abstention",
        },
        "interval_recovery": {
            "n_evaluation_nodes": len(observations.protocol.evaluation_nodes),
            "n_estimated": len(estimates),
            "n_abstained": len(abstentions),
            "abstention_rate": _safe_ratio(
                len(abstentions), len(observations.protocol.evaluation_nodes)
            ),
            "covered": covered,
            "conditional_coverage": _safe_ratio(covered, len(estimates)),
            "mean_interval_width": float(np.mean(widths)) if widths else None,
            "abstention_reasons": {
                node: (
                    submission.node_intervals.get(node).reason
                    if submission.node_intervals.get(node)
                    else "missing_node_submission"
                )
                for node in observations.protocol.evaluation_nodes
                if submission.node_intervals.get(node) is None
                or submission.node_intervals[node].status == "ABSTAIN"
            },
        },
        "held_out_prediction": {
            "n_targets": len(observations.held_out_targets),
            "n_predicted": len(errors),
            "prediction_rate": _safe_ratio(
                len(errors), len(observations.held_out_targets)
            ),
            "mae": float(np.mean(np.abs(errors))) if errors else None,
            "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
        },
        "topology": {
            "candidate_graph": _topology_metrics(
                submission.candidate_edges, truth.true_edges
            ),
            "selected_edges": _topology_metrics(
                submission.selected_edges, truth.true_edges
            ),
        },
        "provenance": {
            "truth_used_by_submission": False,
            "sealed_truth_sha256": _canonical_sha256(truth.to_dict()),
            "public_observations_sha256": _canonical_sha256(observations.to_dict()),
        },
    }
    return StaticTTDGraphScore(
        benchmark_id=truth.benchmark_id,
        method_id=submission.method_id,
        control_id=submission.control_id,
        metrics=metrics,
    )


__all__ = [
    "CandidateGraphControl",
    "CalibrationObservationRow",
    "EdgeKernelTruth",
    "HeldOutPrediction",
    "HeldOutTarget",
    "InputObservationRow",
    "SCHEMA_VERSION",
    "StaticTTDGraphCase",
    "StaticTTDGraphObservations",
    "StaticTTDGraphProtocol",
    "StaticTTDGraphScore",
    "StaticTTDGraphSubmission",
    "StaticTTDGraphTruth",
    "YoungWaterInterval",
    "default_static_ttd_graph_protocol",
    "generate_static_ttd_graph_case",
    "run_all_pre_registered_controls",
    "run_candidate_graph_estimator",
    "run_local_baseline",
    "score_static_ttd_graph_submission",
    "write_static_ttd_graph_case",
]
