"""Truth-sealed virtual benchmark for dynamic graph transit-time distributions.

This module deliberately contains a *forward* synthetic generator and a small,
truth-blind baseline.  It is not a field-data workflow and it does not modify
or condition the local feasible sets in :mod:`hydrosheaf.nuclear.ttd_graph`.

The generator creates a directed recharge--well network with time-varying
edge kernels, local recharge, a bypass pathway, noisy irregular observations,
and a future-time holdout.  The true edge kernels and latent signals are kept
in :class:`DynamicTTDTruth`; an inference routine receives only a separate
:class:`DynamicTTDObservations` instance.  A cryptographic commitment lets an
evaluator verify that those two objects belong to the same frozen case without
placing the truth arrays in the observation bundle.

The included ``recover_dynamic_ttd_baseline`` is intentionally modest: it
uses phase-specific lagged correlations and returns ``ABSTAIN`` whenever the
sampling cannot support enough seasonal-phase estimates.  It is a matched
controlled-synthetic baseline, not a claim that edge TTDs are field-identified.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np


GENERATOR_NAME = "hydrosheaf.ttd_graph_dynamic"
GENERATOR_VERSION = "1.0.0"
SOURCE_NODE = "R"
NETWORK_NODES = ("A", "B", "C")
TRUE_EDGES = (("R", "A"), ("A", "B"), ("B", "C"), ("A", "C"))
SCENARIOS = (
    "nominal",
    "sparse_sampling",
    "wrong_forcing",
    "unmodelled_local_recharge",
    "wrong_topology",
    "reversed_graph",
    "random_graph",
    "edge_removed_graph",
)


@dataclass(frozen=True)
class DynamicTTDGraphConfig:
    """Parameters for one deterministic controlled-synthetic case.

    ``n_steps`` is the number of reported observation-grid steps.  The
    generator uses an internal warm-up period so an observation at step zero
    is not biased by a zero-padded convolution boundary.
    """

    seed: int = 20_260_917
    scenario: str = "nominal"
    n_steps: int = 208
    step_days: float = 7.0
    season_period_steps: int = 52
    max_lag_steps: int = 18
    observation_probability: float = 0.82
    observation_noise_std: float = 0.035
    missing_block_steps: int = 5
    n_phase_bins: int = 4
    training_fraction: float = 0.70
    min_pairs_per_phase: int = 10
    min_identified_phases: int = 2
    correlation_gate: float = 0.22

    def __post_init__(self) -> None:
        if self.scenario not in SCENARIOS:
            raise ValueError(
                f"scenario must be one of {SCENARIOS}, got {self.scenario!r}"
            )
        if self.n_steps < 48:
            raise ValueError("n_steps must be at least 48 for a future holdout")
        if self.step_days <= 0.0:
            raise ValueError("step_days must be positive")
        if self.season_period_steps < self.n_phase_bins:
            raise ValueError("season_period_steps must be >= n_phase_bins")
        if self.max_lag_steps < 2:
            raise ValueError("max_lag_steps must be at least 2")
        if not 0.0 < self.observation_probability <= 1.0:
            raise ValueError("observation_probability must be in (0, 1]")
        if self.observation_noise_std < 0.0:
            raise ValueError("observation_noise_std must be non-negative")
        if self.missing_block_steps < 0:
            raise ValueError("missing_block_steps must be non-negative")
        if not 2 <= self.n_phase_bins <= 12:
            raise ValueError("n_phase_bins must be between 2 and 12")
        if not 0.45 <= self.training_fraction < 0.9:
            raise ValueError("training_fraction must be in [0.45, 0.9)")
        if self.min_pairs_per_phase < 3:
            raise ValueError("min_pairs_per_phase must be at least 3")
        if not 1 <= self.min_identified_phases <= self.n_phase_bins:
            raise ValueError("min_identified_phases is outside phase-bin range")
        if not -1.0 <= self.correlation_gate <= 1.0:
            raise ValueError("correlation_gate must be in [-1, 1]")


@dataclass(frozen=True)
class DynamicTTDTruth:
    """Sealed latent quantities for a virtual case.

    This object is only for generator-side storage and evaluator-side scoring.
    It must not be supplied to a recovery method.
    """

    case_id: str
    scenario: str
    time_steps: np.ndarray
    true_forcing: np.ndarray
    node_signals: Mapping[str, np.ndarray]
    edge_kernels: Mapping[str, np.ndarray]
    local_recharge_signals: Mapping[str, np.ndarray]
    local_recharge_fractions: Mapping[str, np.ndarray]
    direct_path_fraction: np.ndarray
    true_edges: tuple[tuple[str, str], ...]
    expected_identifiability: Mapping[str, bool]
    truth_commitment: str


@dataclass(frozen=True)
class DynamicTTDObservations:
    """Observable input to an inference method, with no true TTD arrays."""

    case_id: str
    scenario: str
    time_steps: np.ndarray
    forcing: np.ndarray
    node_observations: Mapping[str, np.ndarray]
    node_observation_masks: Mapping[str, np.ndarray]
    candidate_edges: tuple[tuple[str, str], ...]
    local_recharge_inputs: Mapping[str, np.ndarray]
    local_recharge_available: bool
    training_end_step: int
    # This one-way commitment is safe for an inference-visible object.  Its
    # name deliberately avoids the ``truth_`` prefix checked by the shared
    # truth-blindness guard.
    sealed_case_commitment: str
    observation_digest: str
    declared_stressors: tuple[str, ...]
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class DynamicTTDRecovery:
    """Result of the truth-blind seasonal lag baseline for one directed edge."""

    edge_id: str
    status: str
    reason: str
    mean_lag_steps: float | None
    phase_lag_steps: tuple[float | None, ...]
    phase_correlations: tuple[float | None, ...]
    phase_pair_counts: tuple[int, ...]
    training_end_step: int
    forecast_n: int
    forecast_rmse: float | None
    forecast_r2: float | None


@dataclass(frozen=True)
class DynamicTTDEvaluation:
    """Controlled-synthetic score calculated only after recovery is frozen."""

    case_id: str
    truth_commitment_verified: bool
    per_edge: Mapping[str, Mapping[str, Any]]
    summary: Mapping[str, Any]
    scope_note: str


@dataclass(frozen=True)
class DynamicTTDScenarioRun:
    """One generated stress case, its truth-blind baseline, and its score."""

    observations: DynamicTTDObservations
    manifest: Mapping[str, Any]
    recoveries: tuple[DynamicTTDRecovery, ...]
    evaluation: DynamicTTDEvaluation


@dataclass(frozen=True)
class _ScenarioControl:
    observation_probability_multiplier: float
    extra_missing_block_steps: int
    forcing_shift_steps: int
    forcing_amplitude_multiplier: float
    local_recharge_multiplier: float
    direct_path_multiplier: float
    local_recharge_available: bool
    candidate_edges: tuple[tuple[str, str], ...]
    stressors: tuple[str, ...]


def edge_id(edge: tuple[str, str]) -> str:
    """Return the canonical directed-edge identifier used by this benchmark."""

    return f"{edge[0]}->{edge[1]}"


def _scenario_control(scenario: str) -> _ScenarioControl:
    """Return declared stress controls without looking at generated truth."""

    controls = {
        "nominal": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.0,
            local_recharge_available=True,
            candidate_edges=TRUE_EDGES,
            stressors=("irregular_observations", "dynamic_kernels"),
        ),
        "sparse_sampling": _ScenarioControl(
            observation_probability_multiplier=0.30,
            extra_missing_block_steps=18,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.0,
            local_recharge_available=True,
            candidate_edges=TRUE_EDGES,
            stressors=(
                "sparse_irregular_observations",
                "long_missing_blocks",
                "dynamic_kernels",
            ),
        ),
        "wrong_forcing": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=5,
            forcing_amplitude_multiplier=0.78,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.0,
            local_recharge_available=True,
            candidate_edges=TRUE_EDGES,
            stressors=("mis-specified_recharge_forcing", "dynamic_kernels"),
        ),
        "unmodelled_local_recharge": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=2.65,
            direct_path_multiplier=1.0,
            local_recharge_available=False,
            candidate_edges=TRUE_EDGES,
            stressors=(
                "unobserved_time_varying_local_recharge",
                "dynamic_mixing",
            ),
        ),
        "wrong_topology": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.75,
            local_recharge_available=True,
            candidate_edges=(("R", "A"), ("A", "B"), ("B", "C")),
            stressors=("omitted_bypass_edge", "dynamic_mixing"),
        ),
        "reversed_graph": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.0,
            local_recharge_available=True,
            candidate_edges=tuple((target, source) for source, target in TRUE_EDGES),
            stressors=("reversed_candidate_edge_directions", "dynamic_kernels"),
        ),
        "random_graph": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.0,
            local_recharge_available=True,
            candidate_edges=(("R", "C"), ("C", "A"), ("C", "B")),
            stressors=("randomized_candidate_graph", "dynamic_kernels"),
        ),
        "edge_removed_graph": _ScenarioControl(
            observation_probability_multiplier=1.0,
            extra_missing_block_steps=0,
            forcing_shift_steps=0,
            forcing_amplitude_multiplier=1.0,
            local_recharge_multiplier=1.0,
            direct_path_multiplier=1.0,
            local_recharge_available=True,
            candidate_edges=tuple(edge for edge in TRUE_EDGES if edge != ("B", "C")),
            stressors=("edge_removed_candidate_graph", "dynamic_kernels"),
        ),
    }
    return controls[scenario]


def _readonly(values: Any, dtype: Any | None = float) -> np.ndarray:
    array = np.array(values, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _readonly_mapping(
    values: Mapping[str, Any], dtype: Any | None = float
) -> Mapping[str, np.ndarray]:
    return MappingProxyType({key: _readonly(value, dtype) for key, value in values.items()})


def _mapping_proxy(values: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(values))


def _config_digest(config: DynamicTTDGraphConfig) -> str:
    payload = json.dumps(
        asdict(config), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _update_array_digest(digest: Any, name: str, values: np.ndarray) -> None:
    array = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    digest.update(name.encode("utf-8"))
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.tobytes())


def _truth_digest(truth: DynamicTTDTruth) -> str:
    digest = sha256()
    digest.update(GENERATOR_NAME.encode("utf-8"))
    digest.update(GENERATOR_VERSION.encode("utf-8"))
    digest.update(truth.case_id.encode("utf-8"))
    digest.update(truth.scenario.encode("utf-8"))
    _update_array_digest(digest, "time_steps", truth.time_steps)
    _update_array_digest(digest, "true_forcing", truth.true_forcing)
    for name in sorted(truth.node_signals):
        _update_array_digest(digest, f"node:{name}", truth.node_signals[name])
    for name in sorted(truth.edge_kernels):
        _update_array_digest(digest, f"kernel:{name}", truth.edge_kernels[name])
    for name in sorted(truth.local_recharge_signals):
        _update_array_digest(
            digest, f"local_signal:{name}", truth.local_recharge_signals[name]
        )
        _update_array_digest(
            digest,
            f"local_fraction:{name}",
            truth.local_recharge_fractions[name],
        )
    _update_array_digest(digest, "direct_path_fraction", truth.direct_path_fraction)
    digest.update(json.dumps(truth.true_edges).encode("utf-8"))
    digest.update(
        json.dumps(dict(sorted(truth.expected_identifiability.items()))).encode("utf-8")
    )
    return digest.hexdigest()


def _observation_digest(observations: DynamicTTDObservations) -> str:
    digest = sha256()
    digest.update(GENERATOR_NAME.encode("utf-8"))
    digest.update(observations.case_id.encode("utf-8"))
    _update_array_digest(digest, "time_steps", observations.time_steps)
    _update_array_digest(digest, "forcing", observations.forcing)
    for name in sorted(observations.node_observations):
        _update_array_digest(
            digest, f"node_observation:{name}", observations.node_observations[name]
        )
        _update_array_digest(
            digest,
            f"node_observation_mask:{name}",
            observations.node_observation_masks[name].astype(float),
        )
    for name in sorted(observations.local_recharge_inputs):
        _update_array_digest(
            digest, f"local_input:{name}", observations.local_recharge_inputs[name]
        )
    digest.update(json.dumps(observations.candidate_edges).encode("utf-8"))
    digest.update(str(observations.training_end_step).encode("ascii"))
    digest.update(observations.sealed_case_commitment.encode("ascii"))
    return digest.hexdigest()


def _event_forcing(
    time_steps: np.ndarray, rng: np.random.Generator, period: int
) -> np.ndarray:
    """Create a seasonal, event-bearing conservative-tracer forcing signal."""

    phase = 2.0 * np.pi * time_steps / float(period)
    forcing = 0.95 * np.sin(phase) + 0.42 * np.cos(2.0 * phase - 0.35)
    # Short positive and negative events break the seasonal lag aliases that
    # make a periodic forcing alone weakly informative about a TTD shape.
    centres = rng.integers(3, len(time_steps) - 3, size=max(8, len(time_steps) // 20))
    for centre in centres:
        width = float(rng.uniform(0.7, 2.4))
        amplitude = float(rng.uniform(0.45, 1.05) * rng.choice((-1.0, 1.0)))
        forcing += amplitude * np.exp(
            -0.5 * ((np.arange(len(time_steps)) - centre) / width) ** 2
        )
    forcing += rng.normal(0.0, 0.045, size=len(time_steps))
    return forcing


def _local_recharge_signal(
    time_steps: np.ndarray, phase_shift: float, rng: np.random.Generator, period: int
) -> np.ndarray:
    phase = 2.0 * np.pi * time_steps / float(period) + phase_shift
    result = 0.35 + 0.62 * np.sin(phase) - 0.19 * np.cos(2.0 * phase)
    # Local inputs have independent event timing, so they cannot be reduced to
    # a fixed offset of the regional recharge tracer.
    centres = rng.integers(2, len(time_steps) - 2, size=max(5, len(time_steps) // 32))
    for centre in centres:
        result += float(rng.uniform(-0.45, 0.45)) * np.exp(
            -0.5 * ((np.arange(len(time_steps)) - centre) / 1.4) ** 2
        )
    return result + rng.normal(0.0, 0.035, size=len(time_steps))


def _dynamic_kernel(
    time_steps: np.ndarray,
    *,
    max_lag_steps: int,
    period: int,
    centre_lag: float,
    lag_amplitude: float,
    width: float,
    phase_shift: float,
) -> np.ndarray:
    """Make normalized causal ``h_e(tau, t)`` rows indexed by output time."""

    lags = np.arange(max_lag_steps + 1, dtype=float)
    phase = 2.0 * np.pi * time_steps / float(period) + phase_shift
    centres = centre_lag + lag_amplitude * np.sin(phase)
    widths = np.maximum(width + 0.22 * np.cos(phase - 0.4), 0.55)
    kernels = np.empty((len(time_steps), len(lags)), dtype=float)
    for row, (centre, spread) in enumerate(zip(centres, widths)):
        weights = np.exp(-0.5 * ((lags - centre) / spread) ** 2)
        # The edge retains a negligible zero-age tail but no acausal support;
        # this helps distinguish an instantaneous mixture from transport.
        weights[0] *= 0.03
        kernels[row] = weights / weights.sum()
    return kernels


def _apply_output_indexed_kernel(signal: np.ndarray, kernels: np.ndarray) -> np.ndarray:
    """Apply a causal time-varying kernel ``h(tau, t)`` to a signal.

    Each row of ``kernels`` is used at output time ``t`` and only accesses
    ``signal[t - tau]``.  The generator uses a warm-up region before retaining
    observations, so the leading zero padding cannot influence a reported
    virtual observation.
    """

    if signal.ndim != 1 or kernels.ndim != 2 or len(signal) != len(kernels):
        raise ValueError("signal and kernels must have compatible one-dimensional time axes")
    output = np.zeros(len(signal), dtype=float)
    for step in range(len(signal)):
        n_lags = min(step, kernels.shape[1] - 1) + 1
        output[step] = float(
            np.dot(kernels[step, :n_lags], signal[step - np.arange(n_lags)])
        )
    return output


def _sample_mask(
    rng: np.random.Generator,
    n_steps: int,
    probability: float,
    period: int,
    missing_block_steps: int,
    node_offset: int,
) -> np.ndarray:
    """Generate non-periodic observation times and one or more missing blocks."""

    phase = 2.0 * np.pi * (np.arange(n_steps) + node_offset) / float(period)
    probabilities = np.clip(probability * (0.80 + 0.20 * np.sin(phase)), 0.02, 1.0)
    mask = rng.random(n_steps) < probabilities
    # A random start and different offset per node make sampling irregular and
    # avoid an artificial common observation calendar.
    total_block = missing_block_steps
    if total_block:
        start_upper = max(1, n_steps - total_block)
        start = int(rng.integers(0, start_upper))
        mask[start : start + total_block] = False
    # Do not force a regular cadence just to meet a convenience target.  The
    # baseline has to abstain when this realised calendar is insufficient.
    return mask


def _phase_index(
    time_steps: np.ndarray, period: int, n_phase_bins: int
) -> np.ndarray:
    within_season = np.mod(time_steps, period)
    phase = np.floor(n_phase_bins * within_season / float(period)).astype(int)
    return np.clip(phase, 0, n_phase_bins - 1)


def _phase_kernel_means(
    kernels: np.ndarray, time_steps: np.ndarray, period: int, n_phase_bins: int
) -> tuple[float, ...]:
    lags = np.arange(kernels.shape[1], dtype=float)
    means = kernels @ lags
    phases = _phase_index(time_steps, period, n_phase_bins)
    return tuple(float(np.mean(means[phases == phase])) for phase in range(n_phase_bins))


def _expected_identifiability(scenario: str) -> Mapping[str, bool]:
    """Pre-declare the target-level information regime of a virtual case.

    These labels are a controlled-synthetic scoring convention, not an oracle
    claim that a real aquifer is identifiable.  They are intentionally
    conservative: the mixed downstream targets are labelled unresolved unless
    a future inference method explicitly models the needed mixing inputs.
    """

    expected = {edge_id(edge): False for edge in TRUE_EDGES}
    if scenario == "nominal":
        expected["R->A"] = True
        expected["A->B"] = True
    elif scenario == "wrong_forcing":
        expected["A->B"] = True
    elif scenario == "wrong_topology":
        expected["R->A"] = True
        expected["A->B"] = True
    return _mapping_proxy(expected)


def generate_dynamic_ttd_case(
    config: DynamicTTDGraphConfig | None = None,
) -> tuple[DynamicTTDTruth, DynamicTTDObservations, Mapping[str, Any]]:
    """Generate one deterministic, truth-separated dynamic graph-TTD case.

    Returns
    -------
    truth, observations, manifest
        ``truth`` is retained by a benchmark evaluator; a recovery algorithm
        should receive only ``observations``.  The manifest contains hashes and
        declared protocol details, never the latent kernel or node-signal data.
    """

    config = config or DynamicTTDGraphConfig()
    control = _scenario_control(config.scenario)
    rng = np.random.default_rng(config.seed)
    config_hash = _config_digest(config)
    case_id = f"ttd-graph-dynamic-{config.scenario}-{config.seed}-{config_hash[:12]}"

    # Run two maximum-lag windows before the first reported point.  This is
    # generator-only state and never appears in the observation bundle.
    warm_up = 2 * config.max_lag_steps
    full_steps = np.arange(-warm_up, config.n_steps, dtype=float)
    retained = slice(warm_up, None)
    full_forcing = _event_forcing(full_steps, rng, config.season_period_steps)

    kernel_ra = _dynamic_kernel(
        full_steps,
        max_lag_steps=config.max_lag_steps,
        period=config.season_period_steps,
        centre_lag=2.7,
        lag_amplitude=0.75,
        width=1.0,
        phase_shift=0.2,
    )
    kernel_ab = _dynamic_kernel(
        full_steps,
        max_lag_steps=config.max_lag_steps,
        period=config.season_period_steps,
        centre_lag=5.4,
        lag_amplitude=1.55,
        width=1.35,
        phase_shift=1.1,
    )
    kernel_bc = _dynamic_kernel(
        full_steps,
        max_lag_steps=config.max_lag_steps,
        period=config.season_period_steps,
        centre_lag=8.2,
        lag_amplitude=2.00,
        width=1.70,
        phase_shift=-0.7,
    )
    kernel_ac = _dynamic_kernel(
        full_steps,
        max_lag_steps=config.max_lag_steps,
        period=config.season_period_steps,
        centre_lag=3.4,
        lag_amplitude=0.90,
        width=1.15,
        phase_shift=0.55,
    )

    local_a = _local_recharge_signal(
        full_steps, 1.6, rng, config.season_period_steps
    )
    local_b = _local_recharge_signal(
        full_steps, -1.1, rng, config.season_period_steps
    )
    local_c = _local_recharge_signal(
        full_steps, 2.35, rng, config.season_period_steps
    )
    phase = 2.0 * np.pi * full_steps / float(config.season_period_steps)
    local_a_fraction = np.clip(
        control.local_recharge_multiplier * (0.06 + 0.035 * np.sin(phase + 0.2)),
        0.01,
        0.72,
    )
    local_b_fraction = np.clip(
        control.local_recharge_multiplier * (0.10 + 0.065 * np.sin(phase - 0.5)),
        0.01,
        0.72,
    )
    local_c_fraction = np.clip(
        control.local_recharge_multiplier * (0.08 + 0.055 * np.sin(phase + 0.9)),
        0.01,
        0.72,
    )
    direct_path_fraction = np.clip(
        control.direct_path_multiplier * (0.24 + 0.12 * np.sin(phase - 0.3)),
        0.04,
        0.75,
    )

    # The forward model is a causal, output-time-indexed dynamic operator.
    # Local fractions and bypass mixing vary at the receiving time step.
    regional_a = _apply_output_indexed_kernel(full_forcing, kernel_ra)
    signal_a = (1.0 - local_a_fraction) * regional_a + local_a_fraction * local_a
    advected_b = _apply_output_indexed_kernel(signal_a, kernel_ab)
    signal_b = (1.0 - local_b_fraction) * advected_b + local_b_fraction * local_b
    via_b = _apply_output_indexed_kernel(signal_b, kernel_bc)
    bypass_a = _apply_output_indexed_kernel(signal_a, kernel_ac)
    advected_c = (1.0 - direct_path_fraction) * via_b + direct_path_fraction * bypass_a
    signal_c = (1.0 - local_c_fraction) * advected_c + local_c_fraction * local_c

    time_steps = np.arange(config.n_steps, dtype=float)
    truth = DynamicTTDTruth(
        case_id=case_id,
        scenario=config.scenario,
        time_steps=_readonly(time_steps),
        true_forcing=_readonly(full_forcing[retained]),
        node_signals=_readonly_mapping(
            {"A": signal_a[retained], "B": signal_b[retained], "C": signal_c[retained]}
        ),
        edge_kernels=_readonly_mapping(
            {
                "R->A": kernel_ra[retained],
                "A->B": kernel_ab[retained],
                "B->C": kernel_bc[retained],
                "A->C": kernel_ac[retained],
            }
        ),
        local_recharge_signals=_readonly_mapping(
            {"A": local_a[retained], "B": local_b[retained], "C": local_c[retained]}
        ),
        local_recharge_fractions=_readonly_mapping(
            {
                "A": local_a_fraction[retained],
                "B": local_b_fraction[retained],
                "C": local_c_fraction[retained],
            }
        ),
        direct_path_fraction=_readonly(direct_path_fraction[retained]),
        true_edges=TRUE_EDGES,
        expected_identifiability=_expected_identifiability(config.scenario),
        truth_commitment="",
    )
    truth = replace(truth, truth_commitment=_truth_digest(truth))

    observed_forcing = truth.true_forcing.copy()
    if control.forcing_shift_steps:
        observed_forcing = np.roll(observed_forcing, control.forcing_shift_steps)
        observed_forcing = control.forcing_amplitude_multiplier * observed_forcing
        observed_forcing += 0.22 * np.sin(
            2.0 * np.pi * (time_steps + 3.0) / float(config.season_period_steps)
        )
    observed_forcing += rng.normal(
        0.0, config.observation_noise_std * 0.45, size=config.n_steps
    )

    observation_probability = (
        config.observation_probability * control.observation_probability_multiplier
    )
    missing_block = config.missing_block_steps + control.extra_missing_block_steps
    node_masks: dict[str, np.ndarray] = {}
    node_observations: dict[str, np.ndarray] = {}
    for offset, node in enumerate(NETWORK_NODES):
        mask = _sample_mask(
            rng,
            config.n_steps,
            observation_probability,
            config.season_period_steps,
            missing_block,
            node_offset=7 * (offset + 1),
        )
        noisy = truth.node_signals[node] + rng.normal(
            0.0, config.observation_noise_std, size=config.n_steps
        )
        node_masks[node] = mask
        node_observations[node] = np.where(mask, noisy, np.nan)

    local_inputs: dict[str, np.ndarray] = {}
    if control.local_recharge_available:
        for node in NETWORK_NODES:
            # These inputs are intentionally imperfect measurements.  The
            # baseline below does not use them, leaving a clear hook for a
            # future source/mixing-aware method to test against the same case.
            local_inputs[node] = truth.local_recharge_signals[node] + rng.normal(
                0.0, config.observation_noise_std * 1.5, size=config.n_steps
            )

    training_end_step = int(np.floor(config.training_fraction * config.n_steps))
    observations = DynamicTTDObservations(
        case_id=case_id,
        scenario=config.scenario,
        time_steps=_readonly(time_steps),
        forcing=_readonly(observed_forcing),
        node_observations=_readonly_mapping(node_observations),
        node_observation_masks=_readonly_mapping(node_masks, dtype=bool),
        candidate_edges=control.candidate_edges,
        local_recharge_inputs=_readonly_mapping(local_inputs),
        local_recharge_available=control.local_recharge_available,
        training_end_step=training_end_step,
        sealed_case_commitment=truth.truth_commitment,
        observation_digest="",
        declared_stressors=control.stressors,
        metadata=_mapping_proxy(
            {
                "generator": GENERATOR_NAME,
                "generator_version": GENERATOR_VERSION,
                "validation_scope": "controlled_synthetic_only",
                "dynamic_operator": "causal output-time-indexed h(tau, t)",
                "max_lag_steps": config.max_lag_steps,
                "n_phase_bins": config.n_phase_bins,
                "season_period_steps": config.season_period_steps,
                "min_pairs_per_phase": config.min_pairs_per_phase,
                "min_identified_phases": config.min_identified_phases,
                "correlation_gate": config.correlation_gate,
                "local_recharge_available": control.local_recharge_available,
                "heldout_target_rule": "destination observations at t >= training_end_step are not used for fitting",
                "edge_truth_visible_to_inference": False,
            }
        ),
    )
    observations = replace(observations, observation_digest=_observation_digest(observations))
    manifest: Mapping[str, Any] = {
        "schema_version": "1.0",
        "generator": GENERATOR_NAME,
        "generator_version": GENERATOR_VERSION,
        "case_id": case_id,
        "seed": config.seed,
        "scenario": config.scenario,
        "configuration_digest": config_hash,
        "n_steps": config.n_steps,
        "step_days": config.step_days,
        "season_period_steps": config.season_period_steps,
        "max_lag_steps": config.max_lag_steps,
        "training_end_step": training_end_step,
        "candidate_edges": [edge_id(edge) for edge in observations.candidate_edges],
        "declared_stressors": list(control.stressors),
        "sealed_case_commitment": truth.truth_commitment,
        "observation_digest": observations.observation_digest,
        "sealed_payload_in_manifest": False,
        "validation_scope": "controlled_synthetic_only",
        "claim_boundary": (
            "This case evaluates recovery and abstention against sealed virtual truth; "
            "it is not field validation or evidence of field superiority."
        ),
    }
    return truth, observations, manifest


def _source_series_for_edge(
    observations: DynamicTTDObservations, edge: tuple[str, str]
) -> np.ndarray:
    source, _ = edge
    if source == SOURCE_NODE:
        return observations.forcing
    if source not in observations.node_observations:
        raise ValueError(f"candidate source node {source!r} has no observed series")
    return observations.node_observations[source]


def _best_phase_lag(
    source: np.ndarray,
    target_training: np.ndarray,
    phase_at_time: np.ndarray,
    phase: int,
    max_lag_steps: int,
    min_pairs: int,
) -> tuple[float | None, float | None, int, float | None, float | None]:
    """Fit one phase's lag and affine response using training data only."""

    best: tuple[float, float, int, float, float] | None = None
    train_end = len(target_training)
    for lag in range(1, max_lag_steps + 1):
        target_idx = np.arange(lag, train_end)
        use = phase_at_time[target_idx] == phase
        if not np.any(use):
            continue
        target_idx = target_idx[use]
        source_values = source[target_idx - lag]
        target_values = target_training[target_idx]
        finite = np.isfinite(source_values) & np.isfinite(target_values)
        source_values = source_values[finite]
        target_values = target_values[finite]
        n_pairs = len(source_values)
        if n_pairs < min_pairs:
            continue
        source_std = float(np.std(source_values))
        target_std = float(np.std(target_values))
        if source_std <= 1.0e-10 or target_std <= 1.0e-10:
            continue
        correlation = float(np.corrcoef(source_values, target_values)[0, 1])
        if not np.isfinite(correlation):
            continue
        design = np.column_stack((source_values, np.ones(n_pairs)))
        slope, intercept = np.linalg.lstsq(design, target_values, rcond=None)[0]
        candidate = (correlation, float(lag), n_pairs, float(slope), float(intercept))
        if best is None or candidate[0] > best[0]:
            best = candidate
    if best is None:
        return None, None, 0, None, None
    correlation, lag, n_pairs, slope, intercept = best
    return lag, correlation, n_pairs, slope, intercept


def _forecast_metrics(
    source: np.ndarray,
    target: np.ndarray,
    phase_at_time: np.ndarray,
    training_end_step: int,
    phase_lags: Sequence[float | None],
    phase_models: Sequence[tuple[float | None, float | None]],
) -> tuple[int, float | None, float | None]:
    predictions: list[float] = []
    observed: list[float] = []
    for step in range(training_end_step, len(target)):
        phase = int(phase_at_time[step])
        lag = phase_lags[phase]
        slope, intercept = phase_models[phase]
        if lag is None or slope is None or intercept is None:
            continue
        source_step = step - int(round(lag))
        if source_step < 0 or not np.isfinite(source[source_step]):
            continue
        if not np.isfinite(target[step]):
            continue
        predictions.append(float(slope * source[source_step] + intercept))
        observed.append(float(target[step]))
    n_forecast = len(observed)
    if n_forecast == 0:
        return 0, None, None
    errors = np.asarray(predictions) - np.asarray(observed)
    rmse = float(np.sqrt(np.mean(errors**2)))
    if n_forecast < 2:
        return n_forecast, rmse, None
    centered = np.asarray(observed) - float(np.mean(observed))
    total = float(np.sum(centered**2))
    r2 = None if total <= 1.0e-12 else float(1.0 - np.sum(errors**2) / total)
    return n_forecast, rmse, r2


def recover_dynamic_ttd_baseline(
    observations: DynamicTTDObservations,
    *,
    max_lag_steps: int | None = None,
    n_phase_bins: int | None = None,
    season_period_steps: int | None = None,
    min_pairs_per_phase: int | None = None,
    min_identified_phases: int | None = None,
    correlation_gate: float | None = None,
) -> tuple[DynamicTTDRecovery, ...]:
    """Run a truth-blind phase-specific lag baseline on candidate edges.

    Destination observations from the future holdout are explicitly replaced
    with ``NaN`` during fitting.  They are used only after the baseline has
    selected its lag/affine model, for reported future-time prediction metrics.
    """

    metadata = observations.metadata
    max_lag_steps = int(
        max_lag_steps if max_lag_steps is not None else metadata.get("max_lag_steps", 18)
    )
    n_phase_bins = int(
        n_phase_bins if n_phase_bins is not None else metadata.get("n_phase_bins", 4)
    )
    season_period_steps = int(
        season_period_steps
        if season_period_steps is not None
        else metadata.get("season_period_steps", 52)
    )
    min_pairs_per_phase = int(
        min_pairs_per_phase
        if min_pairs_per_phase is not None
        else metadata.get("min_pairs_per_phase", 10)
    )
    min_identified_phases = int(
        min_identified_phases
        if min_identified_phases is not None
        else metadata.get("min_identified_phases", 2)
    )
    correlation_gate = float(
        correlation_gate
        if correlation_gate is not None
        else metadata.get("correlation_gate", 0.22)
    )
    if max_lag_steps < 1:
        raise ValueError("max_lag_steps must be positive")
    if not 1 <= min_identified_phases <= n_phase_bins:
        raise ValueError("min_identified_phases is outside phase-bin range")

    phase_at_time = _phase_index(
        observations.time_steps, season_period_steps, n_phase_bins
    )
    recoveries: list[DynamicTTDRecovery] = []
    for edge in observations.candidate_edges:
        source, target_node = edge
        name = edge_id(edge)
        if target_node not in observations.node_observations:
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason=f"candidate target node {target_node!r} has no observed series",
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(None for _ in range(n_phase_bins)),
                    phase_correlations=tuple(None for _ in range(n_phase_bins)),
                    phase_pair_counts=tuple(0 for _ in range(n_phase_bins)),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue
        source_series = _source_series_for_edge(observations, edge)
        target_series = observations.node_observations[target_node]
        target_training = np.array(target_series, copy=True)
        target_training[observations.training_end_step :] = np.nan

        phase_lags: list[float | None] = []
        phase_correlations: list[float | None] = []
        phase_pairs: list[int] = []
        phase_models: list[tuple[float | None, float | None]] = []
        for phase in range(n_phase_bins):
            lag, corr, count, slope, intercept = _best_phase_lag(
                source_series,
                target_training,
                phase_at_time,
                phase,
                max_lag_steps,
                min_pairs_per_phase,
            )
            if corr is None or corr < correlation_gate:
                phase_lags.append(None)
                phase_correlations.append(corr)
                phase_pairs.append(count)
                phase_models.append((None, None))
            else:
                phase_lags.append(lag)
                phase_correlations.append(corr)
                phase_pairs.append(count)
                phase_models.append((slope, intercept))

        n_identified = sum(lag is not None for lag in phase_lags)
        if n_identified < min_identified_phases:
            recoveries.append(
                DynamicTTDRecovery(
                    edge_id=name,
                    status="ABSTAIN",
                    reason=(
                        f"only {n_identified}/{n_phase_bins} seasonal phases passed "
                        f"the {min_pairs_per_phase}-pair and correlation gates"
                    ),
                    mean_lag_steps=None,
                    phase_lag_steps=tuple(phase_lags),
                    phase_correlations=tuple(phase_correlations),
                    phase_pair_counts=tuple(phase_pairs),
                    training_end_step=observations.training_end_step,
                    forecast_n=0,
                    forecast_rmse=None,
                    forecast_r2=None,
                )
            )
            continue
        forecast_n, forecast_rmse, forecast_r2 = _forecast_metrics(
            source_series,
            target_series,
            phase_at_time,
            observations.training_end_step,
            phase_lags,
            phase_models,
        )
        retained_lags = [lag for lag in phase_lags if lag is not None]
        recoveries.append(
            DynamicTTDRecovery(
                edge_id=name,
                status="RECOVERED",
                reason=(
                    "phase-specific lag estimates passed declared sampling and "
                    "correlation gates; controlled-synthetic scoring remains required"
                ),
                mean_lag_steps=float(np.mean(retained_lags)),
                phase_lag_steps=tuple(phase_lags),
                phase_correlations=tuple(phase_correlations),
                phase_pair_counts=tuple(phase_pairs),
                training_end_step=observations.training_end_step,
                forecast_n=forecast_n,
                forecast_rmse=forecast_rmse,
                forecast_r2=forecast_r2,
            )
        )
    return tuple(recoveries)


def evaluate_dynamic_ttd_recovery(
    truth: DynamicTTDTruth,
    observations: DynamicTTDObservations,
    recoveries: Sequence[DynamicTTDRecovery],
    *,
    season_period_steps: int = 52,
    n_phase_bins: int = 4,
    lag_tolerance_steps: float = 2.5,
) -> DynamicTTDEvaluation:
    """Score frozen recoveries against sealed truth after inference completes.

    The evaluator deliberately separates an expected-identifiability label
    from numerical error.  A point estimate on a declared unresolved case is
    counted as an unsupported estimate, even when it happens to be near the
    hidden answer; an abstention is counted separately from a missed recovery.
    """

    if truth.case_id != observations.case_id:
        raise ValueError("truth and observations belong to different benchmark cases")
    commitment_verified = (
        truth.truth_commitment == observations.sealed_case_commitment
        and truth.truth_commitment == _truth_digest(replace(truth, truth_commitment=""))
    )
    if not commitment_verified:
        raise ValueError("truth commitment does not verify against observations")
    if lag_tolerance_steps <= 0.0:
        raise ValueError("lag_tolerance_steps must be positive")

    truth_phase_means = {
        name: _phase_kernel_means(
            kernel,
            truth.time_steps,
            season_period_steps,
            n_phase_bins,
        )
        for name, kernel in truth.edge_kernels.items()
    }
    per_edge: dict[str, Mapping[str, Any]] = {}
    justified = false_abstentions = correct_abstentions = unsupported = misses = 0
    forecast_r2_values: list[float] = []
    forecast_n_total = 0
    for recovery in recoveries:
        expected = bool(truth.expected_identifiability.get(recovery.edge_id, False))
        true_phase = truth_phase_means.get(recovery.edge_id)
        absolute_errors: list[float] = []
        if true_phase is not None:
            for estimated, actual in zip(recovery.phase_lag_steps, true_phase):
                if estimated is not None:
                    absolute_errors.append(abs(float(estimated) - float(actual)))
        mean_abs_error = (
            float(np.mean(absolute_errors)) if absolute_errors else None
        )
        if expected and recovery.status == "ABSTAIN":
            classification = "FALSE_ABSTENTION"
            false_abstentions += 1
        elif not expected and recovery.status == "ABSTAIN":
            classification = "CORRECT_ABSTENTION"
            correct_abstentions += 1
        elif not expected:
            classification = "UNSUPPORTED_POINT_ESTIMATE"
            unsupported += 1
        elif mean_abs_error is not None and mean_abs_error <= lag_tolerance_steps:
            classification = "IDENTIFIABLE_RECOVERY_WITHIN_TOLERANCE"
            justified += 1
        else:
            classification = "IDENTIFIABLE_RECOVERY_MISS"
            misses += 1
        if recovery.forecast_r2 is not None:
            forecast_r2_values.append(recovery.forecast_r2)
        forecast_n_total += recovery.forecast_n
        per_edge[recovery.edge_id] = _mapping_proxy(
            {
                "expected_identifiable_in_virtual_protocol": expected,
                "status": recovery.status,
                "classification": classification,
                "mean_absolute_phase_lag_error_steps": mean_abs_error,
                "forecast_n": recovery.forecast_n,
                "forecast_rmse": recovery.forecast_rmse,
                "forecast_r2": recovery.forecast_r2,
            }
        )
    summary = _mapping_proxy(
        {
            "n_recoveries": len(recoveries),
            "justified_recoveries": justified,
            "identifiable_recovery_misses": misses,
            "false_abstentions": false_abstentions,
            "correct_abstentions": correct_abstentions,
            "unsupported_point_estimates": unsupported,
            "heldout_forecast_n": forecast_n_total,
            "mean_heldout_forecast_r2": (
                float(np.mean(forecast_r2_values)) if forecast_r2_values else None
            ),
        }
    )
    return DynamicTTDEvaluation(
        case_id=truth.case_id,
        truth_commitment_verified=True,
        per_edge=_mapping_proxy(per_edge),
        summary=summary,
        scope_note=(
            "Scores are controlled-synthetic recovery and abstention evidence only; "
            "they do not validate a field groundwater network or establish field superiority."
        ),
    )


def run_dynamic_ttd_stress_suite(
    config: DynamicTTDGraphConfig | None = None,
) -> Mapping[str, DynamicTTDScenarioRun]:
    """Generate and score all declared dynamic-TTD stress scenarios.

    Each scenario has a deterministic seed derived from the supplied base seed
    and runs independently.  This prevents a preceding scenario from changing
    the synthetic draw of a later one.
    """

    base = config or DynamicTTDGraphConfig()
    suite: dict[str, DynamicTTDScenarioRun] = {}
    for index, scenario in enumerate(SCENARIOS):
        scenario_config = replace(base, scenario=scenario, seed=base.seed + index)
        truth, observations, manifest = generate_dynamic_ttd_case(scenario_config)
        recoveries = recover_dynamic_ttd_baseline(
            observations,
            max_lag_steps=scenario_config.max_lag_steps,
            n_phase_bins=scenario_config.n_phase_bins,
            season_period_steps=scenario_config.season_period_steps,
            min_pairs_per_phase=scenario_config.min_pairs_per_phase,
            min_identified_phases=scenario_config.min_identified_phases,
            correlation_gate=scenario_config.correlation_gate,
        )
        evaluation = evaluate_dynamic_ttd_recovery(
            truth,
            observations,
            recoveries,
            season_period_steps=scenario_config.season_period_steps,
            n_phase_bins=scenario_config.n_phase_bins,
        )
        suite[scenario] = DynamicTTDScenarioRun(
            observations=observations,
            manifest=manifest,
            recoveries=recoveries,
            evaluation=evaluation,
        )
    return _mapping_proxy(suite)


__all__ = [
    "DynamicTTDGraphConfig",
    "DynamicTTDTruth",
    "DynamicTTDObservations",
    "DynamicTTDRecovery",
    "DynamicTTDEvaluation",
    "DynamicTTDScenarioRun",
    "GENERATOR_NAME",
    "GENERATOR_VERSION",
    "SCENARIOS",
    "edge_id",
    "generate_dynamic_ttd_case",
    "recover_dynamic_ttd_baseline",
    "evaluate_dynamic_ttd_recovery",
    "run_dynamic_ttd_stress_suite",
]
