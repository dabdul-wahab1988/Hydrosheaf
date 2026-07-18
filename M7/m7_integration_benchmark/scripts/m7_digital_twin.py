"""Operational synthetic-aquifer twin used by the optional M7 extension.

This module is deliberately separate from ``run_m7_integration_benchmark.py``.
The locked M7 integration benchmark asks whether age, topology, and chemistry
jointly reject false edges when joint truth is known.  This extension asks a
different question: can a graph-constrained model be updated sequentially from
sparse observations and make calibrated, genuinely out-of-sample forecasts?

The implementation is a controlled *synthetic* operational twin, not a field
digital twin.  It has:

* a transient hidden truth with graph transport, heterogeneous coefficients,
  nonlinear recharge response, an unmodelled preferential connection, and
  process noise;
* a deliberately simpler linear forecast model (anti-inverse-crime);
* sparse asynchronous head, chloride, and silica monitoring;
* an ensemble Kalman filter (EnKF) for sequential state updating;
* prequential forecasts issued before future observations are assimilated; and
* open-loop, oracle-persistence, seasonal, wrong-topology, shuffled-observation,
  and sensor-dropout controls.

Only NumPy and pandas are required.  The code is deterministic for a fixed seed.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


VARIABLES = ("head_m", "chloride_mmol_l", "silica_mmol_l")
N_VARIABLES = len(VARIABLES)


@dataclass(frozen=True)
class TwinConfig:
    """Locked numerical protocol for the operational-twin benchmark."""

    n_nodes: int = 24
    n_steps: int = 72
    spinup_steps: int = 6
    forecast_start: int = 42
    ensemble_size: int = 80
    n_replicates: int = 24
    horizons: tuple[int, ...] = (1, 3, 6)
    seed: int = 20260718
    inflation: float = 1.025
    dropout_fraction: float = 0.50
    interval: float = 0.90

    def validate(self) -> None:
        if self.n_nodes != 24:
            raise ValueError("The locked M7 topology requires n_nodes=24.")
        if self.n_steps <= max(self.horizons) + self.forecast_start:
            raise ValueError("n_steps must leave observations after every horizon.")
        if self.spinup_steps < 1 or self.forecast_start <= self.spinup_steps:
            raise ValueError("forecast_start must follow a non-empty spin-up period.")
        if self.ensemble_size < 12:
            raise ValueError("ensemble_size must be at least 12.")
        if self.n_replicates < 1:
            raise ValueError("n_replicates must be positive.")
        if not (0.0 <= self.dropout_fraction < 1.0):
            raise ValueError("dropout_fraction must lie in [0, 1).")
        if not (0.0 < self.interval < 1.0):
            raise ValueError("interval must lie in (0, 1).")


@dataclass(frozen=True)
class DynamicModel:
    """Parameters of either the hidden truth or an operational twin."""

    name: str
    upstream_weights: np.ndarray
    retention: np.ndarray
    coupling: np.ndarray
    input_coefficients: np.ndarray
    process_sd: np.ndarray
    node_sensitivity: np.ndarray
    nonlinear_recharge: float = 0.0
    hidden_pulse: bool = False


@dataclass(frozen=True)
class ObservationBatch:
    """Observations available at one operational update cycle."""

    time: int
    flat_indices: np.ndarray
    values: np.ndarray
    standard_deviations: np.ndarray


@dataclass(frozen=True)
class ExperimentData:
    """Known truth, forcings, and sparse observations for one replicate."""

    truth_normalized: np.ndarray
    truth_physical: np.ndarray
    realized_forcing: np.ndarray
    forecast_forcing: np.ndarray
    observations: tuple[ObservationBatch, ...]
    monitor_nodes: tuple[int, ...]


@dataclass(frozen=True)
class FilterRun:
    """Sequential ensemble histories and update diagnostics."""

    method: str
    posterior_ensembles: tuple[np.ndarray, ...]
    diagnostics: pd.DataFrame


def locked_edges() -> tuple[tuple[int, int], ...]:
    """True-edge topology from the locked M7 integration benchmark."""

    return (
        (0, 1),
        (1, 2),
        (2, 6),
        (3, 5),
        (4, 5),
        (5, 6),
        (6, 8),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 12),
        (11, 13),
        (12, 13),
        (13, 15),
        (14, 15),
        (15, 16),
        (16, 17),
        (17, 19),
        (18, 19),
        (19, 20),
        (20, 22),
        (21, 23),
        (22, 23),
    )


def _upstream_matrix(
    n_nodes: int,
    edges: Sequence[tuple[int, int]],
    edge_weights: Mapping[tuple[int, int], float] | None = None,
) -> np.ndarray:
    """Return a row-normalized upstream mixing matrix."""

    matrix = np.zeros((n_nodes, n_nodes), dtype=float)
    weights = edge_weights or {}
    for upstream, downstream in edges:
        matrix[downstream, upstream] += float(weights.get((upstream, downstream), 1.0))
    for node in range(n_nodes):
        total = float(matrix[node].sum())
        if total == 0.0:
            matrix[node, node] = 1.0
        else:
            matrix[node] /= total
    return matrix


def physical_baseline(n_nodes: int) -> np.ndarray:
    """Node-specific physical baseline for the three state variables."""

    progress = np.linspace(0.0, 1.0, n_nodes)
    baseline = np.empty((n_nodes, N_VARIABLES), dtype=float)
    baseline[:, 0] = 78.0 - 28.0 * progress
    baseline[:, 1] = 0.20 + 0.42 * progress
    baseline[:, 2] = 0.14 + 0.38 * progress
    return baseline


def physical_scale(n_nodes: int) -> np.ndarray:
    """Scale mapping dimensionless filter states into physical units."""

    return np.tile(np.array([3.0, 0.10, 0.08], dtype=float), (n_nodes, 1))


def to_physical(states: np.ndarray) -> np.ndarray:
    """Convert normalized state arrays (..., node, variable) to physical units."""

    return physical_baseline(states.shape[-2]) + physical_scale(states.shape[-2]) * states


def build_models(config: TwinConfig) -> tuple[DynamicModel, DynamicModel, DynamicModel]:
    """Construct hidden truth, nominal twin, and wrong-topology control models."""

    config.validate()
    edges = locked_edges()
    nominal_weights = _upstream_matrix(config.n_nodes, edges)

    # Truth includes an unmodelled preferential connection from the second branch
    # into the lower main stem.  Its weights differ from the operational model.
    truth_edges = (*edges, (10, 16))
    truth_edge_weights = {(10, 16): 0.30, (15, 16): 0.70}
    truth_weights = _upstream_matrix(config.n_nodes, truth_edges, truth_edge_weights)

    # The wrong-topology control preserves a directed, branching, down-gradient
    # network but misroutes both branch merges and the lower stem. It is a
    # structural control rather than a one-edge perturbation.
    wrong_edges = (
        (0, 1),
        (1, 2),
        (2, 5),
        (3, 5),
        (4, 5),
        (5, 6),
        (6, 8),
        (7, 9),
        (8, 9),
        (9, 10),
        (10, 12),
        (11, 13),
        (12, 13),
        (13, 15),
        (14, 16),
        (15, 16),
        (16, 17),
        (17, 19),
        (18, 20),
        (19, 20),
        (20, 22),
        (21, 22),
        (22, 23),
    )
    wrong_weights = _upstream_matrix(config.n_nodes, wrong_edges)

    progress = np.linspace(0.0, 1.0, config.n_nodes)[:, None]
    truth_sensitivity = np.hstack(
        (
            0.80 + 0.35 * (1.0 - progress),
            0.75 + 0.45 * progress,
            0.70 + 0.55 * progress,
        )
    )
    twin_sensitivity = np.ones((config.n_nodes, N_VARIABLES), dtype=float)

    truth = DynamicModel(
        name="hidden_truth",
        upstream_weights=truth_weights,
        retention=np.array([0.77, 0.84, 0.80]),
        coupling=np.array([0.14, 0.10, 0.13]),
        input_coefficients=np.array(
            [
                [0.24, -0.18, 0.02],  # head: recharge, pumping, temperature
                [-0.08, 0.035, 0.055],  # chloride: dilution, concentration, evap
                [0.055, 0.018, 0.045],  # silica: weathering/reaction response
            ]
        ),
        process_sd=np.array([0.045, 0.035, 0.035]),
        node_sensitivity=truth_sensitivity,
        nonlinear_recharge=0.030,
        hidden_pulse=True,
    )
    nominal = DynamicModel(
        name="updated_twin",
        upstream_weights=nominal_weights,
        retention=np.array([0.79, 0.86, 0.82]),
        coupling=np.array([0.11, 0.08, 0.10]),
        input_coefficients=np.array(
            [
                [0.21, -0.15, 0.015],
                [-0.065, 0.025, 0.045],
                [0.045, 0.012, 0.035],
            ]
        ),
        process_sd=np.array([0.060, 0.045, 0.045]),
        node_sensitivity=twin_sensitivity,
    )
    wrong = replace(nominal, name="wrong_topology_updated", upstream_weights=wrong_weights)
    return truth, nominal, wrong


def generate_forcing(config: TwinConfig, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Generate realized and operationally forecast forcing.

    Columns are standardized recharge, pumping, and temperature.  The operational
    forecast knows the seasonal climatology and declared pumping schedule but not
    future stochastic recharge anomalies.
    """

    time = np.arange(config.n_steps, dtype=float)
    seasonal_recharge = np.sin(2.0 * np.pi * (time - 1.0) / 12.0)
    temperature = np.sin(2.0 * np.pi * (time - 4.0) / 12.0)
    pumping = np.zeros(config.n_steps, dtype=float)
    ramp_start = min(48, config.n_steps)
    ramp_stop = min(60, config.n_steps)
    if ramp_stop > ramp_start:
        pumping[ramp_start:ramp_stop] = np.linspace(
            0.25, 1.0, ramp_stop - ramp_start
        )
    pumping[ramp_stop:] = 1.0

    anomaly = np.zeros(config.n_steps, dtype=float)
    for step in range(1, config.n_steps):
        anomaly[step] = 0.62 * anomaly[step - 1] + rng.normal(0.0, 0.38)
    # A wet-season anomaly is not encoded in the simplified twin.
    pulse = np.array([0.45, 0.80, 0.55, 0.25])
    pulse_stop = min(55, config.n_steps)
    if pulse_stop > 51:
        anomaly[51:pulse_stop] += pulse[: pulse_stop - 51]

    realized = np.column_stack(
        (seasonal_recharge + anomaly, pumping, temperature)
    )
    forecast = np.column_stack((seasonal_recharge, pumping, temperature))
    return realized, forecast


def advance_ensemble(
    ensemble: np.ndarray,
    forcing: np.ndarray,
    model: DynamicModel,
    rng: np.random.Generator,
    *,
    time: int,
    stochastic: bool = True,
) -> np.ndarray:
    """Advance an ensemble one monthly cycle."""

    upstream = np.einsum("vu,euk->evk", model.upstream_weights, ensemble)
    next_state = (
        model.retention[None, None, :] * ensemble
        + model.coupling[None, None, :] * upstream
    )
    forced = forcing @ model.input_coefficients.T
    next_state += model.node_sensitivity[None, :, :] * forced[None, None, :]

    if model.nonlinear_recharge:
        vulnerability = np.linspace(1.20, 0.70, ensemble.shape[1])[None, :]
        recharge_term = model.nonlinear_recharge * float(forcing[0] ** 2)
        next_state[:, :, 0] += recharge_term * vulnerability
        next_state[:, :, 2] += 0.65 * recharge_term * vulnerability

    if model.hidden_pulse and 52 <= time <= 56:
        # Unmodelled local weathering/contamination pulse: the operational model
        # can learn it only through observations, never from its transition law.
        next_state[:, 16:23, 1] += 0.10
        next_state[:, 16:23, 2] += 0.16

    if stochastic:
        noise = rng.normal(
            0.0,
            model.process_sd[None, None, :],
            size=next_state.shape,
        )
        next_state += noise
    return np.clip(next_state, -5.0, 5.0)


def observation_schedule(config: TwinConfig) -> tuple[int, ...]:
    """Monitoring network selected before simulating any data."""

    del config
    return (1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)


def observation_noise() -> np.ndarray:
    """One-standard-deviation observation errors in physical units."""

    return np.array([0.30, 0.018, 0.014], dtype=float)


def _make_observation_batches(
    truth_physical: np.ndarray,
    config: TwinConfig,
    rng: np.random.Generator,
) -> tuple[ObservationBatch, ...]:
    monitors = observation_schedule(config)
    noise = observation_noise()
    batches: list[ObservationBatch] = []
    for time in range(config.n_steps):
        indices: list[int] = []
        values: list[float] = []
        standard_deviations: list[float] = []
        for node in monitors:
            # Heads are available monthly; laboratory chemistry is staggered.
            variables: Iterable[int] = (0,)
            if time % 2 == node % 2:
                variables = (*variables, 1)
            if time % 3 == node % 3:
                variables = (*variables, 2)
            for variable in variables:
                indices.append(node * N_VARIABLES + variable)
                values.append(
                    float(truth_physical[time, node, variable] + rng.normal(0.0, noise[variable]))
                )
                standard_deviations.append(float(noise[variable]))
        batches.append(
            ObservationBatch(
                time=time,
                flat_indices=np.asarray(indices, dtype=int),
                values=np.asarray(values, dtype=float),
                standard_deviations=np.asarray(standard_deviations, dtype=float),
            )
        )
    return tuple(batches)


def simulate_experiment(config: TwinConfig, seed: int) -> ExperimentData:
    """Create one known-truth operational experiment."""

    config.validate()
    rng = np.random.default_rng(seed)
    truth_model, _, _ = build_models(config)
    realized, forecast = generate_forcing(config, rng)
    truth = np.zeros((config.n_steps, config.n_nodes, N_VARIABLES), dtype=float)
    truth[0] = rng.normal(0.0, 0.18, size=(config.n_nodes, N_VARIABLES))
    current = truth[0][None, :, :]
    for time in range(1, config.n_steps):
        current = advance_ensemble(
            current,
            realized[time],
            truth_model,
            rng,
            time=time,
            stochastic=True,
        )
        truth[time] = current[0]
    physical = to_physical(truth)
    observations = _make_observation_batches(physical, config, rng)
    return ExperimentData(
        truth_normalized=truth,
        truth_physical=physical,
        realized_forcing=realized,
        forecast_forcing=forecast,
        observations=observations,
        monitor_nodes=observation_schedule(config),
    )


def initial_ensemble(config: TwinConfig, rng: np.random.Generator) -> np.ndarray:
    """Biased, uncertain prior shared by all dynamic-model comparisons."""

    progress = np.linspace(-0.15, 0.15, config.n_nodes)[:, None]
    bias = np.hstack((progress, -0.6 * progress, 0.8 * progress))
    return bias[None, :, :] + rng.normal(
        0.0,
        np.array([0.55, 0.50, 0.50])[None, None, :],
        size=(config.ensemble_size, config.n_nodes, N_VARIABLES),
    )


def _transform_observations(
    batch: ObservationBatch,
    *,
    mode: str,
    rng: np.random.Generator,
    dropout_fraction: float,
) -> ObservationBatch:
    """Apply predeclared negative-control transformations."""

    indices = batch.flat_indices.copy()
    values = batch.values.copy()
    sd = batch.standard_deviations.copy()
    if mode == "shuffled":
        for variable in range(N_VARIABLES):
            positions = np.flatnonzero(indices % N_VARIABLES == variable)
            if positions.size > 1:
                values[positions] = np.roll(values[positions], 1)
    elif mode == "dropout" and indices.size:
        keep = rng.random(indices.size) >= dropout_fraction
        # Keep at least one head observation so every cycle remains operational.
        if not np.any(keep):
            keep[int(np.flatnonzero(indices % N_VARIABLES == 0)[0])] = True
        indices, values, sd = indices[keep], values[keep], sd[keep]
    return ObservationBatch(batch.time, indices, values, sd)


def enkf_update(
    forecast: np.ndarray,
    batch: ObservationBatch,
    config: TwinConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    """Perturbed-observation EnKF update in normalized state coordinates."""

    if batch.flat_indices.size == 0:
        return forecast.copy()
    n_ensemble = forecast.shape[0]
    flat = forecast.reshape(n_ensemble, -1)
    base = physical_baseline(config.n_nodes).reshape(-1)[batch.flat_indices]
    scale = physical_scale(config.n_nodes).reshape(-1)[batch.flat_indices]
    simulated = base[None, :] + flat[:, batch.flat_indices] * scale[None, :]

    x_anomaly = flat - flat.mean(axis=0, keepdims=True)
    y_anomaly = simulated - simulated.mean(axis=0, keepdims=True)
    denominator = max(n_ensemble - 1, 1)
    cross_covariance = (x_anomaly.T @ y_anomaly) / denominator
    observation_covariance = (y_anomaly.T @ y_anomaly) / denominator
    observation_covariance += np.diag(batch.standard_deviations**2)
    gain = cross_covariance @ np.linalg.pinv(observation_covariance, rcond=1e-10)

    perturbed = batch.values[None, :] + rng.normal(
        0.0,
        batch.standard_deviations[None, :],
        size=(n_ensemble, batch.values.size),
    )
    updated = flat + (perturbed - simulated) @ gain.T
    mean = updated.mean(axis=0, keepdims=True)
    updated = mean + config.inflation * (updated - mean)
    return np.clip(updated.reshape(forecast.shape), -5.0, 5.0)


def _state_diagnostic(
    method: str,
    time: int,
    prior: np.ndarray,
    posterior: np.ndarray,
    truth: np.ndarray,
    n_observations: int,
    interval: float,
) -> dict[str, float | int | str]:
    truth_physical = to_physical(truth[None, :, :])[0]
    prior_physical = to_physical(prior)
    posterior_physical = to_physical(posterior)
    alpha = (1.0 - interval) / 2.0
    lower = np.quantile(posterior_physical, alpha, axis=0)
    upper = np.quantile(posterior_physical, 1.0 - alpha, axis=0)
    return {
        "method": method,
        "time": time,
        "n_observations": n_observations,
        "prior_rmse": float(np.sqrt(np.mean((prior_physical.mean(axis=0) - truth_physical) ** 2))),
        "posterior_rmse": float(
            np.sqrt(np.mean((posterior_physical.mean(axis=0) - truth_physical) ** 2))
        ),
        "posterior_spread": float(np.mean(np.std(posterior_physical, axis=0, ddof=1))),
        "posterior_coverage": float(np.mean((truth_physical >= lower) & (truth_physical <= upper))),
    }


def run_filter(
    data: ExperimentData,
    model: DynamicModel,
    config: TwinConfig,
    *,
    method: str,
    seed: int,
    observation_mode: str = "normal",
    stop_assimilation_at: int | None = None,
) -> FilterRun:
    """Run sequential forecasts and updates without accessing future observations."""

    ensemble = initial_ensemble(config, np.random.default_rng(seed))
    histories: list[np.ndarray] = []
    diagnostics: list[dict[str, float | int | str]] = []
    for time in range(config.n_steps):
        prior = ensemble.copy()
        if time > 0:
            prior = advance_ensemble(
                ensemble,
                data.realized_forcing[time],
                model,
                np.random.default_rng(seed + 10_000 + time),
                time=time,
                stochastic=True,
            )
        can_assimilate = stop_assimilation_at is None or time <= stop_assimilation_at
        if observation_mode == "none" or not can_assimilate:
            posterior = prior.copy()
            n_observations = 0
        else:
            batch = _transform_observations(
                data.observations[time],
                mode=observation_mode,
                rng=np.random.default_rng(seed + 20_000 + time),
                dropout_fraction=config.dropout_fraction,
            )
            posterior = enkf_update(
                prior,
                batch,
                config,
                np.random.default_rng(seed + 30_000 + time),
            )
            n_observations = int(batch.values.size)
        diagnostics.append(
            _state_diagnostic(
                method,
                time,
                prior,
                posterior,
                data.truth_normalized[time],
                n_observations,
                config.interval,
            )
        )
        histories.append(posterior.copy())
        ensemble = posterior
    return FilterRun(method, tuple(histories), pd.DataFrame(diagnostics))


def forecast_ensemble(
    posterior: np.ndarray,
    data: ExperimentData,
    model: DynamicModel,
    *,
    origin: int,
    horizon: int,
    seed: int,
) -> np.ndarray:
    """Issue a forecast using no observations later than ``origin``."""

    rng = np.random.default_rng(seed)
    ensemble = posterior.copy()
    for target in range(origin + 1, origin + horizon + 1):
        ensemble = advance_ensemble(
            ensemble,
            data.forecast_forcing[target],
            model,
            rng,
            time=target,
            stochastic=True,
        )
    return ensemble


def ensemble_crps(samples: np.ndarray, truth: np.ndarray) -> float:
    """Mean continuous ranked probability score for ensemble predictions."""

    samples = np.asarray(samples, dtype=float)
    truth = np.asarray(truth, dtype=float)
    first = np.mean(np.abs(samples - truth[None, :]), axis=0)
    ordered = np.sort(samples, axis=0)
    ensemble_size = ordered.shape[0]
    coefficients = 2 * np.arange(ensemble_size) - ensemble_size + 1
    pairwise = (
        2.0
        * np.sum(coefficients[:, None] * ordered, axis=0)
        / float(ensemble_size**2)
    )
    return float(np.mean(first - 0.5 * pairwise))


def _metric_row(
    prediction: np.ndarray,
    truth: np.ndarray,
    *,
    replicate: int,
    method: str,
    origin: int,
    horizon: int,
    variable: int,
    domain: str,
    interval: float,
    raw_prediction: np.ndarray | None = None,
    calibration_factor: float = 1.0,
) -> dict[str, float | int | str]:
    """Metrics for one ensemble, variable, forecast origin, and node domain."""

    mean = prediction.mean(axis=0)
    error = mean - truth
    alpha = (1.0 - interval) / 2.0
    lower = np.quantile(prediction, alpha, axis=0)
    upper = np.quantile(prediction, 1.0 - alpha, axis=0)
    raw_coverage = np.nan
    if raw_prediction is not None:
        raw_lower = np.quantile(raw_prediction, alpha, axis=0)
        raw_upper = np.quantile(raw_prediction, 1.0 - alpha, axis=0)
        raw_coverage = float(
            np.mean((truth >= raw_lower) & (truth <= raw_upper))
        )
    return {
        "replicate": replicate,
        "method": method,
        "origin": origin,
        "horizon": horizon,
        "variable": VARIABLES[variable],
        "domain": domain,
        "n_targets": int(truth.size),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "mae": float(np.mean(np.abs(error))),
        "bias": float(np.mean(error)),
        "coverage90": float(np.mean((truth >= lower) & (truth <= upper))),
        "raw_coverage90": raw_coverage,
        "interval_width90": float(np.mean(upper - lower)),
        "calibration_factor": float(calibration_factor),
        "crps": ensemble_crps(prediction, truth),
    }


def calibrate_forecast_spread(
    data: ExperimentData,
    run: FilterRun,
    model: DynamicModel,
    config: TwinConfig,
    *,
    seed: int,
    evaluation_origins: Sequence[int] | None = None,
) -> dict[tuple[int, int, int], float]:
    """Calibrate spread from monitoring residuals available at issue time.

    At each evaluation origin, the calibration pool contains only forecasts whose
    verifying observations have already arrived. This is empirical rolling split
    calibration, not a claim of exchangeable conformal coverage: groundwater time
    series are dependent.
    """

    alpha = (1.0 - config.interval) / 2.0
    if evaluation_origins is None:
        evaluation_origins = tuple(
            range(
                config.forecast_start,
                config.n_steps - max(config.horizons),
                3,
            )
        )
    first_origin = max(config.spinup_steps + max(config.horizons), 12)
    factors: dict[tuple[int, int, int], float] = {}
    for evaluation_origin in evaluation_origins:
        ratios: dict[tuple[int, int], list[float]] = {
            (horizon, variable): []
            for horizon in config.horizons
            for variable in range(N_VARIABLES)
        }
        for calibration_origin in range(first_origin, evaluation_origin, 3):
            for horizon in config.horizons:
                target = calibration_origin + horizon
                if target > evaluation_origin:
                    continue
                forecast = to_physical(
                    forecast_ensemble(
                        run.posterior_ensembles[calibration_origin],
                        data,
                        model,
                        origin=calibration_origin,
                        horizon=horizon,
                        seed=seed + calibration_origin * 101 + horizon,
                    )
                ).reshape(config.ensemble_size, -1)
                batch = data.observations[target]
                for variable in range(N_VARIABLES):
                    positions = np.flatnonzero(
                        batch.flat_indices % N_VARIABLES == variable
                    )
                    if positions.size == 0:
                        continue
                    indices = batch.flat_indices[positions]
                    predicted = forecast[:, indices]
                    mean = predicted.mean(axis=0)
                    lower = np.quantile(predicted, alpha, axis=0)
                    upper = np.quantile(predicted, 1.0 - alpha, axis=0)
                    half_width = np.maximum((upper - lower) / 2.0, 1e-10)
                    standardized = (
                        np.abs(batch.values[positions] - mean) / half_width
                    )
                    ratios[(horizon, variable)].extend(standardized.tolist())

        for key, values in ratios.items():
            horizon, variable = key
            if not values:
                factors[(evaluation_origin, horizon, variable)] = 1.0
                continue
            factors[(evaluation_origin, horizon, variable)] = max(
                1.0,
                float(
                    np.quantile(
                        np.asarray(values),
                        config.interval,
                        method="higher",
                    )
                )
            )
    return factors


def evaluate_filter_forecasts(
    data: ExperimentData,
    run: FilterRun,
    model: DynamicModel,
    config: TwinConfig,
    *,
    replicate: int,
    seed: int,
) -> pd.DataFrame:
    """Evaluate prequential forecasts at predeclared origins and horizons."""

    monitor = np.asarray(data.monitor_nodes, dtype=int)
    unmonitored = np.asarray(
        [node for node in range(config.n_nodes) if node not in data.monitor_nodes],
        dtype=int,
    )
    domains = {
        "all": np.arange(config.n_nodes, dtype=int),
        "monitored": monitor,
        "unmonitored": unmonitored,
    }
    rows: list[dict[str, float | int | str]] = []
    max_horizon = max(config.horizons)
    origins = tuple(range(config.forecast_start, config.n_steps - max_horizon, 3))
    calibration = calibrate_forecast_spread(
        data,
        run,
        model,
        config,
        seed=seed,
        evaluation_origins=origins,
    )
    for origin in origins:
        for horizon in config.horizons:
            target = origin + horizon
            forecast = to_physical(
                forecast_ensemble(
                    run.posterior_ensembles[origin],
                    data,
                    model,
                    origin=origin,
                    horizon=horizon,
                    seed=seed + origin * 101 + horizon,
                )
            )
            truth = data.truth_physical[target]
            for variable in range(N_VARIABLES):
                factor = calibration[(origin, horizon, variable)]
                raw_variable_forecast = forecast[:, :, variable]
                forecast_mean = raw_variable_forecast.mean(axis=0, keepdims=True)
                calibrated_variable_forecast = forecast_mean + factor * (
                    raw_variable_forecast - forecast_mean
                )
                for domain, nodes in domains.items():
                    rows.append(
                        _metric_row(
                            calibrated_variable_forecast[:, nodes],
                            truth[nodes, variable],
                            replicate=replicate,
                            method=run.method,
                            origin=origin,
                            horizon=horizon,
                            variable=variable,
                            domain=domain,
                            interval=config.interval,
                            raw_prediction=raw_variable_forecast[:, nodes],
                            calibration_factor=factor,
                        )
                    )
    return pd.DataFrame(rows)


def evaluate_deterministic_baselines(
    data: ExperimentData,
    config: TwinConfig,
    *,
    replicate: int,
) -> pd.DataFrame:
    """Conservative oracle-persistence and seasonal-climatology comparisons."""

    monitor = np.asarray(data.monitor_nodes, dtype=int)
    unmonitored = np.asarray(
        [node for node in range(config.n_nodes) if node not in data.monitor_nodes],
        dtype=int,
    )
    domains = {
        "all": np.arange(config.n_nodes, dtype=int),
        "monitored": monitor,
        "unmonitored": unmonitored,
    }
    rows: list[dict[str, float | int | str]] = []
    max_horizon = max(config.horizons)
    origins = range(config.forecast_start, config.n_steps - max_horizon, 3)
    training = data.truth_physical[: config.forecast_start]
    for origin in origins:
        for horizon in config.horizons:
            target = origin + horizon
            target_month = target % 12
            month_indices = np.arange(target_month, config.forecast_start, 12)
            climatology = training[month_indices].mean(axis=0)
            baseline_predictions = {
                # Oracle persistence knows the exact current synthetic state.  It
                # is intentionally stronger than an operational last-observation
                # baseline and therefore conservative for claims of skill.
                "oracle_persistence": data.truth_physical[origin],
                "oracle_seasonal_climatology": climatology,
            }
            truth = data.truth_physical[target]
            for method, deterministic in baseline_predictions.items():
                for variable in range(N_VARIABLES):
                    for domain, nodes in domains.items():
                        prediction = deterministic[nodes, variable][None, :]
                        row = _metric_row(
                            prediction,
                            truth[nodes, variable],
                            replicate=replicate,
                            method=method,
                            origin=origin,
                            horizon=horizon,
                            variable=variable,
                            domain=domain,
                            interval=config.interval,
                        )
                        row["coverage90"] = np.nan
                        row["raw_coverage90"] = np.nan
                        row["interval_width90"] = np.nan
                        row["calibration_factor"] = np.nan
                        rows.append(row)
    return pd.DataFrame(rows)


def representative_timeseries(
    data: ExperimentData,
    updated: FilterRun,
    open_loop: FilterRun,
    config: TwinConfig,
    *,
    replicate: int,
) -> pd.DataFrame:
    """Long-form state record for visual and audit inspection."""

    observation_lookup: dict[tuple[int, int], float] = {}
    for batch in data.observations:
        for flat_index, value in zip(batch.flat_indices, batch.values):
            observation_lookup[(batch.time, int(flat_index))] = float(value)

    alpha = (1.0 - config.interval) / 2.0
    rows: list[dict[str, float | int | str]] = []
    for time in range(config.n_steps):
        updated_physical = to_physical(updated.posterior_ensembles[time])
        open_physical = to_physical(open_loop.posterior_ensembles[time])
        for node in range(config.n_nodes):
            for variable in range(N_VARIABLES):
                flat_index = node * N_VARIABLES + variable
                rows.append(
                    {
                        "replicate": replicate,
                        "time": time,
                        "node": f"N{node:02d}",
                        "variable": VARIABLES[variable],
                        "truth": data.truth_physical[time, node, variable],
                        "observation": observation_lookup.get((time, flat_index), np.nan),
                        "updated_mean": updated_physical[:, node, variable].mean(),
                        "updated_q05": np.quantile(
                            updated_physical[:, node, variable], alpha
                        ),
                        "updated_q95": np.quantile(
                            updated_physical[:, node, variable], 1.0 - alpha
                        ),
                        "open_loop_mean": open_physical[:, node, variable].mean(),
                        "forecast_evaluation_period": int(time >= config.forecast_start),
                    }
                )
    return pd.DataFrame(rows)


def run_replicate(
    config: TwinConfig,
    replicate: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """Run all predeclared methods on one shared truth/observation realization."""

    replicate_seed = config.seed + replicate * 10_003
    data = simulate_experiment(config, replicate_seed)
    _, nominal, wrong = build_models(config)
    method_specs = (
        ("updated_twin", nominal, "normal"),
        ("open_loop", nominal, "none"),
        ("wrong_topology_updated", wrong, "normal"),
        ("shuffled_observation_updated", nominal, "shuffled"),
        ("sensor_dropout_updated", nominal, "dropout"),
    )
    runs: dict[str, FilterRun] = {}
    metrics: list[pd.DataFrame] = []
    diagnostics: list[pd.DataFrame] = []
    for method, model, observation_mode in method_specs:
        run = run_filter(
            data,
            model,
            config,
            method=method,
            seed=replicate_seed + 1_000,
            observation_mode=observation_mode,
        )
        run_diag = run.diagnostics.copy()
        run_diag.insert(0, "replicate", replicate)
        diagnostics.append(run_diag)
        metrics.append(
            evaluate_filter_forecasts(
                data,
                run,
                model,
                config,
                replicate=replicate,
                seed=replicate_seed + 5_000,
            )
        )
        runs[method] = run
    metrics.append(evaluate_deterministic_baselines(data, config, replicate=replicate))
    timeseries = None
    if replicate == 0:
        timeseries = representative_timeseries(
            data,
            runs["updated_twin"],
            runs["open_loop"],
            config,
            replicate=replicate,
        )
    return pd.concat(metrics, ignore_index=True), pd.concat(diagnostics, ignore_index=True), timeseries


def summarize_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    """Aggregate prequential metrics without treating origins as replicates."""

    replicate_level = (
        metrics.groupby(
            ["replicate", "method", "horizon", "variable", "domain"],
            as_index=False,
        )
        .agg(
            rmse=("rmse", "mean"),
            mae=("mae", "mean"),
            bias=("bias", "mean"),
            coverage90=("coverage90", "mean"),
            raw_coverage90=("raw_coverage90", "mean"),
            interval_width90=("interval_width90", "mean"),
            calibration_factor=("calibration_factor", "mean"),
            crps=("crps", "mean"),
        )
    )
    rows: list[dict[str, float | int | str]] = []
    grouping = ["method", "horizon", "variable", "domain"]
    for keys, group in replicate_level.groupby(grouping, sort=True):
        method, horizon, variable, domain = keys
        row: dict[str, float | int | str] = {
            "method": method,
            "horizon": int(horizon),
            "variable": variable,
            "domain": domain,
            "n_replicates": int(group["replicate"].nunique()),
        }
        for metric in (
            "rmse",
            "mae",
            "bias",
            "coverage90",
            "raw_coverage90",
            "interval_width90",
            "calibration_factor",
            "crps",
        ):
            values = group[metric].dropna().to_numpy(dtype=float)
            row[f"{metric}_mean"] = float(np.mean(values)) if values.size else np.nan
            row[f"{metric}_sd"] = (
                float(np.std(values, ddof=1)) if values.size > 1 else np.nan
            )
            row[f"{metric}_q025"] = (
                float(np.quantile(values, 0.025)) if values.size else np.nan
            )
            row[f"{metric}_q975"] = (
                float(np.quantile(values, 0.975)) if values.size else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def paired_skill(metrics: pd.DataFrame) -> pd.DataFrame:
    """Replicate-level paired skill against declared comparators."""

    replicate_level = (
        metrics.query("domain == 'all'")
        .groupby(["replicate", "method", "horizon", "variable"], as_index=False)
        .agg(rmse=("rmse", "mean"), crps=("crps", "mean"))
    )
    updated = replicate_level.query("method == 'updated_twin'")
    comparators = (
        "open_loop",
        "oracle_persistence",
        "oracle_seasonal_climatology",
        "wrong_topology_updated",
        "shuffled_observation_updated",
        "sensor_dropout_updated",
    )
    rows: list[dict[str, float | int | str]] = []
    for comparator in comparators:
        other = replicate_level.query("method == @comparator")
        merged = updated.merge(
            other,
            on=["replicate", "horizon", "variable"],
            suffixes=("_updated", "_comparator"),
        )
        for (horizon, variable), group in merged.groupby(["horizon", "variable"]):
            skill = 1.0 - group["rmse_updated"] / np.maximum(
                group["rmse_comparator"], 1e-12
            )
            difference = group["rmse_comparator"] - group["rmse_updated"]
            rows.append(
                {
                    "comparator": comparator,
                    "horizon": int(horizon),
                    "variable": variable,
                    "n_replicates": int(len(group)),
                    "mean_rmse_skill": float(skill.mean()),
                    "skill_q025": float(skill.quantile(0.025)),
                    "skill_q975": float(skill.quantile(0.975)),
                    "mean_paired_rmse_reduction": float(difference.mean()),
                    "reduction_q025": float(difference.quantile(0.025)),
                    "reduction_q975": float(difference.quantile(0.975)),
                    "fraction_replicates_better": float((difference > 0.0).mean()),
                }
            )
    return pd.DataFrame(rows)


def config_as_dict(config: TwinConfig) -> dict[str, object]:
    """JSON-friendly representation of the locked protocol."""

    result = asdict(config)
    result["horizons"] = list(config.horizons)
    result["variables"] = list(VARIABLES)
    result["edges"] = [list(edge) for edge in locked_edges()]
    result["monitor_nodes"] = list(observation_schedule(config))
    result["observation_sd"] = observation_noise().tolist()
    return result
