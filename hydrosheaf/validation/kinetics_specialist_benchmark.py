"""Truth-blind synthetic benchmark for the M8 kinetic specialist.

This module is deliberately independent of the programme runner.  It tests the
kinetic specialist's inverse-facing behaviour with several generator families,
transport/front observations, noise, missingness, predictive intervals and a
selective-risk gate.  The synthetic observation model uses the same structural
quantity exposed by the PHREEQC adapter, ``k * surface_area``.  Consequently a
rate constant and a surface area cannot be separated from chemistry/transport
observations alone; the benchmark reports that fact and abstains from a
separate ``k``/``A`` claim unless an independent area measurement is present.

The fixed effective-rate comparator is a competence-matched diagnostic only. It
is not a PHREEQC inverse solver and must not be used as evidence of field or
general specialist superiority.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np


BENCHMARK_VERSION = "m8-kinetics-synthetic-v1"
COMPARATOR_STATUS = "DIAGNOSTIC_ONLY"
COMPETENCE_MATCHED_COMPARATOR_STATUS = "COMPETENCE_MATCHED_DIAGNOSTIC_ONLY"


@dataclass(frozen=True)
class KineticRegime:
    """Independent synthetic generator specification."""

    name: str
    generator_family: str
    k_per_day: float
    surface_area: float
    velocity_m_per_day: float
    front_delay_days: float
    distance_m: float = 45.0


DEFAULT_REGIMES: Tuple[KineticRegime, ...] = (
    KineticRegime(
        "transport_limited",
        "exponential_front",
        0.060,
        2.20,
        1.05,
        5.0,
    ),
    KineticRegime(
        "reaction_limited",
        "stretched_exponential",
        0.028,
        3.40,
        1.55,
        7.0,
    ),
    KineticRegime(
        "mixed_kinetics",
        "two_stage_reaction",
        0.082,
        1.45,
        0.78,
        12.0,
    ),
)


@dataclass(frozen=True)
class KineticsBenchmarkConfig:
    """Hashable, deterministic design and scoring contract."""

    seed: int = 20260802
    cases_per_regime: int = 4
    n_time_points: int = 12
    time_start_days: float = 4.0
    time_end_days: float = 72.0
    front_noise_sd_m: float = 0.08
    concentration_noise_sd: float = 0.025
    missingness_rate: float = 0.15
    surface_area_measurement_probability: float = 0.50
    surface_area_measurement_sd: float = 0.12
    confidence_level: float = 0.90
    ka_grid_min_per_day: float = 1.0e-4
    ka_grid_max_per_day: float = 0.50
    ka_grid_size: int = 320
    min_concentration_points: int = 3
    max_relative_interval_width: float = 1.50
    interval_discrepancy_inflation: float = 2.50
    development_fraction: float = 0.50
    regimes: Tuple[KineticRegime, ...] = DEFAULT_REGIMES
    benchmark_version: str = BENCHMARK_VERSION

    def __post_init__(self) -> None:
        if self.cases_per_regime < 1 or self.n_time_points < 5:
            raise ValueError("cases_per_regime must be >=1 and n_time_points >=5")
        if self.time_end_days <= self.time_start_days:
            raise ValueError("time_end_days must exceed time_start_days")
        if not 0.0 <= self.missingness_rate < 1.0:
            raise ValueError("missingness_rate must be in [0, 1)")
        if not 0.0 <= self.surface_area_measurement_probability <= 1.0:
            raise ValueError("surface_area_measurement_probability must be in [0, 1]")
        if self.front_noise_sd_m <= 0 or self.concentration_noise_sd <= 0:
            raise ValueError("noise scales must be positive")
        if not 0.0 < self.confidence_level < 1.0:
            raise ValueError("confidence_level must be in (0, 1)")
        if self.ka_grid_min_per_day <= 0 or self.ka_grid_max_per_day <= self.ka_grid_min_per_day:
            raise ValueError("invalid effective-rate grid")
        if self.ka_grid_size < 20:
            raise ValueError("ka_grid_size must be >=20")
        if self.interval_discrepancy_inflation < 1.0:
            raise ValueError("interval_discrepancy_inflation must be >=1")
        if not 0.0 < self.development_fraction < 1.0:
            raise ValueError("development_fraction must be strictly between 0 and 1")

    def as_dict(self) -> Dict[str, Any]:
        value = asdict(self)
        value["regimes"] = [asdict(regime) for regime in self.regimes]
        return value

    @property
    def config_hash(self) -> str:
        payload = json.dumps(self.as_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class KineticsObservation:
    """Public, truth-blind observation record for one synthetic case."""

    case_id: str
    regime: str
    split: str
    times_days: np.ndarray
    front_m: np.ndarray
    front_observed: np.ndarray
    concentration: np.ndarray
    concentration_observed: np.ndarray
    surface_area_measurement: Optional[float]

    def public_dict(self) -> Dict[str, Any]:
        return {
            "case_id": self.case_id,
            "regime": self.regime,
            "split": self.split,
            "times_days": self.times_days.tolist(),
            "front_m": [None if not bool(mask) else float(value) for value, mask in zip(self.front_m, self.front_observed)],
            "concentration": [None if not bool(mask) else float(value) for value, mask in zip(self.concentration, self.concentration_observed)],
            "surface_area_measurement": self.surface_area_measurement,
            "truth_blind": True,
        }


@dataclass(frozen=True)
class _KineticsTruth:
    case_id: str
    regime: str
    generator_family: str
    k_per_day: float
    surface_area: float
    velocity_m_per_day: float
    front_delay_days: float
    effective_rate_per_day: float
    true_front_m: np.ndarray
    true_concentration: np.ndarray


@dataclass
class KineticsBenchmarkDataset:
    """Dataset with public observations and private scoring truth."""

    observations: Tuple[KineticsObservation, ...]
    config_hash: str
    benchmark_version: str
    development_case_ids: Tuple[str, ...]
    locked_case_ids: Tuple[str, ...]
    _truth: Tuple[_KineticsTruth, ...] = field(repr=False)

    def truth_blind_view(self) -> Tuple[Dict[str, Any], ...]:
        """Return observations only; no latent parameter or noiseless trajectory."""

        return tuple(observation.public_dict() for observation in self.observations)

    def _truth_by_case(self) -> Dict[str, _KineticsTruth]:
        return {truth.case_id: truth for truth in self._truth}

    def as_dict(self) -> Dict[str, Any]:
        return {
            "benchmark_version": self.benchmark_version,
            "config_hash": self.config_hash,
            "n_cases": len(self.observations),
            "development_case_ids": list(self.development_case_ids),
            "locked_case_ids": list(self.locked_case_ids),
            "observations": list(self.truth_blind_view()),
            "truth_released_for_scoring": False,
        }


@dataclass(frozen=True)
class KineticsPrediction:
    """Truth-free inverse and forward prediction for one case."""

    case_id: str
    model: str
    effective_rate_per_day: Optional[float]
    effective_rate_interval: Optional[Tuple[float, float]]
    k_per_day: Optional[float]
    k_interval: Optional[Tuple[float, float]]
    surface_area: Optional[float]
    surface_area_interval: Optional[Tuple[float, float]]
    velocity_m_per_day: Optional[float]
    front_delay_days: Optional[float]
    predicted_front_m: np.ndarray
    predicted_concentration: np.ndarray
    front_interval: Tuple[np.ndarray, np.ndarray]
    concentration_interval: Tuple[np.ndarray, np.ndarray]
    abstain: bool
    parameter_abstain: bool
    k_a_identifiable: bool
    identifiability: str
    n_front_observations: int
    n_concentration_observations: int
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "case_id": self.case_id,
            "model": self.model,
            "effective_rate_per_day": self.effective_rate_per_day,
            "effective_rate_interval": self.effective_rate_interval,
            "k_per_day": self.k_per_day,
            "k_interval": self.k_interval,
            "surface_area": self.surface_area,
            "surface_area_interval": self.surface_area_interval,
            "velocity_m_per_day": self.velocity_m_per_day,
            "front_delay_days": self.front_delay_days,
            "predicted_front_m": self.predicted_front_m.tolist(),
            "predicted_concentration": self.predicted_concentration.tolist(),
            "front_interval": [self.front_interval[0].tolist(), self.front_interval[1].tolist()],
            "concentration_interval": [self.concentration_interval[0].tolist(), self.concentration_interval[1].tolist()],
            "abstain": self.abstain,
            "parameter_abstain": self.parameter_abstain,
            "k_a_identifiable": self.k_a_identifiable,
            "identifiability": self.identifiability,
            "n_front_observations": self.n_front_observations,
            "n_concentration_observations": self.n_concentration_observations,
            "diagnostics": dict(self.diagnostics),
            "truth_blind": True,
        }


@dataclass(frozen=True)
class KineticsBenchmarkScore:
    """Aggregated score with explicit selective and identifiability fields."""

    model: str
    n_cases: int
    n_predictive_cases: int
    abstention_rate: float
    parameter_abstention_rate: float
    predictive_rmse: float
    front_rmse: float
    concentration_rmse: float
    effective_rate_log_rmse: float
    k_log_rmse_identified: Optional[float]
    surface_area_log_rmse_identified: Optional[float]
    effective_rate_interval_coverage: float
    front_interval_coverage: float
    concentration_interval_coverage: float
    selective_effective_rate_log_rmse: Optional[float]
    false_commitment_rate: float
    identified_case_rate: float
    identifiability_counts: Mapping[str, int]
    status: str

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class KineticsBenchmarkReport:
    """Complete benchmark output and claim boundary."""

    benchmark_version: str
    config: Mapping[str, Any]
    config_hash: str
    dataset: KineticsBenchmarkDataset
    specialist_predictions: Tuple[KineticsPrediction, ...]
    comparator_predictions: Tuple[KineticsPrediction, ...]
    interval_calibrator: "KineticsIntervalCalibrator"
    comparator_interval_calibrator: "KineticsIntervalCalibrator"
    specialist_score: KineticsBenchmarkScore
    comparator_score: KineticsBenchmarkScore
    comparator_status: str
    claim_status: str
    claim_boundary: str

    def as_dict(self) -> Dict[str, Any]:
        return {
            "benchmark_version": self.benchmark_version,
            "config": dict(self.config),
            "config_hash": self.config_hash,
            "dataset": self.dataset.as_dict(),
            "specialist_predictions": [prediction.as_dict() for prediction in self.specialist_predictions],
            "comparator_predictions": [prediction.as_dict() for prediction in self.comparator_predictions],
            "interval_calibrator": self.interval_calibrator.as_dict(),
            "comparator_interval_calibrator": self.comparator_interval_calibrator.as_dict(),
            "specialist_score": self.specialist_score.as_dict(),
            "comparator_score": self.comparator_score.as_dict(),
            "comparator_status": self.comparator_status,
            "claim_status": self.claim_status,
            "claim_boundary": self.claim_boundary,
            "truth_released_for_scoring": False,
        }


@dataclass(frozen=True)
class KineticsIntervalCalibrator:
    """Development-only conformal scale for predictive intervals.

    The calibrator stores only residual-derived scales.  It is fitted using
    development-case truth and then applied to locked-case predictions before
    locked scoring.  No locked truth is needed by ``apply``.
    """

    fit_scope: str
    confidence_level: float
    development_case_ids: Tuple[str, ...]
    locked_case_ids: Tuple[str, ...]
    effective_rate_scale: float
    front_scale: float
    concentration_scale: float
    max_relative_interval_width: float
    calibration_status: str

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def apply(self, predictions: Sequence[KineticsPrediction]) -> Tuple[KineticsPrediction, ...]:
        calibrated: List[KineticsPrediction] = []
        for prediction in predictions:
            effective_interval = prediction.effective_rate_interval
            if prediction.effective_rate_per_day is not None and effective_interval is not None:
                center = prediction.effective_rate_per_day
                lower = max(effective_interval[0], 1.0e-12)
                upper = max(effective_interval[1], center)
                half_width = max(math.log(center / lower), math.log(upper / center))
                effective_interval = (
                    center / math.exp(half_width * self.effective_rate_scale),
                    center * math.exp(half_width * self.effective_rate_scale),
                )
            front_center = prediction.predicted_front_m
            front_radius = (prediction.front_interval[1] - prediction.front_interval[0]) / 2.0
            front_interval = (
                front_center - front_radius * self.front_scale,
                front_center + front_radius * self.front_scale,
            )
            concentration_center = prediction.predicted_concentration
            concentration_radius = (
                prediction.concentration_interval[1] - prediction.concentration_interval[0]
            ) / 2.0
            concentration_interval = (
                np.clip(concentration_center - concentration_radius * self.concentration_scale, 0.0, 1.0),
                np.clip(concentration_center + concentration_radius * self.concentration_scale, 0.0, 1.0),
            )
            calibrated_relative_width = float("inf")
            if prediction.effective_rate_per_day is not None and effective_interval is not None:
                calibrated_relative_width = (
                    effective_interval[1] - effective_interval[0]
                ) / prediction.effective_rate_per_day
            abstain = prediction.abstain or calibrated_relative_width > self.max_relative_interval_width
            k_interval = prediction.k_interval
            if (
                prediction.k_a_identifiable
                and effective_interval is not None
                and prediction.surface_area_interval is not None
            ):
                area_lower = max(prediction.surface_area_interval[0], 1.0e-12)
                area_upper = max(prediction.surface_area_interval[1], area_lower)
                k_interval = (
                    effective_interval[0] / area_upper,
                    effective_interval[1] / area_lower,
                )
            calibrated.append(
                replace(
                    prediction,
                    effective_rate_interval=effective_interval,
                    k_interval=k_interval,
                    front_interval=front_interval,
                    concentration_interval=concentration_interval,
                    abstain=abstain,
                    diagnostics={
                        **dict(prediction.diagnostics),
                        "interval_calibration_fit_scope": self.fit_scope,
                        "effective_rate_calibration_scale": self.effective_rate_scale,
                        "front_calibration_scale": self.front_scale,
                        "concentration_calibration_scale": self.concentration_scale,
                        "post_calibration_effective_rate_relative_interval_width": calibrated_relative_width,
                    },
                )
            )
        return tuple(calibrated)


def _regime_map(config: KineticsBenchmarkConfig) -> Dict[str, KineticRegime]:
    return {regime.name: regime for regime in config.regimes}


def _normal_quantile(confidence_level: float) -> float:
    # Acklam-free fixed approximation is enough for the declared 0.80--0.99
    # benchmark interval range and avoids introducing a scipy-only path.
    table = {0.80: 1.2816, 0.90: 1.6449, 0.95: 1.9600, 0.99: 2.5758}
    nearest = min(table, key=lambda value: abs(value - confidence_level))
    return table[nearest]


def _progress(family: str, effective_rate: float, elapsed_days: np.ndarray) -> np.ndarray:
    elapsed = np.maximum(np.asarray(elapsed_days, dtype=float), 0.0)
    x = np.maximum(effective_rate * elapsed, 0.0)
    if family == "stretched_exponential":
        return 1.0 - np.exp(-(x ** 0.88))
    if family == "two_stage_reaction":
        first = 1.0 - np.exp(-x)
        second = 0.72 + 0.28 * (1.0 - np.exp(-elapsed / 22.0))
        return first * second
    if family != "exponential_front":
        raise ValueError(f"Unknown kinetic generator family: {family}")
    return 1.0 - np.exp(-x)


def _forward(
    family: str,
    effective_rate: float,
    velocity: float,
    delay: float,
    times_days: np.ndarray,
    distance_m: float,
) -> Tuple[np.ndarray, np.ndarray]:
    elapsed = np.maximum(np.asarray(times_days, dtype=float) - delay, 0.0)
    front = np.minimum(distance_m, velocity * elapsed)
    concentration = _progress(family, effective_rate, elapsed)
    return front, concentration


def generate_kinetics_dataset(config: KineticsBenchmarkConfig = KineticsBenchmarkConfig()) -> KineticsBenchmarkDataset:
    """Generate independent cases with hidden truth and public noisy observations."""

    times = np.linspace(config.time_start_days, config.time_end_days, config.n_time_points)
    observations: List[KineticsObservation] = []
    truth_records: List[_KineticsTruth] = []
    development_case_ids: List[str] = []
    locked_case_ids: List[str] = []
    n_development_cases = max(
        1,
        min(
            config.cases_per_regime - 1,
            int(math.ceil(config.cases_per_regime * config.development_fraction)),
        ),
    )

    for regime_index, regime in enumerate(config.regimes):
        for case_index in range(config.cases_per_regime):
            # Independent SeedSequence streams prevent one regime's number of
            # cases or missingness pattern from changing another regime.
            seed_sequence = np.random.SeedSequence([config.seed, regime_index, case_index])
            rng = np.random.default_rng(seed_sequence)
            front_true, concentration_true = _forward(
                regime.generator_family,
                regime.k_per_day * regime.surface_area,
                regime.velocity_m_per_day,
                regime.front_delay_days,
                times,
                regime.distance_m,
            )
            front_observed = rng.random(config.n_time_points) >= config.missingness_rate
            concentration_observed = rng.random(config.n_time_points) >= config.missingness_rate
            front_observed[0] = True
            concentration_observed[-1] = True
            front_noisy = front_true + rng.normal(0.0, config.front_noise_sd_m, config.n_time_points)
            concentration_noisy = np.clip(
                concentration_true + rng.normal(0.0, config.concentration_noise_sd, config.n_time_points),
                0.0,
                1.0,
            )
            area_measurement: Optional[float] = None
            if rng.random() < config.surface_area_measurement_probability:
                area_measurement = max(
                    1.0e-6,
                    float(regime.surface_area + rng.normal(0.0, config.surface_area_measurement_sd)),
                )
                if rng.random() < config.missingness_rate:
                    area_measurement = None
            case_id = f"{regime.name}-{case_index + 1:02d}"
            split = "development" if case_index < n_development_cases else "locked_test"
            if split == "development":
                development_case_ids.append(case_id)
            else:
                locked_case_ids.append(case_id)
            observations.append(
                KineticsObservation(
                    case_id=case_id,
                    regime=regime.name,
                    split=split,
                    times_days=times.copy(),
                    front_m=np.where(front_observed, front_noisy, np.nan),
                    front_observed=front_observed.copy(),
                    concentration=np.where(concentration_observed, concentration_noisy, np.nan),
                    concentration_observed=concentration_observed.copy(),
                    surface_area_measurement=area_measurement,
                )
            )
            truth_records.append(
                _KineticsTruth(
                    case_id=case_id,
                    regime=regime.name,
                    generator_family=regime.generator_family,
                    k_per_day=regime.k_per_day,
                    surface_area=regime.surface_area,
                    velocity_m_per_day=regime.velocity_m_per_day,
                    front_delay_days=regime.front_delay_days,
                    effective_rate_per_day=regime.k_per_day * regime.surface_area,
                    true_front_m=front_true,
                    true_concentration=concentration_true,
                )
            )

    return KineticsBenchmarkDataset(
        observations=tuple(observations),
        config_hash=config.config_hash,
        benchmark_version=config.benchmark_version,
        development_case_ids=tuple(development_case_ids),
        locked_case_ids=tuple(locked_case_ids),
        _truth=tuple(truth_records),
    )


def _estimate_front(
    observation: KineticsObservation,
    distance_m: float = 45.0,
) -> Tuple[Optional[float], Optional[float], float]:
    mask = observation.front_observed & np.isfinite(observation.front_m)
    times = observation.times_days[mask]
    fronts = observation.front_m[mask]
    if len(times) < 2:
        return None, None, float("nan")
    best: Tuple[float, float, float] = (float("inf"), 0.0, 0.0)
    # Fit the actual front observation model: zero before the front arrives and
    # a distance cap after arrival.  This avoids the systematic early-delay
    # bias caused by treating every observation as a positive straight line.
    delay_grid = np.linspace(0.0, max(float(np.max(times)), 1.0), 100)
    velocity_grid = np.linspace(0.05, max(3.0, float(np.max(fronts)) / max(float(np.max(times)), 1.0) * 3.0), 140)
    for delay in delay_grid:
        elapsed = np.maximum(times - delay, 0.0)
        for velocity in velocity_grid:
            predicted = np.minimum(distance_m, velocity * elapsed)
            residual = float(np.mean((fronts - predicted) ** 2))
            if residual < best[0]:
                best = (residual, float(velocity), float(delay))
    return best[1], best[2], best[0]


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = np.maximum(weights[order], 0.0)
    cumulative = np.cumsum(sorted_weights)
    if cumulative[-1] <= 0:
        return float(sorted_values[len(sorted_values) // 2])
    return float(sorted_values[np.searchsorted(cumulative, quantile * cumulative[-1], side="left")])


def _estimate_effective_rate(
    observation: KineticsObservation,
    delay: Optional[float],
    config: KineticsBenchmarkConfig,
) -> Tuple[Optional[float], Optional[Tuple[float, float]], float, int]:
    mask = observation.concentration_observed & np.isfinite(observation.concentration)
    times = observation.times_days[mask]
    concentrations = observation.concentration[mask]
    if delay is None or len(times) < config.min_concentration_points:
        return None, None, float("nan"), len(times)
    elapsed = np.maximum(times - delay, 0.0)
    grid = np.geomspace(config.ka_grid_min_per_day, config.ka_grid_max_per_day, config.ka_grid_size)
    residual_sse = np.array(
        [
            float(np.sum((concentrations - _progress("exponential_front", rate, elapsed)) ** 2))
            for rate in grid
        ]
    )
    best_index = int(np.argmin(residual_sse))
    best_rate = float(grid[best_index])
    variance = max(config.concentration_noise_sd ** 2, float(residual_sse[best_index]) / max(len(times), 1))
    # The common exponential inverse model is intentionally tested against
    # held-out generator families.  Inflate profile intervals by a frozen
    # discrepancy factor so interval coverage measures both observation noise
    # and declared structural mismatch instead of reporting false precision.
    variance *= config.interval_discrepancy_inflation ** 2
    weights = np.exp(-0.5 * (residual_sse - residual_sse[best_index]) / variance)
    interval = (
        _weighted_quantile(grid, weights, (1.0 - config.confidence_level) / 2.0),
        _weighted_quantile(grid, weights, 1.0 - (1.0 - config.confidence_level) / 2.0),
    )
    return best_rate, interval, float(residual_sse[best_index] / len(times)), len(times)


def _prediction(
    observation: KineticsObservation,
    config: KineticsBenchmarkConfig,
    model: str,
    effective_rate: Optional[float],
    effective_interval: Optional[Tuple[float, float]],
    velocity: Optional[float],
    delay: Optional[float],
    concentration_residual: float,
    front_residual: float,
    area: Optional[float],
    area_interval: Optional[Tuple[float, float]],
    use_area_measurement: bool,
    diagnostics: Mapping[str, Any],
) -> KineticsPrediction:
    z = _normal_quantile(config.confidence_level)
    times = observation.times_days
    safe_rate = float(effective_rate if effective_rate is not None else config.ka_grid_min_per_day)
    safe_velocity = float(velocity if velocity is not None else 0.0)
    safe_delay = float(delay if delay is not None else 0.0)
    regime = next(regime for regime in config.regimes if regime.name == observation.regime)
    predicted_front, predicted_concentration = _forward(
        "exponential_front", safe_rate, safe_velocity, safe_delay, times, regime.distance_m
    )
    front_sd = math.sqrt(config.front_noise_sd_m ** 2 + max(front_residual, 0.0))
    concentration_sd = math.sqrt(config.concentration_noise_sd ** 2 + max(concentration_residual, 0.0))
    front_interval = (predicted_front - z * front_sd, predicted_front + z * front_sd)
    concentration_interval = (
        np.clip(predicted_concentration - z * concentration_sd, 0.0, 1.0),
        np.clip(predicted_concentration + z * concentration_sd, 0.0, 1.0),
    )
    n_front = int(np.sum(observation.front_observed & np.isfinite(observation.front_m)))
    n_concentration = int(np.sum(observation.concentration_observed & np.isfinite(observation.concentration)))
    identifiable = bool(use_area_measurement and effective_rate is not None)
    if identifiable and area is not None and area_interval is not None and effective_interval is not None:
        k_value = float(effective_rate / area)
        k_interval = (float(effective_interval[0] / area_interval[1]), float(effective_interval[1] / area_interval[0]))
        identifiability = "identified_with_independent_surface_area_measurement"
    else:
        k_value = None
        k_interval = None
        identifiability = "structurally_confounded_k_times_A"
    relative_width = float("inf")
    if effective_rate is not None and effective_interval is not None and effective_rate > 0:
        relative_width = (effective_interval[1] - effective_interval[0]) / effective_rate
    abstain = bool(
        effective_rate is None
        or n_concentration < config.min_concentration_points
        or relative_width > config.max_relative_interval_width
    )
    return KineticsPrediction(
        case_id=observation.case_id,
        model=model,
        effective_rate_per_day=effective_rate,
        effective_rate_interval=effective_interval,
        k_per_day=k_value,
        k_interval=k_interval,
        surface_area=area if identifiable else None,
        surface_area_interval=area_interval if identifiable else None,
        velocity_m_per_day=velocity,
        front_delay_days=delay,
        predicted_front_m=predicted_front,
        predicted_concentration=predicted_concentration,
        front_interval=front_interval,
        concentration_interval=concentration_interval,
        abstain=abstain,
        parameter_abstain=not identifiable,
        k_a_identifiable=identifiable,
        identifiability=identifiability,
        n_front_observations=n_front,
        n_concentration_observations=n_concentration,
        diagnostics={**dict(diagnostics), "effective_rate_relative_interval_width": relative_width},
    )


def fit_kinetics_specialist(
    dataset: KineticsBenchmarkDataset,
    config: KineticsBenchmarkConfig = KineticsBenchmarkConfig(),
) -> Tuple[KineticsPrediction, ...]:
    """Fit the HydroSheaf synthetic specialist without accessing hidden truth."""

    predictions: List[KineticsPrediction] = []
    for observation in dataset.observations:
        regime = next(regime for regime in config.regimes if regime.name == observation.regime)
        velocity, delay, front_residual = _estimate_front(observation, regime.distance_m)
        effective_rate, effective_interval, concentration_residual, n_concentration = _estimate_effective_rate(
            observation, delay, config
        )
        area = observation.surface_area_measurement
        area_interval = None
        if area is not None:
            z = _normal_quantile(config.confidence_level)
            half_width = z * config.surface_area_measurement_sd
            area_interval = (max(1.0e-6, area - half_width), area + half_width)
        predictions.append(
            _prediction(
                observation,
                config,
                "hydrosheaf_kinetics",
                effective_rate,
                effective_interval,
                velocity,
                delay,
                concentration_residual,
                front_residual,
                area,
                area_interval,
                area is not None,
                {"fit_method": "front_profile_plus_effective_rate_grid", "n_concentration": n_concentration},
            )
        )
    return tuple(predictions)


def fit_effective_rate_comparator(
    dataset: KineticsBenchmarkDataset,
    config: KineticsBenchmarkConfig = KineticsBenchmarkConfig(),
) -> Tuple[KineticsPrediction, ...]:
    """Fit a fixed effective-rate comparator with matched observations.

    The comparator uses the same front profile and effective-rate grid as the
    specialist, but has no separate ``k`` or surface-area channel.  This makes
    the comparison competence-matched for effective-rate prediction while
    preventing an invalid separate-parameter claim.
    """

    predictions: List[KineticsPrediction] = []
    for observation in dataset.observations:
        regime = next(regime for regime in config.regimes if regime.name == observation.regime)
        velocity, delay, front_residual = _estimate_front(observation, regime.distance_m)
        effective_rate, effective_interval, concentration_residual, n_concentration = _estimate_effective_rate(
            observation, delay, config
        )
        predictions.append(
            _prediction(
                observation,
                config,
                "fixed_effective_rate_comparator",
                effective_rate,
                effective_interval,
                velocity,
                delay,
                concentration_residual,
                front_residual,
                None,
                None,
                False,
                {
                    "fit_method": "matched_effective_rate_profile_grid",
                    "competence_matched": True,
                    "n_concentration": n_concentration,
                    "separate_k_a": False,
                },
            )
        )
    return tuple(predictions)


def _conformal_scale(scores: Sequence[float], confidence_level: float) -> float:
    """Return a finite-sample split-conformal multiplicative scale."""

    finite_scores = np.asarray([score for score in scores if np.isfinite(score)], dtype=float)
    if finite_scores.size == 0:
        return 1.0
    ordered = np.sort(np.maximum(finite_scores, 0.0))
    rank = min(
        len(ordered) - 1,
        max(0, int(math.ceil((len(ordered) + 1) * confidence_level)) - 1),
    )
    return max(1.0, float(ordered[rank]))


def fit_development_interval_calibrator(
    dataset: KineticsBenchmarkDataset,
    predictions: Sequence[KineticsPrediction],
    config: KineticsBenchmarkConfig = KineticsBenchmarkConfig(),
) -> KineticsIntervalCalibrator:
    """Fit interval scales on development truth only.

    The locked cases are never inspected here.  Calibration scores are
    nonconformity ratios relative to each raw interval's half-width, so the
    resulting scales can be applied to predictions from either specialist or
    comparator before locked-test scoring.
    """

    prediction_by_case = {prediction.case_id: prediction for prediction in predictions}
    truth_by_case = dataset._truth_by_case()
    expected = set(truth_by_case)
    if set(prediction_by_case) != expected:
        raise ValueError("predictions must contain exactly one record for every benchmark case")
    if not dataset.development_case_ids or not dataset.locked_case_ids:
        raise ValueError("development calibration requires both development and locked cases")

    effective_scores: List[float] = []
    front_scores: List[float] = []
    concentration_scores: List[float] = []
    for case_id in dataset.development_case_ids:
        prediction = prediction_by_case[case_id]
        truth = truth_by_case[case_id]
        if prediction.effective_rate_per_day is not None and prediction.effective_rate_interval is not None:
            center = prediction.effective_rate_per_day
            lower = max(prediction.effective_rate_interval[0], 1.0e-12)
            upper = max(prediction.effective_rate_interval[1], center)
            half_width = max(math.log(center / lower), math.log(upper / center), 1.0e-12)
            effective_scores.append(abs(math.log(truth.effective_rate_per_day / center)) / half_width)
        front_radius = np.maximum(
            (prediction.front_interval[1] - prediction.front_interval[0]) / 2.0,
            1.0e-12,
        )
        front_scores.append(float(np.max(np.abs(truth.true_front_m - prediction.predicted_front_m) / front_radius)))
        concentration_radius = np.maximum(
            (prediction.concentration_interval[1] - prediction.concentration_interval[0]) / 2.0,
            1.0e-12,
        )
        concentration_scores.append(
            float(
                np.max(
                    np.abs(truth.true_concentration - prediction.predicted_concentration)
                    / concentration_radius
                )
            )
        )
    return KineticsIntervalCalibrator(
        fit_scope="development_only",
        confidence_level=config.confidence_level,
        development_case_ids=dataset.development_case_ids,
        locked_case_ids=dataset.locked_case_ids,
        effective_rate_scale=_conformal_scale(effective_scores, config.confidence_level),
        front_scale=_conformal_scale(front_scores, config.confidence_level),
        concentration_scale=_conformal_scale(concentration_scores, config.confidence_level),
        max_relative_interval_width=config.max_relative_interval_width,
        calibration_status="development_truth_only_locked_apply",
    )


def _log_error(estimate: Optional[float], truth: float) -> Optional[float]:
    if estimate is None or estimate <= 0 or truth <= 0:
        return None
    return float(math.log(estimate / truth))


def _coverage(interval: Optional[Tuple[float, float]], truth: float) -> Optional[bool]:
    if interval is None:
        return None
    return bool(interval[0] <= truth <= interval[1])


def score_kinetics_predictions(
    dataset: KineticsBenchmarkDataset,
    predictions: Sequence[KineticsPrediction],
    config: KineticsBenchmarkConfig = KineticsBenchmarkConfig(),
    case_ids: Optional[Sequence[str]] = None,
) -> KineticsBenchmarkScore:
    """Score predictions against private truth after the truth-blind fit.

    ``case_ids`` is used by the benchmark runner to report locked-test metrics;
    development truth is never pooled into that final score.
    """

    truth_by_case = dataset._truth_by_case()
    prediction_by_case = {prediction.case_id: prediction for prediction in predictions}
    if set(prediction_by_case) != set(truth_by_case):
        raise ValueError("predictions must contain exactly one record for every benchmark case")
    selected_case_ids = tuple(case_ids) if case_ids is not None else tuple(truth_by_case)
    if not selected_case_ids or not set(selected_case_ids).issubset(truth_by_case):
        raise ValueError("case_ids must select at least one known benchmark case")
    selected_predictions = [prediction_by_case[case_id] for case_id in selected_case_ids]
    all_front_errors: List[float] = []
    all_concentration_errors: List[float] = []
    all_predictive_errors: List[float] = []
    ka_errors: List[float] = []
    k_errors: List[float] = []
    area_errors: List[float] = []
    ka_coverages: List[bool] = []
    front_coverages: List[bool] = []
    concentration_coverages: List[bool] = []
    selective_errors: List[float] = []
    false_commitments = 0
    identified = 0
    counts: Dict[str, int] = {}

    for observation in dataset.observations:
        if observation.case_id not in selected_case_ids:
            continue
        truth = truth_by_case[observation.case_id]
        prediction = prediction_by_case[observation.case_id]
        counts[prediction.identifiability] = counts.get(prediction.identifiability, 0) + 1
        front_error = prediction.predicted_front_m - truth.true_front_m
        concentration_error = prediction.predicted_concentration - truth.true_concentration
        all_front_errors.extend(front_error.tolist())
        all_concentration_errors.extend(concentration_error.tolist())
        all_predictive_errors.extend(np.concatenate([front_error, concentration_error]).tolist())
        ka_error = _log_error(prediction.effective_rate_per_day, truth.effective_rate_per_day)
        if ka_error is not None:
            ka_errors.append(ka_error)
            if not prediction.abstain:
                selective_errors.append(ka_error)
        coverage = _coverage(prediction.effective_rate_interval, truth.effective_rate_per_day)
        if coverage is not None:
            ka_coverages.append(coverage)
        front_coverages.extend(
            bool(low <= truth_value <= high)
            for low, high, truth_value in zip(
                prediction.front_interval[0], prediction.front_interval[1], truth.true_front_m
            )
        )
        concentration_coverages.extend(
            bool(low <= truth_value <= high)
            for low, high, truth_value in zip(
                prediction.concentration_interval[0], prediction.concentration_interval[1], truth.true_concentration
            )
        )
        if prediction.k_a_identifiable:
            identified += 1
            k_error = _log_error(prediction.k_per_day, truth.k_per_day)
            area_error = _log_error(prediction.surface_area, truth.surface_area)
            if k_error is not None:
                k_errors.append(k_error)
            if area_error is not None:
                area_errors.append(area_error)
        elif not prediction.parameter_abstain:
            false_commitments += 1

    n_cases = len(selected_case_ids)
    n_predictive = sum(not prediction_by_case[case_id].abstain for case_id in selected_case_ids)
    status = "SCORED" if all_predictive_errors else "ABSTAIN_NO_FINITE_PREDICTIONS"
    return KineticsBenchmarkScore(
        model=predictions[0].model if predictions else "unknown",
        n_cases=n_cases,
        n_predictive_cases=n_predictive,
        abstention_rate=1.0 - n_predictive / max(n_cases, 1),
        parameter_abstention_rate=sum(prediction.parameter_abstain for prediction in selected_predictions) / max(n_cases, 1),
        predictive_rmse=float(np.sqrt(np.mean(np.square(all_predictive_errors)))) if all_predictive_errors else float("nan"),
        front_rmse=float(np.sqrt(np.mean(np.square(all_front_errors)))) if all_front_errors else float("nan"),
        concentration_rmse=float(np.sqrt(np.mean(np.square(all_concentration_errors)))) if all_concentration_errors else float("nan"),
        effective_rate_log_rmse=float(np.sqrt(np.mean(np.square(ka_errors)))) if ka_errors else float("nan"),
        k_log_rmse_identified=float(np.sqrt(np.mean(np.square(k_errors)))) if k_errors else None,
        surface_area_log_rmse_identified=float(np.sqrt(np.mean(np.square(area_errors)))) if area_errors else None,
        effective_rate_interval_coverage=float(np.mean(ka_coverages)) if ka_coverages else float("nan"),
        front_interval_coverage=float(np.mean(front_coverages)) if front_coverages else float("nan"),
        concentration_interval_coverage=float(np.mean(concentration_coverages)) if concentration_coverages else float("nan"),
        selective_effective_rate_log_rmse=float(np.sqrt(np.mean(np.square(selective_errors)))) if selective_errors else None,
        false_commitment_rate=false_commitments / max(n_cases, 1),
        identified_case_rate=identified / max(n_cases, 1),
        identifiability_counts=counts,
        status=status,
    )


def run_kinetics_specialist_benchmark(
    config: KineticsBenchmarkConfig = KineticsBenchmarkConfig(),
) -> KineticsBenchmarkReport:
    """Run the complete deterministic M8 benchmark without releasing truth."""

    dataset = generate_kinetics_dataset(config)
    raw_specialist_predictions = fit_kinetics_specialist(dataset, config)
    raw_comparator_predictions = fit_effective_rate_comparator(dataset, config)
    interval_calibrator = fit_development_interval_calibrator(
        dataset, raw_specialist_predictions, config
    )
    comparator_interval_calibrator = fit_development_interval_calibrator(
        dataset, raw_comparator_predictions, config
    )
    specialist_predictions = interval_calibrator.apply(raw_specialist_predictions)
    comparator_predictions = comparator_interval_calibrator.apply(raw_comparator_predictions)
    specialist_score = score_kinetics_predictions(
        dataset, specialist_predictions, config, case_ids=dataset.locked_case_ids
    )
    comparator_score = score_kinetics_predictions(
        dataset, comparator_predictions, config, case_ids=dataset.locked_case_ids
    )
    return KineticsBenchmarkReport(
        benchmark_version=config.benchmark_version,
        config=config.as_dict(),
        config_hash=config.config_hash,
        dataset=dataset,
        specialist_predictions=specialist_predictions,
        comparator_predictions=comparator_predictions,
        interval_calibrator=interval_calibrator,
        comparator_interval_calibrator=comparator_interval_calibrator,
        specialist_score=specialist_score,
        comparator_score=comparator_score,
        comparator_status=COMPETENCE_MATCHED_COMPARATOR_STATUS,
        claim_status="ABSTAIN_COMPONENT_ONLY",
        claim_boundary=(
            "This is a deterministic controlled-synthetic M8 kinetics diagnostic. "
            "It can report effective-rate prediction, parameter recovery, interval "
            "coverage and selective risk under the declared generators. It does not "
            "establish PHREEQC equivalence, field validity, or general superiority; "
            "separate k and surface area are claimable only when an independent "
            "surface-area measurement is observed."
        ),
    )


__all__ = [
    "BENCHMARK_VERSION",
    "COMPARATOR_STATUS",
    "COMPETENCE_MATCHED_COMPARATOR_STATUS",
    "DEFAULT_REGIMES",
    "KineticRegime",
    "KineticsBenchmarkConfig",
    "KineticsObservation",
    "KineticsBenchmarkDataset",
    "KineticsPrediction",
    "KineticsBenchmarkScore",
    "KineticsBenchmarkReport",
    "KineticsIntervalCalibrator",
    "generate_kinetics_dataset",
    "fit_kinetics_specialist",
    "fit_effective_rate_comparator",
    "fit_development_interval_calibrator",
    "score_kinetics_predictions",
    "run_kinetics_specialist_benchmark",
]
