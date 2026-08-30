"""Independent competence-matched groundwater-age specialist baseline.

This module is intentionally self-contained.  It imports no HydroSheaf
inference, posterior, candidate-generation, or calibration code.  It uses
only the declared tracer-age observation channel and a frozen set of compact
public screening histories plus radioactive-decay constants.

The baseline evaluates a deterministic candidate universe of bounded
exponential and gamma transit-time distributions.  Candidate likelihoods are
combined with a simple Gaussian measurement model, and the resulting mixture
of age distributions supplies a point estimate and a model-based interval.
The interval is not claimed to be calibrated.  Unsupported observations,
large tracer disagreement, multimodality, or a weakly identified interval
produce an explicit ``abstain`` decision.

The public integration point is ``AgeCompetentBaseline.predict`` (or the
``predict_age_baseline`` convenience function).  It accepts the same blind
shape used by the programme runner::

    {"tracer_age_history": {"nodes": {"site": {
        "tritium_TU": 2.5,
        "argon39_pmc": 95.0,
        "sample_date": 2025.5,
    }}}}

The returned records use ``mean_age_years``, ``age_95_low``, and
``age_95_high`` so they can be adapted to the existing specialist scoring
contract without changing this module or the runner.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
from itertools import product
import json
import math

import numpy as np


ABSTAIN = "abstain"
SELECT = "select"
AGE_OUTPUT = "age_interval"

_BASELINE_NAME = "independent_competence_matched_age_specialist"
_BASELINE_VERSION = "1.2"
_EVIDENCE_CHANNEL = "tracer_age_history"

# These compact histories are the fixed screening curves used by the public
# multi-tracer age examples.  They are deliberately embedded as immutable
# tuples so this comparator does not call a HydroSheaf history builder or
# inherit a site-specific/posterior state.
FIXED_PUBLIC_ATMOSPHERIC_HISTORIES: Mapping[
    str, tuple[tuple[float, ...], tuple[float, ...]]
] = {
    "cfc11_pptv": (
        (1940.0, 1950.0, 1960.0, 1970.0, 1980.0, 1990.0, 2000.0, 2010.0, 2020.0, 2025.0),
        (0.0, 2.0, 25.0, 135.0, 225.0, 260.0, 256.0, 240.0, 226.0, 220.0),
    ),
    "cfc12_pptv": (
        (1940.0, 1950.0, 1960.0, 1970.0, 1980.0, 1990.0, 2000.0, 2010.0, 2020.0, 2025.0),
        (0.0, 8.0, 125.0, 285.0, 430.0, 525.0, 540.0, 528.0, 505.0, 496.0),
    ),
    "cfc113_pptv": (
        (1955.0, 1965.0, 1975.0, 1985.0, 1995.0, 2005.0, 2015.0, 2025.0),
        (0.0, 1.0, 12.0, 55.0, 85.0, 82.0, 74.0, 68.0),
    ),
    "sf6_pptv": (
        (1960.0, 1970.0, 1980.0, 1990.0, 2000.0, 2010.0, 2020.0, 2025.0),
        (0.0, 0.03, 0.55, 2.1, 4.2, 7.0, 10.2, 12.0),
    ),
}

_DECAY_TRACERS: Mapping[str, tuple[float, float]] = {
    # field -> (modern reference, half-life in years)
    "tritium_TU": (8.0, 12.32),
    "argon39_pmc": (100.0, 269.0),
    "c14_pmc": (100.0, 5730.0),
}
_GAS_TRACERS = frozenset(FIXED_PUBLIC_ATMOSPHERIC_HISTORIES)
_GAS_PREFIXES = {
    "cfc11_pptv": "cfc11",
    "cfc12_pptv": "cfc12",
    "cfc113_pptv": "cfc113",
    "sf6_pptv": "sf6",
}
_SUPPORTED_TRACERS = (
    "tritium_TU",
    "argon39_pmc",
    "c14_pmc",
    "cfc11_pptv",
    "cfc12_pptv",
    "cfc113_pptv",
    "sf6_pptv",
)
_DEFAULT_SIGMA = {
    "tritium_TU": 0.15,
    "argon39_pmc": 2.0,
    "c14_pmc": 0.8,
    "cfc11_pptv": 0.25,
    "cfc12_pptv": 0.25,
    "cfc113_pptv": 0.15,
    "sf6_pptv": 0.05,
}
_SIGMA_FIELDS = {
    "tritium_TU": "tritium_sigma_TU",
    "argon39_pmc": "argon39_sigma_pmc",
    "c14_pmc": "c14_sigma_pmc",
    "cfc11_pptv": "cfc11_sigma_pptv",
    "cfc12_pptv": "cfc12_sigma_pptv",
    "cfc113_pptv": "cfc113_sigma_pptv",
    "sf6_pptv": "sf6_sigma_pptv",
}
_FORBIDDEN_TRUTH_KEYS = frozenset(
    {
        "age",
        "age_years",
        "age_predictions",
        "candidate_edges",
        "edge",
        "pathline_rows",
        "posterior",
        "posterior_age",
        "process",
        "reference_age",
        "travel_age_years",
        "true_age",
        "true_age_years",
        "true_ages_years",
        "true_edges",
        "truth",
        "truth_age",
    }
)


def _finite(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _fingerprint(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _grid(maximum: float, step: float) -> tuple[float, ...]:
    values = [float(index * step) for index in range(int(math.floor(maximum / step)) + 1)]
    if not values or values[-1] < maximum - 1.0e-10:
        values.append(float(maximum))
    else:
        values[-1] = float(maximum)
    return tuple(round(value, 10) for value in values)


@dataclass(frozen=True)
class AgeBaselineConfig:
    """Frozen controls for the independent candidate universe and abstention."""

    max_age_years: float = 120.0
    candidate_step_years: float = 2.0
    quadrature_step_years: float = 0.5
    gamma_shapes: tuple[float, ...] = (2.0, 4.0)
    min_supported_tracers: int = 2
    ambiguity_age_gap_years: float = 40.0
    max_interval_width_fraction: float = 0.85
    min_peak_separation_years: float = 8.0
    min_secondary_peak_fraction: float = 0.25
    min_tracer_log_likelihood_range: float = 2.0
    excess_air_grid_fraction: tuple[float, ...] = (0.0, 0.005, 0.01, 0.02)
    degradation_grid_fraction: tuple[float, ...] = (0.0, 0.25, 0.50, 0.75)

    def __post_init__(self) -> None:
        maximum = _finite(self.max_age_years)
        candidate_step = _finite(self.candidate_step_years)
        quadrature_step = _finite(self.quadrature_step_years)
        ambiguity_gap = _finite(self.ambiguity_age_gap_years)
        width_fraction = _finite(self.max_interval_width_fraction)
        peak_separation = _finite(self.min_peak_separation_years)
        peak_fraction = _finite(self.min_secondary_peak_fraction)
        tracer_information = _finite(self.min_tracer_log_likelihood_range)
        if maximum is None or maximum <= 0.0:
            raise ValueError("max_age_years must be finite and positive.")
        if candidate_step is None or candidate_step <= 0.0 or candidate_step > maximum:
            raise ValueError("candidate_step_years must lie in (0, max_age_years].")
        if quadrature_step is None or quadrature_step <= 0.0:
            raise ValueError("quadrature_step_years must be finite and positive.")
        if ambiguity_gap is None or ambiguity_gap < 0.0:
            raise ValueError("ambiguity_age_gap_years must be finite and non-negative.")
        if width_fraction is None or not 0.0 < width_fraction <= 1.0:
            raise ValueError("max_interval_width_fraction must lie in (0, 1].")
        if peak_separation is None or peak_separation < 0.0:
            raise ValueError("min_peak_separation_years must be non-negative.")
        if peak_fraction is None or not 0.0 < peak_fraction < 1.0:
            raise ValueError("min_secondary_peak_fraction must lie in (0, 1).")
        if tracer_information is None or tracer_information < 0.0:
            raise ValueError("min_tracer_log_likelihood_range must be finite and non-negative.")
        if int(self.min_supported_tracers) < 1:
            raise ValueError("min_supported_tracers must be at least one.")

        excess_air_grid = tuple(float(value) for value in self.excess_air_grid_fraction)
        degradation_grid = tuple(float(value) for value in self.degradation_grid_fraction)
        if (
            not excess_air_grid
            or any(not math.isfinite(value) or value < 0.0 for value in excess_air_grid)
            or len(set(excess_air_grid)) != len(excess_air_grid)
        ):
            raise ValueError("excess_air_grid_fraction must contain unique finite non-negative values.")
        if (
            not degradation_grid
            or any(not math.isfinite(value) or not 0.0 <= value < 1.0 for value in degradation_grid)
            or len(set(degradation_grid)) != len(degradation_grid)
        ):
            raise ValueError("degradation_grid_fraction must contain unique values in [0, 1).")

        shapes = tuple(float(shape) for shape in self.gamma_shapes)
        if any(not math.isfinite(shape) or shape <= 0.0 or math.isclose(shape, 1.0) for shape in shapes):
            raise ValueError("gamma_shapes must be finite positive values other than one.")
        if len(set(shapes)) != len(shapes):
            raise ValueError("gamma_shapes must be unique.")
        object.__setattr__(self, "max_age_years", float(maximum))
        object.__setattr__(self, "candidate_step_years", float(candidate_step))
        object.__setattr__(self, "quadrature_step_years", float(quadrature_step))
        object.__setattr__(self, "gamma_shapes", shapes)
        object.__setattr__(self, "min_supported_tracers", int(self.min_supported_tracers))
        object.__setattr__(self, "ambiguity_age_gap_years", float(ambiguity_gap))
        object.__setattr__(self, "max_interval_width_fraction", float(width_fraction))
        object.__setattr__(self, "min_peak_separation_years", float(peak_separation))
        object.__setattr__(self, "min_secondary_peak_fraction", float(peak_fraction))
        object.__setattr__(self, "min_tracer_log_likelihood_range", float(tracer_information))
        object.__setattr__(self, "excess_air_grid_fraction", tuple(sorted(excess_air_grid)))
        object.__setattr__(self, "degradation_grid_fraction", tuple(sorted(degradation_grid)))

    def to_dict(self) -> dict[str, object]:
        return {
            "max_age_years": self.max_age_years,
            "candidate_step_years": self.candidate_step_years,
            "quadrature_step_years": self.quadrature_step_years,
            "gamma_shapes": list(self.gamma_shapes),
            "min_supported_tracers": self.min_supported_tracers,
            "ambiguity_age_gap_years": self.ambiguity_age_gap_years,
            "max_interval_width_fraction": self.max_interval_width_fraction,
            "min_peak_separation_years": self.min_peak_separation_years,
            "min_secondary_peak_fraction": self.min_secondary_peak_fraction,
            "min_tracer_log_likelihood_range": self.min_tracer_log_likelihood_range,
            "excess_air_grid_fraction": list(self.excess_air_grid_fraction),
            "degradation_grid_fraction": list(self.degradation_grid_fraction),
        }


@dataclass(frozen=True)
class AgeCandidate:
    """One bounded exponential/gamma transit-time candidate."""

    candidate_id: str
    family: str
    mean_age_years: float
    shape: float
    scale_years: float

    def to_dict(self) -> dict[str, object]:
        return {
            "candidate_id": self.candidate_id,
            "family": self.family,
            "mean_age_years": self.mean_age_years,
            "shape": self.shape,
            "scale_years": self.scale_years,
        }


@dataclass(frozen=True)
class AgeCandidateUniverse:
    """Deterministic, truth-blind candidate set used by the baseline."""

    age_grid_years: tuple[float, ...]
    candidates: tuple[AgeCandidate, ...]
    config: AgeBaselineConfig

    @property
    def candidate_hash(self) -> str:
        return _fingerprint(
            {
                "age_grid_years": list(self.age_grid_years),
                "candidates": [candidate.to_dict() for candidate in self.candidates],
                "config": self.config.to_dict(),
            }
        )

    def to_audit_record(self) -> dict[str, object]:
        return {
            "algorithm": "bounded_exponential_gamma_candidate_grid_v1",
            "candidate_count": len(self.candidates),
            "age_grid_count": len(self.age_grid_years),
            "age_min_years": self.age_grid_years[0],
            "age_max_years": self.age_grid_years[-1],
            "families": sorted({candidate.family for candidate in self.candidates}),
            "candidate_hash": self.candidate_hash,
            "truth_blind": True,
            "config": self.config.to_dict(),
        }


def generate_age_candidates(config: AgeBaselineConfig | None = None) -> AgeCandidateUniverse:
    """Enumerate the independent bounded exponential/gamma candidate universe."""

    cfg = config or AgeBaselineConfig()
    means = _grid(cfg.max_age_years, cfg.candidate_step_years)
    families: tuple[tuple[str, float], ...] = (
        ("exponential", 1.0),
        *(("gamma", shape) for shape in cfg.gamma_shapes),
    )
    candidates: list[AgeCandidate] = []
    for family, shape in families:
        for mean in means:
            scale = float(mean / shape) if mean > 0.0 else 0.0
            candidate_id = f"{family}:shape={shape:g}:mean={mean:g}"
            candidates.append(
                AgeCandidate(
                    candidate_id=candidate_id,
                    family=family,
                    mean_age_years=float(mean),
                    shape=float(shape),
                    scale_years=scale,
                )
            )
    age_grid = _grid(cfg.max_age_years, cfg.quadrature_step_years)
    return AgeCandidateUniverse(age_grid, tuple(candidates), cfg)


def enumerate_age_candidates(config: AgeBaselineConfig | None = None) -> AgeCandidateUniverse:
    """Alias with an explicit enumeration name for audit and test callers."""

    return generate_age_candidates(config)


def _assert_truth_blind(value: object, path: tuple[str, ...] = ()) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            key_text = str(key)
            lowered = key_text.lower()
            if (
                lowered in _FORBIDDEN_TRUTH_KEYS
                or lowered.startswith("true_")
                or lowered.startswith("truth_")
                or "posterior" in lowered
            ):
                dotted = ".".join((*path, key_text))
                raise ValueError(f"Truth/reference field is forbidden: {dotted}")
            _assert_truth_blind(item, (*path, key_text))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _assert_truth_blind(item, (*path, str(index)))


def _extract_nodes(observations: Mapping[str, object]) -> dict[str, Mapping[str, object]]:
    _assert_truth_blind(observations)
    channel: object = observations
    if "tracer_age_history" in observations:
        channel = observations["tracer_age_history"]
    elif "tracer_age" in observations:
        channel = observations["tracer_age"]
    if isinstance(channel, Mapping) and "nodes" in channel:
        channel = channel["nodes"]
    if not isinstance(channel, Mapping):
        raise TypeError("Age observations must contain a mapping of nodes.")

    nodes: dict[str, Mapping[str, object]] = {}
    for raw_target, raw_features in channel.items():
        target = str(raw_target).strip()
        if not target:
            raise ValueError("Age observation targets must be non-empty.")
        if target in nodes:
            raise ValueError(f"Duplicate age observation target: {target!r}")
        if not isinstance(raw_features, Mapping):
            raise TypeError(f"Age features for {target!r} must be a mapping.")
        nodes[target] = raw_features
    return dict(sorted(nodes.items()))


@dataclass(frozen=True)
class _TracerObservation:
    field: str
    value: float
    sigma: float
    sample_year: float | None
    measurement_kind: str = "atmospheric_equivalent"
    solubility_per_pptv: float = 1.0


def _collect_tracers(features: Mapping[str, object]) -> tuple[tuple[_TracerObservation, ...], tuple[str, ...]]:
    usable: list[_TracerObservation] = []
    unsupported: list[str] = []
    sample_year = _finite(features.get("sample_year", features.get("sample_date")))
    for field in _SUPPORTED_TRACERS:
        if field not in features:
            continue
        value = _finite(features.get(field))
        if value is None or value < 0.0:
            unsupported.append(field)
            continue
        if field in _GAS_TRACERS and sample_year is None:
            unsupported.append(f"{field}:missing_sample_year")
            continue
        sigma_field = _SIGMA_FIELDS[field]
        supplied_sigma = features.get(sigma_field)
        if supplied_sigma is None:
            sigma = float(_DEFAULT_SIGMA[field])
        else:
            sigma_value = _finite(supplied_sigma)
            if sigma_value is None or sigma_value <= 0.0:
                unsupported.append(f"{field}:invalid_sigma")
                continue
            sigma = float(sigma_value)
        sigma = max(sigma, 0.08 * max(abs(value), _DEFAULT_SIGMA[field]))
        usable.append(_TracerObservation(field, float(value), sigma, sample_year))
    # When a dissolved-gas observation and an explicit solubility are
    # available, retain the measurement in its native units.  This permits a
    # truth-blind nuisance profile over excess air instead of silently
    # treating the dissolved concentration as an atmospheric equivalent.
    for field, prefix in _GAS_PREFIXES.items():
        if field in features:
            continue
        dissolved = _finite(features.get(f"{prefix}_dissolved"))
        solubility = _finite(features.get(f"{prefix}_solubility_per_pptv"))
        if dissolved is None or solubility is None:
            if dissolved is not None or solubility is not None:
                unsupported.append(f"{prefix}:missing_or_invalid_solubility_schema")
            continue
        if dissolved < 0.0 or solubility <= 0.0 or sample_year is None:
            unsupported.append(
                f"{prefix}:invalid_dissolved_gas_schema"
                if dissolved < 0.0 or solubility <= 0.0
                else f"{prefix}_dissolved:missing_sample_year"
            )
            continue
        sigma_raw = features.get(f"{prefix}_dissolved_sigma")
        sigma_value = _finite(sigma_raw) if sigma_raw is not None else None
        if sigma_raw is not None and (sigma_value is None or sigma_value <= 0.0):
            unsupported.append(f"{prefix}_dissolved:invalid_sigma")
            continue
        sigma = float(sigma_value if sigma_value is not None else _DEFAULT_SIGMA[field] * solubility)
        sigma = max(sigma, 0.08 * max(abs(dissolved), _DEFAULT_SIGMA[field] * solubility))
        usable.append(
            _TracerObservation(
                field,
                float(dissolved),
                sigma,
                sample_year,
                measurement_kind="dissolved_gas",
                solubility_per_pptv=float(solubility),
            )
        )
    return tuple(usable), tuple(sorted(set(unsupported)))


def _optional_fraction(
    features: Mapping[str, object],
    prefix: str,
    suffix: str,
) -> float | None:
    """Read an explicitly declared nuisance fraction, never infer one."""

    for key in (f"{prefix}_{suffix}", suffix):
        if key not in features:
            continue
        value = _finite(features.get(key))
        if value is None:
            return None
        return value
    return None


def _correction_scenarios(
    features: Mapping[str, object],
    tracers: tuple[_TracerObservation, ...],
    config: AgeBaselineConfig,
) -> tuple[dict[str, object], ...]:
    """Build a small, declared nuisance universe for gas corrections.

    No correction is profiled unless the observation schema contains either a
    dissolved-gas/solubility pair or an explicit correction fraction.  An
    explicit fraction is treated as one alternative alongside the no-
    correction control.  This makes the correction assumption visible rather
    than silently choosing it from the held-out truth.
    """

    gas_tracers = [tracer for tracer in tracers if tracer.field in _GAS_PREFIXES]
    if not gas_tracers:
        return (
            {
                "excess_air_fraction": 0.0,
                "degradation_fraction": 0.0,
                "history_model": "public_atmospheric",
            },
        )
    explicit_nuisance = any(
        "excess_air_fraction" in str(key)
        or "degradation_fraction" in str(key)
        or str(key).endswith("_dissolved")
        or str(key).endswith("_solubility_per_pptv")
        for key in features
    )
    # Direct synthetic gas observations are allowed to come from a different
    # recharge-history family.  Profile the alternatives only when the input
    # does not already declare a site-specific excess-air/degradation schema;
    # declared nuisance tests retain their exact two-/four-scenario contract.
    history_models = (
        ("public_atmospheric",)
        if explicit_nuisance
        else (
            "public_atmospheric",
            "exponential_surrogate",
            "two_component_mixture",
            "logistic_recharge",
        )
    )
    excess_options = {0.0}
    common_degradation_options = {0.0}
    common_degradation_declared = "degradation_fraction" in features
    per_tracer_degradation: dict[str, tuple[float, ...]] = {}
    for tracer in gas_tracers:
        prefix = _GAS_PREFIXES[tracer.field]
        supplied = _optional_fraction(features, prefix, "excess_air_fraction")
        if tracer.measurement_kind == "dissolved_gas":
            excess_options.update(
                config.excess_air_grid_fraction if supplied is None else (max(0.0, supplied),)
            )
        elif supplied is not None:
            # A direct gas field may still be accompanied by an explicit
            # correction declaration from a preprocessing step.  Preserve a
            # no-correction control and the declared alternative; do not infer
            # a correction from the gas concentration itself.
            excess_options.add(max(0.0, supplied))
        supplied_degradation = _optional_fraction(features, prefix, "degradation_fraction")
        if supplied_degradation is not None:
            if not 0.0 <= supplied_degradation < 1.0:
                continue
            if f"{prefix}_degradation_fraction" in features:
                per_tracer_degradation[tracer.field] = (0.0, float(supplied_degradation))
            else:
                common_degradation_options.add(float(supplied_degradation))
        elif common_degradation_declared:
            common_degradation_options.update(config.degradation_grid_fraction)
    if common_degradation_declared:
        per_tracer_degradation = {}
    tracer_keys = tuple(sorted(per_tracer_degradation))
    tracer_options = (
        tuple(product(*(per_tracer_degradation[key] for key in tracer_keys)))
        if tracer_keys
        else ((),)
    )
    scenarios: list[dict[str, object]] = []
    for history_model in history_models:
        for excess_air in sorted(excess_options):
            for common_degradation in sorted(common_degradation_options):
                for tracer_values in tracer_options:
                    scenario: dict[str, object] = {
                        "excess_air_fraction": float(excess_air),
                        "degradation_fraction": float(common_degradation),
                        "history_model": history_model,
                    }
                    scenario.update(
                        {
                            f"{key}_degradation_fraction": float(value)
                            for key, value in zip(tracer_keys, tracer_values)
                        }
                    )
                    scenarios.append(scenario)
    return tuple(scenarios)


def _candidate_mass(candidate: AgeCandidate, age_grid: np.ndarray) -> np.ndarray:
    if candidate.mean_age_years <= 0.0:
        mass = np.zeros(age_grid.size, dtype=float)
        mass[0] = 1.0
        return mass
    scale = candidate.scale_years
    if candidate.shape == 1.0:
        log_density = -age_grid / scale - math.log(scale)
    else:
        log_density = np.full(age_grid.shape, -np.inf, dtype=float)
        positive = age_grid > 0.0
        log_density[positive] = (
            (candidate.shape - 1.0) * np.log(age_grid[positive])
            - age_grid[positive] / scale
            - math.lgamma(candidate.shape)
            - candidate.shape * math.log(scale)
        )
    finite = np.isfinite(log_density)
    if not np.any(finite):
        raise RuntimeError(f"Candidate {candidate.candidate_id} has no finite density.")
    mass = np.zeros_like(age_grid, dtype=float)
    maximum = float(np.max(log_density[finite]))
    mass[finite] = np.exp(log_density[finite] - maximum)
    total = float(mass.sum())
    if not math.isfinite(total) or total <= 0.0:
        raise RuntimeError(f"Candidate {candidate.candidate_id} has invalid mass.")
    return mass / total


def _tracer_forward_values(
    tracer: _TracerObservation,
    age_grid: np.ndarray,
    *,
    excess_air_fraction: float = 0.0,
    degradation_fraction: float = 0.0,
    tracer_degradation_fraction: float | None = None,
    history_model: str = "public_atmospheric",
) -> np.ndarray:
    if tracer.field in _DECAY_TRACERS:
        reference, half_life = _DECAY_TRACERS[tracer.field]
        ages = np.asarray(age_grid, dtype=float)
        if history_model == "exponential_surrogate":
            values = {
                "tritium_TU": 0.15 + 6.4 * np.exp(-ages / 17.0) + 0.10 * np.cos(ages / 8.0),
                "argon39_pmc": 97.0 * np.exp(-ages / 330.0),
                "c14_pmc": 96.0 * np.exp(-math.log(2.0) * ages / 5730.0),
            }[tracer.field]
        elif history_model == "two_component_mixture":
            young_fraction = np.clip(0.72 - 0.004 * ages, 0.18, 0.72)
            young = np.exp(-ages / 16.0)
            old = np.exp(-ages / 270.0)
            values = {
                "tritium_TU": 0.12 + 6.0 * (young_fraction * young + 0.025 * old),
                "argon39_pmc": 96.0 * (young_fraction * young + (1.0 - young_fraction) * old),
                "c14_pmc": 98.0 * (0.35 * young + 0.65 * old),
            }[tracer.field]
        else:
            values = reference * np.exp(-math.log(2.0) * ages / half_life)
    else:
        ages = np.asarray(age_grid, dtype=float)
        if history_model == "public_atmospheric":
            years, values_raw = FIXED_PUBLIC_ATMOSPHERIC_HISTORIES[tracer.field]
            values = np.interp(
                float(tracer.sample_year) - ages,
                np.asarray(years, dtype=float),
                np.asarray(values_raw, dtype=float),
            )
        elif history_model == "exponential_surrogate":
            cfc12 = 520.0 * np.exp(-ages / 38.0) / (1.0 + 0.006 * ages)
            values = {
                "cfc11_pptv": 0.72 * cfc12,
                "cfc12_pptv": cfc12,
                "cfc113_pptv": 0.16 * cfc12,
                "sf6_pptv": 9.0 * np.exp(-ages / 42.0),
            }[tracer.field]
        elif history_model == "two_component_mixture":
            young_fraction = np.clip(0.72 - 0.004 * ages, 0.18, 0.72)
            young = np.exp(-ages / 16.0)
            old = np.exp(-ages / 270.0)
            values = {
                "cfc11_pptv": 0.68 * 520.0 * (0.78 * young + 0.22 * old),
                "cfc12_pptv": 520.0 * (0.78 * young + 0.22 * old),
                "cfc113_pptv": 0.15 * 520.0 * (0.78 * young + 0.22 * old),
                "sf6_pptv": 10.0 * (0.86 * young + 0.14 * old),
            }[tracer.field]
        elif history_model == "logistic_recharge":
            recharge_year = float(tracer.sample_year) - ages
            values = {
                "cfc11_pptv": 270.0 / (1.0 + np.exp(-(recharge_year - 1974.0) / 6.8)),
                "cfc12_pptv": 540.0 / (1.0 + np.exp(-(recharge_year - 1976.0) / 8.0)),
                "cfc113_pptv": 82.0 / (1.0 + np.exp(-(recharge_year - 1983.0) / 5.4)),
                "sf6_pptv": 11.0 / (1.0 + np.exp(-(recharge_year - 2004.0) / 10.0)),
            }[tracer.field]
        else:
            raise ValueError(f"Unsupported atmospheric history model: {history_model}")
        effective_degradation = (
            degradation_fraction
            if tracer_degradation_fraction is None
            else tracer_degradation_fraction
        )
        values = values * (1.0 - effective_degradation)
        if tracer.measurement_kind == "dissolved_gas":
            values = values * tracer.solubility_per_pptv * (1.0 + excess_air_fraction)
    return np.asarray(values, dtype=float)


def _expected_tracer(
    tracer: _TracerObservation,
    age_grid: np.ndarray,
    mass: np.ndarray,
    *,
    excess_air_fraction: float = 0.0,
    degradation_fraction: float = 0.0,
    tracer_degradation_fraction: float | None = None,
    history_model: str = "public_atmospheric",
) -> float:
    values = _tracer_forward_values(
        tracer,
        age_grid,
        excess_air_fraction=excess_air_fraction,
        degradation_fraction=degradation_fraction,
        tracer_degradation_fraction=tracer_degradation_fraction,
        history_model=history_model,
    )
    return float(np.dot(mass, values))


def _log_tracer_likelihood(
    tracer: _TracerObservation,
    age_grid: np.ndarray,
    mass: np.ndarray,
    scenario: Mapping[str, object],
) -> float:
    tracer_degradation = scenario.get(f"{tracer.field}_degradation_fraction")
    expected = _expected_tracer(
        tracer,
        age_grid,
        mass,
        excess_air_fraction=float(scenario.get("excess_air_fraction", 0.0)),
        degradation_fraction=float(scenario.get("degradation_fraction", 0.0)),
        tracer_degradation_fraction=(
            None if tracer_degradation is None else float(tracer_degradation)
        ),
        history_model=str(scenario.get("history_model", "public_atmospheric")),
    )
    residual = (tracer.value - expected) / tracer.sigma
    return float(-0.5 * residual * residual - math.log(tracer.sigma))


def _log_marginal_tracer_likelihood(
    tracer: _TracerObservation,
    age_grid: np.ndarray,
    mass: np.ndarray,
    scenarios: tuple[Mapping[str, object], ...],
) -> float:
    values = np.asarray(
        [_log_tracer_likelihood(tracer, age_grid, mass, scenario) for scenario in scenarios],
        dtype=float,
    )
    maximum = float(np.max(values))
    return float(maximum + math.log(float(np.exp(values - maximum).sum())) - math.log(len(values)))


def _log_marginal_joint_likelihood(
    tracers: tuple[_TracerObservation, ...],
    age_grid: np.ndarray,
    mass: np.ndarray,
    scenarios: tuple[Mapping[str, object], ...],
) -> float:
    """Integrate one shared history/nuisance scenario for all tracers.

    Marginalising each tracer independently permits one gas tracer to select
    one atmospheric history while another selects a different history.  That
    is incoherent for a sample-level nuisance model, so the finite-model
    average is taken over a common scenario for all tracers.
    """

    values = np.asarray(
        [
            sum(_log_tracer_likelihood(tracer, age_grid, mass, scenario) for tracer in tracers)
            for scenario in scenarios
        ],
        dtype=float,
    )
    maximum = float(np.max(values))
    return float(maximum + math.log(float(np.exp(values - maximum).sum())) - math.log(len(values)))


def _logsumexp_axis(values: np.ndarray, axis: int) -> np.ndarray:
    """Small NumPy-only log-sum-exp helper for the age likelihood matrix."""

    maximum = np.max(values, axis=axis, keepdims=True)
    result = maximum + np.log(np.exp(values - maximum).sum(axis=axis, keepdims=True))
    return np.squeeze(result, axis=axis)


def _softmax(log_values: np.ndarray) -> np.ndarray:
    if not np.isfinite(log_values).all():
        raise ValueError("Age candidate likelihoods are non-finite.")
    maximum = float(np.max(log_values))
    weights = np.exp(log_values - maximum)
    total = float(weights.sum())
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("Age candidate likelihoods have no positive mass.")
    return weights / total


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    order = np.argsort(values)
    ordered_values = values[order]
    ordered_weights = weights[order]
    cumulative = np.cumsum(ordered_weights)
    index = int(np.searchsorted(cumulative, float(quantile), side="left"))
    return float(ordered_values[min(index, ordered_values.size - 1)])


def _local_peaks(mass: np.ndarray) -> tuple[int, ...]:
    if mass.size < 3:
        return tuple()
    peaks = [
        index
        for index in range(1, mass.size - 1)
        if mass[index] >= mass[index - 1] and mass[index] >= mass[index + 1]
    ]
    return tuple(peaks)


class AgeCompetentBaseline:
    """Independent bounded-distribution age comparator."""

    name = _BASELINE_NAME
    version = _BASELINE_VERSION
    output_kind = AGE_OUTPUT

    def __init__(self, config: AgeBaselineConfig | None = None) -> None:
        self.config = config or AgeBaselineConfig()
        self.universe = generate_age_candidates(self.config)
        age_grid = np.asarray(self.universe.age_grid_years, dtype=float)
        self._age_grid = age_grid
        self._masses = tuple(
            _candidate_mass(candidate, age_grid) for candidate in self.universe.candidates
        )

    def enumerate_candidates(self) -> AgeCandidateUniverse:
        """Return the fixed candidate universe without reading observations."""

        return self.universe

    def to_audit_record(self) -> dict[str, object]:
        return {
            "name": self.name,
            "version": self.version,
            "family": "age_specialist",
            "output_kind": self.output_kind,
            "input_channel": _EVIDENCE_CHANNEL,
            "supported_tracers": list(_SUPPORTED_TRACERS),
            "history_source": "embedded_fixed_public_screening_histories",
            "nuisance_alternatives": {
                "supported_schema": {
                    "dissolved_gas": "<gas>_dissolved plus <gas>_solubility_per_pptv",
                    "excess_air": "<gas>_excess_air_fraction or excess_air_fraction",
                    "degradation": "<gas>_degradation_fraction or degradation_fraction",
                },
                "profile_rule": "declared_schema_only; no correction is inferred when fields are absent",
                "excess_air_grid_fraction": list(self.config.excess_air_grid_fraction),
                "degradation_grid_fraction": list(self.config.degradation_grid_fraction),
            },
            "decay_constants": {
                field: {"reference": reference, "half_life_years": half_life}
                for field, (reference, half_life) in _DECAY_TRACERS.items()
            },
            "candidate_universe": self.universe.to_audit_record(),
            "uncertainty": {
                "type": "candidate_weighted_bounded_age_distribution_quantiles",
                "calibrated": False,
                "interval": "model_based_95_percentile_not_calibrated",
            },
            "abstention": {
                "minimum_supported_tracers": self.config.min_supported_tracers,
                "rules": [
                    "unsupported_or_missing_tracers",
                    "tracer_age_disagreement",
                    "multimodal_age_distribution",
                    "weakly_identified_interval",
                ],
            },
            "truth_blind": True,
            "fingerprint": _fingerprint(
                {
                    "name": self.name,
                    "version": self.version,
                    "candidate_hash": self.universe.candidate_hash,
                    "supported_tracers": list(_SUPPORTED_TRACERS),
                }
            ),
        }

    def _abstention(
        self,
        *,
        reason: str,
        tracers_used: tuple[str, ...] = (),
        unsupported: tuple[str, ...] = (),
    ) -> dict[str, object]:
        maximum = self.config.max_age_years
        return {
            "mean_age_years": maximum / 2.0,
            "age_95_low": 0.0,
            "age_95_high": maximum,
            "uncertainty_years": maximum / 2.0,
            "decision": ABSTAIN,
            "reason": reason,
            "evidence_channel": _EVIDENCE_CHANNEL,
            "tracers_used": list(tracers_used),
            "unsupported_fields": list(unsupported),
            "candidate_count": len(self.universe.candidates),
            "interval_kind": "bounded_abstention_interval_not_an_estimate",
            "calibrated": False,
            "truth_blind": True,
        }

    def _predict_node(self, features: Mapping[str, object]) -> dict[str, object]:
        tracers, unsupported = _collect_tracers(features)
        tracer_names = tuple(sorted(tracer.field for tracer in tracers))
        if len(tracers) < self.config.min_supported_tracers:
            return self._abstention(
                reason="insufficient_supported_age_tracers",
                tracers_used=tracer_names,
                unsupported=unsupported,
            )

        correction_scenarios = _correction_scenarios(features, tracers, self.config)
        mass_matrix = np.asarray(self._masses, dtype=float)
        scenario_count = len(correction_scenarios)
        joint_log_likelihood_matrix = np.zeros(
            (len(self.universe.candidates), scenario_count),
            dtype=float,
        )
        tracer_log_likelihoods: dict[str, np.ndarray] = {}
        tracer_modes: dict[str, float] = {}
        tracer_information: dict[str, float] = {}
        informative_tracer_modes: dict[str, float] = {}
        for tracer in tracers:
            scenario_log_likelihoods = []
            for scenario in correction_scenarios:
                tracer_degradation = scenario.get(f"{tracer.field}_degradation_fraction")
                expected = mass_matrix @ _tracer_forward_values(
                    tracer,
                    self._age_grid,
                    excess_air_fraction=float(scenario.get("excess_air_fraction", 0.0)),
                    degradation_fraction=float(scenario.get("degradation_fraction", 0.0)),
                    tracer_degradation_fraction=(
                        None if tracer_degradation is None else float(tracer_degradation)
                    ),
                    history_model=str(scenario.get("history_model", "public_atmospheric")),
                )
                residual = (tracer.value - expected) / tracer.sigma
                scenario_log_likelihoods.append(
                    -0.5 * residual * residual - math.log(tracer.sigma)
                )
            single_log_likelihood_by_scenario = np.stack(scenario_log_likelihoods, axis=1)
            tracer_log_likelihoods[tracer.field] = single_log_likelihood_by_scenario
            joint_log_likelihood_matrix += single_log_likelihood_by_scenario
            single_log_likelihood = _logsumexp_axis(
                single_log_likelihood_by_scenario,
                axis=1,
            ) - math.log(scenario_count)
            information_range = float(
                np.max(single_log_likelihood) - np.min(single_log_likelihood)
            )
            tracer_information[tracer.field] = information_range
            single_weights = _softmax(single_log_likelihood)
            mode_index = int(np.argmax(single_weights))
            tracer_modes[tracer.field] = self.universe.candidates[mode_index].mean_age_years
            if information_range >= self.config.min_tracer_log_likelihood_range:
                informative_tracer_modes[tracer.field] = tracer_modes[tracer.field]

        candidate_log_likelihood = _logsumexp_axis(
            joint_log_likelihood_matrix,
            axis=1,
        ) - math.log(scenario_count)

        try:
            candidate_weights = _softmax(candidate_log_likelihood)
        except ValueError:
            return self._abstention(
                reason="non_finite_age_likelihood",
                tracers_used=tracer_names,
                unsupported=unsupported,
            )

        age_mass = np.sum(
            np.asarray(candidate_weights, dtype=float)[:, None]
            * np.asarray(self._masses, dtype=float),
            axis=0,
        )
        age_mass = age_mass / float(age_mass.sum())
        mean_age = float(np.dot(self._age_grid, age_mass))
        variance = float(np.dot((self._age_grid - mean_age) ** 2, age_mass))
        low = _weighted_quantile(self._age_grid, age_mass, 0.025)
        high = _weighted_quantile(self._age_grid, age_mass, 0.975)
        interval_width = high - low

        peaks = _local_peaks(age_mass)
        peak_threshold = self.config.min_secondary_peak_fraction * float(age_mass.max())
        substantial_peaks = [index for index in peaks if age_mass[index] >= peak_threshold]
        multimodal = False
        if len(substantial_peaks) >= 2:
            for left_index, left in enumerate(substantial_peaks[:-1]):
                if any(
                    self._age_grid[right] - self._age_grid[left]
                    >= self.config.min_peak_separation_years
                    for right in substantial_peaks[left_index + 1 :]
                ):
                    multimodal = True
                    break

        best_index = int(np.argmax(candidate_weights))
        scenario_weights = np.zeros(len(correction_scenarios), dtype=float)
        joint_log_weights = (
            joint_log_likelihood_matrix - math.log(scenario_count)
        ).reshape(-1)
        if np.isfinite(joint_log_weights).all():
            joint_weights = _softmax(joint_log_weights)
            for index in range(len(self.universe.candidates)):
                start = index * scenario_count
                stop = start + scenario_count
                scenario_weights += candidate_weights[index] * (
                    joint_weights[start:stop] / max(float(joint_weights[start:stop].sum()), 1.0e-300)
                )
        best_scenario_index = int(np.argmax(scenario_weights))
        best_scenario = correction_scenarios[best_scenario_index]
        # Use one coherent scenario for the disagreement diagnostic.  The
        # marginal per-tracer modes above remain in the audit record, but they
        # are not allowed to manufacture disagreement by selecting different
        # atmospheric histories for the same sample.
        coherent_tracer_modes: dict[str, float] = {}
        coherent_tracer_ages: dict[str, float] = {}
        coherent_tracer_information: dict[str, float] = {}
        informative_tracer_modes = {}
        informative_tracer_ages: dict[str, float] = {}
        candidate_mean_ages = np.asarray(
            [candidate.mean_age_years for candidate in self.universe.candidates],
            dtype=float,
        )
        for tracer in tracers:
            coherent_log_likelihood = tracer_log_likelihoods[tracer.field][:, best_scenario_index]
            information_range = float(
                np.max(coherent_log_likelihood) - np.min(coherent_log_likelihood)
            )
            coherent_tracer_information[tracer.field] = information_range
            mode_index = int(np.argmax(coherent_log_likelihood))
            coherent_tracer_modes[tracer.field] = self.universe.candidates[mode_index].mean_age_years
            coherent_weights = _softmax(coherent_log_likelihood)
            coherent_tracer_ages[tracer.field] = float(
                np.dot(coherent_weights, candidate_mean_ages)
            )
            if information_range >= self.config.min_tracer_log_likelihood_range:
                informative_tracer_modes[tracer.field] = coherent_tracer_modes[tracer.field]
                informative_tracer_ages[tracer.field] = coherent_tracer_ages[tracer.field]

        tracer_mode_values = tuple(informative_tracer_ages.values())
        tracer_age_spread = (
            max(tracer_mode_values) - min(tracer_mode_values)
            if tracer_mode_values
            else self.config.max_age_years
        )
        reasons: list[str] = []
        if tracer_age_spread > self.config.ambiguity_age_gap_years:
            reasons.append("tracer_age_disagreement")
        if multimodal:
            reasons.append("multimodal_age_distribution")
        if interval_width > self.config.max_interval_width_fraction * self.config.max_age_years:
            reasons.append("weakly_identified_interval")
        degradation_support_by_tracer = {
            tracer.field: float(
                sum(
                    weight
                    * float(
                        scenario.get(
                            f"{tracer.field}_degradation_fraction",
                            scenario.get("degradation_fraction", 0.0),
                        )
                    )
                    for weight, scenario in zip(scenario_weights, correction_scenarios)
                )
            )
            for tracer in tracers
            if tracer.field in _GAS_PREFIXES
        }
        family_support: dict[str, float] = {}
        for candidate, weight in zip(self.universe.candidates, candidate_weights):
            family_support[candidate.family] = family_support.get(candidate.family, 0.0) + float(weight)
        return {
            "mean_age_years": float(min(max(mean_age, 0.0), self.config.max_age_years)),
            "age_95_low": float(min(max(low, 0.0), self.config.max_age_years)),
            "age_95_high": float(min(max(high, 0.0), self.config.max_age_years)),
            "uncertainty_years": float(math.sqrt(max(variance, 0.0))),
            "decision": ABSTAIN if reasons else SELECT,
            "reason": ";".join(reasons) if reasons else "bounded_exponential_gamma_likelihood",
            "evidence_channel": _EVIDENCE_CHANNEL,
            "tracers_used": list(tracer_names),
            "unsupported_fields": list(unsupported),
            "candidate_count": len(self.universe.candidates),
            "candidate_hash": self.universe.candidate_hash,
            "best_candidate_id": self.universe.candidates[best_index].candidate_id,
            "best_candidate_family": self.universe.candidates[best_index].family,
            "family_support": family_support,
            "correction_model": (
                "profiled_excess_air_and_degradation"
                if any(
                    key in features
                    for key in (
                        "excess_air_fraction",
                        "degradation_fraction",
                        "sf6_excess_air_fraction",
                        "sf6_degradation_fraction",
                    )
                )
                else "profiled_atmospheric_history_and_gas_nuisance"
                if len(correction_scenarios) > 1
                else "none_declared"
            ),
            "correction_scenario_count": len(correction_scenarios),
            "correction_scenarios": [dict(scenario) for scenario in correction_scenarios],
            "correction_support": {
                "excess_air_fraction": float(
                    sum(
                        weight * scenario["excess_air_fraction"]
                        for weight, scenario in zip(scenario_weights, correction_scenarios)
                    )
                ),
                "degradation_fraction": float(
                    sum(
                        weight * scenario["degradation_fraction"]
                        for weight, scenario in zip(scenario_weights, correction_scenarios)
                    )
                ),
                "best_excess_air_fraction": best_scenario["excess_air_fraction"],
                "best_degradation_fraction": best_scenario["degradation_fraction"],
                "degradation_fraction_by_tracer": degradation_support_by_tracer,
                "best_degradation_fraction_by_tracer": {
                    tracer.field: float(
                        best_scenario.get(
                            f"{tracer.field}_degradation_fraction",
                            best_scenario.get("degradation_fraction", 0.0),
                        )
                    )
                    for tracer in tracers
                    if tracer.field in _GAS_PREFIXES
                },
                "history_model_support": {
                    history_model: float(
                        sum(
                            weight
                            for weight, candidate_scenario in zip(
                                scenario_weights, correction_scenarios
                            )
                            if str(
                                candidate_scenario.get(
                                    "history_model", "public_atmospheric"
                                )
                            )
                            == history_model
                        )
                    )
                    for history_model in sorted(
                        {
                            str(
                                candidate_scenario.get(
                                    "history_model", "public_atmospheric"
                                )
                            )
                            for candidate_scenario in correction_scenarios
                        }
                    )
                },
                "best_history_model": str(best_scenario.get("history_model", "public_atmospheric")),
            },
            "tracer_mode_ages_years": tracer_modes,
            "informative_tracer_mode_ages_years": informative_tracer_modes,
            "tracer_information_log_likelihood_range": tracer_information,
            "coherent_tracer_mode_ages_years": coherent_tracer_modes,
            "coherent_tracer_mean_ages_years": coherent_tracer_ages,
            "coherent_tracer_information_log_likelihood_range": coherent_tracer_information,
            "informative_tracers": sorted(informative_tracer_modes),
            "informative_tracer_ages_years": informative_tracer_ages,
            "weak_tracers": sorted(set(tracer_names) - set(informative_tracer_modes)),
            "tracer_age_spread_years": float(tracer_age_spread),
            "multimodal": multimodal,
            "interval_width_years": float(interval_width),
            "interval_kind": "model_based_95_percentile_not_calibrated",
            "calibrated": False,
            "truth_blind": True,
        }

    def predict(self, observations: Mapping[str, object]) -> dict[str, dict[str, object]]:
        """Predict bounded age distributions from declared observations only."""

        nodes = _extract_nodes(observations)
        return {target: self._predict_node(features) for target, features in nodes.items()}


def predict_age_baseline(
    observations: Mapping[str, object],
    *,
    config: AgeBaselineConfig | None = None,
) -> dict[str, dict[str, object]]:
    """Convenience wrapper for the independent age specialist."""

    return AgeCompetentBaseline(config).predict(observations)


__all__ = [
    "ABSTAIN",
    "AGE_OUTPUT",
    "SELECT",
    "FIXED_PUBLIC_ATMOSPHERIC_HISTORIES",
    "AgeBaselineConfig",
    "AgeCandidate",
    "AgeCandidateUniverse",
    "AgeCompetentBaseline",
    "enumerate_age_candidates",
    "generate_age_candidates",
    "predict_age_baseline",
]
