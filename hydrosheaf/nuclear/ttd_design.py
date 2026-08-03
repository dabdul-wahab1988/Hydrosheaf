"""Probability-gated measurement design for transit-time distributions.

This adapter connects groundwater tracer response kernels to HydroSheaf's
generic Bayesian active-learning engine.  An identified set alone does not
define probabilities, so the adapter abstains unless the caller supplies a
probability-bearing TTD ensemble and declares where those probabilities came
from.  In particular, optimizer witnesses are not assigned invented uniform
weights.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Optional, Sequence

import numpy as np

from ..calibration.bayesian_active_learning import (
    AcquisitionConfig,
    MeasurementOption,
    PredictiveScenario,
    rank_measurement_options,
    select_measurement_batch,
)
from .joint_lpm import tracer_response_kernel
from .ttd_identified import AgeFunctional

_PROBABILITY_SEMANTICS = {
    "posterior",
    "declared_prior",
    "calibrated_particle_ensemble",
}


@dataclass(frozen=True)
class TtdHypothesisEnsemble:
    """Discrete TTD hypotheses with optional, provenance-bearing probabilities."""

    hypothesis_ids: Sequence[str]
    age_grid_years: Sequence[float]
    masses: Sequence[Sequence[float]]
    probabilities: Optional[Sequence[float]] = None
    probability_semantics: Optional[str] = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        ids = tuple(str(value) for value in self.hypothesis_ids)
        ages = np.asarray(self.age_grid_years, dtype=float).copy()
        masses = np.asarray(self.masses, dtype=float).copy()
        if len(ids) < 2 or len(set(ids)) != len(ids) or any(not item for item in ids):
            raise ValueError("hypothesis_ids must contain at least two unique names.")
        if ages.ndim != 1 or ages.size == 0 or not np.all(np.isfinite(ages)):
            raise ValueError("age_grid_years must be a finite one-dimensional array.")
        if np.any(ages < 0.0) or np.any(np.diff(ages) <= 0.0):
            raise ValueError(
                "age_grid_years must be non-negative and strictly increasing."
            )
        if masses.shape != (len(ids), ages.size):
            raise ValueError("masses must have one age-bin row per hypothesis.")
        if not np.all(np.isfinite(masses)) or np.any(masses < -1.0e-10):
            raise ValueError("masses must be finite and non-negative.")
        if not np.allclose(masses.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-8):
            raise ValueError("each TTD hypothesis must have unit probability mass.")
        probabilities: Optional[np.ndarray] = None
        if self.probabilities is not None:
            probabilities = np.asarray(self.probabilities, dtype=float).copy()
            if probabilities.shape != (len(ids),):
                raise ValueError("probabilities must match hypothesis_ids.")
            if not np.all(np.isfinite(probabilities)) or np.any(probabilities < 0.0):
                raise ValueError("probabilities must be finite and non-negative.")
            if float(probabilities.sum()) <= 0.0:
                raise ValueError("probabilities must have positive total mass.")
            probabilities /= float(probabilities.sum())
            probabilities.setflags(write=False)
        ages.setflags(write=False)
        masses.setflags(write=False)
        object.__setattr__(self, "hypothesis_ids", ids)
        object.__setattr__(self, "age_grid_years", ages)
        object.__setattr__(self, "masses", masses)
        object.__setattr__(self, "probabilities", probabilities)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class TracerDesignScenario:
    """Declared sensitivity scenario for one candidate measurement."""

    name: str
    response_multiplier: float = 1.0
    response_offset: float = 0.0
    standard_deviation_multiplier: float = 1.0
    weight: float = 1.0

    def __post_init__(self) -> None:
        values = (
            self.response_multiplier,
            self.response_offset,
            self.standard_deviation_multiplier,
            self.weight,
        )
        if not str(self.name).strip() or not all(
            math.isfinite(float(v)) for v in values
        ):
            raise ValueError("design scenario names and values must be finite.")
        if self.standard_deviation_multiplier <= 0.0 or self.weight <= 0.0:
            raise ValueError(
                "scenario uncertainty multipliers and weights must be positive."
            )


@dataclass(frozen=True)
class CandidateTracerMeasurement:
    """A tracer, well/target, and sampling time that could actually be measured."""

    option_id: str
    tracer: str
    target_id: str
    sample_year: float
    standard_deviation: float
    cost: float = 1.0
    feasible: bool = True
    scenarios: Sequence[TracerDesignScenario] = ()
    kernel_kwargs: Mapping[str, Any] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not str(self.option_id).strip() or not str(self.tracer).strip():
            raise ValueError("option_id and tracer must be non-empty.")
        if not str(self.target_id).strip():
            raise ValueError("target_id must be non-empty.")
        if not math.isfinite(float(self.sample_year)):
            raise ValueError("sample_year must be finite.")
        if (
            not math.isfinite(float(self.standard_deviation))
            or self.standard_deviation <= 0
        ):
            raise ValueError("standard_deviation must be finite and positive.")
        if not math.isfinite(float(self.cost)) or self.cost <= 0:
            raise ValueError("cost must be finite and positive.")
        names = [scenario.name for scenario in self.scenarios]
        if len(set(names)) != len(names):
            raise ValueError("candidate scenario names must be unique.")
        object.__setattr__(self, "scenarios", tuple(self.scenarios))
        object.__setattr__(self, "kernel_kwargs", dict(self.kernel_kwargs))
        object.__setattr__(self, "metadata", dict(self.metadata))


def _probability_gate(ensemble: TtdHypothesisEnsemble) -> Optional[dict[str, Any]]:
    if ensemble.probabilities is None or ensemble.probability_semantics is None:
        return {
            "status": "ABSTAIN",
            "reason": "no_probability_model",
            "selected_option_id": None,
            "rankings": [],
        }
    semantics = str(ensemble.probability_semantics).strip().lower()
    if semantics not in _PROBABILITY_SEMANTICS:
        return {
            "status": "ABSTAIN",
            "reason": "unsupported_probability_semantics",
            "declared_probability_semantics": ensemble.probability_semantics,
            "supported_probability_semantics": sorted(_PROBABILITY_SEMANTICS),
            "selected_option_id": None,
            "rankings": [],
        }
    return None


def build_ttd_measurement_options(
    ensemble: TtdHypothesisEnsemble,
    candidates: Sequence[CandidateTracerMeasurement],
) -> list[MeasurementOption]:
    """Generate candidate-specific forward predictions from the TTD ensemble."""
    options: list[MeasurementOption] = []
    for candidate in candidates:
        kernel = tracer_response_kernel(
            candidate.tracer,
            ensemble.age_grid_years,
            candidate.sample_year,
            **candidate.kernel_kwargs,
        )
        nominal_means = np.asarray(ensemble.masses) @ kernel
        declared_scenarios = candidate.scenarios or (
            TracerDesignScenario(name="nominal"),
        )
        scenarios = [
            PredictiveScenario(
                name=scenario.name,
                means=(
                    scenario.response_multiplier * nominal_means
                    + scenario.response_offset
                ),
                standard_deviations=(
                    float(candidate.standard_deviation)
                    * float(scenario.standard_deviation_multiplier)
                ),
                weight=float(scenario.weight),
            )
            for scenario in declared_scenarios
        ]
        options.append(
            MeasurementOption(
                option_id=candidate.option_id,
                measurement_type="groundwater_tracer_concentration",
                target_id=candidate.target_id,
                cost=float(candidate.cost),
                scenarios=scenarios,
                feasible=bool(candidate.feasible),
                metadata={
                    **candidate.metadata,
                    "tracer": candidate.tracer,
                    "sample_year": float(candidate.sample_year),
                    "probability_semantics": ensemble.probability_semantics,
                },
            )
        )
    return options


def _decision_values(
    ensemble: TtdHypothesisEnsemble,
    functionals: Optional[Sequence[AgeFunctional]],
) -> Optional[np.ndarray]:
    if functionals is None:
        return None
    columns: list[np.ndarray] = []
    for functional in functionals:
        coefficients = np.asarray(functional.coefficients, dtype=float)
        if coefficients.shape != np.asarray(ensemble.age_grid_years).shape:
            raise ValueError(
                f"Functional {functional.name!r} does not match the ensemble age grid."
            )
        columns.append(np.asarray(ensemble.masses) @ coefficients)
    return np.column_stack(columns) if columns else None


def rank_ttd_measurements(
    ensemble: TtdHypothesisEnsemble,
    candidates: Sequence[CandidateTracerMeasurement],
    *,
    decision_functionals: Optional[Sequence[AgeFunctional]] = None,
    config: Optional[AcquisitionConfig] = None,
) -> dict[str, Any]:
    """Rank tracer/time/well actions, or abstain when probabilities are absent."""
    gate = _probability_gate(ensemble)
    if gate is not None:
        return gate
    result = rank_measurement_options(
        ensemble.hypothesis_ids,
        np.asarray(ensemble.probabilities),
        build_ttd_measurement_options(ensemble, candidates),
        decision_values=_decision_values(ensemble, decision_functionals),
        config=config,
    )
    result["probability_semantics"] = ensemble.probability_semantics
    result["groundwater_forward_operator"] = "tracer_response_kernel"
    return result


def select_ttd_measurement_batch(
    ensemble: TtdHypothesisEnsemble,
    candidates: Sequence[CandidateTracerMeasurement],
    *,
    batch_size: int,
    decision_functionals: Optional[Sequence[AgeFunctional]] = None,
    config: Optional[AcquisitionConfig] = None,
) -> dict[str, Any]:
    """Select a non-redundant batch under the same probability-provenance gate."""
    gate = _probability_gate(ensemble)
    if gate is not None:
        gate["selected_option_ids"] = []
        return gate
    if decision_functionals is not None:
        raise ValueError(
            "The current batch engine supports joint information gain but not "
            "decision-weighted batch utility; omit decision_functionals."
        )
    result = select_measurement_batch(
        ensemble.hypothesis_ids,
        np.asarray(ensemble.probabilities),
        build_ttd_measurement_options(ensemble, candidates),
        batch_size=batch_size,
        config=config,
    )
    result["probability_semantics"] = ensemble.probability_semantics
    result["groundwater_forward_operator"] = "tracer_response_kernel"
    return result


__all__ = [
    "CandidateTracerMeasurement",
    "TracerDesignScenario",
    "TtdHypothesisEnsemble",
    "build_ttd_measurement_options",
    "rank_ttd_measurements",
    "select_ttd_measurement_batch",
]
