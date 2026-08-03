"""Set-valued transit-time-distribution inference from tracer observations.

The inferential object in this module is an age-bin probability-mass vector,
not a selected lumped-parameter-model family.  Tracer observations define a
convex feasible set through linear response constraints.  Linear programmes
then return sharp lower and upper bounds, with witness TTDs, for quantities the
data can support.  Infeasible or overly wide results are reported explicitly;
the implementation never converts numerical convergence into identifiability.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Any, Mapping, Optional, Sequence

import numpy as np
from scipy.optimize import linprog

from .input_history import InputHistory
from .joint_lpm import (
    TracerFitObservation,
    build_lpm_tracer_observations,
    tracer_response_kernel,
)

_SUPPORTED_INTERVAL_LIKELIHOODS = {
    "gaussian",
    "upper_censored",
    "lower_censored",
}


def _readonly_vector(values: Sequence[float], *, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float).copy()
    if vector.ndim != 1 or vector.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values.")
    vector.setflags(write=False)
    return vector


@dataclass(frozen=True)
class TracerConstraint:
    """One tracer response row and its declared observation uncertainty."""

    tracer: str
    response: Sequence[float]
    observed: float
    sigma: float
    likelihood: str = "gaussian"
    units: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "response", _readonly_vector(self.response, name="response")
        )
        if not str(self.tracer).strip():
            raise ValueError("tracer must be non-empty.")
        if not math.isfinite(float(self.observed)):
            raise ValueError("observed must be finite.")
        if not math.isfinite(float(self.sigma)) or float(self.sigma) <= 0.0:
            raise ValueError("sigma must be finite and strictly positive.")
        likelihood = str(self.likelihood).strip().lower()
        if likelihood not in _SUPPORTED_INTERVAL_LIKELIHOODS:
            raise ValueError(
                f"Unsupported interval likelihood {self.likelihood!r}; supported "
                f"values are {sorted(_SUPPORTED_INTERVAL_LIKELIHOODS)}."
            )
        object.__setattr__(self, "likelihood", likelihood)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class AgeFunctional:
    """A linear, decision-relevant summary of age-bin probability mass."""

    name: str
    coefficients: Sequence[float]
    maximum_reportable_width: Optional[float] = None
    units: str = "fraction"

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "coefficients",
            _readonly_vector(self.coefficients, name="coefficients"),
        )
        if not str(self.name).strip():
            raise ValueError("functional name must be non-empty.")
        threshold = self.maximum_reportable_width
        if threshold is not None and (
            not math.isfinite(float(threshold)) or float(threshold) < 0.0
        ):
            raise ValueError(
                "maximum_reportable_width must be finite and non-negative."
            )


@dataclass(frozen=True)
class IdentifiedBound:
    """Sharp feasible bounds and optimizer witnesses for one functional."""

    name: str
    lower: float
    upper: float
    width: float
    status: str
    units: str
    maximum_reportable_width: Optional[float]
    lower_witness: np.ndarray
    upper_witness: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "lower_witness",
            _readonly_vector(self.lower_witness, name="lower_witness"),
        )
        object.__setattr__(
            self,
            "upper_witness",
            _readonly_vector(self.upper_witness, name="upper_witness"),
        )

    def to_dict(self, *, include_witnesses: bool = True) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "name": self.name,
            "lower": float(self.lower),
            "upper": float(self.upper),
            "width": float(self.width),
            "status": self.status,
            "units": self.units,
            "maximum_reportable_width": self.maximum_reportable_width,
        }
        if include_witnesses:
            payload["lower_witness"] = self.lower_witness.tolist()
            payload["upper_witness"] = self.upper_witness.tolist()
        return payload


@dataclass(frozen=True)
class TtdEvidenceReport:
    """Auditable result for one declared TTD identified-set scenario."""

    status: str
    scenario_name: str
    age_grid_years: np.ndarray
    constraints: tuple[TracerConstraint, ...]
    bounds: Mapping[str, IdentifiedBound]
    response_rank: int
    normalization_augmented_rank: int
    nullity: int
    sigma_multiplier: float
    feasibility_tolerance: float
    feasible_witness: Optional[np.ndarray]
    excluded_observations: tuple[Mapping[str, Any], ...] = ()
    assumptions: tuple[str, ...] = ()
    provenance: Mapping[str, Any] = field(default_factory=dict)
    message: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "age_grid_years",
            _readonly_vector(self.age_grid_years, name="age_grid_years"),
        )
        if self.feasible_witness is not None:
            object.__setattr__(
                self,
                "feasible_witness",
                _readonly_vector(self.feasible_witness, name="feasible_witness"),
            )
        object.__setattr__(self, "constraints", tuple(self.constraints))
        object.__setattr__(self, "bounds", dict(self.bounds))
        object.__setattr__(
            self,
            "excluded_observations",
            tuple(dict(item) for item in self.excluded_observations),
        )
        object.__setattr__(self, "assumptions", tuple(self.assumptions))
        object.__setattr__(self, "provenance", dict(self.provenance))

    def to_dict(self, *, include_witnesses: bool = True) -> dict[str, Any]:
        return {
            "status": self.status,
            "scenario_name": self.scenario_name,
            "age_grid_years": self.age_grid_years.tolist(),
            "bounds": {
                name: bound.to_dict(include_witnesses=include_witnesses)
                for name, bound in self.bounds.items()
            },
            "response_rank": int(self.response_rank),
            "normalization_augmented_rank": int(self.normalization_augmented_rank),
            "nullity": int(self.nullity),
            "sigma_multiplier": float(self.sigma_multiplier),
            "feasibility_tolerance": float(self.feasibility_tolerance),
            "feasible_witness": (
                self.feasible_witness.tolist()
                if include_witnesses and self.feasible_witness is not None
                else None
            ),
            "tracers_used": [constraint.tracer for constraint in self.constraints],
            "excluded_observations": [
                dict(item) for item in self.excluded_observations
            ],
            "assumptions": list(self.assumptions),
            "provenance": dict(self.provenance),
            "message": self.message,
        }


def standard_age_functionals(
    age_grid_years: Sequence[float],
    *,
    cutoffs_years: Sequence[float] = (70.0, 11700.0),
    maximum_fraction_widths: Optional[Mapping[str, float]] = None,
    include_mean_age: bool = False,
    maximum_mean_age_width_years: Optional[float] = None,
) -> tuple[AgeFunctional, ...]:
    """Build the M3 age fractions without imposing hidden reportability gates."""
    ages = _readonly_vector(age_grid_years, name="age_grid_years")
    cutoffs = np.asarray(cutoffs_years, dtype=float)
    if (
        cutoffs.shape != (2,)
        or not np.all(np.isfinite(cutoffs))
        or cutoffs[0] < 0.0
        or cutoffs[1] <= cutoffs[0]
    ):
        raise ValueError(
            "cutoffs_years must contain two finite, increasing, non-negative ages."
        )
    young_cutoff, old_cutoff = (float(value) for value in cutoffs)
    thresholds = dict(maximum_fraction_widths or {})
    definitions = {
        "anthropocene": ages <= young_cutoff,
        "holocene": (ages > young_cutoff) & (ages <= old_cutoff),
        "pleistocene": ages > old_cutoff,
    }
    functionals = [
        AgeFunctional(
            name=name,
            coefficients=mask.astype(float),
            maximum_reportable_width=thresholds.get(name),
            units="fraction",
        )
        for name, mask in definitions.items()
    ]
    if include_mean_age:
        functionals.append(
            AgeFunctional(
                name="mean_age_years",
                coefficients=ages,
                maximum_reportable_width=maximum_mean_age_width_years,
                units="years",
            )
        )
    return tuple(functionals)


def build_tracer_response_matrix(
    tracers: Sequence[str],
    age_grid_years: Sequence[float],
    sample_year: float,
    **kernel_kwargs: Any,
) -> np.ndarray:
    """Stack one linear tracer-response row per requested tracer."""
    rows = [
        tracer_response_kernel(
            tracer,
            age_grid_years,
            sample_year,
            **kernel_kwargs,
        )
        for tracer in tracers
    ]
    if not rows:
        return np.empty((0, len(age_grid_years)), dtype=float)
    return np.vstack(rows)


def build_tracer_constraints(
    observations: Mapping[str, Any] | Sequence[TracerFitObservation],
    age_grid_years: Sequence[float],
    sample_year: float,
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    use_helium4: bool = False,
    initial_c14_pmc: float = 100.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
) -> tuple[tuple[TracerConstraint, ...], tuple[Mapping[str, Any], ...]]:
    """Translate existing M3 observations into auditable linear constraints.

    Robust contaminated-mixture observations do not have equivalent linear
    interval semantics and are therefore excluded explicitly rather than being
    silently treated as Gaussian evidence.
    """
    if isinstance(observations, Mapping):
        fit_observations = build_lpm_tracer_observations(
            observations, use_helium4=use_helium4
        )
    else:
        fit_observations = list(observations)

    constraints: list[TracerConstraint] = []
    excluded: list[Mapping[str, Any]] = []
    for observation in fit_observations:
        likelihood = str(observation.likelihood).strip().lower()
        if likelihood not in _SUPPORTED_INTERVAL_LIKELIHOODS:
            excluded.append(
                {
                    "tracer": observation.tracer,
                    "likelihood": likelihood,
                    "reason": "no_declared_linear_interval_semantics",
                }
            )
            continue
        response = tracer_response_kernel(
            observation.tracer,
            age_grid_years,
            sample_year,
            histories=histories,
            initial_c14_pmc=initial_c14_pmc,
            helium4_background_ccpg=helium4_background_ccpg,
            helium4_accumulation_rate_ccpg_per_year=(
                helium4_accumulation_rate_ccpg_per_year
            ),
            prediction_scale_factors=prediction_scale_factors,
        )
        constraints.append(
            TracerConstraint(
                tracer=observation.tracer,
                response=response,
                observed=observation.value,
                sigma=observation.sigma,
                likelihood=likelihood,
                units=observation.units,
                metadata={
                    "note": observation.note,
                    "legacy_fit_weight": float(observation.weight),
                },
            )
        )
    return tuple(constraints), tuple(excluded)


def _compile_constraints(
    constraints: Sequence[TracerConstraint],
    n_age_bins: int,
    sigma_multiplier: float,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    rows: list[np.ndarray] = []
    limits: list[float] = []
    for constraint in constraints:
        response = np.asarray(constraint.response, dtype=float)
        if response.shape != (n_age_bins,):
            raise ValueError(
                f"Constraint {constraint.tracer!r} has {response.size} response "
                f"values for {n_age_bins} age bins."
            )
        margin = float(sigma_multiplier) * float(constraint.sigma)
        if constraint.likelihood in {"gaussian", "upper_censored"}:
            rows.append(response)
            limits.append(float(constraint.observed) + margin)
        if constraint.likelihood in {"gaussian", "lower_censored"}:
            rows.append(-response)
            limits.append(-(float(constraint.observed) - margin))
    if not rows:
        return None, None
    return np.vstack(rows), np.asarray(limits, dtype=float)


def _audit_witness(
    witness: np.ndarray,
    a_ub: Optional[np.ndarray],
    b_ub: Optional[np.ndarray],
    tolerance: float,
) -> None:
    if not np.all(np.isfinite(witness)):
        raise RuntimeError("Linear-program witness contains non-finite values.")
    if float(np.min(witness)) < -tolerance:
        raise RuntimeError("Linear-program witness violates non-negativity.")
    if abs(float(np.sum(witness)) - 1.0) > tolerance:
        raise RuntimeError("Linear-program witness violates unit probability mass.")
    if a_ub is not None and b_ub is not None:
        violation = float(np.max(a_ub @ witness - b_ub))
        if violation > tolerance:
            raise RuntimeError(
                "Linear-program witness violates tracer constraints by "
                f"{violation:.3g}."
            )


def _solve_lp(
    objective: np.ndarray,
    a_ub: Optional[np.ndarray],
    b_ub: Optional[np.ndarray],
    tolerance: float,
) -> Any:
    result = linprog(
        objective,
        A_ub=a_ub,
        b_ub=b_ub,
        A_eq=np.ones((1, objective.size), dtype=float),
        b_eq=np.ones(1, dtype=float),
        bounds=(0.0, None),
        method="highs",
    )
    if result.success:
        _audit_witness(np.asarray(result.x, dtype=float), a_ub, b_ub, tolerance)
    return result


def solve_ttd_identified_set(
    age_grid_years: Sequence[float],
    constraints: Sequence[TracerConstraint],
    functionals: Sequence[AgeFunctional],
    *,
    sigma_multiplier: float = 1.96,
    feasibility_tolerance: float = 1.0e-7,
    scenario_name: str = "declared",
    excluded_observations: Sequence[Mapping[str, Any]] = (),
    assumptions: Sequence[str] = (),
    provenance: Optional[Mapping[str, Any]] = None,
) -> TtdEvidenceReport:
    """Compute sharp functional bounds over the tracer-consistent TTD set."""
    ages = _readonly_vector(age_grid_years, name="age_grid_years")
    if np.any(np.diff(ages) <= 0.0):
        raise ValueError("age_grid_years must be strictly increasing.")
    if not math.isfinite(float(sigma_multiplier)) or sigma_multiplier < 0.0:
        raise ValueError("sigma_multiplier must be finite and non-negative.")
    if not math.isfinite(float(feasibility_tolerance)) or feasibility_tolerance <= 0:
        raise ValueError("feasibility_tolerance must be finite and positive.")
    constraint_tuple = tuple(constraints)
    functional_tuple = tuple(functionals)
    names = [functional.name for functional in functional_tuple]
    if len(set(names)) != len(names):
        raise ValueError("functional names must be unique.")
    for functional in functional_tuple:
        if np.asarray(functional.coefficients).shape != ages.shape:
            raise ValueError(
                f"Functional {functional.name!r} does not match the age grid."
            )

    a_ub, b_ub = _compile_constraints(
        constraint_tuple, ages.size, float(sigma_multiplier)
    )
    response_matrix = (
        np.vstack([constraint.response for constraint in constraint_tuple])
        if constraint_tuple
        else np.empty((0, ages.size), dtype=float)
    )
    response_rank = int(np.linalg.matrix_rank(response_matrix))
    augmented = np.vstack([np.ones(ages.size, dtype=float), response_matrix])
    augmented_rank = int(np.linalg.matrix_rank(augmented))
    nullity = max(0, int(ages.size - augmented_rank))

    feasibility = _solve_lp(
        np.zeros(ages.size, dtype=float),
        a_ub,
        b_ub,
        float(feasibility_tolerance),
    )
    if not feasibility.success:
        return TtdEvidenceReport(
            status="ABSTAIN",
            scenario_name=scenario_name,
            age_grid_years=ages,
            constraints=constraint_tuple,
            bounds={},
            response_rank=response_rank,
            normalization_augmented_rank=augmented_rank,
            nullity=nullity,
            sigma_multiplier=float(sigma_multiplier),
            feasibility_tolerance=float(feasibility_tolerance),
            feasible_witness=None,
            excluded_observations=tuple(excluded_observations),
            assumptions=tuple(assumptions),
            provenance=dict(provenance or {}),
            message=(
                "No non-negative unit-mass TTD satisfies the declared tracer "
                f"intervals: {feasibility.message}"
            ),
        )

    bounds: dict[str, IdentifiedBound] = {}
    for functional in functional_tuple:
        coefficients = np.asarray(functional.coefficients, dtype=float)
        lower_result = _solve_lp(coefficients, a_ub, b_ub, float(feasibility_tolerance))
        upper_result = _solve_lp(
            -coefficients, a_ub, b_ub, float(feasibility_tolerance)
        )
        if not lower_result.success or not upper_result.success:
            raise RuntimeError(
                f"Functional optimization failed for {functional.name!r}: "
                f"lower={lower_result.message}; upper={upper_result.message}."
            )
        lower = float(lower_result.fun)
        upper = float(-upper_result.fun)
        if upper < lower and lower - upper <= feasibility_tolerance:
            lower, upper = upper, lower
        width = max(0.0, upper - lower)
        threshold = functional.maximum_reportable_width
        if threshold is None:
            bound_status = "BOUNDED"
        elif width <= float(threshold) + float(feasibility_tolerance):
            bound_status = "IDENTIFIED"
        else:
            bound_status = "UNRESOLVED"
        bounds[functional.name] = IdentifiedBound(
            name=functional.name,
            lower=lower,
            upper=upper,
            width=width,
            status=bound_status,
            units=functional.units,
            maximum_reportable_width=threshold,
            lower_witness=np.asarray(lower_result.x, dtype=float),
            upper_witness=np.asarray(upper_result.x, dtype=float),
        )

    identified_count = sum(bound.status == "IDENTIFIED" for bound in bounds.values())
    if bounds and identified_count == len(bounds):
        overall_status = "IDENTIFIED"
    elif identified_count > 0:
        overall_status = "PARTIALLY_IDENTIFIED"
    else:
        overall_status = "ABSTAIN"
    return TtdEvidenceReport(
        status=overall_status,
        scenario_name=scenario_name,
        age_grid_years=ages,
        constraints=constraint_tuple,
        bounds=bounds,
        response_rank=response_rank,
        normalization_augmented_rank=augmented_rank,
        nullity=nullity,
        sigma_multiplier=float(sigma_multiplier),
        feasibility_tolerance=float(feasibility_tolerance),
        feasible_witness=np.asarray(feasibility.x, dtype=float),
        excluded_observations=tuple(excluded_observations),
        assumptions=tuple(assumptions),
        provenance=dict(provenance or {}),
        message=(
            "All requested functionals meet their predeclared width gates."
            if overall_status == "IDENTIFIED"
            else (
                "The tracer-consistent set is feasible, but at least one requested "
                "functional lacks or exceeds a predeclared width gate."
                if overall_status == "PARTIALLY_IDENTIFIED"
                else "The tracer-consistent set is feasible, but no requested "
                "functional meets a predeclared width gate."
            )
        ),
    )


def infer_ttd_evidence(
    observations: Mapping[str, Any] | Sequence[TracerFitObservation],
    sample_year: float,
    age_grid_years: Sequence[float],
    *,
    functionals: Optional[Sequence[AgeFunctional]] = None,
    maximum_fraction_widths: Optional[Mapping[str, float]] = None,
    histories: Optional[Mapping[str, InputHistory]] = None,
    use_helium4: bool = False,
    initial_c14_pmc: float = 100.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
    sigma_multiplier: float = 1.96,
    feasibility_tolerance: float = 1.0e-7,
    scenario_name: str = "declared",
    assumptions: Sequence[str] = (),
    provenance: Optional[Mapping[str, Any]] = None,
) -> TtdEvidenceReport:
    """High-level bridge from the existing M3 observation contract."""
    constraints, excluded = build_tracer_constraints(
        observations,
        age_grid_years,
        sample_year,
        histories=histories,
        use_helium4=use_helium4,
        initial_c14_pmc=initial_c14_pmc,
        helium4_background_ccpg=helium4_background_ccpg,
        helium4_accumulation_rate_ccpg_per_year=(
            helium4_accumulation_rate_ccpg_per_year
        ),
        prediction_scale_factors=prediction_scale_factors,
    )
    requested = (
        tuple(functionals)
        if functionals is not None
        else standard_age_functionals(
            age_grid_years,
            maximum_fraction_widths=maximum_fraction_widths,
        )
    )
    if not constraints:
        return TtdEvidenceReport(
            status="ABSTAIN",
            scenario_name=scenario_name,
            age_grid_years=np.asarray(age_grid_years, dtype=float),
            constraints=(),
            bounds={},
            response_rank=0,
            normalization_augmented_rank=1,
            nullity=max(0, len(age_grid_years) - 1),
            sigma_multiplier=float(sigma_multiplier),
            feasibility_tolerance=float(feasibility_tolerance),
            feasible_witness=None,
            excluded_observations=excluded,
            assumptions=tuple(assumptions),
            provenance=dict(provenance or {}),
            message="No observation with declared linear interval semantics is available.",
        )
    return solve_ttd_identified_set(
        age_grid_years,
        constraints,
        requested,
        sigma_multiplier=sigma_multiplier,
        feasibility_tolerance=feasibility_tolerance,
        scenario_name=scenario_name,
        excluded_observations=excluded,
        assumptions=assumptions,
        provenance=provenance,
    )


def bound_linear_functional(
    report: TtdEvidenceReport,
    functional: AgeFunctional,
) -> IdentifiedBound:
    """Bound a new quantity without adding it to the fitted evidence."""
    if report.feasible_witness is None:
        raise ValueError(
            "Cannot bound a functional from an infeasible evidence report."
        )
    refreshed = solve_ttd_identified_set(
        report.age_grid_years,
        report.constraints,
        [functional],
        sigma_multiplier=report.sigma_multiplier,
        feasibility_tolerance=report.feasibility_tolerance,
        scenario_name=report.scenario_name,
        excluded_observations=report.excluded_observations,
        assumptions=report.assumptions,
        provenance=report.provenance,
    )
    return refreshed.bounds[functional.name]


def predict_tracer_bounds(
    report: TtdEvidenceReport,
    tracer: str,
    sample_year: float,
    *,
    name: Optional[str] = None,
    maximum_reportable_width: Optional[float] = None,
    histories: Optional[Mapping[str, InputHistory]] = None,
    **kernel_kwargs: Any,
) -> IdentifiedBound:
    """Return held-out predictive bounds without conditioning on that tracer."""
    response = tracer_response_kernel(
        tracer,
        report.age_grid_years,
        sample_year,
        histories=histories,
        **kernel_kwargs,
    )
    return bound_linear_functional(
        report,
        AgeFunctional(
            name=name or f"predicted_{tracer}",
            coefficients=response,
            maximum_reportable_width=maximum_reportable_width,
            units="tracer_native",
        ),
    )


def assess_held_out_tracer(
    report: TtdEvidenceReport,
    observation: TracerFitObservation,
    sample_year: float,
    *,
    sigma_multiplier: Optional[float] = None,
    histories: Optional[Mapping[str, InputHistory]] = None,
    **kernel_kwargs: Any,
) -> dict[str, Any]:
    """Check held-out concentration compatibility without refitting the TTD set."""
    prediction = predict_tracer_bounds(
        report,
        observation.tracer,
        sample_year,
        histories=histories,
        **kernel_kwargs,
    )
    multiplier = (
        report.sigma_multiplier if sigma_multiplier is None else float(sigma_multiplier)
    )
    margin = multiplier * float(observation.sigma)
    likelihood = str(observation.likelihood).strip().lower()
    if likelihood == "gaussian":
        observed_lower = float(observation.value) - margin
        observed_upper = float(observation.value) + margin
        compatible = (
            prediction.upper >= observed_lower and prediction.lower <= observed_upper
        )
    elif likelihood == "upper_censored":
        observed_lower = float("-inf")
        observed_upper = float(observation.value) + margin
        compatible = prediction.lower <= observed_upper
    elif likelihood == "lower_censored":
        observed_lower = float(observation.value) - margin
        observed_upper = float("inf")
        compatible = prediction.upper >= observed_lower
    else:
        return {
            "status": "ABSTAIN",
            "tracer": observation.tracer,
            "reason": "no_declared_linear_interval_semantics",
        }
    return {
        "status": "COMPATIBLE" if compatible else "INCOMPATIBLE",
        "tracer": observation.tracer,
        "prediction_lower": prediction.lower,
        "prediction_upper": prediction.upper,
        "observed_interval_lower": observed_lower,
        "observed_interval_upper": observed_upper,
        "sigma_multiplier": multiplier,
        "conditioned_on_held_out_observation": False,
    }


__all__ = [
    "AgeFunctional",
    "IdentifiedBound",
    "TracerConstraint",
    "TtdEvidenceReport",
    "assess_held_out_tracer",
    "bound_linear_functional",
    "build_tracer_constraints",
    "build_tracer_response_matrix",
    "infer_ttd_evidence",
    "predict_tracer_bounds",
    "solve_ttd_identified_set",
    "standard_age_functionals",
]
