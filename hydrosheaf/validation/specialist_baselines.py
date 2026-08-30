"""Truth-blind specialist outputs for the programme benchmark.

The topology baselines in :mod:`hydrosheaf.validation.baselines` operate on
candidate edges.  Age and reaction specialists have different output spaces,
so they are kept in this small companion module rather than forcing every
baseline into an edge-probability interface.

These are deliberately fixed, auditable diagnostic comparators.  They do not
reuse HydroSheaf's posterior or sheaf code, do not receive generator truth,
and do not claim to be a production PHREEQC inverse.  The registry metadata
records that limitation so a benchmark cannot silently promote these rules to
a competence-matched superiority comparator.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from collections.abc import Callable, Iterable, Mapping
from typing import Any

from hydrosheaf.nuclear.multi_tracer import ATMOSPHERIC_TRACER_HISTORIES

from .age_competent_baseline import AgeCompetentBaseline
from .baselines import assert_truth_blind_observations
from .metrics import classification_metrics, interval_coverage, regression_metrics
from .reaction_competent_baseline import ReactionCompetentBaseline


AGE_OUTPUT = "age_interval"
REACTION_OUTPUT = "reaction_family"
ABSTAIN = "abstain"
SELECT = "select"

_AGE_MAX_YEARS = 200.0
_TRITIUM_MODERN_TU = 8.0
_TRITIUM_HALF_LIFE_YEARS = 12.32
_ARGON39_MODERN_PMC = 100.0
_ARGON39_HALF_LIFE_YEARS = 269.0
_AGE_Z95 = 1.959963984540054
_ATMOSPHERIC_AGE_MAX_YEARS = 120.0
_ATMOSPHERIC_AGE_GRID_STEP_YEARS = 0.5
_ATMOSPHERIC_TRACERS = {
    "cfc11_pptv": "CFC11",
    "cfc12_pptv": "CFC12",
    "cfc113_pptv": "CFC113",
    "sf6_pptv": "SF6",
}
_REACTION_FAMILIES = (
    "carbonate",
    "silicate_exchange",
    "sulfate_reduction",
    "iron_reduction",
    "denitrification",
    "sulfate_source",
    "other_redox",
    "other",
    "none",
)


def _finite(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, float(value)))


def _jsonable(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True, default=str))


def _fingerprint(record: Mapping[str, Any]) -> str:
    payload = json.dumps(record, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SpecialistBaselineSpec:
    """A non-topology baseline with a declared output schema."""

    name: str
    version: str
    family: str
    output_kind: str
    input_channels: tuple[str, ...]
    tuning: Mapping[str, Any]
    uncertainty: Mapping[str, Any]
    abstention: Mapping[str, Any]
    cost: Mapping[str, Any]
    predictor: Callable[[Mapping[str, Any]], Mapping[str, Mapping[str, Any]]]
    description: str = ""
    control: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name or not self.version or not self.family:
            raise ValueError("Specialist baseline identity fields are required.")
        if self.output_kind not in {AGE_OUTPUT, REACTION_OUTPUT}:
            raise ValueError(f"Unsupported specialist output kind: {self.output_kind!r}")
        if not self.input_channels:
            raise ValueError("Specialist baselines require input channels.")
        for section in (self.tuning, self.uncertainty, self.abstention, self.cost):
            if not isinstance(section, Mapping):
                raise TypeError("Specialist baseline metadata must be mappings.")

    def predict(self, observations: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
        """Predict from observed channels only and validate output shape."""

        assert_truth_blind_observations(observations)
        output = self.predictor(observations)
        if not isinstance(output, Mapping):
            raise TypeError(f"{self.name} predictor must return a mapping.")
        normalised: dict[str, dict[str, Any]] = {}
        for target, record in output.items():
            if not isinstance(record, Mapping):
                raise TypeError(f"{self.name} output for {target!r} is not a mapping.")
            normalised[str(target)] = dict(record)
        return normalised

    def to_audit_record(self) -> dict[str, Any]:
        record: dict[str, Any] = {
            "name": self.name,
            "version": self.version,
            "family": self.family,
            "output_kind": self.output_kind,
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


class SpecialistBaselineRegistry:
    """Name-addressed registry for age and reaction specialist comparators."""

    def __init__(self, specs: Iterable[SpecialistBaselineSpec] = ()) -> None:
        self._specs: dict[str, SpecialistBaselineSpec] = {}
        for spec in specs:
            self.register(spec)

    def register(self, spec: SpecialistBaselineSpec) -> SpecialistBaselineSpec:
        if spec.name in self._specs:
            raise ValueError(f"Baseline already registered: {spec.name}")
        # Validate metadata and fingerprint at registration time.
        spec.to_audit_record()
        self._specs[spec.name] = spec
        return spec

    def get(self, name: str) -> SpecialistBaselineSpec:
        try:
            return self._specs[name]
        except KeyError as exc:
            raise KeyError(f"Unknown specialist baseline: {name}") from exc

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._specs))

    def specs(self) -> tuple[SpecialistBaselineSpec, ...]:
        return tuple(self._specs[name] for name in self.names())

    def audit_table(self) -> tuple[dict[str, Any], ...]:
        return tuple(spec.to_audit_record() for spec in self.specs())


def local_tracer_age_baseline_spec() -> SpecialistBaselineSpec:
    """Return a fixed local decay baseline using observed age tracers.

    Tritium and argon-39 are converted independently with fixed modern
    reference activities and physical half-lives.  Available estimates are
    precision-weighted; missing tracers produce an explicit wide interval and
    abstention.  No input-history lookup, generator truth, or HydroSheaf age
    posterior is used.
    """

    return SpecialistBaselineSpec(
        name="local_tracer_decay_age",
        version="1.0",
        family="local_tracer_age",
        output_kind=AGE_OUTPUT,
        input_channels=("tracer_age",),
        tuning={
            "reference_tritium_TU": _TRITIUM_MODERN_TU,
            "reference_argon39_pmc": _ARGON39_MODERN_PMC,
            "max_age_years": _AGE_MAX_YEARS,
            "half_life_rule": "physical_decay_only",
            "tuning_data": "none_fixed_constants",
        },
        uncertainty={
            "type": "fixed_measurement_plus_forward_model_interval",
            "calibrated": False,
            "age_sigma_rule": "max(5 years, 0.20 * age)",
        },
        abstention={
            "rule": "abstain when no declared age tracer is finite",
            "missing_evidence_decision": ABSTAIN,
        },
        cost={
            "false_commitment": 1.0,
            "abstain": 0.10,
            "measurement": 2.0,
        },
        predictor=_predict_local_tracer_age,
        description=(
            "Fixed local tracer-decay comparator; diagnostic only and not a "
            "full atmospheric-history inverse."
        ),
        control={"truth_blind": True, "candidate_universe_required": False},
    )


def atmospheric_history_multitracer_age_baseline_spec() -> SpecialistBaselineSpec:
    """Return a fixed atmospheric-history, multi-tracer age comparator.

    This is a stronger specialist diagnostic than the local decay control. It
    evaluates a bounded recharge-age grid against fixed atmospheric histories
    for CFCs/SF6 and physical decay curves for tritium/argon-39. It is still
    deliberately simpler than a site-specific LPM inverse: recharge mixing,
    excess air, degradation, and local atmospheric histories are represented as
    explicit limitations and trigger conservative interpretation.
    """

    return SpecialistBaselineSpec(
        name="multitracer_atmospheric_history_age",
        version="1.0",
        family="multitracer_age",
        output_kind=AGE_OUTPUT,
        input_channels=("tracer_age_history",),
        tuning={
            "max_age_years": _ATMOSPHERIC_AGE_MAX_YEARS,
            "grid_step_years": _ATMOSPHERIC_AGE_GRID_STEP_YEARS,
            "atmospheric_histories": sorted(ATMOSPHERIC_TRACER_HISTORIES),
            "physical_decay_tracers": ["tritium_TU", "argon39_pmc"],
            "tuning_data": "fixed_public_screening_histories_and_half_lives",
        },
        uncertainty={
            "type": "discrete_age_grid_likelihood_interval",
            "calibrated": False,
            "measurement_sigma_rule": "declared_sigma_or_max(relative_error, absolute_floor)",
            "multimodality_reported": True,
        },
        abstention={
            "rule": "abstain when no supported tracer is finite",
            "missing_evidence_decision": ABSTAIN,
            "model_disagreement": "reported_when_posterior_is_multimodal",
        },
        cost={
            "false_commitment": 1.0,
            "abstain": 0.10,
            "measurement": 2.5,
        },
        predictor=_predict_atmospheric_history_multitracer_age,
        description=(
            "Fixed multi-tracer atmospheric-history screening comparator with "
            "an explicit age grid; not a site-specific LPM, excess-air, or "
            "degradation-corrected inverse."
        ),
        control={
            "truth_blind": True,
            "candidate_universe_required": False,
            "null_control": "local_tracer_decay_age",
        },
    )


def local_thermodynamic_reaction_baseline_spec() -> SpecialistBaselineSpec:
    """Return a fixed local chemistry rule comparator.

    The comparator uses only concentration changes across a candidate edge and
    emits a probability distribution over declared reaction families.  It is
    intentionally marked as a rule-based thermodynamic diagnostic, not as a
    PHREEQC inverse or as a replacement for a calibrated reaction model.
    """

    return SpecialistBaselineSpec(
        name="local_thermodynamic_reaction_rules",
        version="1.0",
        family="reaction_family",
        output_kind=REACTION_OUTPUT,
        input_channels=("reaction_chemistry",),
        tuning={
            "delta_thresholds": {
                "major_ion": 0.05,
                "redox": 0.05,
                "iron": 0.01,
            },
            "softmax_temperature": 1.0,
            "tuning_data": "none_fixed_rules",
            "candidate_universe": "caller_supplied_observed_edges",
        },
        uncertainty={
            "type": "rule_score_softmax",
            "calibrated": False,
            "probabilities": "diagnostic_not_calibrated",
        },
        abstention={
            "rule": "abstain below 0.60 maximum family probability or missing pair",
            "select_threshold": 0.60,
            "missing_evidence_decision": ABSTAIN,
        },
        cost={
            "false_commitment": 1.0,
            "abstain": 0.10,
            "measurement": 3.0,
        },
        predictor=_predict_local_reaction,
        description=(
            "Fixed local concentration-difference rules with a softmax over "
            "reaction families; not a full PHREEQC inverse."
        ),
        control={"truth_blind": True, "candidate_universe_required": True},
    )


def stoichiometric_reaction_baseline_spec() -> SpecialistBaselineSpec:
    """Return a constrained stoichiometric reaction-family comparator.

    The predictor projects observed concentration changes onto fixed signed
    reaction templates with non-negative extents and residual diagnostics. It
    is independent of HydroSheaf's reaction posterior and is stronger than a
    one-ion threshold rule, while remaining explicitly short of a nonlinear
    PHREEQC inverse because activities, mineral saturation, and kinetic rates
    are not solved here.
    """

    return SpecialistBaselineSpec(
        name="stoichiometric_reaction_inverse",
        version="1.0",
        family="stoichiometric_reaction",
        output_kind=REACTION_OUTPUT,
        input_channels=("reaction_chemistry",),
        tuning={
            "template_library": "fixed_signed_ion_balance_v1",
            "extent_constraint": "non_negative_projection",
            "residual_scale": "declared_ion_measurement_scales",
            "null_reaction": "zero_change_template_with_complexity_penalty",
            "reaction_complexity_penalty": 0.25,
            "softmax_temperature": 1.0,
            "candidate_universe": "caller_supplied_observed_edges",
            "tuning_data": "none_fixed_templates",
        },
        uncertainty={
            "type": "residual_softmax_diagnostic",
            "calibrated": False,
            "probabilities": "requires_development_temperature_calibration",
            "residual_reported": True,
        },
        abstention={
            "rule": "abstain below 0.60 maximum probability or insufficient paired ions",
            "select_threshold": 0.60,
            "missing_evidence_decision": ABSTAIN,
        },
        cost={
            "false_commitment": 1.0,
            "abstain": 0.10,
            "measurement": 3.0,
        },
        predictor=_predict_stoichiometric_reaction,
        description=(
            "Fixed non-negative stoichiometric-template comparator with residual "
            "and evidence diagnostics; not a coupled nonlinear PHREEQC inverse."
        ),
        control={
            "truth_blind": True,
            "candidate_universe_required": True,
            "null_control": "local_thermodynamic_reaction_rules",
        },
    )


def competence_matched_age_baseline_spec() -> SpecialistBaselineSpec:
    """Return the independent bounded-distribution age specialist.

    The candidate universe and forward observation model live in the separate
    ``age_competent_baseline`` module.  This adapter only exposes that model
    through the programme's registry contract; it never supplies truth or a
    HydroSheaf posterior.
    """

    baseline = AgeCompetentBaseline()
    audit = baseline.to_audit_record()
    return SpecialistBaselineSpec(
        name=baseline.name,
        version=baseline.version,
        family="age_specialist",
        output_kind=AGE_OUTPUT,
        input_channels=("tracer_age_history",),
        tuning={
            "candidate_generator": audit["candidate_universe"],
            "history_source": audit["history_source"],
            "decay_constants": audit["decay_constants"],
            "tuning_data": "fixed_public_histories_and_decay_constants",
        },
        uncertainty=dict(audit["uncertainty"]),
        abstention=dict(audit["abstention"]),
        cost={
            "false_commitment": 1.0,
            "abstain": 0.10,
            "measurement": 2.5,
        },
        predictor=baseline.predict,
        description=(
            "Independent bounded exponential/gamma age-distribution specialist "
            "with fixed tracer histories; not a field-calibrated LPM."
        ),
        control={
            "truth_blind": True,
            "candidate_universe_required": False,
            "candidate_universe_scope": "all eligible observed nodes",
            "candidate_universe_hash": audit["candidate_universe"]["candidate_hash"],
            "uses_hydrosheaf_posterior": False,
        },
    )


def competence_matched_reaction_baseline_spec() -> SpecialistBaselineSpec:
    """Return the independent fixed-template reaction specialist.

    This adapter preserves the baseline's explicit null class, non-negative
    extent projection, residual diagnostics, and uncalibrated-probability
    warning.  It is a competence-matched validation comparator only when the
    benchmark declares the observed edge universe independently.
    """

    baseline = ReactionCompetentBaseline()
    audit = baseline.to_audit_record()
    return SpecialistBaselineSpec(
        name=baseline.name,
        version=baseline.version,
        family="reaction_family",
        output_kind=REACTION_OUTPUT,
        input_channels=("reaction_chemistry",),
        tuning=dict(audit["tuning"]),
        uncertainty=dict(audit["uncertainty"]),
        abstention=dict(audit["abstention"]),
        cost={
            "false_commitment": 1.0,
            "abstain": 0.10,
            "measurement": 3.0,
        },
        predictor=baseline.predict,
        description=(
            "Independent constrained reaction-family candidate generator with "
            "an explicit null model; not a coupled PHREEQC inverse."
        ),
        control={
            "truth_blind": True,
            "candidate_universe_required": True,
            "candidate_universe_scope": "independent observed chemistry edges",
            "uses_hydrosheaf_posterior": False,
            "uses_phreeqc_outputs_as_truth": False,
        },
    )


def default_specialist_baseline_registry() -> SpecialistBaselineRegistry:
    """Return the default age and reaction specialist registry."""

    return SpecialistBaselineRegistry(
        (
            local_tracer_age_baseline_spec(),
            atmospheric_history_multitracer_age_baseline_spec(),
            competence_matched_age_baseline_spec(),
            local_thermodynamic_reaction_baseline_spec(),
            stoichiometric_reaction_baseline_spec(),
            competence_matched_reaction_baseline_spec(),
        )
    )


def score_age_baseline_outputs(
    true_ages_years: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]] | None,
) -> dict[str, Any]:
    """Score specialist age intervals after the truth-release boundary."""

    if not predictions:
        return {"status": "not_available", "n": 0}
    observed: list[float] = []
    predicted: list[float] = []
    lower: list[float] = []
    upper: list[float] = []
    for target, truth in true_ages_years.items():
        row = predictions.get(str(target))
        if not isinstance(row, Mapping):
            continue
        values = tuple(_finite(row.get(field)) for field in ("mean_age_years", "age_95_low", "age_95_high"))
        if any(value is None for value in values):
            continue
        estimate, low, high = values
        assert estimate is not None and low is not None and high is not None
        observed.append(float(truth))
        predicted.append(estimate)
        lower.append(low)
        upper.append(high)
    if not observed:
        return {"status": "no_comparable_outputs", "n": 0}
    return {
        "status": "scored",
        "n": len(observed),
        "point": dict(regression_metrics(observed, predicted)),
        "interval": dict(interval_coverage(observed, lower, upper)),
        "n_abstain": sum(
            str(predictions.get(str(target), {}).get("decision", "")) == ABSTAIN
            for target in true_ages_years
        ),
    }


def score_reaction_baseline_outputs(
    true_processes: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]] | None,
    *,
    candidate_edge_ids: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Score fixed reaction-family predictions on a declared edge universe."""

    truth_by_id = {str(edge_id): value for edge_id, value in true_processes.items()}
    candidate_ids = tuple(
        dict.fromkeys(
            str(edge_id)
            for edge_id in (
                candidate_edge_ids
                if candidate_edge_ids is not None
                else truth_by_id.keys()
            )
        )
    )
    if not predictions:
        return {
            "status": "not_available",
            "n": 0,
            "n_truth": len(candidate_ids),
            "n_missing_outputs": len(candidate_ids),
            "outputs_complete": False,
        }
    expected: list[str] = []
    predicted: list[str] = []
    log_losses: list[float] = []
    brier_scores: list[float] = []
    n_missing = 0
    for edge_id in candidate_ids:
        truth = truth_by_id.get(edge_id, "none")
        row = predictions.get(str(edge_id))
        if not isinstance(row, Mapping):
            n_missing += 1
            continue
        expected_family = _reaction_family(str(truth))
        family = str(row.get("family", "none"))
        expected.append(expected_family)
        predicted.append(family)
        probabilities = row.get("probabilities")
        if isinstance(probabilities, Mapping) and probabilities:
            clean_probabilities = {
                str(key): max(0.0, float(value))
                for key, value in probabilities.items()
                if _finite(value) is not None
            }
            total = sum(clean_probabilities.values())
            if total > 0.0:
                clean_probabilities = {
                    key: value / total for key, value in clean_probabilities.items()
                }
                log_losses.append(
                    -math.log(max(clean_probabilities.get(expected_family, 0.0), 1.0e-12))
                )
                brier_scores.append(
                    sum(
                        (
                            clean_probabilities.get(key, 0.0)
                            - float(key == expected_family)
                        )
                        ** 2
                        for key in set(clean_probabilities) | {expected_family}
                    )
                )
    if not expected:
        return {
            "status": "no_comparable_outputs",
            "n": 0,
            "n_truth": len(candidate_ids),
            "n_missing_outputs": n_missing,
            "outputs_complete": n_missing == 0 and bool(candidate_ids),
        }
    metrics = dict(classification_metrics(expected, predicted))
    metrics["unresolved_rate"] = predicted.count("none") / len(predicted)
    metrics["multiclass_log_loss"] = (
        sum(log_losses) / len(log_losses) if log_losses else float("nan")
    )
    metrics["multiclass_brier"] = (
        sum(brier_scores) / len(brier_scores) if brier_scores else float("nan")
    )
    return {
        "status": "scored",
        "n": len(expected),
        "n_truth": len(candidate_ids),
        "n_missing_outputs": n_missing,
        "outputs_complete": n_missing == 0,
        "metrics": metrics,
        "expected_families": sorted(set(expected)),
        "predicted_families": sorted(set(predicted)),
    }


def _predict_local_tracer_age(
    observations: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    channel = observations.get("tracer_age", {})
    nodes = channel.get("nodes", {}) if isinstance(channel, Mapping) else {}
    if not isinstance(nodes, Mapping):
        return {}
    output: dict[str, dict[str, Any]] = {}
    for target, raw in nodes.items():
        features = raw if isinstance(raw, Mapping) else {}
        estimates: list[tuple[float, float, str]] = []
        tritium = _finite(features.get("tritium_TU"))
        if tritium is not None and tritium >= 0.0:
            estimate = _decay_age(tritium, _TRITIUM_MODERN_TU, _TRITIUM_HALF_LIFE_YEARS)
            estimates.append((estimate, max(5.0, 0.20 * estimate), "tritium_TU"))
        argon39 = _finite(features.get("argon39_pmc"))
        if argon39 is not None and argon39 >= 0.0:
            estimate = _decay_age(argon39, _ARGON39_MODERN_PMC, _ARGON39_HALF_LIFE_YEARS)
            estimates.append((estimate, max(8.0, 0.25 * estimate), "argon39_pmc"))
        if not estimates:
            output[str(target)] = {
                "mean_age_years": _AGE_MAX_YEARS / 2.0,
                "age_95_low": 0.0,
                "age_95_high": _AGE_MAX_YEARS,
                "uncertainty_years": _AGE_MAX_YEARS / 2.0,
                "decision": ABSTAIN,
                "reason": "missing_declared_age_tracer",
                "evidence_channel": "tracer_age",
            }
            continue
        weights = [1.0 / (sigma * sigma) for _, sigma, _ in estimates]
        estimate = sum(weight * value for weight, (value, _, _) in zip(weights, estimates)) / sum(weights)
        sigma = math.sqrt(1.0 / sum(weights))
        output[str(target)] = {
            "mean_age_years": _clamp(estimate, 0.0, _AGE_MAX_YEARS),
            "age_95_low": _clamp(estimate - _AGE_Z95 * sigma, 0.0, _AGE_MAX_YEARS),
            "age_95_high": _clamp(estimate + _AGE_Z95 * sigma, 0.0, _AGE_MAX_YEARS),
            "uncertainty_years": sigma,
            "decision": SELECT,
            "reason": "local_physical_decay_conversion",
            "evidence_channel": "tracer_age",
            "tracers_used": [name for _, _, name in estimates],
        }
    return output


def _predict_atmospheric_history_multitracer_age(
    observations: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    channel = observations.get("tracer_age_history", {})
    nodes = channel.get("nodes", {}) if isinstance(channel, Mapping) else {}
    if not isinstance(nodes, Mapping):
        return {}
    output: dict[str, dict[str, Any]] = {}
    for target, raw in nodes.items():
        features = raw if isinstance(raw, Mapping) else {}
        sample_year = _finite(features.get("sample_year")) or 2025.0
        ages = [
            index * _ATMOSPHERIC_AGE_GRID_STEP_YEARS
            for index in range(
                int(_ATMOSPHERIC_AGE_MAX_YEARS / _ATMOSPHERIC_AGE_GRID_STEP_YEARS) + 1
            )
        ]
        log_likelihood: list[float] = []
        tracers_used: set[str] = set()
        for age in ages:
            terms: list[float] = []
            for field_name, history_name in _ATMOSPHERIC_TRACERS.items():
                observed = _finite(features.get(field_name))
                if observed is None or observed < 0.0:
                    continue
                recharge_year = sample_year - age
                history = ATMOSPHERIC_TRACER_HISTORIES.get(history_name)
                if not history:
                    continue
                predicted = _history_value(history, recharge_year)
                sigma = _age_measurement_sigma(features, field_name, observed)
                terms.append(-0.5 * ((observed - predicted) / sigma) ** 2)
                tracers_used.add(field_name)
            tritium = _finite(features.get("tritium_TU"))
            if tritium is not None and tritium >= 0.0:
                predicted = _TRITIUM_MODERN_TU * math.exp(
                    -math.log(2.0) * age / _TRITIUM_HALF_LIFE_YEARS
                )
                sigma = _age_measurement_sigma(features, "tritium_TU", tritium)
                terms.append(-0.5 * ((tritium - predicted) / sigma) ** 2)
                tracers_used.add("tritium_TU")
            argon39 = _finite(features.get("argon39_pmc"))
            if argon39 is not None and argon39 >= 0.0:
                predicted = _ARGON39_MODERN_PMC * math.exp(
                    -math.log(2.0) * age / _ARGON39_HALF_LIFE_YEARS
                )
                sigma = _age_measurement_sigma(features, "argon39_pmc", argon39)
                terms.append(-0.5 * ((argon39 - predicted) / sigma) ** 2)
                tracers_used.add("argon39_pmc")
            if not terms:
                log_likelihood.append(float("nan"))
            else:
                log_likelihood.append(sum(terms))
        if not tracers_used or not any(math.isfinite(value) for value in log_likelihood):
            output[str(target)] = _missing_age_output(
                reason="missing_supported_atmospheric_or_decay_tracer",
                evidence_channel="tracer_age_history",
            )
            continue
        finite_likelihood = [
            value if math.isfinite(value) else -float("inf") for value in log_likelihood
        ]
        maximum = max(finite_likelihood)
        weights = [
            math.exp(value - maximum) if math.isfinite(value) else 0.0
            for value in finite_likelihood
        ]
        total = sum(weights)
        if total <= 0.0 or not math.isfinite(total):
            output[str(target)] = _missing_age_output(
                reason="non_finite_age_likelihood",
                evidence_channel="tracer_age_history",
            )
            continue
        weights = [value / total for value in weights]
        mean = sum(age * weight for age, weight in zip(ages, weights))
        variance = sum((age - mean) ** 2 * weight for age, weight in zip(ages, weights))
        mode_index = max(range(len(weights)), key=weights.__getitem__)
        low = _posterior_quantile(ages, weights, 0.025)
        high = _posterior_quantile(ages, weights, 0.975)
        mode_count = _count_posterior_modes(weights)
        multimodal = mode_count > 1
        output[str(target)] = {
            "mean_age_years": _clamp(mean, 0.0, _ATMOSPHERIC_AGE_MAX_YEARS),
            "age_95_low": _clamp(low, 0.0, _ATMOSPHERIC_AGE_MAX_YEARS),
            "age_95_high": _clamp(high, 0.0, _ATMOSPHERIC_AGE_MAX_YEARS),
            "uncertainty_years": math.sqrt(max(0.0, variance)),
            "mode_age_years": ages[mode_index],
            "posterior_mode_count": mode_count,
            "multimodal": multimodal,
            "decision": ABSTAIN if multimodal else SELECT,
            "reason": (
                "multimodal_atmospheric_history_posterior"
                if multimodal
                else "multitracer_atmospheric_history_likelihood"
            ),
            "evidence_channel": "tracer_age_history",
            "tracers_used": sorted(tracers_used),
            "sample_year": sample_year,
        }
    return output


def _missing_age_output(*, reason: str, evidence_channel: str) -> dict[str, Any]:
    return {
        "mean_age_years": _ATMOSPHERIC_AGE_MAX_YEARS / 2.0,
        "age_95_low": 0.0,
        "age_95_high": _ATMOSPHERIC_AGE_MAX_YEARS,
        "uncertainty_years": _ATMOSPHERIC_AGE_MAX_YEARS / 2.0,
        "decision": ABSTAIN,
        "reason": reason,
        "evidence_channel": evidence_channel,
    }


def _history_value(history: Mapping[str, Any], year: float) -> float:
    years = [float(value) for value in history.get("years", ())]
    values = [float(value) for value in history.get("values", ())]
    if not years or len(years) != len(values):
        return 0.0
    if year <= years[0]:
        return values[0]
    if year >= years[-1]:
        return values[-1]
    for left in range(len(years) - 1):
        if years[left] <= year <= years[left + 1]:
            fraction = (year - years[left]) / (years[left + 1] - years[left])
            return values[left] + fraction * (values[left + 1] - values[left])
    return values[-1]


def _age_measurement_sigma(
    features: Mapping[str, Any],
    field_name: str,
    observed: float,
) -> float:
    sigma_fields = {
        "tritium_TU": "tritium_sigma_TU",
        "argon39_pmc": "argon39_sigma_pmc",
        "cfc11_pptv": "cfc11_sigma_pptv",
        "cfc12_pptv": "cfc12_sigma_pptv",
        "cfc113_pptv": "cfc113_sigma_pptv",
        "sf6_pptv": "sf6_sigma_pptv",
    }
    supplied = _finite(features.get(sigma_fields.get(field_name, "")))
    floor = {
        "tritium_TU": 0.15,
        "argon39_pmc": 2.0,
        "cfc11_pptv": 0.25,
        "cfc12_pptv": 0.25,
        "cfc113_pptv": 0.15,
        "sf6_pptv": 0.05,
    }.get(field_name, 0.05)
    return max(float(supplied or 0.0), floor, 0.08 * max(abs(observed), floor))


def _posterior_quantile(ages: list[float], weights: list[float], quantile: float) -> float:
    pairs = sorted(zip(ages, weights), key=lambda item: item[0])
    cumulative = 0.0
    for age, weight in pairs:
        cumulative += weight
        if cumulative >= quantile:
            return float(age)
    return float(pairs[-1][0])


def _count_posterior_modes(weights: list[float]) -> int:
    if not weights:
        return 0
    maximum = max(weights)
    if maximum <= 0.0:
        return 0
    return sum(
        1
        for index in range(1, len(weights) - 1)
        if weights[index] > weights[index - 1]
        and weights[index] >= weights[index + 1]
        and weights[index] >= 0.05 * maximum
    ) + int(weights[0] >= 0.05 * maximum and weights[0] > weights[1])


def _predict_local_reaction(
    observations: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    channel = observations.get("reaction_chemistry", {})
    edges = channel.get("edges", {}) if isinstance(channel, Mapping) else {}
    if not isinstance(edges, Mapping):
        return {}
    output: dict[str, dict[str, Any]] = {}
    for edge_id, raw in edges.items():
        features = raw if isinstance(raw, Mapping) else {}
        upstream = features.get("upstream", {})
        downstream = features.get("downstream", {})
        if not isinstance(upstream, Mapping) or not isinstance(downstream, Mapping):
            output[str(edge_id)] = _missing_reaction_output()
            continue
        delta = {
            key: (_finite(downstream.get(key)) - _finite(upstream.get(key)))
            if _finite(downstream.get(key)) is not None and _finite(upstream.get(key)) is not None
            else None
            for key in ("Ca", "Mg", "Na", "HCO3", "SO4", "NO3", "Fe", "pH")
        }
        if not any(value is not None for value in delta.values()):
            output[str(edge_id)] = _missing_reaction_output()
            continue
        scores = {family: 0.0 for family in _REACTION_FAMILIES}
        scores["other"] = 0.25
        if _negative(delta.get("NO3"), 0.05):
            scores["denitrification"] += 4.0
        if _negative(delta.get("SO4"), 0.05):
            scores["sulfate_reduction"] += 4.0
        if _positive(delta.get("Fe"), 0.01):
            scores["iron_reduction"] += 3.5
        if _positive(delta.get("Ca"), 0.05) and _positive(delta.get("HCO3"), 0.05):
            scores["carbonate"] += 4.0
        if _positive(delta.get("Na"), 0.05) and _negative(delta.get("Mg"), 0.05):
            scores["silicate_exchange"] += 4.0
        if _positive(delta.get("SO4"), 0.05):
            scores["sulfate_source"] += 3.0
        if _negative(delta.get("pH"), 0.10):
            scores["other_redox"] += 1.0
        probabilities = _softmax(scores)
        family = max(probabilities, key=probabilities.get)
        probability = probabilities[family]
        output[str(edge_id)] = {
            "family": family,
            "probabilities": probabilities,
            "logits": scores,
            "raw_scores": scores,
            "probability": probability,
            "decision": SELECT if probability >= 0.60 else ABSTAIN,
            "reason": "local_concentration_difference_rules",
            "evidence_channel": "reaction_chemistry",
        }
    return output


_STOICHIOMETRIC_TEMPLATES: dict[str, dict[str, float]] = {
    "carbonate": {"Ca": 1.0, "HCO3": 1.0},
    "silicate_exchange": {"Na": 1.0, "Mg": -1.0, "Ca": -0.5},
    "sulfate_reduction": {"SO4": -1.0, "HCO3": 1.0},
    "iron_reduction": {"Fe": 1.0, "HCO3": 1.0},
    "denitrification": {"NO3": -1.0, "HCO3": 1.0},
    "sulfate_source": {"SO4": 1.0},
    "other_redox": {"pH": -1.0},
}
_STOICHIOMETRIC_SCALES = {
    "Ca": 0.05,
    "Mg": 0.05,
    "Na": 0.05,
    "HCO3": 0.05,
    "SO4": 0.05,
    "NO3": 0.05,
    "Fe": 0.01,
    "pH": 0.10,
}


def _predict_stoichiometric_reaction(
    observations: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    channel = observations.get("reaction_chemistry", {})
    edges = channel.get("edges", {}) if isinstance(channel, Mapping) else {}
    if not isinstance(edges, Mapping):
        return {}
    output: dict[str, dict[str, Any]] = {}
    for edge_id, raw in edges.items():
        features = raw if isinstance(raw, Mapping) else {}
        upstream = features.get("upstream", {})
        downstream = features.get("downstream", {})
        if not isinstance(upstream, Mapping) or not isinstance(downstream, Mapping):
            output[str(edge_id)] = _missing_reaction_output(
                reason="invalid_paired_chemistry"
            )
            continue
        delta: dict[str, float] = {}
        for ion, scale in _STOICHIOMETRIC_SCALES.items():
            source = _finite(upstream.get(ion))
            target = _finite(downstream.get(ion))
            if source is not None and target is not None:
                delta[ion] = (target - source) / scale
        if not delta:
            output[str(edge_id)] = _missing_reaction_output(
                reason="missing_paired_chemistry"
            )
            continue
        scores = {family: -2.0 for family in _REACTION_FAMILIES}
        extents: dict[str, float] = {}
        residuals: dict[str, float] = {}
        fields_used: dict[str, list[str]] = {}
        for family, template in _STOICHIOMETRIC_TEMPLATES.items():
            available = [ion for ion in template if ion in delta]
            if not available:
                continue
            template_values = [float(template[ion]) for ion in available]
            observed_values = [float(delta[ion]) for ion in available]
            denominator = sum(value * value for value in template_values)
            extent = max(
                0.0,
                sum(observed * basis for observed, basis in zip(observed_values, template_values))
                / denominator,
            ) if denominator > 0.0 else 0.0
            residual = math.sqrt(
                sum(
                    (observed - extent * basis) ** 2
                    for observed, basis in zip(observed_values, template_values)
                )
                / len(available)
            )
            extents[family] = extent
            residuals[family] = residual
            fields_used[family] = available
            # A reaction explanation pays a fixed complexity cost. This
            # makes the explicit zero-change/null model win when paired
            # chemistry is unchanged, rather than assigning a reaction family
            # merely because every reaction template has zero residual.
            scores[family] = (
                -residual
                + 0.15 * math.log1p(extent)
                - 0.25
            )
        delta_norm = math.sqrt(
            sum(value * value for value in delta.values()) / len(delta)
        )
        scores["none"] = -delta_norm
        scores["other"] = -0.75
        probabilities = _softmax(scores)
        family = max(probabilities, key=probabilities.get)
        probability = probabilities[family]
        output[str(edge_id)] = {
            "family": family,
            "probabilities": probabilities,
            "logits": scores,
            "raw_scores": scores,
            "probability": probability,
            "decision": SELECT if probability >= 0.60 else ABSTAIN,
            "reason": "nonnegative_stoichiometric_template_projection",
            "evidence_channel": "reaction_chemistry",
            "extent": extents.get(family, 0.0),
            "residual_norm": residuals.get(family),
            "template_fields_used": fields_used.get(family, []),
        }
    return output


def _missing_reaction_output(*, reason: str = "missing_paired_chemistry") -> dict[str, Any]:
    probability = 1.0 / len(_REACTION_FAMILIES)
    return {
        "family": "none",
        "probabilities": {family: probability for family in _REACTION_FAMILIES},
        "probability": probability,
        "decision": ABSTAIN,
        "reason": reason,
        "evidence_channel": "reaction_chemistry",
    }


def _decay_age(value: float, modern: float, half_life: float) -> float:
    ratio = _clamp(value / modern, 1.0e-6, 1.0)
    age = -math.log(ratio) * half_life / math.log(2.0)
    return _clamp(age, 0.0, _AGE_MAX_YEARS)


def _positive(value: float | None, threshold: float) -> bool:
    return value is not None and value >= threshold


def _negative(value: float | None, threshold: float) -> bool:
    return value is not None and value <= -threshold


def _softmax(scores: Mapping[str, float]) -> dict[str, float]:
    maximum = max(scores.values())
    exponentials = {key: math.exp(value - maximum) for key, value in scores.items()}
    total = sum(exponentials.values())
    return {key: value / total for key, value in exponentials.items()}


def _reaction_family(label: str) -> str:
    text = str(label or "").lower()
    if text in {"", "none", "no_reaction", "null", "missing"}:
        return "none"
    if any(token in text for token in ("calcite", "dolomite", "carbonate")):
        return "carbonate"
    if any(token in text for token in ("albite", "anorthite", "feldspar", "silicate", "exch")):
        return "silicate_exchange"
    if "sulfate_reduction" in text:
        return "sulfate_reduction"
    if "iron_reduction" in text:
        return "iron_reduction"
    if "denit" in text:
        return "denitrification"
    if "pyrite" in text:
        return "other_redox"
    if "gypsum" in text or "so4" in text:
        return "sulfate_source"
    return "other"


def normalise_reaction_family(label: str) -> str:
    """Map generator/process labels to the specialist scoring taxonomy."""

    return _reaction_family(label)


__all__ = [
    "ABSTAIN",
    "AGE_OUTPUT",
    "REACTION_OUTPUT",
    "SELECT",
    "SpecialistBaselineRegistry",
    "SpecialistBaselineSpec",
    "default_specialist_baseline_registry",
    "competence_matched_age_baseline_spec",
    "competence_matched_reaction_baseline_spec",
    "atmospheric_history_multitracer_age_baseline_spec",
    "local_thermodynamic_reaction_baseline_spec",
    "local_tracer_age_baseline_spec",
    "normalise_reaction_family",
    "score_age_baseline_outputs",
    "score_reaction_baseline_outputs",
    "stoichiometric_reaction_baseline_spec",
]
