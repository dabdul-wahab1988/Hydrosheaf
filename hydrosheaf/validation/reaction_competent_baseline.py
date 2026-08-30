"""Independent, truth-blind reaction-family competence baseline.

This module is deliberately small and self-contained.  It treats the observed
change between the two ends of a reaction-chemistry edge as a vector and
compares that vector with fixed, signed reaction-family signatures.  The
signatures are dimensionless observation-space heuristics, not balanced
geochemical equations and not outputs from PHREEQC or HydroSheaf inference.

The baseline is intended for validation and model-disagreement reporting.  It
does not estimate a posterior, fit a HydroSheaf model, or claim that a chosen
family is chemically identified.  Its probabilities are deterministic
softmax rankings and are explicitly marked as uncalibrated.  A non-reaction
candidate is always considered when enough paired chemistry is available and
is assigned an explicit complexity penalty (zero by default, while reaction
families carry a positive penalty).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

REACTION_OUTPUT: Final[str] = "reaction_family"
ABSTAIN: Final[str] = "abstain"
SELECT: Final[str] = "select"
NULL_FAMILY: Final[str] = "none"
EVIDENCE_CHANNEL: Final[str] = "reaction_chemistry"

_DEFAULT_ION_SCALES: Final[dict[str, float]] = {
    "Ca": 0.05,
    "Mg": 0.05,
    "Na": 0.05,
    "HCO3": 0.05,
    "SO4": 0.05,
    "NO3": 0.05,
    "Fe": 0.01,
    "pH": 0.10,
}

_CONCENTRATION_FIELDS: Final[frozenset[str]] = frozenset(
    {"Ca", "Mg", "Na", "HCO3", "SO4", "NO3", "Fe"}
)

_FORBIDDEN_TRUTH_KEYS: Final[frozenset[str]] = frozenset(
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


def _finite_number(value: Any) -> float | None:
    """Return a finite float, rejecting booleans and non-numeric values."""

    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _walk_forbidden_keys(value: Any, path: tuple[str, ...] = ()) -> None:
    """Reject labels/reference fields before any candidate is generated."""

    if isinstance(value, Mapping):
        for raw_key, child in value.items():
            key = str(raw_key).strip().lower()
            if key in _FORBIDDEN_TRUTH_KEYS or key.startswith(("truth_", "true_")):
                location = ".".join(path + (str(raw_key),))
                raise ValueError(
                    f"Truth/reference field is forbidden for this baseline: {location}"
                )
            _walk_forbidden_keys(child, path + (str(raw_key),))
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for index, child in enumerate(value):
            _walk_forbidden_keys(child, path + (str(index),))


@dataclass(frozen=True)
class ReactionTemplate:
    """A fixed normalized observation-space signature for one reaction family.

    ``stoichiometry`` accepts either a mapping or an iterable of pairs.  The
    values are signed coefficients in normalized observation space.  They are
    intentionally not interpreted as molar stoichiometric coefficients.
    """

    family: str
    stoichiometry: Mapping[str, float] | Iterable[tuple[str, float]]
    description: str = ""
    complexity_penalty: float = 0.25
    minimum_fields: int | None = None

    def __post_init__(self) -> None:
        family = str(self.family).strip()
        if not family or family == NULL_FAMILY:
            raise ValueError("ReactionTemplate.family must be non-empty and not 'none'")

        items = (
            self.stoichiometry.items()
            if isinstance(self.stoichiometry, Mapping)
            else self.stoichiometry
        )
        canonical: list[tuple[str, float]] = []
        seen: set[str] = set()
        for raw_field, raw_coefficient in items:
            field_name = str(raw_field).strip()
            coefficient = _finite_number(raw_coefficient)
            if not field_name:
                raise ValueError("ReactionTemplate fields must be non-empty")
            if field_name in seen:
                raise ValueError(f"Duplicate reaction field: {field_name}")
            if coefficient is None or coefficient == 0.0:
                raise ValueError(
                    f"Reaction coefficient for {field_name!r} must be finite/non-zero"
                )
            seen.add(field_name)
            canonical.append((field_name, coefficient))

        canonical.sort(key=lambda item: item[0])
        if not canonical:
            raise ValueError(
                "ReactionTemplate.stoichiometry must contain a coefficient"
            )

        penalty = _finite_number(self.complexity_penalty)
        if penalty is None or penalty < 0.0:
            raise ValueError(
                "ReactionTemplate.complexity_penalty must be finite and non-negative"
            )

        minimum_fields = self.minimum_fields
        if minimum_fields is None:
            minimum_fields = 2 if len(canonical) > 1 else 1
        if isinstance(minimum_fields, bool) or not isinstance(minimum_fields, int):
            raise TypeError("ReactionTemplate.minimum_fields must be an integer")
        if minimum_fields < 1 or minimum_fields > len(canonical):
            raise ValueError("ReactionTemplate.minimum_fields must fit its signature")

        object.__setattr__(self, "family", family)
        object.__setattr__(self, "stoichiometry", tuple(canonical))
        object.__setattr__(self, "complexity_penalty", penalty)
        object.__setattr__(self, "minimum_fields", minimum_fields)

    @property
    def fields(self) -> tuple[str, ...]:
        return tuple(field_name for field_name, _ in self.stoichiometry)

    def to_dict(self) -> dict[str, Any]:
        return {
            "family": self.family,
            "stoichiometry": {
                name: coefficient for name, coefficient in self.stoichiometry
            },
            "description": self.description,
            "complexity_penalty": self.complexity_penalty,
            "minimum_fields": self.minimum_fields,
        }


DEFAULT_REACTION_TEMPLATES: Final[tuple[ReactionTemplate, ...]] = (
    ReactionTemplate(
        "carbonate",
        {"Ca": 1.0, "HCO3": 1.0},
        "co-movement of calcium and bicarbonate",
    ),
    ReactionTemplate(
        "silicate_exchange",
        {"Na": 1.0, "Mg": -1.0, "Ca": -0.5},
        "sodium gain opposed by divalent-cation loss",
    ),
    ReactionTemplate(
        "sulfate_reduction",
        {"SO4": -1.0, "HCO3": 1.0},
        "sulfate loss with bicarbonate gain",
    ),
    ReactionTemplate(
        "iron_reduction",
        {"Fe": 1.0, "HCO3": 1.0},
        "iron gain with bicarbonate gain",
    ),
    ReactionTemplate(
        "denitrification",
        {"NO3": -1.0, "HCO3": 1.0},
        "nitrate loss with bicarbonate gain",
    ),
    ReactionTemplate(
        "sulfate_source",
        {"SO4": 1.0},
        "sulfate-only source signal",
        minimum_fields=1,
    ),
    ReactionTemplate(
        "other_redox",
        {"pH": -1.0},
        "pH-only redox signal",
        minimum_fields=1,
    ),
)


@dataclass(frozen=True)
class ReactionBaselineConfig:
    """Deterministic thresholds and fixed conventions for the baseline."""

    templates: tuple[ReactionTemplate, ...] = field(
        default_factory=lambda: DEFAULT_REACTION_TEMPLATES
    )
    ion_scales: Mapping[str, float] = field(
        default_factory=lambda: dict(_DEFAULT_ION_SCALES)
    )
    null_complexity_penalty: float = 0.0
    probability_temperature: float = 0.50
    selection_threshold: float = 0.60
    selection_margin: float = 0.10
    null_tolerance: float = 0.10
    min_paired_fields: int = 2
    max_residual_norm: float = 2.50
    confidence_z: float = 1.959963984540054

    def __post_init__(self) -> None:
        templates = tuple(self.templates)
        families = [template.family for template in templates]
        if len(families) != len(set(families)):
            raise ValueError("Reaction template families must be unique")
        if NULL_FAMILY in families:
            raise ValueError(
                "The null family is supplied by the baseline, not a template"
            )

        canonical_scales: dict[str, float] = {}
        for raw_name, raw_scale in self.ion_scales.items():
            name = str(raw_name).strip()
            scale = _finite_number(raw_scale)
            if not name or scale is None or scale <= 0.0:
                raise ValueError("ion_scales must contain finite positive values")
            canonical_scales[name] = scale

        penalty = _finite_number(self.null_complexity_penalty)
        temperature = _finite_number(self.probability_temperature)
        threshold = _finite_number(self.selection_threshold)
        margin = _finite_number(self.selection_margin)
        tolerance = _finite_number(self.null_tolerance)
        max_residual = _finite_number(self.max_residual_norm)
        confidence_z = _finite_number(self.confidence_z)
        if penalty is None or penalty < 0.0:
            raise ValueError("null_complexity_penalty must be finite and non-negative")
        if temperature is None or temperature <= 0.0:
            raise ValueError("probability_temperature must be positive")
        if threshold is None or not 0.0 <= threshold <= 1.0:
            raise ValueError("selection_threshold must lie in [0, 1]")
        if margin is None or not 0.0 <= margin <= 1.0:
            raise ValueError("selection_margin must lie in [0, 1]")
        if tolerance is None or tolerance < 0.0:
            raise ValueError("null_tolerance must be non-negative")
        if max_residual is None or max_residual < 0.0:
            raise ValueError("max_residual_norm must be non-negative")
        if confidence_z is None or confidence_z < 0.0:
            raise ValueError("confidence_z must be non-negative")
        if (
            isinstance(self.min_paired_fields, bool)
            or not isinstance(self.min_paired_fields, int)
            or self.min_paired_fields < 1
        ):
            raise ValueError("min_paired_fields must be a positive integer")

        object.__setattr__(self, "templates", templates)
        object.__setattr__(self, "ion_scales", dict(sorted(canonical_scales.items())))
        object.__setattr__(self, "null_complexity_penalty", penalty)
        object.__setattr__(self, "probability_temperature", temperature)
        object.__setattr__(self, "selection_threshold", threshold)
        object.__setattr__(self, "selection_margin", margin)
        object.__setattr__(self, "null_tolerance", tolerance)
        object.__setattr__(self, "max_residual_norm", max_residual)
        object.__setattr__(self, "confidence_z", confidence_z)


@dataclass(frozen=True)
class ReactionCandidate:
    """One reaction-family explanation for one observed edge."""

    family: str
    extent: float
    residual_norm: float
    penalized_loss: float
    complexity_penalty: float
    fields_used: tuple[str, ...]
    residuals: Mapping[str, float]
    extent_std: float | None = None
    extent_ci95: tuple[float, float] | None = None
    description: str = ""

    @property
    def score(self) -> float:
        return -self.penalized_loss

    @property
    def is_null(self) -> bool:
        return self.family == NULL_FAMILY

    def to_dict(self) -> dict[str, Any]:
        uncertainty: dict[str, Any] = {
            "extent_std": self.extent_std,
            "extent_ci95": list(self.extent_ci95)
            if self.extent_ci95 is not None
            else None,
            "calibrated": False,
            "method": "residual_linearization",
        }
        return {
            "family": self.family,
            "extent": self.extent,
            "residual_norm": self.residual_norm,
            "penalized_loss": self.penalized_loss,
            "complexity_penalty": self.complexity_penalty,
            "score": self.score,
            "fields_used": list(self.fields_used),
            "residuals": dict(self.residuals),
            "uncertainty": uncertainty,
            "description": self.description,
            "is_null": self.is_null,
        }


@dataclass(frozen=True)
class _EdgeData:
    paired_values: tuple[tuple[str, float], ...]
    normalized_delta: tuple[tuple[str, float], ...]
    missing_fields: tuple[str, ...]
    unsupported_fields: tuple[str, ...]
    invalid_fields: tuple[str, ...]
    status: str | None = None

    @property
    def observed_delta_norm(self) -> float | None:
        if not self.normalized_delta:
            return None
        return math.sqrt(
            sum(delta * delta for _, delta in self.normalized_delta)
            / len(self.normalized_delta)
        )


def _edge_items(observations: Mapping[str, Any]) -> list[tuple[str, Mapping[str, Any]]]:
    channel = observations.get(EVIDENCE_CHANNEL)
    if not isinstance(channel, Mapping):
        return []
    raw_edges = channel.get("edges")
    if isinstance(raw_edges, Mapping):
        result: list[tuple[str, Mapping[str, Any]]] = []
        for raw_id, raw_edge in raw_edges.items():
            if isinstance(raw_edge, Mapping):
                result.append((str(raw_id), raw_edge))
        return result
    if isinstance(raw_edges, Sequence) and not isinstance(
        raw_edges, (str, bytes, bytearray)
    ):
        result = []
        for index, raw_edge in enumerate(raw_edges):
            if not isinstance(raw_edge, Mapping):
                continue
            raw_id = raw_edge.get("edge_id", raw_edge.get("id", f"edge_{index}"))
            result.append((str(raw_id), raw_edge))
        return result
    return []


def _side(
    edge: Mapping[str, Any], primary: str, fallback: str
) -> Mapping[str, Any] | None:
    value = edge.get(primary, edge.get(fallback))
    return value if isinstance(value, Mapping) else None


def _extract_edge_data(
    edge: Mapping[str, Any], config: ReactionBaselineConfig
) -> _EdgeData:
    upstream = _side(edge, "upstream", "source")
    downstream = _side(edge, "downstream", "target_values")
    if downstream is None:
        downstream = _side(edge, "downstream", "destination")
    if upstream is None or downstream is None:
        return _EdgeData((), (), (), (), (), "malformed_edge")

    paired: list[tuple[str, float]] = []
    normalized: list[tuple[str, float]] = []
    missing: list[str] = []
    invalid: list[str] = []
    for field_name, scale in config.ion_scales.items():
        if field_name not in upstream or field_name not in downstream:
            missing.append(field_name)
            continue
        before = _finite_number(upstream[field_name])
        after = _finite_number(downstream[field_name])
        if before is None or after is None:
            invalid.append(field_name)
            continue
        if field_name in _CONCENTRATION_FIELDS and (before < 0.0 or after < 0.0):
            invalid.append(field_name)
            continue
        paired.append((field_name, after - before))
        normalized.append((field_name, (after - before) / scale))

    supported_names = set(config.ion_scales)
    supplied_names = {str(name) for name in upstream} | {
        str(name) for name in downstream
    }
    unsupported = sorted(supplied_names - supported_names)
    status: str | None = None
    if invalid:
        status = "invalid_chemistry_values"
    elif not normalized:
        status = (
            "unsupported_chemistry_fields" if unsupported else "missing_paired_ions"
        )
    elif len(normalized) < config.min_paired_fields:
        status = "insufficient_paired_ions"
    return _EdgeData(
        tuple(paired),
        tuple(normalized),
        tuple(sorted(missing)),
        tuple(unsupported),
        tuple(sorted(invalid)),
        status,
    )


def _candidate_for_template(
    template: ReactionTemplate,
    edge_data: _EdgeData,
    confidence_z: float,
) -> ReactionCandidate | None:
    delta = dict(edge_data.normalized_delta)
    fields = tuple(field_name for field_name in template.fields if field_name in delta)
    if len(fields) < int(template.minimum_fields):
        return None
    coefficients = {
        field_name: dict(template.stoichiometry)[field_name] for field_name in fields
    }
    denominator = sum(coefficients[field_name] ** 2 for field_name in fields)
    if denominator <= 0.0:
        return None
    extent = max(
        0.0,
        sum(coefficients[field_name] * delta[field_name] for field_name in fields)
        / denominator,
    )
    residuals = {
        field_name: delta[field_name] - extent * coefficients[field_name]
        for field_name in fields
    }
    residual_norm = math.sqrt(
        sum(residual * residual for residual in residuals.values()) / len(residuals)
    )
    sse = sum(residual * residual for residual in residuals.values())
    extent_std: float | None = None
    extent_ci95: tuple[float, float] | None = None
    if len(fields) >= 2:
        residual_variance = sse / max(len(fields) - 1, 1)
        extent_std = math.sqrt(max(0.0, residual_variance / denominator))
        extent_ci95 = (
            max(0.0, extent - confidence_z * extent_std),
            extent + confidence_z * extent_std,
        )
    return ReactionCandidate(
        family=template.family,
        extent=extent,
        residual_norm=residual_norm,
        penalized_loss=residual_norm + template.complexity_penalty,
        complexity_penalty=template.complexity_penalty,
        fields_used=fields,
        residuals=residuals,
        extent_std=extent_std,
        extent_ci95=extent_ci95,
        description=template.description,
    )


def _null_candidate(
    edge_data: _EdgeData, config: ReactionBaselineConfig
) -> ReactionCandidate:
    residuals = dict(edge_data.normalized_delta)
    residual_norm = math.sqrt(
        sum(residual * residual for residual in residuals.values()) / len(residuals)
    )
    return ReactionCandidate(
        family=NULL_FAMILY,
        extent=0.0,
        residual_norm=residual_norm,
        penalized_loss=residual_norm + config.null_complexity_penalty,
        complexity_penalty=config.null_complexity_penalty,
        fields_used=tuple(name for name, _ in edge_data.normalized_delta),
        residuals=residuals,
        description="no reaction change",
    )


def _softmax(
    candidates: Sequence[ReactionCandidate], temperature: float
) -> dict[str, float]:
    if not candidates:
        return {}
    logits = [-candidate.penalized_loss / temperature for candidate in candidates]
    maximum = max(logits)
    weights = [math.exp(logit - maximum) for logit in logits]
    total = sum(weights)
    return {
        candidate.family: weight / total
        for candidate, weight in zip(candidates, weights, strict=True)
    }


class ReactionCompetentBaseline:
    """Independent reaction-family comparator for paired chemistry edges."""

    name: Final[str] = "reaction_competent_baseline"
    version: Final[str] = "0.1.0"

    def __init__(self, config: ReactionBaselineConfig | None = None) -> None:
        self.config = config or ReactionBaselineConfig()

    def generate_candidates(
        self, observations: Mapping[str, Any]
    ) -> dict[str, tuple[ReactionCandidate, ...]]:
        """Generate null and reaction-family candidates without selecting one."""

        _walk_forbidden_keys(observations)
        result: dict[str, tuple[ReactionCandidate, ...]] = {}
        for edge_id, edge in _edge_items(observations):
            edge_data = _extract_edge_data(edge, self.config)
            if edge_data.status is not None:
                continue
            candidates: list[ReactionCandidate] = [
                _null_candidate(edge_data, self.config)
            ]
            for template in self.config.templates:
                candidate = _candidate_for_template(
                    template, edge_data, self.config.confidence_z
                )
                if candidate is not None:
                    candidates.append(candidate)
            result[edge_id] = tuple(candidates)
        return result

    def predict(self, observations: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
        """Return selected explanations plus auditable uncertainty diagnostics."""

        _walk_forbidden_keys(observations)
        result: dict[str, dict[str, Any]] = {}
        channel_present = isinstance(observations.get(EVIDENCE_CHANNEL), Mapping)
        for edge_id, edge in _edge_items(observations):
            edge_data = _extract_edge_data(edge, self.config)
            if edge_data.status is not None:
                result[edge_id] = self._abstention_record(edge_id, edge_data)
                continue

            candidates = [_null_candidate(edge_data, self.config)]
            for template in self.config.templates:
                candidate = _candidate_for_template(
                    template, edge_data, self.config.confidence_z
                )
                if candidate is not None:
                    candidates.append(candidate)
            candidates.sort(
                key=lambda candidate: (candidate.penalized_loss, candidate.family)
            )
            probabilities = _softmax(candidates, self.config.probability_temperature)
            best = candidates[0]
            second_probability = (
                probabilities[candidates[1].family] if len(candidates) > 1 else 0.0
            )
            observed_norm = edge_data.observed_delta_norm

            if (
                observed_norm is not None
                and observed_norm <= self.config.null_tolerance
            ):
                probabilities = {candidate.family: 0.0 for candidate in candidates}
                probabilities[NULL_FAMILY] = 1.0
                best = next(candidate for candidate in candidates if candidate.is_null)
                decision = SELECT
                reason = "no_reaction_supported_by_delta_tolerance"
            elif best.residual_norm > self.config.max_residual_norm:
                decision = ABSTAIN
                reason = "high_residual"
            elif (
                probabilities[best.family] < self.config.selection_threshold
                or probabilities[best.family] - second_probability
                < self.config.selection_margin
            ):
                decision = ABSTAIN
                reason = "ambiguous_reaction_explanation"
            else:
                decision = SELECT
                reason = (
                    "null_candidate_selected_after_complexity_penalty"
                    if best.is_null
                    else "reaction_family_selected"
                )

            accepted_family = best.family if decision == SELECT else NULL_FAMILY
            accepted_extent = best.extent if decision == SELECT else None
            accepted_uncertainty = (
                {
                    "extent_std": best.extent_std,
                    "extent_ci95": (
                        list(best.extent_ci95) if best.extent_ci95 is not None else None
                    ),
                    "calibrated": False,
                    "selection_probability_calibrated": False,
                }
                if decision == SELECT
                else {
                    "extent_std": None,
                    "extent_ci95": None,
                    "calibrated": False,
                    "selection_probability_calibrated": False,
                }
            )
            candidate_records = []
            for candidate in candidates:
                record = candidate.to_dict()
                record["probability"] = probabilities[candidate.family]
                candidate_records.append(record)
            result[edge_id] = {
                "edge_id": edge_id,
                "family": accepted_family,
                "candidate_family": best.family,
                "probability": probabilities[best.family],
                "probabilities": probabilities,
                "decision": decision,
                "reason": reason,
                "evidence_channel": EVIDENCE_CHANNEL,
                "extent": accepted_extent,
                "residual_norm": best.residual_norm,
                "uncertainty": accepted_uncertainty,
                "diagnostics": {
                    "observed_delta_norm": observed_norm,
                    "paired_fields": [name for name, _ in edge_data.normalized_delta],
                    "n_paired_fields": len(edge_data.normalized_delta),
                    "candidate_count": len(candidates),
                    "null_tolerance": self.config.null_tolerance,
                    "null_complexity_penalty": self.config.null_complexity_penalty,
                    "reaction_complexity_penalties": {
                        candidate.family: candidate.complexity_penalty
                        for candidate in candidates
                        if not candidate.is_null
                    },
                    "selection_probability_is_uncalibrated": True,
                },
                "candidate_explanations": candidate_records,
            }

        if not channel_present:
            return {}
        return result

    def _abstention_record(self, edge_id: str, edge_data: _EdgeData) -> dict[str, Any]:
        reason = edge_data.status or "unsupported_reaction_chemistry"
        return {
            "edge_id": edge_id,
            "family": NULL_FAMILY,
            "candidate_family": None,
            "probability": 0.0,
            "probabilities": {NULL_FAMILY: 1.0},
            "decision": ABSTAIN,
            "reason": reason,
            "evidence_channel": EVIDENCE_CHANNEL,
            "extent": None,
            "residual_norm": None,
            "uncertainty": {
                "extent_std": None,
                "extent_ci95": None,
                "calibrated": False,
                "selection_probability_calibrated": False,
            },
            "diagnostics": {
                "paired_fields": [name for name, _ in edge_data.normalized_delta],
                "n_paired_fields": len(edge_data.normalized_delta),
                "missing_fields": list(edge_data.missing_fields),
                "unsupported_fields": list(edge_data.unsupported_fields),
                "invalid_fields": list(edge_data.invalid_fields),
                "selection_probability_is_uncalibrated": True,
            },
            "candidate_explanations": [],
        }

    def to_audit_record(self) -> dict[str, Any]:
        """Describe provenance and limitations for a validation audit."""

        return {
            "name": self.name,
            "version": self.version,
            "family": "reaction_family",
            "output_kind": REACTION_OUTPUT,
            "input_channels": [EVIDENCE_CHANNEL],
            "tuning": {
                "templates": [template.to_dict() for template in self.config.templates],
                "ion_scales": dict(self.config.ion_scales),
                "probability_temperature": self.config.probability_temperature,
                "selection_threshold": self.config.selection_threshold,
                "selection_margin": self.config.selection_margin,
                "null_tolerance": self.config.null_tolerance,
                "confidence_z": self.config.confidence_z,
            },
            "uncertainty": {
                "method": "residual_linearization",
                "calibrated": False,
                "selection_probability_calibrated": False,
                "abstain_when_unsupported_or_ambiguous": True,
            },
            "abstention": {
                "unsupported_input": True,
                "invalid_input": True,
                "insufficient_paired_fields": True,
                "high_residual": True,
                "ambiguous_candidate_set": True,
            },
            "cost": {"complexity_penalty_is_explicit": True},
            "control": {
                "truth_blind": True,
                "uses_hydrosheaf_inference": False,
                "uses_hydrosheaf_posterior": False,
                "uses_phreeqc_outputs_as_truth": False,
                "signature_semantics": "fixed normalized observation-space heuristics",
            },
        }


def reaction_competent_baseline(
    config: ReactionBaselineConfig | None = None,
) -> ReactionCompetentBaseline:
    """Construct the independent baseline with optional deterministic config."""

    return ReactionCompetentBaseline(config=config)


def generate_reaction_candidates(
    observations: Mapping[str, Any],
    *,
    config: ReactionBaselineConfig | None = None,
) -> dict[str, tuple[ReactionCandidate, ...]]:
    """Functional wrapper for candidate generation."""

    return reaction_competent_baseline(config).generate_candidates(observations)


def predict_reaction_families(
    observations: Mapping[str, Any],
    *,
    config: ReactionBaselineConfig | None = None,
) -> dict[str, dict[str, Any]]:
    """Functional wrapper for validation predictions and diagnostics."""

    return reaction_competent_baseline(config).predict(observations)


__all__ = [
    "ABSTAIN",
    "DEFAULT_REACTION_TEMPLATES",
    "EVIDENCE_CHANNEL",
    "NULL_FAMILY",
    "REACTION_OUTPUT",
    "SELECT",
    "ReactionBaselineConfig",
    "ReactionCandidate",
    "ReactionCompetentBaseline",
    "ReactionTemplate",
    "generate_reaction_candidates",
    "predict_reaction_families",
    "reaction_competent_baseline",
]


