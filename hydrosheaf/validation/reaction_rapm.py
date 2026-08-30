"""Development-fitted RAPM-style reaction-family evidence model.

This module adds a regularized, case-blocked reaction scorer without replacing
the fixed-template comparator.  Reaction-family labels are used only by
``fit`` on development cases.  ``predict`` accepts observed paired chemistry
only and therefore remains truth-blind for locked evaluation.

The model is deliberately modest: it is a multiclass ridge evidence model
with a separate on/off family-effect diagnostic.  It is not a thermodynamic
equilibrium solver, a PHREEQC inverse, or a causal estimator.  Its purpose is
to learn stable, adjusted evidence weights from development generators while
retaining an explicit null/mixing control and conservative abstention.  The
module also provides a case-held-out temperature layer and a truth-release
performance report; neither turns the model into a mechanistic inverse.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
import hashlib
import json
import math
from typing import Any, Final

import numpy as np


REACTION_OUTPUT: Final[str] = "reaction_family"
ABSTAIN: Final[str] = "abstain"
SELECT: Final[str] = "select"
NULL_FAMILY: Final[str] = "none"
EVIDENCE_CHANNEL: Final[str] = "reaction_chemistry"

REACTION_FIELDS: Final[tuple[str, ...]] = (
    "Ca",
    "Cl",
    "F",
    "Fe",
    "HCO3",
    "K",
    "Mg",
    "Na",
    "NO3",
    "PO4",
    "SO4",
    "SiO2",
    "pH",
)

DEFAULT_REACTION_CLASSES: Final[tuple[str, ...]] = (
    "carbonate",
    "silicate_exchange",
    "sulfate_reduction",
    "iron_reduction",
    "denitrification",
    "sulfate_source",
    "other_redox",
    "other",
    NULL_FAMILY,
    "mixing",
)

TEMPLATE_EVIDENCE_FAMILIES: Final[tuple[str, ...]] = (
    "carbonate",
    "silicate_exchange",
    "sulfate_reduction",
    "iron_reduction",
    "denitrification",
    "sulfate_source",
    "other_redox",
    NULL_FAMILY,
    "mixing",
)

TEMPLATE_EVIDENCE_LIBRARY: Final[dict[str, dict[str, float]]] = {
    "carbonate": {"Ca": 1.0, "HCO3": 1.0},
    "silicate_exchange": {"Na": 1.0, "Mg": -1.0, "Ca": -0.5},
    "sulfate_reduction": {"SO4": -1.0, "HCO3": 1.0},
    "iron_reduction": {"Fe": 1.0, "HCO3": 1.0},
    "denitrification": {"NO3": -1.0, "HCO3": 1.0},
    "sulfate_source": {"SO4": 1.0},
    "other_redox": {"pH": -1.0},
}

DEFAULT_ION_SCALES: Final[dict[str, float]] = {
    "Ca": 0.05,
    "Cl": 0.05,
    "F": 0.01,
    "Fe": 0.01,
    "HCO3": 0.05,
    "K": 0.05,
    "Mg": 0.05,
    "Na": 0.05,
    "NO3": 0.05,
    "PO4": 0.005,
    "SO4": 0.05,
    "SiO2": 0.05,
    "pH": 0.10,
}

_CONCENTRATION_FIELDS: Final[frozenset[str]] = frozenset(
    {"Ca", "Cl", "F", "Fe", "HCO3", "K", "Mg", "Na", "NO3", "PO4", "SO4", "SiO2"}
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


def _finite(value: object, *, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be numeric, got {value!r}.")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric, got {value!r}.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite, got {value!r}.")
    return number


def _optional_finite(value: object) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _fingerprint(record: Mapping[str, Any]) -> str:
    payload = json.dumps(record, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _walk_forbidden_keys(value: Any, path: tuple[str, ...] = ()) -> None:
    """Reject truth/reference fields before locked prediction."""

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


def _normalise_family(label: object) -> str:
    """Map generator labels into the locked reaction taxonomy.

    Pure mixing is deliberately collapsed into the explicit null/mixing class.
    This prevents the model from being rewarded for calling transport mixing a
    chemical reaction.
    """

    text = str(label or "").strip().lower()
    if text in {"", "none", "no_reaction", "null", "missing", "source"}:
        return NULL_FAMILY
    if "mixing" in text or text == "mix":
        return "mixing"
    if any(token in text for token in ("calcite", "dolomite", "carbonate")):
        return "carbonate"
    if any(
        token in text
        for token in ("albite", "anorthite", "feldspar", "silicate", "exchange")
    ):
        return "silicate_exchange"
    if "sulfate_reduction" in text:
        return "sulfate_reduction"
    if "iron_reduction" in text:
        return "iron_reduction"
    if "denit" in text:
        return "denitrification"
    if "pyrite" in text:
        return "other_redox"
    if any(
        token in text
        for token in ("sulfate_source", "sulfate source", "gypsum", "so4")
    ):
        return "sulfate_source"
    return "other"


def normalise_rapm_reaction_family(label: object) -> str:
    """Return the RAPM taxonomy label used for development/locked scoring."""

    return _normalise_family(label)


def _softmax(logits: np.ndarray, temperature: float) -> np.ndarray:
    scaled = np.asarray(logits, dtype=float) / float(temperature)
    shifted = scaled - float(np.max(scaled))
    weights = np.exp(shifted)
    total = float(np.sum(weights))
    if not math.isfinite(total) or total <= 0.0:
        return np.full(len(scaled), 1.0 / len(scaled), dtype=float)
    return weights / total


def _mapping_softmax(
    logits: Mapping[str, object],
    classes: Sequence[str],
    *,
    temperature: float = 1.0,
    class_offsets: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """Return a finite probability vector for a sparse logit mapping."""

    selected: dict[str, float] = {}
    for family in classes:
        value = _optional_finite(logits.get(family, 0.0))
        offset = 0.0 if class_offsets is None else float(class_offsets.get(str(family), 0.0))
        selected[str(family)] = (0.0 if value is None else value) + offset
    scaled = {family: value / float(temperature) for family, value in selected.items()}
    maximum = max(scaled.values())
    weights = {family: math.exp(value - maximum) for family, value in scaled.items()}
    total = sum(weights.values())
    if not math.isfinite(total) or total <= 0.0:
        return {family: 1.0 / len(classes) for family in classes}
    return {family: value / total for family, value in weights.items()}


def _edge_items(observations: Mapping[str, Any]) -> list[tuple[str, Mapping[str, Any]]]:
    channel = observations.get(EVIDENCE_CHANNEL)
    if not isinstance(channel, Mapping):
        return []
    raw_edges = channel.get("edges")
    if isinstance(raw_edges, Mapping):
        return [
            (str(edge_id), edge)
            for edge_id, edge in raw_edges.items()
            if isinstance(edge, Mapping)
        ]
    if isinstance(raw_edges, Sequence) and not isinstance(
        raw_edges, (str, bytes, bytearray)
    ):
        result: list[tuple[str, Mapping[str, Any]]] = []
        for index, edge in enumerate(raw_edges):
            if not isinstance(edge, Mapping):
                continue
            edge_id = edge.get("edge_id", edge.get("id", f"edge_{index}"))
            result.append((str(edge_id), edge))
        return result
    return []


def _side(
    edge: Mapping[str, Any], primary: str, fallback: str
) -> Mapping[str, Any] | None:
    value = edge.get(primary, edge.get(fallback))
    return value if isinstance(value, Mapping) else None


@dataclass(frozen=True)
class ReactionRAPMConfig:
    """Predeclared fitting, weighting, and abstention conventions."""

    ion_scales: Mapping[str, float] = field(
        default_factory=lambda: dict(DEFAULT_ION_SCALES)
    )
    classes: tuple[str, ...] = DEFAULT_REACTION_CLASSES
    ridge_lambda_grid: tuple[float, ...] = (0.01, 0.1, 1.0, 10.0)
    on_off_weight_grid: tuple[float, ...] = (0.0, 0.25, 0.5, 1.0)
    default_ridge_lambda: float = 1.0
    default_on_off_weight: float = 0.25
    on_off_shrinkage: float = 4.0
    probability_temperature: float = 1.0
    selection_threshold: float = 0.60
    selection_margin: float = 0.10
    min_paired_fields: int = 2

    def __post_init__(self) -> None:
        scales: dict[str, float] = {}
        for raw_name, raw_value in self.ion_scales.items():
            name = str(raw_name).strip()
            value = _finite(raw_value, name=f"ion_scales[{name}]")
            if not name or value <= 0.0:
                raise ValueError("ion_scales must contain finite positive values")
            scales[name] = value
        classes = tuple(dict.fromkeys(str(item).strip() for item in self.classes))
        if NULL_FAMILY not in classes or len(classes) < 2:
            raise ValueError("classes must contain at least two families and 'none'")
        if any(not item for item in classes):
            raise ValueError("classes must not contain empty names")
        for name, values in (
            ("ridge_lambda_grid", self.ridge_lambda_grid),
            ("on_off_weight_grid", self.on_off_weight_grid),
        ):
            if not tuple(values):
                raise ValueError(f"{name} must not be empty")
            if any(_finite(value, name=name) < 0.0 for value in values):
                raise ValueError(f"{name} must contain non-negative values")
        for name, value in (
            ("default_ridge_lambda", self.default_ridge_lambda),
            ("default_on_off_weight", self.default_on_off_weight),
            ("on_off_shrinkage", self.on_off_shrinkage),
            ("probability_temperature", self.probability_temperature),
            ("selection_threshold", self.selection_threshold),
            ("selection_margin", self.selection_margin),
        ):
            if _finite(value, name=name) < 0.0:
                raise ValueError(f"{name} must be non-negative")
        if self.probability_temperature <= 0.0:
            raise ValueError("probability_temperature must be positive")
        if not 0.0 < self.selection_threshold < 1.0:
            raise ValueError("selection_threshold must lie in (0, 1)")
        if not 0.0 <= self.selection_margin < 1.0:
            raise ValueError("selection_margin must lie in [0, 1)")
        if (
            isinstance(self.min_paired_fields, bool)
            or not isinstance(self.min_paired_fields, int)
            or self.min_paired_fields < 1
        ):
            raise ValueError("min_paired_fields must be a positive integer")
        object.__setattr__(self, "ion_scales", dict(sorted(scales.items())))
        object.__setattr__(self, "classes", classes)
        object.__setattr__(
            self,
            "ridge_lambda_grid",
            tuple(sorted({_finite(value, name="ridge_lambda_grid") for value in self.ridge_lambda_grid})),
        )
        object.__setattr__(
            self,
            "on_off_weight_grid",
            tuple(sorted({_finite(value, name="on_off_weight_grid") for value in self.on_off_weight_grid})),
        )


@dataclass(frozen=True)
class ReactionRAPMTrainingExample:
    """One labelled development edge used only by :meth:`ReactionRAPM.fit`."""

    case_id: str
    edge_id: str
    upstream: Mapping[str, Any]
    downstream: Mapping[str, Any]
    truth_family: str

    def __post_init__(self) -> None:
        for name in ("case_id", "edge_id", "truth_family"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} is required")
        if not isinstance(self.upstream, Mapping) or not isinstance(self.downstream, Mapping):
            raise TypeError("Training chemistry sides must be mappings")


@dataclass(frozen=True)
class ReactionRAPMCalibrationExample:
    """One case-held-out RAPM logit vector used for development calibration.

    The ``cross_fitted`` flag is deliberately mandatory evidence in the data
    contract.  Calibration on in-sample logits would make the probability
    layer look better than it is and is therefore rejected by the calibrator.
    Truth is accepted here only for development fitting; it is never accepted
    by :meth:`ReactionRAPM.predict` or :meth:`ReactionRAPMCalibrator.apply`.
    """

    case_id: str
    edge_id: str
    truth_family: str
    logits: Mapping[str, Any]
    cross_fitted: bool = True

    def __post_init__(self) -> None:
        for name in ("case_id", "edge_id", "truth_family"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} is required")
        if not self.cross_fitted:
            raise ValueError("Reaction calibration requires case-held-out logits")
        if not isinstance(self.logits, Mapping) or not self.logits:
            raise ValueError("Reaction calibration requires non-empty logits")
        clean = {
            str(key): _finite(value, name=f"logit[{key}]")
            for key, value in self.logits.items()
        }
        object.__setattr__(self, "case_id", str(self.case_id))
        object.__setattr__(self, "edge_id", str(self.edge_id))
        object.__setattr__(self, "truth_family", _normalise_family(self.truth_family))
        object.__setattr__(self, "logits", clean)


@dataclass(frozen=True)
class _FeatureRow:
    case_id: str
    edge_id: str
    family: str
    values: tuple[float, ...]
    paired_fields: tuple[str, ...]
    missing_fields: tuple[str, ...]
    unsupported_fields: tuple[str, ...]
    invalid_fields: tuple[str, ...]
    status: str | None = None


def _feature_names(config: ReactionRAPMConfig) -> tuple[str, ...]:
    names = [f"delta_{field_name}" for field_name in config.ion_scales]
    names.extend(f"present_{field_name}" for field_name in config.ion_scales)
    names.append("delta_norm")
    names.extend(f"template_score_{family}" for family in TEMPLATE_EVIDENCE_FAMILIES)
    return tuple(names)


def _template_evidence_scores(delta: Mapping[str, float]) -> dict[str, float]:
    """Create truth-blind signed stoichiometric evidence covariates.

    These are not equilibrium calculations and are not treated as reaction
    truth.  They expose fixed signed observation-space patterns as regularized
    covariates for RAPM, allowing the learned on/off layer to adjust them
    across generator families.
    """

    scores = {family: -2.0 for family in TEMPLATE_EVIDENCE_FAMILIES}
    for family, template in TEMPLATE_EVIDENCE_LIBRARY.items():
        available = [ion for ion in template if ion in delta]
        if not available:
            continue
        basis = np.asarray([template[ion] for ion in available], dtype=float)
        observed = np.asarray([delta[ion] for ion in available], dtype=float)
        denominator = float(np.dot(basis, basis))
        extent = max(0.0, float(np.dot(observed, basis)) / denominator) if denominator else 0.0
        residual = float(np.sqrt(np.mean((observed - extent * basis) ** 2)))
        scores[family] = -residual + 0.15 * math.log1p(extent) - 0.25
    delta_norm = math.sqrt(
        sum(value * value for value in delta.values()) / max(len(delta), 1)
    )
    scores[NULL_FAMILY] = -delta_norm
    scores["mixing"] = -0.5 * delta_norm - 0.05
    return scores


def _make_feature_row(
    *,
    case_id: str,
    edge_id: str,
    upstream: Mapping[str, Any],
    downstream: Mapping[str, Any],
    family: str,
    config: ReactionRAPMConfig,
) -> _FeatureRow:
    deltas: list[float] = []
    delta_by_field: dict[str, float] = {}
    presence: list[float] = []
    paired: list[str] = []
    missing: list[str] = []
    invalid: list[str] = []
    for field_name, scale in config.ion_scales.items():
        if field_name not in upstream or field_name not in downstream:
            deltas.append(0.0)
            presence.append(0.0)
            missing.append(field_name)
            continue
        before = _optional_finite(upstream[field_name])
        after = _optional_finite(downstream[field_name])
        if before is None or after is None:
            deltas.append(0.0)
            presence.append(0.0)
            invalid.append(field_name)
            continue
        if field_name in _CONCENTRATION_FIELDS and (before < 0.0 or after < 0.0):
            deltas.append(0.0)
            presence.append(0.0)
            invalid.append(field_name)
            continue
        delta = (after - before) / scale
        deltas.append(delta)
        delta_by_field[field_name] = delta
        presence.append(1.0)
        paired.append(field_name)
    supplied = {str(name) for name in upstream} | {str(name) for name in downstream}
    unsupported = sorted(supplied - set(config.ion_scales))
    status: str | None = None
    if invalid:
        status = "invalid_chemistry_values"
    elif len(paired) < config.min_paired_fields:
        status = (
            "unsupported_chemistry_fields"
            if not paired and unsupported
            else "insufficient_paired_fields"
        )
    magnitude = (
        math.sqrt(sum(value * value for value in deltas) / len(paired))
        if paired
        else 0.0
    )
    template_scores = _template_evidence_scores(delta_by_field)
    values = tuple(
        deltas
        + presence
        + [magnitude]
        + [template_scores[family] for family in TEMPLATE_EVIDENCE_FAMILIES]
    )
    return _FeatureRow(
        case_id=str(case_id),
        edge_id=str(edge_id),
        family=_normalise_family(family),
        values=values,
        paired_fields=tuple(paired),
        missing_fields=tuple(missing),
        unsupported_fields=tuple(unsupported),
        invalid_fields=tuple(invalid),
        status=status,
    )


def _training_rows(
    examples: Iterable[ReactionRAPMTrainingExample],
    config: ReactionRAPMConfig,
) -> tuple[tuple[_FeatureRow, ...], int]:
    rows: list[_FeatureRow] = []
    skipped = 0
    for example in examples:
        row = _make_feature_row(
            case_id=example.case_id,
            edge_id=example.edge_id,
            upstream=example.upstream,
            downstream=example.downstream,
            family=example.truth_family,
            config=config,
        )
        if row.status is not None:
            skipped += 1
            continue
        rows.append(row)
    return tuple(rows), skipped


def _standardise(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.mean(values, axis=0)
    scale = np.std(values, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-8), scale, 1.0)
    return (values - mean) / scale, mean, scale


def _case_class_weights(rows: Sequence[_FeatureRow]) -> np.ndarray:
    case_counts = Counter(row.case_id for row in rows)
    class_counts = Counter(row.family for row in rows)
    class_count = max(len(class_counts), 1)
    weights = np.asarray(
        [
            1.0
            / float(case_counts[row.case_id])
            * float(len(rows))
            / float(class_count * class_counts[row.family])
            for row in rows
        ],
        dtype=float,
    )
    mean = float(np.mean(weights))
    return weights / mean if mean > 0.0 else np.ones(len(rows), dtype=float)


def _fit_parameters(
    rows: Sequence[_FeatureRow],
    config: ReactionRAPMConfig,
    ridge_lambda: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray([row.values for row in rows], dtype=float)
    scaled, mean, scale = _standardise(values)
    design = np.column_stack([np.ones(len(rows), dtype=float), scaled])
    classes = config.classes
    targets = np.zeros((len(rows), len(classes)), dtype=float)
    class_index = {name: index for index, name in enumerate(classes)}
    for row_index, row in enumerate(rows):
        targets[row_index, class_index.get(row.family, class_index["other"])] = 1.0
    weights = _case_class_weights(rows)
    weighted_design = design * weights[:, None]
    left = design.T @ weighted_design
    penalty = np.eye(design.shape[1], dtype=float) * float(ridge_lambda)
    penalty[0, 0] = 0.0
    right = design.T @ (targets * weights[:, None])
    coefficients = np.linalg.solve(left + penalty + np.eye(left.shape[0]) * 1.0e-10, right)

    effects = np.zeros((len(classes), scaled.shape[1]), dtype=float)
    shrinkage = float(config.on_off_shrinkage)
    for class_index_value, family in enumerate(classes):
        on = np.asarray([row.family == family for row in rows], dtype=bool)
        off = ~on
        if not np.any(on) or not np.any(off):
            continue
        on_weights = weights * on
        off_weights = weights * off
        on_mean = np.sum(scaled * on_weights[:, None], axis=0) / max(float(np.sum(on_weights)), 1.0e-12)
        off_mean = np.sum(scaled * off_weights[:, None], axis=0) / max(float(np.sum(off_weights)), 1.0e-12)
        on_shrink = float(np.sum(on_weights)) / (float(np.sum(on_weights)) + shrinkage)
        off_shrink = float(np.sum(off_weights)) / (float(np.sum(off_weights)) + shrinkage)
        effects[class_index_value] = (on_mean - off_mean) * min(on_shrink, off_shrink)
    return coefficients, effects, mean, scale, weights


def _logits_for_values(
    values: np.ndarray,
    *,
    coefficients: np.ndarray,
    effects: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
    on_off_weight: float,
) -> np.ndarray:
    scaled = (values - mean) / scale
    design = np.column_stack([np.ones(len(values), dtype=float), scaled])
    return design @ coefficients + float(on_off_weight) * (scaled @ effects.T)


def _case_blocked_cv(
    rows: Sequence[_FeatureRow],
    config: ReactionRAPMConfig,
) -> tuple[float, float, float, int]:
    case_ids = tuple(sorted({row.case_id for row in rows}))
    if len(case_ids) < 2:
        return (
            float("nan"),
            float(config.default_ridge_lambda),
            float(config.default_on_off_weight),
            0,
        )
    best: tuple[float, float, float] | None = None
    for ridge_lambda in config.ridge_lambda_grid:
        for on_off_weight in config.on_off_weight_grid:
            losses: list[float] = []
            for held_out in case_ids:
                train = tuple(row for row in rows if row.case_id != held_out)
                test = tuple(row for row in rows if row.case_id == held_out)
                if not train or not test:
                    continue
                coefficients, effects, mean, scale, _ = _fit_parameters(
                    train, config, ridge_lambda
                )
                values = np.asarray([row.values for row in test], dtype=float)
                logits = _logits_for_values(
                    values,
                    coefficients=coefficients,
                    effects=effects,
                    mean=mean,
                    scale=scale,
                    on_off_weight=on_off_weight,
                )
                class_index = {name: index for index, name in enumerate(config.classes)}
                row_losses = []
                for row, row_logits in zip(test, logits, strict=True):
                    probabilities = _softmax(row_logits, config.probability_temperature)
                    probability = probabilities[class_index.get(row.family, class_index["other"])]
                    row_losses.append(-math.log(max(float(probability), 1.0e-12)))
                losses.append(float(np.mean(row_losses)))
            if not losses:
                continue
            candidate = (float(np.mean(losses)), float(ridge_lambda), float(on_off_weight))
            if best is None or candidate < best:
                best = candidate
    if best is None:
        return (
            float("nan"),
            float(config.default_ridge_lambda),
            float(config.default_on_off_weight),
            0,
        )
    return best[0], best[1], best[2], len(case_ids)


def _design_diagnostics(
    rows: Sequence[_FeatureRow],
) -> tuple[int, tuple[float, ...], float | None, float | None]:
    """Summarise rank and signature coherence without interpreting causes."""

    values = np.asarray([row.values for row in rows], dtype=float)
    scaled, _mean, _scale = _standardise(values)
    singular_values = np.linalg.svd(scaled, compute_uv=False)
    rank = int(np.linalg.matrix_rank(scaled))
    positive = singular_values[singular_values > 1.0e-10]
    condition_number = (
        float(positive[0] / positive[-1])
        if positive.size and positive[-1] > 0.0
        else None
    )
    coherence: float | None = None
    variable_columns = np.std(scaled, axis=0) > 1.0e-8
    if int(np.sum(variable_columns)) > 1 and scaled.shape[0] > 1:
        correlation = np.corrcoef(scaled[:, variable_columns], rowvar=False)
        off_diagonal = np.abs(correlation[np.triu_indices_from(correlation, k=1)])
        finite = off_diagonal[np.isfinite(off_diagonal)]
        if finite.size:
            coherence = float(np.max(finite))
    return (
        rank,
        tuple(float(value) for value in singular_values),
        condition_number,
        coherence,
    )


def _parameter_stability(
    rows: Sequence[_FeatureRow],
    config: ReactionRAPMConfig,
    *,
    ridge_lambda: float,
    full_coefficients: np.ndarray,
    full_effects: np.ndarray,
) -> dict[str, float | int | None]:
    """Compare fixed-hyperparameter fits across held-out development cases."""

    relative_coefficients: list[float] = []
    relative_effects: list[float] = []
    for held_out in sorted({row.case_id for row in rows}):
        train = tuple(row for row in rows if row.case_id != held_out)
        if len({row.case_id for row in train}) < 2:
            continue
        coefficients, effects, _mean, _scale, _weights = _fit_parameters(
            train,
            config,
            ridge_lambda,
        )
        coefficient_scale = max(float(np.linalg.norm(full_coefficients)), 1.0e-12)
        effect_scale = max(float(np.linalg.norm(full_effects)), 1.0e-12)
        relative_coefficients.append(
            float(np.linalg.norm(coefficients - full_coefficients) / coefficient_scale)
        )
        relative_effects.append(
            float(np.linalg.norm(effects - full_effects) / effect_scale)
        )

    if not relative_coefficients:
        return {
            "fold_count": 0,
            "coefficient_relative_l2_mean": None,
            "coefficient_relative_l2_max": None,
            "on_off_relative_l2_mean": None,
            "on_off_relative_l2_max": None,
        }
    return {
        "fold_count": len(relative_coefficients),
        "coefficient_relative_l2_mean": float(np.mean(relative_coefficients)),
        "coefficient_relative_l2_max": float(np.max(relative_coefficients)),
        "on_off_relative_l2_mean": float(np.mean(relative_effects)),
        "on_off_relative_l2_max": float(np.max(relative_effects)),
    }


@dataclass(frozen=True)
class ReactionRAPM:
    """Frozen development-fitted regularized reaction evidence model."""

    name: Final[str] = "reaction_rapm"
    version: Final[str] = "0.3.0"
    config: ReactionRAPMConfig = field(default_factory=ReactionRAPMConfig)
    feature_names: tuple[str, ...] = ()
    feature_means: tuple[float, ...] = ()
    feature_scales: tuple[float, ...] = ()
    coefficients: tuple[tuple[float, ...], ...] = ()
    on_off_effects: tuple[tuple[float, ...], ...] = ()
    ridge_lambda: float = 1.0
    on_off_weight: float = 0.25
    cv_log_loss: float | None = None
    cv_case_count: int = 0
    fit_case_ids: tuple[str, ...] = ()
    fit_record_count: int = 0
    skipped_training_count: int = 0
    class_counts: Mapping[str, int] = field(default_factory=dict)
    design_rank: int = 0
    design_singular_values: tuple[float, ...] = ()
    design_condition_number: float | None = None
    signature_coherence: float | None = None
    parameter_stability: Mapping[str, float | int | None] = field(default_factory=dict)

    @classmethod
    def fit(
        cls,
        examples: Iterable[ReactionRAPMTrainingExample],
        *,
        config: ReactionRAPMConfig | None = None,
        phase: str = "development",
    ) -> "ReactionRAPM":
        if str(phase) != "development":
            raise ValueError("Reaction RAPM fitting is restricted to development data")
        config = config or ReactionRAPMConfig()
        rows, skipped = _training_rows(examples, config)
        if len(rows) < 2:
            raise ValueError("Reaction RAPM requires at least two valid development edges")
        if len({row.case_id for row in rows}) < 2:
            raise ValueError("Reaction RAPM requires at least two development cases")
        cv_loss, ridge_lambda, on_off_weight, cv_case_count = _case_blocked_cv(rows, config)
        coefficients, effects, mean, scale, _ = _fit_parameters(
            rows, config, ridge_lambda
        )
        design_rank, singular_values, condition_number, signature_coherence = (
            _design_diagnostics(rows)
        )
        parameter_stability = _parameter_stability(
            rows,
            config,
            ridge_lambda=ridge_lambda,
            full_coefficients=coefficients,
            full_effects=effects,
        )
        counts = Counter(row.family for row in rows)
        return cls(
            config=config,
            feature_names=_feature_names(config),
            feature_means=tuple(float(value) for value in mean),
            feature_scales=tuple(float(value) for value in scale),
            coefficients=tuple(tuple(float(value) for value in row) for row in coefficients),
            on_off_effects=tuple(tuple(float(value) for value in row) for row in effects),
            ridge_lambda=float(ridge_lambda),
            on_off_weight=float(on_off_weight),
            cv_log_loss=None if not math.isfinite(cv_loss) else float(cv_loss),
            cv_case_count=int(cv_case_count),
            fit_case_ids=tuple(sorted({row.case_id for row in rows})),
            fit_record_count=len(rows),
            skipped_training_count=skipped,
            class_counts=dict(sorted(counts.items())),
            design_rank=design_rank,
            design_singular_values=singular_values,
            design_condition_number=condition_number,
            signature_coherence=signature_coherence,
            parameter_stability=parameter_stability,
        )

    def _predict_row(self, row: _FeatureRow) -> dict[str, Any]:
        if row.status is not None:
            return self._abstention_record(row)
        values = np.asarray([row.values], dtype=float)
        coefficients = np.asarray(self.coefficients, dtype=float)
        effects = np.asarray(self.on_off_effects, dtype=float)
        mean = np.asarray(self.feature_means, dtype=float)
        scale = np.asarray(self.feature_scales, dtype=float)
        logits = _logits_for_values(
            values,
            coefficients=coefficients,
            effects=effects,
            mean=mean,
            scale=scale,
            on_off_weight=self.on_off_weight,
        )[0]
        probabilities_array = _softmax(logits, self.config.probability_temperature)
        probabilities = {
            family: float(probabilities_array[index])
            for index, family in enumerate(self.config.classes)
        }
        ranking = sorted(
            self.config.classes,
            key=lambda family: (-probabilities[family], family),
        )
        best = ranking[0]
        second = probabilities[ranking[1]] if len(ranking) > 1 else 0.0
        if probabilities[best] < self.config.selection_threshold:
            decision = ABSTAIN
            reason = "top_probability_below_threshold"
        elif probabilities[best] - second < self.config.selection_margin:
            decision = ABSTAIN
            reason = "adjusted_family_margin_below_threshold"
        else:
            decision = SELECT
            reason = (
                "null_mixing_control_selected"
                if best == NULL_FAMILY
                else "regularized_adjusted_reaction_family_selected"
            )
        return {
            "edge_id": row.edge_id,
            "family": best if decision == SELECT else NULL_FAMILY,
            "candidate_family": best,
            "probability": probabilities[best],
            "probabilities": probabilities,
            "logits": {
                family: float(logits[index])
                for index, family in enumerate(self.config.classes)
            },
            "decision": decision,
            "reason": reason,
            "evidence_channel": EVIDENCE_CHANNEL,
            "uncertainty": {
                "calibrated": False,
                "selection_probability_calibrated": False,
                "method": "case_blocked_ridge_plus_on_off_effects",
            },
            "diagnostics": {
                "paired_fields": list(row.paired_fields),
                "n_paired_fields": len(row.paired_fields),
                "missing_fields": list(row.missing_fields),
                "unsupported_fields": list(row.unsupported_fields),
                "invalid_fields": list(row.invalid_fields),
                "ridge_lambda": self.ridge_lambda,
                "on_off_weight": self.on_off_weight,
                "cv_log_loss": self.cv_log_loss,
                "cv_case_count": self.cv_case_count,
                "feature_names": list(self.feature_names),
                "on_off_scores": {
                    family: float(
                        np.dot(
                            (np.asarray(row.values) - mean) / scale,
                            effects[index],
                        )
                    )
                    for index, family in enumerate(self.config.classes)
                },
            },
        }

    def predict(self, observations: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
        """Predict from paired chemistry without accessing truth."""

        _walk_forbidden_keys(observations)
        if not isinstance(observations.get(EVIDENCE_CHANNEL), Mapping):
            return {}
        result: dict[str, dict[str, Any]] = {}
        for edge_id, edge in _edge_items(observations):
            upstream = _side(edge, "upstream", "source")
            downstream = _side(edge, "downstream", "target_values")
            if downstream is None:
                downstream = _side(edge, "downstream", "destination")
            if upstream is None or downstream is None:
                row = _FeatureRow(
                    case_id="prediction",
                    edge_id=edge_id,
                    family=NULL_FAMILY,
                    values=tuple(0.0 for _ in self.feature_names),
                    paired_fields=(),
                    missing_fields=(),
                    unsupported_fields=(),
                    invalid_fields=(),
                    status="malformed_edge",
                )
            else:
                row = _make_feature_row(
                    case_id="prediction",
                    edge_id=edge_id,
                    upstream=upstream,
                    downstream=downstream,
                    family=NULL_FAMILY,
                    config=self.config,
                )
            result[edge_id] = self._predict_row(row)
        return result

    @staticmethod
    def _abstention_record(row: _FeatureRow) -> dict[str, Any]:
        return {
            "edge_id": row.edge_id,
            "family": NULL_FAMILY,
            "candidate_family": None,
            "probability": 0.0,
            "probabilities": {NULL_FAMILY: 1.0},
            "logits": {NULL_FAMILY: 0.0},
            "decision": ABSTAIN,
            "reason": row.status or "unsupported_reaction_chemistry",
            "evidence_channel": EVIDENCE_CHANNEL,
            "uncertainty": {
                "calibrated": False,
                "selection_probability_calibrated": False,
            },
            "diagnostics": {
                "paired_fields": list(row.paired_fields),
                "n_paired_fields": len(row.paired_fields),
                "missing_fields": list(row.missing_fields),
                "unsupported_fields": list(row.unsupported_fields),
                "invalid_fields": list(row.invalid_fields),
            },
        }

    def to_audit_record(self) -> dict[str, Any]:
        record: dict[str, Any] = {
            "name": self.name,
            "version": self.version,
            "family": "reaction_family",
            "output_kind": REACTION_OUTPUT,
            "input_channels": [EVIDENCE_CHANNEL],
            "description": (
                "Development-fitted case-blocked ridge reaction evidence model "
                "with a separate on/off family-effect diagnostic."
            ),
            "fit_scope": "development_only",
            "features": list(self.feature_names),
            "template_evidence": {
                "families": list(TEMPLATE_EVIDENCE_FAMILIES),
                "library": {
                    family: dict(template)
                    for family, template in TEMPLATE_EVIDENCE_LIBRARY.items()
                },
                "role": "truth_blind_observation_space_covariates_for_RAPM",
                "not_a_phreeqc_inverse": True,
            },
            "classes": list(self.config.classes),
            "selected_hyperparameters": {
                "ridge_lambda": self.ridge_lambda,
                "on_off_weight": self.on_off_weight,
                "probability_temperature": self.config.probability_temperature,
                "selection_threshold": self.config.selection_threshold,
                "selection_margin": self.config.selection_margin,
            },
            "training": {
                "case_count": len(self.fit_case_ids),
                "record_count": self.fit_record_count,
                "skipped_training_count": self.skipped_training_count,
                "class_counts": dict(self.class_counts),
                "case_blocked_cv_case_count": self.cv_case_count,
                "case_blocked_cv_log_loss": self.cv_log_loss,
                "truth_used_for": "development_fit_only",
            },
            "identifiability": {
                "design_rank": self.design_rank,
                "design_feature_count": len(self.feature_names),
                "design_singular_values": list(self.design_singular_values),
                "design_condition_number": self.design_condition_number,
                "maximum_absolute_feature_coherence": self.signature_coherence,
                "parameter_stability_across_case_folds": dict(self.parameter_stability),
                "interpretation": (
                    "ridge stability is not unique chemical-mechanism identification; "
                    "overlapping signatures require abstention or an equivalence-class claim"
                ),
            },
            "uncertainty": {
                "calibrated": False,
                "selection_probability_calibrated": False,
                "calibration_contract": "apply frozen development-only temperature scaling before locked claim scoring",
            },
            "abstention": {
                "unsupported_input": True,
                "invalid_input": True,
                "insufficient_paired_fields": True,
                "top_probability_below_threshold": True,
                "adjusted_family_margin_below_threshold": True,
            },
            "control": {
                "truth_blind_prediction": True,
                "uses_hydrosheaf_posterior": False,
                "uses_phreeqc_outputs_as_truth": False,
                "candidate_universe_scope": "caller_supplied_observed_chemistry_edges",
                "null_no_reaction_class": NULL_FAMILY,
                "mixing_class": "mixing",
                "signature_semantics": "learned normalized observation-space evidence, not balanced geochemical equations",
            },
            "parameter_fingerprint": _fingerprint(
                {
                    "feature_means": list(self.feature_means),
                    "feature_scales": list(self.feature_scales),
                    "coefficients": [list(row) for row in self.coefficients],
                    "on_off_effects": [list(row) for row in self.on_off_effects],
                    "ridge_lambda": self.ridge_lambda,
                    "on_off_weight": self.on_off_weight,
                }
            ),
        }
        record["fingerprint"] = _fingerprint(record)
        return record

    def to_dict(self) -> dict[str, Any]:
        """Serialize the frozen model and its provenance."""

        return {
            "audit": self.to_audit_record(),
            "config": {
                "ion_scales": dict(self.config.ion_scales),
                "classes": list(self.config.classes),
                "ridge_lambda_grid": list(self.config.ridge_lambda_grid),
                "on_off_weight_grid": list(self.config.on_off_weight_grid),
                "default_ridge_lambda": self.config.default_ridge_lambda,
                "default_on_off_weight": self.config.default_on_off_weight,
                "on_off_shrinkage": self.config.on_off_shrinkage,
                "probability_temperature": self.config.probability_temperature,
                "selection_threshold": self.config.selection_threshold,
                "selection_margin": self.config.selection_margin,
                "min_paired_fields": self.config.min_paired_fields,
            },
            "feature_names": list(self.feature_names),
            "feature_means": list(self.feature_means),
            "feature_scales": list(self.feature_scales),
            "coefficients": [list(row) for row in self.coefficients],
            "on_off_effects": [list(row) for row in self.on_off_effects],
            "ridge_lambda": self.ridge_lambda,
            "on_off_weight": self.on_off_weight,
            "cv_log_loss": self.cv_log_loss,
            "cv_case_count": self.cv_case_count,
            "fit_case_ids": list(self.fit_case_ids),
            "fit_record_count": self.fit_record_count,
            "skipped_training_count": self.skipped_training_count,
            "class_counts": dict(self.class_counts),
            "design_rank": self.design_rank,
            "design_singular_values": list(self.design_singular_values),
            "design_condition_number": self.design_condition_number,
            "signature_coherence": self.signature_coherence,
            "parameter_stability": dict(self.parameter_stability),
        }


@dataclass(frozen=True)
class ReactionRAPMCalibrator:
    """Frozen development-only temperature and selective-risk layer.

    Temperature is fitted on case-held-out development logits by multiclass
    log loss.  Optional acceptance-rule tuning is performed only on those
    cross-fitted development logits against a predeclared coverage/risk target.
    Applying this object therefore needs only raw truth-blind records.
    """

    temperature: float
    classes: tuple[str, ...]
    decision_threshold: float
    decision_margin: float
    fit_count: int
    case_ids: tuple[str, ...]
    fit_scope: str = "development_only"
    calibration_kind: str = "case_blocked_temperature_scaling"
    provenance: Mapping[str, Any] = field(default_factory=dict)
    class_offsets: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        temperature = _finite(self.temperature, name="temperature")
        threshold = _finite(self.decision_threshold, name="decision_threshold")
        margin = _finite(self.decision_margin, name="decision_margin")
        if temperature <= 0.0:
            raise ValueError("temperature must be positive")
        if not 0.0 < threshold < 1.0:
            raise ValueError("decision_threshold must lie in (0, 1)")
        if not 0.0 <= margin < 1.0:
            raise ValueError("decision_margin must lie in [0, 1)")
        classes = tuple(dict.fromkeys(str(item).strip() for item in self.classes))
        if len(classes) < 2 or any(not item for item in classes):
            raise ValueError("Reaction calibration requires at least two classes")
        if int(self.fit_count) < 1 or not self.case_ids:
            raise ValueError("Reaction calibration requires observations and cases")
        if str(self.fit_scope) != "development_only":
            raise ValueError("Reaction calibration must be fitted on development data only")
        object.__setattr__(self, "temperature", temperature)
        object.__setattr__(self, "decision_threshold", threshold)
        object.__setattr__(self, "decision_margin", margin)
        object.__setattr__(self, "classes", classes)
        object.__setattr__(self, "fit_count", int(self.fit_count))
        object.__setattr__(self, "case_ids", tuple(sorted(set(str(item) for item in self.case_ids))))
        offsets = tuple(float(value) for value in self.class_offsets)
        if not offsets:
            offsets = tuple(0.0 for _ in classes)
        if len(offsets) != len(classes) or any(not math.isfinite(value) for value in offsets):
            raise ValueError("class_offsets must contain one finite value per class")
        mean_offset = float(np.mean(offsets))
        object.__setattr__(self, "class_offsets", tuple(value - mean_offset for value in offsets))

    @classmethod
    def fit(
        cls,
        examples: Iterable[ReactionRAPMCalibrationExample],
        *,
        classes: Sequence[str] = DEFAULT_REACTION_CLASSES,
        decision_threshold: float = 0.60,
        decision_margin: float = 0.10,
        tune_selection_rule: bool = False,
        target_coverage: float = 0.25,
        max_selective_risk: float = 0.40,
        phase: str = "development",
    ) -> "ReactionRAPMCalibrator":
        if str(phase) != "development":
            raise ValueError("Reaction calibration may only fit development data")
        rows = tuple(examples)
        if not rows:
            raise ValueError("At least one reaction calibration example is required")
        if len({row.case_id for row in rows}) < 2:
            raise ValueError("Reaction calibration requires at least two cases")
        class_names = tuple(
            dict.fromkeys(
                [str(item) for item in classes]
                + [row.truth_family for row in rows]
                + [str(key) for row in rows for key in row.logits]
            )
        )
        by_case: dict[str, list[ReactionRAPMCalibrationExample]] = {}
        for row in rows:
            by_case.setdefault(row.case_id, []).append(row)

        def loss(
            temperature: float,
            class_offsets: Mapping[str, float] | None = None,
        ) -> float:
            case_losses: list[float] = []
            for case_rows in by_case.values():
                row_losses: list[float] = []
                for row in case_rows:
                    probabilities = _mapping_softmax(
                        row.logits,
                        class_names,
                        temperature=temperature,
                        class_offsets=class_offsets,
                    )
                    row_losses.append(
                        -math.log(max(probabilities.get(row.truth_family, 0.0), 1.0e-12))
                    )
                if row_losses:
                    case_losses.append(float(np.mean(row_losses)))
            return float(np.mean(case_losses)) if case_losses else float("inf")

        candidate_temperatures = tuple(
            math.exp(-2.5 + index * (5.0 / 160.0)) for index in range(161)
        )
        temperature = min(candidate_temperatures, key=lambda value: (loss(value), value))
        offset_values = np.zeros(len(class_names), dtype=float)
        offset_grid = tuple(float(value) for value in np.linspace(-1.5, 1.5, 31))

        def offset_objective(values: np.ndarray) -> float:
            centered = values - float(np.mean(values))
            offsets = {
                family: float(centered[index])
                for index, family in enumerate(class_names)
            }
            return loss(temperature, offsets) + 0.02 * float(np.mean(centered * centered))

        # A small, deterministic vector-scaling layer corrects class-prior
        # distortions left by temperature-only calibration.  The offsets are
        # fit on case-held-out development logits and regularised so that this
        # remains a probability calibration layer rather than a second RAPM
        # predictor.
        for _ in range(6):
            previous_objective = offset_objective(offset_values)
            for index in range(len(class_names)):
                candidates = []
                for value in offset_grid:
                    trial = offset_values.copy()
                    trial[index] = value
                    candidates.append((offset_objective(trial), value))
                best_objective, best_value = min(candidates, key=lambda item: (item[0], item[1]))
                if best_objective <= previous_objective + 1.0e-12:
                    offset_values[index] = best_value
                    previous_objective = best_objective
            if np.max(np.abs(offset_values - np.mean(offset_values))) < 1.0e-8:
                break
        offset_values -= float(np.mean(offset_values))
        selected_offsets = {
            family: float(offset_values[index])
            for index, family in enumerate(class_names)
        }
        selected_threshold = float(decision_threshold)
        selected_margin = float(decision_margin)
        selection_tuning: dict[str, Any] = {
            "enabled": bool(tune_selection_rule),
            "target_coverage": float(target_coverage),
            "max_selective_risk": float(max_selective_risk),
        }
        if tune_selection_rule:
            if not 0.0 < float(target_coverage) <= 1.0:
                raise ValueError("target_coverage must lie in (0, 1]")
            if not 0.0 <= float(max_selective_risk) <= 1.0:
                raise ValueError("max_selective_risk must lie in [0, 1]")
            threshold_grid = tuple(
                sorted(
                    {
                        float(decision_threshold),
                        0.35,
                        0.40,
                        0.45,
                        0.50,
                        0.55,
                        0.65,
                        0.70,
                    }
                )
            )
            margin_grid = tuple(sorted({float(decision_margin), 0.0, 0.05, 0.15, 0.20}))
            candidates: list[tuple[float, float, float, float, int]] = []
            for threshold in threshold_grid:
                if not 0.0 < threshold < 1.0:
                    continue
                for margin in margin_grid:
                    if not 0.0 <= margin < 1.0:
                        continue
                    accepted = 0
                    correct = 0
                    for row in rows:
                        probabilities = _mapping_softmax(
                            row.logits,
                            class_names,
                            temperature=temperature,
                            class_offsets=selected_offsets,
                        )
                        ranking = sorted(
                            class_names,
                            key=lambda family: (-probabilities[family], family),
                        )
                        best = ranking[0]
                        second = probabilities[ranking[1]] if len(ranking) > 1 else 0.0
                        if probabilities[best] >= threshold and probabilities[best] - second >= margin:
                            accepted += 1
                            correct += int(best == row.truth_family)
                    coverage = accepted / len(rows) if rows else 0.0
                    risk = 1.0 - correct / accepted if accepted else 1.0
                    candidates.append((coverage, risk, threshold, margin, accepted))
            eligible = [
                candidate
                for candidate in candidates
                if candidate[0] >= float(target_coverage)
                and candidate[1] <= float(max_selective_risk)
            ]
            chosen = max(
                eligible or candidates,
                key=lambda candidate: (
                    candidate[0] if eligible else -candidate[1],
                    -candidate[1] if eligible else candidate[0],
                    candidate[2],
                    candidate[3],
                ),
            )
            selected_threshold = chosen[2]
            selected_margin = chosen[3]
            selection_tuning.update(
                {
                    "selected_development_coverage": chosen[0],
                    "selected_development_selective_risk": chosen[1],
                    "selected_threshold": selected_threshold,
                    "selected_margin": selected_margin,
                    "target_met": bool(eligible),
                    "candidate_count": len(candidates),
                }
            )
        boundary = temperature in {
            candidate_temperatures[0],
            candidate_temperatures[-1],
        }
        calibration_kind = (
            "case_blocked_temperature_and_selective_rule_calibration"
            if tune_selection_rule
            else "case_blocked_temperature_scaling"
        )
        return cls(
            temperature=float(temperature),
            classes=class_names,
            decision_threshold=float(selected_threshold),
            decision_margin=float(selected_margin),
            fit_count=len(rows),
            case_ids=tuple(sorted(by_case)),
            calibration_kind=calibration_kind,
            class_offsets=tuple(float(value) for value in offset_values),
            provenance={
                "fit_rule": "case_equal_weight_multiclass_log_score_grid",
                "truth_used_for": "development_fit_only",
                "threshold_tuning": dict(selection_tuning),
                "margin_tuning": dict(selection_tuning),
                "temperature_grid": [
                    float(candidate_temperatures[0]),
                    float(candidate_temperatures[-1]),
                    len(candidate_temperatures),
                ],
                "grid_boundary_hit": boundary,
                "class_offset_grid": [
                    float(offset_grid[0]),
                    float(offset_grid[-1]),
                    len(offset_grid),
                ],
                "class_offset_regularization": 0.02,
                "requires_model_audit": (
                    "retain RAPM design rank, coherence, condition number, and "
                    "parameter-stability diagnostics alongside this calibrator"
                ),
            },
        )

    def apply(
        self,
        predictions: Mapping[str, Mapping[str, Any]],
    ) -> dict[str, dict[str, Any]]:
        """Calibrate and re-apply fixed abstention rules without truth."""

        _walk_forbidden_keys(predictions)
        output: dict[str, dict[str, Any]] = {}
        for edge_id, raw in predictions.items():
            if not isinstance(raw, Mapping):
                raise TypeError(f"Reaction prediction {edge_id!r} must be a mapping")
            logits = raw.get("logits")
            input_kind = "logits"
            if not isinstance(logits, Mapping) or not logits:
                probabilities = raw.get("probabilities")
                if isinstance(probabilities, Mapping) and probabilities:
                    logits = {
                        str(key): math.log(max(_finite(value, name="probability"), 1.0e-12))
                        for key, value in probabilities.items()
                    }
                    input_kind = "probabilities_as_logits"
                else:
                    updated = dict(raw)
                    updated.update(
                        {
                            "family": NULL_FAMILY,
                            "candidate_family": None,
                            "decision": ABSTAIN,
                            "calibration_status": "abstained_not_calibrated",
                            "calibration_reason": "missing_probability_or_logit_vector",
                        }
                    )
                    output[str(edge_id)] = updated
                    continue
            probabilities = _mapping_softmax(
                logits,
                self.classes,
                temperature=self.temperature,
                class_offsets={
                    family: self.class_offsets[index]
                    for index, family in enumerate(self.classes)
                },
            )
            ranking = sorted(self.classes, key=lambda family: (-probabilities[family], family))
            best = ranking[0]
            second = probabilities[ranking[1]] if len(ranking) > 1 else 0.0
            if probabilities[best] < self.decision_threshold:
                decision = ABSTAIN
                reason = "calibrated_top_probability_below_threshold"
            elif probabilities[best] - second < self.decision_margin:
                decision = ABSTAIN
                reason = "calibrated_family_margin_below_threshold"
            else:
                decision = SELECT
                reason = (
                    "calibrated_null_or_mixing_selected"
                    if best in {NULL_FAMILY, "mixing"}
                    else "calibrated_reaction_family_selected"
                )
            updated = dict(raw)
            updated.update(
                {
                    "family": best if decision == SELECT else NULL_FAMILY,
                    "candidate_family": best,
                    "probability": probabilities[best],
                    "probabilities": probabilities,
                    "decision": decision,
                    "reason": reason,
                    "calibration_status": "development_fitted",
                    "calibration_temperature": self.temperature,
                    "calibration_input": input_kind,
                }
            )
            uncertainty = dict(raw.get("uncertainty", {}))
            uncertainty.update(
                {
                    "calibrated": True,
                    "selection_probability_calibrated": True,
                    "method": self.calibration_kind,
                }
            )
            updated["uncertainty"] = uncertainty
            diagnostics = dict(raw.get("diagnostics", {}))
            diagnostics["calibration_temperature"] = self.temperature
            diagnostics["calibration_class_offsets"] = dict(
                zip(self.classes, self.class_offsets, strict=True)
            )
            diagnostics["calibration_threshold_fixed"] = self.decision_threshold
            diagnostics["calibration_margin_fixed"] = self.decision_margin
            updated["diagnostics"] = diagnostics
            if bool(self.provenance.get("grid_boundary_hit")):
                updated["calibration_warning"] = "temperature_grid_boundary_hit"
            output[str(edge_id)] = updated
        return output

    def to_audit_record(self) -> dict[str, Any]:
        record: dict[str, Any] = {
            "name": "reaction_rapm_calibrator",
            "version": "0.1.0",
            "fit_scope": self.fit_scope,
            "calibration_kind": self.calibration_kind,
            "temperature": self.temperature,
            "classes": list(self.classes),
            "decision_threshold": self.decision_threshold,
            "decision_margin": self.decision_margin,
            "class_offsets": list(self.class_offsets),
            "fit_count": self.fit_count,
            "case_ids": list(self.case_ids),
            "provenance": dict(self.provenance),
            "control": {
                "truth_blind_apply": True,
                "truth_used_for": "development_fit_only",
                "threshold_tuned_on_locked_truth": False,
                "uses_phreeqc_inverse": False,
            },
        }
        record["fingerprint"] = _fingerprint(record)
        return record

    def to_dict(self) -> dict[str, Any]:
        return self.to_audit_record()


def cross_fitted_reaction_rapm_calibration_examples(
    examples: Iterable[ReactionRAPMTrainingExample],
    *,
    config: ReactionRAPMConfig | None = None,
) -> tuple[ReactionRAPMCalibrationExample, ...]:
    """Generate case-held-out logits for calibration on development cases.

    At least three development cases are required: each held-out case needs a
    training set containing two or more cases so the RAPM fitting contract is
    preserved.  The returned truth labels are fit-only records and must not be
    passed to locked prediction.
    """

    config = config or ReactionRAPMConfig()
    source = tuple(examples)
    case_ids = tuple(sorted({str(example.case_id) for example in source}))
    if len(case_ids) < 3:
        raise ValueError("Cross-fitted reaction calibration requires at least three cases")
    result: list[ReactionRAPMCalibrationExample] = []
    for held_out in case_ids:
        train = tuple(example for example in source if str(example.case_id) != held_out)
        model = ReactionRAPM.fit(train, config=config, phase="development")
        for example in source:
            if str(example.case_id) != held_out:
                continue
            row = _make_feature_row(
                case_id=held_out,
                edge_id=str(example.edge_id),
                upstream=example.upstream,
                downstream=example.downstream,
                family=NULL_FAMILY,
                config=config,
            )
            if row.status is not None:
                continue
            prediction = model._predict_row(row)
            result.append(
                ReactionRAPMCalibrationExample(
                    case_id=held_out,
                    edge_id=str(example.edge_id),
                    truth_family=example.truth_family,
                    logits=prediction["logits"],
                    cross_fitted=True,
                )
            )
    if not result:
        raise ValueError("No valid case-held-out reaction calibration examples were generated")
    return tuple(result)


def training_examples_from_observations(
    observations: Mapping[str, Any],
    truth_by_edge: Mapping[str, object],
    *,
    case_id: str,
) -> tuple[ReactionRAPMTrainingExample, ...]:
    """Build labelled development examples from a blind observation mapping.

    The truth mapping is a separate fit-only argument.  It must never be
    passed to :meth:`ReactionRAPM.predict`; this separation is what permits
    locked inference to remain truth-blind.
    """

    _walk_forbidden_keys(observations)
    examples: list[ReactionRAPMTrainingExample] = []
    for edge_id, edge in _edge_items(observations):
        upstream = _side(edge, "upstream", "source")
        downstream = _side(edge, "downstream", "target_values")
        if downstream is None:
            downstream = _side(edge, "downstream", "destination")
        if upstream is None or downstream is None:
            continue
        examples.append(
            ReactionRAPMTrainingExample(
                case_id=str(case_id),
                edge_id=edge_id,
                upstream=dict(upstream),
                downstream=dict(downstream),
                truth_family=str(truth_by_edge.get(edge_id, NULL_FAMILY)),
            )
        )
    return tuple(examples)


def _prediction_probabilities(
    row: Mapping[str, Any],
) -> dict[str, float] | None:
    raw = row.get("probabilities")
    if not isinstance(raw, Mapping) or not raw:
        return None
    numeric: dict[str, float] = {}
    for key, value in raw.items():
        finite = _optional_finite(value)
        if finite is not None and finite >= 0.0:
            numeric[str(key)] = finite
    total = sum(numeric.values())
    if not numeric or total <= 0.0:
        return None
    return {key: value / total for key, value in numeric.items()}


def _classwise_reliability(
    records: Sequence[Mapping[str, Any]],
    classes: Sequence[str],
    *,
    n_bins: int = 10,
) -> dict[str, Any]:
    """Return one-vs-rest reliability bins for every reaction family."""

    result: dict[str, Any] = {}
    for family in classes:
        members = [
            (
                float(record["probabilities"].get(family, 0.0)),
                float(record["expected"] == family),
            )
            for record in records
            if isinstance(record.get("probabilities"), Mapping)
        ]
        bins: list[dict[str, Any]] = []
        weighted_error = 0.0
        for index in range(n_bins):
            lower = index / n_bins
            upper = (index + 1) / n_bins
            selected = [
                (probability, observed)
                for probability, observed in members
                if lower <= probability < upper
                or (index == n_bins - 1 and probability <= upper)
            ]
            if not selected:
                continue
            mean_probability = float(np.mean([item[0] for item in selected]))
            observed_frequency = float(np.mean([item[1] for item in selected]))
            gap = abs(mean_probability - observed_frequency)
            weighted_error += len(selected) / max(len(members), 1) * gap
            bins.append(
                {
                    "lower": lower,
                    "upper": upper,
                    "n": len(selected),
                    "mean_probability": mean_probability,
                    "observed_frequency": observed_frequency,
                    "absolute_gap": gap,
                }
            )
        result[family] = {
            "n": len(members),
            "ece": weighted_error if members else float("nan"),
            "bins": bins,
        }
    return result


def _partition_performance(
    records: Sequence[Mapping[str, Any]],
    *,
    partition: str,
) -> dict[str, Any]:
    selected = [record for record in records if record.get("partition") == partition]
    accepted = [record for record in selected if bool(record.get("accepted"))]
    correct = [record for record in accepted if record.get("predicted") == record.get("expected")]
    if partition == "decoy":
        false_commitments = [
            record
            for record in accepted
            if str(record.get("predicted")) != NULL_FAMILY
        ]
    else:
        false_commitments = [record for record in accepted if not record.get("correct")]
    return {
        "n": len(selected),
        "n_accepted": len(accepted),
        "coverage": len(accepted) / len(selected) if selected else float("nan"),
        "selective_accuracy": len(correct) / len(accepted) if accepted else float("nan"),
        "selective_risk": 1.0 - len(correct) / len(accepted) if accepted else float("nan"),
        "false_commitment_count": len(false_commitments),
        "false_commitment_rate": (
            len(false_commitments) / len(selected) if selected else float("nan")
        ),
        "non_null_commitment_rate": (
            sum(
                str(record.get("predicted")) != NULL_FAMILY
                for record in accepted
            )
            / len(selected)
            if selected
            else float("nan")
        ),
    }


def _mixing_none_equivalence(
    records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    equivalent = {NULL_FAMILY, "mixing"}
    accepted = [record for record in records if bool(record.get("accepted"))]
    correct = [
        record
        for record in accepted
        if (record.get("expected") in equivalent)
        == (record.get("predicted") in equivalent)
    ]
    raw_confusion = {
        "truth_none_pred_mixing": sum(
            record.get("expected") == NULL_FAMILY
            and record.get("predicted") == "mixing"
            and bool(record.get("accepted"))
            for record in records
        ),
        "truth_mixing_pred_none": sum(
            record.get("expected") == "mixing"
            and record.get("predicted") == NULL_FAMILY
            and bool(record.get("accepted"))
            for record in records
        ),
    }
    return {
        "equivalence_class": "none_or_mixing",
        "n": len(records),
        "n_accepted": len(accepted),
        "coverage": len(accepted) / len(records) if records else float("nan"),
        "equivalent_correct": len(correct),
        "selective_accuracy": len(correct) / len(accepted) if accepted else float("nan"),
        "selective_risk": 1.0 - len(correct) / len(accepted) if accepted else float("nan"),
        "raw_confusion": raw_confusion,
        "interpretation": (
            "none and mixing are reported separately in the primary taxonomy but "
            "also scored as one observationally equivalent no-chemical-change class"
        ),
    }


def score_reaction_rapm_outputs(
    true_processes: Mapping[str, object],
    predictions: Mapping[str, Mapping[str, object]] | None,
    *,
    candidate_edge_ids: Iterable[str] | None = None,
    model_audit: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Score RAPM predictions after truth release on the declared edge universe.

    When ``candidate_edge_ids`` is supplied, edges absent from ``true_processes``
    are explicit no-reaction decoys.  This prevents the null class and
    abstention from being evaluated only on true labelled edges.  When
    ``model_audit`` is supplied, its fitted rank/coherence/conditioning and
    stability diagnostics are copied into the performance record.
    """

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
    reference_ids = set(candidate_ids).intersection(truth_by_id)
    decoy_ids = set(candidate_ids).difference(truth_by_id)

    if not predictions:
        return {
            "status": "not_available",
            "n": 0,
            "n_truth": len(candidate_ids),
            "n_reference_edges": len(reference_ids),
            "n_decoy_edges": len(decoy_ids),
            "n_missing_outputs": len(candidate_ids),
        }
    expected: list[str] = []
    predicted: list[str] = []
    log_losses: list[float] = []
    brier_terms: list[float] = []
    n_abstain = 0
    n_missing = 0
    accepted_correct = 0
    unconditional_correct = 0
    false_commitments = 0
    family_summary: dict[str, dict[str, int]] = {}
    scored_records: list[dict[str, Any]] = []
    observed_classes: set[str] = set(DEFAULT_REACTION_CLASSES)
    for edge_id in candidate_ids:
        truth = truth_by_id.get(edge_id, NULL_FAMILY)
        is_reference = edge_id in reference_ids
        row = predictions.get(str(edge_id))
        if not isinstance(row, Mapping):
            n_missing += 1
            expected_family = _normalise_family(truth)
            expected.append(expected_family)
            predicted.append(NULL_FAMILY)
            summary = family_summary.setdefault(
                expected_family,
                {"n": 0, "n_accepted": 0, "n_correct": 0},
            )
            summary["n"] += 1
            scored_records.append(
                {
                    "edge_id": edge_id,
                    "expected": expected_family,
                    "predicted": NULL_FAMILY,
                    "accepted": False,
                    "correct": False,
                    "missing_output": True,
                    "partition": "reference" if is_reference else "decoy",
                    "probabilities": {},
                    "missing_fields": (),
                    "invalid_fields": (),
                    "unsupported_fields": (),
                    "n_paired_fields": 0,
                }
            )
            continue
        expected_family = _normalise_family(truth)
        predicted_family = str(row.get("family", NULL_FAMILY))
        decision = str(row.get("decision", ABSTAIN))
        accepted = decision != ABSTAIN
        probabilities = _prediction_probabilities(row)
        if probabilities is not None:
            observed_classes.update(probabilities)
        expected.append(expected_family)
        predicted.append(predicted_family)
        summary = family_summary.setdefault(
            expected_family,
            {"n": 0, "n_accepted": 0, "n_correct": 0},
        )
        summary["n"] += 1
        if decision == ABSTAIN:
            n_abstain += 1
        if accepted:
            summary["n_accepted"] += 1
            if predicted_family == expected_family:
                accepted_correct += 1
                unconditional_correct += 1
                summary["n_correct"] += 1
            if not is_reference and predicted_family != NULL_FAMILY:
                false_commitments += 1
        if probabilities is not None:
            p_true = probabilities.get(expected_family, 0.0)
            log_losses.append(-math.log(max(p_true, 1.0e-12)))
            brier_terms.append(
                sum(
                    (probabilities.get(key, 0.0) - float(key == expected_family)) ** 2
                    for key in set(probabilities) | {expected_family}
                )
            )
        diagnostics = row.get("diagnostics", {})
        diagnostics = diagnostics if isinstance(diagnostics, Mapping) else {}
        scored_records.append(
            {
                "edge_id": edge_id,
                "expected": expected_family,
                "predicted": predicted_family,
                "accepted": accepted,
                "correct": predicted_family == expected_family,
                "partition": "reference" if is_reference else "decoy",
                "probabilities": probabilities or {},
                "missing_fields": tuple(diagnostics.get("missing_fields", ()))
                if isinstance(diagnostics.get("missing_fields", ()), Sequence)
                and not isinstance(diagnostics.get("missing_fields", ()), (str, bytes, bytearray))
                else (),
                "invalid_fields": tuple(diagnostics.get("invalid_fields", ()))
                if isinstance(diagnostics.get("invalid_fields", ()), Sequence)
                and not isinstance(diagnostics.get("invalid_fields", ()), (str, bytes, bytearray))
                else (),
                "unsupported_fields": tuple(diagnostics.get("unsupported_fields", ()))
                if isinstance(diagnostics.get("unsupported_fields", ()), Sequence)
                and not isinstance(diagnostics.get("unsupported_fields", ()), (str, bytes, bytearray))
                else (),
                "n_paired_fields": int(diagnostics.get("n_paired_fields", 0) or 0),
            }
        )
    if not expected:
        return {
            "status": "no_comparable_outputs",
            "n": 0,
            "n_truth": len(candidate_ids),
            "n_reference_edges": len(reference_ids),
            "n_decoy_edges": len(decoy_ids),
            "n_missing_outputs": n_missing,
        }
    available = len(candidate_ids)
    accepted = sum(1 for record in scored_records if bool(record.get("accepted")))
    missingness_records = [
        record
        for record in scored_records
        if record["missing_fields"] or record["invalid_fields"] or record["unsupported_fields"]
    ]
    by_paired_fields: dict[str, int] = {}
    for record in scored_records:
        key = str(record["n_paired_fields"])
        by_paired_fields[key] = by_paired_fields.get(key, 0) + 1
    identifiability = {
        "status": "not_supplied",
        "interpretation": (
            "pass the fitted RAPM audit record to attach the actual rank, "
            "coherence, conditioning, and parameter-stability diagnostics"
        ),
    }
    if isinstance(model_audit, Mapping):
        supplied_identifiability = model_audit.get("identifiability")
        if isinstance(supplied_identifiability, Mapping):
            identifiability = dict(supplied_identifiability)
            identifiability["status"] = "supplied_from_model_audit"
    return {
        "status": "scored",
        "n": available,
        "n_truth": len(candidate_ids),
        "n_reference_edges": len(reference_ids),
        "n_decoy_edges": len(decoy_ids),
        "n_missing_outputs": n_missing,
        "n_abstain": n_abstain,
        "n_accepted": accepted,
        "coverage": accepted / available if available else float("nan"),
        "coverage_including_abstention": (
            accepted / len(candidate_ids)
            if candidate_ids
            else float("nan")
        ),
        "selective_accuracy": (
            accepted_correct / accepted if accepted else float("nan")
        ),
        "selective_risk": (
            1.0 - accepted_correct / accepted if accepted else float("nan")
        ),
        "accuracy": unconditional_correct / len(candidate_ids),
        "unconditional_error": 1.0 - unconditional_correct / len(candidate_ids),
        "false_commitment_rate": (
            false_commitments / len(decoy_ids) if decoy_ids else float("nan")
        ),
        "abstention_rate": (n_abstain + n_missing) / available if available else float("nan"),
        "available_output_count": available - n_missing,
        "outputs_complete": n_missing == 0,
        "multiclass_log_loss": (
            sum(log_losses) / len(log_losses) if log_losses else float("nan")
        ),
        "multiclass_brier": (
            sum(brier_terms) / len(brier_terms) if brier_terms else float("nan")
        ),
        "unresolved_rate": predicted.count(NULL_FAMILY) / len(predicted),
        "by_expected_family": {
            family: {
                **values,
                "coverage": (
                    values["n_accepted"] / values["n"] if values["n"] else float("nan")
                ),
                "selective_accuracy": (
                    values["n_correct"] / values["n_accepted"]
                    if values["n_accepted"]
                    else float("nan")
                ),
            }
            for family, values in sorted(family_summary.items())
        },
        "expected_families": sorted(set(expected)),
        "predicted_families": sorted(set(predicted)),
        "reference_performance": _partition_performance(
            scored_records, partition="reference"
        ),
        "decoy_performance": _partition_performance(
            scored_records, partition="decoy"
        ),
        "classwise_reliability": _classwise_reliability(
            scored_records, tuple(sorted(observed_classes))
        ),
        "missingness": {
            "n_with_input_missing_or_invalid_fields": len(missingness_records),
            "rate_with_input_missing_or_invalid_fields": (
                len(missingness_records) / available if available else float("nan")
            ),
            "n_missing_outputs": n_missing,
            "by_paired_field_count": dict(sorted(by_paired_fields.items())),
            "interpretation": (
                "missing output is separated from incomplete chemistry; neither is "
                "silently scored as a correct abstention"
            ),
        },
        "mixing_none_equivalence": _mixing_none_equivalence(scored_records),
        "identifiability_diagnostics_contract": {
            "required_from_model_audit": [
                "design_rank",
                "design_singular_values",
                "design_condition_number",
                "maximum_absolute_feature_coherence",
                "parameter_stability_across_case_folds",
            ],
            "interpretation": (
                "performance metrics do not establish unique chemical mechanism "
                "identification when signatures are rank-deficient or coherent"
            ),
        },
        "identifiability": identifiability,
    }


def fit_reaction_rapm(
    examples: Iterable[ReactionRAPMTrainingExample],
    *,
    config: ReactionRAPMConfig | None = None,
    phase: str = "development",
) -> ReactionRAPM:
    return ReactionRAPM.fit(examples, config=config, phase=phase)


def fit_reaction_rapm_calibrator(
    examples: Iterable[ReactionRAPMCalibrationExample],
    *,
    classes: Sequence[str] = DEFAULT_REACTION_CLASSES,
    decision_threshold: float = 0.60,
    decision_margin: float = 0.10,
    tune_selection_rule: bool = False,
    target_coverage: float = 0.25,
    max_selective_risk: float = 0.40,
    phase: str = "development",
) -> ReactionRAPMCalibrator:
    return ReactionRAPMCalibrator.fit(
        examples,
        classes=classes,
        decision_threshold=decision_threshold,
        decision_margin=decision_margin,
        tune_selection_rule=tune_selection_rule,
        target_coverage=target_coverage,
        max_selective_risk=max_selective_risk,
        phase=phase,
    )


__all__ = [
    "ABSTAIN",
    "DEFAULT_REACTION_CLASSES",
    "DEFAULT_ION_SCALES",
    "EVIDENCE_CHANNEL",
    "NULL_FAMILY",
    "REACTION_FIELDS",
    "REACTION_OUTPUT",
    "ReactionRAPM",
    "ReactionRAPMCalibrationExample",
    "ReactionRAPMCalibrator",
    "ReactionRAPMConfig",
    "ReactionRAPMTrainingExample",
    "SELECT",
    "cross_fitted_reaction_rapm_calibration_examples",
    "fit_reaction_rapm",
    "fit_reaction_rapm_calibrator",
    "normalise_rapm_reaction_family",
    "score_reaction_rapm_outputs",
    "training_examples_from_observations",
]
