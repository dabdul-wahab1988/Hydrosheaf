"""Null-aware, calibrated topology evidence.

This module is intentionally disjoint from the existing topology posterior
implementation.  It provides a small, auditable component for comparing
flow-supporting evidence with evidence that is also plausible under a null
(no-flow) explanation.

Only explicitly allow-listed fields are read from observation and edge
metadata.  In particular, labels, truth/reference fields, and arbitrary
numeric metadata are never used as features.  A row can therefore be built
before a synthetic truth record is attached without changing its features.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
import math
import re
from typing import Any, Optional

import numpy as np
from scipy.optimize import minimize

from ..config import Config
from ..null_models import compute_null_penalty
from ..null_models.chemistry import chemistry_null_score
from ..null_models.endmembers import endmember_null_score
from ..null_models.lithology import lithology_null_score


FLOW_FEATURE_NAMES = (
    "head_drop_m",
    "age_increment_years",
    "age_direction_support",
    "chemistry_consistency",
    "isotope_consistency",
    "transport_support",
)

NULL_FEATURE_NAMES = (
    "null_chemistry_similarity",
    "null_isotope_similarity",
    "common_lithology",
    "shared_recharge",
    "spatial_proximity",
    "common_source",
    "null_score",
)

# ``null_score`` is retained as an auditable derived field for compatibility,
# but it is not an independent observation channel.  Including it in the
# logistic design would double-count the component null evidence used to
# construct it.
_MODEL_NULL_FEATURE_NAMES = tuple(
    name for name in NULL_FEATURE_NAMES if name != "null_score"
)

DECISION_PRESENT = "PRESENT"
DECISION_ABSENT = "ABSENT"
DECISION_ABSTAIN = "ABSTAIN"

CALIBRATION_UNFITTED = "UNFITTED"
CALIBRATION_FITTED = "FITTED"
CALIBRATION_FITTED_MISSING_EVIDENCE = "FITTED_MISSING_EVIDENCE"

_DEPLOYABLE_CALIBRATION_SCOPES = frozenset(
    {"held_out_calibration", "deployment"}
)

_ALL_FEATURE_NAMES = FLOW_FEATURE_NAMES + _MODEL_NULL_FEATURE_NAMES

# These names are deliberately conservative.  A field is rejected if its
# normalized name contains one of these truth/label markers, even when it is
# nested inside an otherwise accepted metadata mapping.
_TRUTH_FIELD_NAMES = {
    "actual",
    "connectivity_truth",
    "flow_truth",
    "ground_truth",
    "is_connected",
    "is_true_flow",
    "label",
    "observed_present",
    "present",
    "reference",
    "reference_flow",
    "target",
    "true_flow",
    "truth",
}
_TRUTH_FIELD_PATTERN = re.compile(
    r"(?:^|[_-])(actual|truth|true|label|target|reference|ground)(?:$|[_-])"
)

_FLOW_ALIASES = {
    "head_drop_m": (
        "head_drop_m",
        "delta_head",
        "delta_h",
        "hydraulic_head_drop",
        "head_gradient",
    ),
    "age_increment_years": (
        "age_increment_years",
        "age_delta_years",
        "delta_age_years",
        "travel_time_years",
    ),
    "age_direction_support": (
        "age_direction_support",
        "age_consistency",
        "age_support",
    ),
    "chemistry_consistency": (
        "chemistry_consistency",
        "chemistry_support",
        "chemistry_score",
    ),
    "isotope_consistency": (
        "isotope_consistency",
        "isotope_support",
        "isotope_score",
    ),
    "transport_support": (
        "transport_support",
        "transport_score",
        "flow_support",
        "flow_score",
    ),
}

_NULL_ALIASES = {
    "null_chemistry_similarity": (
        "null_chemistry_similarity",
        "chemistry_similarity",
        "null_chemistry_score",
    ),
    "null_isotope_similarity": (
        "null_isotope_similarity",
        "isotope_similarity",
        "null_endmember_score",
    ),
    "common_lithology": (
        "common_lithology",
        "lithology_similarity",
        "null_lithology_score",
    ),
    "shared_recharge": (
        "shared_recharge",
        "shared_recharge_score",
        "common_recharge",
    ),
    "spatial_proximity": (
        "spatial_proximity",
        "spatial_similarity",
        "null_spatial_score",
    ),
    "common_source": (
        "common_source",
        "common_source_score",
        "null_anthropogenic_score",
    ),
    "null_score": (
        "null_score",
        "no_flow_score",
        "null_model_score",
    ),
}

_OBSERVATION_ID_KEYS = ("site_id", "node_id", "sample_id", "id")
_UPSTREAM_KEYS = ("u", "upstream", "source", "from", "from_node")
_DOWNSTREAM_KEYS = ("v", "downstream", "target_node", "to", "to_node")
_HEAD_KEYS = ("head_meas", "hydraulic_head", "head", "water_level")
_AGE_KEYS = ("mean_age_years", "age_years", "age")
_LITHOLOGY_KEYS = ("aquifer_layer", "aquifer_unit", "lithology")


def _normalise_key(value: Any) -> str:
    """Return a stable lower-case metadata key representation."""

    return str(value).strip().lower().replace(" ", "_")


def _is_truth_like(value: Any) -> bool:
    key = _normalise_key(value)
    return key in _TRUTH_FIELD_NAMES or bool(_TRUTH_FIELD_PATTERN.search(key))


def _as_finite_float(value: Any) -> Optional[float]:
    if isinstance(value, (bool, np.bool_)):
        return float(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        number = float(value)
        return number if math.isfinite(number) else None
    return None


def _as_mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(values))


def _metadata_sources(metadata: Mapping[str, Any], scope: str) -> list[tuple[str, Mapping[str, Any]]]:
    """Return accepted metadata namespaces in deterministic precedence order."""

    sources: list[tuple[str, Mapping[str, Any]]] = []
    nested_key = f"{scope}_features"
    nested = metadata.get(nested_key)
    if isinstance(nested, Mapping):
        sources.append((nested_key, nested))

    # A generic feature mapping is accepted only for explicitly named keys.
    generic = metadata.get("features")
    if isinstance(generic, Mapping):
        sources.append(("features", generic))

    sources.append(("edge", metadata))
    return sources


def _lookup_explicit(
    metadata: Mapping[str, Any],
    aliases: Sequence[str],
    scope: str,
) -> tuple[Optional[float], Optional[str]]:
    aliases_normalized = {_normalise_key(alias) for alias in aliases}
    for source_name, source in _metadata_sources(metadata, scope):
        for key, value in source.items():
            normalized = _normalise_key(key)
            if normalized not in aliases_normalized or _is_truth_like(normalized):
                continue
            number = _as_finite_float(value)
            if number is not None:
                return number, f"{source_name}.{key}"
    return None, None


def _lookup_value(metadata: Mapping[str, Any], keys: Sequence[str]) -> tuple[Any, Optional[str]]:
    wanted = {_normalise_key(key) for key in keys}
    for key, value in metadata.items():
        normalized = _normalise_key(key)
        if normalized in wanted and not _is_truth_like(normalized):
            return value, str(key)
    return None, None


def _edge_parts(edge: Any, index: int) -> tuple[str, Any, Any, Mapping[str, Any]]:
    if isinstance(edge, Mapping):
        edge_id = edge.get("edge_id", edge.get("id", f"edge-{index}"))
        upstream = next((edge[key] for key in _UPSTREAM_KEYS if key in edge), None)
        downstream = next((edge[key] for key in _DOWNSTREAM_KEYS if key in edge), None)
        attrs = dict(edge)
        nested_attrs = edge.get("attrs", edge.get("metadata"))
        if isinstance(nested_attrs, Mapping):
            attrs.update(nested_attrs)
        return str(edge_id), upstream, downstream, attrs

    if isinstance(edge, (tuple, list)):
        if len(edge) == 2:
            upstream, downstream = edge
            return f"{upstream}->{downstream}", upstream, downstream, {}
        if len(edge) >= 3:
            edge_id, upstream, downstream = edge[:3]
            attrs = edge[3] if len(edge) > 3 and isinstance(edge[3], Mapping) else {}
            return str(edge_id), upstream, downstream, attrs

    edge_id = getattr(edge, "edge_id", getattr(edge, "id", f"edge-{index}"))
    upstream = next((getattr(edge, key, None) for key in _UPSTREAM_KEYS if hasattr(edge, key)), None)
    downstream = next((getattr(edge, key, None) for key in _DOWNSTREAM_KEYS if hasattr(edge, key)), None)
    attrs = getattr(edge, "attrs", getattr(edge, "metadata", {}))
    return str(edge_id), upstream, downstream, _as_mapping(attrs)


def _observation_table(observations: Any) -> dict[str, Mapping[str, Any]]:
    if observations is None:
        return {}
    if isinstance(observations, Mapping):
        table: dict[str, Mapping[str, Any]] = {}
        for key, value in observations.items():
            if isinstance(value, Mapping):
                table[str(key)] = value
        return table

    table = {}
    if isinstance(observations, (str, bytes)):
        return table
    for index, observation in enumerate(observations):
        if not isinstance(observation, Mapping):
            continue
        identifier = next(
            (observation[key] for key in _OBSERVATION_ID_KEYS if key in observation),
            f"observation-{index}",
        )
        table[str(identifier)] = observation
    return table


def _observation_for(table: Mapping[str, Mapping[str, Any]], identifier: Any) -> Mapping[str, Any]:
    if identifier is None:
        return {}
    return table.get(str(identifier), {})


def _derive_head_drop(
    upstream: Mapping[str, Any],
    downstream: Mapping[str, Any],
    keys: Sequence[str],
) -> tuple[Optional[float], Optional[str]]:
    upstream_value, upstream_key = _lookup_value(upstream, keys)
    downstream_value, downstream_key = _lookup_value(downstream, keys)
    upstream_number = _as_finite_float(upstream_value)
    downstream_number = _as_finite_float(downstream_value)
    if upstream_number is None or downstream_number is None:
        return None, None
    return upstream_number - downstream_number, f"observation.{upstream_key}-{downstream_key}"


def _derive_age_increment(
    upstream: Mapping[str, Any],
    downstream: Mapping[str, Any],
    keys: Sequence[str],
) -> tuple[Optional[float], Optional[str]]:
    upstream_value, upstream_key = _lookup_value(upstream, keys)
    downstream_value, downstream_key = _lookup_value(downstream, keys)
    upstream_number = _as_finite_float(upstream_value)
    downstream_number = _as_finite_float(downstream_value)
    if upstream_number is None or downstream_number is None:
        return None, None
    return downstream_number - upstream_number, f"observation.{downstream_key}-{upstream_key}"


def _derive_common_lithology(
    upstream: Mapping[str, Any],
    downstream: Mapping[str, Any],
) -> tuple[Optional[float], Optional[str]]:
    for key in _LITHOLOGY_KEYS:
        upstream_value, upstream_key = _lookup_value(upstream, (key,))
        downstream_value, downstream_key = _lookup_value(downstream, (key,))
        if upstream_value is not None and downstream_value is not None:
            return float(upstream_value == downstream_value), f"observation.{downstream_key}=={upstream_key}"
    return None, None


def _paired_numeric_count(
    upstream: Mapping[str, Any],
    downstream: Mapping[str, Any],
    keys: Sequence[str],
) -> int:
    count = 0
    for key in keys:
        first, _ = _lookup_value(upstream, (key,))
        second, _ = _lookup_value(downstream, (key,))
        if _as_finite_float(first) is not None and _as_finite_float(second) is not None:
            count += 1
    return count


def _derive_null_features(
    upstream: Mapping[str, Any],
    downstream: Mapping[str, Any],
    config: Config,
) -> tuple[dict[str, float], list[str]]:
    """Derive null features from declared observations only.

    Existing heuristic null models remain useful as feature generators, but
    their outputs are not treated as calibrated probabilities.  The fitted
    logistic layer performs the held-out calibration and records its scope.
    """
    null: dict[str, float] = {}
    source_fields: list[str] = []

    chemistry_count = _paired_numeric_count(
        upstream,
        downstream,
        tuple(getattr(config, "ion_order", ())),
    )
    if chemistry_count >= 2:
        chemistry_score, _ = chemistry_null_score(upstream, downstream, config)
        null["null_chemistry_similarity"] = float(chemistry_score)
        source_fields.append("null_model.chemistry")

    isotope_keys = (
        getattr(config, "isotope_d18o_key", "18O"),
        getattr(config, "isotope_d2h_key", "2H"),
    )
    if any(_is_truth_like(key) for key in isotope_keys):
        raise ValueError("configured isotope keys must be truth-blind observation fields")
    if _paired_numeric_count(upstream, downstream, isotope_keys) == 2:
        d18o_u, _ = _lookup_value(upstream, (isotope_keys[0],))
        d18o_v, _ = _lookup_value(downstream, (isotope_keys[0],))
        d2h_u, _ = _lookup_value(upstream, (isotope_keys[1],))
        d2h_v, _ = _lookup_value(downstream, (isotope_keys[1],))
        d18o_u = _as_finite_float(d18o_u)
        d18o_v = _as_finite_float(d18o_v)
        d2h_u = _as_finite_float(d2h_u)
        d2h_v = _as_finite_float(d2h_v)
        assert d18o_u is not None and d18o_v is not None
        assert d2h_u is not None and d2h_v is not None
        isotope_distance = math.sqrt(
            ((d18o_u - d18o_v) / 1.0) ** 2
            + ((d2h_u - d2h_v) / 8.0) ** 2
        )
        null["null_isotope_similarity"] = float(
            math.exp(-0.5 * isotope_distance**2)
        )
        endmember_score, endmember_flags = endmember_null_score(
            upstream, downstream, config
        )
        null["shared_recharge"] = float(
            "null_shared_recharge" in endmember_flags
        )
        source_fields.append("null_model.isotopes")
        if endmember_score > 0.0:
            source_fields.append("null_model.endmembers")

    lithology_score, _ = lithology_null_score(upstream, downstream, config)
    if (
        _paired_numeric_count(
            upstream,
            downstream,
            (getattr(config, "layer_key", "aquifer_layer"),),
        ) > 0
        or any(
            upstream.get(key) is not None and downstream.get(key) is not None
            for key in _LITHOLOGY_KEYS
        )
    ):
        null["common_lithology"] = float(lithology_score)
        source_fields.append("null_model.lithology")

    lat_u = _as_finite_float(upstream.get("lat"))
    lon_u = _as_finite_float(upstream.get("lon"))
    lat_v = _as_finite_float(downstream.get("lat"))
    lon_v = _as_finite_float(downstream.get("lon"))
    if None not in (lat_u, lon_u, lat_v, lon_v):
        distance_degrees = math.hypot(lat_u - lat_v, lon_u - lon_v)  # type: ignore[operator]
        null["spatial_proximity"] = float(math.exp(-distance_degrees / 0.01))
        source_fields.append("observation.lat_lon")

    if null:
        combined_score, combined_flags = compute_null_penalty(
            upstream, downstream, config
        )
        null["null_score"] = float(combined_score)
        if "null_common_anthropogenic" in combined_flags:
            null["common_source"] = 1.0
        source_fields.append("null_model.combined")

    return null, source_fields


@dataclass(frozen=True)
class NullAwareFeatureRow:
    """One metadata-derived candidate-edge feature row.

    Feature values are separated by evidential role.  ``missing_channels``
    contains names such as ``flow:head_drop_m`` and
    ``null:null_chemistry_similarity`` so an abstention can be audited without
    inspecting the model internals.
    """

    edge_id: str
    flow_features: Mapping[str, float] = field(default_factory=dict)
    null_features: Mapping[str, float] = field(default_factory=dict)
    missing_channels: tuple[str, ...] = ()
    source_fields: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        flow = {
            name: float(value)
            for name, value in self.flow_features.items()
            if name in FLOW_FEATURE_NAMES and _as_finite_float(value) is not None
        }
        null = {
            name: float(value)
            for name, value in self.null_features.items()
            if name in NULL_FEATURE_NAMES and _as_finite_float(value) is not None
        }
        missing = list(self.missing_channels)
        missing.extend(f"flow:{name}" for name in FLOW_FEATURE_NAMES if name not in flow)
        missing.extend(f"null:{name}" for name in NULL_FEATURE_NAMES if name not in null)
        object.__setattr__(self, "flow_features", flow)
        object.__setattr__(self, "null_features", null)
        object.__setattr__(self, "missing_channels", _unique(missing))
        object.__setattr__(self, "source_fields", _unique(self.source_fields))

    def to_dict(self) -> dict[str, Any]:
        return {
            "edge_id": self.edge_id,
            "flow_features": dict(self.flow_features),
            "null_features": dict(self.null_features),
            "missing_channels": list(self.missing_channels),
            "source_fields": list(self.source_fields),
        }


# A short alias is useful for callers that do not need the longer class name.
FeatureRow = NullAwareFeatureRow


def build_feature_rows(
    candidate_edges: Iterable[Any],
    observations: Any = None,
    config: Optional[Config] = None,
) -> list[NullAwareFeatureRow]:
    """Build explicit flow/null feature rows for every candidate edge.

    ``candidate_edges`` may contain mappings, ``(u, v)`` pairs,
    ``(edge_id, u, v[, attrs])`` tuples, or objects with ``edge_id``, ``u`` /``v``
    and optional ``attrs``/``metadata`` attributes.  Observation-derived
    differences are used only for the explicitly named head and age fields.
    No truth-like or arbitrary numeric fields are inspected.
    """

    table = _observation_table(observations)
    feature_config = config or Config()
    rows: list[NullAwareFeatureRow] = []
    for index, edge in enumerate(candidate_edges):
        edge_id, upstream_id, downstream_id, metadata = _edge_parts(edge, index)
        upstream = _observation_for(table, upstream_id)
        downstream = _observation_for(table, downstream_id)

        flow: dict[str, float] = {}
        null: dict[str, float] = {}
        source_fields: list[str] = []

        for name, aliases in _FLOW_ALIASES.items():
            value, source = _lookup_explicit(metadata, aliases, "flow")
            if value is not None:
                flow[name] = value
                if source is not None:
                    source_fields.append(source)

        for name, aliases in _NULL_ALIASES.items():
            value, source = _lookup_explicit(metadata, aliases, "null")
            if value is not None:
                null[name] = value
                if source is not None:
                    source_fields.append(source)

        if "head_drop_m" not in flow:
            value, source = _derive_head_drop(upstream, downstream, _HEAD_KEYS)
            if value is not None:
                flow["head_drop_m"] = value
                if source is not None:
                    source_fields.append(source)

        if "age_increment_years" not in flow:
            value, source = _derive_age_increment(upstream, downstream, _AGE_KEYS)
            if value is not None:
                flow["age_increment_years"] = value
                if source is not None:
                    source_fields.append(source)

        if "common_lithology" not in null:
            value, source = _derive_common_lithology(upstream, downstream)
            if value is not None:
                null["common_lithology"] = value
                if source is not None:
                    source_fields.append(source)

        if upstream and downstream:
            derived_null, derived_sources = _derive_null_features(
                upstream,
                downstream,
                feature_config,
            )
            for name, value in derived_null.items():
                null.setdefault(name, value)
            source_fields.extend(derived_sources)

        rows.append(
            NullAwareFeatureRow(
                edge_id=edge_id,
                flow_features=flow,
                null_features=null,
                source_fields=tuple(source_fields),
            )
        )
    return rows


def build_null_aware_feature_rows(
    candidate_edges: Iterable[Any],
    observations: Any = None,
    config: Optional[Config] = None,
) -> list[NullAwareFeatureRow]:
    """Explicitly named alias for :func:`build_feature_rows`."""

    return build_feature_rows(candidate_edges, observations, config=config)


def _coerce_feature_row(row: Any, index: int = 0) -> NullAwareFeatureRow:
    if isinstance(row, NullAwareFeatureRow):
        return row
    if isinstance(row, Mapping):
        return NullAwareFeatureRow(
            edge_id=str(row.get("edge_id", row.get("id", f"row-{index}"))),
            flow_features=_as_mapping(row.get("flow_features", {})),
            null_features=_as_mapping(row.get("null_features", {})),
            missing_channels=tuple(row.get("missing_channels", ())),
            source_fields=tuple(row.get("source_fields", ())),
        )
    raise TypeError("rows must be NullAwareFeatureRow instances or mappings")


def _coerce_rows(rows: Iterable[Any]) -> list[NullAwareFeatureRow]:
    if isinstance(rows, (str, bytes)):
        raise TypeError("rows must be an iterable of feature rows")
    return [_coerce_feature_row(row, index) for index, row in enumerate(rows)]


def _coerce_labels(labels: Any, rows: Sequence[NullAwareFeatureRow]) -> np.ndarray:
    if isinstance(labels, Mapping):
        values = []
        for row in rows:
            if row.edge_id not in labels:
                raise ValueError(f"missing label for edge {row.edge_id!r}")
            values.append(labels[row.edge_id])
    else:
        values = list(labels)
    if len(values) != len(rows):
        raise ValueError("labels and rows must have the same length")
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or not np.all(np.isfinite(result)) or not np.all(np.isin(result, (0.0, 1.0))):
        raise ValueError("labels must be finite binary values")
    return result


def _sigmoid(values: np.ndarray | float) -> np.ndarray | float:
    array = np.asarray(values, dtype=float)
    clipped = np.clip(array, -709.0, 709.0)
    result = np.where(clipped >= 0.0, 1.0 / (1.0 + np.exp(-clipped)), np.exp(clipped) / (1.0 + np.exp(clipped)))
    if np.ndim(values) == 0:
        return float(result)
    return result


def calibration_diagnostics(
    labels: Sequence[float] | np.ndarray,
    probabilities: Sequence[float] | np.ndarray,
    n_bins: int = 10,
) -> dict[str, Any]:
    """Return Brier score, log loss, and equal-width expected calibration error."""

    if n_bins < 2:
        raise ValueError("n_bins must be at least 2")
    y = np.asarray(labels, dtype=float)
    p = np.asarray(probabilities, dtype=float)
    if y.ndim != 1 or p.ndim != 1 or len(y) != len(p) or len(y) == 0:
        raise ValueError("labels and probabilities must be non-empty one-dimensional arrays of equal length")
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(p)):
        raise ValueError("labels and probabilities must be finite")
    if not np.all(np.isin(y, (0.0, 1.0))):
        raise ValueError("labels must be binary")
    if np.any((p < 0.0) | (p > 1.0)):
        raise ValueError("probabilities must be bounded in [0, 1]")

    clipped = np.clip(p, 1e-15, 1.0 - 1e-15)
    brier = float(np.mean((p - y) ** 2))
    log_loss = float(np.mean(-(y * np.log(clipped) + (1.0 - y) * np.log1p(-clipped))))

    # Include p == 1 in the last bin, while leaving p == 0 in the first.
    bins = np.minimum((p * n_bins).astype(int), n_bins - 1)
    reliability: list[dict[str, float | int]] = []
    ece = 0.0
    for bin_index in range(n_bins):
        mask = bins == bin_index
        if not np.any(mask):
            continue
        mean_probability = float(np.mean(p[mask]))
        mean_label = float(np.mean(y[mask]))
        count = int(np.sum(mask))
        weight = count / len(y)
        ece += weight * abs(mean_probability - mean_label)
        reliability.append(
            {
                "bin": bin_index,
                "count": count,
                "mean_probability": mean_probability,
                "mean_label": mean_label,
            }
        )

    return {
        "n": int(len(y)),
        "n_bins": int(n_bins),
        "brier": brier,
        "log_loss": log_loss,
        "ece": float(ece),
        "reliability": reliability,
    }


class CalibrationNotFittedError(RuntimeError):
    """Raised when a prediction is requested from an unfitted calibrator."""


class NullAwareLogisticCalibrator:
    """Deterministic L2-regularized logistic calibrator.

    The feature order, standardization values, coefficients, and optimizer
    settings are all retained in :meth:`to_dict`, making the fitted mapping
    inspectable and reproducible.  Missing numeric features are mean-imputed
    for the mathematical model, but only calibration-supported missingness
    patterns may receive a deployment probability.
    """

    def __init__(
        self,
        l2: float = 1.0,
        max_iter: int = 500,
        tolerance: float = 1e-9,
        fit_scope: str = "development",
        min_samples: int = 4,
        feature_names: Sequence[str] = _ALL_FEATURE_NAMES,
        calibration_provenance: Optional[Mapping[str, Any]] = None,
    ) -> None:
        if l2 < 0.0 or not math.isfinite(float(l2)):
            raise ValueError("l2 must be finite and non-negative")
        if max_iter < 1:
            raise ValueError("max_iter must be positive")
        if tolerance <= 0.0 or not math.isfinite(float(tolerance)):
            raise ValueError("tolerance must be finite and positive")
        if not str(fit_scope).strip():
            raise ValueError("fit_scope must be a non-empty string")
        if min_samples < 2:
            raise ValueError("min_samples must be at least 2")
        names = tuple(feature_names)
        if names != _ALL_FEATURE_NAMES:
            raise ValueError("feature_names must contain the canonical flow and null features in order")

        self.l2 = float(l2)
        self.max_iter = int(max_iter)
        self.tolerance = float(tolerance)
        self.fit_scope = str(fit_scope).strip()
        self.min_samples = int(min_samples)
        self.feature_names = names
        self.intercept_: Optional[float] = None
        self.coefficients_: Optional[np.ndarray] = None
        self.feature_means_: Optional[np.ndarray] = None
        self.feature_scales_: Optional[np.ndarray] = None
        self.converged_: Optional[bool] = None
        self.n_iter_: Optional[int] = None
        self.objective_: Optional[float] = None
        self.n_samples_: Optional[int] = None
        self.training_diagnostics_: Optional[dict[str, Any]] = None
        self.missingness_patterns_: Optional[tuple[tuple[str, ...], ...]] = None
        self.calibration_provenance = (
            dict(calibration_provenance)
            if isinstance(calibration_provenance, Mapping)
            else None
        )

    @property
    def fitted(self) -> bool:
        return self.coefficients_ is not None and self.intercept_ is not None

    @property
    def calibration_status(self) -> str:
        return CALIBRATION_FITTED if self.fitted else CALIBRATION_UNFITTED

    @property
    def deployment_eligible(self) -> bool:
        """Whether the artifact satisfies the public deployment contract.

        A fitted object is not automatically a deployable calibrator.  The
        public path requires a held-out/deployment scope, convergence, and
        explicit provenance asserting that the calibration cases were
        independent of development fitting.  This is metadata enforcement,
        not a substitute for an honest split; the provenance is therefore
        retained in the serialized artifact for audit.
        """

        provenance = self.calibration_provenance
        required = ("generator_id", "split_id", "dataset_hash")
        return bool(
            self.fitted
            and self.converged_ is True
            and self.fit_scope in _DEPLOYABLE_CALIBRATION_SCOPES
            and isinstance(provenance, Mapping)
            and provenance.get("independent") is True
            and all(str(provenance.get(key, "")).strip() for key in required)
        )

    @staticmethod
    def _raw_matrix(rows: Sequence[NullAwareFeatureRow]) -> np.ndarray:
        matrix = np.full((len(rows), len(_ALL_FEATURE_NAMES)), np.nan, dtype=float)
        for row_index, row in enumerate(rows):
            values = {**row.flow_features, **row.null_features}
            for column, name in enumerate(_ALL_FEATURE_NAMES):
                value = values.get(name)
                number = _as_finite_float(value)
                if number is not None:
                    matrix[row_index, column] = number
        return matrix

    def _matrix(self, rows: Sequence[NullAwareFeatureRow], fit: bool = False) -> np.ndarray:
        raw = self._raw_matrix(rows)
        if fit:
            if raw.shape[0] == 0:
                raise ValueError("at least one feature row is required")
            observed = np.isfinite(raw)
            counts = np.sum(observed, axis=0)
            safe_raw = np.where(observed, raw, 0.0)
            means = np.divide(
                np.sum(safe_raw, axis=0),
                counts,
                out=np.zeros(raw.shape[1], dtype=float),
                where=counts > 0,
            )
            deviations = np.where(observed, raw - means, 0.0)
            scales = np.sqrt(
                np.divide(
                    np.sum(deviations**2, axis=0),
                    counts,
                    out=np.ones(raw.shape[1], dtype=float),
                    where=counts > 0,
                )
            )
            scales = np.where(np.isfinite(scales) & (scales > 1e-12), scales, 1.0)
            self.feature_means_ = means
            self.feature_scales_ = scales
        if self.feature_means_ is None or self.feature_scales_ is None:
            raise CalibrationNotFittedError("calibrator has not been fitted")
        filled = np.where(np.isfinite(raw), raw, self.feature_means_)
        return (filled - self.feature_means_) / self.feature_scales_

    def fit(self, rows: Iterable[Any], labels: Any) -> "NullAwareLogisticCalibrator":
        feature_rows = _coerce_rows(rows)
        if not feature_rows:
            raise ValueError("at least one feature row is required")
        y = _coerce_labels(labels, feature_rows)
        if len(feature_rows) < self.min_samples:
            raise ValueError(
                f"at least {self.min_samples} feature rows are required for "
                "calibration"
            )
        if len(np.unique(y)) < 2:
            raise ValueError("both classes are required for calibration")
        x = self._matrix(feature_rows, fit=True)
        design = np.column_stack((np.ones(len(x)), x))
        initial = np.zeros(design.shape[1], dtype=float)

        def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
            logit = design @ parameters
            loss = float(np.mean(np.logaddexp(0.0, logit) - y * logit))
            loss += 0.5 * self.l2 * float(np.dot(parameters[1:], parameters[1:]))
            residual = np.asarray(_sigmoid(logit), dtype=float) - y
            gradient = design.T @ residual / len(y)
            gradient[1:] += self.l2 * parameters[1:]
            return loss, gradient

        result = minimize(
            lambda parameters: objective(parameters)[0],
            initial,
            jac=lambda parameters: objective(parameters)[1],
            method="L-BFGS-B",
            options={
                "maxiter": self.max_iter,
                "ftol": self.tolerance,
                "gtol": self.tolerance,
                "maxls": 50,
            },
        )
        self.intercept_ = float(result.x[0])
        self.coefficients_ = np.asarray(result.x[1:], dtype=float)
        self.converged_ = bool(result.success)
        self.n_iter_ = int(getattr(result, "nit", 0))
        self.objective_ = float(result.fun)
        self.n_samples_ = int(len(y))
        self.missingness_patterns_ = tuple(
            sorted(
                {
                    tuple(sorted(row.missing_channels))
                    for row in feature_rows
                }
            )
        )
        self.training_diagnostics_ = calibration_diagnostics(y, self.predict(feature_rows))
        return self

    def supports_missingness(self, row: NullAwareFeatureRow) -> bool:
        """Return whether ``row`` has a pattern represented during fitting."""

        if self.missingness_patterns_ is None:
            return False
        pattern = tuple(sorted(row.missing_channels))
        return pattern in self.missingness_patterns_

    def _require_fitted(self) -> None:
        if not self.fitted:
            raise CalibrationNotFittedError("calibrator has not been fitted")

    def _explain_rows(self, rows: Sequence[NullAwareFeatureRow]) -> list[dict[str, Any]]:
        self._require_fitted()
        x = self._matrix(rows)
        assert self.coefficients_ is not None
        assert self.intercept_ is not None
        flow_count = len(FLOW_FEATURE_NAMES)
        flow_logits = x[:, :flow_count] @ self.coefficients_[:flow_count]
        # null_logit is intentionally defined as evidence for the no-flow
        # explanation.  A positive null feature with a learned negative
        # coefficient consequently increases null_logit.
        null_logits = -(x[:, flow_count:] @ self.coefficients_[flow_count:])
        logits = self.intercept_ + flow_logits - null_logits
        probabilities = np.asarray(_sigmoid(logits), dtype=float)
        flow_probabilities = np.asarray(_sigmoid(flow_logits), dtype=float)
        # The calibrated no-flow probability is the complement of the same
        # binary flow probability.  The separate null logit remains an
        # explanatory diagnostic, not a second probability used for gating.
        null_scores = 1.0 - probabilities
        null_explanation_scores = np.asarray(_sigmoid(null_logits), dtype=float)
        explanations: list[dict[str, Any]] = []
        for index, row in enumerate(rows):
            explanations.append(
                {
                    "edge_id": row.edge_id,
                    "probability": float(np.clip(probabilities[index], 0.0, 1.0)),
                    "flow_support_probability": float(np.clip(flow_probabilities[index], 0.0, 1.0)),
                    "null_score": float(np.clip(null_scores[index], 0.0, 1.0)),
                    "null_explanation_score": float(
                        np.clip(null_explanation_scores[index], 0.0, 1.0)
                    ),
                    "logit": float(logits[index]),
                    "flow_support_logit": float(flow_logits[index]),
                    "null_logit": float(null_logits[index]),
                }
            )
        return explanations

    def predict(self, rows: Iterable[Any]) -> np.ndarray:
        feature_rows = _coerce_rows(rows)
        explanations = self._explain_rows(feature_rows)
        return np.asarray([item["probability"] for item in explanations], dtype=float)

    def predict_proba(self, rows: Iterable[Any]) -> np.ndarray:
        return self.predict(rows)

    def explain(self, rows: Iterable[Any]) -> list[dict[str, Any]]:
        return self._explain_rows(_coerce_rows(rows))

    def diagnostics(
        self,
        rows: Iterable[Any],
        labels: Any,
        n_bins: int = 10,
    ) -> dict[str, Any]:
        feature_rows = _coerce_rows(rows)
        y = _coerce_labels(labels, feature_rows)
        return calibration_diagnostics(y, self.predict(feature_rows), n_bins=n_bins)

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "calibration_status": self.calibration_status,
            "l2": self.l2,
            "max_iter": self.max_iter,
            "tolerance": self.tolerance,
            "calibration_fit_scope": self.fit_scope,
                "minimum_calibration_samples": self.min_samples,
                "feature_names": list(self.feature_names),
                "calibration_provenance": self.calibration_provenance,
                "deployment_eligible": self.deployment_eligible,
            }
        if not self.fitted:
            return result
        assert self.coefficients_ is not None
        assert self.intercept_ is not None
        assert self.feature_means_ is not None
        assert self.feature_scales_ is not None
        result.update(
            {
                "intercept": self.intercept_,
                "coefficients": self.coefficients_.tolist(),
                "feature_means": self.feature_means_.tolist(),
                "feature_scales": self.feature_scales_.tolist(),
                "converged": self.converged_,
                "n_iter": self.n_iter_,
                "objective": self.objective_,
                "n_samples": self.n_samples_,
                "training_diagnostics": self.training_diagnostics_,
                "calibration_missingness_patterns": [
                    list(pattern)
                    for pattern in (self.missingness_patterns_ or ())
                ],
                "optimizer": "scipy.optimize.minimize:L-BFGS-B",
            }
        )
        return result

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "NullAwareLogisticCalibrator":
        """Restore a fitted calibrator from an auditable ``to_dict`` payload."""

        if not isinstance(payload, Mapping):
            raise TypeError("calibrator payload must be a mapping")
        model = cls(
            l2=float(payload.get("l2", 1.0)),
            max_iter=int(payload.get("max_iter", 500)),
            tolerance=float(payload.get("tolerance", 1e-9)),
            fit_scope=str(payload.get("calibration_fit_scope", "development")),
            min_samples=int(payload.get("minimum_calibration_samples", 4)),
            feature_names=tuple(payload.get("feature_names", _ALL_FEATURE_NAMES)),
            calibration_provenance=payload.get("calibration_provenance"),
        )
        if payload.get("calibration_status") != CALIBRATION_FITTED:
            raise ValueError("payload does not contain a fitted calibrator")

        coefficients = np.asarray(payload.get("coefficients"), dtype=float)
        means = np.asarray(payload.get("feature_means"), dtype=float)
        scales = np.asarray(payload.get("feature_scales"), dtype=float)
        if (
            coefficients.shape != (len(_ALL_FEATURE_NAMES),)
            or means.shape != (len(_ALL_FEATURE_NAMES),)
            or scales.shape != (len(_ALL_FEATURE_NAMES),)
            or not np.all(np.isfinite(coefficients))
            or not np.all(np.isfinite(means))
            or not np.all(np.isfinite(scales))
            or np.any(scales <= 0.0)
        ):
            raise ValueError("payload contains invalid fitted coefficient arrays")
        intercept = _as_finite_float(payload.get("intercept"))
        if intercept is None:
            raise ValueError("payload contains an invalid fitted intercept")

        model.intercept_ = intercept
        model.coefficients_ = coefficients
        model.feature_means_ = means
        model.feature_scales_ = scales
        model.converged_ = bool(payload.get("converged", True))
        model.n_iter_ = int(payload.get("n_iter", 0))
        model.objective_ = _as_finite_float(payload.get("objective"))
        model.n_samples_ = int(payload.get("n_samples", 0))
        model.training_diagnostics_ = payload.get("training_diagnostics")
        raw_patterns = payload.get("calibration_missingness_patterns", [])
        if not isinstance(raw_patterns, Sequence) or isinstance(raw_patterns, (str, bytes)):
            raise ValueError("payload contains invalid calibration missingness patterns")
        model.missingness_patterns_ = tuple(
            tuple(sorted(str(item) for item in pattern))
            for pattern in raw_patterns
            if isinstance(pattern, Sequence) and not isinstance(pattern, (str, bytes))
        )
        return model


class NullAwareTopologyScorer:
    """Issue one auditable topology record per candidate edge."""

    def __init__(
        self,
        calibrator: NullAwareLogisticCalibrator,
        present_threshold: float = 0.75,
        absent_threshold: float = 0.25,
        require_flow_channel: bool = True,
        require_null_channel: bool = True,
    ) -> None:
        if not isinstance(calibrator, NullAwareLogisticCalibrator):
            raise TypeError("calibrator must be a NullAwareLogisticCalibrator")
        if not 0.0 <= absent_threshold < present_threshold <= 1.0:
            raise ValueError("thresholds must satisfy 0 <= absent < present <= 1")
        self.calibrator = calibrator
        self.present_threshold = float(present_threshold)
        self.absent_threshold = float(absent_threshold)
        self.require_flow_channel = bool(require_flow_channel)
        self.require_null_channel = bool(require_null_channel)

    def _has_required_evidence(self, row: NullAwareFeatureRow) -> bool:
        flow_present = bool(row.flow_features)
        null_present = bool(row.null_features)
        return (
            (not self.require_flow_channel or flow_present)
            and (not self.require_null_channel or null_present)
            and self.calibrator.supports_missingness(row)
        )

    def _unfitted_record(self, row: NullAwareFeatureRow) -> dict[str, Any]:
        return {
            "edge_id": row.edge_id,
            "decision": DECISION_ABSTAIN,
            "flow_probability": None,
            "flow_support_probability": None,
            "null_score": None,
            "null_explanation_score": None,
            "missing_channels": list(row.missing_channels),
            "calibration_status": CALIBRATION_UNFITTED,
            "calibration_fit_scope": None,
            "reason": "calibrator is unfitted; no topology decision issued",
            "flow_support_logit": None,
            "null_logit": None,
            "flow_features": dict(row.flow_features),
            "null_features": dict(row.null_features),
            "source_fields": list(row.source_fields),
        }

    def score_rows(self, rows: Iterable[Any]) -> list[dict[str, Any]]:
        feature_rows = _coerce_rows(rows)
        if not self.calibrator.fitted:
            return [self._unfitted_record(row) for row in feature_rows]

        explanations = self.calibrator.explain(feature_rows)
        records: list[dict[str, Any]] = []
        for row, explanation in zip(feature_rows, explanations):
            evidence_complete = self._has_required_evidence(row)
            if not evidence_complete:
                decision = DECISION_ABSTAIN
                status = CALIBRATION_FITTED_MISSING_EVIDENCE
                if not self.calibrator.supports_missingness(row):
                    reason = "missingness pattern was not represented during calibration"
                else:
                    reason = "required flow-support or null evidence channel is missing"
            else:
                probability = explanation["probability"]
                if probability >= self.present_threshold:
                    decision = DECISION_PRESENT
                    reason = "calibrated flow probability meets the PRESENT threshold"
                elif probability <= self.absent_threshold:
                    decision = DECISION_ABSENT
                    reason = "calibrated flow probability meets the ABSENT threshold"
                else:
                    decision = DECISION_ABSTAIN
                    reason = "calibrated flow probability is between decision thresholds"
                status = CALIBRATION_FITTED

            records.append(
                {
                    "edge_id": row.edge_id,
                    "decision": decision,
                    "flow_probability": (
                        float(np.clip(explanation["probability"], 0.0, 1.0))
                        if evidence_complete
                        else None
                    ),
                    "flow_support_probability": (
                        float(np.clip(explanation["flow_support_probability"], 0.0, 1.0))
                        if evidence_complete
                        else None
                    ),
                    "null_score": (
                        float(np.clip(explanation["null_score"], 0.0, 1.0))
                        if evidence_complete
                        else None
                    ),
                    "null_explanation_score": (
                        float(
                            np.clip(
                                explanation["null_explanation_score"], 0.0, 1.0
                            )
                        )
                        if evidence_complete
                        else None
                    ),
                    "missing_channels": list(row.missing_channels),
                    "calibration_status": status,
                    "calibration_fit_scope": self.calibrator.fit_scope,
                    "reason": reason,
                    "flow_support_logit": explanation["flow_support_logit"],
                    "null_logit": explanation["null_logit"],
                    "flow_features": dict(row.flow_features),
                    "null_features": dict(row.null_features),
                    "source_fields": list(row.source_fields),
                }
            )
        return records

    def score_edges(
        self,
        candidate_edges: Iterable[Any],
        observations: Any = None,
        config: Optional[Config] = None,
    ) -> list[dict[str, Any]]:
        return self.score_rows(build_feature_rows(candidate_edges, observations, config=config))

    def score(
        self,
        candidate_edges: Iterable[Any],
        observations: Any = None,
        config: Optional[Config] = None,
    ) -> list[dict[str, Any]]:
        """Alias for :meth:`score_edges` for callers treating this as a scorer."""

        return self.score_edges(candidate_edges, observations, config=config)


def score_null_aware_topology(
    candidate_edges: Iterable[Any],
    observations: Any,
    calibrator: NullAwareLogisticCalibrator,
    **scorer_options: Any,
) -> list[dict[str, Any]]:
    """Convenience function returning one auditable record per candidate edge."""

    config = scorer_options.pop("config", None)
    scorer = NullAwareTopologyScorer(calibrator, **scorer_options)
    return scorer.score_edges(candidate_edges, observations, config=config)


__all__ = [
    "CALIBRATION_FITTED",
    "CALIBRATION_FITTED_MISSING_EVIDENCE",
    "CALIBRATION_UNFITTED",
    "CalibrationNotFittedError",
    "DECISION_ABSENT",
    "DECISION_ABSTAIN",
    "DECISION_PRESENT",
    "FLOW_FEATURE_NAMES",
    "NULL_FEATURE_NAMES",
    "FeatureRow",
    "NullAwareFeatureRow",
    "NullAwareLogisticCalibrator",
    "NullAwareTopologyScorer",
    "build_feature_rows",
    "build_null_aware_feature_rows",
    "calibration_diagnostics",
    "score_null_aware_topology",
]
