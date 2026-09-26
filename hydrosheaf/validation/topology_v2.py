"""Evidence-gated topology-v2 candidate generation and calibration.

This module is intentionally separate from ``hydrosheaf.graph.build`` and
``hydrosheaf.inference.topology_posterior``.  The older path is useful for
legacy inference and diagnostics, but its candidate graph can apply downhill
and nearest-neighbour gates before the evidence model sees an edge.  That is
not a defensible candidate universe for a held-out topology benchmark.

The v2 contract is:

* candidate generation is truth-blind and records every physical rejection;
* direction, distance, screen, aquifer, capacity, and geology are soft
  evidence unless a caller explicitly opts into a hard admissibility rule;
* the candidate prior is an uncalibrated screening quantity, never a final
  posterior probability;
* calibration is fit on whole-case development splits and is deployable only
  with explicit independent provenance;
* missingness is represented and unsupported missingness produces ABSTAIN;
* a selected inference threshold must be strictly inside (0, 1).

The implementation uses NumPy/SciPy only so that the scientific contract is
available in the base package.  Scikit-learn remains useful for optional
cross-checks, but is not required by this module.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
import hashlib
import json
import math
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, ndtr

from .baselines import assert_truth_blind_observations


TOPOLOGY_V2_VERSION = "topology_v2.1.0"
TOPOLOGY_V2_ALGORITHM = "truth_blind_all_pairs_soft_admissibility_v2"

DECISION_PRESENT = "PRESENT"
DECISION_ABSENT = "ABSENT"
DECISION_ABSTAIN = "ABSTAIN"

GEOMETRY_FEATURES = (
    "distance_m",
    "log_distance_m",
    "weak_geometry_prior",
)
ELEVATION_FEATURES = (
    "elevation_delta",
    "elevation_direction_probability",
)
HEAD_FEATURES = (
    "head_delta_m",
    "head_sigma_delta_m",
    "head_z_score",
)
HEAD_GRADIENT_FEATURES = (
    "direction_probability",
    "gradient_m_per_km",
)
PATH_FEATURES = (
    "path_reachability",
    "path_relative_support",
    "path_log_support",
    "path_steps",
    "path_transition_efficiency",
    "path_activity_support",
    "path_endpoint_sink_support",
)
SCREEN_FEATURES = (
    "vertical_separation_m",
    "screen_overlap_fraction",
    "screen_compatibility",
    "aquifer_compatibility",
)
CAPACITY_FEATURES = ("hydraulic_capacity_proxy",)
GEOLOGY_FEATURES = ("geology_similarity",)
SOURCE_SINK_FEATURES = (
    "source_pumping_strength",
    "target_pumping_sink_strength",
    "source_recharge_strength",
    "target_river_sink_strength",
    "target_boundary_sink_strength",
)
CHEMISTRY_TRACER_FEATURES = (
    "chemistry_similarity",
    "isotope_similarity",
    "tracer_direction_support",
)

MODEL_FEATURE_SETS: dict[str, tuple[str, ...]] = {
    "A_geometry": GEOMETRY_FEATURES,
    "M4_A_sparse": GEOMETRY_FEATURES + ELEVATION_FEATURES,
    "B_geometry_head": GEOMETRY_FEATURES + HEAD_FEATURES,
    "C_geometry_head_gradient": (
        GEOMETRY_FEATURES + HEAD_FEATURES + HEAD_GRADIENT_FEATURES
    ),
    "M4_B_archive_informed": (
        GEOMETRY_FEATURES
        + HEAD_FEATURES
        + HEAD_GRADIENT_FEATURES
        + SOURCE_SINK_FEATURES
    ),
    "M4_C_path_aware": (
        GEOMETRY_FEATURES
        + HEAD_FEATURES
        + HEAD_GRADIENT_FEATURES
        + SOURCE_SINK_FEATURES
        + PATH_FEATURES
    ),
    "D_plus_screen": (
        GEOMETRY_FEATURES
        + HEAD_FEATURES
        + HEAD_GRADIENT_FEATURES
        + SCREEN_FEATURES
    ),
    "E_plus_hydraulic_capacity": (
        GEOMETRY_FEATURES
        + HEAD_FEATURES
        + HEAD_GRADIENT_FEATURES
        + SCREEN_FEATURES
        + CAPACITY_FEATURES
    ),
    "F_plus_geology": (
        GEOMETRY_FEATURES
        + HEAD_FEATURES
        + HEAD_GRADIENT_FEATURES
        + SCREEN_FEATURES
        + CAPACITY_FEATURES
        + GEOLOGY_FEATURES
    ),
    "G_chemistry_tracer_exploratory": (
        GEOMETRY_FEATURES
        + HEAD_FEATURES
        + HEAD_GRADIENT_FEATURES
        + SCREEN_FEATURES
        + CAPACITY_FEATURES
        + GEOLOGY_FEATURES
        + CHEMISTRY_TRACER_FEATURES
    ),
}

_HEAD_KEYS = (
    "head_evidence_value",
    "head_meas",
    "hydraulic_head",
    "head",
    "water_level",
)
_HEAD_SIGMA_KEYS = (
    "head_sigma_m",
    "head_measurement_sigma_m",
    "head_evidence_sigma_m",
)
_ELEVATION_KEYS = ("elevation_m", "elevation", "z_m", "z")
_GEOMETRY_KEYS = ("x_m", "y_m", "lat", "lon")
_SCREEN_TOP_KEYS = ("screen_top_m", "screen_depth", "screen_top")
_SCREEN_BOTTOM_KEYS = ("screen_bottom_m", "well_depth", "screen_bottom")
_AQUIFER_KEYS = ("aquifer_unit", "aquifer_layer", "lithology")
_GEOLOGY_KEYS = (
    "geology_symbol",
    "geology_code",
    "geology_unit",
    "geologic_unit",
    "mapped_geology",
    "aquifer_unit",
)
_CAPACITY_KEYS = (
    "transmissivity_m2_per_day",
    "transmissivity",
    "hydraulic_conductivity_m_per_s",
    "hydraulic_conductivity",
    "conductivity",
)
_THICKNESS_KEYS = ("aquifer_thickness_m", "thickness_m", "screen_thickness_m")
_CHEMISTRY_KEYS = (
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
    "Fe",
    "PO4",
    "SiO2",
    "pH",
)
_ISOTOPE_KEYS = (
    "18O",
    "d18O",
    "delta18O",
    "delta_18O",
    "2H",
    "d2H",
    "delta2H",
    "delta_2H",
)
_AGE_KEYS = ("mean_age_years", "age_years", "age")
_TOPOLOGY_V2_TRUTH_NAMES = {
    "actual",
    "edge_label",
    "heldout_truth",
    "is_connected",
    "is_true_flow",
    "label",
    "observed_present",
    "reference",
    "reference_edges",
    "target",
    "true_edges",
    "truth_edges",
}


def _finite(value: object) -> float | None:
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _normalise_key(value: object) -> str:
    return str(value).strip().lower().replace(" ", "_")


def _jsonable(value: object) -> object:
    return json.loads(json.dumps(value, sort_keys=True, default=str))


def _hash_payload(value: object) -> str:
    payload = json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"), default=str
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, float(value)))


def _lookup(row: Mapping[str, object], keys: Sequence[str]) -> tuple[object, str | None]:
    normalised = {_normalise_key(key): (value, str(key)) for key, value in row.items()}
    for key in keys:
        found = normalised.get(_normalise_key(key))
        if found is not None:
            return found
    return None, None


def _lookup_number(
    row: Mapping[str, object], keys: Sequence[str]
) -> tuple[float | None, str | None]:
    value, source = _lookup(row, keys)
    return _finite(value), source


def _lookup_text(
    row: Mapping[str, object], keys: Sequence[str]
) -> tuple[str | None, str | None]:
    value, source = _lookup(row, keys)
    return _text(value), source


def _normalise_observations(
    observations: Iterable[Mapping[str, object]] | Mapping[str, Mapping[str, object]],
) -> tuple[dict[str, object], ...]:
    if isinstance(observations, Mapping):
        rows: list[dict[str, object]] = []
        for key, value in observations.items():
            row = dict(value)
            row.setdefault("site_id", key)
            rows.append(row)
    else:
        rows = [dict(row) for row in observations]
    assert_truth_blind_observations(rows)

    def visit(value: object, path: tuple[str, ...]) -> None:
        if isinstance(value, Mapping):
            for key, nested in value.items():
                key_text = str(key)
                lowered = _normalise_key(key_text)
                if (
                    lowered in _TOPOLOGY_V2_TRUTH_NAMES
                    or lowered.endswith("_truth")
                    or lowered.endswith("_label")
                    or "reference" in lowered
                ):
                    dotted = ".".join((*path, key_text))
                    raise ValueError(f"Truth/reference field is forbidden: {dotted}")
                visit(nested, (*path, key_text))
        elif isinstance(value, (list, tuple)):
            for index, nested in enumerate(value):
                visit(nested, (*path, str(index)))

    visit(rows, ())

    normalised: list[dict[str, object]] = []
    for row in rows:
        raw_id = row.get("site_id", row.get("sample_id", row.get("node_id")))
        site_id = _text(raw_id)
        if site_id is None:
            raise ValueError("Every topology-v2 observation needs site_id or sample_id.")
        row["site_id"] = site_id
        normalised.append(row)
    normalised.sort(key=lambda row: str(row["site_id"]))
    ids = [str(row["site_id"]) for row in normalised]
    if len(ids) != len(set(ids)):
        raise ValueError("Topology-v2 observations contain duplicate site IDs.")
    return tuple(normalised)


def _coordinate_system(rows: Sequence[Mapping[str, object]]) -> str:
    cartesian = all(
        _finite(row.get("x_m")) is not None and _finite(row.get("y_m")) is not None
        for row in rows
    )
    geographic = all(
        _finite(row.get("lat")) is not None and _finite(row.get("lon")) is not None
        for row in rows
    )
    if cartesian:
        return "x_m_y_m"
    if geographic:
        return "lon_lat_local_tangent_m"
    raise ValueError(
        "Topology-v2 requires complete x_m/y_m or complete lat/lon coordinates "
        "for every observation; mixed or missing coordinate systems are not silently combined."
    )


def _project_coordinates(
    rows: Sequence[Mapping[str, object]], coordinate_system: str
) -> dict[str, tuple[float, float]]:
    if coordinate_system == "x_m_y_m":
        return {
            str(row["site_id"]): (float(row["x_m"]), float(row["y_m"]))
            for row in rows
        }
    lat0 = float(np.mean([float(row["lat"]) for row in rows]))
    lon0 = float(np.mean([float(row["lon"]) for row in rows]))
    cos_lat = max(0.1, math.cos(math.radians(lat0)))
    metres_lat = 110_540.0
    metres_lon = 111_320.0 * cos_lat
    return {
        str(row["site_id"]): (
            (float(row["lon"]) - lon0) * metres_lon,
            (float(row["lat"]) - lat0) * metres_lat,
        )
        for row in rows
    }


def _head_info(row: Mapping[str, object], default_sigma_m: float) -> tuple[float | None, float, str | None]:
    value, source = _lookup_number(row, _HEAD_KEYS)
    sigma, sigma_source = _lookup_number(row, _HEAD_SIGMA_KEYS)
    if sigma is None or sigma <= 0.0:
        sigma = float(default_sigma_m)
        sigma_source = None
    if value is None:
        return None, float(sigma), None
    return float(value), float(sigma), source or sigma_source


def _signed_log_context(
    row: Mapping[str, object], keys: Sequence[str], *, positive: bool
) -> tuple[float | None, str | None]:
    """Return a graded non-negative source/sink strength from a signed field."""

    value, source = _lookup_number(row, keys)
    if value is None:
        return None, source
    magnitude = value if positive else -value
    return float(math.log1p(max(0.0, magnitude))), source


def _screen_interval(row: Mapping[str, object]) -> tuple[float, float, str | None, str | None] | None:
    top, top_source = _lookup_number(row, _SCREEN_TOP_KEYS)
    bottom, bottom_source = _lookup_number(row, _SCREEN_BOTTOM_KEYS)
    if top is None or bottom is None:
        return None
    if bottom < top:
        top, bottom = bottom, top
        top_source, bottom_source = bottom_source, top_source
    if math.isclose(bottom, top):
        return None
    return float(top), float(bottom), top_source, bottom_source


def _interval_features(
    source: Mapping[str, object], target: Mapping[str, object]
) -> tuple[float | None, float | None, float | None, tuple[str, ...]]:
    source_interval = _screen_interval(source)
    target_interval = _screen_interval(target)
    if source_interval is None or target_interval is None:
        return None, None, None, ()
    source_top, source_bottom, source_top_key, source_bottom_key = source_interval
    target_top, target_bottom, target_top_key, target_bottom_key = target_interval
    overlap = max(0.0, min(source_bottom, target_bottom) - max(source_top, target_top))
    source_length = source_bottom - source_top
    target_length = target_bottom - target_top
    overlap_fraction = overlap / max(min(source_length, target_length), 1e-12)
    source_mid = 0.5 * (source_top + source_bottom)
    target_mid = 0.5 * (target_top + target_bottom)
    separation = abs(source_mid - target_mid)
    compatibility = 1.0 if overlap > 0.0 else 0.0
    return (
        float(separation),
        float(overlap_fraction),
        float(compatibility),
        tuple(
            key
            for key in (source_top_key, source_bottom_key, target_top_key, target_bottom_key)
            if key is not None
        ),
    )


def _compatibility(
    source: Mapping[str, object], target: Mapping[str, object], keys: Sequence[str]
) -> tuple[float | None, tuple[str, ...]]:
    source_value, source_key = _lookup_text(source, keys)
    target_value, target_key = _lookup_text(target, keys)
    if source_value is None or target_value is None:
        return None, tuple(key for key in (source_key, target_key) if key is not None)
    return (
        1.0 if _normalise_key(source_value) == _normalise_key(target_value) else 0.0,
        tuple(key for key in (source_key, target_key) if key is not None),
    )


def _capacity_value(row: Mapping[str, object]) -> tuple[float | None, tuple[str, ...]]:
    direct, direct_key = _lookup_number(row, _CAPACITY_KEYS)
    thickness, thickness_key = _lookup_number(row, _THICKNESS_KEYS)
    if direct is None or direct <= 0.0:
        return None, tuple(key for key in (direct_key, thickness_key) if key is not None)
    raw = direct * thickness if thickness is not None and thickness > 0.0 else direct
    if raw <= 0.0 or not math.isfinite(raw):
        return None, tuple(key for key in (direct_key, thickness_key) if key is not None)
    return float(math.log10(raw)), tuple(
        key for key in (direct_key, thickness_key) if key is not None
    )


def _similarity(
    source: Mapping[str, object],
    target: Mapping[str, object],
    keys: Sequence[str],
) -> tuple[float | None, tuple[str, ...]]:
    differences: list[float] = []
    used: list[str] = []
    source_norm = {_normalise_key(key): (value, str(key)) for key, value in source.items()}
    target_norm = {_normalise_key(key): (value, str(key)) for key, value in target.items()}
    for key in keys:
        left = source_norm.get(_normalise_key(key))
        right = target_norm.get(_normalise_key(key))
        if left is None or right is None:
            continue
        left_value = _finite(left[0])
        right_value = _finite(right[0])
        if left_value is None or right_value is None:
            continue
        scale = abs(left_value) + abs(right_value) + 1e-12
        differences.append(abs(left_value - right_value) / scale)
        used.extend((left[1], right[1]))
    if not differences:
        return None, tuple(dict.fromkeys(used))
    return float(1.0 / (1.0 + float(np.mean(differences)))), tuple(dict.fromkeys(used))


@dataclass(frozen=True)
class TopologyV2Config:
    """Candidate-generation and evidence defaults for topology-v2."""

    max_distance_m: float | None = None
    max_neighbors: int | None = None
    hard_direction: bool = False
    hard_screen_compatibility: bool = False
    hard_aquifer_compatibility: bool = False
    default_head_sigma_m: float = 0.10
    head_tie_tolerance_m: float = 0.05
    geometry_scale_m: float | None = None

    def __post_init__(self) -> None:
        if self.max_distance_m is not None and (
            not math.isfinite(float(self.max_distance_m)) or float(self.max_distance_m) <= 0.0
        ):
            raise ValueError("max_distance_m must be finite and positive when supplied.")
        if self.max_neighbors is not None and int(self.max_neighbors) < 1:
            raise ValueError("max_neighbors must be at least one when supplied.")
        if not math.isfinite(float(self.default_head_sigma_m)) or self.default_head_sigma_m <= 0.0:
            raise ValueError("default_head_sigma_m must be finite and positive.")
        if not math.isfinite(float(self.head_tie_tolerance_m)) or self.head_tie_tolerance_m < 0.0:
            raise ValueError("head_tie_tolerance_m must be finite and non-negative.")
        if self.geometry_scale_m is not None and (
            not math.isfinite(float(self.geometry_scale_m)) or float(self.geometry_scale_m) <= 0.0
        ):
            raise ValueError("geometry_scale_m must be finite and positive when supplied.")

    def to_dict(self) -> dict[str, object]:
        return {
            "max_distance_m": self.max_distance_m,
            "max_neighbors": self.max_neighbors,
            "hard_direction": self.hard_direction,
            "hard_screen_compatibility": self.hard_screen_compatibility,
            "hard_aquifer_compatibility": self.hard_aquifer_compatibility,
            "default_head_sigma_m": self.default_head_sigma_m,
            "head_tie_tolerance_m": self.head_tie_tolerance_m,
            "geometry_scale_m": self.geometry_scale_m,
        }


@dataclass(frozen=True)
class TopologyV2CandidateEdge:
    u: str
    v: str
    attrs: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        u = str(self.u).strip()
        v = str(self.v).strip()
        if not u or not v or u == v:
            raise ValueError("Topology-v2 candidate edges need distinct non-empty endpoints.")
        object.__setattr__(self, "u", u)
        object.__setattr__(self, "v", v)
        object.__setattr__(self, "attrs", dict(self.attrs))

    @property
    def edge_id(self) -> str:
        return f"{self.u}->{self.v}"

    def to_dict(self) -> dict[str, object]:
        return {"edge_id": self.edge_id, "u": self.u, "v": self.v, "attrs": _jsonable(self.attrs)}


@dataclass(frozen=True)
class TopologyV2CandidateUniverse:
    edges: tuple[TopologyV2CandidateEdge, ...]
    nodes: tuple[str, ...]
    algorithm: str
    version: str
    coordinate_system: str
    input_channels: tuple[str, ...]
    parameters: Mapping[str, object]
    input_hash: str
    candidate_hash: str
    rejection_counts: Mapping[str, int] = field(default_factory=dict)
    truth_blind: bool = True

    def __post_init__(self) -> None:
        ordered = tuple(sorted(self.edges, key=lambda edge: edge.edge_id))
        if ordered != tuple(self.edges):
            raise ValueError("Topology-v2 candidate edges must be deterministically sorted.")
        if not self.truth_blind:
            raise ValueError("Topology-v2 candidate universes must be truth-blind.")

    @property
    def edge_ids(self) -> tuple[str, ...]:
        return tuple(edge.edge_id for edge in self.edges)

    def candidate_graph_recall(self, reference_edges: Iterable[object]) -> float | None:
        reference = _normalise_edge_set(reference_edges)
        if not reference:
            return None
        retained = set(self.edge_ids)
        return float(len(reference & retained) / len(reference))

    def to_dict(self) -> dict[str, object]:
        return {
            "algorithm": self.algorithm,
            "version": self.version,
            "coordinate_system": self.coordinate_system,
            "input_channels": list(self.input_channels),
            "parameters": _jsonable(self.parameters),
            "input_hash": self.input_hash,
            "candidate_hash": self.candidate_hash,
            "n_nodes": len(self.nodes),
            "n_candidate_edges": len(self.edges),
            "rejection_counts": dict(self.rejection_counts),
            "truth_blind": self.truth_blind,
            "candidate_edges": [edge.to_dict() for edge in self.edges],
        }


def _pair_attributes(
    source: Mapping[str, object],
    target: Mapping[str, object],
    coordinates: Mapping[str, tuple[float, float]],
    config: TopologyV2Config,
    geometry_scale_m: float,
) -> dict[str, object]:
    u = str(source["site_id"])
    v = str(target["site_id"])
    x1, y1 = coordinates[u]
    x2, y2 = coordinates[v]
    distance_m = math.hypot(x1 - x2, y1 - y2)
    source_head, source_sigma, source_head_key = _head_info(source, config.default_head_sigma_m)
    target_head, target_sigma, target_head_key = _head_info(target, config.default_head_sigma_m)
    if source_head is not None and target_head is not None:
        delta = source_head - target_head
        sigma_delta = math.sqrt(source_sigma**2 + target_sigma**2)
        z_score = delta / sigma_delta
        direction_probability = float(ndtr(z_score))
        gradient = delta / (distance_m / 1000.0) if distance_m > 0.0 else None
    else:
        delta = None
        sigma_delta = None
        z_score = None
        direction_probability = None
        gradient = None
    source_elevation, source_elevation_key = _lookup_number(source, _ELEVATION_KEYS)
    target_elevation, target_elevation_key = _lookup_number(target, _ELEVATION_KEYS)
    if source_elevation is not None and target_elevation is not None:
        elevation_delta = source_elevation - target_elevation
        elevation_scale = max(abs(elevation_delta) * 0.25, 1e-6)
        elevation_direction_probability = float(ndtr(elevation_delta / elevation_scale))
    else:
        elevation_delta = None
        elevation_direction_probability = None
    screen_separation, overlap_fraction, screen_compatibility, screen_keys = _interval_features(
        source, target
    )
    aquifer_compatibility, aquifer_keys = _compatibility(source, target, _AQUIFER_KEYS)
    geometry_prior = 1.0 / (1.0 + distance_m / max(geometry_scale_m, 1e-12))
    direction_for_prior = 0.5 if direction_probability is None else direction_probability
    compatibility_values = [
        value for value in (screen_compatibility, aquifer_compatibility) if value is not None
    ]
    compatibility_prior = (
        float(np.mean(compatibility_values)) if compatibility_values else 0.5
    )
    candidate_prior = _clamp(
        0.50 * direction_for_prior + 0.30 * geometry_prior + 0.20 * compatibility_prior
    )
    return {
        "candidate_generator": TOPOLOGY_V2_ALGORITHM,
        "coordinate_system": "projected_local_m",
        "distance_m": float(distance_m),
        "weak_geometry_prior": float(geometry_prior),
        "candidate_prior_probability": float(candidate_prior),
        "head_source_field": source_head_key,
        "head_target_field": target_head_key,
        "head_delta_m": delta,
        "head_sigma_delta_m": sigma_delta,
        "head_z_score": z_score,
        "direction_probability": direction_probability,
        "gradient_m_per_km": gradient,
        "elevation_source_field": source_elevation_key,
        "elevation_target_field": target_elevation_key,
        "elevation_delta": elevation_delta,
        "elevation_direction_probability": elevation_direction_probability,
        "screen_overlap_fraction": overlap_fraction,
        "screen_compatibility": screen_compatibility,
        "vertical_separation_m": screen_separation,
        "aquifer_compatibility": aquifer_compatibility,
        "screen_source_fields": list(screen_keys),
        "aquifer_source_fields": list(aquifer_keys),
        "hard_direction_gate": config.hard_direction,
        "hard_screen_gate": config.hard_screen_compatibility,
        "hard_aquifer_gate": config.hard_aquifer_compatibility,
    }


def generate_topology_v2_candidate_universe(
    observations: Iterable[Mapping[str, object]] | Mapping[str, Mapping[str, object]],
    *,
    config: TopologyV2Config | None = None,
) -> TopologyV2CandidateUniverse:
    """Generate a deterministic directed candidate graph without truth.

    With the default configuration every directed pair is retained.  A caller
    may impose distance, neighbour, direction, screen, or aquifer gates, but
    each such rejection is counted and the resulting candidate recall must be
    reported against a reference graph later.  In particular, a downhill gate
    is never implicit.
    """

    resolved = config or TopologyV2Config()
    rows = _normalise_observations(observations)
    if len(rows) < 2:
        raise ValueError("Topology-v2 candidate generation requires at least two observations.")
    coordinate_system = _coordinate_system(rows)
    coordinates = _project_coordinates(rows, coordinate_system)
    by_id = {str(row["site_id"]): row for row in rows}
    pair_distances = [
        math.hypot(
            coordinates[left][0] - coordinates[right][0],
            coordinates[left][1] - coordinates[right][1],
        )
        for index, left in enumerate(sorted(by_id))
        for right in sorted(by_id)[index + 1 :]
    ]
    positive_distances = sorted(distance for distance in pair_distances if distance > 0.0)
    if resolved.geometry_scale_m is not None:
        geometry_scale_m = float(resolved.geometry_scale_m)
    elif positive_distances:
        middle = len(positive_distances) // 2
        geometry_scale_m = (
            positive_distances[middle]
            if len(positive_distances) % 2
            else 0.5 * (positive_distances[middle - 1] + positive_distances[middle])
        )
    else:
        geometry_scale_m = 1.0

    candidate_targets: dict[str, list[str]] = {site_id: [] for site_id in by_id}
    for source_id in sorted(by_id):
        ranked = sorted(
            (
                math.hypot(
                    coordinates[source_id][0] - coordinates[target_id][0],
                    coordinates[source_id][1] - coordinates[target_id][1],
                ),
                target_id,
            )
            for target_id in by_id
            if target_id != source_id
        )
        if resolved.max_neighbors is not None:
            ranked = ranked[: int(resolved.max_neighbors)]
        candidate_targets[source_id] = [target_id for _, target_id in ranked]

    rejection_counts: dict[str, int] = {}
    edges: list[TopologyV2CandidateEdge] = []
    for source_id in sorted(by_id):
        source = by_id[source_id]
        for target_id in candidate_targets[source_id]:
            target = by_id[target_id]
            attrs = _pair_attributes(source, target, coordinates, resolved, geometry_scale_m)
            reasons: list[str] = []
            distance_m = float(attrs["distance_m"])
            # Distinct wells/cells may share a projected coordinate (for
            # example, vertically stacked cells).  They are not self-loops
            # and must remain in the default candidate universe; distance is
            # then a neutral/tied geometry feature rather than a rejection.
            if resolved.max_distance_m is not None and distance_m > float(resolved.max_distance_m):
                reasons.append("max_distance")
            direction = attrs.get("direction_probability")
            if resolved.hard_direction and direction is not None:
                delta = _finite(attrs.get("head_delta_m"))
                if delta is not None and delta <= float(resolved.head_tie_tolerance_m):
                    reasons.append("hard_direction")
            screen = attrs.get("screen_compatibility")
            if resolved.hard_screen_compatibility and screen is not None and float(screen) <= 0.0:
                reasons.append("hard_screen_compatibility")
            aquifer = attrs.get("aquifer_compatibility")
            if resolved.hard_aquifer_compatibility and aquifer is not None and float(aquifer) <= 0.0:
                reasons.append("hard_aquifer_compatibility")
            if reasons:
                for reason in reasons:
                    rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
                continue
            attrs["physical_admissibility"] = "soft_or_validated"
            edges.append(TopologyV2CandidateEdge(source_id, target_id, attrs))
    edges.sort(key=lambda edge: edge.edge_id)
    input_channels = tuple(
        key
        for key in (
            "site_id",
            "sample_id",
            "x_m",
            "y_m",
            "lat",
            "lon",
            "head_meas",
            "hydraulic_head",
            "head_sigma_m",
            "elevation_m",
            "elevation",
            "well_rate",
            "river_leakage",
            "recharge",
            "head_boundary_flux",
            "screen_depth",
            "well_depth",
            "aquifer_unit",
            "aquifer_layer",
        )
        if any(key in row for row in rows)
    )
    parameters = {
        **resolved.to_dict(),
        "geometry_scale_m_resolved": geometry_scale_m,
        "default_candidate_policy": "all_directed_pairs_when_no_explicit_cap",
        "direction_is_soft_by_default": not resolved.hard_direction,
        "candidate_prior_is_not_calibrated_posterior": True,
    }
    return TopologyV2CandidateUniverse(
        edges=tuple(edges),
        nodes=tuple(sorted(by_id)),
        algorithm=TOPOLOGY_V2_ALGORITHM,
        version=TOPOLOGY_V2_VERSION,
        coordinate_system=coordinate_system,
        input_channels=input_channels,
        parameters=parameters,
        input_hash=_hash_payload(rows),
        candidate_hash=_hash_payload([edge.to_dict() for edge in edges]),
        rejection_counts=dict(sorted(rejection_counts.items())),
    )


@dataclass(frozen=True)
class TopologyV2FeatureRow:
    edge_id: str
    case_id: str
    u: str
    v: str
    features: Mapping[str, float | None]
    missing_features: tuple[str, ...]
    source_fields: tuple[str, ...]
    candidate_prior_probability: float | None = None

    def __post_init__(self) -> None:
        if self.edge_id != f"{self.u}->{self.v}":
            raise ValueError("Topology-v2 feature edge_id must match u->v.")
        if tuple(sorted(self.missing_features)) != tuple(self.missing_features):
            raise ValueError("Topology-v2 missing feature names must be sorted.")

    def to_dict(self) -> dict[str, object]:
        return {
            "edge_id": self.edge_id,
            "case_id": self.case_id,
            "u": self.u,
            "v": self.v,
            "features": _jsonable(self.features),
            "missing_features": list(self.missing_features),
            "source_fields": list(self.source_fields),
            "candidate_prior_probability": self.candidate_prior_probability,
        }


def _feature_values(
    edge: TopologyV2CandidateEdge,
    source: Mapping[str, object],
    target: Mapping[str, object],
    config: TopologyV2Config,
) -> tuple[dict[str, float | None], tuple[str, ...]]:
    attrs = dict(edge.attrs)
    source_fields: list[str] = []
    for key in (
        "head_source_field",
        "head_target_field",
        "elevation_source_field",
        "elevation_target_field",
        "screen_source_fields",
        "aquifer_source_fields",
    ):
        value = attrs.get(key)
        if isinstance(value, str) and value:
            source_fields.append(value)
        elif isinstance(value, (list, tuple)):
            source_fields.extend(str(item) for item in value if str(item))

    source_capacity, source_capacity_fields = _capacity_value(source)
    target_capacity, target_capacity_fields = _capacity_value(target)
    capacity_values = [
        value for value in (source_capacity, target_capacity) if value is not None
    ]
    capacity = float(np.mean(capacity_values)) if capacity_values else None
    source_fields.extend(source_capacity_fields)
    source_fields.extend(target_capacity_fields)
    geology, geology_fields = _compatibility(source, target, _GEOLOGY_KEYS)
    source_fields.extend(geology_fields)
    chemistry, chemistry_fields = _similarity(source, target, _CHEMISTRY_KEYS)
    isotopes, isotope_fields = _similarity(source, target, _ISOTOPE_KEYS)
    source_fields.extend(chemistry_fields)
    source_fields.extend(isotope_fields)

    source_age, source_age_key = _lookup_number(source, _AGE_KEYS)
    target_age, target_age_key = _lookup_number(target, _AGE_KEYS)
    if source_age is not None and target_age is not None:
        age_delta = target_age - source_age
        tracer_direction_support = float(ndtr(age_delta / max(1.0, abs(age_delta) * 0.25)))
    else:
        tracer_direction_support = None
    source_fields.extend(key for key in (source_age_key, target_age_key) if key is not None)

    source_pumping, source_pumping_key = _signed_log_context(
        source, ("well_rate", "pumping_rate", "well_flux"), positive=True
    )
    target_pumping_sink, target_pumping_key = _signed_log_context(
        target, ("well_rate", "pumping_rate", "well_flux"), positive=False
    )
    source_recharge, source_recharge_key = _signed_log_context(
        source, ("recharge", "recharge_rate", "recharge_flux"), positive=True
    )
    target_river_sink, target_river_key = _signed_log_context(
        target, ("river_leakage", "river_flux"), positive=False
    )
    target_boundary_sink, target_boundary_key = _signed_log_context(
        target, ("head_boundary_flux", "boundary_flux", "constant_head_flux"), positive=False
    )
    source_fields.extend(
        key
        for key in (
            source_pumping_key,
            target_pumping_key,
            source_recharge_key,
            target_river_key,
            target_boundary_key,
        )
        if key is not None
    )

    distance = _finite(attrs.get("distance_m"))
    features: dict[str, float | None] = {
        "distance_m": distance,
        "log_distance_m": math.log1p(distance) if distance is not None else None,
        "weak_geometry_prior": _finite(attrs.get("weak_geometry_prior")),
        "head_delta_m": _finite(attrs.get("head_delta_m")),
        "head_sigma_delta_m": _finite(attrs.get("head_sigma_delta_m")),
        "head_z_score": _finite(attrs.get("head_z_score")),
        "direction_probability": _finite(attrs.get("direction_probability")),
        "gradient_m_per_km": _finite(attrs.get("gradient_m_per_km")),
        "elevation_delta": _finite(attrs.get("elevation_delta")),
        "elevation_direction_probability": _finite(
            attrs.get("elevation_direction_probability")
        ),
        "vertical_separation_m": _finite(attrs.get("vertical_separation_m")),
        "screen_overlap_fraction": _finite(attrs.get("screen_overlap_fraction")),
        "screen_compatibility": _finite(attrs.get("screen_compatibility")),
        "aquifer_compatibility": _finite(attrs.get("aquifer_compatibility")),
        "hydraulic_capacity_proxy": capacity,
        "geology_similarity": geology,
        "chemistry_similarity": chemistry,
        "isotope_similarity": isotopes,
        "tracer_direction_support": tracer_direction_support,
        "source_pumping_strength": source_pumping,
        "target_pumping_sink_strength": target_pumping_sink,
        "source_recharge_strength": source_recharge,
        "target_river_sink_strength": target_river_sink,
        "target_boundary_sink_strength": target_boundary_sink,
        "path_reachability": _finite(attrs.get("path_reachability")),
        "path_relative_support": _finite(attrs.get("path_relative_support")),
        "path_log_support": _finite(attrs.get("path_log_support")),
        "path_steps": _finite(attrs.get("path_steps")),
        "path_transition_efficiency": _finite(
            attrs.get("path_transition_efficiency")
        ),
        "path_activity_support": _finite(attrs.get("path_activity_support")),
        "path_endpoint_sink_support": _finite(
            attrs.get("path_endpoint_sink_support")
        ),
    }
    missing = tuple(sorted(name for name, value in features.items() if value is None))
    return features, tuple(dict.fromkeys(str(item) for item in source_fields))


def build_topology_v2_feature_rows(
    candidate_universe: TopologyV2CandidateUniverse,
    observations: Iterable[Mapping[str, object]] | Mapping[str, Mapping[str, object]],
    *,
    case_id: str = "",
    config: TopologyV2Config | None = None,
) -> tuple[TopologyV2FeatureRow, ...]:
    """Build all v2 feature rows from blind observations only."""

    resolved = config or TopologyV2Config()
    rows = _normalise_observations(observations)
    by_id = {str(row["site_id"]): row for row in rows}
    missing_nodes = set(candidate_universe.nodes) - set(by_id)
    if missing_nodes:
        raise ValueError(f"Observations are missing candidate nodes: {sorted(missing_nodes)}")
    feature_rows: list[TopologyV2FeatureRow] = []
    for edge in candidate_universe.edges:
        values, source_fields = _feature_values(
            edge, by_id[edge.u], by_id[edge.v], resolved
        )
        feature_rows.append(
            TopologyV2FeatureRow(
                edge_id=edge.edge_id,
                case_id=str(case_id),
                u=edge.u,
                v=edge.v,
                features=values,
                missing_features=tuple(sorted(name for name, value in values.items() if value is None)),
                source_fields=tuple(sorted(source_fields)),
                candidate_prior_probability=_finite(edge.attrs.get("candidate_prior_probability")),
            )
        )
    return tuple(feature_rows)


def _row_features(row: TopologyV2FeatureRow | Mapping[str, object]) -> Mapping[str, object]:
    if isinstance(row, TopologyV2FeatureRow):
        return row.features
    features = row.get("features") if isinstance(row, Mapping) else None
    if not isinstance(features, Mapping):
        raise TypeError("Topology-v2 rows must expose a features mapping.")
    return features


def _row_missing(row: TopologyV2FeatureRow | Mapping[str, object], feature_names: Sequence[str]) -> tuple[str, ...]:
    if isinstance(row, TopologyV2FeatureRow):
        return tuple(sorted(set(row.missing_features) & set(feature_names)))
    return tuple(
        sorted(
            name
            for name in feature_names
            if _finite(_row_features(row).get(name)) is None
        )
    )


class TopologyV2LogisticCalibrator:
    """Deterministic, schema-explicit logistic probability calibrator."""

    calibration_version = "topology_v2_logistic_calibrator_v1"

    def __init__(self, *, feature_names: Sequence[str], l2: float = 1.0, max_iter: int = 500):
        names = tuple(str(name) for name in feature_names)
        if not names or len(names) != len(set(names)):
            raise ValueError("A topology-v2 calibrator needs unique feature names.")
        if float(l2) < 0.0 or not math.isfinite(float(l2)):
            raise ValueError("l2 must be finite and non-negative.")
        if int(max_iter) < 20:
            raise ValueError("max_iter must be at least 20.")
        self.feature_names = names
        self.l2 = float(l2)
        self.max_iter = int(max_iter)
        self.means: np.ndarray | None = None
        self.scales: np.ndarray | None = None
        self.coefficients: np.ndarray | None = None
        self.intercept: float | None = None
        self.fitted = False
        self.converged = False
        self.scope = "unfitted"
        self.independent = False
        self.generator_id = ""
        self.split_id = ""
        self.dataset_hash = ""
        self.fit_n = 0
        self.fit_positive_count = 0
        self.supported_missing_patterns: tuple[tuple[str, ...], ...] = ()
        self.constant_features: tuple[str, ...] = ()
        self.optimization_message = ""

    @property
    def deployment_eligible(self) -> bool:
        return bool(
            self.fitted
            and self.converged
            and self.scope in {"held_out_calibration", "deployment"}
            and self.independent
            and self.generator_id
            and self.split_id
            and self.dataset_hash
        )

    def _design_matrix(
        self, rows: Sequence[TopologyV2FeatureRow | Mapping[str, object]], *, fit: bool = False
    ) -> np.ndarray:
        matrix = np.full((len(rows), len(self.feature_names)), np.nan, dtype=float)
        for row_index, row in enumerate(rows):
            features = _row_features(row)
            for column_index, name in enumerate(self.feature_names):
                matrix[row_index, column_index] = _finite(features.get(name))
        if fit:
            means = np.zeros(matrix.shape[1], dtype=float)
            raw_scales = np.zeros(matrix.shape[1], dtype=float)
            for column_index in range(matrix.shape[1]):
                finite_values = matrix[np.isfinite(matrix[:, column_index]), column_index]
                if len(finite_values):
                    means[column_index] = float(np.mean(finite_values))
                    raw_scales[column_index] = float(np.std(finite_values))
            scales = np.where(
                np.isfinite(raw_scales) & (raw_scales > 1e-12),
                raw_scales,
                1.0,
            )
            self.means = means.astype(float)
            self.scales = scales.astype(float)
            self.constant_features = tuple(
                self.feature_names[index]
                for index, scale in enumerate(raw_scales)
                if not math.isfinite(float(scale)) or float(scale) <= 1e-12
            )
        if self.means is None or self.scales is None:
            raise RuntimeError("Topology-v2 calibrator has no imputation statistics.")
        matrix = np.where(np.isfinite(matrix), matrix, self.means[None, :])
        return (matrix - self.means[None, :]) / self.scales[None, :]

    def fit(
        self,
        rows: Sequence[TopologyV2FeatureRow | Mapping[str, object]],
        labels: Sequence[int | float],
        *,
        scope: str,
        independent: bool,
        generator_id: str,
        split_id: str,
        dataset_hash: str,
    ) -> "TopologyV2LogisticCalibrator":
        if not rows:
            raise ValueError("Cannot fit topology-v2 calibration with no rows.")
        y = np.asarray(labels, dtype=float)
        if len(y) != len(rows):
            raise ValueError("Topology-v2 calibration labels and rows have different lengths.")
        if not np.all(np.isfinite(y)) or not np.all(np.isin(y, [0.0, 1.0])):
            raise ValueError("Topology-v2 labels must be finite binary values.")
        if len(np.unique(y)) < 2:
            raise ValueError("Topology-v2 calibration needs both positive and negative labels.")
        X = self._design_matrix(rows, fit=True)
        self.supported_missing_patterns = tuple(
            sorted({
                _row_missing(row, self.feature_names)
                for row in rows
            })
        )
        n_features = X.shape[1]

        def objective(theta: np.ndarray) -> tuple[float, np.ndarray]:
            intercept = float(theta[0])
            beta = theta[1:]
            probabilities = expit(intercept + X @ beta)
            probabilities = np.clip(probabilities, 1e-12, 1.0 - 1e-12)
            loss = float(
                -np.sum(y * np.log(probabilities) + (1.0 - y) * np.log1p(-probabilities))
                / len(y)
                + 0.5 * self.l2 * np.sum(beta**2)
            )
            residual = probabilities - y
            gradient = np.concatenate(
                (
                    np.asarray([float(np.mean(residual))]),
                    (X.T @ residual) / len(y) + self.l2 * beta,
                )
            )
            return loss, gradient

        result = minimize(
            lambda theta: objective(theta)[0],
            np.zeros(n_features + 1, dtype=float),
            jac=lambda theta: objective(theta)[1],
            method="L-BFGS-B",
            options={"maxiter": self.max_iter, "ftol": 1e-12, "gtol": 1e-9, "maxls": 50},
        )
        self.coefficients = np.asarray(result.x[1:], dtype=float)
        self.intercept = float(result.x[0])
        self.fitted = True
        self.converged = bool(result.success and np.all(np.isfinite(result.x)))
        self.scope = str(scope)
        self.independent = bool(independent)
        self.generator_id = str(generator_id)
        self.split_id = str(split_id)
        self.dataset_hash = str(dataset_hash)
        self.fit_n = int(len(y))
        self.fit_positive_count = int(np.sum(y == 1.0))
        self.optimization_message = str(result.message)
        return self

    def predict_proba(
        self, rows: Sequence[TopologyV2FeatureRow | Mapping[str, object]]
    ) -> np.ndarray:
        if not self.fitted or self.coefficients is None or self.intercept is None:
            raise RuntimeError("Topology-v2 calibrator is not fitted.")
        X = self._design_matrix(rows)
        return np.asarray(expit(self.intercept + X @ self.coefficients), dtype=float)

    def to_dict(self) -> dict[str, object]:
        return {
            "calibration_version": self.calibration_version,
            "feature_names": list(self.feature_names),
            "l2": self.l2,
            "max_iter": self.max_iter,
            "means": None if self.means is None else self.means.tolist(),
            "scales": None if self.scales is None else self.scales.tolist(),
            "coefficients": None if self.coefficients is None else self.coefficients.tolist(),
            "intercept": self.intercept,
            "fitted": self.fitted,
            "converged": self.converged,
            "deployment_eligible": self.deployment_eligible,
            "scope": self.scope,
            "independent": self.independent,
            "generator_id": self.generator_id,
            "split_id": self.split_id,
            "dataset_hash": self.dataset_hash,
            "fit_n": self.fit_n,
            "fit_positive_count": self.fit_positive_count,
            "supported_missing_patterns": [list(pattern) for pattern in self.supported_missing_patterns],
            "constant_features": list(self.constant_features),
            "optimization_message": self.optimization_message,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "TopologyV2LogisticCalibrator":
        calibrator = cls(
            feature_names=tuple(str(name) for name in payload.get("feature_names", ())),
            l2=float(payload.get("l2", 1.0)),
            max_iter=int(payload.get("max_iter", 500)),
        )
        for attribute in (
            "means",
            "scales",
            "coefficients",
        ):
            value = payload.get(attribute)
            if value is not None:
                setattr(calibrator, attribute, np.asarray(value, dtype=float))
        calibrator.intercept = (
            None if payload.get("intercept") is None else float(payload["intercept"])
        )
        calibrator.fitted = bool(payload.get("fitted", False))
        calibrator.converged = bool(payload.get("converged", False))
        calibrator.scope = str(payload.get("scope", "unfitted"))
        calibrator.independent = bool(payload.get("independent", False))
        calibrator.generator_id = str(payload.get("generator_id", ""))
        calibrator.split_id = str(payload.get("split_id", ""))
        calibrator.dataset_hash = str(payload.get("dataset_hash", ""))
        calibrator.fit_n = int(payload.get("fit_n", 0))
        calibrator.fit_positive_count = int(payload.get("fit_positive_count", 0))
        calibrator.supported_missing_patterns = tuple(
            tuple(sorted(str(name) for name in pattern))
            for pattern in payload.get("supported_missing_patterns", ())
        )
        calibrator.constant_features = tuple(
            str(name) for name in payload.get("constant_features", ())
        )
        calibrator.optimization_message = str(payload.get("optimization_message", ""))
        return calibrator


@dataclass(frozen=True)
class TopologyV2ThresholdPolicy:
    present_threshold: float = 0.75
    absent_threshold: float = 0.25
    tuning_rule: str = "validation_max_recall_subject_to_fdr"
    target_fdr: float = 0.25
    minimum_gap: float = 0.10

    def __post_init__(self) -> None:
        if not (0.0 < float(self.absent_threshold) < float(self.present_threshold) < 1.0):
            raise ValueError(
                "Topology-v2 selected thresholds must satisfy "
                "0 < absent_threshold < present_threshold < 1."
            )
        if not (0.0 <= float(self.target_fdr) <= 1.0):
            raise ValueError("target_fdr must be between zero and one.")
        if not (0.0 <= float(self.minimum_gap) < 1.0):
            raise ValueError("minimum_gap must be in [0, 1).")

    def to_dict(self) -> dict[str, object]:
        return {
            "present_threshold": self.present_threshold,
            "absent_threshold": self.absent_threshold,
            "tuning_rule": self.tuning_rule,
            "target_fdr": self.target_fdr,
            "minimum_gap": self.minimum_gap,
        }


class TopologyV2Scorer:
    """Apply a frozen calibrated model and emit PRESENT/ABSENT/ABSTAIN."""

    def __init__(
        self,
        calibrator: TopologyV2LogisticCalibrator,
        *,
        policy: TopologyV2ThresholdPolicy | None = None,
        require_deployment_eligibility: bool = True,
    ):
        self.calibrator = calibrator
        self.policy = policy or TopologyV2ThresholdPolicy()
        self.require_deployment_eligibility = bool(require_deployment_eligibility)

    def score_rows(
        self, rows: Sequence[TopologyV2FeatureRow | Mapping[str, object]]
    ) -> tuple[dict[str, object], ...]:
        if self.require_deployment_eligibility and not self.calibrator.deployment_eligible:
            return tuple(
                {
                    "edge_id": getattr(row, "edge_id", None),
                    "case_id": getattr(row, "case_id", None),
                    "probability": None,
                    "decision": DECISION_ABSTAIN,
                    "missing_features": list(_row_missing(row, self.calibrator.feature_names)),
                    "abstention_reason": "calibration_not_deployment_eligible",
                }
                for row in rows
            )
        probabilities = self.calibrator.predict_proba(rows)
        records: list[dict[str, object]] = []
        supported_patterns = set(self.calibrator.supported_missing_patterns)
        for row, probability in zip(rows, probabilities):
            missing = _row_missing(row, self.calibrator.feature_names)
            edge_id = getattr(row, "edge_id", None)
            case_id = getattr(row, "case_id", None)
            if isinstance(row, Mapping):
                edge_id = row.get("edge_id", edge_id)
                case_id = row.get("case_id", case_id)
            if tuple(missing) not in supported_patterns:
                records.append(
                    {
                        "edge_id": edge_id,
                        "case_id": case_id,
                        "probability": None,
                        "decision": DECISION_ABSTAIN,
                        "missing_features": list(missing),
                        "abstention_reason": "unsupported_missingness_pattern",
                    }
                )
                continue
            value = float(probability)
            if value >= self.policy.present_threshold:
                decision = DECISION_PRESENT
            elif value <= self.policy.absent_threshold:
                decision = DECISION_ABSENT
            else:
                decision = DECISION_ABSTAIN
            records.append(
                {
                    "edge_id": edge_id,
                    "case_id": case_id,
                    "probability": value,
                    "decision": decision,
                    "missing_features": list(missing),
                    "abstention_reason": "probability_between_thresholds"
                    if decision == DECISION_ABSTAIN
                    else None,
                }
            )
        return tuple(records)


def _normalise_edge_set(edges: Iterable[object]) -> set[str]:
    result: set[str] = set()
    for edge in edges:
        if isinstance(edge, str) and "->" in edge:
            left, right = edge.split("->", 1)
        elif isinstance(edge, Mapping):
            left = edge.get("u", edge.get("source"))
            right = edge.get("v", edge.get("target"))
            if left is None or right is None:
                edge_id = edge.get("edge_id")
                if edge_id is None or "->" not in str(edge_id):
                    continue
                left, right = str(edge_id).split("->", 1)
        elif isinstance(edge, (tuple, list)) and len(edge) >= 2:
            left, right = edge[0], edge[1]
        else:
            left = getattr(edge, "u", None)
            right = getattr(edge, "v", None)
        if left is not None and right is not None and str(left) != str(right):
            result.add(f"{str(left)}->{str(right)}")
    return result


def _average_precision(labels: np.ndarray, probabilities: np.ndarray) -> float | None:
    positives = int(np.sum(labels == 1.0))
    if positives == 0:
        return None
    order = np.argsort(-probabilities, kind="mergesort")
    sorted_labels = labels[order]
    cumulative = np.cumsum(sorted_labels == 1.0)
    ranks = np.arange(1, len(labels) + 1, dtype=float)
    precision = cumulative / ranks
    return float(np.sum(precision[sorted_labels == 1.0]) / positives)


def _roc_auc(labels: np.ndarray, probabilities: np.ndarray) -> float | None:
    positives = labels == 1.0
    negatives = labels == 0.0
    n_positive = int(np.sum(positives))
    n_negative = int(np.sum(negatives))
    if n_positive == 0 or n_negative == 0:
        return None
    order = np.argsort(probabilities, kind="mergesort")
    sorted_probabilities = probabilities[order]
    ranks = np.empty(len(probabilities), dtype=float)
    index = 0
    while index < len(order):
        end = index + 1
        while end < len(order) and sorted_probabilities[end] == sorted_probabilities[index]:
            end += 1
        ranks[order[index:end]] = 0.5 * (index + 1 + end)
        index = end
    rank_sum = float(np.sum(ranks[positives]))
    return float((rank_sum - n_positive * (n_positive + 1) / 2.0) / (n_positive * n_negative))


def _reliability(
    labels: np.ndarray, probabilities: np.ndarray, n_bins: int = 10
) -> tuple[dict[str, object], ...]:
    bins: list[dict[str, object]] = []
    for index in range(int(n_bins)):
        lower = index / float(n_bins)
        upper = (index + 1) / float(n_bins)
        mask = (probabilities >= lower) & (
            probabilities < upper if index < n_bins - 1 else probabilities <= upper
        )
        count = int(np.sum(mask))
        bins.append(
            {
                "bin_lower": lower,
                "bin_upper": upper,
                "count": count,
                "mean_probability": float(np.mean(probabilities[mask])) if count else None,
                "observed_fraction": float(np.mean(labels[mask])) if count else None,
            }
        )
    return tuple(bins)


def _mcc(tp: int, tn: int, fp: int, fn: int) -> float | None:
    denominator = math.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator <= 0.0:
        return None
    return float((tp * tn - fp * fn) / denominator)


def topology_v2_metrics(
    labels: Sequence[int | float],
    probabilities: Sequence[float],
    *,
    policy: TopologyV2ThresholdPolicy | None = None,
    candidate_recall: float | None = None,
    n_reliability_bins: int = 10,
) -> dict[str, object]:
    """Return probability, tri-state, and candidate-universe metrics."""

    resolved = policy or TopologyV2ThresholdPolicy()
    y = np.asarray(labels, dtype=float)
    p = np.asarray(probabilities, dtype=float)
    if len(y) != len(p) or len(y) == 0:
        raise ValueError("Topology-v2 metric labels and probabilities must be non-empty and aligned.")
    if not np.all(np.isin(y, [0.0, 1.0])) or not np.all(np.isfinite(p)):
        raise ValueError("Topology-v2 metrics require finite binary labels and probabilities.")
    p = np.clip(p, 0.0, 1.0)
    present = p >= resolved.present_threshold
    absent = p <= resolved.absent_threshold
    abstain = ~(present | absent)
    predicted_present = present
    predicted_absent = absent
    tp = int(np.sum(predicted_present & (y == 1.0)))
    fp = int(np.sum(predicted_present & (y == 0.0)))
    fn = int(np.sum(predicted_absent & (y == 1.0)))
    tn = int(np.sum(predicted_absent & (y == 0.0)))
    abstain_positive = int(np.sum(abstain & (y == 1.0)))
    abstain_negative = int(np.sum(abstain & (y == 0.0)))
    n = int(len(y))
    positive_count = int(np.sum(y == 1.0))
    negative_count = int(np.sum(y == 0.0))
    predicted_present_count = int(np.sum(predicted_present))
    resolved_count = int(np.sum(~abstain))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / positive_count if positive_count else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    fdr = fp / (tp + fp) if tp + fp else 0.0
    selective_precision = tp / (tp + fp) if tp + fp else 0.0
    selective_recall = tp / (tp + fn) if tp + fn else 0.0
    brier = float(np.mean((p - y) ** 2))
    log_loss = float(
        -np.mean(y * np.log(np.clip(p, 1e-12, 1.0)) + (1.0 - y) * np.log(np.clip(1.0 - p, 1e-12, 1.0)))
    )
    reliability = _reliability(y, p, n_bins=n_reliability_bins)
    ece = float(
        sum(
            (float(item["count"]) / n)
            * abs(float(item["mean_probability"]) - float(item["observed_fraction"]))
            for item in reliability
            if item["count"]
        )
    )
    return {
        "n_candidates": n,
        "n_positive": positive_count,
        "n_negative": negative_count,
        "pr_auc": _average_precision(y, p),
        "roc_auc": _roc_auc(y, p),
        "brier": brier,
        "log_loss": log_loss,
        "ece": ece,
        "reliability": [dict(item) for item in reliability],
        "present_threshold": resolved.present_threshold,
        "absent_threshold": resolved.absent_threshold,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "abstain_positive": abstain_positive,
        "abstain_negative": abstain_negative,
        "abstain_count": int(np.sum(abstain)),
        "abstention_rate": float(np.mean(abstain)),
        "coverage": float(resolved_count / n),
        "predicted_present_count": predicted_present_count,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "fdr": fdr,
        "selective_precision": selective_precision,
        "selective_recall": selective_recall,
        "mcc_resolved": _mcc(tp, tn, fp, fn),
        "candidate_graph_recall": candidate_recall,
    }


def tune_topology_thresholds(
    labels: Sequence[int | float],
    probabilities: Sequence[float],
    *,
    target_fdr: float = 0.25,
    minimum_gap: float = 0.10,
    grid_size: int = 19,
) -> tuple[TopologyV2ThresholdPolicy, tuple[dict[str, object], ...]]:
    """Tune a tri-state policy on validation data only.

    The primary selection rule is maximum recall subject to the declared FDR
    target.  An all-abstain/all-absent policy is not accepted as a meaningful
    selected topology when a non-empty PRESENT operating point exists.  If the
    FDR constraint is infeasible for non-empty selections, the fallback
    minimizes FDR and then maximizes recall, and records that the target was
    infeasible.  Thresholds of exactly zero or one are never in the grid.
    """

    if not (0.0 <= float(target_fdr) <= 1.0):
        raise ValueError("target_fdr must be between zero and one.")
    if not (0.0 <= float(minimum_gap) < 1.0):
        raise ValueError("minimum_gap must be in [0, 1).")
    grid = np.linspace(0.05, 0.95, int(grid_size))
    candidates: list[dict[str, object]] = []
    for present in grid:
        for absent in grid:
            if absent >= present or present - absent + 1e-12 < float(minimum_gap):
                continue
            policy = TopologyV2ThresholdPolicy(
                present_threshold=float(present),
                absent_threshold=float(absent),
                target_fdr=float(target_fdr),
                minimum_gap=float(minimum_gap),
            )
            metrics = topology_v2_metrics(labels, probabilities, policy=policy)
            candidates.append(
                {
                    "policy": policy.to_dict(),
                    "precision": metrics["precision"],
                    "recall": metrics["recall"],
                    "fdr": metrics["fdr"],
                    "f1": metrics["f1"],
                    "coverage": metrics["coverage"],
                    "abstention_rate": metrics["abstention_rate"],
                    "predicted_present_count": metrics["predicted_present_count"],
                    "meets_fdr_target": float(metrics["fdr"]) <= float(target_fdr),
                }
            )
    if not candidates:
        raise ValueError("Threshold grid produced no valid tri-state policies.")
    nonempty = [item for item in candidates if int(item["predicted_present_count"]) > 0]
    feasible = [
        item
        for item in nonempty
        if bool(item["meets_fdr_target"])
    ]
    if feasible:
        pool = feasible
        selected = max(
            pool,
            key=lambda item: (
                float(item["recall"]),
                float(item["precision"]),
                float(item["coverage"]),
                -float(item["abstention_rate"]),
            ),
        )
    elif nonempty:
        pool = nonempty
        selected = min(
            pool,
            key=lambda item: (
                float(item["fdr"]),
                -float(item["recall"]),
                -float(item["precision"]),
            ),
        )
    else:
        # This is only possible when every validation probability is below the
        # interior grid.  Keep the result explicit rather than manufacturing a
        # selected edge.
        selected = min(
            candidates,
            key=lambda item: (
                float(item["fdr"]),
                -float(item["recall"]),
                -float(item["precision"]),
            ),
        )
    selected_policy = TopologyV2ThresholdPolicy(**dict(selected["policy"]))
    return selected_policy, tuple(candidates)


def bootstrap_case_metric_ci(
    labels: Sequence[int | float],
    probabilities: Sequence[float],
    case_ids: Sequence[str],
    *,
    metric: str = "pr_auc",
    policy: TopologyV2ThresholdPolicy | None = None,
    n_bootstrap: int = 1000,
    seed: int = 20260922,
) -> dict[str, object]:
    """Bootstrap a metric by whole case, never by individual edge."""

    y = np.asarray(labels, dtype=float)
    p = np.asarray(probabilities, dtype=float)
    cases = np.asarray([str(case) for case in case_ids], dtype=object)
    unique_cases = tuple(sorted(set(str(case) for case in cases)))
    if len(y) != len(p) or len(y) != len(cases) or not unique_cases:
        raise ValueError("Case bootstrap inputs must be aligned and non-empty.")
    rng = np.random.default_rng(int(seed))
    values: list[float] = []
    for _ in range(int(n_bootstrap)):
        sampled_cases = rng.choice(unique_cases, size=len(unique_cases), replace=True)
        mask = np.concatenate([np.flatnonzero(cases == case) for case in sampled_cases])
        sample_metrics = topology_v2_metrics(y[mask], p[mask], policy=policy)
        value = sample_metrics.get(metric)
        if value is not None and math.isfinite(float(value)):
            values.append(float(value))
    if not values:
        return {
            "metric": metric,
            "n_cases": len(unique_cases),
            "n_bootstrap": int(n_bootstrap),
            "resampling_unit": "whole_case",
            "estimate": None,
            "ci95_low": None,
            "ci95_high": None,
        }
    values_array = np.asarray(values, dtype=float)
    return {
        "metric": metric,
        "n_cases": len(unique_cases),
        "n_bootstrap": int(n_bootstrap),
        "n_valid_bootstrap": len(values),
        "resampling_unit": "whole_case",
        "estimate": float(np.mean(values_array)),
        "ci95_low": float(np.quantile(values_array, 0.025)),
        "ci95_high": float(np.quantile(values_array, 0.975)),
    }


__all__ = [
    "CAPACITY_FEATURES",
    "CHEMISTRY_TRACER_FEATURES",
    "DECISION_ABSENT",
    "DECISION_ABSTAIN",
    "DECISION_PRESENT",
    "ELEVATION_FEATURES",
    "GEOMETRY_FEATURES",
    "GEOLOGY_FEATURES",
    "HEAD_FEATURES",
    "HEAD_GRADIENT_FEATURES",
    "PATH_FEATURES",
    "MODEL_FEATURE_SETS",
    "SCREEN_FEATURES",
    "SOURCE_SINK_FEATURES",
    "TOPOLOGY_V2_ALGORITHM",
    "TOPOLOGY_V2_VERSION",
    "TopologyV2CandidateEdge",
    "TopologyV2CandidateUniverse",
    "TopologyV2Config",
    "TopologyV2FeatureRow",
    "TopologyV2LogisticCalibrator",
    "TopologyV2Scorer",
    "TopologyV2ThresholdPolicy",
    "bootstrap_case_metric_ci",
    "build_topology_v2_feature_rows",
    "generate_topology_v2_candidate_universe",
    "topology_v2_metrics",
    "tune_topology_thresholds",
]
