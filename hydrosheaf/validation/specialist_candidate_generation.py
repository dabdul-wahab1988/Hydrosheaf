"""Independent candidate generation for specialist comparator runs.

The HydroSheaf inference path is deliberately not imported here.  This module
constructs a small, deterministic directed candidate graph from observed site
geometry and measured heads only.  It exists to prevent a comparator from
being evaluated only on edges that HydroSheaf already proposed.

The graph is a k-nearest-neighbour geometry graph.  Each undirected neighbour
pair is oriented from higher observed head to lower observed head.  Tied or
missing heads retain both directions so that the generator does not silently
turn missing evidence into a directional claim.  The resulting support value
is a raw geometry/head heuristic, not a calibrated probability and not a
HydroSheaf posterior.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

from .baselines import assert_truth_blind_observations


_ALGORITHM = "independent_geometry_head_knn_v1"
_VERSION = "1.0"
_EARTH_METRES_PER_DEGREE_LAT = 110_540.0
_EARTH_METRES_PER_DEGREE_LON = 111_320.0


def _finite(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _jsonable(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True, default=str))


def _hash_payload(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, float(value)))


@dataclass(frozen=True)
class IndependentCandidateEdge:
    """One edge generated without HydroSheaf candidate or truth fields."""

    u: str
    v: str
    attrs: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        source = str(self.u).strip()
        target = str(self.v).strip()
        if not source or not target:
            raise ValueError("Independent candidate edges require non-empty node IDs.")
        if source == target:
            raise ValueError("Independent candidate edges cannot be self-loops.")
        object.__setattr__(self, "u", source)
        object.__setattr__(self, "v", target)
        object.__setattr__(self, "attrs", dict(self.attrs))

    @property
    def edge(self) -> tuple[str, str]:
        return self.u, self.v

    def to_audit_record(self) -> dict[str, object]:
        return {"edge": [self.u, self.v], "attrs": _jsonable(self.attrs)}


@dataclass(frozen=True)
class IndependentCandidateUniverse:
    """Auditable candidate graph and its blind-input provenance."""

    edges: tuple[IndependentCandidateEdge, ...]
    nodes: tuple[str, ...]
    coordinate_nodes: tuple[str, ...]
    dropped_nodes: tuple[str, ...]
    algorithm: str
    version: str
    input_channels: tuple[str, ...]
    parameters: Mapping[str, object]
    input_hash: str
    candidate_hash: str
    coordinate_system: str
    truth_blind: bool = True
    truth_fields_seen: tuple[str, ...] = ()
    source: str = "hydrosheaf.validation.specialist_candidate_generation"

    def __post_init__(self) -> None:
        if not self.truth_blind:
            raise ValueError("Independent candidate generation must be truth-blind.")
        if self.truth_fields_seen:
            raise ValueError("Truth-blind candidate generation cannot record truth fields.")
        if not self.algorithm or not self.version or not self.input_hash:
            raise ValueError("Independent candidate provenance is incomplete.")
        expected_edges = tuple(sorted(self.edges, key=lambda edge: edge.edge))
        if expected_edges != self.edges:
            raise ValueError("Independent candidate edges must be sorted deterministically.")
        pairs = [edge.edge for edge in self.edges]
        if len(pairs) != len(set(pairs)):
            raise ValueError("Independent candidate universe contains duplicate edges.")

    def to_audit_record(self) -> dict[str, object]:
        return {
            "algorithm": self.algorithm,
            "version": self.version,
            "source": self.source,
            "input_channels": list(self.input_channels),
            "parameters": _jsonable(self.parameters),
            "coordinate_system": self.coordinate_system,
            "nodes": list(self.nodes),
            "coordinate_nodes": list(self.coordinate_nodes),
            "dropped_nodes": list(self.dropped_nodes),
            "edge_count": len(self.edges),
            "edges": [edge.to_audit_record() for edge in self.edges],
            "input_hash": self.input_hash,
            "candidate_hash": self.candidate_hash,
            "truth_blind": self.truth_blind,
            "truth_fields_seen": list(self.truth_fields_seen),
        }


@dataclass(frozen=True)
class _Node:
    site_id: str
    x_m: float
    y_m: float
    head_m: float | None
    coordinate_source: str


def generate_independent_candidate_universe(
    observations: Iterable[Mapping[str, object]],
    *,
    max_neighbors: int = 4,
    max_distance_km: float | None = None,
    head_tie_tolerance_m: float = 0.10,
) -> IndependentCandidateUniverse:
    """Generate a deterministic candidate universe from blind observations.

    Only ``site_id``/``sample_id``, coordinates, and the observed head fields
    are used.  Generator truth, HydroSheaf edge attributes, posterior fields,
    and process labels are rejected or ignored at this boundary.  The input
    hash is computed from the declared fields after sorting by site ID, so the
    result is invariant to observation-row order.
    """

    if int(max_neighbors) < 1:
        raise ValueError("max_neighbors must be at least one.")
    max_neighbors = int(max_neighbors)
    if max_distance_km is not None and (
        not math.isfinite(float(max_distance_km)) or float(max_distance_km) <= 0.0
    ):
        raise ValueError("max_distance_km must be finite and positive when supplied.")
    if not math.isfinite(float(head_tie_tolerance_m)) or float(head_tie_tolerance_m) < 0.0:
        raise ValueError("head_tie_tolerance_m must be finite and non-negative.")
    rows = [dict(row) for row in observations]
    assert_truth_blind_observations(rows)
    if not rows:
        return _empty_universe(
            max_neighbors=max_neighbors,
            max_distance_km=max_distance_km,
            head_tie_tolerance_m=float(head_tie_tolerance_m),
        )

    raw_records: list[dict[str, object]] = []
    for row in rows:
        raw_site = row.get("site_id", row.get("sample_id"))
        site_id = str(raw_site).strip() if raw_site is not None else ""
        if not site_id:
            raise ValueError("Every blind observation needs site_id or sample_id.")
        x_m = _finite(row.get("x_m"))
        y_m = _finite(row.get("y_m"))
        lat = _finite(row.get("lat"))
        lon = _finite(row.get("lon"))
        if (x_m is None) != (y_m is None):
            raise ValueError(f"Site {site_id!r} has only one Cartesian coordinate.")
        if (lat is None) != (lon is None):
            raise ValueError(f"Site {site_id!r} has only one geographic coordinate.")
        if x_m is None and (lat is None or lon is None):
            coordinate_source = "missing"
        elif x_m is not None:
            coordinate_source = "x_m_y_m"
        else:
            coordinate_source = "lon_lat"
        head = _finite(row.get("head_meas"))
        head_source = "head_meas"
        if head is None:
            head = _finite(row.get("hydraulic_head"))
            head_source = "hydraulic_head" if head is not None else "missing"
        raw_records.append(
            {
                "site_id": site_id,
                "x_m": x_m,
                "y_m": y_m,
                "lat": lat,
                "lon": lon,
                "head_meas": head,
                "coordinate_source": coordinate_source,
                "head_source": head_source,
            }
        )
    raw_records.sort(key=lambda item: str(item["site_id"]))
    site_ids = [str(item["site_id"]) for item in raw_records]
    if len(site_ids) != len(set(site_ids)):
        raise ValueError("Blind observations contain duplicate site IDs.")

    input_hash = _hash_payload(raw_records)
    coordinate_sources = {
        str(item["coordinate_source"])
        for item in raw_records
        if str(item["coordinate_source"]) != "missing"
    }
    if len(coordinate_sources) > 1:
        raise ValueError(
            "Independent candidate generation cannot combine Cartesian and "
            "geographic coordinates without a declared transform."
        )
    coordinate_system = next(iter(coordinate_sources), "none")
    nodes = _project_nodes(raw_records, coordinate_system)
    by_id = {node.site_id: node for node in nodes}
    dropped_nodes = tuple(
        sorted(site_id for site_id in site_ids if site_id not in by_id)
    )
    candidate_pairs = _nearest_pairs(
        tuple(nodes),
        max_neighbors=max_neighbors,
        max_distance_km=max_distance_km,
    )
    distance_scale = _distance_scale(candidate_pairs, by_id)
    edges: list[IndependentCandidateEdge] = []
    for left, right in sorted(candidate_pairs):
        left_node = by_id[left]
        right_node = by_id[right]
        distance_m = math.hypot(
            left_node.x_m - right_node.x_m,
            left_node.y_m - right_node.y_m,
        )
        left_head = left_node.head_m
        right_head = right_node.head_m
        head_drop = (
            left_head - right_head
            if left_head is not None and right_head is not None
            else None
        )
        if head_drop is None or abs(head_drop) <= float(head_tie_tolerance_m):
            directions = ((left, right), (right, left))
            reason = (
                "missing_observed_head_both_directions"
                if head_drop is None
                else "head_tie_both_directions"
            )
        elif head_drop > 0.0:
            directions = ((left, right),)
            reason = "higher_observed_head_to_lower"
        else:
            directions = ((right, left),)
            reason = "higher_observed_head_to_lower"
        for source, target in directions:
            source_node = by_id[source]
            target_node = by_id[target]
            directed_drop = (
                source_node.head_m - target_node.head_m
                if source_node.head_m is not None and target_node.head_m is not None
                else None
            )
            distance_km = distance_m / 1000.0
            distance_score = 1.0 / (1.0 + distance_m / distance_scale)
            if directed_drop is None:
                directional_score = 0.5
            else:
                scale = max(0.10, float(head_tie_tolerance_m))
                directional_score = 1.0 / (1.0 + math.exp(-directed_drop / scale))
            independent_support = _clamp(
                0.5 * distance_score + 0.5 * directional_score,
                0.0,
                1.0,
            )
            attrs: dict[str, object] = {
                "candidate_generator": _ALGORITHM,
                "direction_rule": reason,
                "distance_km": distance_km,
                "independent_support": independent_support,
            }
            if directed_drop is not None:
                attrs["head_drop"] = directed_drop
                attrs["head_gradient_m_per_km"] = (
                    directed_drop / distance_km if distance_km > 0.0 else directed_drop
                )
            edges.append(IndependentCandidateEdge(source, target, attrs))
    edges.sort(key=lambda edge: edge.edge)
    candidate_hash = _hash_payload([edge.to_audit_record() for edge in edges])
    return IndependentCandidateUniverse(
        edges=tuple(edges),
        nodes=tuple(site_ids),
        coordinate_nodes=tuple(sorted(by_id)),
        dropped_nodes=dropped_nodes,
        algorithm=_ALGORITHM,
        version=_VERSION,
        input_channels=("site_id", "x_m", "y_m", "lat", "lon", "head_meas"),
        parameters={
            "max_neighbors": max_neighbors,
            "max_distance_km": max_distance_km,
            "head_tie_tolerance_m": float(head_tie_tolerance_m),
        },
        input_hash=input_hash,
        candidate_hash=candidate_hash,
        coordinate_system=coordinate_system,
    )


def _empty_universe(
    *,
    max_neighbors: int,
    max_distance_km: float | None,
    head_tie_tolerance_m: float,
) -> IndependentCandidateUniverse:
    input_hash = _hash_payload([])
    return IndependentCandidateUniverse(
        edges=(),
        nodes=(),
        coordinate_nodes=(),
        dropped_nodes=(),
        algorithm=_ALGORITHM,
        version=_VERSION,
        input_channels=("site_id", "x_m", "y_m", "lat", "lon", "head_meas"),
        parameters={
            "max_neighbors": max_neighbors,
            "max_distance_km": max_distance_km,
            "head_tie_tolerance_m": head_tie_tolerance_m,
        },
        input_hash=input_hash,
        candidate_hash=_hash_payload([]),
        coordinate_system="none",
    )


def _project_nodes(
    records: Sequence[Mapping[str, object]],
    coordinate_system: str,
) -> tuple[_Node, ...]:
    if coordinate_system == "x_m_y_m":
        return tuple(
            _Node(
                str(record["site_id"]),
                float(record["x_m"]),
                float(record["y_m"]),
                _finite(record.get("head_meas")),
                "x_m_y_m",
            )
            for record in records
            if record.get("x_m") is not None and record.get("y_m") is not None
        )
    if coordinate_system == "lon_lat":
        geo_records = [
            record
            for record in records
            if record.get("lon") is not None and record.get("lat") is not None
        ]
        if not geo_records:
            return ()
        lat0 = sum(float(record["lat"]) for record in geo_records) / len(geo_records)
        lon0 = sum(float(record["lon"]) for record in geo_records) / len(geo_records)
        cos_lat = max(0.1, math.cos(math.radians(lat0)))
        return tuple(
            _Node(
                str(record["site_id"]),
                (float(record["lon"]) - lon0)
                * _EARTH_METRES_PER_DEGREE_LON
                * cos_lat,
                (float(record["lat"]) - lat0) * _EARTH_METRES_PER_DEGREE_LAT,
                _finite(record.get("head_meas")),
                "lon_lat",
            )
            for record in geo_records
        )
    return ()


def _nearest_pairs(
    nodes: Sequence[_Node],
    *,
    max_neighbors: int,
    max_distance_km: float | None,
) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for node in nodes:
        neighbours: list[tuple[float, str]] = []
        for other in nodes:
            if node.site_id == other.site_id:
                continue
            distance_m = math.hypot(node.x_m - other.x_m, node.y_m - other.y_m)
            if max_distance_km is not None and distance_m > float(max_distance_km) * 1000.0:
                continue
            neighbours.append((distance_m, other.site_id))
        neighbours.sort(key=lambda item: (item[0], item[1]))
        for _, other_id in neighbours[:max_neighbors]:
            pairs.add(tuple(sorted((node.site_id, other_id))))
    return pairs


def _distance_scale(
    pairs: Iterable[tuple[str, str]],
    by_id: Mapping[str, _Node],
) -> float:
    distances = [
        math.hypot(by_id[left].x_m - by_id[right].x_m, by_id[left].y_m - by_id[right].y_m)
        for left, right in pairs
    ]
    positive = sorted(distance for distance in distances if distance > 0.0)
    if not positive:
        return 1.0
    middle = len(positive) // 2
    if len(positive) % 2:
        return positive[middle]
    return 0.5 * (positive[middle - 1] + positive[middle])


__all__ = [
    "IndependentCandidateEdge",
    "IndependentCandidateUniverse",
    "generate_independent_candidate_universe",
]
