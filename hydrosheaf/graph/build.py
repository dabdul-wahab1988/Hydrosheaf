"""Graph construction helpers."""

import math
from typing import Iterable, List, Mapping, Optional, Tuple, Union

from .types import Edge


EdgeInput = Union[Tuple[str, str], Tuple[str, str, str], Edge]


def build_edges(edges: Iterable[EdgeInput]) -> List[Edge]:
    built: List[Edge] = []
    for entry in edges:
        if isinstance(entry, Edge):
            built.append(entry)
            continue
        if len(entry) == 2:
            u, v = entry
            edge_id = f"{u}->{v}"
        elif len(entry) == 3:
            edge_id, u, v = entry
        else:
            raise ValueError("Edge entry must be (u, v), (edge_id, u, v), or Edge.")
        built.append(Edge(edge_id=edge_id, u=u, v=v))
    return built


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    radius_km = 6371.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    d_phi = math.radians(lat2 - lat1)
    d_lambda = math.radians(lon2 - lon1)
    a = math.sin(d_phi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2.0) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return radius_km * c


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def _head_estimate(
    sample: Mapping[str, object],
    head_key: str,
    dtw_key: str,
    elevation_key: str,
    sigma_meas: float,
    sigma_dtw: float,
    sigma_elev: float,
    sigma_topo: float,
) -> Tuple[Optional[float], Optional[float], str]:
    if sample.get(head_key) not in (None, ""):
        return float(sample[head_key]), sigma_meas, "A"
    if sample.get(dtw_key) not in (None, "") and sample.get(elevation_key) not in (None, ""):
        head = float(sample[elevation_key]) - float(sample[dtw_key])
        sigma = math.sqrt(sigma_elev**2 + sigma_dtw**2)
        return head, sigma, "B"
    if sample.get(elevation_key) not in (None, ""):
        return float(sample[elevation_key]), sigma_topo, "C"
    return None, None, "missing"


def infer_edges_probabilistic(
    samples: Iterable[Mapping[str, object]],
    radius_km: float,
    max_neighbors: int,
    p_min: float,
    head_key: str = "head_meas",
    dtw_key: str = "dtw",
    elevation_key: str = "elevation",
    aquifer_key: str = "aquifer_unit",
    screen_depth_key: str = "screen_depth",
    well_depth_key: str = "well_depth",
    sigma_meas: float = 0.5,
    sigma_dtw: float = 1.0,
    sigma_elev: float = 1.0,
    sigma_topo: float = 10.0,
    gradient_min: float = 1e-4,
    depth_mismatch: float = 20.0,
) -> List[Edge]:
    samples_list = list(samples)
    node_rows: List[Tuple[str, float, float, Mapping[str, object], float, float, str]] = []
    for sample in samples_list:
        site_id = sample.get("site_id")
        if site_id is None:
            continue
        if sample.get("lat") in (None, "") or sample.get("lon") in (None, ""):
            continue
        head, sigma, tier = _head_estimate(
            sample,
            head_key,
            dtw_key,
            elevation_key,
            sigma_meas,
            sigma_dtw,
            sigma_elev,
            sigma_topo,
        )
        if head is None or sigma is None:
            continue
        node_rows.append(
            (
                str(site_id),
                float(sample["lat"]),
                float(sample["lon"]),
                sample,
                float(head),
                float(sigma),
                tier,
            )
        )

    edges: List[Edge] = []
    edge_ids = set()
    for node_id, lat, lon, sample, head, sigma, tier in node_rows:
        candidates: List[Tuple[float, float, str, Mapping[str, object], float, float, str, List[str]]] = []
        for other_id, o_lat, o_lon, other_sample, other_head, other_sigma, other_tier in node_rows:
            if other_id == node_id:
                continue
            if aquifer_key and sample.get(aquifer_key) and other_sample.get(aquifer_key):
                if sample.get(aquifer_key) != other_sample.get(aquifer_key):
                    continue
            distance = _haversine_km(lat, lon, o_lat, o_lon)
            if radius_km and distance > radius_km:
                continue
            if distance == 0:
                continue

            delta_h = head - other_head
            sigma_delta = math.sqrt(sigma**2 + other_sigma**2)
            if sigma_delta == 0:
                p_uv = 0.5 if delta_h == 0 else (1.0 if delta_h > 0 else 0.0)
            else:
                p_uv = _normal_cdf(delta_h / sigma_delta)

            gradient = abs(delta_h) / (distance * 1000.0)
            flags: List[str] = []
            if gradient_min > 0 and gradient < gradient_min:
                flags.append("flat_gradient")
                factor = gradient / gradient_min
                p_uv = 0.5 + (p_uv - 0.5) * factor

            depth_a = sample.get(screen_depth_key) or sample.get(well_depth_key)
            depth_b = other_sample.get(screen_depth_key) or other_sample.get(well_depth_key)
            if depth_a not in (None, "") and depth_b not in (None, ""):
                if abs(float(depth_a) - float(depth_b)) > depth_mismatch:
                    flags.append("depth_mismatch")

            if p_uv < p_min:
                continue

            candidates.append(
                (
                    p_uv,
                    distance,
                    other_id,
                    other_sample,
                    delta_h,
                    sigma_delta,
                    f"{tier}/{other_tier}",
                    flags,
                )
            )

        candidates.sort(key=lambda item: (-item[0], item[1]))
        if max_neighbors > 0:
            candidates = candidates[:max_neighbors]

        for p_uv, distance, other_id, other_sample, delta_h, sigma_delta, tier_pair, flags in candidates:
            edge_id = f"{node_id}->{other_id}"
            if edge_id in edge_ids:
                continue
            edges.append(
                Edge(
                    edge_id=edge_id,
                    u=node_id,
                    v=other_id,
                    attrs={
                        "distance_km": distance,
                        "delta_h": delta_h,
                        "sigma_delta_h": sigma_delta,
                        "p_uv": p_uv,
                        "edge_confidence": p_uv,
                        "source_tier": tier_pair,
                        "flags": ",".join(flags),
                    },
                )
            )
            edge_ids.add(edge_id)

    return edges


def infer_edges_from_coordinates(
    samples: Iterable[Mapping[str, object]],
    max_neighbors: int = 1,
    allow_uphill: bool = False,
    flow_to_key: str = "flow_to",
    head_key: str = "hydraulic_head",
    elevation_key: str = "elevation",
) -> List[Edge]:
    samples_list = list(samples)
    sample_map = {str(sample.get("site_id")): sample for sample in samples_list if sample.get("site_id")}

    nodes: List[Tuple[str, float, float, Optional[float]]] = []
    for sample in samples_list:
        site_id = sample.get("site_id")
        if site_id is None:
            continue
        if "lat" not in sample or "lon" not in sample:
            continue
        elevation = sample.get(head_key)
        if elevation is None:
            elevation = sample.get(elevation_key)
        nodes.append((str(site_id), float(sample["lat"]), float(sample["lon"]), elevation))

    node_lookup = {node_id: (lat, lon, elevation) for node_id, lat, lon, elevation in nodes}
    edges: List[Edge] = []
    edge_keys = set()

    for node_id, lat, lon, elevation in nodes:
        sample = sample_map.get(node_id)
        if sample and flow_to_key in sample:
            target = str(sample[flow_to_key])
            if target in node_lookup:
                edge_id = f"{node_id}->{target}"
                if edge_id not in edge_keys:
                    edges.append(Edge(edge_id=edge_id, u=node_id, v=target))
                    edge_keys.add(edge_id)
                continue

        candidates: List[Tuple[float, str]] = []
        for other_id, o_lat, o_lon, o_elev in nodes:
            if other_id == node_id:
                continue
            if elevation is not None and o_elev is not None and not allow_uphill:
                if o_elev >= elevation:
                    continue
            distance = _haversine_km(lat, lon, o_lat, o_lon)
            candidates.append((distance, other_id))
        candidates.sort(key=lambda item: item[0])
        for _, target in candidates[:max_neighbors]:
            edge_id = f"{node_id}->{target}"
            if edge_id in edge_keys:
                continue
            edges.append(Edge(edge_id=edge_id, u=node_id, v=target))
            edge_keys.add(edge_id)

    return edges
