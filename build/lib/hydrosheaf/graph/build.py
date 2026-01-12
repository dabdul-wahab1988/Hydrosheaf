"""Graph construction helpers."""

import math
from typing import Dict, Iterable, List, Mapping, Optional, Tuple, Union

from .types import Edge
from .head_inference import (
    infer_head_heuristic,
    infer_heads_bayesian_linear,
    infer_heads_bayesian_mcmc,
)


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
    head_inference: str = "heuristic",
    dtw_prior_mu: float = 5.0,
    dtw_prior_sigma: float = 5.0,
    head_prior_mu: float = 0.0,
    head_prior_sigma: float = 1000.0,
    mcmc_draws: int = 1000,
    mcmc_chains: int = 2,
    mcmc_target_accept: float = 0.9,
    mcmc_warmup_fraction: float = 0.5,
) -> List[Edge]:
    samples_list = list(samples)
    node_rows: List[Tuple[str, float, float, Mapping[str, object], float, float, str, Optional[int]]] = []
    head_cov = None
    if head_inference in {"bayesian", "bayesian_mcmc"}:
        sample_lookup: Dict[str, Mapping[str, object]] = {}
        for sample in samples_list:
            site_id = sample.get("site_id")
            if site_id is None:
                continue
            if sample.get("lat") in (None, "") or sample.get("lon") in (None, ""):
                continue
            sample_lookup[str(site_id)] = sample

        if head_inference == "bayesian_mcmc":
            posterior = infer_heads_bayesian_mcmc(
                list(sample_lookup.values()),
                node_id_key="site_id",
                head_key=head_key,
                dtw_key=dtw_key,
                elevation_key=elevation_key,
                sigma_meas=sigma_meas,
                sigma_dtw=sigma_dtw,
                sigma_elev=sigma_elev,
                sigma_topo=sigma_topo,
                dtw_prior_mu=dtw_prior_mu,
                dtw_prior_sigma=dtw_prior_sigma,
                head_prior_mu=head_prior_mu,
                head_prior_sigma=head_prior_sigma,
                mcmc_draws=mcmc_draws,
                mcmc_chains=mcmc_chains,
                mcmc_target_accept=mcmc_target_accept,
                mcmc_warmup_fraction=mcmc_warmup_fraction,
            )
        else:
            posterior = infer_heads_bayesian_linear(
                list(sample_lookup.values()),
                node_id_key="site_id",
                head_key=head_key,
                dtw_key=dtw_key,
                elevation_key=elevation_key,
                sigma_meas=sigma_meas,
                sigma_dtw=sigma_dtw,
                sigma_elev=sigma_elev,
                sigma_topo=sigma_topo,
                dtw_prior_mu=dtw_prior_mu,
                dtw_prior_sigma=dtw_prior_sigma,
                head_prior_mu=head_prior_mu,
                head_prior_sigma=head_prior_sigma,
            )
        idx_map = posterior.index()
        head_cov = posterior.head_cov
        for node_id in posterior.node_ids:
            sample = sample_lookup.get(node_id)
            if sample is None:
                continue
            idx = idx_map[node_id]
            head_mean = float(posterior.head_mean[idx])
            head_sigma = float(math.sqrt(float(head_cov[idx, idx])) if head_cov is not None else sigma_topo)
            node_rows.append(
                (
                    node_id,
                    float(sample["lat"]),
                    float(sample["lon"]),
                    sample,
                    head_mean,
                    head_sigma,
                    posterior.tiers[idx],
                    idx,
                )
            )
    else:
        for sample in samples_list:
            site_id = sample.get("site_id")
            if site_id is None:
                continue
            if sample.get("lat") in (None, "") or sample.get("lon") in (None, ""):
                continue
            head, sigma, tier = infer_head_heuristic(
                sample,
                head_key=head_key,
                dtw_key=dtw_key,
                elevation_key=elevation_key,
                sigma_meas=sigma_meas,
                sigma_dtw=sigma_dtw,
                sigma_elev=sigma_elev,
                sigma_topo=sigma_topo,
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
                    None,
                )
            )

    edges: List[Edge] = []
    edge_ids = set()
    for node_id, lat, lon, sample, head, sigma, tier, node_idx in node_rows:
        candidates: List[Tuple[float, float, str, Mapping[str, object], float, float, str, List[str]]] = []
        for other_id, o_lat, o_lon, other_sample, other_head, other_sigma, other_tier, other_idx in node_rows:
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
            if head_cov is not None and node_idx is not None and other_idx is not None:
                var_delta = float(head_cov[node_idx, node_idx] + head_cov[other_idx, other_idx] - 2.0 * head_cov[node_idx, other_idx])
                sigma_delta = math.sqrt(max(var_delta, 0.0))
            else:
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
