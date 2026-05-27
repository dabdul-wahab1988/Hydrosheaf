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
    a = (
        math.sin(d_phi / 2.0) ** 2
        + math.cos(phi1) * math.cos(phi2) * math.sin(d_lambda / 2.0) ** 2
    )
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
    node_rows: List[
        Tuple[str, float, float, Mapping[str, object], float, float, str, Optional[int]]
    ] = []
    head_cov = None
    if head_inference in {"bayesian", "bayesian_mcmc"}:
        sample_lookup: Dict[str, Mapping[str, object]] = {}
        for sample in samples_list:
            site_id = sample.get("site_id")
            if site_id is None:
                continue
            
            # Safe float conversion for lat/lon
            lat_v = sample.get("lat")
            lon_v = sample.get("lon")
            if lat_v is None or lon_v is None:
                continue
                
            try:
                lat_val = float(lat_v) # type: ignore
                lon_val = float(lon_v) # type: ignore
            except (ValueError, TypeError):
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
            head_sigma = float(
                math.sqrt(float(head_cov[idx, idx]))
                if head_cov is not None
                else sigma_topo
            )
            
            # lat/lon already validated above
            lat_v = sample.get("lat")
            lon_v = sample.get("lon")
            if lat_v is None or lon_v is None:
                continue
            
            # Runtime validation
            if not isinstance(lat_v, (int, float, str)) or not isinstance(lon_v, (int, float, str)):
                continue

            try:
                lat_val = float(f"{lat_v}")
                lon_val = float(f"{lon_v}")
            except (ValueError, TypeError):
                continue
            
            node_rows.append(
                (
                    node_id,
                    lat_val,
                    lon_val,
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
            
            lat_v = sample.get("lat")
            lon_v = sample.get("lon")
            if lat_v is None or lon_v is None:
                continue
                
            try:
                lat_val = float(lat_v) # type: ignore
                lon_val = float(lon_v) # type: ignore
            except (ValueError, TypeError):
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
                    lat_val,
                    lon_val,
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
        candidates: List[
            Tuple[float, float, str, Mapping[str, object], float, float, str, List[str], bool]
        ] = []
        for (
            other_id,
            o_lat,
            o_lon,
            other_sample,
            other_head,
            other_sigma,
            other_tier,
            other_idx,
        ) in node_rows:
            if other_id == node_id:
                continue
            if (
                aquifer_key
                and sample.get(aquifer_key)
                and other_sample.get(aquifer_key)
            ):
                if sample.get(aquifer_key) != other_sample.get(aquifer_key):
                    continue
            distance = _haversine_km(lat, lon, o_lat, o_lon)
            if radius_km and distance > radius_km:
                continue
            if distance == 0:
                continue

            delta_h = head - other_head
            if head_cov is not None and node_idx is not None and other_idx is not None:
                var_delta = float(
                    head_cov[node_idx, node_idx]
                    + head_cov[other_idx, other_idx]
                    - 2.0 * head_cov[node_idx, other_idx]
                )
                sigma_delta = math.sqrt(max(var_delta, 0.0))
            else:
                sigma_delta = math.sqrt(sigma**2 + other_sigma**2)
            if sigma_delta == 0:
                p_uv = 0.5 if delta_h == 0 else (1.0 if delta_h > 0 else 0.0)
            else:
                p_uv = _normal_cdf(delta_h / sigma_delta)

            gradient = abs(delta_h) / (distance * 1000.0)
            flags: List[str] = []
            is_lateral = False
            
            # Check for lateral dispersive edges (low gradient but spatially close)
            # These are preserved as type="lateral" even if p_uv is low, 
            # to serve as mixing candidates for dispersion.
            if gradient_min > 0 and gradient < gradient_min:
                if radius_km > 0 and distance < radius_km:
                    # Mark as lateral candidate
                    is_lateral = True
                
                flags.append("flat_gradient")
                factor = gradient / gradient_min
                p_uv = 0.5 + (p_uv - 0.5) * factor

            depth_a = sample.get(screen_depth_key) or sample.get(well_depth_key)
            depth_b = other_sample.get(screen_depth_key) or other_sample.get(
                well_depth_key
            )
            if depth_a is not None and depth_b is not None:
                try:
                    da_val = float(depth_a)
                    db_val = float(depth_b)
                    if abs(da_val - db_val) > depth_mismatch:
                        flags.append("depth_mismatch")
                except (ValueError, TypeError):
                    pass

            if p_uv < p_min and not is_lateral:
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
                    is_lateral,
                )
            )

        candidates.sort(key=lambda item: (-item[0], item[1]))
        if max_neighbors > 0:
            # We want to keep lateral edges regardless of max_neighbors limit for primary flow
            # But simpler to just include them in the pool for now. 
            # Ideally we separate them.
            pass
            # candidates = candidates[:max_neighbors] 
            # NOTE: Lateral edges should not displace primary edges.
            # Logic below handles it.

        primary_count = 0
        
        for (
            p_uv,
            distance,
            other_id,
            other_sample,
            delta_h,
            sigma_delta,
            tier_pair,
            flags,
            is_lateral,
        ) in candidates:
            # Respect max_neighbors only for primary edges
            if not is_lateral:
                if max_neighbors > 0 and primary_count >= max_neighbors:
                    continue
                primary_count += 1
            
            edge_id = f"{node_id}->{other_id}"
            if edge_id in edge_ids:
                continue
            
            edge_type = "lateral" if is_lateral else "primary"
            
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
                    type=edge_type
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
    def _edge_attrs(
        distance_km: float,
        source_head: Optional[float],
        target_head: Optional[float],
        rank_index: int,
        rank_total: int,
        explicit_flow: bool = False,
    ) -> Dict[str, object]:
        attrs: Dict[str, object] = {
            "distance_km": distance_km,
            "head_source": head_key if source_head is not None else elevation_key,
            "allow_uphill": allow_uphill,
        }
        if source_head is not None and target_head is not None:
            delta_h = source_head - target_head
            attrs["delta_h"] = delta_h
            if distance_km > 0:
                attrs["head_gradient"] = delta_h / distance_km

        if explicit_flow:
            attrs["p_uv"] = 1.0
            attrs["edge_confidence"] = 1.0
        elif rank_total > 0:
            prob = max(1e-6, (rank_total - rank_index) / float(rank_total + 1))
            attrs["p_uv"] = prob
            attrs["edge_confidence"] = prob
        return attrs

    samples_list = list(samples)
    sample_map = {
        str(sample.get("site_id")): sample
        for sample in samples_list
        if sample.get("site_id")
    }

    nodes: List[Tuple[str, float, float, Optional[float]]] = []
    for sample in samples_list:
        site_id = sample.get("site_id")
        if site_id is None:
            continue
        
        lat_val = sample.get("lat")
        lon_val = sample.get("lon")
        
        if lat_val is None or lon_val is None:
            continue
            
        try:
            lat_val = float(lat_val) # type: ignore
            lon_val = float(lon_val) # type: ignore
        except (ValueError, TypeError):
            continue
            
        elevation = sample.get(head_key)
        if elevation is None:
            elevation = sample.get(elevation_key)
            
        elev_val: Optional[float] = None
        if elevation is not None:
            try:
                elev_val = float(elevation) # type: ignore
            except (ValueError, TypeError):
                elev_val = None
                
        nodes.append(
            (str(site_id), lat_val, lon_val, elev_val)
        )

    node_lookup = {
        node_id: (lat, lon, elevation) for node_id, lat, lon, elevation in nodes
    }
    edges: List[Edge] = []
    edge_keys = set()

    for node_id, lat, lon, elevation in nodes:
        sample = sample_map.get(node_id)
        if sample and flow_to_key in sample:
            target = str(sample[flow_to_key])
            if target in node_lookup:
                edge_id = f"{node_id}->{target}"
                if edge_id not in edge_keys:
                    o_lat, o_lon, o_elev = node_lookup[target]
                    distance = _haversine_km(lat, lon, o_lat, o_lon)
                    edges.append(
                        Edge(
                            edge_id=edge_id,
                            u=node_id,
                            v=target,
                            attrs=_edge_attrs(
                                distance,
                                elevation,
                                o_elev,
                                rank_index=0,
                                rank_total=1,
                                explicit_flow=True,
                            ),
                        )
                    )
                    edge_keys.add(edge_id)
                continue

        candidates: List[Tuple[float, str, Optional[float]]] = []
        for other_id, o_lat, o_lon, o_elev in nodes:
            if other_id == node_id:
                continue
            if elevation is not None and o_elev is not None and not allow_uphill:
                if o_elev >= elevation:
                    continue
            distance = _haversine_km(lat, lon, o_lat, o_lon)
            candidates.append((distance, other_id, o_elev))
        candidates.sort(key=lambda item: item[0])
        selected_candidates = candidates[:max_neighbors]
        rank_total = len(selected_candidates)
        for rank_index, (distance, target, target_head) in enumerate(selected_candidates):
            edge_id = f"{node_id}->{target}"
            if edge_id in edge_keys:
                continue
            edges.append(
                Edge(
                    edge_id=edge_id,
                    u=node_id,
                    v=target,
                    attrs=_edge_attrs(
                        distance,
                        elevation,
                        target_head,
                        rank_index=rank_index,
                        rank_total=rank_total,
                    ),
                )
            )
            edge_keys.add(edge_id)

    return edges
