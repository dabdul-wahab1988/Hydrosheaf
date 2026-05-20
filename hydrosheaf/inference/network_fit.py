"""Network-level fitting pipeline."""

import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import Config
from ..log import get_logger
from dataclasses import replace

logger = get_logger("inference.network_fit")


from ..data.schema import vector_from_sample
from ..graph.build import (
    EdgeInput,
    build_edges,
    infer_edges_from_coordinates,
    infer_edges_probabilistic,
)
from ..graph3d.build_3d import build_network_3d
from ..graph3d.types_3d import Edge3D
from ..graph.types import Edge
from ..models.reactions import build_reaction_dictionary
from ..phreeqc.constraints import build_edge_bounds
from ..phreeqc.runner import run_phreeqc
from ..models.ec_tds import predict_ec_tds
from ..models.redox import get_redox_constraints
from ..sheaf.topology_refine import refine_edges_with_sheaf
from .edge_fit import EdgeResult, fit_edge
from ..nitrate_source_v2 import infer_node_posteriors
import pandas as pd


def _safe_float(value: object) -> Optional[float]:
    if value is None:
        return None
    try:
        f = float(value)  # type: ignore
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None


def estimate_edge_residence_time_days(
    edge_attrs: Mapping[str, object], config: Config
) -> Tuple[Optional[float], Optional[float]]:
    """Estimate residence time from edge geometry and hydraulic gradient.

    Uses Darcy-based travel time approximation:
      tau = L * porosity / (K * i)
    where i is hydraulic gradient and L is path length.
    
    Returns:
        Tuple of (mean_tau_days, std_tau_days). Standard deviation is calculated
        via first-order Taylor expansion error propagation based on CV parameters.
    """
    # Explicit override (e.g., from MODPATH priors / external physics model).
    for key in (
        "edge_residence_time_days",
        "physics_tau_mean_days",
        "tau_days",
        "residence_time_days",
    ):
        explicit = _safe_float(edge_attrs.get(key))
        if explicit is not None and explicit > 0:
            std_explicit = _safe_float(edge_attrs.get("physics_tau_std_days"))
            return float(explicit), float(std_explicit) if std_explicit is not None else None

    k_m_per_day = float(getattr(config, "residence_time_hydraulic_k", 0.0) or 0.0)
    porosity = float(getattr(config, "residence_time_porosity", 0.0) or 0.0)
    if k_m_per_day <= 0 or not (0.0 < porosity <= 1.0):
        return None, None

    length_m = _safe_float(edge_attrs.get("distance_3d_m"))
    if length_m is None:
        dist_km = _safe_float(edge_attrs.get("distance_km"))
        if dist_km is not None:
            length_m = dist_km * 1000.0
    if length_m is None or length_m <= 0:
        return None, None

    delta_h = _safe_float(edge_attrs.get("delta_h"))
    gradient = None
    if delta_h is not None:
        gradient = abs(delta_h) / length_m
    else:
        h_grad = _safe_float(edge_attrs.get("horizontal_gradient"))
        v_grad = _safe_float(edge_attrs.get("vertical_gradient"))
        if h_grad is not None:
            gradient = abs(h_grad)
        elif v_grad is not None:
            gradient = abs(v_grad)

    if gradient is None or gradient <= 0:
        return None, None

    tau_days = (length_m * porosity) / (k_m_per_day * gradient)
    if tau_days <= 0 or not (tau_days < 1e12):
        return None, None
        
    # Error propagation for tau = L * n / (K * i)
    # CV_tau^2 = CV_n^2 + CV_K^2 + CV_i^2 (assuming L is exact enough)
    cv_k = float(getattr(config, "residence_time_hydraulic_k_cv", 1.0))
    cv_n = float(getattr(config, "residence_time_porosity_cv", 0.2))
    cv_i = float(getattr(config, "residence_time_gradient_cv", 0.5))
    
    cv_tau = (cv_k**2 + cv_n**2 + cv_i**2)**0.5
    std_tau = tau_days * cv_tau

    return float(tau_days), float(std_tau)


def _sample_map(samples: object) -> Dict[str, Mapping[str, object]]:
    if isinstance(samples, Mapping):
        # We need to ensure the values are Mapping[str, object]
        return {str(k): v for k, v in samples.items()}
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, object]] = {}
        for row in samples:
            site_id = row.get("site_id")
            if site_id is None:
                site_id = row.get("SampleID") or row.get("Code")
            if site_id is None:
                raise ValueError("Each sample row must include site_id, SampleID, or Code.")
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def fit_network(
    samples: object,
    edges: Iterable[EdgeInput],
    config: Config,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
    residence_time_overrides: Optional[Mapping[str, float]] = None,
    lateral_neighbors: Optional[Mapping[str, List[str]]] = None,
    pre_si_mask: Optional[Mapping[str, float]] = None,
) -> List[EdgeResult]:
    """Fit inverse geochemical reactions over a directed flow network."""

    sample_map = _sample_map(samples)
    built_edges = build_edges(edges)
    lateral_neighbors = lateral_neighbors or {}

    # 1. Run PHREEQC for all samples to get thermodynamic constraints
    if config.phreeqc_enabled and phreeqc_results is None:
        phreeqc_results = run_phreeqc(sample_map.values(), config)

    # Prepare reaction dictionary once to know active labels
    _, reaction_labels, mineral_mask, _ = build_reaction_dictionary(config)

    results: List[EdgeResult] = []
    for edge in built_edges:
        if edge.u not in sample_map or edge.v not in sample_map:
            continue

        sample_u = sample_map[edge.u]
        sample_v = sample_map[edge.v]

        x_u, sample_u_norm = vector_from_sample(
            sample_u,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if x_u is None:
            continue

        x_v, sample_v_norm = vector_from_sample(
            sample_v,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if x_v is None:
            continue

        # Build edge-specific configuration (handling layer-based priors)
        layer_key = getattr(config, "layer_key", "aquifer_layer")
        layer_value = sample_v.get(layer_key)
        layer_idx = None
        try:
            if layer_value not in (None, ""):
                layer_idx = int(layer_value)
        except (TypeError, ValueError):
            layer_idx = None

        config_edge = config
        if layer_idx is not None:
            minerals_for_layer = getattr(config, "layer_mineral_map", {}).get(layer_idx)
            if minerals_for_layer:
                config_edge = replace(config, active_minerals=list(minerals_for_layer))

        # Build thermodynamic bounds for this edge
        edge_bounds: Dict[str, object] = {}
        if config_edge.phreeqc_enabled and phreeqc_results is not None:
            edge_bounds = build_edge_bounds(
                phreeqc_results,
                [edge],
                reaction_labels,
                mineral_mask,
                config_edge,
            ).get(edge.edge_id, {})

        # Ensure bounds structure exists for redox overrides
        if "lb" not in edge_bounds or not isinstance(edge_bounds["lb"], list):
            edge_bounds["lb"] = [None] * len(reaction_labels)
        if "ub" not in edge_bounds or not isinstance(edge_bounds["ub"], list):
            edge_bounds["ub"] = [None] * len(reaction_labels)
        if "constraints_active" not in edge_bounds or not isinstance(
            edge_bounds["constraints_active"], dict
        ):
            edge_bounds["constraints_active"] = {}

        tau_edge = None
        std_tau_edge = None
        if residence_time_overrides is not None:
            tau_override = residence_time_overrides.get(edge.edge_id)
            try:
                tau_edge = float(tau_override) if tau_override is not None else None
            except (TypeError, ValueError):
                tau_edge = None
        if tau_edge is None:
            tau_edge, std_tau_edge = estimate_edge_residence_time_days(edge.attrs or {}, config_edge)

        # Apply redox overrides
        redox_overrides = get_redox_constraints(sample_v_norm, reaction_labels)
        if redox_overrides:
            lb_list = edge_bounds["lb"]
            ub_list = edge_bounds["ub"]
            assert isinstance(lb_list, list)
            assert isinstance(ub_list, list)

            for i, label in enumerate(reaction_labels):
                if label in redox_overrides:
                    lb_val, ub_val = redox_overrides[label]
                    lb_list[i] = lb_val
                    ub_list[i] = ub_val

            ca_dict = edge_bounds["constraints_active"]
            assert isinstance(ca_dict, dict)
            ca_dict["redox"] = "active"

        # Prepare extra endmembers from lateral neighbors
        extra_endmembers: Dict[str, List[float]] = {}
        neighbors_of_u = lateral_neighbors.get(edge.u, [])
        for neighbor_id in neighbors_of_u:
            if neighbor_id not in sample_map:
                continue
            n_sample = sample_map[neighbor_id]
            x_n, _ = vector_from_sample(
                 n_sample, 
                 config.ion_order, 
                 config.missing_policy, 
                 config.detection_limit_policy
            )
            if x_n is not None:
                extra_endmembers[f"lateral_{neighbor_id}"] = x_n

        result = fit_edge(
            x_u,
            x_v,
            config_edge,
            edge_id=edge.edge_id,
            u=edge.u,
            v=edge.v,
            obs_v=sample_v_norm,
            bounds=edge_bounds,
            obs_u=sample_u_norm,
            residence_time_days=tau_edge,
            residence_time_std_days=std_tau_edge,
            extra_endmembers=extra_endmembers,
            pre_si_mask=pre_si_mask,
        )
        
        if result.transport_model != "mix":
             logger.debug(f"Edge {edge.edge_id}: Selected model '{result.transport_model}' (obj={result.objective_score:.4f})")

        result.edge_residence_time_days = tau_edge
        result.physics_tau_mean_days = tau_edge
        result.physics_tau_std_days = std_tau_edge

        # Copy edge-level spatial attributes
        edge_attrs = edge.attrs or {}

        def _get_float(key: str) -> Optional[float]:
            v = edge_attrs.get(key)
            if v is None:
                return None
            try:
                return float(v)
            except (ValueError, TypeError):
                return None

        def _get_str(key: str) -> Optional[str]:
            v = edge_attrs.get(key)
            if v is None:
                return None
            return str(v)

        result.edge_sheaf_score_local = _get_float("sheaf_score_local")
        result.edge_sheaf_score_global = _get_float("sheaf_score_global")
        result.edge_sheaf_residual = _get_float("sheaf_residual_global")
        result.edge_sheaf_pi_evap = _get_float("sheaf_pi_evap")
        result.edge_sheaf_cost_iso = _get_float("sheaf_cost_iso")
        result.edge_sheaf_cost_cl = _get_float("sheaf_cost_cl")
        result.edge_sheaf_flags = _get_str("sheaf_flags")
        result.edge_distance_km = _get_float("distance_km")
        result.edge_delta_h = _get_float("delta_h")
        result.edge_sigma_delta_h = _get_float("sigma_delta_h")
        result.edge_source_tier = _get_str("source_tier")
        result.edge_flags = _get_str("flags")
        result.edge_distance_3d_m = _get_float("distance_3d_m")
        result.edge_horizontal_distance_m = _get_float("horizontal_distance_m")
        result.edge_vertical_distance_m = _get_float("vertical_distance_m")
        result.edge_type_3d = _get_str("edge_type_3d")
        result.edge_prob_head = _get_float("prob_head")
        result.edge_prob_distance = _get_float("prob_distance")
        result.edge_prob_layer = _get_float("prob_layer")
        result.edge_horizontal_gradient = _get_float("horizontal_gradient")
        result.edge_vertical_gradient = _get_float("vertical_gradient")
        result.physics_source = _get_str("physics_source")
        results.append(result)

    # Nitrate Source Discrimination (v2) Integration
    if config.nitrate_source_enabled:
        df = pd.DataFrame(list(sample_map.values()))
        if "site_id" not in df.columns:
            df["site_id"] = df.index
        df.set_index("site_id", drop=False, inplace=True)
        all_posteriors = infer_node_posteriors(df, results, config=config)
        for result in results:
            node_v = result.v
            post = all_posteriors.get(node_v)
            if post:
                result.nitrate_source_p_manure = post.p_manure
                result.nitrate_source_logit = post.logit_manure
                result.nitrate_source_evidence = post.evidence
                result.nitrate_source_gates = post.gates

    return results


def fit_edges(
    samples: object,
    edges: Iterable[object],
    config: Config,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
) -> List[EdgeResult]:
    return fit_network(samples, edges, config, phreeqc_results=phreeqc_results)


def infer_edges(
    samples: object,
    max_neighbors: int = 1,
    allow_uphill: bool = False,
    head_key: str = "hydraulic_head",
    elevation_key: str = "elevation",
    method: str = "probabilistic",
    config: Optional[Config] = None,
    edge_attr_overrides: Optional[Mapping[str, Mapping[str, object]]] = None,
    layer_definition: Optional[Dict[str, object]] = None,
) -> List[Edge]:
    if isinstance(samples, Mapping):
        samples_iter = list(samples.values())
    elif isinstance(samples, Sequence):
        samples_iter = list(samples)
    else:
        raise TypeError("Unsupported samples input type.")
    if method == "simple":
        return infer_edges_from_coordinates(
            samples_iter,
            max_neighbors=max_neighbors,
            allow_uphill=allow_uphill,
            head_key=head_key,
            elevation_key=elevation_key,
        )
    config = config or Config()

    def _apply_overrides(edges: List[Edge]) -> List[Edge]:
        if not edge_attr_overrides:
            return edges
        for edge in edges:
            override = edge_attr_overrides.get(edge.edge_id)
            if not override:
                continue
            attrs = dict(edge.attrs or {})
            attrs.update(dict(override))
            edge.attrs = attrs
        return edges

    def infer_probabilistic_edges_from_config(config_for_edges: Config) -> List[Edge]:
        if getattr(config_for_edges, "network_3d_enabled", False):
            network = build_network_3d(
                list(samples_iter),
                config_for_edges,
                layer_definition=layer_definition,
                use_haversine=True,
            )
            edges_3d: List[Edge3D] = network.edges
            converted: List[Edge] = []
            for edge3d in edges_3d:
                attrs = {
                    "distance_km": edge3d.distance_3d / 1000.0,
                    "horizontal_distance_m": edge3d.horizontal_distance_m,
                    "vertical_distance_m": edge3d.vertical_distance_m,
                    "distance_3d_m": edge3d.distance_3d,
                    "edge_type_3d": edge3d.edge_type,
                    "horizontal_gradient": edge3d.horizontal_gradient,
                    "vertical_gradient": edge3d.vertical_gradient,
                    "prob_head": edge3d.prob_head,
                    "prob_distance": edge3d.prob_distance,
                    "prob_layer": edge3d.prob_layer,
                    "p_uv": edge3d.prob_combined,
                    "edge_confidence": edge3d.prob_combined,
                }
                converted.append(
                    Edge(edge_id=edge3d.edge_id, u=edge3d.u, v=edge3d.v, attrs=attrs)
                )
            return _apply_overrides(converted)
        inferred = infer_edges_probabilistic(
            samples_iter,
            radius_km=config_for_edges.edge_radius_km,
            max_neighbors=config_for_edges.edge_max_neighbors,
            p_min=config_for_edges.edge_p_min,
            head_key=config_for_edges.edge_head_key,
            dtw_key=config_for_edges.edge_dtw_key,
            elevation_key=config_for_edges.edge_elevation_key,
            aquifer_key=config_for_edges.edge_aquifer_key,
            screen_depth_key=config_for_edges.edge_screen_depth_key,
            well_depth_key=config_for_edges.edge_well_depth_key,
            sigma_meas=config_for_edges.edge_sigma_meas,
            sigma_dtw=config_for_edges.edge_sigma_dtw,
            sigma_elev=config_for_edges.edge_sigma_elev,
            sigma_topo=config_for_edges.edge_sigma_topo,
            gradient_min=config_for_edges.edge_gradient_min,
            depth_mismatch=config_for_edges.edge_depth_mismatch,
            head_inference=getattr(config_for_edges, "edge_head_inference", "heuristic"),
        )
        return _apply_overrides(inferred)

    if method == "probabilistic_map":
        config.validate()
        prior_weight = float(getattr(config, "edge_map_prior_weight", 0.0) or 0.0)
        candidate_multiplier = int(getattr(config, "edge_map_candidate_multiplier", 5) or 5)
        candidate_config = replace(
            config,
            edge_p_min=float(getattr(config, "edge_map_p_min", 0.1)),
            edge_max_neighbors=int(config.edge_max_neighbors) * candidate_multiplier,
        )
        candidate_edges = infer_probabilistic_edges_from_config(candidate_config)
        if not candidate_edges:
            return []
        scoring_config = replace(config, uncertainty_method="none", nitrate_source_enabled=False)
        candidate_results = fit_network(samples_iter, candidate_edges, scoring_config)
        result_by_edge: Dict[str, float] = {res.edge_id: float(res.objective_score) for res in candidate_results}
        scored_by_u: Dict[str, List[tuple]] = {}
        for edge in candidate_edges:
            if edge.edge_id not in result_by_edge: continue
            chem_score = result_by_edge[edge.edge_id]
            attrs = edge.attrs or {}
            p_uv = float(attrs.get("edge_confidence", attrs.get("p_uv", 1.0)))
            map_score = chem_score + prior_weight * (-math.log(max(1e-12, p_uv)))
            scored_by_u.setdefault(edge.u, []).append((map_score, edge))
        selected = []
        for _, scored in scored_by_u.items():
            scored.sort(key=lambda item: item[0])
            for _, edge in scored[: int(config.edge_max_neighbors)]:
                selected.append(edge)
        return selected

    if method == "probabilistic_sheaf":
        candidate_multiplier = int(getattr(config, "edge_map_candidate_multiplier", 5) or 5)
        candidate_config = replace(config, edge_max_neighbors=int(config.edge_max_neighbors) * candidate_multiplier)
        candidate_edges = infer_probabilistic_edges_from_config(candidate_config)
        if not candidate_edges: return []
        return _apply_overrides(refine_edges_with_sheaf(samples_iter, candidate_edges, config))

    return infer_probabilistic_edges_from_config(config)


def summarize_network(results: List[EdgeResult]) -> Dict[str, Any]:
    """Calculate summary statistics for a network of edge results."""
    if not results:
        return {}

    n = len(results)
    avg_objective = sum(r.objective_score for r in results) / n
    avg_anomaly = sum(r.anomaly_norm for r in results) / n

    # Model counts
    models: Dict[str, int] = {}
    for r in results:
        models[r.transport_model] = models.get(r.transport_model, 0) + 1

    # Reaction statistics
    labels = results[0].z_labels
    m = len(labels)
    counts = [0] * m
    sums = [0.0] * m

    for r in results:
        for i, val in enumerate(r.z_extents):
            if abs(val) > 1e-6:
                counts[i] += 1
                sums[i] += val

    summary = {
        "edge_count": n,
        "avg_objective": avg_objective,
        "avg_anomaly_norm": avg_anomaly,
        "transport_counts": models,
        "reaction_labels": labels,
        "reaction_nonzero": counts,
        "reaction_means": [s / max(1, c) for s, c in zip(sums, counts)],
    }

    return summary


def fit_edges(
    samples: object,
    edges: Iterable[object],
    config: Config,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
) -> List[EdgeResult]:
    """Legacy helper for fitting edges."""
    return fit_network(samples, edges, config, phreeqc_results=phreeqc_results)





def edge_process_maps(results: List[EdgeResult]) -> Dict[str, List[Dict[str, object]]]:
    transport_rows: List[Dict[str, object]] = []
    reaction_rows: List[Dict[str, object]] = []
    for result in results:
        transport_row = {"edge_id": result.edge_id, **result.transport_probabilities}
        reaction_row = {
            "edge_id": result.edge_id,
            **{
                label: abs(value)
                for label, value in zip(result.z_labels, result.z_extents)
            },
        }
        transport_rows.append(transport_row)
        reaction_rows.append(reaction_row)
    return {
        "transport_likelihoods": transport_rows,
        "reaction_intensity": reaction_rows,
    }


def predict_node_ec_tds(samples: object, config: Config) -> List[Dict[str, object]]:
    sample_map = _sample_map(samples)
    rows: List[Dict[str, object]] = []
    for site_id, sample in sample_map.items():
        values, sample_norm = vector_from_sample(
            sample,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if values is None:
            continue
        ec_pred, tds_pred = predict_ec_tds(values, config)
        rows.append(
            {
                "site_id": site_id,
                "sample_id": sample_norm.get("sample_id"),
                "ec_pred": ec_pred,
                "tds_pred": tds_pred,
            }
        )
    return rows
