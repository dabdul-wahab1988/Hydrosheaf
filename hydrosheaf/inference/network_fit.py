"""Network-level fitting pipeline."""

from typing import Dict, Iterable, List, Mapping, Optional, Sequence

from ..config import Config
from ..data.schema import vector_from_sample
from ..graph.build import build_edges, infer_edges_from_coordinates, infer_edges_probabilistic
from ..graph.types import Edge
from ..models.reactions import build_reaction_dictionary
from ..phreeqc.constraints import build_edge_bounds
from ..phreeqc.runner import run_phreeqc
from ..models.ec_tds import predict_ec_tds
from ..models.redox import get_redox_constraints
from .edge_fit import EdgeResult, fit_edge
from ..nitrate_source_v2 import infer_node_posteriors
import pandas as pd


def _sample_map(samples: object) -> Dict[str, Mapping[str, float]]:
    if isinstance(samples, Mapping):
        return samples  # Already mapping site_id -> record
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, float]] = {}
        for row in samples:
            site_id = row.get("site_id")
            if site_id is None:
                raise ValueError("Each sample row must include site_id.")
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def fit_network(
    samples: object,
    edges: Iterable[object],
    config: Config,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
) -> List[EdgeResult]:
    sample_map = _sample_map(samples)
    built_edges = build_edges(edges)

    bounds_by_edge: Dict[str, Dict[str, object]] = {}
    if config.phreeqc_enabled:
        if phreeqc_results is None:
            phreeqc_results = run_phreeqc(sample_map.values(), config)
        _, labels, mineral_mask = build_reaction_dictionary(config)
        bounds_by_edge = build_edge_bounds(phreeqc_results, built_edges, labels, mineral_mask, config)

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

        edge_bounds = bounds_by_edge.get(edge.edge_id, {})
        # Apply redox overrides
        _, labels, _ = build_reaction_dictionary(config)
        redox_overrides = get_redox_constraints(sample_v_norm, labels)
        if redox_overrides:
            # Merge with existing bounds
            if "lb" not in edge_bounds: edge_bounds["lb"] = [None] * len(labels)
            if "ub" not in edge_bounds: edge_bounds["ub"] = [None] * len(labels)
            
            for i, label in enumerate(labels):
                if label in redox_overrides:
                    l, u = redox_overrides[label]
                    edge_bounds["lb"][i] = l
                    edge_bounds["ub"][i] = u
            
            edge_bounds.setdefault("constraints_active", {})["redox"] = "active"

        result = fit_edge(
            x_u,
            x_v,
            config,
            edge_id=edge.edge_id,
            u=edge.u,
            v=edge.v,
            obs_v=sample_v_norm,
            bounds=edge_bounds,
            obs_u=sample_u_norm,
        )
        edge_attrs = edge.attrs or {}
        result.edge_confidence = edge_attrs.get("edge_confidence", edge_attrs.get("p_uv"))
        result.edge_distance_km = edge_attrs.get("distance_km")
        result.edge_delta_h = edge_attrs.get("delta_h")
        result.edge_sigma_delta_h = edge_attrs.get("sigma_delta_h")
        result.edge_source_tier = edge_attrs.get("source_tier")
        result.edge_flags = edge_attrs.get("flags")
        results.append(result)

    # Nitrate Source Discrimination (v2) Integration
    if config.nitrate_source_enabled:
        # 1. Convert samples to DataFrame for bulk stats
        df = pd.DataFrame(list(sample_map.values()))
        if "site_id" not in df.columns:
            df["site_id"] = df.index  # Fallback if needed
        df.set_index("site_id", drop=False, inplace=True)
        
        # 2. Run Inference with Overrides
        overrides = {
            "weights": config.nitrate_source_weights,
            "prior_prob": config.nitrate_source_prior,
            "evap_gate_factor": config.nitrate_source_evap_gate,
            "nitrate_source_d_excess_p25": config.nitrate_source_d_excess_p25,
            "nitrate_source_po4_p90": config.nitrate_source_po4_p90,
        }
        node_posteriors = infer_node_posteriors(df, results, overrides)

        # 3. Attach Node Posteriors to Edges (for output)
        for res in results:
            if res.v in node_posteriors:
                nr = node_posteriors[res.v]
                res.nitrate_source_p_manure = nr.p_manure
                res.nitrate_source_logit = nr.logit_score
                res.nitrate_source_evidence = nr.top_evidence
                res.nitrate_source_gates = nr.gating_flags

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
    return infer_edges_probabilistic(
        samples_iter,
        radius_km=config.edge_radius_km,
        max_neighbors=config.edge_max_neighbors,
        p_min=config.edge_p_min,
        head_key=config.edge_head_key,
        dtw_key=config.edge_dtw_key,
        elevation_key=config.edge_elevation_key,
        aquifer_key=config.edge_aquifer_key,
        screen_depth_key=config.edge_screen_depth_key,
        well_depth_key=config.edge_well_depth_key,
        sigma_meas=config.edge_sigma_meas,
        sigma_dtw=config.edge_sigma_dtw,
        sigma_elev=config.edge_sigma_elev,
        sigma_topo=config.edge_sigma_topo,
        gradient_min=config.edge_gradient_min,
        depth_mismatch=config.edge_depth_mismatch,
    )


def summarize_network(results: List[EdgeResult]) -> Dict[str, object]:
    if not results:
        return {
            "edge_count": 0,
            "transport_counts": {},
            "transport_probabilities_mean": {},
            "reaction_means": {},
            "reaction_intensity_mean": {},
            "reaction_nonzero": {},
        }

    transport_counts: Dict[str, int] = {}
    transport_prob_sums: Dict[str, float] = {}
    reaction_sums: Dict[str, float] = {}
    reaction_counts: Dict[str, int] = {}
    reaction_nonzero: Dict[str, int] = {}
    reaction_intensity_sums: Dict[str, float] = {}

    for result in results:
        transport_counts[result.transport_model] = transport_counts.get(result.transport_model, 0) + 1
        for key, value in result.transport_probabilities.items():
            transport_prob_sums[key] = transport_prob_sums.get(key, 0.0) + value
        for label, value in zip(result.z_labels, result.z_extents):
            reaction_sums[label] = reaction_sums.get(label, 0.0) + value
            reaction_counts[label] = reaction_counts.get(label, 0) + 1
            if abs(value) > 1e-9:
                reaction_nonzero[label] = reaction_nonzero.get(label, 0) + 1
            reaction_intensity_sums[label] = reaction_intensity_sums.get(label, 0.0) + abs(value)

    reaction_means = {
        label: reaction_sums[label] / reaction_counts[label]
        for label in reaction_sums
    }
    transport_probabilities_mean = {
        key: transport_prob_sums[key] / len(results)
        for key in transport_prob_sums
    }
    reaction_intensity_mean = {
        label: reaction_intensity_sums[label] / reaction_counts[label]
        for label in reaction_intensity_sums
    }

    return {
        "edge_count": len(results),
        "transport_counts": transport_counts,
        "transport_probabilities_mean": transport_probabilities_mean,
        "reaction_means": reaction_means,
        "reaction_intensity_mean": reaction_intensity_mean,
        "reaction_nonzero": reaction_nonzero,
    }


def edge_process_maps(results: List[EdgeResult]) -> Dict[str, List[Dict[str, object]]]:
    transport_rows: List[Dict[str, object]] = []
    reaction_rows: List[Dict[str, object]] = []
    for result in results:
        transport_row = {"edge_id": result.edge_id, **result.transport_probabilities}
        reaction_row = {
            "edge_id": result.edge_id,
            **{label: abs(value) for label, value in zip(result.z_labels, result.z_extents)},
        }
        transport_rows.append(transport_row)
        reaction_rows.append(reaction_row)
    return {"transport_likelihoods": transport_rows, "reaction_intensity": reaction_rows}


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
