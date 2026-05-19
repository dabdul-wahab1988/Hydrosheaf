"""High-level API helpers for optional Hydrosheaf modules."""

from __future__ import annotations

from datetime import datetime
from statistics import median
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .config import Config
from .data.parsing import (
    extract_sample_decimal_year as _extract_sample_decimal_year,
    sample_list as _sample_list,
)
from .data.schema import parse_numeric
from .data.validation import (
    auto_disable_missing_modules,
    validate_required_inputs,
)
from .graph.build import EdgeInput, build_edges
from .graph.types import Edge
from .inference.edge_fit import EdgeResult
from .inference.network_fit import fit_network
from .log import get_logger
from .physics.priors import PhysicsPrior, apply_physics_priors
from .temporal import TemporalEdgeResult, TemporalNode
from .temporal.temporal_edge_fit import fit_temporal_edge
from .nuclear.nuclides import get_nuclide
from .models.latent import identify_latent_endmembers

from .vadose.contracts import (
    VadoseForcingSample,
    VadoseLinksRow,
    VadoseProfile,
    VadoseRunConfig,
)
from .vadose.run import build_vadose_edge_priors

logger = get_logger("api")


def infer_network_ages_bayesian(*args: Any, **kwargs: Any) -> Any:
    """Lazy proxy to avoid importing nuclear MCMC stack at module import time."""
    from .nuclear.network_aging import infer_network_ages_bayesian as _infer

    return _infer(*args, **kwargs)





def build_vadose_priors(
    profile: VadoseProfile,
    forcing: Sequence[VadoseForcingSample],
    links: Sequence[VadoseLinksRow],
    *,
    config: Optional[VadoseRunConfig] = None,
    water_table_depth_m: Optional[Sequence[Tuple[datetime, float]]] = None,
) -> Tuple[List[PhysicsPrior], Dict[str, Dict[str, object]]]:
    """Run vadose physics and return PhysicsPrior objects plus diagnostics."""
    vadose_priors, details = build_vadose_edge_priors(
        profile,
        forcing,
        links,
        config=config,
        water_table_depth_m=water_table_depth_m,
    )
    physics_priors = [prior.to_physics_prior() for prior in vadose_priors]
    return physics_priors, details


def fit_temporal_edges(
    temporal_nodes: Mapping[str, TemporalNode],
    edges: Iterable[EdgeInput],
    config: Config,
    *,
    hydraulic_params_by_edge: Optional[Mapping[str, Dict[str, float]]] = None,
) -> Tuple[List[TemporalEdgeResult], Dict[str, float]]:
    """Fit temporal edges and return residence time overrides."""
    built_edges = build_edges(edges)
    temporal_results: List[TemporalEdgeResult] = []
    residence_time_overrides: Dict[str, float] = {}
    for edge in built_edges:
        if edge.u not in temporal_nodes or edge.v not in temporal_nodes:
            continue
        hydraulic_params: Optional[Dict[str, float]] = None
        if hydraulic_params_by_edge is not None:
            hydraulic_params = hydraulic_params_by_edge.get(edge.edge_id)
        temporal_result = fit_temporal_edge(
            temporal_nodes[edge.u],
            temporal_nodes[edge.v],
            config,
            edge_id=edge.edge_id,
            hydraulic_params=hydraulic_params,
        )
        temporal_results.append(temporal_result)
        residence_time_overrides[edge.edge_id] = temporal_result.residence_time_days
    return temporal_results, residence_time_overrides


def attach_temporal_results(
    edge_results: List[EdgeResult],
    temporal_results: Sequence[TemporalEdgeResult],
) -> None:
    """Attach temporal summaries to EdgeResult objects in-place."""
    temporal_by_edge = {res.edge_id: res for res in temporal_results}
    for res in edge_results:
        temporal = temporal_by_edge.get(res.edge_id)
        if temporal is None:
            continue
        res.temporal_residence_time_days = temporal.residence_time_days
        res.temporal_residence_time_method = temporal.residence_time_method
        res.temporal_residence_time_uncertainty = temporal.residence_time_uncertainty
        res.temporal_transport_model = temporal.transport_model
        res.temporal_gamma_mean = temporal.gamma_mean
        res.temporal_gamma_std = temporal.gamma_std
        res.temporal_f_mean = temporal.f_mean
        res.temporal_f_std = temporal.f_std
        res.temporal_reaction_extents_mean = list(temporal.reaction_extents_mean or [])
        res.temporal_reaction_extents_std = list(temporal.reaction_extents_std or [])
        res.temporal_total_residual_norm = temporal.total_residual_norm
        res.temporal_n_time_points = len(temporal.timestamps or [])
        res.temporal_residence_time_flags = list(temporal.residence_time_flags or [])
        res.temporal_residence_time_details = dict(
            temporal.residence_time_details or {}
        )


def fit_network_with_priors(
    samples: object,
    edges: Iterable[EdgeInput],
    config: Config,
    *,
    auto_disable_missing: bool = True,
    physics_priors: Optional[Iterable[PhysicsPrior]] = None,
    physics_priors_mode: str = "override",
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
    residence_time_overrides: Optional[Mapping[str, float]] = None,
) -> Tuple[List[EdgeResult], List[Edge]]:
    """Fit a network after applying physics priors to the edge set."""
    if getattr(config, "strict_input_validation", False):
        validate_required_inputs(samples, config)
    elif auto_disable_missing:
        config = auto_disable_missing_modules(samples, config)
    built_edges = build_edges(edges)
    if physics_priors:
        built_edges = apply_physics_priors(
            built_edges, physics_priors, mode=physics_priors_mode
        )
    results = fit_network(
        samples,
        built_edges,
        config,
        phreeqc_results=phreeqc_results,
        residence_time_overrides=residence_time_overrides,
    )
    return results, built_edges


def fit_network_pipeline(
    samples: object,
    edges: Iterable[EdgeInput],
    config: Config,
    *,
    auto_disable_missing: bool = True,
    physics_priors: Optional[Iterable[PhysicsPrior]] = None,
    physics_priors_mode: str = "override",
    temporal_nodes: Optional[Mapping[str, TemporalNode]] = None,
    temporal_hydraulic_params: Optional[Mapping[str, Dict[str, float]]] = None,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
) -> Tuple[List[EdgeResult], Dict[str, object]]:
    """Run a connected pipeline with optional physics priors and temporal fits."""
    if getattr(config, "strict_input_validation", False):
        validate_required_inputs(samples, config)
    elif auto_disable_missing:
        config = auto_disable_missing_modules(samples, config)
        
    # Latent Endmembers (Ultra Upgrade)
    virtual_nodes = []
    if getattr(config, "latent_endmembers_enabled", False):
        try:
            # We assume samples is a list of dicts. If it's a dataframe, this might fail.
            # But the contract usually expects list of mappings.
            sample_list = _sample_list(samples)
            virtual_nodes = identify_latent_endmembers(
                sample_list, 
                config.ion_order, 
                n_endmembers=getattr(config, "latent_endmembers_count", 2)
            )
            # We must inject these into the samples object used by fit_network
            # But fit_network takes 'samples'. If it's a list, we append.
            if isinstance(samples, list):
                samples.extend(virtual_nodes)
            logger.info(f"Injected {len(virtual_nodes)} latent virtual nodes.")
        except Exception as e:
            logger.warning(f"Latent endmember identification failed: {e}")

    built_edges = build_edges(edges)
    
    # Connect Virtual Nodes? 
    # If we added virtual nodes, we should add candidate edges from them to all other nodes.
    if virtual_nodes:
        # We need to know who is 'downstream'. 
        # Heuristic: Virtual nodes are sources. They connect to everyone.
        # But this explodes the graph.
        # Let's rely on the user to run candidate generation again OR we add them here.
        # For 'fit_network_pipeline', edges are usually already provided.
        # So we append edges from virtual nodes to all real nodes?
        # That's O(N_virtual * N_real). Acceptable for small N_virtual.
        sample_list = _sample_list(samples)
        for vn in virtual_nodes:
            uid = vn["site_id"]
            for s in sample_list:
                vid = s.get("site_id")
                if vid and vid != uid:
                    # Add edge with a penalty? Or just add it.
                    # We create a new Edge object
                    from .graph.types import Edge
                    built_edges.append(Edge(edge_id=f"{uid}->{vid}", u=str(uid), v=str(vid)))

    if physics_priors:
        built_edges = apply_physics_priors(
            built_edges, physics_priors, mode=physics_priors_mode
        )

    temporal_results: List[TemporalEdgeResult] = []
    residence_time_overrides: Optional[Dict[str, float]] = None
    graph = None
    if temporal_nodes is not None:

        temporal_results, residence_time_overrides = fit_temporal_edges(
            temporal_nodes,
            built_edges,
            config,
            hydraulic_params_by_edge=temporal_hydraulic_params,
        )

    results = fit_network(
        samples,
        built_edges,
        config,
        phreeqc_results=phreeqc_results,
        residence_time_overrides=residence_time_overrides,
    )
    if temporal_results:
        attach_temporal_results(results, temporal_results)

    # 3. Nuclear Aging (Network-Enhanced Bayesian)
    nuclear_results = None
    if getattr(config, "sheaf_age_enabled", False):
        # Build DiGraph for the aging solver
        import networkx as nx
        graph = nx.DiGraph()
        for edge in built_edges:
            graph.add_edge(edge.u, edge.v, length_m=edge.attrs.get("length_m", 1.0))
        
        # Prepare observations
        node_obs = {}
        observed_sample_years: List[float] = []
        # Nuclear Tracer: default to Tritium if not specified
        tracer_name = getattr(config, "residence_time_tracer", "3H")
        
        sample_list = _sample_list(samples)
        sample_map = {}
        for sample in sample_list:
            node_id = sample.get("site_id") or sample.get("sample_id")
            if node_id is not None:
                sample_map[node_id] = sample
        
        for node_id, sample in sample_map.items():
            val = parse_numeric(sample.get(tracer_name), config.detection_limit_policy)
            if val is not None:
                node_obs[node_id] = val
                sample_year = _extract_sample_decimal_year(sample)
                if sample_year is not None:
                    observed_sample_years.append(sample_year)
        
        if node_obs:
            try:
                if observed_sample_years:
                    sample_date = float(median(observed_sample_years))
                else:
                    sample_date = float(datetime.now().year)
                    logger.warning(
                        "No valid sample dates found for nuclear tracer '%s'; using current year %.1f as fallback.",
                        tracer_name,
                        sample_date,
                    )
                
                nuclide = get_nuclide(tracer_name)
                if nuclide:
                    nuclear_results = infer_network_ages_bayesian(
                        graph,
                        node_obs,
                        {}, # sigmas auto-calculated if empty
                        sample_date,
                        nuclide=nuclide,
                        model_type=getattr(config, "nuclear_model", "PFM")
                    )
            except Exception as exc:
                logger.warning(
                    "Nuclear aging inference failed and was skipped: %s",
                    exc,
                    exc_info=True,
                )

    extras = {
        "edges": built_edges,
        "temporal_results": temporal_results,
        "residence_time_overrides": residence_time_overrides or {},
        "nuclear_results": nuclear_results,
        "graph": locals().get("graph") # pass graph for plotting
    }
    return results, extras
