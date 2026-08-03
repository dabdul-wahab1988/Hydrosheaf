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
from .sheaf.topology_refine import refine_edges_with_sheaf

from .vadose.contracts import (
    VadoseForcingSample,
    VadoseLinksRow,
    VadoseProfile,
    VadoseRunConfig,
)
from .vadose.run import build_vadose_edge_priors

logger = get_logger("api")


class PipelineStageError(RuntimeError):
    """Raised when an explicitly required pipeline stage cannot complete."""

    def __init__(self, stage: str, message: str):
        self.stage = str(stage)
        super().__init__(f"Pipeline stage '{self.stage}' did not complete: {message}")


_PIPELINE_STAGES = (
    "latent_endmembers",
    "temporal",
    "nuclear_age",
    "sheaf_refinement",
    "network_fit",
)


def _stage_record(
    status: str,
    *,
    requested: bool,
    detail: str = "",
    error: Optional[BaseException] = None,
) -> Dict[str, object]:
    record: Dict[str, object] = {
        "status": status,
        "requested": bool(requested),
    }
    if detail:
        record["detail"] = detail
    if error is not None:
        record["error_type"] = type(error).__name__
        record["error"] = str(error)
    return record


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
    sheaf_refinement_enabled: bool = False,
    nuclear_inference_options: Optional[Mapping[str, object]] = None,
    strict_stage_completion: bool = False,
    required_stages: Optional[Iterable[str]] = None,
) -> Tuple[List[EdgeResult], Dict[str, object]]:
    """Run the connected inference pipeline and report explicit stage status.

    ``strict_stage_completion`` converts a requested optional-stage skip or
    failure into :class:`PipelineStageError`. ``sheaf_refinement_enabled`` is
    opt-in for backward compatibility; when age evidence is enabled, nuclear
    posteriors are attached before candidate edges are sheaf-refined.
    """
    requested = {
        "latent_endmembers": bool(
            getattr(config, "latent_endmembers_enabled", False)
        ),
        "temporal": temporal_nodes is not None,
        "nuclear_age": bool(getattr(config, "sheaf_age_enabled", False)),
        "sheaf_refinement": bool(
            sheaf_refinement_enabled
            or getattr(config, "topology_posterior_enabled", False)
        ),
        "network_fit": True,
    }
    if required_stages is None:
        required = {name for name, is_requested in requested.items() if is_requested}
    else:
        required = {str(name) for name in required_stages}
        unknown = required - set(_PIPELINE_STAGES)
        if unknown:
            raise ValueError(f"Unknown required pipeline stages: {sorted(unknown)}")
        not_requested = required - {
            name for name, is_requested in requested.items() if is_requested
        }
        if not_requested:
            raise ValueError(
                "Required pipeline stages were not requested: "
                f"{sorted(not_requested)}"
            )

    stage_status = {
        name: _stage_record("not_requested", requested=is_requested)
        for name, is_requested in requested.items()
    }

    def _fail_stage(stage: str, message: str, exc: Optional[BaseException] = None) -> None:
        stage_status[stage] = _stage_record(
            "failed", requested=True, detail=message, error=exc
        )
        posterior_fail_closed = bool(
            stage == "sheaf_refinement"
            and getattr(config, "topology_posterior_enabled", False)
        )
        if (strict_stage_completion and stage in required) or posterior_fail_closed:
            raise PipelineStageError(stage, message) from exc

    if getattr(config, "strict_input_validation", False):
        validate_required_inputs(samples, config)
    elif auto_disable_missing:
        config = auto_disable_missing_modules(samples, config)

    # Never mutate a caller-owned list when virtual nodes or inferred age fields
    # are added for downstream stages.
    pipeline_samples: object = samples
    if isinstance(samples, list):
        pipeline_samples = [
            dict(sample) if isinstance(sample, Mapping) else sample for sample in samples
        ]

    # Latent Endmembers (Ultra Upgrade)
    virtual_nodes = []
    if requested["latent_endmembers"]:
        try:
            sample_list = _sample_list(pipeline_samples)
            virtual_nodes = identify_latent_endmembers(
                sample_list,
                config.ion_order,
                n_endmembers=getattr(config, "latent_endmembers_count", 2)
            )
            if isinstance(pipeline_samples, list):
                pipeline_samples.extend(virtual_nodes)
            else:
                pipeline_samples = sample_list + list(virtual_nodes)
            stage_status["latent_endmembers"] = _stage_record(
                "completed",
                requested=True,
                detail=f"identified {len(virtual_nodes)} virtual nodes",
            )
            logger.info(f"Injected {len(virtual_nodes)} latent virtual nodes.")
        except Exception as e:
            _fail_stage("latent_endmembers", str(e), e)
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
        sample_list = _sample_list(pipeline_samples)
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
        try:
            temporal_results, residence_time_overrides = fit_temporal_edges(
                temporal_nodes,
                built_edges,
                config,
                hydraulic_params_by_edge=temporal_hydraulic_params,
            )
            stage_status["temporal"] = _stage_record(
                "completed",
                requested=True,
                detail=f"fit {len(temporal_results)} temporal edges",
            )
        except Exception as exc:
            _fail_stage("temporal", str(exc), exc)
            logger.warning("Temporal inference failed and was skipped: %s", exc)

    # 3. Nuclear Aging (Network-Enhanced Bayesian)
    nuclear_results = None
    if requested["nuclear_age"]:
        # Build DiGraph for the aging solver
        import networkx as nx
        graph = nx.DiGraph()
        for edge in built_edges:
            graph.add_nodes_from((edge.u, edge.v))
            # When topology is what the sheaf is about to infer, conditioning
            # node ages on every candidate edge is circular and candidate
            # cycles can invalidate the DAG age model. Infer local ages on an
            # edge-free graph first; retain the legacy graph-conditioned path
            # when sheaf refinement is not requested.
            if not sheaf_refinement_enabled:
                graph.add_edge(
                    edge.u,
                    edge.v,
                    length_m=edge.attrs.get("length_m", 1.0),
                )
        
        # Prepare observations
        node_obs = {}
        observed_sample_years: List[float] = []
        # Nuclear Tracer: default to Tritium if not specified
        tracer_name = getattr(config, "residence_time_tracer", "3H")
        
        sample_list = _sample_list(pipeline_samples)
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
                    inference_options = dict(nuclear_inference_options or {})
                    reserved_options = {
                        "graph",
                        "node_observations",
                        "node_sigmas",
                        "sample_date",
                        "nuclide",
                        "model_type",
                    }
                    invalid_options = reserved_options & set(inference_options)
                    if invalid_options:
                        raise ValueError(
                            "nuclear_inference_options cannot override: "
                            f"{sorted(invalid_options)}"
                        )
                    nuclear_results = infer_network_ages_bayesian(
                        graph,
                        node_obs,
                        {}, # sigmas auto-calculated if empty
                        sample_date,
                        nuclide=nuclide,
                        model_type=getattr(config, "nuclear_model", "PFM"),
                        **inference_options,
                    )
                    if not isinstance(nuclear_results, Mapping) or not any(
                        node_id in nuclear_results for node_id in node_obs
                    ):
                        raise RuntimeError(
                            "age inference returned no node posterior results"
                        )
                    stage_status["nuclear_age"] = _stage_record(
                        "completed",
                        requested=True,
                        detail=f"inferred ages for {len(node_obs)} observed nodes",
                    )
                else:
                    raise RuntimeError(
                        f"unknown residence-time tracer '{tracer_name}'"
                    )
            except Exception as exc:
                _fail_stage("nuclear_age", str(exc), exc)
                logger.warning(
                    "Nuclear aging inference failed and was skipped: %s",
                    exc,
                    exc_info=True,
                )
        else:
            message = f"no numeric observations for tracer '{tracer_name}'"
            stage_status["nuclear_age"] = _stage_record(
                "skipped", requested=True, detail=message
            )
            if strict_stage_completion and "nuclear_age" in required:
                raise PipelineStageError("nuclear_age", message)

    if requested["sheaf_refinement"]:
        try:
            sheaf_samples = [
                dict(row) if isinstance(row, Mapping) else row
                for row in _sample_list(pipeline_samples)
            ]
            if requested["nuclear_age"]:
                if nuclear_results is None:
                    raise RuntimeError(
                        "age evidence was requested but nuclear inference did not complete"
                    )
                for row in sheaf_samples:
                    if not isinstance(row, dict):
                        continue
                    node_id = row.get("site_id") or row.get("sample_id")
                    posterior = (
                        nuclear_results.get(node_id)
                        if isinstance(nuclear_results, Mapping)
                        else None
                    )
                    if not isinstance(posterior, Mapping):
                        continue
                    if "mean_age_years" in posterior:
                        row["mean_age_years"] = posterior["mean_age_years"]
                    if "std_age_years" in posterior:
                        row["mean_age_std_years"] = posterior["std_age_years"]
                    if "tracer_identifiable" in posterior:
                        row["tracer_identifiable"] = posterior[
                            "tracer_identifiable"
                        ]
            refined_edges = refine_edges_with_sheaf(
                sheaf_samples,
                built_edges,
                config,
            )
            if not refined_edges:
                raise RuntimeError("sheaf refinement returned no edges")
            built_edges = refined_edges
            stage_status["sheaf_refinement"] = _stage_record(
                "completed",
                requested=True,
                detail=f"retained {len(built_edges)} refined edges",
            )
        except Exception as exc:
            _fail_stage("sheaf_refinement", str(exc), exc)
            logger.warning("Sheaf refinement failed: %s", exc)

    try:
        results = fit_network(
            pipeline_samples,
            built_edges,
            config,
            phreeqc_results=phreeqc_results,
            residence_time_overrides=residence_time_overrides,
        )
        stage_status["network_fit"] = _stage_record(
            "completed",
            requested=True,
            detail=f"fit {len(results)} edges",
        )
    except Exception as exc:
        _fail_stage("network_fit", str(exc), exc)
        raise
    if temporal_results:
        attach_temporal_results(results, temporal_results)

    extras = {
        "edges": built_edges,
        "temporal_results": temporal_results,
        "residence_time_overrides": residence_time_overrides or {},
        "nuclear_results": nuclear_results,
        "graph": graph,
        "stage_status": stage_status,
        "strict_stage_completion": bool(strict_stage_completion),
        "required_stages": sorted(required),
        "virtual_nodes": virtual_nodes,
    }
    return results, extras
