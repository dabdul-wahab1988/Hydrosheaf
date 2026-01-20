"""High-level API helpers for optional Hydrosheaf modules."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from .config import Config
from .data.schema import parse_numeric
from .graph.build import EdgeInput, build_edges
from .graph.types import Edge
from .inference.edge_fit import EdgeResult
from .inference.network_fit import fit_network
from .physics.priors import PhysicsPrior, apply_physics_priors
from .temporal import TemporalEdgeResult, TemporalNode
from .temporal.temporal_edge_fit import fit_temporal_edge
from .vadose.contracts import (
    VadoseForcingSample,
    VadoseLinksRow,
    VadoseProfile,
    VadoseRunConfig,
)
from .vadose.run import build_vadose_edge_priors


def _sample_list(
    samples: Union[Mapping[str, Any], Sequence[Any]]
) -> List[Mapping[str, Any]]:
    if isinstance(samples, Mapping):
        return list(samples.values())
    if isinstance(samples, Sequence):
        return list(samples)
    raise TypeError("Unsupported samples input type.")


def _any_numeric(
    samples: Sequence[Mapping[str, object]],
    key: str,
    detection_policy: str,
) -> bool:
    for sample in samples:
        if parse_numeric(sample.get(key), detection_policy) is not None:
            return True
    return False


def _any_pair_numeric(
    samples: Sequence[Mapping[str, object]],
    key_a: str,
    key_b: str,
    detection_policy: str,
) -> bool:
    for sample in samples:
        if (
            parse_numeric(sample.get(key_a), detection_policy) is not None
            and parse_numeric(sample.get(key_b), detection_policy) is not None
        ):
            return True
    return False


def _missing_keys(
    samples: Sequence[Mapping[str, object]],
    keys: Sequence[str],
    detection_policy: str,
) -> List[str]:
    missing_ids: List[str] = []
    for sample in samples:
        missing = False
        for key in keys:
            if parse_numeric(sample.get(key), detection_policy) is None:
                missing = True
                break
        if missing:
            sample_id = sample.get("site_id") or sample.get("sample_id") or "unknown"
            missing_ids.append(str(sample_id))
    return missing_ids


def validate_required_inputs(samples: object, config: Config) -> None:
    """Raise if required inputs for enabled modules are missing."""
    sample_list = _sample_list(samples)
    detection_policy = config.detection_limit_policy
    missing_reports: List[str] = []

    if config.phreeqc_enabled:
        required_phreeqc = [
            "pH",
            "Ca",
            "Mg",
            "Na",
            "K",
            "Cl",
            "SO4",
            "NO3",
            "F",
            "HCO3",
        ]
        missing = _missing_keys(sample_list, required_phreeqc, detection_policy)
        if missing:
            missing_reports.append(
                "PHREEQC requires pH and major ions (Ca, Mg, Na, K, Cl, SO4, NO3, F, HCO3) "
                f"for all samples (missing: {missing})"
            )

    if config.isotope_enabled and config.lmwl_defined:
        missing = _missing_keys(
            sample_list,
            [config.isotope_d18o_key, config.isotope_d2h_key],
            detection_policy,
        )
        if missing:
            missing_reports.append(
                "Isotope penalties require both "
                f"{config.isotope_d18o_key} and {config.isotope_d2h_key} "
                f"for all samples (missing: {missing})"
            )

    if config.nitrate_source_enabled:
        missing = _missing_keys(sample_list, ["NO3"], detection_policy)
        if missing:
            missing_reports.append(
                f"Nitrate source requires NO3 for all samples (missing: {missing})"
            )

    if missing_reports:
        raise ValueError("; ".join(missing_reports))


def auto_disable_missing_modules(samples: object, config: Config) -> Config:
    """Disable feature flags when required inputs are missing across samples."""
    sample_list = _sample_list(samples)
    detection_policy = config.detection_limit_policy
    updates: Dict[str, object] = {}

    if config.phreeqc_enabled and not _any_numeric(sample_list, "pH", detection_policy):
        updates["phreeqc_enabled"] = False

    if config.isotope_enabled and not _any_pair_numeric(
        sample_list,
        config.isotope_d18o_key,
        config.isotope_d2h_key,
        detection_policy,
    ):
        updates["isotope_enabled"] = False

    if config.nitrate_source_enabled and not _any_numeric(
        sample_list, "NO3", detection_policy
    ):
        updates["nitrate_source_enabled"] = False

    if updates:
        return replace(config, **updates)
    return config


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
    built_edges = build_edges(edges)
    if physics_priors:
        built_edges = apply_physics_priors(
            built_edges, physics_priors, mode=physics_priors_mode
        )

    temporal_results: List[TemporalEdgeResult] = []
    residence_time_overrides: Optional[Dict[str, float]] = None
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

    extras = {
        "edges": built_edges,
        "temporal_results": temporal_results,
        "residence_time_overrides": residence_time_overrides or {},
    }
    return results, extras
