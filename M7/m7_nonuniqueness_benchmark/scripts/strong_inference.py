"""Truth-blind HydroSheaf inference contract for M7.2."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass
import math
from typing import Dict, Mapping, Optional, Sequence, Tuple

import networkx as nx
import numpy as np

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.network_fit import fit_network, infer_edges
from hydrosheaf.inference.topology_posterior import run_topology_posterior
from hydrosheaf.nuclear.network_aging import (
    TracerObservationSet,
    infer_network_ages_bayesian,
)
from hydrosheaf.nuclear.input_history import InputHistory
from hydrosheaf.nuclear.nuclides import ARGON39
from hydrosheaf.phreeqc.runner import run_phreeqc
from hydrosheaf.sheaf.topology_refine import refine_edges_with_sheaf

from independent_modflow_generator import ION_ORDER

FUSION_FEATURES = (
    "hydraulic_logit",
    "negative_age_cost",
    "negative_chemistry_log_objective",
)


@dataclass(frozen=True)
class StrongInferenceResult:
    edge_features: Tuple[Dict[str, object], ...]
    candidate_edges: Tuple[str, ...]
    age_results: Dict[str, Dict[str, float]]
    age_diagnostics: Dict[str, object]
    phreeqc_diagnostics: Dict[str, object]
    posterior_diagnostics: Dict[str, object]
    module_calls: Dict[str, int]

    def serializable(self) -> Dict[str, object]:
        return asdict(self)


def strong_config(
    *,
    phreeqc_enabled: bool = False,
    measured_ions: Optional[Sequence[str]] = None,
) -> Config:
    selected_ions = set(measured_ions or ION_ORDER)
    unknown = selected_ions - set(ION_ORDER)
    if unknown:
        raise ValueError(f"Unknown chemistry ions: {sorted(unknown)}")
    config = Config(
        ion_order=list(ION_ORDER),
        weights=[1.0 if ion in selected_ions else 0.0 for ion in ION_ORDER],
        conservative_weights=[
            (1.0 if ion == "Cl" else 0.01) if ion in selected_ions else 0.0
            for ion in ION_ORDER
        ],
        measured_ions=[ion for ion in ION_ORDER if ion in selected_ions],
        phreeqc_enabled=phreeqc_enabled,
        missing_policy="impute_zero",
        uncertainty_method="none",
        nitrate_source_enabled=False,
        edge_radius_km=5.0,
        edge_max_neighbors=8,
        edge_p_min=0.01,
        edge_head_key="head_meas",
        edge_sigma_meas=0.20,
        edge_sigma_topo=1.0,
        edge_gradient_min=0.0,
        sheaf_age_enabled=True,
        sheaf_isotope_enabled=False,
        sheaf_cl_enabled=False,
        sheaf_weight_head_prior=0.0,
        sheaf_weight_isotope=0.0,
        sheaf_weight_cl=0.0,
        sheaf_weight_age=1.0,
        sheaf_weight_global=0.0,
        sheaf_age_velocity_m_year=50.0,
        sheaf_age_travel_time_cv=0.75,
        sheaf_age_process_sigma_years=4.0,
        sheaf_age_default_sigma_years=10.0,
        sheaf_max_iter=1,
        sheaf_soft_beta=1.0,
        active_minerals=[
            "calcite",
            "dolomite",
            "gypsum",
            "albite",
            "anorthite",
            "k_feldspar",
        ],
        reaction_processes_enabled=[
            "denitrification",
            "sulfate_reduction",
            "iron_reduction",
        ],
        exchange_enabled=True,
        transport_models_enabled=["evap"],
        allow_signed_reactions=True,
        lambda_l1=0.01,
        topology_posterior_samples=2500,
        topology_posterior_burnin=750,
        topology_posterior_chains=4,
        topology_posterior_beta=1.0e-9,
        topology_posterior_edge_penalty=0.25,
        topology_posterior_min_edges=9,
        topology_posterior_max_out_degree=3,
        topology_posterior_require_acyclic=True,
        topology_posterior_initialization_steps=250,
        topology_posterior_gibbs_probability=0.70,
        topology_posterior_updates_per_sample=16,
    )
    config.validate()
    return config


def _bounded_probability(value: float) -> float:
    return min(1.0 - 1e-6, max(1e-6, float(value)))


def _reaction_family(label: Optional[str]) -> str:
    if not label:
        return "none"
    lower = label.lower()
    if any(token in lower for token in ("calcite", "dolomite", "carbonate")):
        return "carbonate"
    if any(
        token in lower
        for token in (
            "albite",
            "anorthite",
            "feldspar",
            "silicate",
            "exch",
        )
    ):
        return "silicate_exchange"
    if "sulfate_reduction" in lower:
        return "sulfate_reduction"
    if "iron_reduction" in lower:
        return "iron_reduction"
    if "denit" in lower:
        return "denitrification"
    if "pyrite" in lower:
        return "other_redox"
    if "gypsum" in lower or "so4" in lower:
        return "sulfate_source"
    return "other"


def _fusion_probability(
    row: Mapping[str, object],
    fusion_model: Mapping[str, object],
) -> float:
    feature_names = tuple(fusion_model.get("feature_names", FUSION_FEATURES))
    means = np.asarray(fusion_model.get("means", [0.0] * len(feature_names)), float)
    scales = np.asarray(fusion_model.get("scales", [1.0] * len(feature_names)), float)
    coefficients = np.asarray(fusion_model["coefficients"], float)
    vector = np.asarray([float(row[name]) for name in feature_names], float)
    standardized = (vector - means) / np.where(scales > 0.0, scales, 1.0)
    linear = float(fusion_model["intercept"]) + float(standardized @ coefficients)
    probability = 1.0 / (1.0 + math.exp(-max(-40.0, min(40.0, linear))))
    if fusion_model.get("kind") == "age_compatibility_gate":
        threshold = float(fusion_model["age_cost_max"])
        if float(row["age_cost"]) > threshold:
            probability = min(
                probability,
                float(fusion_model.get("incompatible_probability", 1.0e-6)),
            )
    return probability


def _feasible_initial_edges(
    edges: Sequence[Edge],
    probability_by_id: Mapping[str, float],
    n_edges: int,
) -> list[Edge]:
    graph = nx.DiGraph()
    selected: list[Edge] = []
    for edge in sorted(
        edges,
        key=lambda item: probability_by_id[item.edge_id],
        reverse=True,
    ):
        if graph.has_node(edge.u) and int(graph.out_degree[edge.u]) >= 3:
            continue
        graph.add_edge(edge.u, edge.v)
        if not nx.is_directed_acyclic_graph(graph):
            graph.remove_edge(edge.u, edge.v)
            continue
        selected.append(edge)
        if len(selected) >= n_edges:
            return selected
    raise RuntimeError("Could not construct a feasible topology-posterior start.")


def run_strong_inference(
    observations: Sequence[Mapping[str, object]],
    seed: int,
    *,
    fusion_model: Optional[Mapping[str, object]] = None,
    run_posterior: bool = False,
    age_draws: int = 500,
    age_chains: int = 4,
    topology_samples: Optional[int] = None,
    topology_updates_per_sample: Optional[int] = None,
    age_travel_cost_weight: Optional[float] = None,
    chemistry_ions: Optional[Sequence[str]] = None,
) -> StrongInferenceResult:
    """Run HydroSheaf using observables only."""

    rows = [dict(row) for row in observations]
    forbidden = {
        "true_age",
        "true_edges",
        "true_process",
        "true_age_hidden_from_inference",
    }
    leaked = sorted(forbidden & set().union(*(set(row) for row in rows)))
    if leaked:
        raise ValueError(f"Truth fields are forbidden in inference rows: {leaked}")

    calls = {
        "infer_edges": 0,
        "infer_network_ages_bayesian": 0,
        "refine_edges_with_sheaf": 0,
        "run_phreeqc": 0,
        "fit_network_unconstrained": 0,
        "fit_network_constrained": 0,
        "run_topology_posterior": 0,
    }
    base_config = strong_config(
        phreeqc_enabled=False,
        measured_ions=chemistry_ions,
    )
    if topology_samples is not None:
        base_config.topology_posterior_samples = int(topology_samples)
        base_config.topology_posterior_burnin = min(
            base_config.topology_posterior_burnin,
            max(0, int(topology_samples) // 3),
        )
    if topology_updates_per_sample is not None:
        base_config.topology_posterior_updates_per_sample = int(
            topology_updates_per_sample
        )
    if age_travel_cost_weight is not None:
        base_config.sheaf_age_travel_cost_weight = float(age_travel_cost_weight)
    base_config.validate()
    candidates = infer_edges(rows, method="probabilistic", config=base_config)
    calls["infer_edges"] += 1

    independent_graph = nx.DiGraph()
    independent_graph.add_nodes_from(str(row["site_id"]) for row in rows)
    tritium = {str(row["site_id"]): float(row["tritium_TU"]) for row in rows}
    tritium_sigma = {node: max(0.10, 0.12 * value) for node, value in tritium.items()}
    argon = {str(row["site_id"]): float(row["argon39_pmc"]) for row in rows}
    argon_sigma = {node: 1.8 for node in argon}
    # The synthetic experiment treats the recharge-boundary tracer forcing as
    # known, just as a field inversion would use an independently measured
    # atmospheric/recharge history.  The generator still uses different
    # nonlinear response equations and adds unmodelled process discrepancy.
    steady_tritium_boundary = InputHistory(
        np.asarray([1850.0, 2035.0], dtype=float),
        np.asarray([6.2, 6.2], dtype=float),
    )
    age_results = infer_network_ages_bayesian(
        independent_graph,
        tritium,
        tritium_sigma,
        sample_date=2025.5,
        input_hist=steady_tritium_boundary,
        additional_tracers=[
            TracerObservationSet(
                name="39Ar",
                nuclide=ARGON39,
                observations=argon,
                sigmas=argon_sigma,
                sample_date=2025.5,
                reference_concentration=97.0,
            )
        ],
        n_samples=int(age_draws),
        n_chains=int(age_chains),
        random_seed=int(seed),
        sampler="grid",
        max_age_years=200.0,
        root_age_prior_median_years=20.0,
        root_age_prior_log_sigma=1.2,
    )
    calls["infer_network_ages_bayesian"] += 1
    age_diagnostics = dict(age_results["_diagnostics"])

    age_rows = deepcopy(rows)
    for row in age_rows:
        node_id = str(row["site_id"])
        posterior = age_results[node_id]
        row["mean_age_years"] = posterior["mean_age_years"]
        row["mean_age_std_years"] = posterior["std_age_years"]
        row["tracer_identifiable"] = bool(posterior["tracer_identifiable"])
    age_edges = refine_edges_with_sheaf(
        age_rows,
        deepcopy(candidates),
        base_config,
    )
    calls["refine_edges_with_sheaf"] += 1
    age_by_id = {edge.edge_id: edge for edge in age_edges}

    unconstrained_config = strong_config(
        phreeqc_enabled=False,
        measured_ions=chemistry_ions,
    )
    unconstrained_results = fit_network(
        rows, deepcopy(candidates), unconstrained_config
    )
    calls["fit_network_unconstrained"] += 1
    phreeqc_config = strong_config(
        phreeqc_enabled=True,
        measured_ions=chemistry_ions,
    )
    phreeqc_results = run_phreeqc(rows, phreeqc_config)
    calls["run_phreeqc"] += 1
    constrained_results = fit_network(
        rows,
        deepcopy(candidates),
        phreeqc_config,
        phreeqc_results=phreeqc_results,
    )
    calls["fit_network_constrained"] += 1
    unconstrained_by_id = {result.edge_id: result for result in unconstrained_results}
    constrained_by_id = {result.edge_id: result for result in constrained_results}

    edge_features = []
    for edge in candidates:
        attrs = edge.attrs or {}
        hydraulic_probability = _bounded_probability(
            float(
                attrs.get(
                    "prior_edge_probability",
                    attrs.get("p_uv", attrs.get("edge_confidence", 0.5)),
                )
            )
        )
        age_edge = age_by_id.get(edge.edge_id)
        age_cost = float(
            (age_edge.attrs or {}).get("sheaf_cost_age", 0.0)
            if age_edge is not None
            else 0.0
        )
        unconstrained = unconstrained_by_id[edge.edge_id]
        constrained = constrained_by_id[edge.edge_id]
        chemistry_objective = max(1e-12, float(constrained.objective_score))
        extent_map = {
            str(label): float(extent)
            for label, extent in zip(constrained.z_labels, constrained.z_extents)
        }
        dominant_label = (
            max(extent_map, key=lambda label: abs(extent_map[label]))
            if extent_map
            else None
        )
        unconstrained_extent_map = {
            str(label): float(extent)
            for label, extent in zip(unconstrained.z_labels, unconstrained.z_extents)
        }
        unconstrained_dominant = (
            max(
                unconstrained_extent_map,
                key=lambda label: abs(unconstrained_extent_map[label]),
            )
            if unconstrained_extent_map
            else None
        )
        row: Dict[str, object] = {
            "edge_id": edge.edge_id,
            "u": edge.u,
            "v": edge.v,
            "distance_km": float(attrs.get("distance_km", 0.0)),
            "hydraulic_probability": hydraulic_probability,
            "hydraulic_logit": math.log(
                hydraulic_probability / (1.0 - hydraulic_probability)
            ),
            "age_cost": age_cost,
            # Compress rare, very large reversal costs so a few grossly
            # incompatible candidates cannot determine the entire fusion
            # coefficient.  This monotone transform preserves the ranking.
            "negative_age_cost": -math.log1p(age_cost),
            "chemistry_objective_unconstrained": float(unconstrained.objective_score),
            "chemistry_objective_constrained": chemistry_objective,
            "negative_chemistry_log_objective": -math.log1p(chemistry_objective),
            "phreeqc_objective_change": float(
                constrained.objective_score - unconstrained.objective_score
            ),
            "phreeqc_ok": bool(constrained.phreeqc_ok),
            "thermodynamic_constraints_active_count": int(
                constrained.thermodynamic_constraints_active_count
            ),
            "thermodynamic_bound_hit_count": int(
                constrained.thermodynamic_bound_hit_count
            ),
            "constraints_binding": dict(constrained.constraints_binding),
            "reaction_extents": extent_map,
            "dominant_reaction": dominant_label,
            "dominant_reaction_family": _reaction_family(dominant_label),
            "dominant_reaction_unconstrained": unconstrained_dominant,
            "dominant_reaction_family_unconstrained": _reaction_family(
                unconstrained_dominant
            ),
        }
        if fusion_model is not None:
            row["fusion_probability"] = _fusion_probability(row, fusion_model)
        edge_features.append(row)

    common_ids = sorted(set(unconstrained_by_id) & set(constrained_by_id))
    objective_changes = np.asarray(
        [
            constrained_by_id[edge_id].objective_score
            - unconstrained_by_id[edge_id].objective_score
            for edge_id in common_ids
        ],
        dtype=float,
    )
    phreeqc_diagnostics: Dict[str, object] = {
        "n_samples": len(phreeqc_results),
        "n_successful_samples": int(
            sum(bool(value.get("phreeqc_ok")) for value in phreeqc_results.values())
        ),
        "success_fraction": float(
            np.mean(
                [bool(value.get("phreeqc_ok")) for value in phreeqc_results.values()]
            )
        ),
        "n_candidate_edge_fits": len(common_ids),
        "n_edges_with_active_direction_constraints": int(
            sum(
                result.thermodynamic_constraints_active_count > 0
                for result in constrained_by_id.values()
            )
        ),
        "n_edges_hitting_thermodynamic_bound": int(
            sum(
                result.thermodynamic_bound_hit_count > 0
                for result in constrained_by_id.values()
            )
        ),
        "n_edges_with_material_objective_change": int(
            np.sum(np.abs(objective_changes) > 1e-6)
        ),
        "mean_absolute_objective_change": float(np.mean(np.abs(objective_changes))),
        "max_absolute_objective_change": float(np.max(np.abs(objective_changes))),
    }

    posterior_diagnostics: Dict[str, object] = {}
    if run_posterior:
        if fusion_model is None:
            raise ValueError("A frozen fusion_model is required for posterior runs.")
        probability_by_id = {
            str(row["edge_id"]): float(row["fusion_probability"])
            for row in edge_features
        }
        posterior_universe = deepcopy(candidates)
        for edge in posterior_universe:
            edge.attrs = dict(edge.attrs or {})
            edge.attrs["prior_edge_probability"] = _bounded_probability(
                probability_by_id[edge.edge_id]
            )
        initial = _feasible_initial_edges(
            posterior_universe,
            probability_by_id,
            n_edges=9,
        )
        posterior_diagnostics = run_topology_posterior(
            posterior_universe,
            cost_fn=lambda selected: 0.0,
            config=base_config,
            initial_edges=initial,
            seed=int(seed),
        )
        calls["run_topology_posterior"] += 1

    return StrongInferenceResult(
        edge_features=tuple(edge_features),
        candidate_edges=tuple(edge.edge_id for edge in candidates),
        age_results=age_results,
        age_diagnostics=age_diagnostics,
        phreeqc_diagnostics=phreeqc_diagnostics,
        posterior_diagnostics=posterior_diagnostics,
        module_calls=calls,
    )
