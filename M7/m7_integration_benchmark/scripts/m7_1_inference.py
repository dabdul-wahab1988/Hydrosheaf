"""Truth-blind HydroSheaf execution contract for M7.1.

No truth-model module is imported here.  The public function accepts observations,
a frozen threshold, and a seed; all returned scores are derived from observables
and actual HydroSheaf module outputs.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, replace
from typing import Dict, List, Mapping, Sequence, Tuple

import math
import networkx as nx
import warnings

from hydrosheaf.config import Config
from hydrosheaf.data.minerals import get_mineral_stoich
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.network_fit import fit_network, infer_edges
from hydrosheaf.inference.topology_posterior import infer_topology_map_edges
from hydrosheaf.nuclear.network_aging import infer_network_ages_bayesian
from hydrosheaf.phreeqc.runner import run_phreeqc
from hydrosheaf.sheaf.topology_refine import refine_edges_with_sheaf


TRAINING_IONS = ["Ca", "Na", "K", "HCO3", "Cl", "NO3", "F", "PO4"]
HELD_OUT_IONS = ["Mg", "SO4", "Fe"]


@dataclass(frozen=True)
class BlindInferenceResult:
    edge_scores: Tuple[Dict[str, object], ...]
    selected_edges: Tuple[str, ...]
    candidate_edges: Tuple[str, ...]
    module_calls: Dict[str, int]
    posterior_diagnostics: Dict[str, object]
    bayesian_ages: Dict[str, Dict[str, float]]
    phreeqc_summary: Dict[str, object]

    def serializable(self) -> Dict[str, object]:
        return asdict(self)


def _edge_id(edge: Edge) -> str:
    return str(edge.edge_id)


def _softmax_by_source(
    edges: Sequence[Edge], raw: Mapping[str, float], lower_is_better: bool = False
) -> Dict[str, float]:
    grouped: Dict[str, List[Edge]] = {}
    for edge in edges:
        grouped.setdefault(edge.u, []).append(edge)
    probabilities: Dict[str, float] = {}
    for group in grouped.values():
        values = []
        for edge in group:
            value = float(raw.get(_edge_id(edge), 0.0))
            values.append(-value if lower_is_better else value)
        centre = sum(values) / len(values)
        scale = math.sqrt(sum((value - centre) ** 2 for value in values) / len(values))
        scale = max(scale, 1.0e-6)
        logits = [(value - centre) / scale for value in values]
        maximum = max(logits)
        weights = [math.exp(max(-40.0, min(40.0, value - maximum))) for value in logits]
        total = sum(weights)
        for edge, weight in zip(group, weights):
            probabilities[_edge_id(edge)] = weight / total
    return probabilities


def benchmark_config() -> Config:
    """Frozen M7.1 inference settings, selected without test-truth access."""

    config = Config(
        ion_order=list(TRAINING_IONS),
        weights=[1.0] * len(TRAINING_IONS),
        conservative_weights=[
            0.01 if ion != "Cl" else 1.0 for ion in TRAINING_IONS
        ],
        measured_ions=list(TRAINING_IONS),
        phreeqc_enabled=False,
        missing_policy="impute_zero",
        uncertainty_method="none",
        nitrate_source_enabled=False,
        edge_radius_km=5.5,
        # Candidate generation is intentionally broad.  Selection happens
        # downstream; this setting avoids censoring the immediately downstream
        # true edge merely because a larger head drop exists farther away.
        edge_max_neighbors=12,
        edge_p_min=0.05,
        edge_head_key="head_meas",
        edge_sigma_meas=1.2,
        edge_sigma_topo=2.5,
        edge_gradient_min=0.0,
        sheaf_soft_beta=1.0,
        sheaf_max_iter=2,
        topology_posterior_samples=5000,
        topology_posterior_burnin=1500,
        topology_posterior_chains=3,
        topology_posterior_beta=0.05,
        topology_posterior_min_edges=7,
        topology_posterior_max_out_degree=2,
        topology_posterior_require_acyclic=True,
        topology_posterior_require_weak_connectivity=True,
        topology_posterior_edge_penalty=0.05,
        topology_posterior_invalid_cost=1.0e5,
        active_minerals=["calcite", "albite", "gypsum"],
        transport_models_enabled=["evap"],
        lambda_l1=0.015,
    )
    config.validate()
    return config


def run_blind_inference(
    observations: Sequence[Mapping[str, object]],
    threshold: float,
    seed: int,
    *,
    run_posterior: bool = False,
    run_heavy_modules: bool = False,
    aging_draws: int = 1000,
    aging_chains: int = 4,
    fusion_model: Mapping[str, object] | None = None,
) -> BlindInferenceResult:
    """Run HydroSheaf without accepting any hidden truth fields."""

    rows = [dict(row) for row in observations]
    config = benchmark_config()
    calls = {
        "infer_edges": 0,
        "refine_edges_with_sheaf": 0,
        "fit_network": 0,
        "infer_topology_map_edges": 0,
        "run_phreeqc": 0,
        "infer_network_ages_bayesian": 0,
    }

    candidates = infer_edges(rows, method="probabilistic", config=config)
    calls["infer_edges"] += 1
    # Each ablation receives an independent copy to prevent mutable Edge attrs
    # from contaminating another method.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Recharge date .* exceeds input history start"
        )
        sheaf_edges = refine_edges_with_sheaf(
            rows, deepcopy(candidates), config
        )
        age_only_config = replace(
            config,
            sheaf_weight_head_prior=0.0,
            sheaf_weight_isotope=0.0,
            sheaf_weight_cl=0.0,
            sheaf_weight_age=1.0,
            sheaf_weight_global=0.0,
            sheaf_max_iter=0,
            sheaf_isotope_enabled=False,
            sheaf_cl_enabled=False,
        )
        age_only_edges = refine_edges_with_sheaf(
            rows, deepcopy(candidates), age_only_config
        )
    calls["refine_edges_with_sheaf"] += 2
    fits = fit_network(rows, deepcopy(candidates), config)
    calls["fit_network"] += 1

    fit_by_id = {result.edge_id: float(result.objective_score) for result in fits}
    fit_result_by_id = {result.edge_id: result for result in fits}
    chemistry_probability = _softmax_by_source(candidates, fit_by_id, lower_is_better=True)
    prior_raw = {}
    for edge in candidates:
        attrs = edge.attrs or {}
        probability = float(
            attrs.get("prior_edge_probability", attrs.get("p_uv", 0.5))
        )
        distance_km = float(attrs.get("distance_km", 0.0))
        prior_raw[_edge_id(edge)] = math.log(max(probability, 1.0e-9)) - (
            distance_km / 1.5
        )
    prior_probability = _softmax_by_source(candidates, prior_raw)
    sheaf_raw = {
        _edge_id(edge): float((edge.attrs or {}).get("sheaf_weight", 0.0))
        for edge in sheaf_edges
    }
    sheaf_probability = _softmax_by_source(candidates, sheaf_raw)
    age_only_raw = {
        _edge_id(edge): float((edge.attrs or {}).get("sheaf_weight", 0.0))
        for edge in age_only_edges
    }
    age_only_probability = _softmax_by_source(candidates, age_only_raw)

    edge_scores: List[Dict[str, object]] = []
    for edge in candidates:
        edge_id = _edge_id(edge)
        p_h = min(1.0 - 1.0e-6, max(1.0e-6, prior_probability[edge_id]))
        p_s = min(1.0 - 1.0e-6, max(1.0e-6, sheaf_probability.get(edge_id, 1.0e-6)))
        p_a = min(
            1.0 - 1.0e-6,
            max(1.0e-6, age_only_probability.get(edge_id, 1.0e-6)),
        )
        p_c = min(1.0 - 1.0e-6, max(1.0e-6, chemistry_probability.get(edge_id, 1.0e-6)))
        # Predeclared equal-weight geometric evidence pool.  Source-normalized
        # component probabilities keep a near-certain head direction from
        # swamping the independent age/sheaf and chemistry evidence streams.
        joint = (p_h * p_a * p_c) ** (1.0 / 3.0)
        logistic_joint = None
        if fusion_model is not None:
            coefficients = fusion_model.get("coefficients", {})
            if not isinstance(coefficients, Mapping):
                raise TypeError("fusion_model coefficients must be a mapping.")
            linear = float(fusion_model.get("intercept", 0.0))
            linear += float(coefficients.get("hydraulic_probability", 0.0)) * p_h
            linear += float(coefficients.get("age_only_probability", 0.0)) * p_a
            linear += float(coefficients.get("chemistry_probability", 0.0)) * p_c
            logistic_joint = 1.0 / (
                1.0 + math.exp(-max(-40.0, min(40.0, linear)))
            )
        fit_result = fit_result_by_id.get(edge_id)
        held_out_reaction_delta = {ion: 0.0 for ion in HELD_OUT_IONS}
        reaction_extents: Dict[str, float] = {}
        if fit_result is not None:
            reaction_extents = {
                str(label): float(extent)
                for label, extent in zip(
                    fit_result.z_labels, fit_result.z_extents
                )
            }
            for label, extent in reaction_extents.items():
                try:
                    coefficients = get_mineral_stoich(label)
                except ValueError:
                    coefficients = {
                        "NO3src": {"NO3": 1.0},
                        "denit": {"HCO3": config.denit_kappa, "NO3": -1.0},
                        "CaNa_exch": {"Ca": 1.0, "Na": -2.0},
                        "NaCa_exch": {"Ca": -1.0, "Na": 2.0},
                        "MgNa_exch": {"Mg": 1.0, "Na": -2.0},
                        "NaMg_exch": {"Mg": -1.0, "Na": 2.0},
                    }.get(label, {})
                for ion in HELD_OUT_IONS:
                    held_out_reaction_delta[ion] += extent * float(
                        coefficients.get(ion, 0.0)
                    )
        edge_scores.append(
            {
                "edge_id": edge_id,
                "u": edge.u,
                "v": edge.v,
                "hydraulic_probability": p_h,
                "age_only_probability": p_a,
                "sheaf_multievidence_probability": p_s,
                "chemistry_probability": p_c,
                "joint_probability": joint,
                "joint_logistic_probability": logistic_joint,
                "chemistry_objective": fit_by_id.get(edge_id),
                "reaction_extents": reaction_extents,
                # Strict holdout: neither endpoint's held-out measurements are
                # read here. This is the fitted reaction-only delta.
                "held_out_reaction_delta": held_out_reaction_delta,
            }
        )

    posterior_diagnostics: Dict[str, object] = {}
    if run_posterior and candidates:
        posterior_rank_key = (
            "joint_logistic_probability"
            if fusion_model is not None
            else "joint_probability"
        )
        posterior_items_by_target: Dict[str, List[Dict[str, object]]] = {}
        for item in edge_scores:
            posterior_items_by_target.setdefault(str(item["v"]), []).append(item)
        posterior_ids = {
            str(item["edge_id"])
            for items in posterior_items_by_target.values()
            for item in sorted(
                items,
                key=lambda row: float(row[posterior_rank_key]),
                reverse=True,
            )[:3]
        }
        posterior_universe = [
            deepcopy(edge)
            for edge in candidates
            if _edge_id(edge) in posterior_ids
        ]
        hydraulic_by_id = {
            str(item["edge_id"]): float(item["hydraulic_probability"])
            for item in edge_scores
        }
        for edge in posterior_universe:
            attrs = dict(edge.attrs or {})
            attrs["prior_edge_probability"] = min(
                0.98,
                max(0.02, hydraulic_by_id[_edge_id(edge)]),
            )
            edge.attrs = attrs
        candidates_by_id = {
            _edge_id(edge): edge for edge in posterior_universe
        }
        score_by_pair = {
            (str(item["u"]), str(item["v"])): item
            for item in edge_scores
            if str(item["edge_id"]) in posterior_ids
        }
        ordered_nodes = [
            str(row["site_id"])
            for row in sorted(
                rows,
                key=lambda row: float(row["head_meas"]),
                reverse=True,
            )
        ]
        earlier: set[str] = {ordered_nodes[0]}
        out_degree = {node: 0 for node in ordered_nodes}
        initial_ids: set[str] = set()
        for node in ordered_nodes[1:]:
            feasible = [
                item
                for (u, v), item in score_by_pair.items()
                if v == node and u in earlier and out_degree[u] < 2
            ]
            if not feasible:
                raise RuntimeError(
                    "Could not construct a feasible observable-derived "
                    f"posterior initial graph for node {node}."
                )
            parent_edge = max(
                feasible,
                key=lambda item: float(
                        item[
                        "joint_logistic_probability"
                        if fusion_model is not None
                        else "joint_probability"
                    ]
                ),
            )
            initial_ids.add(str(parent_edge["edge_id"]))
            out_degree[str(parent_edge["u"])] += 1
            earlier.add(node)
        posterior_initial = [
            deepcopy(candidates_by_id[edge_id])
            for edge_id in initial_ids
            if edge_id in candidates_by_id
        ]
        posterior_config = replace(
            config,
            edge_max_neighbors=2,
            topology_posterior_require_root_reachability=True,
            topology_posterior_root_nodes=[ordered_nodes[0]],
        )
        _, posterior_diagnostics = infer_topology_map_edges(
            rows,
            posterior_universe,
            posterior_config,
            initial_edges=posterior_initial,
            max_neighbors=2,
            seed=int(seed),
        )
        calls["infer_topology_map_edges"] += 1

    phreeqc_summary: Dict[str, object] = {}
    bayesian_ages: Dict[str, Dict[str, float]] = {}
    if run_heavy_modules:
        phreeqc_config = replace(config, phreeqc_enabled=True)
        phreeqc = run_phreeqc(rows, phreeqc_config)
        calls["run_phreeqc"] += 1
        ok = sum(bool(value.get("phreeqc_ok")) for value in phreeqc.values())
        phreeqc_summary = {
            "n_samples": len(phreeqc),
            "n_ok": int(ok),
            "success_fraction": float(ok / len(phreeqc)) if phreeqc else 0.0,
            "failure_reasons": sorted(
                {
                    str(value.get("skipped_reason"))
                    for value in phreeqc.values()
                    if not value.get("phreeqc_ok")
                }
            ),
        }
        constrained_fits = fit_network(
            rows,
            deepcopy(candidates),
            phreeqc_config,
            phreeqc_results=phreeqc,
        )
        calls["fit_network"] += 1
        constrained_by_id = {
            result.edge_id: result for result in constrained_fits
        }
        shared_fit_ids = sorted(set(fit_result_by_id) & set(constrained_by_id))
        phreeqc_summary["n_edge_fits_constrained"] = len(shared_fit_ids)
        phreeqc_summary["mean_absolute_objective_change"] = (
            float(
                sum(
                    abs(
                        constrained_by_id[edge_id].objective_score
                        - fit_result_by_id[edge_id].objective_score
                    )
                    for edge_id in shared_fit_ids
                )
                / len(shared_fit_ids)
            )
            if shared_fit_ids
            else None
        )

        graph = nx.DiGraph()
        for row in rows:
            graph.add_node(str(row["site_id"]))
        selection_key = (
            "joint_logistic_probability"
            if fusion_model is not None
            else "joint_probability"
        )
        posterior_map_ids = {
            str(edge_id)
            for edge_id in posterior_diagnostics.get("map_edges", [])
        }
        selected_ids = posterior_map_ids or {
            str(item["edge_id"])
            for item in edge_scores
            if float(item[selection_key]) >= float(threshold)
        }
        phreeqc_summary["aging_graph_source"] = (
            "topology_posterior_map"
            if posterior_map_ids
            else "logistic_threshold_fallback"
        )
        by_id = {_edge_id(edge): edge for edge in candidates}
        for edge_id in selected_ids:
            edge = by_id[edge_id]
            distance_km = float((edge.attrs or {}).get("distance_km", 0.001))
            graph.add_edge(edge.u, edge.v, length_m=max(1.0, distance_km * 1000.0))
        observations_by_node = {
            str(row["site_id"]): float(row["tritium_TU"])
            for row in rows
            if row.get("tritium_TU") is not None
        }
        sigmas = {
            node_id: max(0.10, 0.12 * value)
            for node_id, value in observations_by_node.items()
        }
        bayesian_ages = infer_network_ages_bayesian(
            graph,
            observations_by_node,
            sigmas,
            sample_date=2025.5,
            n_samples=int(aging_draws),
            n_chains=int(aging_chains),
            random_seed=int(seed),
            model_type="PFM",
        )
        calls["infer_network_ages_bayesian"] += 1
        age_feedback_rows = [dict(row) for row in rows]
        for row in age_feedback_rows:
            node_id = str(row["site_id"])
            age_result = bayesian_ages.get(node_id)
            if age_result:
                row["mean_age_years"] = float(age_result["mean_age_years"])
                row["mean_age_std_years"] = float(age_result["std_age_years"])
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Recharge date .* exceeds input history start"
            )
            feedback_edges = refine_edges_with_sheaf(
                age_feedback_rows,
                deepcopy(candidates),
                age_only_config,
            )
        calls["refine_edges_with_sheaf"] += 1
        feedback_ids = {_edge_id(edge) for edge in feedback_edges}
        pre_feedback_ids = {_edge_id(edge) for edge in age_only_edges}
        phreeqc_summary["age_feedback_edge_jaccard"] = (
            len(feedback_ids & pre_feedback_ids)
            / len(feedback_ids | pre_feedback_ids)
            if feedback_ids | pre_feedback_ids
            else 1.0
        )
        phreeqc_summary["age_feedback_changed_edge_count"] = len(
            feedback_ids ^ pre_feedback_ids
        )

    selection_key = (
        "joint_logistic_probability"
        if fusion_model is not None
        else "joint_probability"
    )
    selected = tuple(
        sorted(
            str(item["edge_id"])
            for item in edge_scores
            if float(item[selection_key]) >= float(threshold)
        )
    )
    return BlindInferenceResult(
        edge_scores=tuple(edge_scores),
        selected_edges=selected,
        candidate_edges=tuple(sorted(_edge_id(edge) for edge in candidates)),
        module_calls=calls,
        posterior_diagnostics=posterior_diagnostics,
        bayesian_ages=bayesian_ages,
        phreeqc_summary=phreeqc_summary,
    )
