"""Network-level fitting pipeline."""

import math
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

from ..config import Config
from ..log import get_logger
from dataclasses import replace

logger = get_logger("inference.network_fit")


from ..data.schema import normalize_sample, vector_from_sample
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
from ..uncertainty.bootstrap import bootstrap_edge_fit
from ..uncertainty.propagation import monte_carlo_propagate
from ..reactive_transport.validation import validate_network_forward


def _safe_float(value: object) -> Optional[float]:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def estimate_edge_residence_time_days(
    edge_attrs: Mapping[str, object], config: Config
) -> Optional[float]:
    """Estimate residence time from edge geometry and hydraulic gradient.

    Uses Darcy-based travel time approximation:
      tau = L * porosity / (K * i)
    where i is hydraulic gradient and L is path length.
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
            return float(explicit)

    k_m_per_day = float(getattr(config, "residence_time_hydraulic_k", 0.0) or 0.0)
    porosity = float(getattr(config, "residence_time_porosity", 0.0) or 0.0)
    if k_m_per_day <= 0 or not (0.0 < porosity <= 1.0):
        return None

    length_m = _safe_float(edge_attrs.get("distance_3d_m"))
    if length_m is None:
        dist_km = _safe_float(edge_attrs.get("distance_km"))
        if dist_km is not None:
            length_m = dist_km * 1000.0
    if length_m is None or length_m <= 0:
        return None

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
        return None

    tau_days = (length_m * porosity) / (k_m_per_day * gradient)
    if tau_days <= 0 or not (tau_days < 1e12):
        return None
    return float(tau_days)


def _sample_map(samples: object) -> Dict[str, Mapping[str, object]]:
    if isinstance(samples, Mapping):
        # We need to ensure the values are Mapping[str, object]
        return {str(k): v for k, v in samples.items()}
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, object]] = {}
        for row in samples:
            site_id = row.get("site_id")
            if site_id is None:
                raise ValueError("Each sample row must include site_id.")
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def fit_network(
    samples: object,
    edges: Iterable[EdgeInput],
    config: Config,
    phreeqc_results: Optional[Mapping[str, Mapping[str, object]]] = None,
    residence_time_overrides: Optional[Mapping[str, float]] = None,
) -> List[EdgeResult]:
    sample_map = _sample_map(samples)
    built_edges = build_edges(edges)
    
    logger.info(f"Fitting network with {len(built_edges)} edges and {len(sample_map)} samples.")

    if config.phreeqc_enabled and phreeqc_results is None:
        logger.info("Running global PHREEQC speciation...")
        phreeqc_results = run_phreeqc(sample_map.values(), config)


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

        # Apply layer-specific mineral dictionary if configured (default uses Config.active_minerals).
        config_edge = config
        layer_key = getattr(config, "layer_key", "aquifer_layer")
        layer_value = sample_v.get(layer_key)
        try:
            val = layer_value
            if val is None or val == "":
                layer_idx = None
            else:
                layer_idx = int(val)
        except (TypeError, ValueError):
            layer_idx = None
        if layer_idx is not None:
            minerals_for_layer = getattr(config, "layer_mineral_map", {}).get(layer_idx)
            if minerals_for_layer:
                config_edge = replace(config, active_minerals=list(minerals_for_layer))

        reaction_matrix, reaction_labels, mineral_mask = build_reaction_dictionary(
            config_edge
        )

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
        if residence_time_overrides is not None:
            tau_override = residence_time_overrides.get(edge.edge_id)
            try:
                tau_edge = float(tau_override) if tau_override is not None else None
            except (TypeError, ValueError):
                tau_edge = None
        if tau_edge is None:
            tau_edge = estimate_edge_residence_time_days(edge.attrs or {}, config_edge)

        # Apply redox overrides
        redox_overrides = get_redox_constraints(sample_v_norm, reaction_labels)
        if redox_overrides:
            # Type-safe list modification
            lb_list = edge_bounds["lb"]
            ub_list = edge_bounds["ub"]
            assert isinstance(lb_list, list)
            assert isinstance(ub_list, list)

            for i, label in enumerate(reaction_labels):
                if label in redox_overrides:
                    lb_val, ub_val = redox_overrides[label]
                    lb_list[i] = lb_val
                    ub_list[i] = ub_val

            # Type-safe dict modification
            ca_dict = edge_bounds["constraints_active"]
            assert isinstance(ca_dict, dict)
            ca_dict["redox"] = "active"

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
        )
        
        # Log significant findings (Science Level)
        if result.transport_model != "mix": # Just an example filter
             logger.science(f"Edge {edge.edge_id}: Selected model '{result.transport_model}' (obj={result.objective_score:.4f})")

        result.edge_residence_time_days = tau_edge


        # Optional uncertainty quantification (per-edge)
        if getattr(config_edge, "uncertainty_method", "none") != "none":
            method = getattr(config_edge, "uncertainty_method", "none")
            seed = getattr(config_edge, "uncertainty_seed", None)
            try:
                if method == "bootstrap":
                    uq = bootstrap_edge_fit(
                        x_u,
                        x_v,
                        config_edge,
                        obs_u=sample_u_norm,
                        obs_v=sample_v_norm,
                        bounds=edge_bounds,
                        n_resamples=getattr(config_edge, "bootstrap_n_resamples", 1000),
                        random_state=seed,
                    )
                elif method == "monte_carlo":
                    uq = monte_carlo_propagate(
                        x_u,
                        x_v,
                        config_edge,
                        obs_u=sample_u_norm,
                        obs_v=sample_v_norm,
                        bounds=edge_bounds,
                        input_uncertainty_pct=getattr(
                            config_edge, "input_uncertainty_pct", 5.0
                        ),
                        n_samples=getattr(config_edge, "monte_carlo_n_samples", 1000),
                        random_state=seed,
                    )
                elif method == "bayesian":
                    # Only implemented for evaporation gamma in bayesian_edge_fit.
                    from ..uncertainty.bayesian import bayesian_edge_fit

                    if result.transport_model != "evap":
                        uq = None
                    else:
                        bounds_list = None
                        if edge_bounds:
                            lb = edge_bounds.get("lb")
                            ub = edge_bounds.get("ub")
                        if isinstance(lb, list) and isinstance(ub, list):
                            bounds_list = []
                            for lb_v, ub_v in zip(lb, ub):
                                l_val = float(lb_v) if lb_v is not None else None
                                u_val = float(ub_v) if ub_v is not None else None
                                bounds_list.append((l_val, u_val))
                    uq = bayesian_edge_fit(
                        x_u,
                        x_v,
                        reaction_matrix,
                        reaction_labels,
                        config_edge,
                        n_samples=getattr(config_edge, "bayesian_n_samples", 5000),
                        n_chains=getattr(config_edge, "bayesian_n_chains", 4),
                        target_accept=getattr(
                            config_edge, "bayesian_target_accept", 0.95
                        ),
                        bounds=bounds_list,
                    )
                else:
                    uq = None

                if uq is not None:
                    result.uncertainty = uq
                    result.uncertainty_method = uq.method

                    result.gamma_std = uq.gamma_std
                    result.gamma_ci_low = uq.gamma_ci_low
                    result.gamma_ci_high = uq.gamma_ci_high
                    result.f_std = uq.f_std
                    result.f_ci_low = uq.f_ci_low
                    result.f_ci_high = uq.f_ci_high
                    result.extents_std = uq.extents_std
                    result.extents_ci_low = uq.extents_ci_low
                    result.extents_ci_high = uq.extents_ci_high
                    result.uncertainty_r_hat = uq.r_hat or {}
                    result.uncertainty_ess = uq.ess or {}
            except Exception:
                # Keep baseline result if UQ fails (missing PyMC, numerical issues, etc.)
                pass
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

        result.edge_confidence = _get_float("edge_confidence") or _get_float("p_uv")
        result.edge_map_penalty = _get_float("edge_map_penalty")
        result.edge_map_score = _get_float("edge_map_score")
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
        result.physics_tau_mean_days = _get_float("physics_tau_mean_days")
        result.physics_tau_std_days = _get_float("physics_tau_std_days")
        result.physics_tau_p10_days = _get_float("physics_tau_p10_days")
        result.physics_tau_p90_days = _get_float("physics_tau_p90_days")
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

    # Optional forward reactive transport validation (per-edge annotations)
    if getattr(config, "reactive_transport_validation", False):
        try:
            # Ensure numeric samples for validators (CSV inputs may be strings)
            numeric_samples: Dict[str, Dict[str, float]] = {}
            for site_id, sample in sample_map.items():
                numeric_samples[site_id] = normalize_sample(
                    sample, config.ion_order, config.detection_limit_policy
                )
            summary = validate_network_forward(results, numeric_samples, config)
            for res in results:
                rt = summary.edge_results.get(res.edge_id)
                if rt is None:
                    continue
                res.rt_validated = True
                res.rt_simulator = rt.simulator
                res.rt_residence_time_days = rt.inverse_residence_time_days
                res.rt_rmse = rt.rmse
                res.rt_nse = rt.nse
                res.rt_pbias = rt.pbias
                res.rt_thermodynamic_consistent = rt.thermodynamic_consistent
        except Exception:
            # Keep baseline results if validation fails (missing PHREEQC kinetic, etc.)
            pass

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
            # Infer edges using 3D graph construction; convert Edge3D -> Edge.
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
            head_inference=getattr(
                config_for_edges, "edge_head_inference", "heuristic"
            ),
            dtw_prior_mu=getattr(config_for_edges, "edge_dtw_prior_mu", 5.0),
            dtw_prior_sigma=getattr(config_for_edges, "edge_dtw_prior_sigma", 5.0),
            head_prior_mu=getattr(config_for_edges, "edge_head_prior_mu", 0.0),
            head_prior_sigma=getattr(config_for_edges, "edge_head_prior_sigma", 1000.0),
            mcmc_draws=getattr(config_for_edges, "bayesian_n_samples", 1000),
            mcmc_chains=getattr(config_for_edges, "bayesian_n_chains", 2),
            mcmc_target_accept=getattr(config_for_edges, "bayesian_target_accept", 0.9),
            mcmc_warmup_fraction=getattr(
                config_for_edges, "bayesian_warmup_fraction", 0.5
            ),
        )
        return _apply_overrides(inferred)

    if method == "probabilistic_map":
        config.validate()
        prior_weight = float(getattr(config, "edge_map_prior_weight", 0.0) or 0.0)
        if prior_weight <= 0:
            raise ValueError(
                "--infer-edges-method probabilistic_map requires edge_map_prior_weight > 0."
            )

        candidate_multiplier = int(
            getattr(config, "edge_map_candidate_multiplier", 5) or 5
        )
        candidate_config = replace(
            config,
            edge_p_min=float(getattr(config, "edge_map_p_min", 0.1)),
            edge_max_neighbors=int(config.edge_max_neighbors) * candidate_multiplier,
        )
        candidate_config.validate()
        candidate_edges = infer_probabilistic_edges_from_config(candidate_config)
        if not candidate_edges:
            return []

        # Score candidate edges using chemistry objective, then select top-k per upstream node with MAP weighting.
        scoring_config = replace(
            config,
            uncertainty_method="none",
            reactive_transport_validation=False,
            nitrate_source_enabled=False,
        )
        scoring_config.validate()
        candidate_results = fit_network(
            samples_iter, candidate_edges, scoring_config, phreeqc_results=None
        )
        result_by_edge: Dict[str, float] = {
            res.edge_id: float(res.objective_score) for res in candidate_results
        }

        eps = 1e-12
        scored_by_u: Dict[str, List[tuple]] = {}
        for edge in candidate_edges:
            if edge.edge_id not in result_by_edge:
                continue
            chem_score = result_by_edge[edge.edge_id]
            attrs = edge.attrs or {}
            prior_p = attrs.get("edge_confidence", attrs.get("p_uv"))
            try:
                p_uv = float(prior_p) if prior_p is not None else 1.0
            except (TypeError, ValueError):
                p_uv = 1.0
            p_uv = min(1.0, max(eps, p_uv))
            map_penalty = -math.log(p_uv)
            map_score = chem_score + prior_weight * map_penalty
            attrs["edge_map_penalty"] = map_penalty
            attrs["edge_map_score"] = map_score
            edge.attrs = attrs
            scored_by_u.setdefault(edge.u, []).append((map_score, edge))

        selected: List[Edge] = []
        for _, scored in scored_by_u.items():
            scored.sort(key=lambda item: item[0])
            for _, edge in scored[: int(config.edge_max_neighbors)]:
                selected.append(edge)
        return selected

    if method == "probabilistic_sheaf":
        config.validate()
        candidate_multiplier = int(
            getattr(config, "edge_map_candidate_multiplier", 5) or 5
        )
        candidate_config = replace(
            config,
            edge_p_min=float(getattr(config, "edge_map_p_min", 0.1)),
            edge_max_neighbors=int(config.edge_max_neighbors) * candidate_multiplier,
        )
        candidate_config.validate()
        candidate_edges = infer_probabilistic_edges_from_config(candidate_config)
        if not candidate_edges:
            return []

        refined = refine_edges_with_sheaf(samples_iter, candidate_edges, config)
        return _apply_overrides(refined)

    # Default probabilistic inference (2D or 3D).
    return infer_probabilistic_edges_from_config(config)


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
        transport_counts[result.transport_model] = (
            transport_counts.get(result.transport_model, 0) + 1
        )
        for key, value in result.transport_probabilities.items():
            transport_prob_sums[key] = transport_prob_sums.get(key, 0.0) + value
        for label, value in zip(result.z_labels, result.z_extents):
            reaction_sums[label] = reaction_sums.get(label, 0.0) + value
            reaction_counts[label] = reaction_counts.get(label, 0) + 1
            if abs(value) > 1e-9:
                reaction_nonzero[label] = reaction_nonzero.get(label, 0) + 1
            reaction_intensity_sums[label] = reaction_intensity_sums.get(
                label, 0.0
            ) + abs(value)

    reaction_means = {
        label: reaction_sums[label] / reaction_counts[label] for label in reaction_sums
    }
    transport_probabilities_mean = {
        key: transport_prob_sums[key] / len(results) for key in transport_prob_sums
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
