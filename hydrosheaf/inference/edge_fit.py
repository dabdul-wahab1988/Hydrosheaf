"""Per-edge fitting pipeline."""

import math
from dataclasses import dataclass, field, replace
from typing import Dict, List, Mapping, Optional

from ..config import Config
from ..log import get_logger
from ..data.qc import qc_flags

logger = get_logger("inference.edge_fit")
from ..models.ec_tds import ec_tds_penalty
from ..models.gibbs import gibbs_evaporation_penalty, compute_gibbs_metrics
from ..models.mixing import fit_evaporation
from ..models.reactions import ReactionFit, build_reaction_dictionary, fit_reactions
from ..isotopes import extract_isotopes, isotope_penalty
from ..physics.kinetic_limit import apply_kinetic_penalties


@dataclass
class EdgeResult:
    edge_id: str
    u: str
    v: str
    transport_model: str
    gamma: Optional[float] = None
    f: Optional[float] = None
    endmember_id: Optional[str] = None
    z_extents: List[float] = field(default_factory=list)
    z_labels: List[str] = field(default_factory=list)
    transport_residual_norm: float = 0.0
    anomaly_norm: float = 0.0
    objective_score: float = 0.0
    l1_norm: float = 0.0
    reaction_iterations: int = 0
    reaction_converged: bool = True
    ec_tds_penalty: float = 0.0
    transport_probabilities: Dict[str, float] = field(default_factory=dict)
    candidate_scores: List[Dict[str, object]] = field(default_factory=list)
    constraints_active: Dict[str, str] = field(default_factory=dict)
    si_u: Dict[str, float] = field(default_factory=dict)
    si_v: Dict[str, float] = field(default_factory=dict)
    phreeqc_ok: bool = False
    charge_error: Optional[float] = None
    skipped_reason: Optional[str] = None
    edge_confidence: Optional[float] = None
    edge_map_penalty: Optional[float] = None
    edge_map_score: Optional[float] = None
    edge_sheaf_score_local: Optional[float] = None
    edge_sheaf_score_global: Optional[float] = None
    edge_sheaf_residual: Optional[float] = None
    edge_sheaf_pi_evap: Optional[float] = None
    edge_sheaf_cost_iso: Optional[float] = None
    edge_sheaf_cost_cl: Optional[float] = None
    edge_sheaf_flags: Optional[str] = None
    # Sheaf cohomology diagnostics
    sheaf_h0_dim: Optional[int] = None
    sheaf_h1_dim: Optional[int] = None
    sheaf_obstruction_energy: Optional[float] = None
    sheaf_obstruction_leverage: Optional[float] = None
    sheaf_cycle_obstruction_max: Optional[float] = None
    sheaf_cycle_count: Optional[int] = None
    # Hydraulic Hodge diagnostics
    head_hodge_obstruction_energy: Optional[float] = None
    head_hodge_obstruction_leverage: Optional[float] = None
    head_hodge_cycle_obstruction_max: Optional[float] = None
    head_hodge_cycle_count: Optional[int] = None
    head_hodge_direction_violation_count: Optional[int] = None
    head_hodge_normalized_cost: Optional[float] = None
    head_hodge_reference_distance_km: Optional[float] = None
    head_hodge_standardized_drop: Optional[float] = None
    # Bayesian topology posterior
    posterior_edge_probability: Optional[float] = None
    posterior_edge_log_odds: Optional[float] = None
    posterior_map_selected: Optional[bool] = None
    # Legacy name retained for exported-schema compatibility. The value is the
    # sum of marginal Bernoulli edge entropies, not joint graph entropy.
    posterior_topology_entropy: Optional[float] = None
    posterior_marginal_edge_entropy: Optional[float] = None
    posterior_entropy_definition: Optional[str] = None
    posterior_n_edges_mean: Optional[float] = None
    posterior_n_edges_ci95: Optional[str] = None
    posterior_acceptance_rate: Optional[float] = None
    # Optimal transport plausibility
    ot_total_cost: Optional[float] = None
    ot_score_contribution: Optional[float] = None
    ot_balanced_cost: Optional[float] = None
    ot_creation_mass: Optional[float] = None
    ot_destruction_mass: Optional[float] = None
    ot_conservative_mismatch: Optional[float] = None
    ot_reaction_plausibility: Optional[float] = None
    # Causal discovery layer
    causal_support_score: Optional[float] = None
    causal_confounded_score: Optional[float] = None
    causal_p_value: Optional[float] = None
    causal_method: Optional[str] = None
    causal_n_observations: Optional[int] = None
    causal_status: Optional[str] = None
    edge_distance_km: Optional[float] = None
    edge_delta_h: Optional[float] = None
    edge_sigma_delta_h: Optional[float] = None
    edge_source_tier: Optional[str] = None
    edge_flags: Optional[str] = None
    edge_distance_3d_m: Optional[float] = None
    edge_horizontal_distance_m: Optional[float] = None
    edge_vertical_distance_m: Optional[float] = None
    edge_type_3d: Optional[str] = None
    edge_prob_head: Optional[float] = None
    edge_prob_distance: Optional[float] = None
    edge_prob_layer: Optional[float] = None
    edge_horizontal_gradient: Optional[float] = None
    edge_vertical_gradient: Optional[float] = None
    edge_residence_time_days: Optional[float] = None
    physics_source: Optional[str] = None
    physics_tau_mean_days: Optional[float] = None
    physics_tau_std_days: Optional[float] = None
    physics_tau_p10_days: Optional[float] = None
    physics_tau_p90_days: Optional[float] = None
    isotope_penalty: float = 0.0
    isotope_metrics: Dict[str, float] = field(default_factory=dict)
    isotope_used: bool = False
    gibbs_penalty: float = 0.0
    gibbs_metrics: Dict[str, object] = field(default_factory=dict)
    gibbs_used: bool = False
    isotope_consistency_penalty: float = 0.0
    qc_flags: List[str] = field(default_factory=list)
    nitrate_source_p_manure: Optional[float] = None
    nitrate_source_logit: Optional[float] = None
    nitrate_source_evidence: List[str] = field(default_factory=list)
    nitrate_source_gates: List[str] = field(default_factory=list)
    gamma_std: Optional[float] = None
    gamma_ci_low: Optional[float] = None
    gamma_ci_high: Optional[float] = None
    f_std: Optional[float] = None
    f_ci_low: Optional[float] = None
    f_ci_high: Optional[float] = None
    extents_std: List[float] = field(default_factory=list)
    extents_ci_low: List[float] = field(default_factory=list)
    extents_ci_high: List[float] = field(default_factory=list)
    uncertainty_method: Optional[str] = None
    reaction_fit: Optional[ReactionFit] = None
    residual_vector: List[float] = field(default_factory=list)
    chemistry_r2: Optional[float] = None
    rt_validated: bool = False
    rt_simulator: Optional[str] = None
    rt_residence_time_days: Optional[float] = None
    rt_rmse: Optional[float] = None
    rt_nse: Optional[float] = None
    rt_pbias: Optional[float] = None
    rt_thermodynamic_consistent: Optional[bool] = None
    uncertainty_r_hat: Dict[str, float] = field(default_factory=dict)
    uncertainty_ess: Dict[str, float] = field(default_factory=dict)
    temporal_residence_time_days: Optional[float] = None
    temporal_residence_time_method: Optional[str] = None
    temporal_residence_time_uncertainty: Optional[float] = None
    temporal_transport_model: Optional[str] = None
    temporal_gamma_mean: Optional[float] = None
    temporal_gamma_std: Optional[float] = None
    temporal_f_mean: Optional[float] = None
    temporal_f_std: Optional[float] = None
    temporal_reaction_extents_mean: List[float] = field(default_factory=list)
    temporal_reaction_extents_std: List[float] = field(default_factory=list)
    temporal_total_residual_norm: Optional[float] = None
    temporal_n_time_points: Optional[int] = None
    temporal_residence_time_flags: List[str] = field(default_factory=list)
    temporal_residence_time_details: Dict[str, object] = field(default_factory=dict)
    uncertainty: Optional[object] = None


def _weighted_mean(values: List[float], weights: List[float]) -> float:
    total_weight = sum(weights)
    if total_weight <= 0.0:
        return sum(values) / len(values) if values else 0.0
    return sum(w * v for w, v in zip(weights, values)) / total_weight


def _information_score(residual_norm: float, n_observations: int, n_parameters: int) -> float:
    n = max(1, int(n_observations))
    rss = max(float(residual_norm), 1e-300)
    return n * math.log(rss / n) + 2.0 * max(1, int(n_parameters))


def _bounds_for_bayesian(
    bounds: Optional[Dict[str, object]], n_reactions: int
) -> Optional[List[tuple]]:
    if not bounds:
        return None
    lb_raw = bounds.get("lb")
    ub_raw = bounds.get("ub")
    if not isinstance(lb_raw, list) or not isinstance(ub_raw, list):
        return None
    out = []
    for idx in range(n_reactions):
        lb = lb_raw[idx] if idx < len(lb_raw) else None
        ub = ub_raw[idx] if idx < len(ub_raw) else None
        out.append((lb, ub))
    return out


def _attach_uncertainty(result: EdgeResult, uncertainty: object) -> None:
    result.uncertainty = uncertainty
    method = getattr(uncertainty, "method", None)
    result.uncertainty_method = str(method) if method is not None else None
    result.gamma_std = getattr(uncertainty, "gamma_std", None)
    result.gamma_ci_low = getattr(uncertainty, "gamma_ci_low", None)
    result.gamma_ci_high = getattr(uncertainty, "gamma_ci_high", None)
    result.f_std = getattr(uncertainty, "f_std", None)
    result.f_ci_low = getattr(uncertainty, "f_ci_low", None)
    result.f_ci_high = getattr(uncertainty, "f_ci_high", None)
    result.extents_std = list(getattr(uncertainty, "extents_std", []) or [])
    result.extents_ci_low = list(getattr(uncertainty, "extents_ci_low", []) or [])
    result.extents_ci_high = list(getattr(uncertainty, "extents_ci_high", []) or [])
    result.uncertainty_r_hat = dict(getattr(uncertainty, "r_hat", {}) or {})
    result.uncertainty_ess = dict(getattr(uncertainty, "ess", {}) or {})
    if result.reaction_fit is not None:
        result.reaction_fit.extents_std = result.extents_std
        result.reaction_fit.extents_ci_low = result.extents_ci_low
        result.reaction_fit.extents_ci_high = result.extents_ci_high
        result.reaction_fit.uncertainty_result = uncertainty


def fit_edge(
    x_u: List[float],
    x_v: List[float],
    config: Config,
    edge_id: str = "",
    u: str = "",
    v: str = "",
    obs_v: Optional[Mapping[str, float]] = None,
    bounds: Optional[Dict[str, object]] = None,
    obs_u: Optional[Mapping[str, float]] = None,
    residence_time_days: Optional[float] = None,
    residence_time_std_days: Optional[float] = None,
    extra_endmembers: Optional[Dict[str, List[float]]] = None,
    pre_si_mask: Optional[Mapping[str, float]] = None,
) -> EdgeResult:
    config.validate()

    add_no3src = True
    add_denit = True
    no3src_scale = 1.0
    denit_scale = 1.0
    
    if obs_u is not None:
        no3_val = float(obs_u.get("NO3") or 0.0)
        if no3_val < 0.16:
            add_denit = False
            logger.debug("Logic Gate: Pruning 'denit' (NO3 < 0.16 mmol/L)")
        if no3_val > 0.8:
            no3src_scale = 0.1
            logger.debug("Logic Gate: Forcing 'NO3src' selection (NO3 > 0.8 mmol/L)")

    # Short residence times are treated as unfavorable for denitrification.
    if residence_time_days is not None and residence_time_days < 30.0:
        denit_scale = 10.0
        logger.debug("Topology Redox Proxy: Penalizing denitrification (short residence time < 30d)")
    reaction_matrix, labels, mineral_mask, penalty_scales = build_reaction_dictionary(
        config, pre_si_mask=pre_si_mask, sample=obs_u, dynamic_denit_scale=denit_scale
    )
    signed_mask = [label in config.signed_reaction_labels for label in labels]
    
    lb: Optional[List[float]] = None
    ub: Optional[List[float]] = None
    if bounds:
        lb_raw = bounds.get("lb")
        ub_raw = bounds.get("ub")
        if isinstance(lb_raw, list):
            lb = [float(v) if v is not None else -float("inf") for v in lb_raw]
        if isinstance(ub_raw, list):
            ub = [float(v) if v is not None else float("inf") for v in ub_raw]

    candidates: List[Dict[str, object]] = []
    target_y = [v_val - u_val for v_val, u_val in zip(x_v, x_u)]
    weights_phys = list(config.conservative_weights)
    active_weights = list(config.weights)
    
    if getattr(config, "compositional_weighting", False):
        for i, val in enumerate(x_v):
            if val > 1e-6:
                scale = 1.0 / (val * val)
                active_weights[i] *= scale
                weights_phys[i] *= scale

    if "evap" in config.transport_models_enabled:
        transport_vec = list(x_u)
        gamma_val, _, transport_residual = fit_evaporation(x_u, x_v, weights_phys)
        delta_gamma = gamma_val - 1.0
        transport_contrib = [delta_gamma * val for val in transport_vec]
        chem_target = [y - t for y, t in zip(target_y, transport_contrib)]
        
        lambda_l1 = config.lambda_l1_value()
        chem_fit = fit_reactions(
            chem_target, reaction_matrix, active_weights, lambda_l1,
            max_iter=config.reaction_max_iter, tol=config.reaction_tol,
            signed_mask=signed_mask, lb=lb, ub=ub, penalty_scales=penalty_scales,
            lambda_l2=config.lambda_l2
        )
        
        candidates.append({
            "type": "evap", "end_id": None, "gamma": gamma_val, "f": None,
            "fit": chem_fit, "z": chem_fit.extents, "transport_residual": transport_residual, "labels": labels,
            "transport_k": 1, "chem_target": chem_target
        })

    if "mix" in config.transport_models_enabled:
        mix_sources = list(config.mixing_endmembers.items())
        if extra_endmembers: mix_sources.extend(extra_endmembers.items())

        for end_id, endmember in mix_sources:
            transport_vec = [e - u_val for e, u_val in zip(endmember, x_u)]
            phys_fit = fit_reactions(target_y, [transport_vec], weights_phys, lambda_l1=0.0, signed_mask=[False], lb=[0.0], ub=[1.0])
            f_val = phys_fit.extents[0]
            transport_contrib = [f_val * val for val in transport_vec]
            chem_target = [y - t for y, t in zip(target_y, transport_contrib)]
            
            lambda_l1 = config.lambda_l1_value()
            chem_fit = fit_reactions(
                chem_target, reaction_matrix, active_weights, lambda_l1,
                max_iter=config.reaction_max_iter, tol=config.reaction_tol,
                signed_mask=signed_mask, lb=lb, ub=ub, penalty_scales=penalty_scales
            )
            
            candidates.append({
                "type": "mix", "end_id": end_id, "gamma": None, "f": f_val,
                "fit": chem_fit, "z": chem_fit.extents, "transport_residual": phys_fit.residual_norm, "labels": labels,
                "transport_k": 1, "chem_target": chem_target
            })

    best_result: Optional[EdgeResult] = None
    candidate_entries: List[Dict[str, object]] = []
    
    for cand in candidates:
        transport_model = str(cand["type"])
        chem_fit = cand["fit"] # type: ignore
        gamma_value = cand["gamma"] # type: ignore
        f_value = cand["f"] # type: ignore
        end_id = cand["end_id"] # type: ignore
        
        modeled_x_v = [obs - r for obs, r in zip(x_v, chem_fit.residual)]
        penalty = 0.0
        if obs_v is not None and config.ec_tds_penalty_enabled:
            penalty = ec_tds_penalty(modeled_x_v, obs_v, config)

        iso_penalty = 0.0
        iso_metrics: Dict[str, float] = {}
        iso_used = False
        if config.isotope_enabled and obs_u is not None and obs_v is not None and config.lmwl_defined:
            iso_u = extract_isotopes(obs_u, config.isotope_d18o_key, config.isotope_d2h_key)
            iso_v = extract_isotopes(obs_v, config.isotope_d18o_key, config.isotope_d2h_key)
            iso_end = config.mixing_endmembers_isotopes.get(end_id) if transport_model == "mix" and end_id else None
            if iso_u and iso_v:
                iso_raw, iso_metrics = isotope_penalty(iso_u[0], iso_u[1], iso_v[0], iso_v[1], config.lmwl_a, config.lmwl_b, transport_model, d_excess_weight=config.isotope_d_excess_weight, endmember_iso=iso_end)
                iso_penalty = config.isotope_weight * iso_raw
                iso_used = True

        gibbs_penalty_val = 0.0
        gibbs_metrics_val: Dict[str, object] = {}
        gibbs_used = False
        if config.gibbs_enabled and obs_v is not None:
            gibbs_metrics_val = compute_gibbs_metrics(obs_v)
            raw_gibbs = gibbs_evaporation_penalty(obs_v, tds_precipitation=config.gibbs_tds_precipitation, tds_evaporation=config.gibbs_tds_evaporation)
            gibbs_penalty_val = config.gibbs_weight * (0.3 if iso_used else 1.0) * raw_gibbs
            gibbs_used = True

        iso_consistency_penalty = 0.0
        if transport_model == "evap" and iso_used and obs_u is not None and obs_v is not None:
            cl_idx = config.ion_order.index("Cl") if "Cl" in config.ion_order else -1
            if cl_idx >= 0:
                cl_ratio = (x_v[cl_idx] / x_u[cl_idx]) if x_u[cl_idx] > 0 else 1.0
                if gamma_value is not None:
                    mismatch = abs(cl_ratio - gamma_value)
                    if mismatch > 0.5: iso_consistency_penalty = config.isotope_consistency_weight * mismatch

        # Penalize reactions that are implausibly fast for the available residence time.
        kin_penalty = apply_kinetic_penalties(chem_fit.extents, labels, residence_time_days or residence_time_days, config)

        chem_target_for_score = list(cand.get("chem_target", []))  # type: ignore[arg-type]
        mean_v = _weighted_mean(chem_target_for_score, active_weights)
        sst = sum(
            w * (v_val - mean_v) ** 2
            for w, v_val in zip(active_weights, chem_target_for_score)
        )
        chem_r2 = 1.0 - (chem_fit.residual_norm / sst) if sst > 1e-12 else 0.0

        transport_residual_norm = float(cand.get("transport_residual", 0.0))
        combined_residual_norm = chem_fit.residual_norm + transport_residual_norm
        l1_penalty = config.lambda_l1_value() * chem_fit.l1_norm
        objective = combined_residual_norm + l1_penalty + penalty + iso_penalty + gibbs_penalty_val + iso_consistency_penalty + kin_penalty

        n_params = int(cand.get("transport_k", 1)) + len(
            [extent for extent in chem_fit.extents if abs(float(extent)) > 1e-6]
        ) + 1
        information_score = _information_score(
            combined_residual_norm, len(chem_fit.residual), n_params
        )

        candidate_entries.append({
            "transport_model": transport_model,
            "endmember_id": end_id,
            "objective_score": objective,
            "information_score": information_score,
            "transport_residual_norm": transport_residual_norm,
            "chemistry_residual_norm": chem_fit.residual_norm,
            "combined_residual_norm": combined_residual_norm,
            "n_parameters": n_params,
        })
        
        final_reaction_fit = chem_fit
        constraints_active: Dict[str, str] = {}
        si_u_dict: Dict[str, float] = {}
        si_v_dict: Dict[str, float] = {}
        phreeqc_ok_val = False
        charge_error_val: Optional[float] = None
        skipped_reason_val: Optional[str] = None
        
        if bounds:
            constraints_active = {str(k): str(v) for k, v in bounds.get("constraints_active", {}).items()}
            si_u_dict = {str(k): float(v) for k, v in bounds.get("si_u", {}).items()}
            si_v_dict = {str(k): float(v) for k, v in bounds.get("si_v", {}).items()}
            phreeqc_ok_val = bool(bounds.get("phreeqc_ok", False))
            ce = bounds.get("charge_error")
            if ce is not None: charge_error_val = float(ce)
            sr = bounds.get("skipped_reason")
            if sr is not None: skipped_reason_val = str(sr)

        res_obj = EdgeResult(
            edge_id=edge_id, u=u, v=v, transport_model=transport_model, gamma=gamma_value, f=f_value, endmember_id=end_id,
            z_extents=final_reaction_fit.extents, z_labels=labels, transport_residual_norm=transport_residual_norm,
            anomaly_norm=chem_fit.residual_norm, objective_score=objective, l1_norm=final_reaction_fit.l1_norm,
            reaction_iterations=final_reaction_fit.iterations, reaction_converged=final_reaction_fit.converged,
            ec_tds_penalty=penalty, qc_flags=[], constraints_active=constraints_active, si_u=si_u_dict, si_v=si_v_dict,
            phreeqc_ok=phreeqc_ok_val, charge_error=charge_error_val, skipped_reason=skipped_reason_val,
            isotope_penalty=iso_penalty, isotope_metrics=iso_metrics, isotope_used=iso_used,
            gibbs_penalty=gibbs_penalty_val, gibbs_metrics=gibbs_metrics_val, gibbs_used=gibbs_used,
            isotope_consistency_penalty=iso_consistency_penalty, reaction_fit=final_reaction_fit,
            residual_vector=list(final_reaction_fit.residual), chemistry_r2=chem_r2,
            edge_residence_time_days=residence_time_days,
            physics_tau_mean_days=residence_time_days,
            physics_tau_std_days=residence_time_std_days
        )

        if best_result is None or res_obj.objective_score < best_result.objective_score:
            best_result = res_obj

    if best_result is None: raise RuntimeError("No transport candidates evaluated.")

    scores = [float(e["information_score"]) for e in candidate_entries if e.get("information_score") is not None]
    min_score = min(scores) if scores else 0.0
    weights = [math.exp(-0.5 * (score - min_score)) for score in scores]
    total = sum(weights) or 1.0
    transport_probs: Dict[str, float] = {}
    for entry, weight in zip(candidate_entries, weights):
        t = str(entry["transport_model"])
        transport_probs[t] = transport_probs.get(t, 0.0) + weight / total

    qc = qc_flags(x_v, config.ion_order, config.charge_balance_limit)
    if config.ec_tds_penalty_limit and best_result.ec_tds_penalty > config.ec_tds_penalty_limit:
        qc.append("ec_tds_consistency")
    best_result.qc_flags = qc
    best_result.transport_probabilities = transport_probs
    best_result.candidate_scores = candidate_entries
    uncertainty_method = str(getattr(config, "uncertainty_method", "none") or "none")
    if uncertainty_method != "none":
        uq_config = replace(config, uncertainty_method="none")
        try:
            if uncertainty_method == "bootstrap":
                from ..uncertainty.bootstrap import bootstrap_edge_fit

                uncertainty = bootstrap_edge_fit(
                    x_u,
                    x_v,
                    uq_config,
                    obs_u=obs_u,
                    obs_v=obs_v,
                    bounds=bounds,
                    n_resamples=int(getattr(config, "bootstrap_n_resamples", 1000)),
                )
                _attach_uncertainty(best_result, uncertainty)
            elif uncertainty_method == "monte_carlo":
                from ..uncertainty.propagation import monte_carlo_propagate

                uncertainty = monte_carlo_propagate(
                    x_u,
                    x_v,
                    uq_config,
                    obs_u=obs_u,
                    obs_v=obs_v,
                    bounds=bounds,
                    input_uncertainty_pct=float(getattr(config, "input_uncertainty_pct", 5.0)),
                    n_samples=int(getattr(config, "monte_carlo_n_samples", 1000)),
                )
                _attach_uncertainty(best_result, uncertainty)
            elif uncertainty_method == "bayesian":
                if best_result.transport_model != "evap":
                    logger.warning(
                        "Bayesian edge uncertainty currently supports the evaporation model; skipping edge %s.",
                        edge_id or "<unknown>",
                    )
                else:
                    from ..uncertainty.bayesian import bayesian_edge_fit

                    uncertainty = bayesian_edge_fit(
                        x_u,
                        x_v,
                        reaction_matrix,
                        labels,
                        uq_config,
                        n_samples=int(getattr(config, "bayesian_n_samples", 5000)),
                        n_chains=int(getattr(config, "bayesian_n_chains", 4)),
                        target_accept=float(getattr(config, "bayesian_target_accept", 0.95)),
                        bounds=_bounds_for_bayesian(bounds, len(labels)),
                    )
                    _attach_uncertainty(best_result, uncertainty)
        except Exception as exc:
            logger.warning(
                "Uncertainty quantification failed for edge %s with method %s: %s",
                edge_id or "<unknown>",
                uncertainty_method,
                exc,
            )
    return best_result
