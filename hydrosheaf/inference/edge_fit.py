"""Per-edge fitting pipeline."""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Tuple

from ..config import Config
from ..data.qc import qc_flags
from ..models.ec_tds import ec_tds_penalty
from ..models.gibbs import gibbs_evaporation_penalty, compute_gibbs_metrics
from ..models.reactions import ReactionFit, build_reaction_dictionary, fit_reactions
from ..models.mixing import fit_evaporation, fit_mixing
from ..isotopes import extract_isotopes, isotope_penalty


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

    # Uncertainty fields (populated if uncertainty_method != "none")
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


    # Optional forward-validation fields (reactive transport)
    rt_validated: bool = False
    rt_simulator: Optional[str] = None
    rt_residence_time_days: Optional[float] = None
    rt_rmse: Optional[float] = None
    rt_nse: Optional[float] = None
    rt_pbias: Optional[float] = None
    rt_thermodynamic_consistent: Optional[bool] = None

    # Optional uncertainty diagnostics (Bayesian)
    uncertainty_r_hat: Dict[str, float] = field(default_factory=dict)
    uncertainty_ess: Dict[str, float] = field(default_factory=dict)

    # Optional temporal summaries (when temporal-enabled is used)
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

    # Full uncertainty result object (for plotting posterior ridges)
    uncertainty: Optional[object] = None


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
    extra_endmembers: Optional[Dict[str, List[float]]] = None,
) -> EdgeResult:
    config.validate()

    reaction_matrix, labels, _ = build_reaction_dictionary(config)
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

    # Hierarchical Fitting Strategy (Two-Step)
    # ----------------------------------------
    # Step 1: Fit Transport (f or gamma) using Conservative Tracers ONLY.
    #         Objective: Minimize weighted residual of Cl, Br, Isotopes.
    #         Constraint: No L1 penalty (transport is physical baseline).
    #
    # Step 2: Fit Reactions (z) on the Residual.
    #         Objective: Minimize weighted residual of all ions + L1 penalty.
    #         Constraint: Transport is fixed from Step 1.

    # Common Target: y = Obs_v - Obs_u
    target_y = [v - u_val for v, u_val in zip(x_v, x_u)]

    # Weights for Step 1 (Conservative Physics)
    # We use config.conservative_weights to isolate Cl, Br, etc.
    weights_phys = list(config.conservative_weights)
    
    # Ultra Upgrade: Compositional Weighting
    # If enabled, scale weights by 1/(val + eps) to approximate log-error minimization
    active_weights = list(config.weights)
    if getattr(config, "compositional_weighting", False):
        # Weight ~ 1 / variance. If error is proportional to value (relative error), sigma ~ value.
        # So weight ~ 1 / value^2.
        # Let's use 1/value for robust relative error (approximates Gamma likelihood).
        for i, val in enumerate(x_v):
            if val > 1e-6:
                scale = 1.0 / (val * val)
                active_weights[i] *= scale
                weights_phys[i] *= scale
            else:
                # Keep original weight for trace/zero values to avoid explosion
                pass

    if "evap" in config.transport_models_enabled:
        # Step 1: Physics (Evaporation)
        # Transport Vector: Obs_u (since C_v = C_u * (1 + delta) => C_v - C_u = delta * C_u)
        transport_vec = list(x_u)
        
        # Solve for delta_gamma with NO L1 penalty and conservative weights
        phys_fit = fit_reactions(
            target_y, 
            [transport_vec], 
            weights_phys, 
            lambda_l1=0.0, # No penalty for physics
            max_iter=100,
            signed_mask=[False], # delta_gamma >= 0
            lb=[0.0],
            ub=[999.0]
        )
        
        delta_gamma = phys_fit.extents[0]
        gamma_val = 1.0 + delta_gamma
        
        # Step 2: Chemistry (Reactions)
        # Calculate residual after transport: y_chem = y_total - delta_gamma * Obs_u
        transport_contrib = [delta_gamma * val for val in transport_vec]
        chem_target = [y - t for y, t in zip(target_y, transport_contrib)]
        
        # Ultra Upgrade: Iterative Jacobian
        # If enabled, we refine the reaction matrix based on the current solution
        current_matrix = reaction_matrix
        current_labels = list(labels) # Copy
        
        n_iter = 1
        if getattr(config, "iterative_jacobian_enabled", False):
            n_iter = getattr(config, "iterative_jacobian_max_iter", 3)
            
        chem_fit = None
        
        for _ in range(n_iter):
            # Solve for reactions with L1 penalty
            lambda_l1 = config.lambda_l1_value()
            chem_fit = fit_reactions(
                chem_target,
                current_matrix,
                active_weights,
                lambda_l1,
                max_iter=config.reaction_max_iter,
                tol=config.reaction_tol,
                signed_mask=signed_mask,
                lb=lb,
                ub=ub
            )
            
            # TODO: If we had a live PHREEQC link here, we would:
            # 1. Take chem_fit.extents
            # 2. Run forward model on x_u + transport
            # 3. Calculate new local slopes (Jacobian) for the active reactions
            # 4. Update current_matrix columns
            # For now, we assume the matrix is static or pre-linearized.
            pass
        
        candidates.append({
            "type": "evap",
            "end_id": None,
            "gamma": gamma_val,
            "f": None,
            "fit": chem_fit, # Store chemical fit for residual/stats
            "z": chem_fit.extents,
            "transport_residual": phys_fit.residual_norm, # Metric of how well physics explains Cl
            "labels": current_labels
        })

    if "mix" in config.transport_models_enabled:
        # Gather all endmembers
        mix_sources = list(config.mixing_endmembers.items())
        if extra_endmembers:
            mix_sources.extend(extra_endmembers.items())

        for end_id, endmember in mix_sources:
            # Step 1: Physics (Mixing)
            # Transport Vector: Endmember - Obs_u
            transport_vec = [e - u_val for e, u_val in zip(endmember, x_u)]
            
            # Solve for f with NO L1 penalty and conservative weights
            phys_fit = fit_reactions(
                target_y,
                [transport_vec],
                weights_phys,
                lambda_l1=0.0,
                max_iter=100,
                signed_mask=[False], # f >= 0
                lb=[0.0],
                ub=[1.0]
            )
            
            f_val = phys_fit.extents[0]
            
            # Step 2: Chemistry
            # Residual: y_chem = y_total - f * (End - Obs_u)
            transport_contrib = [f_val * val for val in transport_vec]
            chem_target = [y - t for y, t in zip(target_y, transport_contrib)]
            
            # Ultra Upgrade: Iterative Jacobian (Same logic)
            current_matrix = reaction_matrix
            chem_fit = None
            n_iter = 1
            if getattr(config, "iterative_jacobian_enabled", False):
                n_iter = getattr(config, "iterative_jacobian_max_iter", 3)

            for _ in range(n_iter):
                lambda_l1 = config.lambda_l1_value()
                chem_fit = fit_reactions(
                    chem_target,
                    current_matrix,
                    active_weights,
                    lambda_l1,
                    max_iter=config.reaction_max_iter,
                    tol=config.reaction_tol,
                    signed_mask=signed_mask,
                    lb=lb,
                    ub=ub
                )
            
            candidates.append({
                "type": "mix",
                "end_id": end_id,
                "gamma": None,
                "f": f_val,
                "fit": chem_fit,
                "z": chem_fit.extents,
                "transport_residual": phys_fit.residual_norm,
                "labels": labels
            })

    best_result: Optional[EdgeResult] = None
    candidate_entries: List[Dict[str, object]] = []
    
    for cand in candidates:
        transport_model = str(cand["type"])
        chem_fit = cand["fit"] # type: ignore
        gamma_value = cand["gamma"] # type: ignore
        f_value = cand["f"] # type: ignore
        end_id = cand["end_id"] # type: ignore
        z_vals = cand["z"] # type: ignore
        
        # Reconstruct "Modeled V" for penalties
        # modeled_v = Obs_u + Transport + Reaction
        # chem_fit.residual = (Obs_v - Obs_u - Transport) - Reaction
        # Obs_v - chem_fit.residual = Obs_u + Transport + Reaction
        modeled_x_v = [obs - r for obs, r in zip(x_v, chem_fit.residual)]

        penalty = 0.0
        if obs_v is not None and config.ec_tds_penalty_enabled:
            penalty = ec_tds_penalty(modeled_x_v, obs_v, config)

        iso_penalty = 0.0
        iso_metrics: Dict[str, float] = {}
        iso_used = False
        if (
            config.isotope_enabled
            and obs_u is not None
            and obs_v is not None
            and config.lmwl_defined
        ):
            iso_u = extract_isotopes(
                obs_u, config.isotope_d18o_key, config.isotope_d2h_key
            )
            iso_v = extract_isotopes(
                obs_v, config.isotope_d18o_key, config.isotope_d2h_key
            )
            
            # Handle Endmember Isotopes for Mixing Penalty
            iso_end = None
            if transport_model == "mix" and end_id:
                # Look up endmember isotopes in config
                iso_end = config.mixing_endmembers_isotopes.get(end_id)
                # If not found there, check if they were provided in extra_endmembers_isotopes?
                # For now, just config.
                
                if iso_end is None and config.strict_input_validation:
                    import warnings
                    warnings.warn(f"Endmember '{end_id}' has no isotopic definition. Isotope mixing check bypassed.")

            if iso_u and iso_v:
                iso_raw, iso_metrics = isotope_penalty(
                    iso_u[0],
                    iso_u[1],
                    iso_v[0],
                    iso_v[1],
                    config.lmwl_a,
                    config.lmwl_b,
                    transport_model,
                    d_excess_weight=config.isotope_d_excess_weight,
                    endmember_iso=iso_end
                )
                iso_penalty = config.isotope_weight * iso_raw
                iso_used = True

        # Gibbs fallback/supplement
        gibbs_penalty_val = 0.0
        gibbs_metrics_val: Dict[str, object] = {}
        gibbs_used = False
        if config.gibbs_enabled and obs_v is not None:
            gibbs_metrics_val = compute_gibbs_metrics(obs_v)
            raw_gibbs = gibbs_evaporation_penalty(
                obs_v,
                tds_precipitation=config.gibbs_tds_precipitation,
                tds_evaporation=config.gibbs_tds_evaporation,
            )
            if transport_model == "evap":
                # Penalize evaporation if Gibbs suggests otherwise
                if iso_used:
                    # Isotope available: Gibbs supplements with lower weight
                    gibbs_penalty_val = config.gibbs_weight * 0.3 * raw_gibbs
                else:
                    # No isotope: Gibbs replaces with higher weight
                    gibbs_penalty_val = config.gibbs_weight * raw_gibbs
                gibbs_used = True

        iso_consistency_penalty = 0.0
        if (
            transport_model == "evap"
            and iso_used
            and obs_u is not None
            and obs_v is not None
        ):
            # Cross-check chloride shift vs isotopic enrichment
            cl_idx = -1
            if "Cl" in config.ion_order:
                cl_idx = config.ion_order.index("Cl")

            if cl_idx >= 0:
                cl_u = x_u[cl_idx]
                cl_v = x_v[cl_idx]
                cl_ratio = (cl_v / cl_u) if cl_u > 0 else 1.0
                if gamma_value is not None:
                    mismatch = abs(cl_ratio - gamma_value)
                    if mismatch > 0.5:
                        iso_consistency_penalty = (
                            config.isotope_consistency_weight * mismatch
                        )

        # Calculate Chemistry R2 (weighted)
        mean_v = sum(x_v) / len(x_v) if x_v else 0.0
        sst = sum(w * (v - mean_v) ** 2 for w, v in zip(config.weights, x_v))
        sse = chem_fit.residual_norm
        chem_r2 = 1.0 - (sse / sst) if sst > 1e-12 else 0.0

        objective = (
            chem_fit.residual_norm
            + lambda_l1 * chem_fit.l1_norm
            + penalty
            + iso_penalty
            + gibbs_penalty_val
            + iso_consistency_penalty
        )


        candidate_entries.append(
            {
                "transport_model": transport_model,
                "endmember_id": end_id,
                "objective_score": objective,
                "transport_residual_norm": chem_fit.residual_norm,
            }
        )
        
        # Populate ReactionFit object for result
        # Since chem_fit DOES NOT contain transport parameters, we use it directly.
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
            if ce is not None:
                try:
                    charge_error_val = float(ce)
                except (ValueError, TypeError):
                    pass
            
            sr = bounds.get("skipped_reason")
            if sr is not None:
                skipped_reason_val = str(sr)

        result = EdgeResult(
            edge_id=edge_id,
            u=u,
            v=v,
            transport_model=transport_model,
            gamma=gamma_value,
            f=f_value,
            endmember_id=end_id,
            z_extents=final_reaction_fit.extents,
            z_labels=labels,
            transport_residual_norm=float(cand.get("transport_residual", 0.0)),
            anomaly_norm=chem_fit.residual_norm,
            objective_score=objective,
            l1_norm=final_reaction_fit.l1_norm,
            reaction_iterations=final_reaction_fit.iterations,
            reaction_converged=final_reaction_fit.converged,
            ec_tds_penalty=penalty,
            qc_flags=[],
            constraints_active=constraints_active,
            si_u=si_u_dict,
            si_v=si_v_dict,

                                                        phreeqc_ok=phreeqc_ok_val,

                                                        charge_error=charge_error_val,

                                                        skipped_reason=skipped_reason_val,

                                                        isotope_penalty=iso_penalty,

                                                        isotope_metrics=iso_metrics,

                                                        isotope_used=iso_used,

                                                        gibbs_penalty=gibbs_penalty_val,

                                                        gibbs_metrics=gibbs_metrics_val,

                                                        gibbs_used=gibbs_used,

                                                        isotope_consistency_penalty=iso_consistency_penalty,

            reaction_fit=final_reaction_fit,
            residual_vector=list(final_reaction_fit.residual),
            chemistry_r2=chem_r2,
        )

        if best_result is None or result.objective_score < best_result.objective_score:
            best_result = result

    if best_result is None:
        raise RuntimeError("No transport candidates evaluated.")

    # Cast to float to satisfy mypy
    scores = [float(entry["objective_score"]) for entry in candidate_entries if entry.get("objective_score") is not None]
    if not scores:
        # Fallback if somehow scores are empty but best_result is not None (should not happen)
        scores = [0.0]
        
    min_score = min(scores)
    weights = [math.exp(-(score - min_score)) for score in scores]
    total = sum(weights) or 1.0
    transport_probs: Dict[str, float] = {}
    for entry, weight in zip(candidate_entries, weights):
        transport = str(entry["transport_model"])
        transport_probs[transport] = (
            transport_probs.get(transport, 0.0) + weight / total
        )

    qc = qc_flags(x_v, config.ion_order, config.charge_balance_limit)
    if (
        config.ec_tds_penalty_limit
        and best_result.ec_tds_penalty > config.ec_tds_penalty_limit
    ):
        qc.append("ec_tds_consistency")
    best_result.qc_flags = qc
    best_result.transport_probabilities = transport_probs
    best_result.candidate_scores = candidate_entries
    return best_result
