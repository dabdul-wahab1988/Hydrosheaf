"""Tabular output helpers."""

import json
from typing import List

from ..inference.edge_fit import EdgeResult


def edge_results_table(results: List[EdgeResult]) -> List[dict]:
    rows = []
    label_set = []
    transport_set = []
    for result in results:
        for label in result.z_labels:
            if label not in label_set:
                label_set.append(label)
        for key in result.transport_probabilities:
            if key not in transport_set:
                transport_set.append(key)
    for result in results:
        reaction_map = dict(zip(result.z_labels, result.z_extents))
        reaction_columns = {
            f"reaction_{label}": reaction_map.get(label) for label in label_set
        }
        transport_columns = {
            f"transport_prob_{key}": result.transport_probabilities.get(key)
            for key in transport_set
        }
        rows.append(
            {
                "edge_id": result.edge_id,
                "u": result.u,
                "v": result.v,
                "transport_model": result.transport_model,
                "gamma": result.gamma,
                "f": result.f,
                "endmember_id": result.endmember_id,
                "transport_residual_norm": result.transport_residual_norm,
                "anomaly_norm": result.anomaly_norm,
                "objective_score": result.objective_score,
                "chemistry_r2": result.chemistry_r2,
                "l1_norm": result.l1_norm,

                "ec_tds_penalty": result.ec_tds_penalty,
                "reaction_iterations": result.reaction_iterations,
                "reaction_converged": result.reaction_converged,
                "constraints_active": json.dumps(result.constraints_active),
                "si_u": json.dumps(result.si_u),
                "si_v": json.dumps(result.si_v),
                "phreeqc_ok": result.phreeqc_ok,
                "charge_error": result.charge_error,
                "skipped_reason": result.skipped_reason,
                "edge_confidence": result.edge_confidence,
                "edge_map_penalty": result.edge_map_penalty,
                "edge_map_score": result.edge_map_score,
                "edge_sheaf_score_local": result.edge_sheaf_score_local,
                "edge_sheaf_score_global": result.edge_sheaf_score_global,
                "edge_sheaf_residual": result.edge_sheaf_residual,
                "edge_sheaf_pi_evap": result.edge_sheaf_pi_evap,
                "edge_sheaf_cost_iso": result.edge_sheaf_cost_iso,
                "edge_sheaf_cost_cl": result.edge_sheaf_cost_cl,
                "edge_sheaf_flags": result.edge_sheaf_flags,
                # Sheaf cohomology
                "sheaf_h0_dim": result.sheaf_h0_dim,
                "sheaf_h1_dim": result.sheaf_h1_dim,
                "sheaf_obstruction_energy": result.sheaf_obstruction_energy,
                "sheaf_obstruction_leverage": result.sheaf_obstruction_leverage,
                "sheaf_cycle_obstruction_max": result.sheaf_cycle_obstruction_max,
                "sheaf_cycle_count": result.sheaf_cycle_count,
                # Hydraulic Hodge diagnostics
                "head_hodge_obstruction_energy": result.head_hodge_obstruction_energy,
                "head_hodge_obstruction_leverage": result.head_hodge_obstruction_leverage,
                "head_hodge_cycle_obstruction_max": result.head_hodge_cycle_obstruction_max,
                "head_hodge_cycle_count": result.head_hodge_cycle_count,
                "head_hodge_direction_violation_count": result.head_hodge_direction_violation_count,
                "head_hodge_normalized_cost": result.head_hodge_normalized_cost,
                "head_hodge_reference_distance_km": result.head_hodge_reference_distance_km,
                "head_hodge_standardized_drop": result.head_hodge_standardized_drop,
                # Topology posterior
                "posterior_edge_probability": result.posterior_edge_probability,
                "posterior_edge_log_odds": result.posterior_edge_log_odds,
                "posterior_map_selected": result.posterior_map_selected,
                "posterior_topology_entropy": result.posterior_topology_entropy,
                "posterior_n_edges_mean": result.posterior_n_edges_mean,
                "posterior_n_edges_ci95": result.posterior_n_edges_ci95,
                "posterior_acceptance_rate": result.posterior_acceptance_rate,
                # Optimal transport
                "ot_total_cost": result.ot_total_cost,
                "ot_score_contribution": result.ot_score_contribution,
                "ot_balanced_cost": result.ot_balanced_cost,
                "ot_creation_mass": result.ot_creation_mass,
                "ot_destruction_mass": result.ot_destruction_mass,
                "ot_conservative_mismatch": result.ot_conservative_mismatch,
                "ot_reaction_plausibility": result.ot_reaction_plausibility,
                # Causal discovery
                "causal_support_score": result.causal_support_score,
                "causal_confounded_score": result.causal_confounded_score,
                "causal_p_value": result.causal_p_value,
                "causal_method": result.causal_method,
                "causal_n_observations": result.causal_n_observations,
                "causal_status": result.causal_status,
                "edge_distance_km": result.edge_distance_km,
                "edge_delta_h": result.edge_delta_h,
                "edge_sigma_delta_h": result.edge_sigma_delta_h,
                "edge_source_tier": result.edge_source_tier,
                "edge_flags": result.edge_flags,
                "edge_distance_3d_m": result.edge_distance_3d_m,
                "edge_horizontal_distance_m": result.edge_horizontal_distance_m,
                "edge_vertical_distance_m": result.edge_vertical_distance_m,
                "edge_type_3d": result.edge_type_3d,
                "edge_prob_head": result.edge_prob_head,
                "edge_prob_distance": result.edge_prob_distance,
                "edge_prob_layer": result.edge_prob_layer,
                "edge_horizontal_gradient": result.edge_horizontal_gradient,
                "edge_vertical_gradient": result.edge_vertical_gradient,
                "edge_residence_time_days": result.edge_residence_time_days,
                "physics_source": result.physics_source,
                "physics_tau_mean_days": result.physics_tau_mean_days,
                "physics_tau_std_days": result.physics_tau_std_days,
                "physics_tau_p10_days": result.physics_tau_p10_days,
                "physics_tau_p90_days": result.physics_tau_p90_days,
                "isotope_penalty": result.isotope_penalty,
                "isotope_metrics": json.dumps(result.isotope_metrics),
                "isotope_used": result.isotope_used,
                "gibbs_penalty": result.gibbs_penalty,
                "gibbs_metrics": json.dumps(result.gibbs_metrics),
                "gibbs_used": result.gibbs_used,
                "isotope_consistency_penalty": result.isotope_consistency_penalty,
                "qc_flags": ",".join(result.qc_flags),
                "nitrate_source_p_manure": result.nitrate_source_p_manure,
                "nitrate_source_logit": result.nitrate_source_logit,
                "nitrate_source_evidence": ",".join(result.nitrate_source_evidence),
                "nitrate_source_gates": ",".join(result.nitrate_source_gates),
                "uncertainty_method": result.uncertainty_method,
                "gamma_std": result.gamma_std,
                "gamma_ci_low": result.gamma_ci_low,
                "gamma_ci_high": result.gamma_ci_high,
                "f_std": result.f_std,
                "f_ci_low": result.f_ci_low,
                "f_ci_high": result.f_ci_high,
                "extents_std": json.dumps(result.extents_std),
                "extents_ci_low": json.dumps(result.extents_ci_low),
                "extents_ci_high": json.dumps(result.extents_ci_high),
                "uncertainty_r_hat": json.dumps(result.uncertainty_r_hat),
                "uncertainty_ess": json.dumps(result.uncertainty_ess),
                "temporal_residence_time_days": result.temporal_residence_time_days,
                "temporal_residence_time_method": result.temporal_residence_time_method,
                "temporal_residence_time_uncertainty": result.temporal_residence_time_uncertainty,
                "temporal_transport_model": result.temporal_transport_model,
                "temporal_gamma_mean": result.temporal_gamma_mean,
                "temporal_gamma_std": result.temporal_gamma_std,
                "temporal_f_mean": result.temporal_f_mean,
                "temporal_f_std": result.temporal_f_std,
                "temporal_reaction_extents_mean": json.dumps(
                    result.temporal_reaction_extents_mean
                ),
                "temporal_reaction_extents_std": json.dumps(
                    result.temporal_reaction_extents_std
                ),
                "temporal_total_residual_norm": result.temporal_total_residual_norm,
                "temporal_n_time_points": result.temporal_n_time_points,
                "temporal_residence_time_flags": ",".join(
                    result.temporal_residence_time_flags or []
                ),
                "temporal_residence_time_details": json.dumps(
                    result.temporal_residence_time_details or {}
                ),
                "rt_validated": result.rt_validated,
                "rt_simulator": result.rt_simulator,
                "rt_residence_time_days": result.rt_residence_time_days,
                "rt_rmse": result.rt_rmse,
                "rt_nse": result.rt_nse,
                "rt_pbias": result.rt_pbias,
                "rt_thermodynamic_consistent": result.rt_thermodynamic_consistent,
                **transport_columns,
                **reaction_columns,
            }
        )
    return rows
