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
            f"reaction_{label}": reaction_map.get(label)
            for label in label_set
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
                "edge_distance_km": result.edge_distance_km,
                "edge_delta_h": result.edge_delta_h,
                "edge_sigma_delta_h": result.edge_sigma_delta_h,
                "edge_source_tier": result.edge_source_tier,
                "edge_flags": result.edge_flags,
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
                **transport_columns,
                **reaction_columns,
            }
        )
    return rows
