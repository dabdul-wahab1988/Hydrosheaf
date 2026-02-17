"""Temporal results serialization helpers."""

from __future__ import annotations

from typing import Any, Dict, List

from ..temporal import TemporalEdgeResult


def temporal_edge_result_to_dict(result: TemporalEdgeResult) -> Dict[str, Any]:
    return {
        "edge_id": result.edge_id,
        "u": result.u,
        "v": result.v,
        "residence_time_days": result.residence_time_days,
        "residence_time_method": result.residence_time_method,
        "residence_time_uncertainty": result.residence_time_uncertainty,
        "residence_time_flags": list(result.residence_time_flags or []),
        "residence_time_details": dict(result.residence_time_details or {}),
        "transport_model": result.transport_model,
        "gamma_mean": result.gamma_mean,
        "gamma_std": result.gamma_std,
        "f_mean": result.f_mean,
        "f_std": result.f_std,
        "reaction_extents_mean": result.reaction_extents_mean,
        "reaction_extents_std": result.reaction_extents_std,
        "total_residual_norm": result.total_residual_norm,
        "per_time_residual": result.per_time_residual,
        "timestamps": [ts.isoformat() for ts in result.timestamps],
        "reaction_extents_series": result.reaction_extents_series,
    }


def temporal_edge_results_to_rows(
    results: List[TemporalEdgeResult],
) -> List[Dict[str, Any]]:
    return [temporal_edge_result_to_dict(r) for r in results]
