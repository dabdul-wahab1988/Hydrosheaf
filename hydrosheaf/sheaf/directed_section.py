"""Directed sheaf section solver for chemistry vectors."""

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..config import Config
from ..graph.types import Edge
from ..models.mixing import fit_evaporation, fit_mixing
from ..models.reactions import build_reaction_dictionary, fit_reactions


@dataclass
class DirectedEdgeMap:
    edge: Edge
    alpha: float
    offset: List[float]
    weight: float
    objective: float
    transport_model: str
    endmember_id: Optional[str]
    residual_norm: float


def _edge_confidence(edge: Edge) -> float:
    attrs = edge.attrs or {}
    prior = attrs.get("edge_confidence", attrs.get("p_uv", 1.0))
    p_uv = 1.0
    if isinstance(prior, (int, float, str)):
        try:
            p_uv = float(prior)
        except (TypeError, ValueError):
            p_uv = 1.0
    return min(1.0, max(1e-6, p_uv))


def _reaction_vector(
    pre_residual: Sequence[float], post_residual: Sequence[float]
) -> List[float]:
    return [pre - post for pre, post in zip(pre_residual, post_residual)]


def _fit_edge_map(
    x_u: Sequence[float],
    x_v: Sequence[float],
    config: Config,
) -> Optional[Tuple[float, List[float], float, str, Optional[str]]]:
    transport_weights = getattr(config, "conservative_weights", config.weights)
    candidates: List[
        Tuple[str, Optional[str], Optional[float], Optional[float], List[float], float]
    ] = []

    if "evap" in config.transport_models_enabled:
        gamma, residual, norm = fit_evaporation(x_u, x_v, transport_weights)
        candidates.append(("evap", None, gamma, None, residual, norm))

    if "mix" in config.transport_models_enabled:
        for end_id, endmember in config.mixing_endmembers.items():
            f, residual, norm = fit_mixing(x_u, x_v, endmember, transport_weights)
            candidates.append(("mix", str(end_id), None, f, residual, norm))

    if not candidates:
        return None

    reaction_matrix, labels, _ = build_reaction_dictionary(config)
    signed_mask = [label in config.signed_reaction_labels for label in labels]
    lambda_l1 = config.lambda_l1_value()

    best: Optional[Tuple[float, List[float], float, str, Optional[str]]] = None
    best_objective = float("inf")

    for transport_model, end_id, gamma, f, residual, _ in candidates:
        reaction_fit = fit_reactions(
            list(residual),
            reaction_matrix,
            weights=list(config.weights),
            lambda_l1=lambda_l1,
            max_iter=config.reaction_max_iter,
            tol=config.reaction_tol,
            signed_mask=signed_mask,
        )
        objective = reaction_fit.residual_norm + lambda_l1 * reaction_fit.l1_norm

        if transport_model == "evap":
            alpha = float(gamma if gamma is not None else 1.0)
            offset = [0.0] * len(x_u)
            transport_pred = [alpha * float(x) for x in x_u]
        else:
            f_val = float(f if f is not None else 0.0)
            endmember = config.mixing_endmembers.get(str(end_id), [])
            if len(endmember) != len(x_u):
                continue
            alpha = 1.0 - f_val
            offset = [f_val * float(v) for v in endmember]
            transport_pred = [
                float(u) + f_val * (float(e) - float(u))
                for u, e in zip(x_u, endmember)
            ]

        reaction_pred = _reaction_vector(residual, reaction_fit.residual)
        offset = [o + r for o, r in zip(offset, reaction_pred)]

        if objective < best_objective:
            best_objective = objective
            best = (alpha, offset, objective, transport_model, end_id)

    return best


def build_edge_maps(
    edges: Iterable[Edge],
    node_values: Mapping[str, Sequence[float]],
    config: Config,
    prior_weight: float = 1.0,
) -> Dict[str, DirectedEdgeMap]:
    maps: Dict[str, DirectedEdgeMap] = {}
    for edge in edges:
        x_u = node_values.get(edge.u)
        x_v = node_values.get(edge.v)
        if x_u is None or x_v is None:
            continue
        fit = _fit_edge_map(x_u, x_v, config)
        if fit is None:
            continue
        alpha, offset, objective, transport_model, end_id = fit
        weight = prior_weight * _edge_confidence(edge)
        maps[edge.edge_id] = DirectedEdgeMap(
            edge=edge,
            alpha=alpha,
            offset=offset,
            weight=weight,
            objective=objective,
            transport_model=transport_model,
            endmember_id=end_id,
            residual_norm=objective,
        )
    return maps


def solve_directed_section(
    node_ids: Iterable[str],
    edge_maps: Iterable[DirectedEdgeMap],
    node_obs: Mapping[str, Optional[Sequence[float]]],
    obs_weight: float = 1.0,
    diag_eps: float = 1e-6,
) -> Dict[str, List[float]]:
    nodes = list(node_ids)
    if not nodes:
        return {}
    idx = {node_id: i for i, node_id in enumerate(nodes)}
    n = len(nodes)

    dim = None
    for obs in node_obs.values():
        if obs is not None:
            dim = len(obs)
            break
    if dim is None:
        for edge_map in edge_maps:
            dim = len(edge_map.offset)
            break
    if dim is None:
        return {}

    results: Dict[str, List[float]] = {}
    for node_id in nodes:
        results[node_id] = [0.0] * dim

    for d in range(dim):
        mat = np.zeros((n, n))
        vec = np.zeros(n)

        for edge_map in edge_maps:
            edge = edge_map.edge
            i = idx.get(edge.u)
            j = idx.get(edge.v)
            if i is None or j is None:
                continue
            w = float(edge_map.weight)
            if w <= 0:
                continue
            a = float(edge_map.alpha)
            b = float(edge_map.offset[d])
            mat[i, i] += w * a * a
            mat[j, j] += w
            mat[i, j] -= w * a
            mat[j, i] -= w * a
            vec[i] += -w * a * b
            vec[j] += w * b

        if obs_weight > 0:
            for node_id, obs in node_obs.items():
                if obs is None or d >= len(obs) or obs[d] is None:
                    continue
                i = idx.get(node_id)
                if i is None:
                    continue
                mat[i, i] += obs_weight
                vec[i] += obs_weight * float(obs[d])

        for i in range(n):
            mat[i, i] += diag_eps

        try:
            sol = np.linalg.solve(mat, vec)
        except np.linalg.LinAlgError:
            sol = np.linalg.lstsq(mat, vec, rcond=None)[0]

        for node_id, i in idx.items():
            results[node_id][d] = float(sol[i])

    return results


def compute_edge_section_residuals(
    edge_maps: Mapping[str, DirectedEdgeMap],
    node_values: Mapping[str, Sequence[float]],
    weights: Sequence[float],
) -> Dict[str, float]:
    residuals: Dict[str, float] = {}
    for edge_id, edge_map in edge_maps.items():
        edge = edge_map.edge
        x_u = node_values.get(edge.u)
        x_v = node_values.get(edge.v)
        if x_u is None or x_v is None:
            continue
        res_sq = 0.0
        for w, u_val, v_val, offset in zip(weights, x_u, x_v, edge_map.offset):
            diff = edge_map.alpha * float(u_val) + float(offset) - float(v_val)
            res_sq += float(w) * diff * diff
        residuals[edge_id] = res_sq ** 0.5
    return residuals
