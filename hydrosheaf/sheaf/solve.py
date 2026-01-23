"""Weighted sheaf smoothing solver for isotopes."""

from typing import Dict, Iterable, Mapping, Optional, Tuple

import numpy as np

from ..graph.types import Edge


def solve_isotope_field(
    node_ids: Iterable[str],
    edges: Iterable[Edge],
    node_obs: Mapping[str, Tuple[Optional[float], Optional[float]]],
    edge_weights: Mapping[str, float],
    obs_weight: float = 1.0,
    diag_eps: float = 1e-6,
) -> Dict[str, Tuple[float, float]]:
    nodes = list(node_ids)
    if not nodes:
        return {}
    idx = {node_id: i for i, node_id in enumerate(nodes)}
    n = len(nodes)

    obs_d18o = np.zeros(n)
    obs_d2h = np.zeros(n)
    mask_d18o = np.zeros(n, dtype=bool)
    mask_d2h = np.zeros(n, dtype=bool)

    for node_id, values in node_obs.items():
        if node_id not in idx:
            continue
        d18o, d2h = values
        i = idx[node_id]
        if d18o is not None:
            obs_d18o[i] = float(d18o)
            mask_d18o[i] = True
        if d2h is not None:
            obs_d2h[i] = float(d2h)
            mask_d2h[i] = True

    def _solve_dimension(obs_values: np.ndarray, mask: np.ndarray) -> np.ndarray:
        mat = np.zeros((n, n))
        vec = np.zeros(n)
        if obs_weight > 0:
            for i in range(n):
                if mask[i]:
                    mat[i, i] += obs_weight
                    vec[i] += obs_weight * obs_values[i]
        for edge in edges:
            w = float(edge_weights.get(edge.edge_id, 0.0))
            if w <= 0:
                continue
            i = idx.get(edge.u)
            j = idx.get(edge.v)
            if i is None or j is None:
                continue
            mat[i, i] += w
            mat[j, j] += w
            mat[i, j] -= w
            mat[j, i] -= w
        for i in range(n):
            mat[i, i] += diag_eps
        try:
            return np.linalg.solve(mat, vec)
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(mat, vec, rcond=None)[0]

    sol_d18o = _solve_dimension(obs_d18o, mask_d18o)
    sol_d2h = _solve_dimension(obs_d2h, mask_d2h)

    results: Dict[str, Tuple[float, float]] = {}
    for node_id, i in idx.items():
        results[node_id] = (float(sol_d18o[i]), float(sol_d2h[i]))
    return results


def compute_edge_residuals(
    edges: Iterable[Edge],
    node_values: Mapping[str, Tuple[float, float]],
    sigma_d18o: float,
    sigma_d2h: float,
) -> Dict[str, float]:
    residuals: Dict[str, float] = {}
    sig18 = max(1e-6, float(sigma_d18o))
    sig2h = max(1e-6, float(sigma_d2h))
    for edge in edges:
        u_val = node_values.get(edge.u)
        v_val = node_values.get(edge.v)
        if u_val is None or v_val is None:
            continue
        d18o_u, d2h_u = u_val
        d18o_v, d2h_v = v_val
        res = ((d18o_v - d18o_u) / sig18) ** 2 + ((d2h_v - d2h_u) / sig2h) ** 2
        residuals[edge.edge_id] = float(res) ** 0.5
    return residuals
