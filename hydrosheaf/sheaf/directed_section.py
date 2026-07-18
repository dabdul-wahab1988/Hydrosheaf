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
    prior = attrs.get(
        "sheaf_weight",
        attrs.get(
            "prior_edge_probability",
            attrs.get("edge_confidence", attrs.get("p_uv", 1.0)),
        ),
    )
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
    pre_si_mask: Optional[Mapping[str, float]] = None,
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

    reaction_matrix, labels, mineral_mask, penalty_scales = build_reaction_dictionary(config, pre_si_mask=pre_si_mask)
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
            penalty_scales=penalty_scales,
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
    pre_si_masks: Optional[Mapping[str, Mapping[str, float]]] = None,
) -> Dict[str, DirectedEdgeMap]:
    maps: Dict[str, DirectedEdgeMap] = {}
    for edge in edges:
        x_u = node_values.get(edge.u)
        x_v = node_values.get(edge.v)
        if x_u is None or x_v is None:
            continue
        
        pre_si = None
        if pre_si_masks:
            pre_si = pre_si_masks.get(edge.u)

        fit = _fit_edge_map(x_u, x_v, config, pre_si_mask=pre_si)
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


try:
    from scipy.optimize import lsq_linear, least_squares
    from scipy.sparse import lil_matrix, csr_matrix
except ImportError:
    raise ImportError("Hydrosheaf requires scipy.optimize and scipy.sparse.")

SPECIES_CHARGES = {
    "Ca": 2.0, "Mg": 2.0, "Na": 1.0, "K": 1.0, 
    "Cl": -1.0, "SO4": -2.0, "NO3": -1.0, "HCO3": -1.0, "CO3": -2.0, 
    "F": -1.0, "Br": -1.0, "H": 1.0, "OH": -1.0
}

def solve_coupled_section(
    node_ids: Iterable[str],
    edge_maps: Iterable[DirectedEdgeMap],
    node_obs: Mapping[str, Optional[Sequence[float]]],
    species_names: List[str],
    obs_weight: float = 1.0,
    charge_balance_weight: float = 0.0,
    diag_eps: float = 1e-6,
) -> Dict[str, List[float]]:
    """
    Solve the sheaf section problem with coupled chemical constraints.
    
    Unlike solve_directed_section which solves each species independently,
    this solver minimizes a global non-linear objective that can include
    cross-species constraints like charge balance or equilibrium.
    
    Minimizes:
        || Transport Residuals ||^2 
      + || Observation Residuals ||^2
      + || Charge Balance ||^2 (if weight > 0)
      + Regularization
      
    Parameters
    ----------
    node_ids : Iterable[str]
        List of node IDs
    edge_maps : Iterable[DirectedEdgeMap]
        Edge definitions and transport parameters
    node_obs : Mapping
        Observations per node
    species_names : List[str]
        List of species corresponding to the vector indices (e.g. ['Ca', 'Cl', ...])
    obs_weight : float
        Weight for observations
    charge_balance_weight : float
        Weight for charge balance penalty
    """
    nodes = list(node_ids)
    if not nodes:
        return {}
    
    n_nodes = len(nodes)
    n_species = len(species_names)
    idx = {node_id: i for i, node_id in enumerate(nodes)}
    
    # Identify charges
    charges = np.array([SPECIES_CHARGES.get(s, 0.0) for s in species_names])
    
    # Flattened state vector size
    n_params = n_nodes * n_species
    
    # Pre-process observations
    # obs_flat: map (node_idx, species_idx) -> value
    obs_flat = {}
    for node_id, obs_vec in node_obs.items():
        if obs_vec is None: continue
        u_idx = idx[node_id]
        for s_idx, val in enumerate(obs_vec):
            if val is not None:
                obs_flat[(u_idx, s_idx)] = float(val)

    edge_maps_list = list(edge_maps)

    def residual_function(x):
        # x is shape (n_nodes * n_species)
        # Reshape for easier access: (n_nodes, n_species)
        X = x.reshape((n_nodes, n_species))
        
        residuals = []
        
        # 1. Transport / Edge Residuals
        # w * (alpha * X_u + offset - X_v)
        for em in edge_maps_list:
            u_i = idx.get(em.edge.u)
            v_i = idx.get(em.edge.v)
            if u_i is None or v_i is None: continue
            
            w_sqrt = em.weight ** 0.5
            
            # Vectorized for all species on this edge
            # X[u_i] is (n_species,)
            pred_v = em.alpha * X[u_i] + np.array(em.offset[:n_species])
            res = w_sqrt * (pred_v - X[v_i])
            residuals.extend(res)
            
        # 2. Observation Residuals
        w_obs_sqrt = obs_weight ** 0.5
        for (u_i, s_i), val in obs_flat.items():
            res = w_obs_sqrt * (X[u_i, s_i] - val)
            residuals.append(res)
            
        # 3. Regularization
        # diag_eps * x
        if diag_eps > 0:
            reg_residuals = (diag_eps ** 0.5) * x
            residuals.extend(reg_residuals)
            
        # 4. Charge Balance (Coupling)
        # Sum (z_i * c_i) at each node -> 0
        if charge_balance_weight > 0:
            w_cb = charge_balance_weight ** 0.5
            # Dot product of X rows with charges
            imbalances = X @ charges # (n_nodes,)
            residuals.extend(w_cb * imbalances)
            
        return np.array(residuals)

    # Initial guess: 
    # Use mean of observations or zeros
    x0 = np.zeros(n_params)
    for (u_i, s_i), val in obs_flat.items():
        x0[u_i * n_species + s_i] = val
        
    # Solve
    # Bounds: concentrations >= 0
    res = least_squares(
        residual_function, 
        x0, 
        bounds=(0.0, np.inf), 
        method='trf', 
        ftol=1e-4, 
        xtol=1e-4
    )
    
    # Unpack result
    X_final = res.x.reshape((n_nodes, n_species))
    
    results = {}
    for i, node_id in enumerate(nodes):
        results[node_id] = X_final[i].tolist()
        
    return results


def solve_directed_section(
    node_ids: Iterable[str],
    edge_maps: Iterable[DirectedEdgeMap],
    node_obs: Mapping[str, Optional[Sequence[float]]],
    obs_weight: float = 1.0,
    diag_eps: float = 1e-6,
    non_negative: bool = True,
) -> Dict[str, List[float]]:
    """Solve the sheaf section problem using constrained least squares.

    Minimizes: sum_e w_e * (alpha_e * x_u + offset_e - x_v)^2
             + sum_obs w_obs * (x_obs - obs_val)^2
             + diag_eps * ||x||^2
    Subject to: x >= 0 (if non_negative=True)
    """
    nodes = list(node_ids)
    if not nodes:
        return {}
    idx = {node_id: i for i, node_id in enumerate(nodes)}
    n = len(nodes)

    # Determine dimension of the chemical vector
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

    # Pre-build the graph structure of the design matrix A
    # Rows: 1 per edge, 1 per observation (potential), 1 per node (regularization)
    # Actually, observations might be sparse/missing per dimension, but we can alloc max
    # We will rebuild A per dimension if observations vary, but the edge part is constant.
    # To optimize, we'll build the constant edge part first.

    edge_maps_list = list(edge_maps)
    n_edges = len(edge_maps_list)
    n_obs_max = len(node_obs)
    n_reg = n

    # Total rows estimate
    n_rows_total = n_edges + n_obs_max + n_reg

    results: Dict[str, List[float]] = {node_id: [0.0] * dim for node_id in nodes}

    # Edges contribute to A independently of dimension (alpha is scalar)
    # Residual: sqrt(w) * (alpha * x_u - x_v) ~ -sqrt(w) * offset
    A_edges = lil_matrix((n_edges, n), dtype=float)
    edge_weights_sqrt = []
    
    valid_edges = []
    
    for row_idx, edge_map in enumerate(edge_maps_list):
        u_idx = idx.get(edge_map.edge.u)
        v_idx = idx.get(edge_map.edge.v)
        
        if u_idx is None or v_idx is None:
            continue
            
        w = float(edge_map.weight)
        if w <= 0:
            continue
            
        w_sqrt = w ** 0.5
        edge_weights_sqrt.append(w_sqrt)
        
        # Term: w_sqrt * (alpha * x_u - x_v)
        # Col u: w_sqrt * alpha
        # Col v: -w_sqrt
        A_edges[row_idx, u_idx] = w_sqrt * float(edge_map.alpha)
        A_edges[row_idx, v_idx] = -w_sqrt
        valid_edges.append((row_idx, edge_map))

    # Convert constant parts to CSR for efficiency if possible, 
    # but we need to stack with variable parts.
    # Actually, just use LIL for flexible construction then CSR.
    
    # Regularization part: sqrt(eps) * x_i ~ 0
    reg_val = diag_eps ** 0.5
    
    for d in range(dim):
        # 1. Build A and b for this dimension
        # Rows: Edges + Observations + Regularization
        
        # We need to filter observations that exist for this dimension
        active_obs = []
        if obs_weight > 0:
            for node_id, obs in node_obs.items():
                if obs is not None and d < len(obs) and obs[d] is not None:
                    active_obs.append((idx[node_id], float(obs[d])))
        
        n_active_obs = len(active_obs)
        total_rows = n_edges + n_active_obs + n_reg
        
        A = lil_matrix((total_rows, n), dtype=float)
        b = np.zeros(total_rows)
        
        # Fill Edges
        # Copying from pre-built A_edges.
        A[:n_edges, :] = A_edges
        
        # Populate b vector for edges
        current_valid_idx = 0
        for row_idx, edge_map in valid_edges:
             w_sqrt = edge_weights_sqrt[current_valid_idx]
             offset_val = float(edge_map.offset[d])
             b[row_idx] = -w_sqrt * offset_val
             current_valid_idx += 1
        
        # Observations
        obs_w_sqrt = obs_weight ** 0.5
        obs_start_row = n_edges
        for k, (node_idx, val) in enumerate(active_obs):
            row = obs_start_row + k
            A[row, node_idx] = obs_w_sqrt
            b[row] = obs_w_sqrt * val
            
        # Regularization
        reg_start_row = n_edges + n_active_obs
        for k in range(n):
            row = reg_start_row + k
            A[row, k] = reg_val
            b[row] = 0.0
            
        # Solve
        A_csr = A.tocsr()
        
        if non_negative:
            # lsq_linear handles bounds
            res = lsq_linear(A_csr, b, bounds=(0, np.inf), method='trf', tol=1e-8)
            sol = res.x
        else:
            # Unconstrained
            res = lsq_linear(A_csr, b, bounds=(-np.inf, np.inf), method='trf', tol=1e-8)
            sol = res.x
            
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
