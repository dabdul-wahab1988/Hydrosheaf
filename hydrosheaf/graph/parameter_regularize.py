"""Convex parameter-level graph regularization for groundwater age inference."""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import minimize

def parameter_graph_regularize(
    single_taus, edges, edge_weights, lambda_smooth, *, min_increment_years=0.0
):
    if lambda_smooth <= 0 or not edges:
        return dict(single_taus)
    nodes = list(single_taus.keys())
    node_to_idx = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)
    tau_raw = np.array([single_taus[node] for node in nodes])
    log_tau_raw = np.log(np.maximum(tau_raw, 1e-6))
    edge_i, edge_j, edge_w = [], [], []
    for u, v in edges:
        if u in node_to_idx and v in node_to_idx:
            edge_i.append(node_to_idx[u])
            edge_j.append(node_to_idx[v])
            edge_w.append(edge_weights.get((u, v), 1.0))
    if not edge_i:
        return dict(single_taus)
    edge_i_arr = np.array(edge_i)
    edge_j_arr = np.array(edge_j)
    weights_arr = np.array(edge_w)
    def loss_and_grad(tau_opt):
        tau_opt = np.maximum(tau_opt, 1e-6)
        log_tau_opt = np.log(tau_opt)
        log_diff = log_tau_opt - log_tau_raw
        fidelity_loss = float(np.sum(log_diff**2))
        grad = 2.0 * log_diff / tau_opt
        flow_violation = tau_opt[edge_i_arr] + min_increment_years - tau_opt[edge_j_arr]
        mask = flow_violation > 0.0
        penalty = weights_arr * (flow_violation * mask) ** 2
        flow_loss = lambda_smooth * float(np.sum(penalty))
        grad_i = 2.0 * lambda_smooth * weights_arr * flow_violation * mask
        grad_j = -grad_i
        np.add.at(grad, edge_i_arr, grad_i)
        np.add.at(grad, edge_j_arr, grad_j)
        return fidelity_loss + flow_loss, grad
    bounds = [(1e-3, 100000.0) for _ in range(n)]
    x0 = tau_raw.copy()
    res = minimize(loss_and_grad, x0, method="L-BFGS-B", jac=True, bounds=bounds,
                   options={"disp": False, "ftol": 1e-6, "maxiter": 200})
    return {node: float(res.x[i]) for i, node in enumerate(nodes)}

def _akaike_weights(aiccs):
    if not aiccs:
        return {}
    values = list(aiccs.values())
    min_aic = min(values)
    deltas = {k: v - min_aic for k, v in aiccs.items()}
    exp_terms = {k: math.exp(-0.5 * d) for k, d in deltas.items()}
    total = sum(exp_terms.values())
    if total <= 0 or not math.isfinite(total):
        return {k: 1.0 / len(aiccs) for k in aiccs}
    return {k: v / total for k, v in exp_terms.items()}

def akaike_weight_similarity(aiccs_i, aiccs_j):
    w_i = _akaike_weights(aiccs_i)
    w_j = _akaike_weights(aiccs_j)
    if not w_i or not w_j:
        return 0.0
    all_models = set(w_i) | set(w_j)
    return sum(w_i.get(m, 0.0) * w_j.get(m, 0.0) for m in all_models)

def build_aicc_edge_weights(edges, node_aiccs):
    edge_weights = {}
    for u, v in edges:
        aiccs_u = node_aiccs.get(u, {})
        aiccs_v = node_aiccs.get(v, {})
        sim = akaike_weight_similarity(aiccs_u, aiccs_v)
        edge_weights[(u, v)] = max(sim, 0.01)
    return edge_weights

def reconstruct_age_from_regularized_tau(tau, model_family, secondary_params):
    if not math.isfinite(tau):
        return float("nan")
    family = (model_family or "").upper()
    if family.startswith("BMM-"):
        return float("nan")
    return max(tau, 0.0)