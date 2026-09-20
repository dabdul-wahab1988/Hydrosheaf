"""Causal dynamic graph transport and convolution forward operator.

This module provides the multi-node causal forward model for dynamic groundwater
networks with time-varying edge kernels h_{uv}(tau, t), local recharge r_v(t),
state-dependent mixing fractions pi_{uv}(t), and local recharge fractions rho_v(t):

    g_v(t) = rho_v(t) * r_v(t) + sum_{u in Pa(v)} pi_{uv}(t) * [g_u * h_{uv}](t)

where:
    [g_u * h_{uv}](t) = sum_{tau=0}^{min(t, tau_max)} h_{uv}(tau, t) * g_u(t - tau)
"""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

import numpy as np

from .dynamic_edge_kernel import DynamicEdgeKernel


def edge_id_str(source: str, target: str) -> str:
    """Format standard edge identifier string 'source->target'."""
    return f"{source}->{target}"


def parse_edge_id(edge_id: str) -> tuple[str, str]:
    """Parse 'source->target' string into (source, target) tuple."""
    parts = edge_id.split("->")
    if len(parts) != 2:
        raise ValueError(f"Invalid edge_id {edge_id!r}, expected 'source->target'")
    return parts[0], parts[1]


def topological_sort(nodes: Sequence[str], edges: Sequence[tuple[str, str]]) -> list[str]:
    """Return nodes in topological order or raise ValueError if cyclic."""
    node_set = sorted(set(nodes))
    in_degree = {n: 0 for n in node_set}
    adj: dict[str, list[str]] = {n: [] for n in node_set}

    for u, v in edges:
        if u in in_degree and v in in_degree:
            in_degree[v] += 1
            adj[u].append(v)

    queue = [n for n in node_set if in_degree[n] == 0]
    queue.sort()
    order: list[str] = []

    while queue:
        u = queue.pop(0)
        order.append(u)
        for v in sorted(adj[u]):
            in_degree[v] -= 1
            if in_degree[v] == 0:
                queue.append(v)
                queue.sort()

    if len(order) != len(node_set):
        raise ValueError(f"Graph has cycles or disconnected components without sources: {edges}")
    return order


def validate_candidate_graph(
    nodes: Sequence[str],
    edges: Sequence[tuple[str, str]],
    source_node: str = "R",
) -> tuple[bool, str]:
    """Validate a declared candidate graph without using benchmark labels.

    The regional recharge boundary is an observed source, not an ordinary
    downstream node.  A candidate graph is therefore rejected when it gives
    that boundary an incoming edge, contains a cycle, or disconnects it from
    the candidate network.  This is a structural identifiability gate, not an
    oracle scenario gate.
    """
    node_set = tuple(dict.fromkeys(nodes))
    edge_list = tuple(edges)
    if source_node not in node_set:
        return False, "regional_forcing_boundary_missing_from_candidate_nodes"
    incoming_source = [edge for edge in edge_list if edge[1] == source_node]
    outgoing_source = [edge for edge in edge_list if edge[0] == source_node]
    if incoming_source:
        return False, "regional_forcing_boundary_has_incoming_candidate_edge"
    if not outgoing_source:
        return False, "regional_forcing_boundary_has_no_candidate_outgoing_edge"
    try:
        topological_sort(node_set, edge_list)
    except ValueError as exc:
        return False, f"candidate_graph_not_acyclic:{exc}"
    return True, "candidate_graph_passed_boundary_and_acyclicity_checks"


def convolve_causal_dynamic_kernel(
    source_signal: np.ndarray,
    kernel_matrix: np.ndarray,
    lags: np.ndarray,
    decay_constant_per_time: float = 0.0,
    initial_value: Optional[float] = None,
    step_days: float = 1.0,
) -> np.ndarray:
    """Convolve a 1D source signal with a time-varying causal kernel h(tau, t).

    Parameters
    ----------
    source_signal : np.ndarray
        Source series g_u(t) of length T.
    kernel_matrix : np.ndarray
        Kernel matrix of shape (T, L) where row t is h(tau, t).
    lags : np.ndarray
        Lag grid coordinates of length L, in the same physical units as
        ``step_days`` (the project convention is days).
    decay_constant_per_time : float
        Exponential decay rate lambda. Multiplies kernel by exp(-lambda * tau).
    initial_value : Optional[float]
        Value used for times prior to t=0 (if None, uses source_signal[0]).
    step_days : float
        Sampling interval of the source/output series.  Every lag must map to
        an integer number of source steps; non-grid-aligned lags are rejected.

    Returns
    -------
    np.ndarray
        Convolved output series of length T.
    """
    source = np.asarray(source_signal, dtype=float).reshape(-1)
    lag_arr = np.asarray(lags, dtype=float).reshape(-1)
    kernel = np.asarray(kernel_matrix, dtype=float)
    if step_days <= 0.0:
        raise ValueError("step_days must be positive")
    if source.ndim != 1 or len(source) == 0:
        raise ValueError("source_signal must be a non-empty one-dimensional array")
    if lag_arr.ndim != 1 or len(lag_arr) == 0:
        raise ValueError("lags must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(lag_arr)) or np.any(lag_arr < 0.0):
        raise ValueError("lags must be finite and non-negative")
    if len(lag_arr) > 1 and np.any(np.diff(lag_arr) <= 0.0):
        raise ValueError("lags must be strictly increasing")
    lag_steps_float = lag_arr / float(step_days)
    lag_steps = np.rint(lag_steps_float).astype(int)
    if not np.allclose(lag_steps_float, lag_steps, atol=1e-7, rtol=0.0):
        raise ValueError("lags must align with the declared step_days grid")

    T = len(source)
    L = len(lag_arr)
    if kernel.shape != (T, L):
        raise ValueError(f"kernel_matrix shape {kernel.shape} != (T={T}, L={L})")
    if not np.all(np.isfinite(kernel)) or np.any(kernel < -1e-10):
        raise ValueError("kernel_matrix must be finite and non-negative")

    output = np.zeros(T, dtype=float)
    init_val = float(source[0]) if initial_value is None else float(initial_value)

    # Compute decay factor across lags if decay is non-zero
    if decay_constant_per_time > 0.0:
        decay_weights = np.exp(-decay_constant_per_time * lag_arr)
    else:
        decay_weights = None

    for t in range(T):
        k_row = kernel[t]
        if decay_weights is not None:
            k_row = k_row * decay_weights

        # For discrete integer lag grids:
        # tau is lag in time steps
        val = 0.0
        for tau_idx, lag_step in enumerate(lag_steps):
            weight = k_row[tau_idx]
            if weight <= 0.0:
                continue
            t_past = t - int(lag_step)
            src_val = source[t_past] if t_past >= 0 else init_val
            val += weight * src_val
        output[t] = val

    return output


def simulate_dynamic_graph_transport(
    nodes: Sequence[str],
    edges: Sequence[tuple[str, str]],
    time_grid: np.ndarray,
    dynamic_kernel: DynamicEdgeKernel,
    local_inputs: Mapping[str, np.ndarray],
    local_recharge_fractions: Mapping[str, np.ndarray | float],
    edge_mixing_fractions: Mapping[str, np.ndarray | float],
    decay_constant_per_time: float = 0.0,
    check_simplex: bool = True,
    simplex_tolerance: float = 1e-4,
    step_days: Optional[float] = None,
) -> dict[str, np.ndarray]:
    """Forward simulation of dynamic causal groundwater network transport.

    Parameters
    ----------
    nodes : Sequence[str]
        List of all nodes in the network.
    edges : Sequence[tuple[str, str]]
        List of candidate or true directed edges (u, v).
    time_grid : np.ndarray
        Time coordinates of length T.
    dynamic_kernel : DynamicEdgeKernel
        Kernel container storing h_{uv}(tau, t) for all edges.
    local_inputs : Mapping[str, np.ndarray]
        Mapping from node id to local recharge series r_v(t) of length T.
    local_recharge_fractions : Mapping[str, np.ndarray | float]
        Mapping from node id to local recharge fraction rho_v(t).
    edge_mixing_fractions : Mapping[str, np.ndarray | float]
        Mapping from edge_id ("u->v") to upstream mixing fraction pi_{uv}(t).
    decay_constant_per_time : float
        Tracer radioactive decay rate lambda.
    check_simplex : bool
        If True, validates rho_v(t) + sum_u pi_{uv}(t) = 1.0.
    simplex_tolerance : float
        Tolerance for simplex equality.
    step_days : Optional[float]
        Sampling interval for the dynamic-kernel lag grid.  If omitted, the
        value is read from ``dynamic_kernel.metadata['step_days']`` and then
        defaults to one day.

    Returns
    -------
    dict[str, np.ndarray]
        Simulated node signal series g_v(t) for each node v in nodes.
    """
    T = len(time_grid)
    kernel_step_days = float(
        step_days if step_days is not None else dynamic_kernel.metadata.get("step_days", 1.0)
    )
    topo_nodes = topological_sort(nodes, edges)

    # Pre-parse incoming parents per node
    parents: dict[str, list[str]] = {n: [] for n in nodes}
    for u, v in edges:
        parents[v].append(u)

    # Standardize mixing fractions to shape (T,)
    rho_dict: dict[str, np.ndarray] = {}
    for n in nodes:
        raw_rho = local_recharge_fractions.get(n, 1.0 if not parents[n] else 0.0)
        if isinstance(raw_rho, (int, float, np.floating)):
            rho_dict[n] = np.full(T, float(raw_rho), dtype=float)
        else:
            arr = np.asarray(raw_rho, dtype=float)
            if len(arr) != T:
                raise ValueError(f"local_recharge_fraction for {n} length {len(arr)} != {T}")
            rho_dict[n] = arr

    pi_dict: dict[str, np.ndarray] = {}
    for u, v in edges:
        eid = edge_id_str(u, v)
        raw_pi = edge_mixing_fractions.get(eid, 1.0 / max(len(parents[v]), 1))
        if isinstance(raw_pi, (int, float, np.floating)):
            pi_dict[eid] = np.full(T, float(raw_pi), dtype=float)
        else:
            arr = np.asarray(raw_pi, dtype=float)
            if len(arr) != T:
                raise ValueError(f"edge_mixing_fraction for {eid} length {len(arr)} != {T}")
            pi_dict[eid] = arr

    # Simplex validation
    if check_simplex:
        for v in nodes:
            total_frac = np.array(rho_dict[v], copy=True)
            for u in parents[v]:
                eid = edge_id_str(u, v)
                total_frac += pi_dict[eid]
            max_dev = float(np.max(np.abs(total_frac - 1.0)))
            if max_dev > simplex_tolerance:
                raise ValueError(
                    f"Simplex conservation violation at node {v!r}: max discrepancy is {max_dev:.2e} "
                    f"(tolerance={simplex_tolerance:.2e})."
                )

    # Simulate in topological order
    node_signals: dict[str, np.ndarray] = {}

    for node in topo_nodes:
        up_parents = parents[node]
        g_v = np.zeros(T, dtype=float)

        # 1. Local recharge contribution: rho_v(t) * r_v(t)
        if node in local_inputs:
            r_v = np.asarray(local_inputs[node], dtype=float)
            if len(r_v) != T:
                raise ValueError(f"local input for {node} has length {len(r_v)} != {T}")
            g_v += rho_dict[node] * r_v
        else:
            # If no local input supplied, verify rho is zero (unless pure source)
            if not up_parents and np.any(rho_dict[node] > 0.0):
                raise ValueError(f"Source node {node!r} has no local input series r_v(t).")

        # 2. Upstream parent contributions: sum_u pi_{uv}(t) * [g_u * h_{uv}](t)
        for u in up_parents:
            eid = edge_id_str(u, node)
            g_u = node_signals[u]
            e_idx = dynamic_kernel.edge_index(eid)
            k_mat = dynamic_kernel.values[e_idx]  # shape (T, L)

            convolved = convolve_causal_dynamic_kernel(
                source_signal=g_u,
                kernel_matrix=k_mat,
                lags=dynamic_kernel.lag_grid,
                decay_constant_per_time=decay_constant_per_time,
                step_days=kernel_step_days,
            )
            g_v += pi_dict[eid] * convolved

        node_signals[node] = g_v

    return node_signals
