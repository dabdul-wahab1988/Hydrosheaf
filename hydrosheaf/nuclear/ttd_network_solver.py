"""Regularized convex quadratic programming solver for non-parametric network TTD inversion.

Solves the coupled network objective:
    min_{g_1, ..., g_N} sum_i 1/2 || (A_i g_i - c_i) / sigma_i ||^2
                        + lambda_smooth sum_i || D_2 g_i ||^2
                        + lambda_graph sum_v || g_v - (rho_v r_v + sum_u pi_uv T_uv g_u) ||^2
subject to g_ik >= 0 and sum_k g_ik = 1.0 for all nodes i.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import time
from typing import Any, Dict, Mapping, Optional, Tuple

import networkx as nx
import numpy as np
from scipy.optimize import minimize

from .ttd_diagnostics import (
    CODE_MISSING_FORWARD_SYSTEM,
    CODE_MISSING_MIXING_SPEC,
    CODE_MISSING_TRANSPORT_OPERATOR,
    CODE_NODE_ID_COLLISION,
    CODE_OPTIMIZATION_FAILURE,
    audit_graph_topology,
    audit_node_physical_evidence,
)
from .ttd_grid import TTDGrid, build_mass_aware_curvature_matrix, wasserstein_1d
from .ttd_kernel_builder import MultiTracerForwardSystem
from .ttd_transport import (
    EdgeTransportOperator,
    NodeMixingSpecification,
    build_local_recharge_distribution,
)


@dataclass(frozen=True)
class SingleNodeTTDResult:
    """TTD inversion result for an independent single node (no network coupling)."""

    node_id: str
    status: str  # "ESTIMATED" or "ABSTAIN"
    g: np.ndarray
    young_water_fraction: float
    mean_transit_time_years: float
    chi_squared: float
    reduced_chi_squared: float
    shannon_entropy: float
    chi_squared_per_observation: float = float("nan")
    effective_degrees_of_freedom: float = float("nan")
    abstention_reasons: Tuple[str, ...] = ()
    diagnostics: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class NetworkTTDResult:
    """Joint network TTD inversion result across all nodes."""

    status: str  # "ESTIMATED" or "ABSTAIN"
    node_results: Mapping[str, SingleNodeTTDResult]
    edge_discrepancies_w1: Mapping[Tuple[str, str], float]
    edge_discrepancies_l2: Mapping[Tuple[str, str], float]
    network_sheaf_energy: float
    abstention_reasons: Tuple[str, ...] = ()
    diagnostics: Mapping[str, Any] = field(default_factory=dict)

    def node_g(self, node_id: str) -> Optional[np.ndarray]:
        res = self.node_results.get(node_id)
        return res.g if res is not None else None

    def young_water_fractions(self) -> Dict[str, float]:
        return {nid: res.young_water_fraction for nid, res in self.node_results.items()}

    def mean_transit_times(self) -> Dict[str, float]:
        return {nid: res.mean_transit_time_years for nid, res in self.node_results.items()}


def _abstained_node_result(
    node_id: str,
    grid: TTDGrid,
    reasons: Tuple[str, ...],
    *,
    diagnostics: Optional[Mapping[str, Any]] = None,
) -> SingleNodeTTDResult:
    """Create a consistently shaped abstention record for a node."""
    return SingleNodeTTDResult(
        node_id=str(node_id),
        status="ABSTAIN",
        g=np.zeros(grid.n_bins),
        young_water_fraction=float("nan"),
        mean_transit_time_years=float("nan"),
        chi_squared=float("nan"),
        reduced_chi_squared=float("nan"),
        shannon_entropy=float("nan"),
        chi_squared_per_observation=float("nan"),
        effective_degrees_of_freedom=float("nan"),
        abstention_reasons=tuple(reasons),
        diagnostics=dict(diagnostics or {}),
    )


def _optimizer_is_valid(
    opt: Any,
    x: np.ndarray,
    n_nodes: int,
    k_bins: int,
) -> Tuple[bool, Dict[str, float]]:
    """Validate optimizer success, finiteness, bounds, and simplex constraints."""
    metrics = {
        "optimizer_success": float(bool(getattr(opt, "success", False))),
        "max_simplex_error": float("inf"),
        "max_bound_violation": float("inf"),
    }
    if not bool(getattr(opt, "success", False)):
        return False, metrics
    if x.shape != (n_nodes * k_bins,) or not np.all(np.isfinite(x)):
        return False, metrics

    max_simplex_error = 0.0
    max_bound_violation = 0.0
    for i in range(n_nodes):
        block = x[i * k_bins : (i + 1) * k_bins]
        max_simplex_error = max(max_simplex_error, abs(float(np.sum(block) - 1.0)))
        max_bound_violation = max(max_bound_violation, float(max(0.0, -np.min(block))))
        max_bound_violation = max(max_bound_violation, float(max(0.0, np.max(block) - 1.0)))

    metrics["max_simplex_error"] = max_simplex_error
    metrics["max_bound_violation"] = max_bound_violation
    return max_simplex_error <= 1e-6 and max_bound_violation <= 1e-7, metrics


def _abstained_network_result(
    graph: nx.DiGraph,
    grid: TTDGrid,
    reasons: Tuple[str, ...],
    *,
    diagnostics: Optional[Mapping[str, Any]] = None,
) -> NetworkTTDResult:
    """Create a consistently shaped network-level abstention result."""
    node_results = {
        str(node_id): _abstained_node_result(str(node_id), grid, reasons)
        for node_id in graph.nodes()
    }
    return NetworkTTDResult(
        status="ABSTAIN",
        node_results=node_results,
        edge_discrepancies_w1={},
        edge_discrepancies_l2={},
        network_sheaf_energy=float("nan"),
        abstention_reasons=tuple(reasons),
        diagnostics=dict(diagnostics or {}),
    )


def solve_single_node_ttd(
    system: MultiTracerForwardSystem,
    grid: TTDGrid,
    *,
    lambda_smoothness: float = 0.05,
    max_condition_number: float = 1e8,
    d2_matrix: Optional[np.ndarray] = None,
) -> SingleNodeTTDResult:
    """Solve a regularized non-parametric TTD for one node.

    The result is a regularized point estimate.  ``reduced_chi_squared`` is
    retained for API compatibility but is intentionally ``NaN`` because a
    non-parametric simplex inverse does not have a conventional residual
    degrees-of-freedom calculation.  Use ``chi_squared_per_observation`` for a
    descriptive normalized misfit instead.
    """
    if not np.isfinite(float(lambda_smoothness)) or lambda_smoothness < 0.0:
        return _abstained_node_result(
            system.node_id,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={"message": "lambda_smoothness must be finite and non-negative"},
        )

    gate = audit_node_physical_evidence(system, max_condition_number=max_condition_number)
    if not gate.can_proceed:
        return _abstained_node_result(
            system.node_id,
            grid,
            gate.reason_codes,
            diagnostics={"gate_messages": gate.messages, "gate_metrics": gate.metrics},
        )

    k_bins = grid.n_bins
    d2 = d2_matrix if d2_matrix is not None else build_mass_aware_curvature_matrix(grid)
    d2 = np.asarray(d2, dtype=float)
    if d2.ndim != 2 or d2.shape[1] != k_bins or not np.all(np.isfinite(d2)):
        return _abstained_node_result(
            system.node_id,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={"message": "d2_matrix has an invalid shape or non-finite values"},
        )

    # Data quadratic form: 1/2 g^T (A^T W A) g - (c^T W A) g
    a_mat = system.matrix
    c_obs = system.observations
    sigmas = system.sigmas
    weights = system.weights

    inv_var = weights / (sigmas ** 2)
    h_data = a_mat.T @ (inv_var[:, None] * a_mat)
    f_data = -a_mat.T @ (inv_var * c_obs)

    # Smoothness quadratic form on the density-equivalent representation of
    # the mass vector, including non-uniform-grid quadrature weighting.
    h_smooth = 2.0 * float(lambda_smoothness) * (d2.T @ d2)
    h_total = h_data + h_smooth

    # Regularized objective
    def objective(g: np.ndarray) -> float:
        return float(0.5 * g @ (h_total @ g) + f_data @ g)

    def jacobian(g: np.ndarray) -> np.ndarray:
        return h_total @ g + f_data

    # Initial guess: uniform distribution
    g0 = np.ones(k_bins) / k_bins
    bounds = [(0.0, 1.0) for _ in range(k_bins)]
    constraints = [{"type": "eq", "fun": lambda g: np.sum(g) - 1.0}]

    opt = minimize(
        objective,
        g0,
        jac=jacobian,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 500, "ftol": 1e-9},
    )

    valid, optimizer_metrics = _optimizer_is_valid(opt, np.asarray(opt.x), 1, k_bins)
    if not valid:
        return _abstained_node_result(
            system.node_id,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={
                "iterations": getattr(opt, "nit", None),
                "success": bool(getattr(opt, "success", False)),
                "message": str(getattr(opt, "message", "")),
                **optimizer_metrics,
            },
        )

    g_opt = np.maximum(0.0, opt.x)
    g_sum = float(np.sum(g_opt))
    if g_sum <= 0.0 or not np.isfinite(g_sum):
        return _abstained_node_result(
            system.node_id,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={"message": "optimizer returned zero or non-finite simplex mass"},
        )
    g_opt = g_opt / g_sum

    chi2 = system.chi_squared(g_opt)
    chi2_per_observation = chi2 / max(1, system.n_tracers)

    yf = grid.young_water_fraction(g_opt)
    mtt = grid.mean_transit_time(g_opt)
    ent = float(-np.sum(g_opt * np.log(np.maximum(g_opt, 1e-12))))

    return SingleNodeTTDResult(
        node_id=system.node_id,
        status="ESTIMATED",
        g=g_opt,
        young_water_fraction=yf,
        mean_transit_time_years=mtt,
        chi_squared=chi2,
        reduced_chi_squared=float("nan"),
        shannon_entropy=ent,
        chi_squared_per_observation=chi2_per_observation,
        effective_degrees_of_freedom=float("nan"),
        diagnostics={
            "iterations": opt.nit,
            "success": opt.success,
            "message": opt.message,
            **optimizer_metrics,
            "regularization_semantics": "mass-aware curvature on density-equivalent vector",
            "identifiability": "regularized point estimate; full TTD is not identified by tracer count alone",
        },
    )


def solve_network_ttd(
    graph: nx.DiGraph,
    forward_systems: Mapping[str, MultiTracerForwardSystem],
    transport_operators: Mapping[Tuple[str, str], EdgeTransportOperator],
    grid: TTDGrid,
    *,
    mixing_specs: Optional[Mapping[str, NodeMixingSpecification]] = None,
    lambda_smoothness: float = 0.05,
    lambda_graph: float = 1.0,
    max_condition_number: float = 1e8,
    node_heads: Optional[Mapping[str, float]] = None,
    allow_default_mixing: bool = False,
) -> NetworkTTDResult:
    """Solve the joint regularized network TTD inverse problem across a DAG.

    Parameters
    ----------
    graph : nx.DiGraph
        Directed groundwater flow graph.
    forward_systems : Mapping[str, MultiTracerForwardSystem]
        Multi-tracer linear systems for nodes.
    transport_operators : Mapping[Tuple[str, str], EdgeTransportOperator]
        Edge transit operators T_{uv}.
    grid : TTDGrid
        Discrete age grid.
    mixing_specs : Optional[Mapping[str, NodeMixingSpecification]]
        Pre-declared node mixing models.  Downstream nodes must have an
        explicit specification unless ``allow_default_mixing`` is enabled.
    lambda_smoothness : float
        Weight for within-node curvature regularization.
    lambda_graph : float
        Weight for graph transport continuity regularization.
    max_condition_number : float
        Condition number threshold for numerical gating.
    node_heads : Optional[Mapping[str, float]]
        Hydraulic head measurements for topological gradient gating.
    allow_default_mixing : bool
        Compatibility escape hatch for the historical 20% local / 80%
        upstream prior.  The default is ``False`` because that prior is an
        assumption, not information supplied by the observations.
    """
    t_start = time.perf_counter()

    if (
        not np.isfinite(float(lambda_smoothness))
        or float(lambda_smoothness) < 0.0
        or not np.isfinite(float(lambda_graph))
        or float(lambda_graph) < 0.0
    ):
        return _abstained_network_result(
            graph,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={
                "message": "lambda_smoothness and lambda_graph must be finite and non-negative",
            },
        )

    # Pre-solve topology and physical evidence gating
    topo_gate = audit_graph_topology(graph, node_heads=node_heads)
    if not topo_gate.can_proceed:
        return _abstained_network_result(
            graph,
            grid,
            topo_gate.reason_codes,
            diagnostics={
                "topology_messages": topo_gate.messages,
                "topology_metrics": topo_gate.metrics,
            },
        )

    raw_nodes = list(graph.nodes())
    nodes = sorted(str(n) for n in raw_nodes)
    if len(set(nodes)) != len(raw_nodes):
        return _abstained_network_result(
            graph,
            grid,
            (CODE_NODE_ID_COLLISION,),
            diagnostics={
                "message": "Distinct graph node objects collapse to the same string identifier.",
                "node_ids": nodes,
            },
        )

    raw_node_by_id = {str(node_id): node_id for node_id in raw_nodes}
    n_nodes = len(nodes)
    k_bins = grid.n_bins
    if n_nodes == 0:
        return _abstained_network_result(
            graph,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={"message": "Cannot solve a network with zero nodes."},
        )

    missing_systems = [nid for nid in nodes if nid not in forward_systems]
    if missing_systems:
        return _abstained_network_result(
            graph,
            grid,
            (CODE_MISSING_FORWARD_SYSTEM,),
            diagnostics={
                "missing_nodes": missing_systems,
                "message": "Every graph node requires a validated multi-tracer forward system.",
            },
        )

    # Run the same physical-evidence gate used by the single-node solver for
    # every node before assembling a coupled objective.  A network estimate
    # must not silently drop a contradictory or unidentifiable node.
    node_gate_reports: Dict[str, Dict[str, Any]] = {}
    gate_reasons: list[str] = []
    for nid in nodes:
        system = forward_systems[nid]
        report = audit_node_physical_evidence(
            system,
            max_condition_number=max_condition_number,
        )
        node_gate_reports[nid] = {
            "status": report.status,
            "can_proceed": report.can_proceed,
            "reason_codes": report.reason_codes,
            "messages": report.messages,
            "metrics": report.metrics,
        }
        for reason in report.reason_codes:
            if reason not in gate_reasons:
                gate_reasons.append(reason)

        if system.grid.n_bins != k_bins or not np.allclose(system.grid.taus, grid.taus):
            reason = CODE_OPTIMIZATION_FAILURE
            if reason not in gate_reasons:
                gate_reasons.append(reason)
            node_gate_reports[nid]["reason_codes"] = tuple(
                list(report.reason_codes) + [reason]
            )
            node_gate_reports[nid]["messages"] = tuple(
                list(report.messages)
                + ["Forward-system age grid does not match the network age grid."]
            )

    if gate_reasons:
        return _abstained_network_result(
            graph,
            grid,
            tuple(gate_reasons),
            diagnostics={
                "node_gate_reports": node_gate_reports,
                "message": "At least one node failed the pre-solve physical evidence gate.",
            },
        )

    # Every graph edge needs a declared transport map.  Falling back to an
    # identity map would erase travel time and make the graph penalty appear
    # physically informed when it is not.
    missing_transport = [
        (str(u), str(v))
        for u, v in graph.edges()
        if (str(u), str(v)) not in transport_operators
    ]
    bad_transport_dimensions = [
        (str(u), str(v))
        for u, v in graph.edges()
        if (str(u), str(v)) in transport_operators
        and transport_operators[(str(u), str(v))].n_bins != k_bins
    ]
    if missing_transport:
        return _abstained_network_result(
            graph,
            grid,
            (CODE_MISSING_TRANSPORT_OPERATOR,),
            diagnostics={"missing_edges": missing_transport},
        )
    if bad_transport_dimensions:
        return _abstained_network_result(
            graph,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={
                "bad_transport_dimensions": bad_transport_dimensions,
                "message": "Transport operator dimensions do not match the network age grid.",
            },
        )

    d2 = build_mass_aware_curvature_matrix(grid)
    d2_quad = d2.T @ d2

    node_idx = {nid: idx for idx, nid in enumerate(nodes)}

    # Default recharge distribution: modern infiltration
    default_recharge = build_local_recharge_distribution(grid, mean_recharge_age_years=1.0)

    # Build node mixing parameters.  Roots have only local recharge by
    # definition; every non-root requires a declared mixing model unless the
    # caller explicitly opts into the historical compatibility prior.
    specs: Dict[str, NodeMixingSpecification] = {}
    missing_mixing: list[str] = []
    mixing_assumption = "explicit_node_mixing_specifications"
    for nid in nodes:
        if mixing_specs and nid in mixing_specs:
            specs[nid] = mixing_specs[nid]
        else:
            parents = [str(u) for u in graph.predecessors(raw_node_by_id[nid])]
            if parents:
                if not allow_default_mixing:
                    missing_mixing.append(nid)
                    continue
                eq_weight = 0.8 / len(parents)
                weights = {p: eq_weight for p in parents}
                specs[nid] = NodeMixingSpecification(
                    node_id=nid,
                    local_fraction=0.2,
                    upstream_weights=weights,
                    recharge_distribution=default_recharge,
                )
            else:
                specs[nid] = NodeMixingSpecification(
                    node_id=nid,
                    local_fraction=1.0,
                    upstream_weights={},
                    recharge_distribution=default_recharge,
                )

    if missing_mixing:
        return _abstained_network_result(
            graph,
            grid,
            (CODE_MISSING_MIXING_SPEC,),
            diagnostics={
                "missing_nodes": missing_mixing,
                "allow_default_mixing": bool(allow_default_mixing),
                "message": "Downstream nodes require explicit local and upstream mixing fractions.",
            },
        )

    if allow_default_mixing and any(
        str(u) != str(v) and str(v) in specs and str(u) in specs
        for u, v in graph.edges()
        if str(v) not in (mixing_specs or {})
    ):
        mixing_assumption = "compatibility_default_20_percent_local_80_percent_upstream"

    # Validate the dimensions and identifiers of explicit mixing declarations
    # before their composite distributions enter the objective.
    for nid, spec in specs.items():
        if str(spec.node_id) != nid:
            return _abstained_network_result(
                graph,
                grid,
                (CODE_OPTIMIZATION_FAILURE,),
                diagnostics={
                    "message": f"Mixing specification node_id {spec.node_id!r} does not match {nid!r}.",
                },
            )
        if spec.recharge_distribution.shape != (k_bins,):
            return _abstained_network_result(
                graph,
                grid,
                (CODE_OPTIMIZATION_FAILURE,),
                diagnostics={
                    "message": f"Recharge distribution for {nid!r} does not match the network age grid.",
                },
            )
        parents = {str(u) for u in graph.predecessors(raw_node_by_id[nid])}
        unknown_parents = sorted(set(spec.upstream_weights) - parents)
        if unknown_parents:
            return _abstained_network_result(
                graph,
                grid,
                (CODE_OPTIMIZATION_FAILURE,),
                diagnostics={
                    "node_id": nid,
                    "unknown_mixing_parents": unknown_parents,
                },
            )

    # Precompute per-node data Hessian and linear terms
    h_data_list: list[np.ndarray] = []
    f_data_list: list[np.ndarray] = []
    for nid in nodes:
        sys = forward_systems.get(nid)
        a_mat = sys.matrix
        c_obs = sys.observations
        inv_var = sys.weights / (sys.sigmas ** 2)
        h_data_list.append(a_mat.T @ (inv_var[:, None] * a_mat))
        f_data_list.append(-a_mat.T @ (inv_var * c_obs))

    lam_s = float(lambda_smoothness)
    lam_g = float(lambda_graph)

    # Packed variable x of shape (n_nodes * k_bins,)
    def unpack(x: np.ndarray) -> Dict[str, np.ndarray]:
        return {nid: x[i * k_bins : (i + 1) * k_bins] for i, nid in enumerate(nodes)}

    def objective(x: np.ndarray) -> float:
        total = 0.0
        g_dict = unpack(x)

        # 1. Data misfit and smoothness per node
        for i, nid in enumerate(nodes):
            gi = g_dict[nid]
            h_d = h_data_list[i]
            f_d = f_data_list[i]
            total += 0.5 * gi @ (h_d @ gi) + f_d @ gi
            total += lam_s * gi @ (d2_quad @ gi)

        # 2. Graph transport discrepancy
        if lam_g > 0.0:
            for nid in nodes:
                spec = specs[nid]
                if not spec.upstream_weights:
                    continue
                gi = g_dict[nid]
                # Expected arrival from upstream plus local recharge
                expected = spec.composite_distribution(g_dict, transport_operators)
                diff = gi - expected
                total += 0.5 * lam_g * float(diff @ diff)

        return float(total)

    def jacobian(x: np.ndarray) -> np.ndarray:
        grad = np.zeros_like(x)
        g_dict = unpack(x)

        # 1. Data and smoothness gradient
        for i, nid in enumerate(nodes):
            gi = g_dict[nid]
            h_d = h_data_list[i]
            f_d = f_data_list[i]
            g_slice = slice(i * k_bins, (i + 1) * k_bins)
            grad[g_slice] += h_d @ gi + f_d + 2.0 * lam_s * (d2_quad @ gi)

        # 2. Graph transport gradient
        if lam_g > 0.0:
            for nid in nodes:
                spec = specs[nid]
                if not spec.upstream_weights:
                    continue
                v_idx = node_idx[nid]
                gv = g_dict[nid]
                expected = spec.composite_distribution(g_dict, transport_operators)
                diff = gv - expected

                # Direct derivative w.r.t gv
                grad[v_idx * k_bins : (v_idx + 1) * k_bins] += lam_g * diff

                # Derivative w.r.t upstream parents: -lam_g * (pi_uv * T_uv)^T @ diff
                for parent, weight in spec.upstream_weights.items():
                    u_idx = node_idx[parent]
                    op = transport_operators.get((parent, nid))
                    if op is not None:
                        t_mat = op.matrix
                    else:
                        t_mat = np.eye(k_bins)
                    grad[u_idx * k_bins : (u_idx + 1) * k_bins] -= lam_g * weight * (t_mat.T @ diff)

        return grad

    # Initial guess: uniform mass for each node
    x0 = np.tile(np.ones(k_bins) / k_bins, n_nodes)
    bounds = [(0.0, 1.0) for _ in range(n_nodes * k_bins)]

    # Equality constraints: sum_k g_{ik} = 1.0 for each node i
    constraints = []
    for i in range(n_nodes):
        s_idx = i * k_bins
        e_idx = (i + 1) * k_bins
        constraints.append({"type": "eq", "fun": lambda x, s=s_idx, e=e_idx: np.sum(x[s:e]) - 1.0})

    opt = minimize(
        objective,
        x0,
        jac=jacobian,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 600, "ftol": 1e-9},
    )

    valid, optimizer_metrics = _optimizer_is_valid(
        opt,
        np.asarray(getattr(opt, "x", np.array([])), dtype=float),
        n_nodes,
        k_bins,
    )
    if not valid:
        return _abstained_network_result(
            graph,
            grid,
            (CODE_OPTIMIZATION_FAILURE,),
            diagnostics={
                "iterations": getattr(opt, "nit", None),
                "success": bool(getattr(opt, "success", False)),
                "message": str(getattr(opt, "message", "")),
                "node_gate_reports": node_gate_reports,
                **optimizer_metrics,
            },
        )

    g_final_dict = unpack(np.asarray(opt.x, dtype=float))
    node_results: Dict[str, SingleNodeTTDResult] = {}
    for nid in nodes:
        raw_g = np.maximum(0.0, g_final_dict[nid])
        g_sum = float(np.sum(raw_g))
        if not np.isfinite(g_sum) or g_sum <= 0.0:
            return _abstained_network_result(
                graph,
                grid,
                (CODE_OPTIMIZATION_FAILURE,),
                diagnostics={"message": f"Node {nid!r} has invalid optimized simplex mass."},
            )
        g_norm = raw_g / g_sum
        sys = forward_systems.get(nid)

        chi2 = sys.chi_squared(g_norm)
        chi2_per_observation = chi2 / max(1, sys.n_tracers)

        yf = grid.young_water_fraction(g_norm)
        mtt = grid.mean_transit_time(g_norm)
        ent = float(-np.sum(g_norm * np.log(np.maximum(g_norm, 1e-12))))

        node_results[nid] = SingleNodeTTDResult(
            node_id=nid,
            status="ESTIMATED",
            g=g_norm,
            young_water_fraction=yf,
            mean_transit_time_years=mtt,
            chi_squared=chi2,
            reduced_chi_squared=float("nan"),
            shannon_entropy=ent,
            chi_squared_per_observation=chi2_per_observation,
            effective_degrees_of_freedom=float("nan"),
            diagnostics={
                "regularization_semantics": "mass-aware curvature on density-equivalent vector",
                "identifiability": "regularized point estimate; full TTD is not identified by tracer count alone",
                **optimizer_metrics,
            },
        )

    # Edge transport discrepancy metrics
    edge_w1: Dict[Tuple[str, str], float] = {}
    edge_l2: Dict[Tuple[str, str], float] = {}
    total_sheaf_energy = 0.0

    for u, v in graph.edges():
        edge_key = (str(u), str(v))
        gu = node_results[str(u)].g
        gv = node_results[str(v)].g
        op = transport_operators.get(edge_key)
        if op is not None:
            gu_routed = op.transport(gu)
        else:
            gu_routed = gu

        l2_gap = float(np.linalg.norm(gv - gu_routed))
        w1_gap = wasserstein_1d(gv, gu_routed, grid)

        edge_l2[edge_key] = l2_gap
        edge_w1[edge_key] = w1_gap
        total_sheaf_energy += l2_gap ** 2

    elapsed = time.perf_counter() - t_start

    return NetworkTTDResult(
        status="ESTIMATED",
        node_results=node_results,
        edge_discrepancies_w1=edge_w1,
        edge_discrepancies_l2=edge_l2,
        network_sheaf_energy=total_sheaf_energy,
        abstention_reasons=(),
        diagnostics={
            "elapsed_seconds": elapsed,
            "iterations": opt.nit,
            "optimization_success": opt.success,
            "final_objective": float(opt.fun),
            **optimizer_metrics,
            "node_gate_reports": node_gate_reports,
            "mixing_assumption": mixing_assumption,
            "chi_squared_semantics": "chi_squared_per_observation is descriptive; reduced_chi_squared is NaN for non-parametric simplex fits",
        },
    )
