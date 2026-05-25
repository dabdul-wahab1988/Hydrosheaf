"""Sheaf cohomology diagnostics for global flow-network consistency.

Computes obstruction energies from the affine sheaf section problem:
    sqrt(w_e) * (alpha_e * x_u + offset_e - x_v) ~ 0

The quantity of interest is not just H1 (dimension of obstruction space) but the
*affine obstruction energy*: how much the fitted transport/reaction relations
fail to close globally.  A non-zero obstruction energy flags a physical
inconsistency (e.g. a cycle where chemistry cannot be simultaneously satisfied).
"""

from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
from scipy.optimize import lsq_linear

from ..log import get_logger
from .directed_section import DirectedEdgeMap

logger = get_logger("sheaf.cohomology")


def _build_graph(edge_maps: Iterable[DirectedEdgeMap]) -> nx.DiGraph:
    """Build a NetworkX DiGraph from directed edge maps."""
    g = nx.DiGraph()
    for em in edge_maps:
        g.add_edge(em.edge.u, em.edge.v, edge_map=em)
    return g


def _node_index_map(edge_maps: Iterable[DirectedEdgeMap]) -> Tuple[Dict[str, int], List[str]]:
    """Assign dense indices to nodes."""
    nodes: Set[str] = set()
    for em in edge_maps:
        nodes.add(em.edge.u)
        nodes.add(em.edge.v)
    node_list = sorted(nodes)
    idx = {n: i for i, n in enumerate(node_list)}
    return idx, node_list


def build_coboundary_matrix(
    edge_maps: Iterable[DirectedEdgeMap],
    dim: int,
) -> csr_matrix:
    """Build the sparse coboundary matrix D for the sheaf section problem.

    Each edge e contributes one row per species dimension:
        D_e_row = [..., sqrt(w_e) * alpha_e, ..., -sqrt(w_e), ...]
    where the nonzeros are at columns u and v.

    Parameters
    ----------
    edge_maps : iterable of DirectedEdgeMap
        The directed edge maps from the section solver.
    dim : int
        Number of chemical species dimensions.

    Returns
    -------
    D : csr_matrix, shape (n_edges * dim, n_nodes * dim)
        Sparse coboundary matrix.
    """
    idx, node_list = _node_index_map(edge_maps)
    n_nodes = len(node_list)
    n_edges = len(list(edge_maps))
    n_rows = n_edges * dim
    n_cols = n_nodes * dim

    D = lil_matrix((n_rows, n_cols), dtype=float)

    for row_base, em in enumerate(edge_maps):
        u_i = idx.get(em.edge.u)
        v_i = idx.get(em.edge.v)
        if u_i is None or v_i is None:
            continue

        w_sqrt = float(em.weight) ** 0.5
        alpha = float(em.alpha)

        for d in range(dim):
            row = row_base * dim + d
            col_u = u_i * dim + d
            col_v = v_i * dim + d
            D[row, col_u] = w_sqrt * alpha
            D[row, col_v] = -w_sqrt

    return D.tocsr()


def build_rhs_vector(
    edge_maps: Sequence[DirectedEdgeMap],
    dim: int,
) -> np.ndarray:
    """Build the right-hand side vector b for the affine coboundary equation.

    b_e = -sqrt(w_e) * offset_e[d]
    """
    n_edges = len(edge_maps)
    n_rows = n_edges * dim
    b = np.zeros(n_rows)

    for row_base, em in enumerate(edge_maps):
        w_sqrt = float(em.weight) ** 0.5
        for d in range(min(dim, len(em.offset))):
            row = row_base * dim + d
            b[row] = -w_sqrt * float(em.offset[d])

    return b


def compute_cohomology(
    edge_maps: Sequence[DirectedEdgeMap],
    dim: Optional[int] = None,
) -> Dict[str, object]:
    """Compute sheaf cohomology diagnostics for a set of directed edge maps.

    Parameters
    ----------
    edge_maps : sequence of DirectedEdgeMap
        Edge maps from the section solver.
    dim : int, optional
        Chemical dimension. Auto-detected from the first edge map if not given.

    Returns
    -------
    dict with keys:
        h0_dim : int
            Dimension of global section space (n_vertex_dofs - rank(D)).
        h1_dim : int
            Dimension of first cohomology (n_edge_dofs - rank(D)).
        obstruction_energy : float
            Minimum squared residual of the affine system: min ||D x - b||^2.
        rank_D : int
            Numerical rank of the coboundary matrix.
    """
    edge_maps_list = list(edge_maps)
    if not edge_maps_list:
        return {
            "h0_dim": 0,
            "h1_dim": 0,
            "obstruction_energy": 0.0,
            "rank_D": 0,
        }

    if dim is None:
        dim = len(edge_maps_list[0].offset)

    idx, node_list = _node_index_map(edge_maps_list)
    n_nodes = len(node_list)

    D = build_coboundary_matrix(edge_maps_list, dim)
    b = build_rhs_vector(edge_maps_list, dim)

    # Numerical rank.  Hydrosheaf matrices are small (n_nodes * dim is
    # typically < 2000), so a dense SVD via np.linalg.matrix_rank is safe.
    if D.shape[0] > 0 and D.shape[1] > 0:
        try:
            rank_D = int(np.linalg.matrix_rank(D.toarray()))
        except Exception:
            rank_D = min(D.shape[0], D.shape[1])
    else:
        rank_D = 0

    n_vertex_dofs = n_nodes * dim
    n_edge_dofs = D.shape[0]
    h0_dim = max(0, n_vertex_dofs - rank_D)
    h1_dim = max(0, n_edge_dofs - rank_D)

    # Affine obstruction energy
    if D.shape[0] > 0 and D.shape[1] > 0:
        try:
            res = lsq_linear(D, b, bounds=(-np.inf, np.inf), method="trf", tol=1e-8)
            obstruction_energy = float(np.sum(res.fun ** 2))
        except Exception:
            # Fallback: compute min-norm residual from normal equations
            try:
                x_sol, _, _, _ = np.linalg.lstsq(D.toarray(), b, rcond=None)
                obstruction_energy = float(np.sum((D @ x_sol - b) ** 2))
            except Exception:
                obstruction_energy = float(np.sum(b ** 2))
    else:
        obstruction_energy = 0.0

    return {
        "h0_dim": h0_dim,
        "h1_dim": h1_dim,
        "obstruction_energy": obstruction_energy,
        "rank_D": rank_D,
    }


def compute_edge_leverage(
    edge_maps: Sequence[DirectedEdgeMap],
    dim: Optional[int] = None,
) -> Dict[str, float]:
    """Compute per-edge leverage scores.

    Edge leverage = obstruction_energy(G) - obstruction_energy(G \\ {e}).
    A large positive value means removing this edge significantly reduces
    global obstruction, i.e. this edge is a primary source of inconsistency.
    """
    edge_maps_list = list(edge_maps)
    if dim is None and edge_maps_list:
        dim = len(edge_maps_list[0].offset)

    if dim is None:
        return {}

    base_energy = float(compute_cohomology(edge_maps_list, dim)["obstruction_energy"])
    leverage: Dict[str, float] = {}

    for i, em in enumerate(edge_maps_list):
        reduced = edge_maps_list[:i] + edge_maps_list[i + 1 :]
        reduced_energy = float(compute_cohomology(reduced, dim)["obstruction_energy"])
        leverage[em.edge.edge_id] = base_energy - reduced_energy

    return leverage


def compute_cycle_obstruction(
    edge_maps: Sequence[DirectedEdgeMap],
    dim: Optional[int] = None,
) -> Dict[str, object]:
    """Compute obstruction energy per graph cycle.

    Uses networkx.cycle_basis on the undirected skeleton. For each cycle,
    extract the subgraph of edge maps and compute its obstruction energy.

    Returns
    -------
    dict with keys:
        cycle_obstruction_max : float
            Maximum per-cycle obstruction energy.
        cycle_count : int
            Number of cycles found.
        cycle_energies : List[Tuple[List[str], float]]
            (cycle_node_ids, obstruction_energy) per cycle.
    """
    edge_maps_list = list(edge_maps)
    if dim is None and edge_maps_list:
        dim = len(edge_maps_list[0].offset)

    if dim is None or not edge_maps_list:
        return {"cycle_obstruction_max": 0.0, "cycle_count": 0, "cycle_energies": []}

    # Build undirected graph for cycle detection
    ug = nx.Graph()
    edge_by_nodes: Dict[Tuple[str, str], DirectedEdgeMap] = {}
    for em in edge_maps_list:
        ug.add_edge(em.edge.u, em.edge.v)
        edge_by_nodes[(em.edge.u, em.edge.v)] = em
        edge_by_nodes[(em.edge.v, em.edge.u)] = em

    cycle_energies: List[Tuple[List[str], float]] = []

    try:
        cycles = nx.cycle_basis(ug)
    except Exception:
        cycles = []

    for cycle_nodes in cycles:
        cycle_ems: List[DirectedEdgeMap] = []
        for i in range(len(cycle_nodes)):
            u, v = cycle_nodes[i], cycle_nodes[(i + 1) % len(cycle_nodes)]
            em = edge_by_nodes.get((u, v))
            if em is not None:
                cycle_ems.append(em)

        if cycle_ems:
            coh = compute_cohomology(cycle_ems, dim)
            cycle_energies.append((cycle_nodes, float(coh["obstruction_energy"])))

    max_obstruction = max((e for _, e in cycle_energies), default=0.0)

    return {
        "cycle_obstruction_max": max_obstruction,
        "cycle_count": len(cycles),
        "cycle_energies": cycle_energies,
    }


def attach_cohomology_attrs(
    selected_edges: List[object],
    edge_maps: Mapping[str, DirectedEdgeMap],
    dim: Optional[int] = None,
    compute_leverage: bool = True,
) -> None:
    """Compute cohomology diagnostics and attach per-edge attributes to selected edges.

    Parameters
    ----------
    selected_edges : list of Edge
        Edges that were selected by the sheaf refinement.
    edge_maps : dict mapping edge_id -> DirectedEdgeMap
        All edge maps used in the section solve.
    dim : int, optional
        Chemical dimension.
    compute_leverage : bool
        Whether to compute per-edge leverage scores (O(n_edges^2) SVDs).
    """
    if not selected_edges:
        return

    selected_maps = [edge_maps[e.edge_id] for e in selected_edges if e.edge_id in edge_maps]
    if not selected_maps:
        return

    coh = compute_cohomology(selected_maps, dim)
    cycle_info = compute_cycle_obstruction(selected_maps, dim)

    leverage: Dict[str, float] = {}
    if compute_leverage and len(selected_maps) <= 200:
        leverage = compute_edge_leverage(selected_maps, dim)
    elif compute_leverage:
        logger.info("Skipping per-edge leverage (too many edges: %d).", len(selected_maps))

    for edge in selected_edges:
        attrs = dict(edge.attrs or {})
        attrs["sheaf_h0_dim"] = coh["h0_dim"]
        attrs["sheaf_h1_dim"] = coh["h1_dim"]
        attrs["sheaf_obstruction_energy"] = coh["obstruction_energy"]
        attrs["sheaf_obstruction_leverage"] = leverage.get(edge.edge_id)
        attrs["sheaf_cycle_obstruction_max"] = cycle_info["cycle_obstruction_max"]
        attrs["sheaf_cycle_count"] = cycle_info["cycle_count"]
        edge.attrs = attrs
