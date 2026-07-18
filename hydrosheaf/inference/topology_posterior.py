"""Bayesian topology posterior over candidate flow-network edges.

Defines the posterior:
    log p(G | data) = sum(edge inclusion priors)
                    - beta * sheaf_global_cost(G)
                    - edge_penalty * |E_G|
subject to constraint-violating edges being penalised.

Uses Metropolis-Hastings edge flips with bitmap-based state caching.
"""

import math
import random
from dataclasses import replace
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from ..config import Config
from ..graph.types import Edge
from ..log import get_logger

logger = get_logger("inference.topology_posterior")


def _get_edge_id(e: Any) -> str:
    if hasattr(e, "edge_id"):
        return str(e.edge_id)
    if isinstance(e, tuple):
        if len(e) == 3:
            return str(e[0])
        elif len(e) == 2:
            return f"{e[0]}->{e[1]}"
    return str(e)


def _get_edge_u_v(e: Any) -> Tuple[str, str]:
    if hasattr(e, "u") and hasattr(e, "v"):
        return str(e.u), str(e.v)
    if isinstance(e, tuple):
        if len(e) == 3:
            return str(e[1]), str(e[2])
        elif len(e) == 2:
            return str(e[0]), str(e[1])
    raise ValueError(f"Cannot resolve nodes from edge: {e}")


def _get_edge_attrs(e: Any) -> Dict[str, Any]:
    if hasattr(e, "attrs") and e.attrs is not None:
        return dict(e.attrs)
    return {}


def _get_config_float(config: Any, name: str, default: float) -> float:
    val = getattr(config, name, default)
    if val is None:
        return default
    if isinstance(val, dict):
        for key in ("value", "weight", "val", "default"):
            if key in val:
                try:
                    return float(val[key])
                except (ValueError, TypeError):
                    pass
        for v in val.values():
            try:
                return float(v)
            except (ValueError, TypeError):
                pass
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def _get_config_int(config: Any, name: str, default: int) -> int:
    return int(_get_config_float(config, name, float(default)))


def _edge_prior_probability(edge: Any) -> float:
    attrs = _get_edge_attrs(edge)
    value = attrs.get(
        "prior_edge_probability",
        attrs.get("edge_confidence", attrs.get("p_uv", 0.5)),
    )
    try:
        probability = float(value)
    except (TypeError, ValueError):
        probability = 0.5
    if not math.isfinite(probability):
        probability = 0.5
    return min(1.0 - 1e-12, max(1e-12, probability))


def _sample_map(samples: object) -> Dict[str, Mapping[str, object]]:
    if isinstance(samples, Mapping):
        return {str(k): v for k, v in samples.items()}
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, object]] = {}
        for row in samples:
            if not isinstance(row, Mapping):
                continue
            site_id = row.get("site_id") or row.get("node_id") or row.get("sample_id")
            if site_id is None:
                continue
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def _edges_to_bitset(edges: Sequence[Any], universe: Sequence[Any]) -> Tuple[int, int]:
    """Convert a list of selected edges to an integer bitset.

    Returns (bitset, n_edges).
    """
    edge_index = {_get_edge_id(e): i for i, e in enumerate(universe)}
    bitset = 0
    for e in edges:
        idx = edge_index.get(_get_edge_id(e))
        if idx is not None:
            bitset |= 1 << idx
    return bitset, len(universe)


def _bitset_to_edges(bitset: int, universe: Sequence[Any]) -> List[Any]:
    """Convert a bitset back to a list of edges."""
    return [universe[i] for i in range(len(universe)) if (bitset >> i) & 1]



def _compute_graph_cost(
    edges: Sequence[Any],
    universe: Sequence[Any],
    cost_fn: Callable[[Sequence[Any]], float],
    cache: Dict[str, float],
) -> float:
    """Compute or retrieve the graph cost from cache.

    Uses the universe for consistent bitset encoding so different edge subsets
    always get different cache keys.
    """
    bitset, _ = _edges_to_bitset(edges, universe)
    key = str(bitset)
    if key not in cache:
        cache[key] = cost_fn(list(edges))
    return cache[key]


def _log_posterior(
    edge_set: Sequence[Any],
    universe: Sequence[Any],
    cost_fn: Callable[[Sequence[Any]], float],
    cost_cache: Dict[str, float],
    beta: float,
    edge_penalty: float,
) -> float:
    """Compute unnormalized log-posterior for a given edge set.

    Uses the full Bernoulli log-prior over the universe:
        sum_{e in G} log(p_e) + sum_{e not in G} log(1 - p_e)
    """
    logp = 0.0

    # Full Bernoulli prior over all edges in the universe
    included_ids = {_get_edge_id(e) for e in edge_set}
    for e in universe:
        p = _edge_prior_probability(e)
        if _get_edge_id(e) in included_ids:
            logp += math.log(p)
        else:
            logp += math.log(1.0 - p)

    # Global cost term
    cost = _compute_graph_cost(edge_set, universe, cost_fn, cost_cache)
    logp -= beta * cost

    # Edge count penalty
    logp -= edge_penalty * len(edge_set)

    return logp


def _propose_flip(
    current: Sequence[Any],
    universe: Sequence[Any],
    rng: random.Random,
) -> Tuple[List[Any], str]:
    """Propose a symmetric single-edge toggle from the universe.

    Returns (proposed_edge_list, proposal_label).
    """
    current_set = {_get_edge_id(e) for e in current}
    universe_list = list(universe)
    if not universe_list:
        return list(current), "noop"

    edge = rng.choice(universe_list)
    edge_id = _get_edge_id(edge)
    if edge_id in current_set:
        proposed = [e for e in current if _get_edge_id(e) != edge_id]
        return proposed, f"remove_{edge_id}"

    proposed = list(current) + [edge]
    return proposed, f"add_{edge_id}"


def _propose_swap(
    current: Sequence[Any],
    universe: Sequence[Any],
    rng: random.Random,
) -> Tuple[List[Any], str]:
    """Symmetrically replace one included edge with one excluded edge."""

    current_ids = {_get_edge_id(edge) for edge in current}
    included = list(current)
    excluded = [
        edge for edge in universe if _get_edge_id(edge) not in current_ids
    ]
    if not included or not excluded:
        return _propose_flip(current, universe, rng)
    removed = rng.choice(included)
    added = rng.choice(excluded)
    proposed = [
        edge
        for edge in current
        if _get_edge_id(edge) != _get_edge_id(removed)
    ]
    proposed.append(added)
    return (
        proposed,
        f"swap_{_get_edge_id(removed)}_for_{_get_edge_id(added)}",
    )


def _propose_parent_swap(
    current: Sequence[Any],
    universe: Sequence[Any],
    rng: random.Random,
) -> Tuple[List[Any], str]:
    """Swap an incoming parent while preserving the child's in-degree.

    The proposal is symmetric because the eligible child and its included and
    excluded incoming-edge counts are unchanged by the replacement.
    """

    current_ids = {_get_edge_id(edge) for edge in current}
    included_by_v: Dict[str, List[Any]] = {}
    excluded_by_v: Dict[str, List[Any]] = {}
    for edge in universe:
        _, v = _get_edge_u_v(edge)
        target = (
            included_by_v
            if _get_edge_id(edge) in current_ids
            else excluded_by_v
        )
        target.setdefault(v, []).append(edge)
    eligible = sorted(set(included_by_v) & set(excluded_by_v))
    if not eligible:
        return _propose_swap(current, universe, rng)
    v = rng.choice(eligible)
    removed = rng.choice(included_by_v[v])
    added = rng.choice(excluded_by_v[v])
    proposed = [
        edge
        for edge in current
        if _get_edge_id(edge) != _get_edge_id(removed)
    ]
    proposed.append(added)
    return (
        proposed,
        f"parent_swap_{_get_edge_id(removed)}_for_{_get_edge_id(added)}",
    )


def run_topology_posterior(
    universe: Sequence[Any],
    cost_fn: Callable[[Sequence[Any]], float],
    config: Config,
    initial_edges: Optional[Sequence[Any]] = None,
    seed: int = 42,
) -> Dict[str, Any]:
    """Run Metropolis-Hastings to sample the topology posterior.

    Parameters
    ----------
    universe : sequence of Edge or tuple
        All candidate edges.
    cost_fn : callable
        Function mapping a list of edges to a scalar cost (lower is better).
    config : Config
        Hydrosheaf configuration with topology posterior settings.
    initial_edges : sequence of Edge, optional
        Starting point. Defaults to all edges.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with keys:
        edge_probabilities : Dict[str, float]
            Posterior inclusion probability per edge.
        edge_log_odds : Dict[str, float]
            Log-odds of inclusion per edge.
        map_edges : List[str]
            Edge IDs of the maximum a posteriori graph.
        entropy : float
            Legacy alias for the sum of marginal Bernoulli edge entropies.
        n_edges_mean : float
            Posterior mean number of edges.
        n_edges_ci95 : Tuple[float, float]
            95% credible interval for number of edges.
        acceptance_rate : float
            M-H acceptance rate.
    """
    n_samples = _get_config_int(config, "topology_posterior_samples", 2000)
    n_burnin = _get_config_int(config, "topology_posterior_burnin", 500)
    n_chains = _get_config_int(config, "topology_posterior_chains", 1)
    beta = _get_config_float(config, "topology_posterior_beta", 1.0)
    edge_penalty = _get_config_float(config, "topology_posterior_edge_penalty", 0.0)

    cost_cache: Dict[str, float] = {}
    universe_list = list(universe)
    edge_ids = [_get_edge_id(edge) for edge in universe_list]
    if not universe_list:
        return {
            "edge_probabilities": {},
            "edge_log_odds": {},
            "map_edges": [],
            "entropy": 0.0,
            "marginal_edge_entropy": 0.0,
            "entropy_definition": "sum_of_marginal_bernoulli_edge_entropies",
            "n_edges_mean": 0.0,
            "n_edges_ci95": (0.0, 0.0),
            "acceptance_rate": 0.0,
            "acceptance_rates_by_chain": [],
            "n_chains": n_chains,
            "n_edges_r_hat": None,
            "n_edges_ess": 0.0,
            "edge_r_hat": {},
            "edge_ess": {},
            "initial_edge_ids": [],
            "initial_edge_ids_by_chain": [],
        }

    aggregate_counts: Dict[str, int] = {edge_id: 0 for edge_id in edge_ids}
    chain_counts: List[Dict[str, int]] = []
    chain_edge_traces: List[Dict[str, List[int]]] = []
    chain_n_edges: List[List[int]] = []
    acceptance_rates: List[float] = []
    initial_edge_ids_by_chain: List[List[str]] = []
    best_logp = -math.inf
    best_edges: List[Any] = []
    total_samples = n_samples + n_burnin
    constrained = any(
        (
            _get_config_int(config, "topology_posterior_min_edges", 0) > 0,
            _get_config_int(config, "topology_posterior_max_out_degree", 0) > 0,
            bool(getattr(config, "topology_posterior_require_acyclic", False)),
            bool(
                getattr(
                    config, "topology_posterior_require_weak_connectivity", False
                )
            ),
            bool(
                getattr(
                    config,
                    "topology_posterior_require_root_reachability",
                    False,
                )
            ),
        )
    )
    base_initial = (
        list(initial_edges) if initial_edges is not None else list(universe_list)
    )

    for chain_index in range(n_chains):
        rng = random.Random(seed + 104729 * chain_index)
        if chain_index == 0 or constrained:
            # Empty/prior-random starts can be knowingly invalid under DAG,
            # connectivity, reachability, or degree constraints. Start every
            # constrained chain from the caller's feasible graph and let
            # independent burn-in trajectories disperse it.
            current = list(base_initial)
        elif chain_index % 2 == 1:
            current = []
        else:
            current = [
                edge
                for edge in universe_list
                if rng.random() < _edge_prior_probability(edge)
            ]
        initial_edge_ids_by_chain.append(
            [_get_edge_id(edge) for edge in current]
        )
        current_logp = _log_posterior(
            current, universe_list, cost_fn, cost_cache, beta, edge_penalty
        )
        if current_logp > best_logp:
            best_logp = current_logp
            best_edges = list(current)

        counts = {edge_id: 0 for edge_id in edge_ids}
        edge_trace = {edge_id: [] for edge_id in edge_ids}
        n_edges_trace: List[int] = []
        acceptance_count = 0
        for iteration in range(total_samples):
            if constrained and rng.random() < 0.85:
                proposed_edges, label = _propose_parent_swap(
                    current, universe_list, rng
                )
            else:
                proposed_edges, label = _propose_flip(
                    current, universe_list, rng
                )
            if label == "noop":
                continue
            proposed_logp = _log_posterior(
                proposed_edges,
                universe_list,
                cost_fn,
                cost_cache,
                beta,
                edge_penalty,
            )
            log_ratio = proposed_logp - current_logp
            if log_ratio >= 0.0 or rng.random() < math.exp(log_ratio):
                current = proposed_edges
                current_logp = proposed_logp
                if iteration >= n_burnin:
                    acceptance_count += 1

            if iteration >= n_burnin:
                included = {_get_edge_id(edge) for edge in current}
                for edge_id in edge_ids:
                    present = int(edge_id in included)
                    counts[edge_id] += present
                    edge_trace[edge_id].append(present)
                n_edges_trace.append(len(current))
                if current_logp > best_logp:
                    best_logp = current_logp
                    best_edges = list(current)

        for edge_id, count in counts.items():
            aggregate_counts[edge_id] += count
        chain_counts.append(counts)
        chain_edge_traces.append(edge_trace)
        chain_n_edges.append(n_edges_trace)
        acceptance_rates.append(acceptance_count / max(1, n_samples))

    n_post_samples = max(1, n_samples * n_chains)
    edge_probabilities = {
        edge_id: count / n_post_samples
        for edge_id, count in aggregate_counts.items()
    }
    edge_log_odds = {
        eid: math.log(max(1e-12, p / (1.0 - p + 1e-12))) for eid, p in edge_probabilities.items()
    }

    # Sum of marginal Bernoulli edge entropies. This is not the entropy of the
    # joint graph posterior because dependencies among edges are not represented.
    marginal_edge_entropy = 0.0
    for p in edge_probabilities.values():
        if 0.0 < p < 1.0:
            marginal_edge_entropy -= p * math.log(p) + (1.0 - p) * math.log(1.0 - p)

    flattened_n_edges = [
        value for trace in chain_n_edges for value in trace
    ]
    n_edges_mean = (
        float(np.mean(flattened_n_edges)) if flattened_n_edges else 0.0
    )
    if flattened_n_edges:
        sorted_trace = sorted(flattened_n_edges)
        lo_idx = int(0.025 * len(sorted_trace))
        hi_idx = int(0.975 * len(sorted_trace))
        n_edges_ci95 = (float(sorted_trace[max(0, lo_idx)]), float(sorted_trace[min(len(sorted_trace) - 1, hi_idx)]))
    else:
        n_edges_ci95 = (0.0, 0.0)

    acceptance_rate = (
        float(np.mean(acceptance_rates)) if acceptance_rates else 0.0
    )
    n_edges_r_hat: Optional[float] = None
    n_edges_ess = float(len(flattened_n_edges))
    edge_r_hat: Dict[str, Optional[float]] = {}
    edge_ess: Dict[str, float] = {}
    if n_chains > 1 and n_samples > 1:
        from ..uncertainty.bayesian import compute_ess, compute_r_hat

        n_edge_array = np.asarray(chain_n_edges, dtype=float)
        n_edges_r_hat = compute_r_hat(n_edge_array)
        n_edges_ess = sum(compute_ess(chain) for chain in n_edge_array)
        for edge_id in edge_ids:
            traces = np.asarray(
                [chain[edge_id] for chain in chain_edge_traces],
                dtype=float,
            )
            chain_means = np.mean(traces, axis=1)
            if np.all(np.var(traces, axis=1) == 0.0) and not np.allclose(
                chain_means, chain_means[0]
            ):
                edge_r_hat[edge_id] = None
            else:
                edge_r_hat[edge_id] = compute_r_hat(traces)
            edge_ess[edge_id] = sum(compute_ess(chain) for chain in traces)

    logger.info(
        "Topology posterior: %d edges, %d chains, mean edges=%.1f, "
        "marginal entropy=%.2f, accept=%.3f",
        len(universe_list),
        n_chains,
        n_edges_mean,
        marginal_edge_entropy,
        acceptance_rate,
    )

    return {
        "edge_probabilities": edge_probabilities,
        "edge_log_odds": edge_log_odds,
        "map_edges": [_get_edge_id(e) for e in best_edges],
        "entropy": marginal_edge_entropy,
        "marginal_edge_entropy": marginal_edge_entropy,
        "entropy_definition": "sum_of_marginal_bernoulli_edge_entropies",
        "n_edges_mean": n_edges_mean,
        "n_edges_ci95": n_edges_ci95,
        "acceptance_rate": acceptance_rate,
        "acceptance_rates_by_chain": acceptance_rates,
        "n_chains": n_chains,
        "n_edges_r_hat": n_edges_r_hat,
        "n_edges_ess": float(n_edges_ess),
        "edge_r_hat": edge_r_hat,
        "edge_ess": edge_ess,
        "initial_edge_ids": initial_edge_ids_by_chain[0],
        "initial_edge_ids_by_chain": initial_edge_ids_by_chain,
    }


def make_topology_cost_fn(
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    node_info: Optional[Mapping[str, object]] = None,
    stats: Optional[object] = None,
    reference_distance_km: Optional[float] = None,
    gradient_map: Optional[Mapping[str, tuple[float, float, float]]] = None,
    local_residuals: Optional[Dict[str, float]] = None,
) -> Callable[[Sequence[Any]], float]:
    """Build a cost function that evaluates topology plausibility for a graph.

    The graph cost combines:
      1. Edge-local sheaf scores when chemistry/isotope data exist.
      2. Global sheaf residual energy when chemistry vectors exist.
      3. Projected-gradient alignment cost when *gradient_map* provided.
      4. Local head-plane residuals when *local_residuals* provided.
    """
    from ..data.schema import vector_from_sample
    from ..sheaf.hydraulic_hodge import head_hodge_graph_cost
    from ..sheaf.isotope_metrics import compute_isotope_stats
    from ..sheaf.topology_refine import _build_node_info, _score_candidates
    from ..sheaf.directed_section import build_edge_maps, solve_directed_section

    if stats is None:
        stats = compute_isotope_stats(sample_map.values(), config)
    if node_info is None:
        node_info = _build_node_info(sample_map, stats, config)

    node_vectors = {}
    for node_id, sample in sample_map.items():
        values, _ = vector_from_sample(
            sample,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        if values is not None:
            node_vectors[node_id] = values

    global_weight = _get_config_float(config, "sheaf_weight_global", 1.0)
    hydraulic_weight = _get_config_float(config, "hydraulic_hodge_weight", 1.0)
    invalid_cost = _get_config_float(
        config, "topology_posterior_invalid_cost", 1.0e6
    )
    likelihood_config = replace(config, sheaf_weight_head_prior=0.0)

    def constraint_cost(edges: Sequence[Any]) -> float:
        min_edges = _get_config_int(config, "topology_posterior_min_edges", 0)
        if len(edges) < min_edges:
            return invalid_cost

        require_acyclic = bool(
            getattr(config, "topology_posterior_require_acyclic", False)
        )
        require_connected = bool(
            getattr(config, "topology_posterior_require_weak_connectivity", False)
        )
        require_reachability = bool(
            getattr(config, "topology_posterior_require_root_reachability", False)
        )
        max_out_degree = _get_config_int(
            config, "topology_posterior_max_out_degree", 0
        )
        if not (
            require_acyclic
            or require_connected
            or require_reachability
            or max_out_degree > 0
        ):
            return 0.0

        import networkx as nx

        graph = nx.DiGraph()
        graph.add_nodes_from(str(node_id) for node_id in sample_map)
        graph.add_edges_from(_get_edge_u_v(edge) for edge in edges)
        if require_acyclic and not nx.is_directed_acyclic_graph(graph):
            return invalid_cost
        if require_connected and (
            graph.number_of_nodes() > 1
            and (
                graph.number_of_edges() == 0
                or not nx.is_weakly_connected(graph)
            )
        ):
            return invalid_cost
        if max_out_degree > 0 and any(
            degree > max_out_degree for _, degree in graph.out_degree()
        ):
            return invalid_cost
        if require_reachability:
            roots = [
                str(node)
                for node in getattr(config, "topology_posterior_root_nodes", [])
                if str(node) in graph
            ]
            if not roots:
                return invalid_cost
            reachable = set(roots)
            for root in roots:
                reachable.update(nx.descendants(graph, root))
            if reachable != set(graph.nodes):
                return invalid_cost
        return 0.0

    def cost_fn(edges: Sequence[Any]) -> float:
        invalid = constraint_cost(edges)
        if invalid > 0.0:
            return invalid
        if not edges:
            return 0.0

        # Score likelihood terms without reusing the hydraulic Bernoulli prior.
        scores = _score_candidates(
            edges, node_info, sample_map, stats, likelihood_config
        )
        local_cost = sum(s.local_score for s in scores.values())

        if node_vectors:
            try:
                maps = build_edge_maps(
                    edges, node_vectors, config,
                    prior_weight=0.0,
                )
                node_estimates = solve_directed_section(
                    list(node_vectors.keys()),
                    list(maps.values()),
                    node_vectors,
                    obs_weight=1.0,
                    diag_eps=1e-6,
                )
                from ..sheaf.directed_section import compute_edge_section_residuals

                residuals = compute_edge_section_residuals(
                    maps, node_estimates, config.weights
                )
                global_cost = sum(residuals.values())
            except Exception:
                logger.warning(
                    "Sheaf section solve failed; assigning high cost penalty.",
                    exc_info=True,
                )
                global_cost = 1e6  # large penalty, not zero (avoids biasing MCMC)
        else:
            global_cost = 0.0

        hydraulic_cost = 0.0
        if getattr(config, "hydraulic_hodge_enabled", False):
            edge_list = [edge for edge in edges if isinstance(edge, Edge)]
            if edge_list and len(edge_list) < len(edges):
                logger.warning(
                    "hydraulic_hodge_enabled=True but %d/%d edges are not Edge objects; "
                    "hydraulic cost computed on reduced set.",
                    len(edges) - len(edge_list),
                    len(edges),
                )
            if edge_list:
                hydraulic_cost = head_hodge_graph_cost(
                    edge_list,
                    sample_map,
                    config,
                    reference_distance_km=reference_distance_km,
                    gradient_map=gradient_map,
                    local_residuals=local_residuals,
                )

        return local_cost + global_weight * global_cost + hydraulic_weight * hydraulic_cost

    return cost_fn


def select_posterior_edges(
    candidate_edges: Sequence[Any],
    posterior_result: Dict[str, Any],
    max_neighbors: Optional[int] = None,
    probability_threshold: Optional[float] = None,
) -> List[Any]:
    probs = posterior_result.get("edge_probabilities", {})
    map_edge_ids = set(posterior_result.get("map_edges", []))

    if max_neighbors is not None and max_neighbors > 0 and isinstance(probs, Mapping):
        threshold = 0.0 if probability_threshold is None else float(probability_threshold)
        grouped: Dict[str, List[Any]] = {}
        for edge in candidate_edges:
            prob = float(probs.get(_get_edge_id(edge), 0.0) or 0.0)
            if prob < threshold:
                continue
            u, _ = _get_edge_u_v(edge)
            grouped.setdefault(u, []).append(edge)

        trimmed: List[Any] = []
        for edges in grouped.values():
            ranked = sorted(
                edges,
                key=lambda edge: float(probs.get(_get_edge_id(edge), 0.0) or 0.0),
                reverse=True,
            )
            trimmed.extend(ranked[:max_neighbors])
        if trimmed:
            return trimmed

    selected = [edge for edge in candidate_edges if _get_edge_id(edge) in map_edge_ids]
    if selected:
        return selected

    if isinstance(probs, Mapping):
        threshold = 0.5 if probability_threshold is None else float(probability_threshold)
        return [
            edge
            for edge in candidate_edges
            if float(probs.get(_get_edge_id(edge), 0.0) or 0.0) >= threshold
        ]

    return []


def infer_topology_map_edges(
    samples: object,
    candidate_edges: Sequence[Edge],
    config: Config,
    initial_edges: Optional[Sequence[Edge]] = None,
    max_neighbors: Optional[int] = None,
    probability_threshold: Optional[float] = None,
    seed: int = 42,
    grid_geometry: Optional[object] = None,
    head_map: Optional[Dict[int, float]] = None,
    node_df: Optional[object] = None,
) -> Tuple[List[Edge], Dict[str, Any]]:
    """Run the full topology posterior pipeline.

    Parameters
    ----------
    samples, candidate_edges, config, initial_edges, max_neighbors,
    probability_threshold, seed:
        As in :func:`run_topology_posterior`.
    grid_geometry:
        Optional :class:`~hydrosheaf.physics.modflow_head.GridGeometry` for
        continuous head-gradient computation.
    head_map:
        Optional ``{cell_index: head_value}`` from :func:`~hydrosheaf.physics.modflow_head.parse_fhd`.
    node_df:
        Optional DataFrame with ``node_id`` column (``"cell_{int}"`` convention)
        for joining gradients to nodes.
    """
    sample_map = _sample_map(samples)

    # --- Pre-compute diagnostics used by the posterior ---
    gradient_map: Optional[Dict[str, tuple[float, float, float]]] = None
    local_residuals: Optional[Dict[str, float]] = None

    # Local head-plane residuals can run standalone — no gradient prerequisites.
    if getattr(config, "head_plane_residual_enabled", False):
        from ..sheaf.hydraulic_hodge import (
            compute_local_head_plane_residuals,
            extract_node_heads,
        )

        head_dict = extract_node_heads(sample_map, config)
        n_neighbors = int(
            getattr(config, "head_plane_residual_neighbors", 8) or 8
        )
        local_residuals = compute_local_head_plane_residuals(
            head_dict,
            sample_map,
            n_neighbors=n_neighbors,
        )

    if (
        getattr(config, "projected_gradient_enabled", False)
        and grid_geometry is not None
        and head_map is not None
        and node_df is not None
    ):
        from ..physics.modflow_head import compute_head_gradient, map_gradient_to_nodes

        sigma = _get_config_float(config, "projected_gradient_smoothing_sigma", 1.0)
        cell_gradient = compute_head_gradient(head_map, grid_geometry, sigma=sigma)
        gradient_map = map_gradient_to_nodes(cell_gradient, node_df)

    # --- Apply flow-direction priors ---
    if getattr(config, "projected_gradient_enabled", False) and gradient_map is not None:
        from ..graph.flow_direction import apply_flow_direction_priors

        apply_flow_direction_priors(candidate_edges, sample_map, config, gradient_map)
    elif getattr(config, "steepest_descent_enabled", False):
        from ..graph.flow_direction import apply_steepest_descent_priors

        apply_steepest_descent_priors(candidate_edges, sample_map, config)

    # --- Reference distance for backward compat ---
    reference_distance_km = None
    if getattr(config, "hydraulic_hodge_enabled", False):
        from ..sheaf.hydraulic_hodge import infer_reference_distance_km

        reference_distance_km = infer_reference_distance_km(
            candidate_edges, sample_map, config,
        )

    # --- Build cost function with all diagnostics ---
    cost_fn = make_topology_cost_fn(
        sample_map=sample_map,
        config=config,
        reference_distance_km=reference_distance_km,
        gradient_map=gradient_map,
        local_residuals=local_residuals,
    )
    posterior_result = run_topology_posterior(
        universe=candidate_edges,
        cost_fn=cost_fn,
        config=config,
        initial_edges=initial_edges,
        seed=seed,
    )
    selected = select_posterior_edges(
        candidate_edges=candidate_edges,
        posterior_result=posterior_result,
        max_neighbors=max_neighbors,
        probability_threshold=probability_threshold,
    )
    attach_posterior_attrs(
        selected_edges=selected,
        candidate_edges=candidate_edges,
        posterior_result=posterior_result,
        mode="diagnostic",
    )
    return [edge for edge in selected if isinstance(edge, Edge)], posterior_result


def attach_posterior_attrs(
    selected_edges: List[Any],
    candidate_edges: Sequence[Any],
    posterior_result: Dict[str, Any],
    mode: str = "diagnostic",
) -> None:
    """Attach posterior attributes to selected edges.

    Parameters
    ----------
    selected_edges : list of Edge
        Currently selected edges (may be modified in 'select' mode).
    candidate_edges : sequence of Edge
        All candidate edges in the universe.
    posterior_result : dict
        Output from run_topology_posterior.
    mode : str
        "diagnostic": attach attrs, keep current selection.
        "select": use MAP edges (not yet implemented as replacement).
    """
    probs = posterior_result.get("edge_probabilities", {})
    log_odds = posterior_result.get("edge_log_odds", {})
    map_edges = set(posterior_result.get("map_edges", []))
    entropy = posterior_result.get("entropy", 0.0)
    n_mean = posterior_result.get("n_edges_mean", 0.0)
    ci = posterior_result.get("n_edges_ci95", (0.0, 0.0))
    accept = posterior_result.get("acceptance_rate", 0.0)

    for edge in candidate_edges:
        if isinstance(edge, tuple) or not hasattr(edge, "attrs"):
            continue
        attrs = dict(edge.attrs or {})
        eid = _get_edge_id(edge)
        attrs["posterior_edge_probability"] = probs.get(eid)
        attrs["posterior_edge_log_odds"] = log_odds.get(eid)
        attrs["posterior_map_selected"] = eid in map_edges
        attrs["posterior_topology_entropy"] = entropy
        attrs["posterior_marginal_edge_entropy"] = posterior_result.get(
            "marginal_edge_entropy"
        )
        attrs["posterior_entropy_definition"] = posterior_result.get(
            "entropy_definition"
        )
        attrs["posterior_n_edges_mean"] = n_mean
        attrs["posterior_n_edges_ci95"] = f"{ci[0]:.2f}-{ci[1]:.2f}"
        attrs["posterior_acceptance_rate"] = accept
        edge.attrs = attrs
