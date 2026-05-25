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
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from ..config import Config
from ..graph.types import Edge
from ..log import get_logger

logger = get_logger("inference.topology_posterior")


def _edges_to_bitset(edges: Sequence[Edge], universe: Sequence[Edge]) -> Tuple[int, int]:
    """Convert a list of selected edges to an integer bitset.

    Returns (bitset, n_edges).
    """
    edge_index = {e.edge_id: i for i, e in enumerate(universe)}
    bitset = 0
    for e in edges:
        idx = edge_index.get(e.edge_id)
        if idx is not None:
            bitset |= 1 << idx
    return bitset, len(universe)


def _bitset_to_edges(bitset: int, universe: Sequence[Edge]) -> List[Edge]:
    """Convert a bitset back to a list of edges."""
    return [universe[i] for i in range(len(universe)) if (bitset >> i) & 1]



def _compute_graph_cost(
    edges: Sequence[Edge],
    universe: Sequence[Edge],
    cost_fn: Callable[[Sequence[Edge]], float],
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


def _edge_logit_prior(edge: Edge) -> float:
    """Logit-prior for including an edge: log(p/(1-p)).

    Positive values encourage inclusion; negative discourage it.
    Mathematically equivalent to the full Bernoulli log-prior up to
    a state-independent constant that cancels in MH ratios.
    """
    attrs = edge.attrs or {}
    prior = attrs.get("edge_confidence", attrs.get("p_uv", 0.5))
    p = float(prior) if prior is not None else 0.5
    p = min(1.0 - 1e-12, max(1e-12, p))
    if p >= 1.0 - 1e-12:
        return 20.0  # log(1e12) ~ 27.6, clip for numerical safety
    if p <= 1e-12:
        return -20.0
    return math.log(p / (1.0 - p))


def _log_posterior(
    edge_set: Sequence[Edge],
    universe: Sequence[Edge],
    cost_fn: Callable[[Sequence[Edge]], float],
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
    included_ids = {e.edge_id for e in edge_set}
    for e in universe:
        if e.edge_id in included_ids:
            attrs = e.attrs or {}
            p = float(attrs.get("edge_confidence", attrs.get("p_uv", 0.5)) or 0.5)
            p = min(1.0 - 1e-12, max(1e-12, p))
            logp += math.log(p)
        else:
            attrs = e.attrs or {}
            p = float(attrs.get("edge_confidence", attrs.get("p_uv", 0.5)) or 0.5)
            p = min(1.0 - 1e-12, max(1e-12, p))
            logp += math.log(1.0 - p)

    # Global cost term
    cost = _compute_graph_cost(edge_set, universe, cost_fn, cost_cache)
    logp -= beta * cost

    # Edge count penalty
    logp -= edge_penalty * len(edge_set)

    return logp


def _propose_flip(
    current: Sequence[Edge],
    universe: Sequence[Edge],
    rng: random.Random,
) -> Tuple[List[Edge], str]:
    """Propose adding or removing a random edge from universe.

    Returns (proposed_edge_list, proposal_label).
    """
    current_set = {e.edge_id for e in current}
    universe_list = list(universe)

    # Randomly choose add or remove
    can_remove = len(current_set) > 0
    can_add = len(current_set) < len(universe_list)

    if can_remove and (not can_add or rng.random() < 0.5):
        # Remove a random edge from current
        removal = rng.choice(list(current))
        proposed = [e for e in current if e.edge_id != removal.edge_id]
        return proposed, f"remove_{removal.edge_id}"
    elif can_add:
        # Add a random edge not in current
        absent = [e for e in universe_list if e.edge_id not in current_set]
        addition = rng.choice(absent)
        proposed = list(current) + [addition]
        return proposed, f"add_{addition.edge_id}"
    else:
        return list(current), "noop"


def run_topology_posterior(
    universe: Sequence[Edge],
    cost_fn: Callable[[Sequence[Edge]], float],
    config: Config,
    initial_edges: Optional[Sequence[Edge]] = None,
    seed: int = 42,
) -> Dict[str, Any]:
    """Run Metropolis-Hastings to sample the topology posterior.

    Parameters
    ----------
    universe : sequence of Edge
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
            Posterior entropy over graphs.
        n_edges_mean : float
            Posterior mean number of edges.
        n_edges_ci95 : Tuple[float, float]
            95% credible interval for number of edges.
        acceptance_rate : float
            M-H acceptance rate.
    """
    n_samples = int(getattr(config, "topology_posterior_samples", 2000))
    n_burnin = int(getattr(config, "topology_posterior_burnin", 500))
    beta = float(getattr(config, "topology_posterior_beta", 1.0))
    edge_penalty = float(getattr(config, "topology_posterior_edge_penalty", 0.0))

    rng = random.Random(seed)
    cost_cache: Dict[str, float] = {}

    universe_list = list(universe)
    if not universe_list:
        return {
            "edge_probabilities": {},
            "edge_log_odds": {},
            "map_edges": [],
            "entropy": 0.0,
            "n_edges_mean": 0.0,
            "n_edges_ci95": (0.0, 0.0),
            "acceptance_rate": 0.0,
        }

    # Initialize
    current = list(initial_edges) if initial_edges else list(universe_list)
    current_logp = _log_posterior(current, universe_list, cost_fn, cost_cache, beta, edge_penalty)
    best_logp = current_logp
    best_edges = list(current)

    # Track per-edge inclusion counts
    edge_counts: Dict[str, int] = {e.edge_id: 0 for e in universe_list}
    n_edges_trace: List[int] = []
    acceptance_count = 0

    total_samples = n_samples + n_burnin
    for iteration in range(total_samples):
        proposed_edges, label = _propose_flip(current, universe_list, rng)
        if label == "noop":
            continue

        proposed_logp = _log_posterior(
            proposed_edges, universe_list, cost_fn, cost_cache, beta, edge_penalty
        )

        # Metropolis-Hastings acceptance
        log_ratio = proposed_logp - current_logp
        if log_ratio >= 0.0 or rng.random() < math.exp(log_ratio):
            current = proposed_edges
            current_logp = proposed_logp
            if iteration >= n_burnin:
                acceptance_count += 1

        # Track posterior
        if iteration >= n_burnin:
            for e in current:
                edge_counts[e.edge_id] += 1
            n_edges_trace.append(len(current))

            if current_logp > best_logp:
                best_logp = current_logp
                best_edges = list(current)

    n_post_samples = n_samples if n_samples > 0 else 1
    edge_probabilities = {
        eid: count / n_post_samples for eid, count in edge_counts.items()
    }
    edge_log_odds = {
        eid: math.log(max(1e-12, p / (1.0 - p + 1e-12))) for eid, p in edge_probabilities.items()
    }

    # Entropy over the posterior (approximated from edge marginal probs)
    entropy = 0.0
    for p in edge_probabilities.values():
        if 0.0 < p < 1.0:
            entropy -= p * math.log(p) + (1.0 - p) * math.log(1.0 - p)

    n_edges_mean = float(np.mean(n_edges_trace)) if n_edges_trace else 0.0
    if n_edges_trace:
        sorted_trace = sorted(n_edges_trace)
        lo_idx = int(0.025 * len(sorted_trace))
        hi_idx = int(0.975 * len(sorted_trace))
        n_edges_ci95 = (float(sorted_trace[max(0, lo_idx)]), float(sorted_trace[min(len(sorted_trace) - 1, hi_idx)]))
    else:
        n_edges_ci95 = (0.0, 0.0)

    acceptance_rate = acceptance_count / max(1, n_samples)

    logger.info(
        "Topology posterior: %d edges, mean edges=%.1f, entropy=%.2f, accept=%.3f",
        len(universe_list),
        n_edges_mean,
        entropy,
        acceptance_rate,
    )

    return {
        "edge_probabilities": edge_probabilities,
        "edge_log_odds": edge_log_odds,
        "map_edges": [e.edge_id for e in best_edges],
        "entropy": entropy,
        "n_edges_mean": n_edges_mean,
        "n_edges_ci95": n_edges_ci95,
        "acceptance_rate": acceptance_rate,
    }


def _make_sheaf_cost_fn(
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
    node_info: Mapping[str, object],
    stats: object,
) -> Callable[[Sequence[Edge]], float]:
    """Build a cost function that evaluates sheaf global cost for a set of edges.

    This uses the same scoring as refine_edges_with_sheaf to maintain consistency.
    """
    from ..sheaf.topology_refine import _score_candidates
    from ..sheaf.directed_section import build_edge_maps, solve_directed_section

    def cost_fn(edges: Sequence[Edge]) -> float:
        if not edges:
            return 0.0

        # Score candidates using the same scoring function
        scores = _score_candidates(edges, node_info, sample_map, stats, config)
        local_cost = sum(s.local_score for s in scores.values())

        # Add global section residual if chemistry data is available
        node_vectors = {}
        for node_id, sample in sample_map.items():
            from ..data.schema import vector_from_sample

            values, _ = vector_from_sample(
                sample,
                config.ion_order,
                config.missing_policy,
                config.detection_limit_policy,
            )
            if values is not None:
                node_vectors[node_id] = values

        global_weight = float(getattr(config, "sheaf_weight_global", 1.0))

        if node_vectors:
            try:
                maps = build_edge_maps(
                    edges, node_vectors, config,
                    prior_weight=float(getattr(config, "sheaf_weight_head_prior", 1.0)),
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
                global_cost = 0.0
        else:
            global_cost = 0.0

        return local_cost + global_weight * global_cost

    return cost_fn


def attach_posterior_attrs(
    selected_edges: List[Edge],
    candidate_edges: Sequence[Edge],
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
        attrs = dict(edge.attrs or {})
        attrs["posterior_edge_probability"] = probs.get(edge.edge_id)
        attrs["posterior_edge_log_odds"] = log_odds.get(edge.edge_id)
        attrs["posterior_map_selected"] = edge.edge_id in map_edges
        attrs["posterior_topology_entropy"] = entropy
        attrs["posterior_n_edges_mean"] = n_mean
        attrs["posterior_n_edges_ci95"] = f"{ci[0]:.2f}-{ci[1]:.2f}"
        attrs["posterior_acceptance_rate"] = accept
        edge.attrs = attrs
