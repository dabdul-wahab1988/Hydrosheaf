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
        attrs = _get_edge_attrs(e)
        _p = attrs.get("edge_confidence", attrs.get("p_uv", None))
        p = float(_p) if _p is not None else 0.5
        p = min(1.0 - 1e-12, max(1e-12, p))
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
    """Propose adding or removing a random edge from universe.

    Returns (proposed_edge_list, proposal_label).
    """
    current_set = {_get_edge_id(e) for e in current}
    universe_list = list(universe)

    # Randomly choose add or remove
    can_remove = len(current_set) > 0
    can_add = len(current_set) < len(universe_list)

    if can_remove and (not can_add or rng.random() < 0.5):
        # Remove a random edge from current
        removal = rng.choice(list(current))
        proposed = [e for e in current if _get_edge_id(e) != _get_edge_id(removal)]
        return proposed, f"remove_{_get_edge_id(removal)}"
    elif can_add:
        # Add a random edge not in current
        absent = [e for e in universe_list if _get_edge_id(e) not in current_set]
        addition = rng.choice(absent)
        proposed = list(current) + [addition]
        return proposed, f"add_{_get_edge_id(addition)}"
    else:
        return list(current), "noop"


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
            Posterior entropy over graphs.
        n_edges_mean : float
            Posterior mean number of edges.
        n_edges_ci95 : Tuple[float, float]
            95% credible interval for number of edges.
        acceptance_rate : float
            M-H acceptance rate.
    """
    n_samples = _get_config_int(config, "topology_posterior_samples", 2000)
    n_burnin = _get_config_int(config, "topology_posterior_burnin", 500)
    beta = _get_config_float(config, "topology_posterior_beta", 1.0)
    edge_penalty = _get_config_float(config, "topology_posterior_edge_penalty", 0.0)

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
    edge_counts: Dict[str, int] = {_get_edge_id(e): 0 for e in universe_list}
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
                edge_counts[_get_edge_id(e)] += 1
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
        "map_edges": [_get_edge_id(e) for e in best_edges],
        "entropy": entropy,
        "n_edges_mean": n_edges_mean,
        "n_edges_ci95": n_edges_ci95,
        "acceptance_rate": acceptance_rate,
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

    def cost_fn(edges: Sequence[Any]) -> float:
        if not edges:
            return 0.0

        # Score candidates using the same scoring function
        scores = _score_candidates(edges, node_info, sample_map, stats, config)
        local_cost = sum(s.local_score for s in scores.values())

        if node_vectors:
            try:
                maps = build_edge_maps(
                    edges, node_vectors, config,
                    prior_weight=_get_config_float(config, "sheaf_weight_head_prior", 1.0),
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
        attrs["posterior_n_edges_mean"] = n_mean
        attrs["posterior_n_edges_ci95"] = f"{ci[0]:.2f}-{ci[1]:.2f}"
        attrs["posterior_acceptance_rate"] = accept
        edge.attrs = attrs
