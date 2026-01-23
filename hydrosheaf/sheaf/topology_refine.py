"""Sheaf-inspired topology refinement using isotopes and Cl."""

from dataclasses import dataclass
import math
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import Config
from ..data.schema import parse_numeric
from ..graph.types import Edge
from ..isotopes import extract_isotopes, isotope_penalty
from .isotope_metrics import (
    IsotopeStats,
    compute_evaporation_probability,
    compute_isotope_stats,
    sample_depth_m,
)
from .solve import compute_edge_residuals, solve_isotope_field


@dataclass
class NodeIsotopeInfo:
    node_id: str
    d18o: Optional[float]
    d2h: Optional[float]
    d_excess: Optional[float]
    evap_index: Optional[float]
    cl: Optional[float]
    depth_m: Optional[float]
    p_evap: float


@dataclass
class EdgeSheafScore:
    edge: Edge
    prior_penalty: float
    iso_cost: float
    cl_cost: float
    cons_cost: float
    pi_evap: float
    flags: List[str]

    @property
    def local_score(self) -> float:
        return self.prior_penalty + self.iso_cost + self.cl_cost


def _sample_map(samples: object) -> Dict[str, Mapping[str, object]]:
    if isinstance(samples, Mapping):
        return {str(k): v for k, v in samples.items()}
    if isinstance(samples, Sequence):
        mapping: Dict[str, Mapping[str, object]] = {}
        for row in samples:
            site_id = row.get("site_id") or row.get("node_id") or row.get("sample_id")
            if site_id is None:
                continue
            mapping[str(site_id)] = row
        return mapping
    raise TypeError("Unsupported samples input type.")


def _log_mix_cost(cost_a: float, cost_b: float, weight_b: float) -> float:
    weight_b = min(1.0, max(0.0, weight_b))
    weight_a = 1.0 - weight_b
    if weight_b == 0.0:
        return cost_a
    if weight_a == 0.0:
        return cost_b
    min_cost = min(cost_a, cost_b)
    exp_a = math.exp(-(cost_a - min_cost))
    exp_b = math.exp(-(cost_b - min_cost))
    return -math.log(weight_a * exp_a + weight_b * exp_b) + min_cost


def _edge_prior_penalty(edge: Edge, weight: float) -> float:
    attrs = edge.attrs or {}
    prior = attrs.get("edge_confidence", attrs.get("p_uv", 1.0))
    p_uv = 1.0
    if isinstance(prior, (int, float, str)):
        try:
            p_uv = float(prior)
        except (TypeError, ValueError):
            p_uv = 1.0
    p_uv = min(1.0, max(1e-12, p_uv))
    return weight * (-math.log(p_uv))


def _edge_cl_cost(
    cl_u: Optional[float],
    cl_v: Optional[float],
    pi_evap: float,
) -> Tuple[float, Optional[float]]:
    if cl_u is None or cl_v is None:
        return 0.0, None
    if cl_u <= 0 or cl_v <= 0:
        return 0.0, None
    ratio = cl_v / cl_u
    if ratio <= 0:
        return 0.0, None
    log_ratio = math.log(ratio)
    cons_cost = abs(log_ratio)
    evap_cost = max(0.0, -log_ratio)
    return (1.0 - pi_evap) * cons_cost + pi_evap * evap_cost, ratio


def _edge_iso_cost(
    node_u: NodeIsotopeInfo,
    node_v: NodeIsotopeInfo,
    stats: IsotopeStats,
    config: Config,
) -> Tuple[float, float, float, List[str]]:
    flags: List[str] = []
    if node_u.d18o is None or node_u.d2h is None:
        flags.append("iso_missing_u")
    if node_v.d18o is None or node_v.d2h is None:
        flags.append("iso_missing_v")
    if flags:
        return 0.0, 0.0, 0.0, flags

    assert node_u.d18o is not None
    assert node_u.d2h is not None
    assert node_v.d18o is not None
    assert node_v.d2h is not None

    sigma_d18o = float(getattr(config, "sheaf_iso_sigma_d18o", 0.2))
    sigma_d2h = float(getattr(config, "sheaf_iso_sigma_d2h", 1.0))

    cons_cost = ((node_v.d18o - node_u.d18o) / sigma_d18o) ** 2
    cons_cost += ((node_v.d2h - node_u.d2h) / sigma_d2h) ** 2

    pi_evap = 0.5 * (node_u.p_evap + node_v.p_evap)
    shallow_depth = float(getattr(config, "sheaf_shallow_depth_m", 30.0))
    if node_v.depth_m is not None and shallow_depth > 0 and node_v.depth_m > shallow_depth:
        pi_evap *= shallow_depth / node_v.depth_m
    if node_u.evap_index is not None and node_v.evap_index is not None:
        if abs(node_v.evap_index) <= abs(node_u.evap_index) + 1e-6:
            pi_evap *= 0.5
    if node_u.d_excess is not None and node_v.d_excess is not None:
        if node_v.d_excess > node_u.d_excess + 0.5:
            pi_evap *= 0.5
    pi_evap = min(0.8, max(0.0, pi_evap))
    if pi_evap > 0.3:
        flags.append("evap_candidate")

    evap_penalty, _ = isotope_penalty(
        node_u.d18o,
        node_u.d2h,
        node_v.d18o,
        node_v.d2h,
        config.lmwl_a,
        config.lmwl_b,
        "evap",
        d_excess_weight=config.isotope_d_excess_weight,
    )
    iso_cost = _log_mix_cost(cons_cost, evap_penalty, pi_evap)
    return iso_cost, cons_cost, pi_evap, flags


def _build_node_info(
    sample_map: Mapping[str, Mapping[str, object]],
    stats: IsotopeStats,
    config: Config,
) -> Dict[str, NodeIsotopeInfo]:
    node_info: Dict[str, NodeIsotopeInfo] = {}
    for node_id, sample in sample_map.items():
        isotopes = extract_isotopes(
            sample,
            d18o_key=config.isotope_d18o_key,
            d2h_key=config.isotope_d2h_key,
        )
        d18o = d2h = None
        d_excess = evap_index = None
        if isotopes is not None:
            d18o, d2h = isotopes
            d_excess = d2h - 8.0 * d18o
            evap_index = d2h - (config.lmwl_a + config.lmwl_b * d18o)

        cl_val = parse_numeric(sample.get("Cl"), config.detection_limit_policy)
        depth_m = sample_depth_m(sample, config)
        p_evap = compute_evaporation_probability(evap_index, d_excess, depth_m, stats, config)

        node_info[node_id] = NodeIsotopeInfo(
            node_id=node_id,
            d18o=d18o,
            d2h=d2h,
            d_excess=d_excess,
            evap_index=evap_index,
            cl=cl_val,
            depth_m=depth_m,
            p_evap=p_evap,
        )
    return node_info


def _score_candidates(
    candidates: Iterable[Edge],
    node_info: Mapping[str, NodeIsotopeInfo],
    stats: IsotopeStats,
    config: Config,
) -> Dict[str, EdgeSheafScore]:
    scores: Dict[str, EdgeSheafScore] = {}
    weight_head = float(getattr(config, "sheaf_weight_head_prior", 1.0))
    weight_iso = float(getattr(config, "sheaf_weight_isotope", 1.0))
    weight_cl = float(getattr(config, "sheaf_weight_cl", 0.5))

    for edge in candidates:
        node_u = node_info.get(edge.u)
        node_v = node_info.get(edge.v)
        if node_u is None or node_v is None:
            continue

        prior_penalty = _edge_prior_penalty(edge, weight_head)
        iso_cost = 0.0
        cons_cost = 0.0
        pi_evap = 0.0
        flags: List[str] = []

        if getattr(config, "sheaf_isotope_enabled", True):
            iso_cost, cons_cost, pi_evap, flags = _edge_iso_cost(
                node_u, node_v, stats, config
            )
            iso_cost *= weight_iso

        cl_cost = 0.0
        cl_ratio = None
        if getattr(config, "sheaf_cl_enabled", True):
            cl_cost, cl_ratio = _edge_cl_cost(node_u.cl, node_v.cl, pi_evap)
            cl_cost *= weight_cl
            if cl_ratio is None:
                flags.append("cl_missing")

        scores[edge.edge_id] = EdgeSheafScore(
            edge=edge,
            prior_penalty=prior_penalty,
            iso_cost=iso_cost,
            cl_cost=cl_cost,
            cons_cost=cons_cost,
            pi_evap=pi_evap,
            flags=flags,
        )
    return scores


def _select_by_score(
    scores: Mapping[str, EdgeSheafScore],
    candidate_groups: Mapping[str, List[str]],
    max_neighbors: int,
    global_residuals: Optional[Mapping[str, float]] = None,
    global_weight: float = 0.0,
) -> List[Edge]:
    selected: List[Edge] = []
    for _, edge_ids in candidate_groups.items():
        scored = []
        for edge_id in edge_ids:
            score_obj = scores.get(edge_id)
            if score_obj is None:
                continue
            residual = 0.0
            if global_residuals is not None:
                residual = float(global_residuals.get(edge_id, 0.0))
            total = score_obj.local_score + global_weight * residual
            scored.append((total, score_obj.edge))
        scored.sort(key=lambda item: item[0])
        for _, edge in scored[: max_neighbors or 0]:
            selected.append(edge)
    return selected


def refine_edges_with_sheaf(
    samples: object,
    candidates: Iterable[Edge],
    config: Config,
) -> List[Edge]:
    sample_map = _sample_map(samples)
    candidate_list = list(candidates)
    if not candidate_list:
        return []

    stats = compute_isotope_stats(sample_map.values(), config)
    node_info = _build_node_info(sample_map, stats, config)
    scores = _score_candidates(candidate_list, node_info, stats, config)

    grouped: Dict[str, List[str]] = {}
    for edge in candidate_list:
        grouped.setdefault(edge.u, []).append(edge.edge_id)

    max_neighbors = int(getattr(config, "edge_max_neighbors", 1))
    selected = _select_by_score(scores, grouped, max_neighbors)

    global_weight = float(getattr(config, "sheaf_weight_global", 1.0))
    max_iter = int(getattr(config, "sheaf_max_iter", 3))
    sigma_d18o = float(getattr(config, "sheaf_iso_sigma_d18o", 0.2))
    sigma_d2h = float(getattr(config, "sheaf_iso_sigma_d2h", 1.0))

    has_isotopes = any(
        info.d18o is not None and info.d2h is not None
        for info in node_info.values()
    )
    node_estimates: Dict[str, Tuple[float, float]] = {}

    for _ in range(max_iter if has_isotopes else 0):
        edge_weights: Dict[str, float] = {}
        for edge in selected:
            score_obj = scores.get(edge.edge_id)
            if score_obj is None:
                continue
            weight = max(0.0, 1.0 - score_obj.pi_evap)
            weight = weight / (1.0 + score_obj.cons_cost)
            edge_weights[edge.edge_id] = weight

        node_obs = {
            node_id: (info.d18o, info.d2h) for node_id, info in node_info.items()
        }
        node_estimates = solve_isotope_field(
            node_info.keys(),
            selected,
            node_obs,
            edge_weights,
            obs_weight=1.0,
            diag_eps=1e-6,
        )
        residuals = compute_edge_residuals(
            candidate_list, node_estimates, sigma_d18o, sigma_d2h
        )

        updated = _select_by_score(
            scores, grouped, max_neighbors, residuals, global_weight
        )
        if {edge.edge_id for edge in updated} == {edge.edge_id for edge in selected}:
            selected = updated
            break
        selected = updated

    final_residuals: Dict[str, float] = {}
    if has_isotopes and node_estimates:
        final_residuals = compute_edge_residuals(
            selected, node_estimates, sigma_d18o, sigma_d2h
        )

    for edge in selected:
        score_obj = scores.get(edge.edge_id)
        if score_obj is None:
            continue
        attrs = dict(edge.attrs or {})
        attrs["sheaf_score_local"] = score_obj.local_score
        attrs["sheaf_cost_iso"] = score_obj.iso_cost
        attrs["sheaf_cost_cl"] = score_obj.cl_cost
        attrs["sheaf_pi_evap"] = score_obj.pi_evap
        attrs["sheaf_flags"] = ",".join(score_obj.flags)
        attrs["sheaf_residual_global"] = final_residuals.get(edge.edge_id)
        attrs["sheaf_score_global"] = score_obj.local_score + global_weight * float(
            final_residuals.get(edge.edge_id, 0.0)
        )
        edge.attrs = attrs

    return selected
