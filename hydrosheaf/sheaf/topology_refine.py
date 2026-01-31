"""Sheaf-inspired topology refinement using isotopes and Cl."""

from dataclasses import dataclass
import math
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import Config
from ..log import get_logger
from ..data.schema import parse_numeric, vector_from_sample

logger = get_logger("sheaf.topology_refine")

from ..graph.types import Edge
from ..isotopes import extract_isotopes, isotope_penalty
from .directed_section import (
    build_edge_maps,
    compute_edge_section_residuals,
    solve_directed_section,
)
from .isotope_metrics import (
    IsotopeStats,
    compute_evaporation_probability,
    compute_isotope_stats,
    sample_depth_m,
)

try:
    from hydrosheaf.nuclear.invert import infer_age_from_tracer
    from hydrosheaf.nuclear.nuclides import get_nuclide, TRITIUM
    _NUCLEAR_AVAILABLE = True
except ImportError:
    _NUCLEAR_AVAILABLE = False
    infer_age_from_tracer = None
    get_nuclide = None
    TRITIUM = None


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
    age_years: Optional[float] = None
    age_sigma_years: Optional[float] = None



@dataclass
class EdgeSheafScore:
    edge: Edge
    prior_penalty: float
    iso_cost: float
    cl_cost: float
    age_cost: float
    cons_cost: float
    pi_evap: float
    flags: List[str]

    @property
    def local_score(self) -> float:
        return self.prior_penalty + self.iso_cost + self.cl_cost + self.age_cost



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


def _edge_age_cost(
    node_u: NodeIsotopeInfo,
    node_v: NodeIsotopeInfo,
) -> Tuple[float, List[str]]:
    flags: List[str] = []
    if node_u.age_years is None or node_v.age_years is None:
        return 0.0, flags
        
    diff = node_v.age_years - node_u.age_years
    
    # Calculate uncertainty
    sigma_u = node_u.age_sigma_years or 1.0
    sigma_v = node_v.age_sigma_years or 1.0
    sigma_diff = math.sqrt(sigma_u**2 + sigma_v**2)
    
    # Check if age decreases significantly
    # User said: penalize if A_v < A_u beyond uncertainty
    # A_v - A_u >= -k * sigma -> Cost = 0
    # Otherwise -> Cost = ((A_u - A_v) / sigma)^2
    
    k = 1.0 # Tolerance factor (1 sigma)
    
    if diff >= -k * sigma_diff:
        return 0.0, flags
    
    cost = ((node_u.age_years - node_v.age_years) / sigma_diff) ** 2
    flags.append("age_reversal")
    
    return cost, flags



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

    # Compute evaporation probability pi_evap based on physical proxies
    # Base pi_evap is derived from isotope stats (LMWL departure)
    pi_evap = 0.5 * (node_u.p_evap + node_v.p_evap)
    
    # Physical Heuristic 1: Depth constraint
    # Evaporation occurs at or near the surface. Groundwater at depth (>30m) 
    # is unlikely to experience direct evaporative enrichment unless it was 
    # recharged recently through an open surface water body.
    shallow_depth = float(getattr(config, "sheaf_shallow_depth_m", 30.0))
    if node_v.depth_m is not None and shallow_depth > 0 and node_v.depth_m > shallow_depth:
        pi_evap *= shallow_depth / node_v.depth_m
        
    # Physical Heuristic 2: Thermodynamic Directionality
    # Evaporation always increases enrichment (more positive delta values) 
    # and increases the "Evaporation Index" (departure from LMWL).
    # If the downstream sample is LESS enriched/deviated than upstream, 
    # the process is likely mixing/dilution, not evaporation.
    if node_u.evap_index is not None and node_v.evap_index is not None:
        if abs(node_v.evap_index) <= abs(node_u.evap_index) + 1e-6:
            pi_evap *= 0.5
            
    # Physical Heuristic 3: d-Excess Signature
    # Kinetic fractionation during evaporation uniquely decreases d-excess.
    # If d-excess increases along a flow path, it strongly suggests a different 
    # water source (e.g. mountain-front recharge or different precip regime).
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
    import datetime
    
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

        # Nuclear age inference
        age_years = None
        age_sigma = None
        
        # Check for pre-calculated age
        if "mean_age_years" in sample:
             age_years = parse_numeric(sample["mean_age_years"], config.detection_limit_policy)
        
        # If no pre-calculated age, try to infer from Tritium
        if age_years is None and _NUCLEAR_AVAILABLE and infer_age_from_tracer is not None and get_nuclide is not None:
            t_keys = ["3H", "tritium", "H3", "H3_TU", "tritium_TU", "Tritium"]
            t_val = None
            for k in t_keys:
                if k in sample:
                    t_val = parse_numeric(sample[k], config.detection_limit_policy)
                    if t_val is not None: break
            
            if t_val is not None:
                # Try to get date
                date_val = sample.get("date") or sample.get("timestamp") or sample.get("sample_date")
                year_fraction = None
                if date_val:
                    # Simple parser (assume YYYY-MM-DD or similar string, or year float)
                    try:
                        if isinstance(date_val, (int, float)):
                            if 1900 < date_val < 2100:
                                year_fraction = float(date_val)
                        else:
                            # Assume ISO format YYYY-MM-DD
                            d_str = str(date_val).strip()
                            if len(d_str) >= 4:
                                dt = datetime.datetime.fromisoformat(d_str.replace("Z", "+00:00"))
                                start = datetime.datetime(dt.year, 1, 1, tzinfo=dt.tzinfo)
                                year_fraction = dt.year + (dt - start).total_seconds() / (365.25 * 86400)
                    except Exception:
                        pass
                
                if year_fraction is not None:
                    # Infer
                    sigma_val = max(0.5, t_val * 0.1)
                    try:
                        nuclide = get_nuclide("3H") or TRITIUM
                        if nuclide:
                            res = infer_age_from_tracer(
                                t_val, sigma_val, year_fraction, 
                                nuclide=nuclide, 
                                model="PFM"
                            )
                            age_years = res["tau_map_years"]
                            age_sigma = (res["tau_ci_high_years"] - res["tau_ci_low_years"]) / 4.0
                    except Exception:
                        pass

        node_info[node_id] = NodeIsotopeInfo(
            node_id=node_id,
            d18o=d18o,
            d2h=d2h,
            d_excess=d_excess,
            evap_index=evap_index,
            cl=cl_val,
            depth_m=depth_m,
            p_evap=p_evap,
            age_years=age_years,
            age_sigma_years=age_sigma,
        )
    return node_info



def _build_node_vectors(
    sample_map: Mapping[str, Mapping[str, object]],
    config: Config,
) -> Dict[str, Optional[List[float]]]:
    node_vectors: Dict[str, Optional[List[float]]] = {}
    for node_id, sample in sample_map.items():
        values, _ = vector_from_sample(
            sample,
            config.ion_order,
            config.missing_policy,
            config.detection_limit_policy,
        )
        node_vectors[node_id] = values
    return node_vectors


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
    weight_age = float(getattr(config, "sheaf_weight_age", 2.0))

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
        
        age_cost = 0.0
        if getattr(config, "sheaf_age_enabled", True):
            age_c, age_flags = _edge_age_cost(node_u, node_v)
            age_cost = age_c * weight_age
            if age_flags:
                flags.extend(age_flags)

        scores[edge.edge_id] = EdgeSheafScore(
            edge=edge,
            prior_penalty=prior_penalty,
            iso_cost=iso_cost,
            cl_cost=cl_cost,
            age_cost=age_cost,
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
    soft_beta: float = 2.0,
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
        
        if not scored:
            continue

        # Softmax weighting (lower score is better)
        # We negate scores and multiply by beta (sharpness)
        vals = [-s * soft_beta for s, _ in scored]
        max_val = max(vals)
        exps = [math.exp(v - max_val) for v in vals]
        sum_exps = sum(exps)
        
        for (score, edge), weight_factor in zip(scored, exps):
            # Normalized weight
            prob = weight_factor / sum_exps if sum_exps > 0 else 0.0
            
            # Update edge attributes for the solver
            attrs = edge.attrs or {}
            attrs["sheaf_weight"] = prob
            # Update edge_confidence so build_edge_maps uses it naturally
            attrs["edge_confidence"] = prob
            edge.attrs = attrs
            
            # Keep edge if it has significant weight or if we want to keep all
            # Pruning very small weights helps performance
            if prob > 0.001:
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
    # Higher beta means sharper selection (closer to argmax)
    soft_beta = float(getattr(config, "sheaf_soft_beta", 2.0))

    selected = _select_by_score(
        scores, 
        grouped, 
        max_neighbors, 
        soft_beta=soft_beta
    )

    global_weight = float(getattr(config, "sheaf_weight_global", 1.0))
    max_iter = int(getattr(config, "sheaf_max_iter", 3))

    node_vectors = _build_node_vectors(sample_map, config)
    node_ids = list(
        {node_id for node_id in node_vectors}
        | {edge.u for edge in candidate_list}
        | {edge.v for edge in candidate_list}
    )

    has_chemistry = any(values is not None for values in node_vectors.values())
    node_estimates: Dict[str, List[float]] = {
        node_id: values
        for node_id, values in node_vectors.items()
        if values is not None
    }

    for iter_idx in range(max_iter if has_chemistry else 0):
        logger.info(f"Sheaf Refinement Iteration {iter_idx + 1}/{max_iter}")
        edge_maps = build_edge_maps(
            candidate_list, # Always use all candidates for building map potential
            node_estimates,
            config,
            prior_weight=float(getattr(config, "sheaf_weight_head_prior", 1.0)),
        )
        if not edge_maps:
            break
            
        # Filter edge_maps by 'selected' logic (which now contains all weighted edges > 0.001)
        selected_ids = {e.edge_id for e in selected}
        selected_maps = [m for eid, m in edge_maps.items() if eid in selected_ids]
        
        if not selected_maps:
            break

        node_estimates = solve_directed_section(
            node_ids,
            selected_maps,
            node_vectors,
            obs_weight=1.0,
            diag_eps=1e-6,
        )
        
        # Log global residual energy for tracking convergence
        current_energy = sum(
            (node_estimates[nid][d] - (val[d] if val else 0.0))**2 
            for nid, val in node_vectors.items() 
            if val is not None
            for d in range(len(val))
        )
        logger.science(f"Iter {iter_idx+1}: Global Section Energy = {current_energy:.4f}")

        edge_maps = build_edge_maps(
            candidate_list,
            node_estimates,
            config,
            prior_weight=float(getattr(config, "sheaf_weight_head_prior", 1.0)),
        )
        residuals = compute_edge_section_residuals(
            edge_maps, node_estimates, config.weights
        )

        updated = _select_by_score(
            scores, 
            grouped, 
            max_neighbors, 
            residuals, 
            global_weight,
            soft_beta=soft_beta
        )
        
        selected = updated

    final_residuals: Dict[str, float] = {}
    if has_chemistry and node_estimates:
        edge_maps = build_edge_maps(
            selected,
            node_estimates,
            config,
            prior_weight=float(getattr(config, "sheaf_weight_head_prior", 1.0)),
        )
        final_residuals = compute_edge_section_residuals(
            edge_maps, node_estimates, config.weights
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
