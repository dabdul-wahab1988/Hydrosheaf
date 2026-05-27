"""Sheaf-inspired topology refinement using isotopes and Cl."""

from dataclasses import dataclass, field
import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from ..config import Config
from ..log import get_logger
from ..data.schema import parse_numeric, vector_from_sample

logger = get_logger("sheaf.topology_refine")


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

# Null-model imports (lazy-loaded to avoid circular dependency)
_null_models_imported = False
_compute_null_penalty = None

def _ensure_null_models():
    global _null_models_imported, _compute_null_penalty
    if not _null_models_imported:
        try:
            from hydrosheaf.null_models import compute_null_penalty as _cnp
            _compute_null_penalty = _cnp
        except ImportError:
            _compute_null_penalty = None
        _null_models_imported = True

from ..graph.types import Edge
from ..isotopes import extract_isotopes, isotope_penalty
from ..null_models import compute_null_penalty
from ..validation.evidence import EdgeEvidenceClass, classify_edge_evidence
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
    null_score: float = 0.0
    evidence_class: str = ""
    evidence_reason: str = ""
    ot_cost: float = 0.0
    ot_attrs: Dict[str, float] = field(default_factory=dict)
    causal_attrs: Dict[str, object] = field(default_factory=dict)

    @property
    def local_score(self) -> float:
        return self.prior_penalty + self.iso_cost + self.cl_cost + self.age_cost + self.ot_cost



def _extract_temporal_history(
    node_id: str,
    sample_map: Mapping[str, Mapping[str, object]],
) -> Optional[List[Mapping[str, object]]]:
    """Extract temporal history for a node from the sample map.

    Returns a list of samples for the given node if the sample_map contains
    temporal data (date/timestamp keys or multiple samples per site_id).
    Returns None when only a single static sample is available.
    """
    temporal_keys = {"date", "timestamp", "sample_date", "datetime"}
    sample = sample_map.get(node_id)
    if sample is None:
        return None

    # Check for nested temporal data under a "history" or "replicates" key
    for hist_key in ("history", "replicates", "time_series", "temporal"):
        hist = sample.get(hist_key)
        if isinstance(hist, list) and len(hist) >= 2:
            return hist

    # Check if the sample itself has temporal keys
    has_temporal = any(k in sample for k in temporal_keys)
    if has_temporal:
        return [sample]

    # Check if there are multiple samples with the same site_id prefix
    # (e.g., "well_A_2020-01", "well_A_2020-02")
    prefix = node_id.split("_")[0] if "_" in node_id else node_id
    matching = []
    for sid, s in sample_map.items():
        if sid == node_id or sid.startswith(prefix + "_"):
            if any(k in s for k in temporal_keys):
                matching.append(s)
    if len(matching) >= 2:
        return matching

    return None


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

    sigma_d18o = _get_config_float(config, "sheaf_iso_sigma_d18o", 0.2)
    sigma_d2h = _get_config_float(config, "sheaf_iso_sigma_d2h", 1.0)

    cons_cost = ((node_v.d18o - node_u.d18o) / sigma_d18o) ** 2
    cons_cost += ((node_v.d2h - node_u.d2h) / sigma_d2h) ** 2

    # Compute evaporation probability pi_evap based on physical proxies
    # Base pi_evap is derived from isotope stats (LMWL departure)
    pi_evap = 0.5 * (node_u.p_evap + node_v.p_evap)
    
    # Physical Heuristic 1: Depth constraint
    # Evaporation occurs at or near the surface. Groundwater at depth (>30m) 
    # is unlikely to experience direct evaporative enrichment unless it was 
    # recharged recently through an open surface water body.
    shallow_depth = _get_config_float(config, "sheaf_shallow_depth_m", 30.0)
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
    sample_map: Dict[str, Mapping[str, object]],
    stats: IsotopeStats,
    config: Config,
) -> Dict[str, EdgeSheafScore]:
    scores: Dict[str, EdgeSheafScore] = {}
    weight_head = _get_config_float(config, "sheaf_weight_head_prior", 1.0)
    weight_iso = _get_config_float(config, "sheaf_weight_isotope", 1.0)
    weight_cl = _get_config_float(config, "sheaf_weight_cl", 0.5)
    weight_age = _get_config_float(config, "sheaf_weight_age", 2.0)

    null_enabled = getattr(config, "null_model_enabled", False)
    evidence_enabled = getattr(config, "evidence_ladder_enabled", False)
    # Master gate: assumption_calibration_enabled enables both sub-features
    if getattr(config, "assumption_calibration_enabled", False):
        null_enabled = True
        evidence_enabled = True
    null_weight = _get_config_float(config, "null_model_weight", 0.5)

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

        # Track whether isotope/Cl data is missing for evidence classification
        iso_missing = False
        cl_missing = False

        if getattr(config, "sheaf_isotope_enabled", True):
            iso_cost, cons_cost, pi_evap, flags = _edge_iso_cost(
                node_u, node_v, stats, config
            )
            iso_cost *= weight_iso
            if "iso_missing_u" in flags or "iso_missing_v" in flags:
                iso_missing = True

        cl_cost = 0.0
        cl_ratio = None
        if getattr(config, "sheaf_cl_enabled", True):
            cl_cost, cl_ratio = _edge_cl_cost(node_u.cl, node_v.cl, pi_evap)
            cl_cost *= weight_cl
            if cl_ratio is None:
                flags.append("cl_missing")
                cl_missing = True
        
        age_cost = 0.0
        if getattr(config, "sheaf_age_enabled", True):
            age_c, age_flags = _edge_age_cost(node_u, node_v)
            age_cost = age_c * weight_age
            if age_flags:
                flags.extend(age_flags)

        # === Null-model integration (Phase 0-1) ===
        null_score = 0.0
        edge_evidence_class = ""
        evidence_reason = ""

        if null_enabled:
            sample_a = sample_map.get(edge.u, {})
            sample_b = sample_map.get(edge.v, {})
            try:
                null_score, null_flags = compute_null_penalty(sample_a, sample_b, config)
                # Downgrade chemistry/isotope support proportionally to null plausibility
                # Higher null_score -> larger added cost
                iso_cost += null_score * null_weight * weight_iso
                flags.extend(null_flags)
            except Exception:
                logger.warning(
                    "Null model failed for edge %s; flagging and continuing.",
                    edge.edge_id,
                )
                flags.append("null_error")

        # Optimal transport plausibility
        ot_cost = 0.0
        ot_attrs: Dict[str, float] = {}
        if getattr(config, "ot_enabled", False):
            try:
                from ..models.optimal_transport import compute_unbalanced_ot
                from ..data.schema import vector_from_sample

                sample_u = sample_map.get(edge.u, {})
                sample_v = sample_map.get(edge.v, {})
                x_u, _ = vector_from_sample(
                    sample_u, config.ion_order,
                    config.missing_policy, config.detection_limit_policy,
                )
                x_v, _ = vector_from_sample(
                    sample_v, config.ion_order,
                    config.missing_policy, config.detection_limit_policy,
                )
                if x_u is not None and x_v is not None:
                    ot_result = compute_unbalanced_ot(
                        x_u, x_v, config.ion_order, config
                    )
                    ot_weight_val = _get_config_float(config, "ot_weight", 0.25)
                    raw_ot_total = float(ot_result.get("ot_total_cost", 0.0))
                    ot_cost = ot_weight_val * raw_ot_total
                    ot_attrs = {
                        "ot_total_cost": raw_ot_total,
                        "ot_score_contribution": ot_cost,
                        "ot_balanced_cost": float(ot_result.get("ot_balanced_cost", 0.0)),
                        "ot_creation_mass": float(ot_result.get("ot_creation_mass", 0.0)),
                        "ot_destruction_mass": float(ot_result.get("ot_destruction_mass", 0.0)),
                        "ot_conservative_mismatch": float(ot_result.get("ot_conservative_mismatch", 0.0)),
                        "ot_reaction_plausibility": float(ot_result.get("ot_reaction_plausibility", 0.0)),
                    }
            except Exception:
                logger.debug("OT computation failed for edge %s", edge.edge_id)

        # Causal discovery support
        causal_support = 0.0
        causal_confound = 0.0
        causal_attrs: Dict[str, object] = {}
        if getattr(config, "causal_discovery_enabled", False):
            try:
                from ..causal.discovery import compute_causal_support

                sample_u = sample_map.get(edge.u, {})
                sample_v = sample_map.get(edge.v, {})

                # Build temporal histories from sample_map if available.
                # Samples with the same site_id prefix or containing a date/timestamp
                # key are grouped as time series.
                upstream_history = _extract_temporal_history(
                    edge.u, sample_map
                )
                downstream_history = _extract_temporal_history(
                    edge.v, sample_map
                )

                causal_result = compute_causal_support(
                    sample_u, sample_v,
                    upstream_history=upstream_history,
                    downstream_history=downstream_history,
                    config=config,
                )
                causal_support = float(causal_result.get("causal_support_score", 0.0))
                causal_confound = float(causal_result.get("causal_confounded_score", 0.0))
                causal_status = str(causal_result.get("causal_status", ""))
                causal_weight_val = _get_config_float(config, "causal_weight", 0.25)
                causal_attrs = causal_result
                # Support reduces cost; only apply confound penalty when
                # we have enough data to meaningfully assess it.
                iso_cost -= causal_weight_val * causal_support
                if causal_status != "insufficient_data":
                    iso_cost += causal_weight_val * causal_confound
            except Exception:
                logger.debug("Causal discovery failed for edge %s", edge.edge_id)

        if evidence_enabled:
            if iso_missing or cl_missing:
                edge_evidence_class = EdgeEvidenceClass.AMBIGUOUS.name
                evidence_reason = "missing key isotope or chloride data"
                flags.append("missing_evidence")
            else:
                edge_evidence_class, evidence_reason = classify_edge_evidence(
                    local_score=prior_penalty + iso_cost + cl_cost + age_cost,
                    null_score=null_score,
                    flags=flags,
                    config=config,
                )

        scores[edge.edge_id] = EdgeSheafScore(
            edge=edge,
            prior_penalty=prior_penalty,
            iso_cost=iso_cost,
            cl_cost=cl_cost,
            age_cost=age_cost,
            cons_cost=cons_cost,
            pi_evap=pi_evap,
            flags=flags,
            null_score=null_score,
            evidence_class=edge_evidence_class,
            evidence_reason=evidence_reason,
            ot_cost=ot_cost,
            ot_attrs=ot_attrs,
            causal_attrs=causal_attrs,
        )
    return scores



def _select_by_score(
    scores: Mapping[str, EdgeSheafScore],
    candidate_groups: Mapping[str, List[str]],
    max_neighbors: int,
    global_residuals: Optional[Mapping[str, float]] = None,
    global_weight: float = 0.0,
    soft_beta: float = 2.0,
    selection_penalties: Optional[Mapping[str, float]] = None,
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
            penalty = 0.0
            if selection_penalties is not None:
                penalty = float(selection_penalties.get(edge_id, 0.0))
            total = score_obj.local_score + global_weight * residual + penalty
            scored.append((total, score_obj.edge))
        
        if not scored:
            continue

        # Softmax weighting (lower score is better)
        # We negate scores and multiply by beta (sharpness)
        vals = [-s * soft_beta for s, _ in scored]
        max_val = max(vals)
        exps = [math.exp(v - max_val) for v in vals]
        sum_exps = sum(exps)
        
        weighted_edges: List[Tuple[float, float, Edge]] = []
        for (score, edge), weight_factor in zip(scored, exps):
            # Normalized weight
            prob = weight_factor / sum_exps if sum_exps > 0 else 0.0
            
            # Update edge attributes for the solver
            attrs = edge.attrs or {}
            attrs["sheaf_weight"] = prob
            # Update edge_confidence so build_edge_maps uses it naturally
            attrs["edge_confidence"] = prob
            edge.attrs = attrs

            weighted_edges.append((prob, score, edge))

        keep_n = max(0, int(max_neighbors))
        if keep_n <= 0:
            continue

        weighted_edges.sort(key=lambda item: (-item[0], item[1]))
        for prob, _, edge in weighted_edges[:keep_n]:
            if prob <= 0.0:
                continue
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
    scores = _score_candidates(candidate_list, node_info, sample_map, stats, config)

    grouped: Dict[str, List[str]] = {}
    for edge in candidate_list:
        grouped.setdefault(edge.u, []).append(edge.edge_id)

    max_neighbors = _get_config_int(config, "edge_max_neighbors", 1)
    # Higher beta means sharper selection (closer to argmax)
    soft_beta = _get_config_float(config, "sheaf_soft_beta", 2.0)
    use_hydraulic_hodge = getattr(config, "hydraulic_hodge_enabled", False)

    reference_distance_km = None
    hodge_local_residuals: Optional[Dict[str, float]] = None
    if use_hydraulic_hodge:
        try:
            from .hydraulic_hodge import (
                compute_local_head_plane_residuals,
                extract_node_heads,
                infer_reference_distance_km,
            )

            reference_distance_km = infer_reference_distance_km(
                candidate_list,
                sample_map,
                config,
            )
            if getattr(config, "head_plane_residual_enabled", False):
                _head_map = extract_node_heads(sample_map, config)
                _n_neighbors = _get_config_int(config, "head_plane_residual_neighbors", 8)
                hodge_local_residuals = compute_local_head_plane_residuals(
                    _head_map, sample_map, n_neighbors=_n_neighbors
                )
        except Exception:
            logger.warning(
                "Hydraulic Hodge reference distance failed; continuing without it.",
                exc_info=True,
            )

    selected = _select_by_score(
        scores, 
        grouped, 
        max_neighbors, 
        soft_beta=soft_beta
    )

    global_weight = _get_config_float(config, "sheaf_weight_global", 1.0)
    max_iter = _get_config_int(config, "sheaf_max_iter", 3)

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

    iter_count = max_iter if (has_chemistry or use_hydraulic_hodge) else 0
    for iter_idx in range(iter_count):
        logger.info(f"Sheaf Refinement Iteration {iter_idx + 1}/{max_iter}")
        residuals: Dict[str, float] = {}
        if has_chemistry:
            edge_maps = build_edge_maps(
                candidate_list, # Always use all candidates for building map potential
                node_estimates,
                config,
                prior_weight=_get_config_float(config, "sheaf_weight_head_prior", 1.0),
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
                prior_weight=_get_config_float(config, "sheaf_weight_head_prior", 1.0),
            )
            residuals = compute_edge_section_residuals(
                edge_maps, node_estimates, config.weights
            )

        head_penalties: Dict[str, float] = {}
        if use_hydraulic_hodge and selected:
            try:
                from .hydraulic_hodge import compute_head_hodge_edge_penalties

                head_penalties = compute_head_hodge_edge_penalties(
                    selected,
                    sample_map,
                    config,
                    reference_distance_km=reference_distance_km,
                    local_residuals=hodge_local_residuals,
                )
            except Exception:
                logger.warning(
                    "Hydraulic Hodge leverage feedback failed; continuing.",
                    exc_info=True,
                )

        updated = _select_by_score(
            scores, 
            grouped, 
            max_neighbors, 
            residuals, 
            global_weight,
            soft_beta=soft_beta,
            selection_penalties=head_penalties,
        )
        
        selected = updated

    final_residuals: Dict[str, float] = {}
    if has_chemistry and node_estimates:
        edge_maps = build_edge_maps(
            selected,
            node_estimates,
            config,
            prior_weight=_get_config_float(config, "sheaf_weight_head_prior", 1.0),
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
        # Evidence ladder (Phase 0-1)
        if score_obj.evidence_class:
            attrs["evidence_class"] = score_obj.evidence_class
        attrs["evidence_score"] = 1.0 / (1.0 + score_obj.local_score)
        attrs["null_score"] = score_obj.null_score
        if score_obj.flags:
            attrs["evidence_flags"] = ",".join(score_obj.flags)
        if score_obj.evidence_reason:
            attrs["evidence_reason"] = score_obj.evidence_reason
        # Optimal transport attrs
        if score_obj.ot_attrs:
            for k, v in score_obj.ot_attrs.items():
                attrs[k] = v
        # Causal discovery attrs
        if score_obj.causal_attrs:
            for k, v in score_obj.causal_attrs.items():
                attrs[k] = v
        edge.attrs = attrs

    # Sheaf cohomology diagnostics
    if getattr(config, "sheaf_cohomology_enabled", False):
        try:
            from .cohomology import attach_cohomology_attrs

            all_maps = build_edge_maps(
                selected,
                node_estimates if has_chemistry and node_estimates else {},
                config,
                prior_weight=_get_config_float(config, "sheaf_weight_head_prior", 1.0),
            )
            attach_cohomology_attrs(
                selected,
                all_maps,
                dim=len(config.ion_order) if has_chemistry else None,
                compute_leverage=True,
            )
        except Exception:
            logger.warning("Sheaf cohomology diagnostics failed; continuing.", exc_info=True)

    if use_hydraulic_hodge and selected:
        try:
            from .hydraulic_hodge import attach_head_hodge_attrs

            attach_head_hodge_attrs(
                selected,
                sample_map,
                config,
                reference_distance_km=reference_distance_km,
                compute_leverage=True,
            )
        except Exception:
            logger.warning(
                "Hydraulic Hodge diagnostics failed; continuing.",
                exc_info=True,
            )

    # Bayesian topology posterior
    if getattr(config, "topology_posterior_enabled", False):
        try:
            from ..inference.topology_posterior import (
                attach_posterior_attrs,
                make_topology_cost_fn,
                run_topology_posterior,
            )

            cost_fn = make_topology_cost_fn(
                sample_map=sample_map,
                config=config,
                node_info=node_info,
                stats=stats,
                reference_distance_km=reference_distance_km,
            )
            posterior_result = run_topology_posterior(
                universe=candidate_list,
                cost_fn=cost_fn,
                config=config,
                initial_edges=selected,
            )
            attach_posterior_attrs(
                selected_edges=selected,
                candidate_edges=candidate_list,
                posterior_result=posterior_result,
                mode="diagnostic",
            )
        except Exception:
            logger.warning(
                "Topology posterior sampling failed; continuing.", exc_info=True
            )

    return selected
