"""Per-edge fitting pipeline."""

import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Tuple

from ..config import Config
from ..data.qc import qc_flags
from ..models.ec_tds import ec_tds_penalty
from ..models.gibbs import gibbs_evaporation_penalty, compute_gibbs_metrics
from ..models.reactions import ReactionFit, build_reaction_dictionary, fit_reactions
from ..models.transport import fit_evaporation, fit_mixing
from ..isotopes import extract_isotopes, isotope_penalty


@dataclass
class EdgeResult:
    edge_id: str
    u: str
    v: str
    transport_model: str
    gamma: Optional[float] = None
    f: Optional[float] = None
    endmember_id: Optional[str] = None
    z_extents: List[float] = field(default_factory=list)
    z_labels: List[str] = field(default_factory=list)
    transport_residual_norm: float = 0.0
    anomaly_norm: float = 0.0
    objective_score: float = 0.0
    l1_norm: float = 0.0
    reaction_iterations: int = 0
    reaction_converged: bool = True
    ec_tds_penalty: float = 0.0
    transport_probabilities: Dict[str, float] = field(default_factory=dict)
    candidate_scores: List[Dict[str, object]] = field(default_factory=list)
    constraints_active: Dict[str, str] = field(default_factory=dict)
    si_u: Dict[str, float] = field(default_factory=dict)
    si_v: Dict[str, float] = field(default_factory=dict)
    phreeqc_ok: bool = False
    charge_error: Optional[float] = None
    skipped_reason: Optional[str] = None
    edge_confidence: Optional[float] = None
    edge_distance_km: Optional[float] = None
    edge_delta_h: Optional[float] = None
    edge_sigma_delta_h: Optional[float] = None
    edge_source_tier: Optional[str] = None
    edge_flags: Optional[str] = None
    isotope_penalty: float = 0.0
    isotope_metrics: Dict[str, float] = field(default_factory=dict)
    isotope_used: bool = False
    gibbs_penalty: float = 0.0
    gibbs_metrics: Dict[str, object] = field(default_factory=dict)
    gibbs_used: bool = False
    isotope_consistency_penalty: float = 0.0
    qc_flags: List[str] = field(default_factory=list)
    nitrate_source_p_manure: Optional[float] = None
    nitrate_source_logit: Optional[float] = None
    nitrate_source_evidence: List[str] = field(default_factory=list)
    nitrate_source_gates: List[str] = field(default_factory=list)

    # Uncertainty fields (populated if uncertainty_method != "none")
    gamma_std: Optional[float] = None
    gamma_ci_low: Optional[float] = None
    gamma_ci_high: Optional[float] = None
    f_std: Optional[float] = None
    f_ci_low: Optional[float] = None
    f_ci_high: Optional[float] = None
    extents_std: List[float] = field(default_factory=list)
    extents_ci_low: List[float] = field(default_factory=list)
    extents_ci_high: List[float] = field(default_factory=list)
    uncertainty_method: Optional[str] = None
    reaction_fit: Optional[ReactionFit] = None
    residual_vector: List[float] = field(default_factory=list)


def fit_edge(
    x_u: List[float],
    x_v: List[float],
    config: Config,
    edge_id: str = "",
    u: str = "",
    v: str = "",
    obs_v: Optional[Mapping[str, float]] = None,
    bounds: Optional[Dict[str, object]] = None,
    obs_u: Optional[Mapping[str, float]] = None,
) -> EdgeResult:
    config.validate()

    reaction_matrix, labels, _ = build_reaction_dictionary(config)
    signed_mask = [label in config.signed_reaction_labels for label in labels]
    lb = bounds.get("lb") if bounds else None
    ub = bounds.get("ub") if bounds else None

    candidates: List[Tuple[str, Optional[str], Optional[float], Optional[float], List[float], float]] = []
    if "evap" in config.transport_models_enabled:
        gamma, evap_residual, evap_norm = fit_evaporation(x_u, x_v, config.weights)
        candidates.append(("evap", None, gamma, None, evap_residual, evap_norm))

    if "mix" in config.transport_models_enabled:
        for end_id, endmember in config.mixing_endmembers.items():
            f, mix_residual, mix_norm = fit_mixing(x_u, x_v, endmember, config.weights)
            candidates.append(("mix", end_id, None, f, mix_residual, mix_norm))

    best_result: Optional[EdgeResult] = None
    candidate_entries: List[Dict[str, object]] = []
    for transport_model, end_id, gamma_value, f_value, residual, transport_norm in candidates:
        reaction_fit: ReactionFit = fit_reactions(
            residual,
            reaction_matrix,
            weights=config.weights,
            lambda_l1=config.lambda_l1_value(),
            max_iter=config.reaction_max_iter,
            tol=config.reaction_tol,
            signed_mask=signed_mask,
            lb=lb,
            ub=ub,
        )

        modeled_x_v = [obs - r for obs, r in zip(x_v, reaction_fit.residual)]
        penalty = 0.0
        if obs_v is not None and config.ec_tds_penalty_enabled:
            penalty = ec_tds_penalty(modeled_x_v, obs_v, config)

        iso_penalty = 0.0
        iso_metrics: Dict[str, float] = {}
        iso_used = False
        if config.isotope_enabled and obs_u is not None and obs_v is not None and config.lmwl_defined:
            iso_u = extract_isotopes(obs_u, config.isotope_d18o_key, config.isotope_d2h_key)
            iso_v = extract_isotopes(obs_v, config.isotope_d18o_key, config.isotope_d2h_key)
            if iso_u and iso_v:
                iso_raw, iso_metrics = isotope_penalty(
                    iso_u[0],
                    iso_u[1],
                    iso_v[0],
                    iso_v[1],
                    config.lmwl_a,
                    config.lmwl_b,
                    transport_model,
                    d_excess_weight=config.isotope_d_excess_weight,
                )
                iso_penalty = config.isotope_weight * iso_raw
                iso_used = True

        # Gibbs fallback/supplement
        gibbs_penalty_val = 0.0
        gibbs_metrics_val: Dict[str, object] = {}
        gibbs_used = False
        if config.gibbs_enabled and obs_v is not None:
            gibbs_metrics_val = compute_gibbs_metrics(obs_v)
            raw_gibbs = gibbs_evaporation_penalty(
                obs_v,
                tds_precipitation=config.gibbs_tds_precipitation,
                tds_evaporation=config.gibbs_tds_evaporation,
            )
            if transport_model == "evap":
                # Penalize evaporation if Gibbs suggests otherwise
                if iso_used:
                    # Isotope available: Gibbs supplements with lower weight
                    gibbs_penalty_val = config.gibbs_weight * 0.3 * raw_gibbs
                else:
                    # No isotope: Gibbs replaces with higher weight
                    gibbs_penalty_val = config.gibbs_weight * raw_gibbs
                gibbs_used = True

        iso_consistency_penalty = 0.0
        if transport_model == "evap" and iso_used and obs_u is not None and obs_v is not None:
            # Cross-check chloride shift vs isotopic enrichment
            cl_idx = -1
            if "Cl" in config.ion_order:
                cl_idx = config.ion_order.index("Cl")
            
            if cl_idx >= 0:
                cl_u = x_u[cl_idx]
                cl_v = x_v[cl_idx]
                # If Cl increases much more than isotopes suggest, penalty
                # Or if isotopes show enrichment but Cl stays flat, penalty
                cl_ratio = (cl_v / cl_u) if cl_u > 0 else 1.0
                # Gamma is the concentration factor from biology/physics
                # Isotopes should roughly increase linearly with ln(gamma)
                # Let's use a simpler heuristic: mismatch between cl_ratio and gamma_value
                # If gamma is 1.2 but Cl ratio is 10.0, it's likely dissolution, not evap.
                if gamma_value is not None:
                    mismatch = abs(cl_ratio - gamma_value)
                    # Only penalize if the mismatch is substantial
                    if mismatch > 0.5:
                        iso_consistency_penalty = config.isotope_consistency_weight * mismatch

        objective = (
            reaction_fit.residual_norm
            + config.lambda_l1_value() * reaction_fit.l1_norm
            + penalty
            + iso_penalty
            + gibbs_penalty_val
            + iso_consistency_penalty
        )

        candidate_entries.append(
            {
                "transport_model": transport_model,
                "endmember_id": end_id,
                "objective_score": objective,
                "transport_residual_norm": transport_norm,
            }
        )

        result = EdgeResult(
            edge_id=edge_id,
            u=u,
            v=v,
            transport_model=transport_model,
            gamma=gamma_value,
            f=f_value,
            endmember_id=end_id,
            z_extents=reaction_fit.extents,
            z_labels=labels,
            transport_residual_norm=transport_norm,
            anomaly_norm=reaction_fit.residual_norm,
            objective_score=objective,
            l1_norm=reaction_fit.l1_norm,
            reaction_iterations=reaction_fit.iterations,
            reaction_converged=reaction_fit.converged,
            ec_tds_penalty=penalty,
            qc_flags=[],
            constraints_active=bounds.get("constraints_active", {}) if bounds else {},
            si_u=bounds.get("si_u", {}) if bounds else {},
            si_v=bounds.get("si_v", {}) if bounds else {},
            phreeqc_ok=bounds.get("phreeqc_ok", False) if bounds else False,
            charge_error=bounds.get("charge_error") if bounds else None,
            skipped_reason=bounds.get("skipped_reason") if bounds else None,
            isotope_penalty=iso_penalty,
            isotope_metrics=iso_metrics,
            isotope_used=iso_used,
            gibbs_penalty=gibbs_penalty_val,
            gibbs_metrics=gibbs_metrics_val,
            gibbs_used=gibbs_used,
            isotope_consistency_penalty=iso_consistency_penalty,
        )
        if best_result is None or result.objective_score < best_result.objective_score:
            best_result = result

    if best_result is None:
        raise RuntimeError("No transport candidates evaluated.")

    scores = [entry["objective_score"] for entry in candidate_entries]
    min_score = min(scores)
    weights = [math.exp(-(score - min_score)) for score in scores]
    total = sum(weights) or 1.0
    transport_probs: Dict[str, float] = {}
    for entry, weight in zip(candidate_entries, weights):
        transport = str(entry["transport_model"])
        transport_probs[transport] = transport_probs.get(transport, 0.0) + weight / total

    qc = qc_flags(x_v, config.ion_order, config.charge_balance_limit)
    if config.ec_tds_penalty_limit and best_result.ec_tds_penalty > config.ec_tds_penalty_limit:
        qc.append("ec_tds_consistency")
    best_result.qc_flags = qc
    best_result.transport_probabilities = transport_probs
    best_result.candidate_scores = candidate_entries
    return best_result
