"""
Residence time estimation methods.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from . import TemporalNode

try:
    from hydrosheaf.nuclear.invert import infer_age_from_tracer
    from hydrosheaf.nuclear.nuclides import get_nuclide
    _NUCLEAR_AVAILABLE = True
except ImportError:
    _NUCLEAR_AVAILABLE = False



def estimate_residence_time(
    node_u: TemporalNode,
    node_v: TemporalNode,
    method: str = "cross_correlation",
    tracer_ion: str = "Cl",
    ion_order: Optional[List[str]] = None,
    hydraulic_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float, str]:
    """
    Estimate groundwater travel time between two nodes.

    Parameters
    ----------
    node_u, node_v : TemporalNode
        Upstream and downstream nodes with time-series
    method : str
        "cross_correlation" - lag via signal correlation
        "gradient" - Darcy's law estimate
        "tracer_decay" - radioactive tracer
    tracer_ion : str
        Conservative tracer for correlation (typically "Cl")
    ion_order : Optional[List[str]]
        Ion order to find tracer index
    hydraulic_params : Optional[Dict]
        Required for "gradient" method:
        {"distance_m": float, "K_m_day": float, "gradient": float, "porosity": float}

    Returns
    -------
    Tuple[float, float, str]
        (estimated_tau_days, uncertainty_days, method_used)

    Mathematical Implementation
    ---------------------------
    Cross-correlation:
        1. Extract tracer time-series: u_t = C_u^{tracer}(t), v_t = C_v^{tracer}(t)
        2. Compute normalized cross-correlation:
           r(τ) = Σ_t (u_t - μ_u)(v_{t+τ} - μ_v) / (σ_u σ_v N)
        3. Find τ* = argmax_τ r(τ)
        4. Uncertainty from peak width at r(τ*) - 0.1

    Gradient:
        τ = (distance_m * porosity) / (K_m_day * gradient)

    Tracer decay (e.g., tritium, t_1/2 = 12.32 years):
        τ = -t_1/2 / ln(2) * ln(C_v / C_u)
    """
    tau, uncertainty, used, _, _ = estimate_residence_time_with_details(
        node_u,
        node_v,
        method=method,
        tracer_ion=tracer_ion,
        ion_order=ion_order,
        hydraulic_params=hydraulic_params,
    )
    return tau, uncertainty, used


def estimate_residence_time_with_details(
    node_u: TemporalNode,
    node_v: TemporalNode,
    method: str = "cross_correlation",
    tracer_ion: str = "Cl",
    ion_order: Optional[List[str]] = None,
    hydraulic_params: Optional[Dict[str, float]] = None,
) -> Tuple[float, float, str, Dict[str, object], List[str]]:
    """Estimate residence time and return diagnostics for multi-tracer consensus.

    Returns
    -------
    Tuple[float, float, str, Dict[str, object], List[str]]
        (tau_days, uncertainty_days, method_used, details, flags)
    """
    if method == "cross_correlation":
        return _estimate_residence_time_cross_correlation_consensus(
            node_u,
            node_v,
            tracer_ion=tracer_ion,
            ion_order=ion_order,
            hydraulic_params=hydraulic_params,
        )
    elif method == "bayesian_lag":
        return _estimate_residence_time_bayesian_lag_consensus(
            node_u,
            node_v,
            tracer_ion=tracer_ion,
            ion_order=ion_order,
            hydraulic_params=hydraulic_params,
        )
    elif method == "ttd":
        return _estimate_residence_time_ttd_convolution_consensus(
            node_u,
            node_v,
            tracer_ion=tracer_ion,
            ion_order=ion_order,
            hydraulic_params=hydraulic_params,
        )
    elif method == "gradient":
        tau, uncertainty, used = _estimate_residence_time_gradient(hydraulic_params)
        return (
            tau,
            uncertainty,
            used,
            {"physics": {"tau_days": tau, "uncertainty_days": uncertainty}},
            [],
        )
    elif method == "tracer_decay":
        tau, uncertainty, used = _estimate_residence_time_tracer_decay(
            node_u, node_v, tracer_ion, ion_order
        )
        return tau, uncertainty, used, {}, []
    elif method == "recharge_piston":
        # Placeholder average lag for reporting
        # Actual variable lag logic is handled in align_time_series_recharge
        tau = 0.0
        details = {}
        if hydraulic_params:
            storage = hydraulic_params.get("storage_mm", 500.0)
            avg_recharge = hydraulic_params.get("avg_recharge_mm_day", 2.0)
            if avg_recharge > 0:
                tau = storage / avg_recharge
            # Only store safe summary types
            details = {
                "method": "recharge_piston",
                "storage_mm": float(storage),
                "avg_recharge_mm_day": float(avg_recharge),
                "estimated_avg_lag_days": float(tau),
            }
        return tau, tau * 0.5, "recharge_piston", details, []
    else:
        raise ValueError(f"Unknown method: {method}")


def _parse_tracer_candidates(tracer_ion: str) -> List[str]:
    raw = (tracer_ion or "").strip()
    if not raw:
        return []
    if raw.lower() == "auto":
        return ["Cl", "18O", "2H"]
    if "," in raw:
        return [item.strip() for item in raw.split(",") if item.strip()]
    return [raw]


def _extract_tracer_series(
    node: TemporalNode,
    tracer: str,
    ion_order: Optional[List[str]],
) -> Optional[np.ndarray]:
    if not node.samples:
        return None

    if ion_order is not None and tracer in ion_order:
        try:
            tracer_idx = ion_order.index(tracer)
        except ValueError:
            tracer_idx = None
        if tracer_idx is not None:
            return np.array(
                [float(s.concentrations[tracer_idx]) for s in node.samples], dtype=float
            )

    synonyms: List[str] = []
    if tracer in {"18O", "d18O"}:
        synonyms = ["18O", "d18O"]
    elif tracer in {"2H", "d2H"}:
        synonyms = ["2H", "d2H"]
    elif tracer in {"3H", "tritium", "H3"}:
        synonyms = ["3H", "tritium", "H3", "H-3"]
    elif tracer in {"14C", "C14", "carbon14"}:
        synonyms = ["14C", "C14", "carbon14", "C-14"]
    else:
        synonyms = [tracer]


    values: List[float] = []
    for sample in node.samples:
        if not sample.isotopes:
            return None
        if tracer in sample.isotopes:
            values.append(float(sample.isotopes[tracer]))
            continue
        found = False
        for key in synonyms:
            if key in sample.isotopes:
                values.append(float(sample.isotopes[key]))
                found = True
                break
        if not found:
            return None

    return np.array(values, dtype=float)


def _estimate_residence_time_cross_correlation(
    node_u: TemporalNode,
    node_v: TemporalNode,
    tracer_ion: str,
    ion_order: Optional[List[str]],
) -> Tuple[float, float, str]:
    """
    Estimate residence time via cross-correlation of tracer signal.
    """
    if not node_u.samples or not node_v.samples:
        return 0.0, 0.0, "cross_correlation_failed"

    tau, uncertainty, used, _ = _estimate_residence_time_cross_correlation_detailed(
        node_u, node_v, tracer_ion=tracer_ion, ion_order=ion_order
    )
    return tau, uncertainty, used


def _estimate_residence_time_cross_correlation_detailed(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    tracer_ion: str,
    ion_order: Optional[List[str]],
) -> Tuple[float, float, str, Dict[str, float]]:
    if not node_u.samples or not node_v.samples:
        return 0.0, 0.0, "cross_correlation_failed", {}

    # Use absolute time axis (days since epoch) to avoid origin mismatch between nodes.
    u_times = np.array(
        [s.timestamp.timestamp() / 86400.0 for s in node_u.samples], dtype=float
    )
    v_times = np.array(
        [s.timestamp.timestamp() / 86400.0 for s in node_v.samples], dtype=float
    )

    u_values = _extract_tracer_series(node_u, tracer_ion, ion_order)
    v_values = _extract_tracer_series(node_v, tracer_ion, ion_order)
    if u_values is None or v_values is None:
        return 0.0, 0.0, f"cross_correlation_missing_tracer({tracer_ion})", {}

    if len(u_values) < 3 or len(v_values) < 3:
        return 0.0, 0.0, f"cross_correlation_insufficient_data({tracer_ion})", {}

    u_mean = float(np.mean(u_values))
    v_mean = float(np.mean(v_values))
    u_std = float(np.std(u_values))
    v_std = float(np.std(v_values))
    if u_std < 1e-6 or v_std < 1e-6:
        return (
            0.0,
            0.0,
            f"cross_correlation_no_variation({tracer_ion})",
            {"u_std": u_std, "v_std": v_std},
        )

    u_norm = (u_values - u_mean) / u_std
    v_norm = (v_values - v_mean) / v_std

    max_lag_days = min(365.0, float(v_times[-1] - u_times[0]))
    if max_lag_days <= 0:
        return 0.0, 0.0, f"cross_correlation_no_overlap({tracer_ion})", {}
    lag_step_days = max(1.0, max_lag_days / 100.0)

    lags = np.arange(0.0, max_lag_days, lag_step_days, dtype=float)
    correlations: List[float] = []
    for lag in lags:
        corr_sum = 0.0
        count = 0
        for i, v_t in enumerate(v_times):
            u_t_target = v_t - lag
            if u_t_target < u_times[0] or u_t_target > u_times[-1]:
                continue
            u_interp = float(np.interp(u_t_target, u_times, u_norm))
            corr_sum += u_interp * float(v_norm[i])
            count += 1
        correlations.append(corr_sum / count if count > 0 else 0.0)

    corr_arr = np.array(correlations, dtype=float)
    max_corr = float(np.max(corr_arr)) if corr_arr.size else 0.0
    med_corr = float(np.median(corr_arr)) if corr_arr.size else 0.0
    peak_idx = int(np.argmax(corr_arr)) if corr_arr.size else 0
    peak_lag = float(lags[peak_idx]) if lags.size else 0.0

    positive_corr_mask = (
        corr_arr > 0.1 * max_corr
        if max_corr > 0
        else np.zeros_like(corr_arr, dtype=bool)
    )
    if int(np.sum(positive_corr_mask)) < 2:
        best_lag = peak_lag
        uncertainty = 10.0
    else:
        valid_lags = lags[positive_corr_mask]
        valid_corrs = corr_arr[positive_corr_mask]
        sum_weights = float(np.sum(valid_corrs)) or 1.0
        center_of_mass = float(np.sum(valid_lags * valid_corrs) / sum_weights)
        best_lag = center_of_mass
        variance = float(
            np.sum(((valid_lags - center_of_mass) ** 2) * valid_corrs) / sum_weights
        )
        uncertainty = float(np.sqrt(max(0.0, variance)))

    metrics = {
        "u_std": u_std,
        "v_std": v_std,
        "max_corr": max_corr,
        "median_corr": med_corr,
        "peak_lag": peak_lag,
        "max_lag_days": float(max_lag_days),
    }
    return (
        float(best_lag),
        float(uncertainty),
        f"cross_correlation({tracer_ion})",
        metrics,
    )


def _isotope_evaporation_gate(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    lmwl_intercept: float,
    lmwl_slope: float,
) -> Tuple[float, List[str], Dict[str, float]]:
    d18o_u = _extract_tracer_series(node_u, "18O", None)
    d2h_u = _extract_tracer_series(node_u, "2H", None)
    d18o_v = _extract_tracer_series(node_v, "18O", None)
    d2h_v = _extract_tracer_series(node_v, "2H", None)
    if d18o_u is None or d2h_u is None or d18o_v is None or d2h_v is None:
        return 1.0, ["isotope_missing"], {}

    def _rmse_residual(d18o: np.ndarray, d2h: np.ndarray) -> float:
        resid = d2h - (lmwl_intercept + lmwl_slope * d18o)
        return float(np.sqrt(np.mean(resid**2)))

    def _d_excess_median(d18o: np.ndarray, d2h: np.ndarray) -> float:
        d_excess = d2h - 8.0 * d18o
        return float(np.median(d_excess))

    rmse_u = _rmse_residual(d18o_u, d2h_u)
    rmse_v = _rmse_residual(d18o_v, d2h_v)
    d_ex_u = _d_excess_median(d18o_u, d2h_u)
    d_ex_v = _d_excess_median(d18o_v, d2h_v)

    flags: List[str] = []
    weight = 1.0
    # Heuristic: strong departure from LMWL and/or collapsed d-excess indicates local evaporation fractionation.
    if rmse_u > 6.0 or rmse_v > 6.0:
        flags.append("isotope_lmwl_departure")
        weight *= 0.35
    if d_ex_u < 5.0 or d_ex_v < 5.0:
        flags.append("isotope_low_d_excess")
        weight *= 0.35

    metrics = {
        "lmwl_rmse_u": rmse_u,
        "lmwl_rmse_v": rmse_v,
        "d_excess_u_med": d_ex_u,
        "d_excess_v_med": d_ex_v,
    }
    return weight, flags, metrics


def _chloride_nonconservative_gate(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    ion_order: Optional[List[str]],
) -> Tuple[float, List[str], Dict[str, float]]:
    cl_u = _extract_tracer_series(node_u, "Cl", ion_order)
    cl_v = _extract_tracer_series(node_v, "Cl", ion_order)
    if cl_u is None or cl_v is None:
        return 1.0, ["cl_missing"], {}

    def _stepiness(values: np.ndarray) -> float:
        if values.size < 4:
            return 0.0
        diffs = np.diff(values)
        denom = float(np.std(values)) + 1e-12
        return float(np.percentile(np.abs(diffs), 95) / denom)

    step_u = _stepiness(cl_u)
    step_v = _stepiness(cl_v)
    step_ratio = step_v / (step_u + 1e-6)

    flags: List[str] = []
    weight = 1.0

    if step_v > 3.0 and step_ratio > 2.0:
        flags.append("cl_step_change_downstream")
        weight *= 0.4

    # Optional Na/Cl ratio mismatch check (helps flag NaCl-like point sources or strong ionic perturbations).
    if ion_order is not None and "Na" in ion_order and "Cl" in ion_order:
        na_u = _extract_tracer_series(node_u, "Na", ion_order)
        na_v = _extract_tracer_series(node_v, "Na", ion_order)
        if na_u is not None and na_v is not None:
            ratio_u = float(np.median(na_u / (cl_u + 1e-9)))
            ratio_v = float(np.median(na_v / (cl_v + 1e-9)))
            ratio_change = max(ratio_u, ratio_v) / max(1e-9, min(ratio_u, ratio_v))
            if ratio_change > 2.5 and "cl_step_change_downstream" in flags:
                flags.append("na_cl_ratio_shift")
                weight *= 0.5
            metrics = {
                "cl_step_u": step_u,
                "cl_step_v": step_v,
                "na_cl_ratio_u_med": ratio_u,
                "na_cl_ratio_v_med": ratio_v,
            }
            return weight, flags, metrics

    metrics = {"cl_step_u": step_u, "cl_step_v": step_v}
    return weight, flags, metrics


def _estimate_residence_time_cross_correlation_consensus(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    tracer_ion: str,
    ion_order: Optional[List[str]],
    hydraulic_params: Optional[Dict[str, float]],
) -> Tuple[float, float, str, Dict[str, object], List[str]]:
    tracers = _parse_tracer_candidates(tracer_ion)
    if not tracers:
        tracers = [tracer_ion]

    # Consensus settings (can be overridden via hydraulic_params, if provided).
    tol = float((hydraulic_params or {}).get("agreement_tolerance", 0.4))
    min_peak_corr = float((hydraulic_params or {}).get("min_peak_corr", 0.2))
    max_rel_unc = float((hydraulic_params or {}).get("max_relative_uncertainty", 1.5))
    max_unc_days = float((hydraulic_params or {}).get("max_uncertainty_days", 180.0))

    lmwl_intercept = float((hydraulic_params or {}).get("lmwl_intercept", 10.0))
    lmwl_slope = float((hydraulic_params or {}).get("lmwl_slope", 8.0))

    details: Dict[str, object] = {"candidates": {}}
    flags: List[str] = []

    # Gates
    isotope_weight, isotope_flags, isotope_metrics = _isotope_evaporation_gate(
        node_u, node_v, lmwl_intercept=lmwl_intercept, lmwl_slope=lmwl_slope
    )
    if isotope_flags and isotope_flags != ["isotope_missing"]:
        flags.extend(isotope_flags)
    details["isotope_gate"] = {
        "weight": isotope_weight,
        "metrics": isotope_metrics,
        "flags": isotope_flags,
    }

    cl_weight, cl_flags, cl_metrics = _chloride_nonconservative_gate(
        node_u, node_v, ion_order=ion_order
    )
    if cl_flags and cl_flags != ["cl_missing"]:
        flags.extend(cl_flags)
    details["chloride_gate"] = {
        "weight": cl_weight,
        "metrics": cl_metrics,
        "flags": cl_flags,
    }

    accepted: List[Tuple[str, float, float, float]] = []  # (tracer, tau, unc, weight)
    for tracer in tracers:
        tau, unc, used, metrics = _estimate_residence_time_cross_correlation_detailed(
            node_u, node_v, tracer_ion=tracer, ion_order=ion_order
        )
        candidate_flags: List[str] = []
        ok = used.startswith("cross_correlation(")
        if not ok:
            candidate_flags.append(used)
        max_corr = float(metrics.get("max_corr", 0.0))
        max_lag_days = float(metrics.get("max_lag_days", 0.0))
        rel_unc = (
            float(unc) / (float(tau) + 1e-9)
            if float(tau) > 0
            else float("inf") if float(unc) > 0 else 0.0
        )

        if ok and max_corr < min_peak_corr:
            ok = False
            candidate_flags.append("low_peak_corr")
        if ok and (
            unc > max_unc_days
            or (max_lag_days > 0 and unc > max_lag_days)
            or rel_unc > max_rel_unc
        ):
            ok = False
            candidate_flags.append("high_uncertainty")

        gate_weight = 1.0
        if tracer in {"18O", "2H", "d18O", "d2H"}:
            gate_weight *= isotope_weight
        if tracer == "Cl":
            gate_weight *= cl_weight

        base_weight = (max(0.0, max_corr) + 1e-6) / (float(unc) + 1.0)
        weight = gate_weight * base_weight

        details["candidates"][tracer] = {
            "tau_days": tau,
            "uncertainty_days": unc,
            "max_corr": max_corr,
            "peak_lag": float(metrics.get("peak_lag", 0.0)),
            "accepted": bool(ok),
            "weight": float(weight),
            "flags": candidate_flags,
        }
        if ok and weight > 0:
            accepted.append((tracer, float(tau), float(unc), float(weight)))

    if not accepted:
        # Physics fallback if available
        tau_phy, unc_phy, used_phy = _estimate_residence_time_gradient(hydraulic_params)
        if used_phy == "gradient" and tau_phy > 0:
            flags.append("tau_fallback_physics")
            details["physics"] = {"tau_days": tau_phy, "uncertainty_days": unc_phy}
            return tau_phy, unc_phy, "gradient", details, flags
        flags.append("tau_failed_all_tracers")
        return 0.0, 0.0, "cross_correlation_failed", details, flags

    # Disagreement handling across accepted tracers
    taus = [item[1] for item in accepted]
    tau_min = min(taus)
    tau_max = max(taus)
    rel_spread = (tau_max - tau_min) / max(1e-9, tau_max)
    details["consensus"] = {"rel_spread": float(rel_spread), "tolerance": float(tol)}

    accepted.sort(key=lambda item: item[3], reverse=True)
    best_tracer, best_tau, best_unc, best_w = accepted[0]

    if rel_spread > tol and len(accepted) >= 2:
        flags.append("tau_ambiguous")
        # Prefer the best-weight tracer but widen uncertainty to reflect disagreement.
        spread_days = float(tau_max - tau_min)
        inflated_unc = float(np.sqrt(best_unc**2 + (0.5 * spread_days) ** 2))

        # Physics prior blend if available and strongly inconsistent with best tracer.
        tau_phy, unc_phy, used_phy = _estimate_residence_time_gradient(hydraulic_params)
        if used_phy == "gradient" and tau_phy > 0:
            rel_diff_phy = abs(best_tau - tau_phy) / max(best_tau, tau_phy, 1e-9)
            details["physics"] = {
                "tau_days": tau_phy,
                "uncertainty_days": unc_phy,
                "rel_diff_best": float(rel_diff_phy),
            }
            if rel_diff_phy > 0.5:
                flags.append("physics_prior_blend")
                p_best = 1.0 / max(1e-9, best_unc) ** 2
                p_phy = 1.0 / max(1e-9, float(unc_phy)) ** 2
                blended_tau = (best_tau * p_best + float(tau_phy) * p_phy) / (
                    p_best + p_phy
                )
                blended_unc = float(
                    np.sqrt(1.0 / (p_best + p_phy) + (0.5 * spread_days) ** 2)
                )
                return (
                    float(blended_tau),
                    float(blended_unc),
                    f"cross_correlation_consensus({best_tracer}+physics)",
                    details,
                    flags,
                )

        return (
            best_tau,
            inflated_unc,
            f"cross_correlation_consensus({best_tracer})",
            details,
            flags,
        )

    # Normal case: weighted mean across accepted tracers
    total_w = sum(item[3] for item in accepted) or 1.0
    consensus_tau = sum(item[1] * item[3] for item in accepted) / total_w
    # Uncertainty: weighted RMS of per-tracer uncertainties + inter-tracer spread.
    mean_unc = sum(item[2] * item[3] for item in accepted) / total_w
    spread = (
        float(np.std(np.array(taus, dtype=float), ddof=0)) if len(taus) > 1 else 0.0
    )
    consensus_unc = float(np.sqrt(mean_unc**2 + spread**2))

    details["consensus"].update(
        {
            "tau_days": float(consensus_tau),
            "uncertainty_days": float(consensus_unc),
            "accepted_tracers": [item[0] for item in accepted],
            "weights": {item[0]: float(item[3]) for item in accepted},
        }
    )
    return (
        float(consensus_tau),
        float(consensus_unc),
        "cross_correlation_consensus",
        details,
        flags,
    )


def _build_uniform_grid_overlap(
    node_u: TemporalNode, node_v: TemporalNode, *, dt_days: float
) -> Optional[np.ndarray]:
    if not node_u.samples or not node_v.samples or dt_days <= 0:
        return None
    u_start = node_u.samples[0].timestamp.timestamp() / 86400.0
    u_end = node_u.samples[-1].timestamp.timestamp() / 86400.0
    v_start = node_v.samples[0].timestamp.timestamp() / 86400.0
    v_end = node_v.samples[-1].timestamp.timestamp() / 86400.0
    start = max(u_start, v_start)
    end = min(u_end, v_end)
    if end <= start:
        return None
    n = int(np.floor((end - start) / dt_days)) + 1
    if n < 5:
        return None
    return start + np.arange(n, dtype=float) * float(dt_days)


def _resample_on_grid(
    node: TemporalNode,
    *,
    tracer: str,
    ion_order: Optional[List[str]],
    grid_days: np.ndarray,
) -> Optional[np.ndarray]:
    values = _extract_tracer_series(node, tracer, ion_order)
    if values is None:
        return None
    times = np.array(
        [s.timestamp.timestamp() / 86400.0 for s in node.samples], dtype=float
    )
    if times.size != values.size or times.size < 2:
        return None
    return np.interp(grid_days, times, values).astype(float)


def _toeplitz_lag_matrix(u: np.ndarray, n_lags: int) -> np.ndarray:
    n = int(u.size)
    m = int(n_lags)
    x = np.zeros((n, m), dtype=float)
    for lag in range(m):
        if lag == 0:
            x[:, lag] = u
        else:
            x[lag:, lag] = u[: n - lag]
    return x


def _fit_ttd_nnls(
    u: np.ndarray,
    v: np.ndarray,
    *,
    dt_days: float,
    max_lag_days: float,
    smoothness_lambda: float,
    attenuation_k_grid: np.ndarray,
) -> Tuple[Dict[str, object], List[str]]:
    flags: List[str] = []
    if u.size != v.size or u.size < 10:
        return {"ok": False, "reason": "ttd_insufficient_data"}, [
            "ttd_insufficient_data"
        ]

    v_std = float(np.std(v))
    if v_std < 1e-8:
        return {"ok": False, "reason": "ttd_no_variation"}, ["ttd_no_variation"]

    n_lags = int(np.floor(float(max_lag_days) / float(dt_days))) + 1
    if n_lags < 2:
        return {"ok": False, "reason": "ttd_invalid_lag_grid"}, ["ttd_invalid_lag_grid"]

    lags = np.arange(n_lags, dtype=float) * float(dt_days)
    x_base = _toeplitz_lag_matrix(u, n_lags)

    # Centering to handle intercept b analytically: b = mean(v - Xc)
    ybar = float(np.mean(v))
    y0 = (v - ybar).astype(float)
    xmean = np.mean(x_base, axis=0)
    x0_base = (x_base - xmean).astype(float)

    # Optional smoothness penalty (first differences)
    dmat = None
    lam = float(smoothness_lambda)
    if lam > 0 and n_lags >= 3:
        dmat = np.zeros((n_lags - 1, n_lags), dtype=float)
        for i in range(n_lags - 1):
            dmat[i, i] = -1.0
            dmat[i, i + 1] = 1.0

    best: Optional[Dict[str, object]] = None
    try:
        from scipy.optimize import nnls  # type: ignore
    except Exception:
        nnls = None  # type: ignore

    for k in attenuation_k_grid:
        atten = np.exp(-float(k) * lags)
        x0 = x0_base * atten[None, :]

        a = x0
        b = y0
        if dmat is not None and lam > 0:
            a = np.vstack([x0, np.sqrt(lam) * dmat])
            b = np.concatenate([y0, np.zeros((dmat.shape[0],), dtype=float)])

        if nnls is None:
            # Fallback: projected gradient descent (simple and robust for small lag grids).
            c = np.zeros((n_lags,), dtype=float)
            lr = 1e-3
            for _ in range(2000):
                grad = a.T @ (a @ c - b)
                c = np.maximum(0.0, c - lr * grad)
        else:
            c, _ = nnls(a, b)

        c = np.asarray(c, dtype=float)
        if float(np.sum(c)) <= 0:
            continue

        # Recover intercept for original X (not centered).
        b0 = ybar - float(np.dot(xmean * atten, c))
        yhat = x_base @ (c * atten) + b0

        sse = float(np.sum((v - yhat) ** 2))
        sst = float(np.sum((v - ybar) ** 2)) or 1e-12
        r2 = 1.0 - sse / sst

        if best is None or sse < float(best["sse"]):
            best = {
                "ok": True,
                "k": float(k),
                "c": (c * atten).tolist(),
                "intercept": float(b0),
                "sse": float(sse),
                "r2": float(r2),
                "n_lags": int(n_lags),
                "dt_days": float(dt_days),
                "max_lag_days": float(max_lag_days),
            }

    if best is None:
        return {"ok": False, "reason": "ttd_fit_failed"}, ["ttd_fit_failed"]

    c_eff = np.array(best["c"], dtype=float)
    amp = float(np.sum(c_eff))
    w = c_eff / max(1e-12, amp)
    mean_tau = float(np.sum(lags * w))
    var_tau = float(np.sum(((lags - mean_tau) ** 2) * w))
    std_tau = float(np.sqrt(max(0.0, var_tau)))
    cdf = np.cumsum(w)
    q025 = float(np.interp(0.025, cdf, lags))
    q975 = float(np.interp(0.975, cdf, lags))

    best.update(
        {
            "amp": amp,
            "tau_mean_days": mean_tau,
            "tau_std_days": std_tau,
            "tau_ci_low_days": q025,
            "tau_ci_high_days": q975,
            "weights": w.tolist(),
        }
    )
    return best, flags


def _estimate_residence_time_ttd_convolution_consensus(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    tracer_ion: str,
    ion_order: Optional[List[str]],
    hydraulic_params: Optional[Dict[str, float]],
) -> Tuple[float, float, str, Dict[str, object], List[str]]:
    tracers = _parse_tracer_candidates(tracer_ion) or [tracer_ion]

    dt_days = float((hydraulic_params or {}).get("grid_dt_days", 30.0))
    max_lag_days = float((hydraulic_params or {}).get("max_lag_days", 365.0))
    smoothness_lambda = float((hydraulic_params or {}).get("smoothness_lambda", 0.0))
    min_r2 = float((hydraulic_params or {}).get("ttd_min_r2", 0.2))
    tol = float((hydraulic_params or {}).get("agreement_tolerance", 0.4))
    k_max = float((hydraulic_params or {}).get("attenuation_k_max", 0.02))
    k_steps = int((hydraulic_params or {}).get("attenuation_k_steps", 6))
    k_steps = max(1, k_steps)
    attenuation_k_grid = np.linspace(0.0, max(0.0, k_max), k_steps, dtype=float)

    details: Dict[str, object] = {"candidates": {}}
    flags: List[str] = []

    # Gates reuse
    lmwl_intercept = float((hydraulic_params or {}).get("lmwl_intercept", 10.0))
    lmwl_slope = float((hydraulic_params or {}).get("lmwl_slope", 8.0))
    isotope_weight, isotope_flags, isotope_metrics = _isotope_evaporation_gate(
        node_u, node_v, lmwl_intercept=lmwl_intercept, lmwl_slope=lmwl_slope
    )
    if isotope_flags and isotope_flags != ["isotope_missing"]:
        flags.extend(isotope_flags)
    details["isotope_gate"] = {
        "weight": isotope_weight,
        "metrics": isotope_metrics,
        "flags": isotope_flags,
    }

    cl_weight, cl_flags, cl_metrics = _chloride_nonconservative_gate(
        node_u, node_v, ion_order=ion_order
    )
    if cl_flags and cl_flags != ["cl_missing"]:
        flags.extend(cl_flags)
    details["chloride_gate"] = {
        "weight": cl_weight,
        "metrics": cl_metrics,
        "flags": cl_flags,
    }

    accepted: List[Tuple[str, float, float, float]] = []  # (tracer, tau, unc, weight)
    grid = _build_uniform_grid_overlap(node_u, node_v, dt_days=dt_days)
    if grid is None:
        flags.append("ttd_no_overlap_grid")
        return 0.0, 0.0, "ttd_failed", details, flags

    for tracer in tracers:
        u_res = _resample_on_grid(
            node_u, tracer=tracer, ion_order=ion_order, grid_days=grid
        )
        v_res = _resample_on_grid(
            node_v, tracer=tracer, ion_order=ion_order, grid_days=grid
        )
        if u_res is None or v_res is None:
            details["candidates"][tracer] = {
                "accepted": False,
                "flags": ["missing_tracer"],
            }
            continue

        fit, fit_flags = _fit_ttd_nnls(
            u_res,
            v_res,
            dt_days=dt_days,
            max_lag_days=max_lag_days,
            smoothness_lambda=smoothness_lambda,
            attenuation_k_grid=attenuation_k_grid,
        )
        candidate_flags = list(fit_flags)
        ok = bool(fit.get("ok"))
        r2 = float(fit.get("r2", 0.0)) if ok else 0.0
        tau = float(fit.get("tau_mean_days", 0.0)) if ok else 0.0
        unc = float(fit.get("tau_std_days", 0.0)) if ok else 0.0

        if ok and r2 < min_r2:
            ok = False
            candidate_flags.append("ttd_low_r2")

        gate_weight = 1.0
        if tracer in {"18O", "2H", "d18O", "d2H"}:
            gate_weight *= isotope_weight
        if tracer == "Cl":
            gate_weight *= cl_weight

        weight = gate_weight * max(0.0, r2) / (unc + 1.0)
        details["candidates"][tracer] = {
            "accepted": bool(ok),
            "weight": float(weight),
            "r2": float(r2),
            "tau_mean_days": tau,
            "tau_std_days": unc,
            "tau_ci_low_days": fit.get("tau_ci_low_days"),
            "tau_ci_high_days": fit.get("tau_ci_high_days"),
            "k": fit.get("k"),
            "amp": fit.get("amp"),
            "flags": candidate_flags,
        }
        if ok and weight > 0 and tau > 0:
            accepted.append((tracer, tau, unc, float(weight)))

    if not accepted:
        tau_phy, unc_phy, used_phy = _estimate_residence_time_gradient(hydraulic_params)
        if used_phy == "gradient" and tau_phy > 0:
            flags.append("tau_fallback_physics")
            details["physics"] = {"tau_days": tau_phy, "uncertainty_days": unc_phy}
            return tau_phy, unc_phy, "gradient", details, flags
        flags.append("ttd_failed_all_tracers")
        return 0.0, 0.0, "ttd_failed", details, flags

    taus = [item[1] for item in accepted]
    tau_min = min(taus)
    tau_max = max(taus)
    rel_spread = (tau_max - tau_min) / max(1e-9, tau_max)
    details["consensus"] = {"rel_spread": float(rel_spread), "tolerance": float(tol)}

    accepted.sort(key=lambda item: item[3], reverse=True)
    best_tracer, best_tau, best_unc, _ = accepted[0]

    if rel_spread > tol and len(accepted) >= 2:
        flags.append("tau_ambiguous")
        spread_days = float(tau_max - tau_min)
        inflated_unc = float(np.sqrt(best_unc**2 + (0.5 * spread_days) ** 2))

        tau_phy, unc_phy, used_phy = _estimate_residence_time_gradient(hydraulic_params)
        if used_phy == "gradient" and tau_phy > 0:
            rel_diff_phy = abs(best_tau - tau_phy) / max(best_tau, tau_phy, 1e-9)
            details["physics"] = {
                "tau_days": tau_phy,
                "uncertainty_days": unc_phy,
                "rel_diff_best": float(rel_diff_phy),
            }
            blend_threshold = float(
                (hydraulic_params or {}).get("physics_blend_threshold", 0.5)
            )
            if rel_diff_phy > blend_threshold:
                flags.append("physics_prior_blend")
                p_best = 1.0 / max(1e-9, best_unc) ** 2
                p_phy = 1.0 / max(1e-9, float(unc_phy)) ** 2
                blended_tau = (best_tau * p_best + float(tau_phy) * p_phy) / (
                    p_best + p_phy
                )
                blended_unc = float(
                    np.sqrt(1.0 / (p_best + p_phy) + (0.5 * spread_days) ** 2)
                )
                return (
                    float(blended_tau),
                    float(blended_unc),
                    f"ttd_consensus({best_tracer}+physics)",
                    details,
                    flags,
                )

        return best_tau, inflated_unc, f"ttd_consensus({best_tracer})", details, flags

    total_w = sum(item[3] for item in accepted) or 1.0
    consensus_tau = sum(item[1] * item[3] for item in accepted) / total_w
    mean_unc = sum(item[2] * item[3] for item in accepted) / total_w
    spread = (
        float(np.std(np.array(taus, dtype=float), ddof=0)) if len(taus) > 1 else 0.0
    )
    consensus_unc = float(np.sqrt(mean_unc**2 + spread**2))

    details["consensus"].update(
        {
            "tau_days": float(consensus_tau),
            "uncertainty_days": float(consensus_unc),
            "accepted_tracers": [item[0] for item in accepted],
            "weights": {item[0]: float(item[3]) for item in accepted},
        }
    )
    return float(consensus_tau), float(consensus_unc), "ttd_consensus", details, flags


def _bayes_truncated_normal_prior(
    log_tau: np.ndarray, mu: float, sigma: float
) -> np.ndarray:
    sigma = max(1e-9, float(sigma))
    return -0.5 * ((log_tau - float(mu)) / sigma) ** 2


def _ols_sse_with_intercept(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    """Return (sse, a, b) for y ≈ a*x + b."""
    if x.size != y.size or x.size < 3:
        return float("inf"), 0.0, float(np.mean(y) if y.size else 0.0)
    x_mean = float(np.mean(x))
    y_mean = float(np.mean(y))
    x_var = float(np.sum((x - x_mean) ** 2))
    if x_var < 1e-12:
        b0 = y_mean
        sse = float(np.sum((y - b0) ** 2))
        return sse, 0.0, b0
    cov = float(np.sum((x - x_mean) * (y - y_mean)))
    a = cov / x_var
    b = y_mean - a * x_mean
    resid = y - (a * x + b)
    sse = float(np.sum(resid**2))
    return sse, float(a), float(b)


def _evaluate_bayesian_lag_posterior(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    tracer: str,
    ion_order: Optional[List[str]],
    tau_grid: np.ndarray,
    k_grid: np.ndarray,
    tau_prior_mu: float,
    tau_prior_sigma: float,
    min_pairs: int,
) -> Tuple[Dict[str, object], List[str]]:
    flags: List[str] = []
    u_times = np.array(
        [s.timestamp.timestamp() / 86400.0 for s in node_u.samples], dtype=float
    )
    v_times = np.array(
        [s.timestamp.timestamp() / 86400.0 for s in node_v.samples], dtype=float
    )
    u_values = _extract_tracer_series(node_u, tracer, ion_order)
    v_values = _extract_tracer_series(node_v, tracer, ion_order)
    if u_values is None or v_values is None:
        return {"ok": False, "reason": "missing_tracer"}, ["missing_tracer"]
    if u_values.size < 3 or v_values.size < 3:
        return {"ok": False, "reason": "insufficient_data"}, ["insufficient_data"]

    # For each tau, build aligned pairs (v at t, u at t-tau) by interpolation.
    y = v_values.astype(float)
    y_mean = float(np.mean(y))
    sst = float(np.sum((y - y_mean) ** 2)) or 1e-12

    log_post = np.full((tau_grid.size,), -np.inf, dtype=float)
    best_map: Dict[str, object] = {"sse": float("inf")}
    details_by_k: List[Dict[str, float]] = []

    for k in k_grid:
        sse_list = []
        n_list = []
        for tau in tau_grid:
            # Determine which v times can be paired.
            target = v_times - float(tau)
            mask = (target >= u_times[0]) & (target <= u_times[-1])
            if int(np.sum(mask)) < int(min_pairs):
                sse_list.append(float("inf"))
                n_list.append(int(np.sum(mask)))
                continue
            u_interp = np.interp(target[mask], u_times, u_values).astype(float)
            x = u_interp * float(np.exp(-float(k) * float(tau)))
            y_sub = y[mask]
            sse, a, b = _ols_sse_with_intercept(x, y_sub)
            sse_list.append(sse)
            n_list.append(int(y_sub.size))

            if sse < float(best_map.get("sse", float("inf"))):
                best_map = {
                    "tau_map_days": float(tau),
                    "k_map": float(k),
                    "sse": float(sse),
                    "a": float(a),
                    "b": float(b),
                    "n_pairs": int(y_sub.size),
                }

        sse_arr = np.array(sse_list, dtype=float)
        # Marginal likelihood under Jeffreys prior for (a,b,σ): p(y|tau,k) ∝ SSE^{-(n-p)/2}
        # Use n_pairs(tau) to set exponent; p=2.
        log_like = np.full_like(sse_arr, -np.inf, dtype=float)
        for idx, (sse, n_pairs) in enumerate(zip(sse_arr, n_list)):
            if not np.isfinite(sse) or n_pairs <= 2:
                continue
            dof = max(1, int(n_pairs) - 2)
            log_like[idx] = -0.5 * float(dof) * float(np.log(sse + 1e-12))

        # Add prior on tau (truncated normal on [0, max])
        log_prior = _bayes_truncated_normal_prior(
            tau_grid, tau_prior_mu, tau_prior_sigma
        )
        log_joint = log_like + log_prior

        # Log-sum-exp across k (marginalize k) by accumulating in log_post
        m = np.maximum(log_post, log_joint)
        log_post = m + np.log(np.exp(log_post - m) + np.exp(log_joint - m))

        # track per-k best for details
        finite = np.isfinite(log_joint)
        if np.any(finite):
            idx_best = int(np.argmax(log_joint))
            details_by_k.append(
                {"k": float(k), "tau_map_days": float(tau_grid[idx_best])}
            )

    if not np.any(np.isfinite(log_post)):
        flags.append("bayes_failed")
        return {"ok": False, "reason": "bayes_failed"}, flags

    # Normalize posterior
    log_post = log_post - float(np.max(log_post))
    p = np.exp(log_post)
    total = float(np.sum(p)) or 1.0
    p = p / total

    mean_tau = float(np.sum(tau_grid * p))
    var_tau = float(np.sum(((tau_grid - mean_tau) ** 2) * p))
    std_tau = float(np.sqrt(max(0.0, var_tau)))
    cdf = np.cumsum(p)
    q025 = float(np.interp(0.025, cdf, tau_grid))
    q975 = float(np.interp(0.975, cdf, tau_grid))

    idx_map = int(np.argmax(p))
    tau_map = float(tau_grid[idx_map])
    post_mass_near_map = float(
        np.sum(p[(tau_grid >= max(0.0, tau_map - 2.0)) & (tau_grid <= tau_map + 2.0)])
    )
    if post_mass_near_map < 0.05:
        flags.append("bayes_flat_posterior")

    # R2 at MAP using best_map SSE if it matches tau_map closely; otherwise compute a quick one.
    sse_map = float(best_map.get("sse", float("inf")))
    r2 = 1.0 - (sse_map / sst) if np.isfinite(sse_map) else 0.0

    return (
        {
            "ok": True,
            "tau_mean_days": mean_tau,
            "tau_std_days": std_tau,
            "tau_ci_low_days": q025,
            "tau_ci_high_days": q975,
            "tau_map_days": tau_map,
            "posterior_mass_near_map": post_mass_near_map,
            "r2_map": float(r2),
            "best_map": best_map,
            "per_k_map": details_by_k,
        },
        flags,
    )


def _estimate_residence_time_bayesian_lag_consensus(
    node_u: TemporalNode,
    node_v: TemporalNode,
    *,
    tracer_ion: str,
    ion_order: Optional[List[str]],
    hydraulic_params: Optional[Dict[str, float]],
) -> Tuple[float, float, str, Dict[str, object], List[str]]:
    tracers = _parse_tracer_candidates(tracer_ion) or [tracer_ion]
    details: Dict[str, object] = {"candidates": {}}
    flags: List[str] = []

    # Prior on tau: prefer explicit priors (e.g., from MODPATH), else derive from Darcy/3D physics.
    tau_prior_mu = (hydraulic_params or {}).get("tau_prior_mu_days")
    tau_prior_sigma = (hydraulic_params or {}).get("tau_prior_sigma_days")
    tau_phy, unc_phy, used_phy = _estimate_residence_time_gradient(hydraulic_params)
    if tau_prior_mu is not None and tau_prior_sigma is not None:
        try:
            prior_mu = float(tau_prior_mu)
            prior_sigma = float(tau_prior_sigma)
            details["physics_prior"] = {
                "tau_mu_days": prior_mu,
                "tau_sigma_days": prior_sigma,
                "source": "explicit",
            }
        except (TypeError, ValueError):
            prior_mu = 0.0
            prior_sigma = (
                float((hydraulic_params or {}).get("bayes_lag_max_lag_days", 365.0))
                / 2.0
            )
            flags.append("physics_prior_invalid")
            details["physics_prior"] = {
                "tau_mu_days": prior_mu,
                "tau_sigma_days": prior_sigma,
                "source": "invalid",
            }
    elif used_phy == "gradient" and tau_phy > 0:
        prior_mu = float(tau_phy)
        prior_sigma = float(unc_phy) * float(
            (hydraulic_params or {}).get("bayes_lag_prior_sigma_multiplier", 1.0)
        )
        details["physics_prior"] = {
            "tau_mu_days": prior_mu,
            "tau_sigma_days": prior_sigma,
            "source": "darcy",
        }
    else:
        prior_mu = 0.0
        prior_sigma = (
            float((hydraulic_params or {}).get("bayes_lag_max_lag_days", 365.0)) / 2.0
        )
        flags.append("physics_prior_missing")
        details["physics_prior"] = {
            "tau_mu_days": prior_mu,
            "tau_sigma_days": prior_sigma,
            "source": "none",
        }

    dt = float((hydraulic_params or {}).get("bayes_lag_grid_dt_days", 5.0))
    max_lag = float((hydraulic_params or {}).get("bayes_lag_max_lag_days", 365.0))
    min_pairs = int((hydraulic_params or {}).get("bayes_lag_min_pairs", 5))
    tau_grid = np.arange(0.0, max_lag + 1e-9, dt, dtype=float)

    k_max = float((hydraulic_params or {}).get("attenuation_k_max", 0.02))
    k_steps = int((hydraulic_params or {}).get("attenuation_k_steps", 6))
    k_steps = max(1, k_steps)
    k_grid = np.linspace(0.0, max(0.0, k_max), k_steps, dtype=float)

    # Gates reused
    lmwl_intercept = float((hydraulic_params or {}).get("lmwl_intercept", 10.0))
    lmwl_slope = float((hydraulic_params or {}).get("lmwl_slope", 8.0))
    isotope_weight, isotope_flags, isotope_metrics = _isotope_evaporation_gate(
        node_u, node_v, lmwl_intercept=lmwl_intercept, lmwl_slope=lmwl_slope
    )
    if isotope_flags and isotope_flags != ["isotope_missing"]:
        flags.extend(isotope_flags)
    details["isotope_gate"] = {
        "weight": isotope_weight,
        "metrics": isotope_metrics,
        "flags": isotope_flags,
    }

    cl_weight, cl_flags, cl_metrics = _chloride_nonconservative_gate(
        node_u, node_v, ion_order=ion_order
    )
    if cl_flags and cl_flags != ["cl_missing"]:
        flags.extend(cl_flags)
    details["chloride_gate"] = {
        "weight": cl_weight,
        "metrics": cl_metrics,
        "flags": cl_flags,
    }

    accepted: List[Tuple[str, float, float, float]] = []
    tol = float((hydraulic_params or {}).get("agreement_tolerance", 0.4))
    min_r2 = float((hydraulic_params or {}).get("ttd_min_r2", 0.2))

    for tracer in tracers:
        fit, fit_flags = _evaluate_bayesian_lag_posterior(
            node_u,
            node_v,
            tracer=tracer,
            ion_order=ion_order,
            tau_grid=tau_grid,
            k_grid=k_grid,
            tau_prior_mu=prior_mu,
            tau_prior_sigma=prior_sigma,
            min_pairs=min_pairs,
        )
        candidate_flags = list(fit_flags)
        ok = bool(fit.get("ok"))
        tau = float(fit.get("tau_mean_days", 0.0)) if ok else 0.0
        unc = float(fit.get("tau_std_days", 0.0)) if ok else 0.0
        r2 = float(fit.get("r2_map", 0.0)) if ok else 0.0
        if ok and r2 < min_r2:
            candidate_flags.append("bayes_low_r2")

        gate_weight = 1.0
        if tracer in {"18O", "2H", "d18O", "d2H"}:
            gate_weight *= isotope_weight
        if tracer == "Cl":
            gate_weight *= cl_weight

        weight = gate_weight * max(0.0, r2) / (unc + 1.0)
        details["candidates"][tracer] = {
            **fit,
            "accepted": bool(ok),
            "weight": float(weight),
            "flags": candidate_flags,
        }
        if ok and r2 >= min_r2 and weight > 0 and tau > 0:
            accepted.append((tracer, tau, unc, float(weight)))

    if not accepted:
        if used_phy == "gradient" and tau_phy > 0:
            flags.append("tau_fallback_physics")
            return float(tau_phy), float(unc_phy), "gradient", details, flags
        flags.append("bayes_failed_all_tracers")
        return 0.0, 0.0, "bayesian_lag_failed", details, flags

    taus = [item[1] for item in accepted]
    tau_min = min(taus)
    tau_max = max(taus)
    rel_spread = (tau_max - tau_min) / max(1e-9, tau_max)
    details["consensus"] = {"rel_spread": float(rel_spread), "tolerance": float(tol)}

    accepted.sort(key=lambda item: item[3], reverse=True)
    best_tracer, best_tau, best_unc, _ = accepted[0]
    if rel_spread > tol and len(accepted) >= 2:
        flags.append("tau_ambiguous")
        spread_days = float(tau_max - tau_min)
        inflated_unc = float(np.sqrt(best_unc**2 + (0.5 * spread_days) ** 2))
        return (
            best_tau,
            inflated_unc,
            f"bayesian_lag_consensus({best_tracer})",
            details,
            flags,
        )

    total_w = sum(item[3] for item in accepted) or 1.0
    consensus_tau = sum(item[1] * item[3] for item in accepted) / total_w
    mean_unc = sum(item[2] * item[3] for item in accepted) / total_w
    spread = (
        float(np.std(np.array(taus, dtype=float), ddof=0)) if len(taus) > 1 else 0.0
    )
    consensus_unc = float(np.sqrt(mean_unc**2 + spread**2))
    details["consensus"].update(
        {
            "tau_days": float(consensus_tau),
            "uncertainty_days": float(consensus_unc),
            "accepted_tracers": [item[0] for item in accepted],
            "weights": {item[0]: float(item[3]) for item in accepted},
        }
    )
    return (
        float(consensus_tau),
        float(consensus_unc),
        "bayesian_lag_consensus",
        details,
        flags,
    )


def _estimate_residence_time_gradient(
    hydraulic_params: Optional[Dict[str, float]]
) -> Tuple[float, float, str]:
    """
    Estimate residence time using Darcy's law.

    τ = (distance * porosity) / (K * gradient)
    """
    if not hydraulic_params:
        return 0.0, 0.0, "gradient_no_params"

    distance = hydraulic_params.get("distance_m", 0.0)
    K = hydraulic_params.get("K_m_day", 1.0)
    gradient = hydraulic_params.get("gradient", 0.001)
    porosity = hydraulic_params.get("porosity", 0.2)

    if K <= 0 or gradient <= 0 or distance <= 0:
        return 0.0, 0.0, "gradient_invalid_params"

    # Darcy velocity: v = K * i
    # Pore velocity: v_p = v / n_e
    # Travel time: τ = distance / v_p = distance * porosity / (K * gradient)

    tau_days = (distance * porosity) / (K * gradient)

    # Uncertainty: assume 50% uncertainty on parameters
    uncertainty = tau_days * 0.5

    return float(tau_days), float(uncertainty), "gradient"


def _estimate_residence_time_tracer_decay(
    node_u: TemporalNode,
    node_v: TemporalNode,
    tracer_ion: str,
    ion_order: Optional[List[str]],
) -> Tuple[float, float, str]:
    """
    Estimate residence time using radioactive tracer decay.
    
    Infers Mean Residence Time (MRT) at each node using a lumped parameter model
    (default PFM) and input history, then computes travel time as difference.
    """
    if not _NUCLEAR_AVAILABLE:
        return 0.0, 0.0, "tracer_decay_module_missing"

    if not node_u.samples or not node_v.samples:
        return 0.0, 0.0, "tracer_decay_no_data"

    # Identify nuclide
    nuclide = get_nuclide(tracer_ion)
    if nuclide is None:
         return 0.0, 0.0, f"tracer_decay_unknown_nuclide({tracer_ion})"
    
    # Extract data
    u_vals = _extract_tracer_series(node_u, tracer_ion, ion_order)
    v_vals = _extract_tracer_series(node_v, tracer_ion, ion_order)
    
    if u_vals is None or v_vals is None:
        return 0.0, 0.0, "tracer_decay_missing_values"

    import datetime

    def get_mean_date_val(node, vals):
        if vals.size == 0: return None, None
        ts = np.array([s.timestamp.timestamp() for s in node.samples])
        # Use mean timestamp and value
        mean_ts = np.mean(ts)
        mean_val = np.mean(vals)
        
        # Convert to fractional year
        dt = datetime.datetime.fromtimestamp(mean_ts)
        # simplistic fractional year
        start_of_year = datetime.datetime(dt.year, 1, 1)
        year_fraction = (dt - start_of_year).total_seconds() / (365.25 * 24 * 3600)
        year = dt.year + year_fraction
        return year, mean_val

    u_year, u_conc = get_mean_date_val(node_u, u_vals)
    v_year, v_conc = get_mean_date_val(node_v, v_vals)
    
    if u_year is None or v_year is None:
        return 0.0, 0.0, "tracer_decay_no_dates"

    # Infer ages
    # Default error approx: 10% or 0.5 units
    u_sigma = max(0.5, u_conc * 0.1)
    v_sigma = max(0.5, v_conc * 0.1)
    
    # Infer MRT at each node
    # Using PFM is robust baseline as requested
    res_u = infer_age_from_tracer(u_conc, u_sigma, u_year, nuclide, model="PFM") 
    res_v = infer_age_from_tracer(v_conc, v_sigma, v_year, nuclide, model="PFM")
    
    age_u = res_u["tau_map_years"]
    age_v = res_v["tau_map_years"]
    
    # Travel time = difference in ages
    tau_years = age_v - age_u
    
    # Uncertainty propagation
    sigma_age_u = (res_u["tau_ci_high_years"] - res_u["tau_ci_low_years"]) / 4.0
    sigma_age_v = (res_v["tau_ci_high_years"] - res_v["tau_ci_low_years"]) / 4.0
    sigma_tau_years = np.sqrt(sigma_age_u**2 + sigma_age_v**2)
    
    # Convert to days
    tau_days = tau_years * 365.25
    sigma_tau_days = sigma_tau_years * 365.25
    
    method_str = f"tracer_decay({nuclide.symbol})"
    flags = []
    if res_u.get("multimodal") or res_v.get("multimodal"):
        flags.append("multimodal")
    
    if tau_days < 0:
        # Negative travel time implies downstream is younger than upstream
        # This conflicts with flow direction unless mixing with recent water occurred.
        tau_days = 0.0
        flags.append("neg_age_diff")
        
    if flags:
        method_str += f"_{'_'.join(flags)}"

    return float(tau_days), float(sigma_tau_days), method_str



def compute_residence_time_from_velocity(
    distance_m: float, velocity_m_day: float, porosity: float = 0.2
) -> float:
    """
    Simple residence time calculation from known velocity.

    Parameters
    ----------
    distance_m : float
        Distance along flow path
    velocity_m_day : float
        Darcy velocity in m/day
    porosity : float
        Effective porosity

    Returns
    -------
    float
        Residence time in days
    """
    if velocity_m_day <= 0:
        return 0.0

    pore_velocity = velocity_m_day / porosity
    tau_days = distance_m / pore_velocity

    return tau_days
