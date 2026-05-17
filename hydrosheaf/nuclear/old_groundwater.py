from __future__ import annotations

import json
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


def _finite_float(value: Any) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return out if math.isfinite(out) else float("nan")


def _valid_pmc(value: Any) -> bool:
    number = _finite_float(value)
    return math.isfinite(number) and 0.0 < number <= 130.0


def _parse_candidate_sequence(value: Any) -> list[Any]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        parsed = [token.strip() for token in text.split("|") if token.strip()]
    if isinstance(parsed, list):
        return parsed
    return [parsed]


def aggregate_c14_correction_candidates(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for sample_id, group in df.groupby("SampleID", dropna=False):
        plausible = group.copy()
        if "Model" in plausible.columns:
            model_text = plausible["Model"].astype(str).str.upper()
            plausible = plausible[~model_text.str.contains("NO 14C|NO MAJOR ION", regex=True, na=False)]
        corrected_series = plausible["Corrected_14C_sample_pmC"] if "Corrected_14C_sample_pmC" in plausible.columns else pd.Series(dtype=float)
        corrected_values = [
            float(value)
            for value in pd.to_numeric(corrected_series, errors="coerce").tolist()
            if _valid_pmc(value)
        ]
        a0_series = plausible["Corrected_Ao_pmC"] if "Corrected_Ao_pmC" in plausible.columns else pd.Series(dtype=float)
        a0_values = [
            float(value)
            for value in pd.to_numeric(a0_series, errors="coerce").tolist()
            if _valid_pmc(value)
        ]
        models = [str(value).strip() for value in plausible.get("Model", pd.Series(dtype=object)).tolist() if str(value).strip()]
        rows.append(
            {
                "SampleID": sample_id,
                "c14_candidate_count": int(max(len(corrected_values), len(a0_values), len(models))),
                "c14_candidate_corrected_pmc_json": json.dumps(corrected_values),
                "c14_candidate_a0_pmc_json": json.dumps(a0_values),
                "c14_candidate_models_json": json.dumps(models),
            }
        )
    return pd.DataFrame(rows)


def initial_c14_activity(observations: Mapping[str, Any]) -> float:
    corrected_a0 = _finite_float(observations.get("corrected_a0_pmc"))
    if _valid_pmc(corrected_a0):
        return corrected_a0
    return 100.0


def _effective_c14_candidates(
    candidate_corrected_pmcs: Any,
    candidate_initial_pmcs: Any,
    candidate_models: Any,
) -> tuple[list[float], list[float], list[str]]:
    corrected = [_finite_float(value) for value in _parse_candidate_sequence(candidate_corrected_pmcs)]
    corrected = [value for value in corrected if _valid_pmc(value)]
    initial = [_finite_float(value) for value in _parse_candidate_sequence(candidate_initial_pmcs)]
    initial = [value for value in initial if _valid_pmc(value)]
    models = [str(value).strip() for value in _parse_candidate_sequence(candidate_models) if str(value).strip()]
    return corrected, initial, models


def prepare_c14_observation(
    observations: Mapping[str, Any],
    *,
    mode: str = "selected",
    candidate_corrected_pmcs: Any = None,
    candidate_initial_pmcs: Any = None,
    candidate_models: Any = None,
) -> tuple[dict[str, Any], float, dict[str, Any]]:
    out = dict(observations)
    corrected_candidates, initial_candidates, models = _effective_c14_candidates(
        candidate_corrected_pmcs,
        candidate_initial_pmcs,
        candidate_models,
    )
    diagnostics: dict[str, Any] = {
        "c14_strategy": mode,
        "c14_candidate_count": int(max(len(corrected_candidates), len(initial_candidates), len(models))),
        "c14_candidate_models": "|".join(models),
        "c14_candidate_spread_pmc": float(np.ptp(corrected_candidates)) if len(corrected_candidates) >= 2 else 0.0,
    }

    if mode == "raw":
        initial = 100.0
        out["c14_initial_pmc"] = initial
        diagnostics["c14_effective_source"] = "raw"
        diagnostics["c14_effective_pmc"] = _finite_float(out.get("c14_pmc"))
        return out, initial, diagnostics

    if mode == "fixed_a0":
        initial = 100.0
        out["c14_initial_pmc"] = initial
        diagnostics["c14_effective_source"] = "fixed_a0"
        diagnostics["c14_effective_pmc"] = _finite_float(out.get("c14_pmc"))
        return out, initial, diagnostics

    if mode == "ensemble" and corrected_candidates:
        weights = []
        for m in models:
            m_upper = m.upper()
            if "PHREEQC" in m_upper: w = 8.0
            elif "NETPATH" in m_upper: w = 7.0
            elif "REV" in m_upper or "F&G" in m_upper: w = 5.0
            elif "TAMERS" in m_upper: w = 3.0
            elif "SCALING FACTOR" in m_upper: w = 1.0
            else: w = 0.5
            weights.append(w)
        if sum(weights[:len(corrected_candidates)]) > 0:
            out["c14_pmc"] = float(np.average(corrected_candidates, weights=weights[:len(corrected_candidates)]))
            initial = float(np.average(initial_candidates, weights=weights[:len(initial_candidates)])) if initial_candidates else initial_c14_activity(out)
        else:
            out["c14_pmc"] = float(np.median(corrected_candidates))
            initial = float(np.median(initial_candidates)) if initial_candidates else initial_c14_activity(out)

        out["c14_initial_pmc"] = initial
        diagnostics["c14_effective_source"] = "ensemble_weighted"
        diagnostics["c14_effective_pmc"] = float(out["c14_pmc"])
        diagnostics["c14_effective_a0_pmc"] = initial
        return out, initial, diagnostics

    corrected_c14 = _finite_float(out.get("corrected_c14_pmc"))
    if _valid_pmc(corrected_c14):
        out["c14_pmc"] = corrected_c14
        initial = 100.0
        out["c14_initial_pmc"] = initial
        diagnostics["c14_effective_source"] = "selected_corrected"
        diagnostics["c14_effective_pmc"] = corrected_c14
        return out, initial, diagnostics

    initial = initial_c14_activity(out)
    out["c14_initial_pmc"] = initial
    diagnostics["c14_effective_source"] = "raw_fallback"
    diagnostics["c14_effective_pmc"] = _finite_float(out.get("c14_pmc"))
    diagnostics["c14_effective_a0_pmc"] = initial
    return out, initial, diagnostics


def calibrated_he4_uncertainty_fraction(observations: Mapping[str, Any]) -> float:
    calibration_fields = [
        _finite_float(observations.get("he4_u_ppm")),
        _finite_float(observations.get("he4_th_ppm")),
        _finite_float(observations.get("he4_porosity")),
        _finite_float(observations.get("he4_bulk_density")),
    ]
    present = sum(math.isfinite(value) for value in calibration_fields)
    if present >= 4:
        return 0.25
    if present >= 2:
        return 0.35
    return 0.50


def apply_he4_uncertainty_mode(
    observations: Mapping[str, Any],
    *,
    mode: str = "calibrated",
) -> tuple[dict[str, Any], dict[str, Any]]:
    out = dict(observations)
    diagnostics = {
        "he4_uncertainty_mode": mode,
        "he4_rate_uncertainty_fraction": float("nan"),
        "he4_sigma_effective_ccpg": _finite_float(out.get("he4_sigma_ccpg")),
    }
    if mode != "calibrated_uncertainty":
        return out, diagnostics

    he4 = _finite_float(out.get("he4_ccpg"))
    background = _finite_float(out.get("he4_background_ccpg"))
    sigma = _finite_float(out.get("he4_sigma_ccpg"))
    if not math.isfinite(background):
        background = 0.0
    if not math.isfinite(sigma):
        sigma = 0.0
    if not math.isfinite(he4) or he4 <= background:
        diagnostics["he4_rate_uncertainty_fraction"] = calibrated_he4_uncertainty_fraction(out)
        return out, diagnostics

    rate_fraction = calibrated_he4_uncertainty_fraction(out)
    excess = max(he4 - background, 0.0)
    sigma_effective = max(sigma, excess * rate_fraction, 1.0e-12)
    out["he4_sigma_ccpg"] = sigma_effective
    diagnostics["he4_rate_uncertainty_fraction"] = rate_fraction
    diagnostics["he4_sigma_effective_ccpg"] = sigma_effective
    return out, diagnostics


def _c14_apparent_age_years(c14_pmc: Any, initial_pmc: Any) -> float:
    corrected = _finite_float(c14_pmc)
    initial = _finite_float(initial_pmc)
    if not _valid_pmc(corrected) or not _valid_pmc(initial) or corrected > initial:
        return float("nan")
    return -8033.0 * math.log(corrected / initial)


def _he4_apparent_age_years(he4_ccpg: Any, background_ccpg: Any, accumulation_rate: Any) -> float:
    he4 = _finite_float(he4_ccpg)
    background = _finite_float(background_ccpg)
    rate = _finite_float(accumulation_rate)
    if not math.isfinite(background):
        background = 0.0
    if not math.isfinite(he4) or not math.isfinite(rate) or rate <= 0 or he4 <= background:
        return float("nan")
    return (he4 - background) / rate


def diagnose_old_groundwater_constraints(
    observations: Mapping[str, Any],
    *,
    c14_initial_pmc: float | None = None,
) -> dict[str, Any]:
    initial = c14_initial_pmc if c14_initial_pmc is not None else observations.get("c14_initial_pmc")
    c14_age = _c14_apparent_age_years(observations.get("c14_pmc"), initial)
    he4_age = _he4_apparent_age_years(
        observations.get("he4_ccpg"),
        observations.get("he4_background_ccpg"),
        observations.get("he4_accumulation_rate_ccpg_per_year"),
    )
    has_c14 = math.isfinite(c14_age)
    has_he4 = math.isfinite(he4_age)
    if has_c14 and has_he4:
        log_gap = abs(math.log10(max(c14_age, 1.0)) - math.log10(max(he4_age, 1.0)))
        if log_gap <= 0.30:
            status = "agreement"
        elif log_gap <= 0.70:
            status = "tension"
        else:
            status = "conflict"
        case = "14C+4He"
        active = "14C|4He"
    elif has_c14:
        log_gap = float("nan")
        status = "single_tracer"
        case = "14C_only"
        active = "14C"
    elif has_he4:
        log_gap = float("nan")
        status = "single_tracer"
        case = "4He_only"
        active = "4He"
    else:
        log_gap = float("nan")
        status = "none"
        case = "none"
        active = ""
    return {
        "old_groundwater_case": case,
        "old_groundwater_active_constraints": active,
        "old_groundwater_constraint_count": int(has_c14) + int(has_he4),
        "old_groundwater_status": status,
        "old_groundwater_apparent_c14_age": c14_age,
        "old_groundwater_apparent_he4_age": he4_age,
        "old_groundwater_constraint_gap_log10": log_gap,
    }
