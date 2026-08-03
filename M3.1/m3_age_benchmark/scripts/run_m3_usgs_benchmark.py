"""
USGS Groundwater Age Benchmark for M3.
Evaluates Hydrosheaf's multi-tracer performance on National and Regional datasets.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple
from functools import lru_cache
import math
import json

import pandas as pd
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.log import get_logger
from hydrosheaf.nuclear.joint_lpm import (
    SUPPORTED_LPM_MODELS,
    age_fraction_predictions,
    compute_gated_bma_age,
    fit_lpm_model,
    fit_lpm_models,
)
from hydrosheaf.nuclear.dissolved_gas import fit_and_correct_dissolved_gases
from hydrosheaf.nuclear.multi_tracer import (
    calculate_tracer_reliability_weights,
    infer_multi_tracer_age,
)
from hydrosheaf.nuclear.tracer_inputs import (
    SiteInputContext,
    build_site_tracer_histories,
    site_input_history_metadata,
)
from hydrosheaf.nuclear.old_groundwater import (
    OldGroundwaterPrior,
    aggregate_c14_correction_candidates,
    apply_he4_uncertainty_mode,
    build_old_groundwater_priors,
    diagnose_old_groundwater_constraints,
    prepare_c14_observation,
    _lookup_oldwater_prior,
)

logger = get_logger("m3.usgs_benchmark")

# Benchmark Paths
BENCHMARK_DIR = Path(__file__).resolve().parents[1]
M2_USGS_DATA = REPO_ROOT / "M2" / "m2_benchmark" / "external" / "usgs_age" / "input" / "DataForNationalGroundwaterAge_1_1"
RESULT_DIR = BENCHMARK_DIR / "results"
QA_DIR = BENCHMARK_DIR / "docs"

M3_DEFAULT_AGE_STEPS = 90
FIT_ADEQUACY_MAX_STANDARDIZED_RMSE = 2.0
FIT_IDENTIFIABILITY_MAX_LOG10_AGE_SPAN = 0.5

def _make_site_context(row: Mapping[str, Any]) -> SiteInputContext:
    return SiteInputContext(
        site_id=str(row.get("site_id", "")),
        sample_year=float(row.get("sample_year", 2020.0)),
        latitude=_finite_float(row.get("lat")),
        longitude=_finite_float(row.get("lon")),
        study_unit=str(row.get("StudyUnit", "")),
        aquifer_group=str(row.get("AqGroup", "")),
        recharge_temperature_c=_finite_float(row.get("recharge_temperature_c")),
        elevation_m=_finite_float(row.get("elevation_m")),
    )


def _get_site_histories(row: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, str]]:
    ctx = _make_site_context(row)
    histories = build_site_tracer_histories(ctx)
    meta = site_input_history_metadata(ctx)
    return histories, meta

def _finite_float(val: Any) -> float | None:
    try:
        f = float(val)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None

def _decimal_year(date_str: Any) -> float | None:
    if pd.isna(date_str): return None
    try:
        dt = pd.to_datetime(date_str)
        return dt.year + (dt.dayofyear - 1) / 365.25
    except Exception:
        return None

def _parse_age(val):
    if pd.isna(val): return np.nan
    s = str(val).strip()
    if not s: return np.nan
    part = s.split()[0].split('(')[0]
    try: return float(part)
    except ValueError: return np.nan

def _log10_abs_error(est, ref):
    return abs(math.log10(est) - math.log10(ref)) if est > 0 and ref > 0 else np.nan

def _age_class(reference_age: float) -> str:
    if not math.isfinite(reference_age): return "unknown"
    if reference_age <= 50: return "modern_le_50"
    if reference_age <= 1000: return "intermediate_50_1000"
    if reference_age <= 30000: return "old_1000_30000"
    return "very_old_gt_30000"

def _deduplicate_c14_rows(df_c14):
    if "Corrected_Ao_pmC" in df_c14.columns:
        # Prefer rows with corrected values, then prefer PHREEQC model
        df_c14 = df_c14.assign(
            has_corr=df_c14["Corrected_Ao_pmC"].notna(),
            is_phreeqc=(df_c14["Model"].str.contains("PHREEQC", na=False))
        )
        df_c14 = df_c14.sort_values(
            ["SampleID", "has_corr", "is_phreeqc"], 
            ascending=[True, False, False]
        )
        df_c14 = df_c14.drop(columns=["has_corr", "is_phreeqc"])
    return df_c14.drop_duplicates("SampleID").reset_index(drop=True)

def _fit_age(fit):
    if hasattr(fit, "age_years"):
        return fit.age_years
    if hasattr(fit, "parameters"):
        params = fit.parameters
        if "mean_age_years" in params:
            return params["mean_age_years"]
        if "mean_age_1_years" in params and "mean_age_2_years" in params:
            fraction = float(params.get("binary_fraction", 0.5))
            return fraction * params["mean_age_1_years"] + (1.0 - fraction) * params["mean_age_2_years"]
    return np.nan


def _fit_identifiability(
    fit,
    *,
    reported_model_supported: bool | None = True,
) -> tuple[bool, str]:
    """Return a reference-free reportability decision for a fitted scalar age."""
    if not bool(getattr(fit, "converged", False)):
        return False, "fit_not_converged"

    n_tracers = int(getattr(fit, "n_tracers", 0) or 0)
    n_parameters = int(getattr(fit, "effective_n_params", 0) or 0)
    reasons: list[str] = []
    if reported_model_supported is False:
        reasons.append("reported_model_not_supported_no_fallback_age")
    if n_tracers <= n_parameters:
        reasons.append("nonpositive_residual_degrees_of_freedom")

    fit_rmse = _finite_float(getattr(fit, "rmse_standardized", np.nan))
    if fit_rmse is None or fit_rmse > FIT_ADEQUACY_MAX_STANDARDIZED_RMSE:
        reasons.append("gross_standardized_misfit")

    profile_span = _finite_float(getattr(fit, "age_profile_log10_span", np.nan))
    if profile_span is not None and profile_span > FIT_IDENTIFIABILITY_MAX_LOG10_AGE_SPAN:
        reasons.append("near_optimal_age_profile_spans_more_than_factor_3")

    return not reasons, "|".join(reasons)

def _nonidentifiability_reason(obs, lpm_tracers_mod):
    if "3He" in str(lpm_tracers_mod) and "trit" in str(lpm_tracers_mod).lower():
        if math.isnan(obs.get("tritium_TU", np.nan)):
            return "3H/3He requires finite positive tritium"
    return None


def _gas_likelihood_counts(obs: Mapping[str, Any]) -> dict[str, int]:
    """Count how many gas observations are marked with each likelihood type."""
    from hydrosheaf.nuclear.joint_lpm import build_lpm_tracer_observations
    fit_obs = build_lpm_tracer_observations(obs)
    counts = {"gaussian": 0, "upper_censored": 0, "lower_censored": 0, "contaminated_mixture": 0}
    for o in fit_obs:
        if o.tracer in {"SF6", "CFC11", "CFC12", "CFC113"}:
            counts[o.likelihood] = counts.get(o.likelihood, 0) + 1
    return counts

def _supported_reported_model(value: Any) -> str:
    model = str(value or "").strip().upper()
    return model if model in SUPPORTED_LPM_MODELS else "GA"

def _reported_tracer_tokens(value: Any) -> set[str]:
    """Map USGS reported tracer text into Hydrosheaf keys."""
    text = str(value or "").strip()
    tokens = set()
    if not text:
        return tokens
    
    mapping = {
        "3H": {"tritium_TU"},
        "3HE(TRIT)": {"tritium_TU", "he3_trit_TU"},
        "3H/3HE": {"tritium_TU", "he3_trit_TU"},
        "SF6": {"sf6_pptv"},
        "CFC-11": {"cfc11_pptv"},
        "CFC11": {"cfc11_pptv"},
        "CFC-12": {"cfc12_pptv"},
        "CFC12": {"cfc12_pptv"},
        "CFC-113": {"cfc113_pptv"},
        "CFC113": {"cfc113_pptv"},
        "14C": {"c14_pmc"},
        "4HE": {"he4_ccpg"},
        "39AR": {"ar39_pmc"},
        "85KR": {"kr85_pptv"},
    }
    
    # Split by common delimiters
    parts = [p.strip().upper() for p in text.replace(",", "|").replace(";", "|").split("|")]
    for p in parts:
        if p in mapping:
            tokens.update(mapping[p])
    return tokens


_REPORTED_LPM_TRACER_NAME_COLUMNS = (
    "LPM_Tracer_Name_01",
    "LPM_Tracer_Name_02",
    "LPM_Tracer_Name_03",
    "LPM_Tracer_Name_04",
    "LPM_Tracer_Name_05",
    "LPM_Tracer_Name_06",
    "LPM_Tracer_Name_07",
    # The source header repeats Name_09 and omits Name_08. Pandas mangles the
    # second duplicate as .1; these positional pairings preserve slots 8-9.
    "LPM_Tracer_Name_09",
    "LPM_Tracer_Name_09.1",
    "LPM_Tracer_Name_10",
)


def _canonical_reported_scale_tracer(value: Any) -> str:
    key = str(value or "").strip().upper().replace("-", "")
    return {
        "3H": "3H",
        "3HE(TRIT)": "3H/3He",
        "14C": "14C",
        "SF6": "SF6",
        "CFC11": "CFC11",
        "CFC12": "CFC12",
        "CFC113": "CFC113",
        "39AR": "39Ar",
        "85KR": "85Kr",
    }.get(key, "")


def _reported_lpm_metadata(df_lpm: pd.DataFrame) -> pd.DataFrame:
    """Extract auditable input-scaling and optimization metadata from Table 4."""
    if df_lpm.empty or "SampleID" not in df_lpm:
        return pd.DataFrame(columns=["SampleID"])

    work = df_lpm.copy()
    reported_age = pd.to_numeric(work.get("Rpt_TotAge_yrs"), errors="coerce")
    work["_valid_reported_age"] = reported_age.notna() & (reported_age > 0)
    work = work.sort_values("_valid_reported_age", ascending=False).drop_duplicates("SampleID")

    records: list[dict[str, Any]] = []
    for _, source in work.iterrows():
        scales: dict[str, float] = {}
        for slot, name_column in enumerate(_REPORTED_LPM_TRACER_NAME_COLUMNS, start=1):
            tracer = _canonical_reported_scale_tracer(source.get(name_column))
            scale = _finite_float(source.get(f"LPM_ScaleFact_{slot:02d}"))
            # Table 4 records zero for 4He because it uses a production-rate
            # model, not a recharge-input scale. Unsupported/zero scales are
            # deliberately omitted.
            if tracer and scale is not None and scale > 0:
                scales[tracer] = float(scale)

        records.append(
            {
                "SampleID": source["SampleID"],
                "reported_lpm_prediction_scale_factors_json": json.dumps(scales, sort_keys=True),
                "reported_lpm_optimization_parameters": source.get("LPM_Parms", ""),
                "reported_lpm_initial_age_c1_years": source.get("LPM_I_Age_C1"),
                "reported_lpm_initial_model_param1_c1": source.get("LPM_I_ModParm1_C1"),
                "reported_lpm_initial_model_param2_c1": source.get("LPM_I_ModParm2_C1"),
                "reported_lpm_initial_fraction_c1": source.get("LPM_I_Fraction_C1"),
                "reported_lpm_initial_age_c2_years": source.get("LPM_I_Age_C2"),
                "reported_lpm_initial_model_param1_c2": source.get("LPM_I_ModParm1_C2"),
                "reported_lpm_initial_model_param2_c2": source.get("LPM_I_ModParm2_C2"),
            }
        )
    return pd.DataFrame.from_records(records)


def _reported_prediction_scale_factors(row: Mapping[str, Any]) -> dict[str, float]:
    raw = row.get("reported_lpm_prediction_scale_factors_json", "{}")
    try:
        parsed = json.loads(raw) if isinstance(raw, str) else dict(raw or {})
    except (TypeError, ValueError, json.JSONDecodeError):
        return {}
    out: dict[str, float] = {}
    for tracer, value in parsed.items():
        scale = _finite_float(value)
        if scale is not None and scale > 0:
            out[str(tracer)] = float(scale)
    return out


def _reported_fit_configuration(
    row: Mapping[str, Any], model: str
) -> tuple[dict[str, float], tuple[str, ...]]:
    """Translate compatible Table 4 initial/fixed parameters into fit controls.

    Only DM parameterizations with a direct Hydrosheaf equivalent are used.
    The published *initial* configuration and optimization declaration are
    inputs; reported final fitted values and total age are never used here.
    """
    tokens = {
        token.strip().lower()
        for token in str(row.get("reported_lpm_optimization_parameters", "")).split(",")
        if token.strip()
    }
    model_key = str(model).upper()

    if model_key == "DM":
        age = _finite_float(row.get("reported_lpm_initial_age_c1_years"))
        dispersion = _finite_float(row.get("reported_lpm_initial_model_param1_c1"))
        if age is None or age <= 0 or dispersion is None or dispersion <= 0:
            return {}, ()
        template = {"mean_age_years": age, "dispersion": dispersion}
        free = []
        if "mean age" in tokens:
            free.append("mean_age_years")
        if "model parm 1" in tokens:
            free.append("dispersion")
        return template, tuple(free)

    if model_key == "BMM-DM-DM":
        values = {
            "mean_age_1_years": _finite_float(row.get("reported_lpm_initial_age_c1_years")),
            "dispersion": _finite_float(row.get("reported_lpm_initial_model_param1_c1")),
            "binary_fraction": _finite_float(row.get("reported_lpm_initial_fraction_c1")),
            "mean_age_2_years": _finite_float(row.get("reported_lpm_initial_age_c2_years")),
            "dispersion_2": _finite_float(row.get("reported_lpm_initial_model_param1_c2")),
        }
        if any(value is None for value in values.values()):
            return {}, ()
        template = {key: float(value) for key, value in values.items() if value is not None}
        if (
            template["mean_age_1_years"] <= 0
            or template["mean_age_2_years"] <= 0
            or template["dispersion"] <= 0
            or template["dispersion_2"] <= 0
            or not 0 < template["binary_fraction"] < 1
        ):
            return {}, ()
        free = []
        if "mean age" in tokens:
            free.append("mean_age_1_years")
        if "model parm 1" in tokens:
            free.append("dispersion")
        if "fraction" in tokens:
            free.append("binary_fraction")
        if "2nd mean age" in tokens:
            free.append("mean_age_2_years")
        return template, tuple(free)

    return {}, ()

def _apply_reported_tracer_mask(obs: dict[str, Any], lpm_tracers_mod: Any) -> dict[str, Any]:
    """Set non-reported tracer values and sigmas to None."""
    reported = _reported_tracer_tokens(lpm_tracers_mod)
    if not reported:
        return obs
        
    all_tracers = {
        "tritium_TU", "he3_trit_TU", "sf6_pptv", 
        "cfc11_pptv", "cfc12_pptv", "cfc113_pptv", 
        "c14_pmc", "he4_ccpg", "ar39_pmc", "kr85_pptv"
    }
    
    result = obs.copy()
    for tracer in all_tracers:
        if tracer not in reported:
            result[tracer] = None
            sigma_key = tracer.replace("_TU", "_sigma_TU").replace("_pptv", "_sigma_pptv").replace("_pmc", "_sigma_pmc").replace("_ccpg", "_sigma_ccpg")
            if sigma_key in result:
                result[sigma_key] = None
    return result

def _reported_uztt_years(row: Mapping[str, Any]) -> float:
    """Return reported UZ travel time, defaulting to 0.0."""
    val = _finite_float(row.get("Rpt_UZtt_yrs"))
    return max(0.0, val) if val is not None else 0.0

def _apply_age_target_mode(est_saturated_age: float, row: Mapping[str, Any], factors: Mapping[str, Any]) -> tuple[float, dict[str, Any]]:
    """Apply UZ travel-time corrections based on age_target_mode and uztt_mode."""
    uztt = _reported_uztt_years(row)
    uztt_mode = factors.get("uztt_mode", "ignore")
    age_target_mode = factors.get("age_target_mode", "saturated_only")
    
    est_total = est_saturated_age
    if uztt_mode == "add_reported":
        est_total = est_saturated_age + uztt
        
    diag = {
        "est_age_saturated_years": est_saturated_age,
        "est_age_total_years": est_total,
        "reported_uztt_years": uztt,
        "age_target_mode": age_target_mode,
        "uztt_mode": uztt_mode,
    }
    
    # If comparing to total age, the target age is Rpt_TotAge_yrs
    # If comparing to saturated age, the reference is Rpt_TotAge_yrs - Rpt_UZtt_yrs
    ref_total = _finite_float(row.get("Rpt_TotAge_yrs"))
    if age_target_mode == "reported_total":
        diag["est_age_multi"] = est_total
        if ref_total is not None:
            diag["reference_age_years"] = ref_total
    else:
        diag["est_age_multi"] = est_saturated_age
        if ref_total is not None:
            diag["reference_age_years"] = max(0.1, ref_total - uztt)

    return diag["est_age_multi"], diag

def _safe_join(values: Iterable[Any]) -> str:
    return "|".join(str(value) for value in values if str(value).strip())

def _ratio_within(est: Any, ref: Any, factor: float) -> bool:
    est_value = _finite_float(est)
    ref_value = _finite_float(ref)
    if est_value is None or ref_value is None or est_value <= 0 or ref_value <= 0:
        return False
    ratio = est_value / ref_value
    return (1.0 / factor) <= ratio <= factor

def _estimate_3h_age(obs: Mapping[str, Any], sample_year: float) -> float:
    tritium_obs = {
        "tritium_TU": obs.get("tritium_TU"),
        "tritium_sigma_TU": obs.get("tritium_sigma_TU"),
        "he3_trit_TU": obs.get("he3_trit_TU"),
        "he3_trit_sigma_TU": obs.get("he3_trit_sigma_TU"),
    }
    estimate = infer_multi_tracer_age(tritium_obs, sample_year=sample_year, use_helium4=False)
    age = _finite_float(estimate.get("age_years"))
    return float(age) if age is not None else np.nan

def _modern_fraction_proxy(age_years: Any) -> float:
    age = _finite_float(age_years)
    if age is None or age < 0:
        return np.nan
    return float(1.0 / (1.0 + age / 50.0))

def _proxy_age_coherence(obs: Mapping[str, Any], sample_year: float) -> tuple[float, str, str]:
    young_obs = {
        key: obs.get(key)
        for key in (
            "tritium_TU",
            "tritium_sigma_TU",
            "he3_trit_TU",
            "he3_trit_sigma_TU",
            "sf6_pptv",
            "sf6_sigma_pptv",
            "cfc11_pptv",
            "cfc11_sigma_pptv",
            "cfc12_pptv",
            "cfc12_sigma_pptv",
            "cfc113_pptv",
            "cfc113_sigma_pptv",
        )
    }
    estimate = infer_multi_tracer_age(young_obs, sample_year=sample_year, use_helium4=False)
    rows = [
        item for item in estimate.get("estimates", [])
        if _finite_float(item.get("age_years")) is not None and _finite_float(item.get("age_years")) >= 0
    ]
    if not rows:
        return np.nan, "", ""
    ages = [max(float(item["age_years"]), 0.1) for item in rows]
    names = [str(item.get("tracer", "")) for item in rows]
    values = [f"{name}:{age:.3g}" for name, age in zip(names, ages)]
    coherence = float(np.std(np.log10(ages))) if len(ages) >= 2 else 0.0
    return coherence, _safe_join(names), _safe_join(values)

def _apply_tracer_set(obs: Mapping[str, Any], tracer_set: str) -> dict[str, Any]:
    out = dict(obs)
    if tracer_set == "young_only":
        for key in ("c14_pmc", "c14_sigma_pmc", "he4_ccpg", "he4_sigma_ccpg"):
            out[key] = None
    elif tracer_set == "old_only":
        for key in (
            "tritium_TU",
            "tritium_sigma_TU",
            "he3_trit_TU",
            "he3_trit_sigma_TU",
            "sf6_pptv",
            "sf6_sigma_pptv",
            "cfc11_pptv",
            "cfc11_sigma_pptv",
            "cfc12_pptv",
            "cfc12_sigma_pptv",
            "cfc113_pptv",
            "cfc113_sigma_pptv",
        ):
            out[key] = None
    return out

def _apply_design_factors(obs, row, factors):
    out = dict(obs)
    factors = dict(factors) if factors else {}
    gas_mode = factors.get("gas_correction_mode", "corrected")
    factors["factor_gas_correction_mode"] = gas_mode
    
    gas_keys = ["sf6_pptv", "cfc11_pptv", "cfc12_pptv", "cfc113_pptv"]
    if gas_mode in {"corrected", "usgs_dgm"}:
        corrected_count = 0
        for k in gas_keys:
            dgm_key = f"dgm_{k}"
            dgm_value = _finite_float(row.get(dgm_key))
            if dgm_value is None:
                continue
            out[k] = dgm_value
            corrected_count += 1
            sigma_k = k.replace("_pptv", "_sigma_pptv")
            dgm_sigma = _finite_float(row.get(f"dgm_{sigma_k}"))
            if dgm_sigma is not None and dgm_sigma > 0:
                out[sigma_k] = dgm_sigma
        factors["gas_correction_input_source"] = (
            "usgs_table6_dgm" if corrected_count else "reported_atmospheric_equivalent"
        )
        factors["gas_correction_comparison_supported"] = True
    elif gas_mode == "raw":
        raw_count = 0
        for k in gas_keys:
            raw_key = f"raw_{k}"
            raw_value = _finite_float(row.get(raw_key))
            out[k] = raw_value
            if raw_value is not None:
                raw_count += 1
            sigma_k = k.replace("_pptv", "_sigma_pptv")
            raw_sigma_key = f"raw_{sigma_k}"
            raw_sigma = _finite_float(row.get(raw_sigma_key))
            out[sigma_k] = raw_sigma
        factors["gas_correction_input_source"] = "independent_raw_atmospheric_equivalent"
        factors["gas_correction_comparison_supported"] = raw_count > 0
        if raw_count == 0:
            factors["scenario_supported"] = False
            factors["scenario_support_reason"] = "raw atmospheric-equivalent gas observations unavailable"
    
    he4_mode = factors.get("he4_mode")
    if he4_mode == "disabled":
        out["he4_ccpg"] = None
        out["he4_sigma_ccpg"] = None
        
    return out, factors

def _choose_screened_gas_result(res_corr, res_raw):
    corr_coh = res_corr.get("young_gas_proxy_coherence", float('inf'))
    raw_coh = res_raw.get("young_gas_proxy_coherence", float('inf'))
    
    # Priority 1: Proxy Coherence (scientific consistency)
    if corr_coh < raw_coh - 0.1:
        return "usgs_dgm", "corrected preferred (proxy coherence)"
    if raw_coh < corr_coh - 0.05:
        return "raw", "raw preferred (proxy coherence)"
        
    # Priority 2: Fit Pathology
    corr_obj = res_corr.get("fit_objective", float('inf'))
    raw_obj = res_raw.get("fit_objective", float('inf'))
    if corr_obj > 100 and raw_obj < corr_obj:
        return "raw", "raw preferred (corrected pathological)"

    # Priority 3: Information Theory
    if res_raw.get("fit_aic", float('inf')) < res_corr.get("fit_aic", float('inf')) - 2.0:
         return "raw", "raw preferred (AIC)"
         
    return "usgs_dgm", "corrected preferred (default)"

def _has_screenable_gas_difference(row):
    for gas in ["sf6", "cfc11", "cfc12", "cfc113"]:
        k = f"{gas}_pptv"
        raw_k = f"raw_{k}"
        corrected = _finite_float(row.get(f"dgm_{k}"))
        if corrected is None:
            corrected = _finite_float(row.get(k))
        raw = _finite_float(row.get(raw_k))
        if corrected is not None and raw is not None and abs(corrected - raw) > 1e-5:
            return True
    return False

def load_benchmark_dataset(sources: str | None = None):
    """Return the analysis table for the benchmark.

    ``sources`` (or the ``M3_DATA_SOURCES`` environment variable) selects which
    of the releases declared in Section 2.2 are used:

    ``"national"`` - the national public-supply release only (historical default).
    ``"all"``      - all three releases, harmonised and de-duplicated by
                     :mod:`m3_data_sources`.

    Rows carry ``source_release``, ``has_coordinates`` and ``has_reported_lpm``
    so that scenarios which cannot run on a given release (graph families need
    coordinates, which the Western release does not publish; strict TracerLPM
    parity needs the reported LPM fields, which MRVA does not provide) are
    reported on their true support rather than silently averaged together.
    """
    import os

    national = load_usgs_national_dataset()
    choice = (sources or os.environ.get("M3_DATA_SOURCES", "national")).lower()
    if choice in ("national", "", "none"):
        national = national.copy()
        national["source_release"] = "national_public_supply"
        national["has_coordinates"] = national["lat"].notna() & national["lon"].notna()
        national["has_reported_lpm"] = True
        return national

    import m3_data_sources as _ds

    return _ds.load_combined_dataset(national)


def load_usgs_national_dataset():
    tracers_path = M2_USGS_DATA / "Table_3_Tracers.txt"
    ages_path = M2_USGS_DATA / "Table_2_Ages.txt"
    sites_path = M2_USGS_DATA / "Table_1_Sites.txt"
    c14_path = M2_USGS_DATA / "Table_7_Carbon14.txt"
    lpm_path = M2_USGS_DATA / "Table_4_LPM_ModOut.txt"
    diss_gas_path = M2_USGS_DATA / "Table_5_DissGas_ModOut.txt"
    tracer_output_path = M2_USGS_DATA / "Table_6_TrcOut.txt"
    
    if not tracers_path.exists(): return pd.DataFrame()
    
    df_tracers = pd.read_csv(tracers_path, sep="\t", na_values="na", low_memory=False).drop_duplicates("SampleID")
    df_ages = pd.read_csv(ages_path, sep="\t", na_values="na", low_memory=False).drop_duplicates("SampleID")
    df_sites = pd.read_csv(sites_path, sep="\t", na_values="na", low_memory=False).drop_duplicates("SampleID")
    
    df = pd.merge(df_tracers, df_ages, on="SampleID", how="inner")
    df = pd.merge(df, df_sites[["SampleID", "LatDD83", "LongDD83", "Depth_m", "StudyUnit", "AqGroup"]], on="SampleID", how="inner")
    
    if c14_path.exists():
        df_c14 = pd.read_csv(c14_path, sep="\t", na_values="na", encoding="latin-1", low_memory=False)
        df = pd.merge(df, _deduplicate_c14_rows(df_c14), on="SampleID", how="left", suffixes=("", "_c14"))
        df = pd.merge(df, aggregate_c14_correction_candidates(df_c14), on="SampleID", how="left")

    if lpm_path.exists():
        df_lpm = pd.read_csv(lpm_path, sep="\t", na_values="na", encoding="latin-1", low_memory=False)
        he4_columns = [
            col for col in (
                "SampleID",
                "LPM_He_SolnRate_ccpgpyr",
                "LPM_U_ppm",
                "LPM_Th_ppm",
                "LPM_Porosity",
                "LPM_BulkDensity",
            )
            if col in df_lpm.columns
        ]
        df = pd.merge(
            df,
            df_lpm[he4_columns].drop_duplicates("SampleID"),
            on="SampleID",
            how="left",
        )
        df = pd.merge(
            df,
            _reported_lpm_metadata(df_lpm),
            on="SampleID",
            how="left",
        )
        
    if diss_gas_path.exists():
        df_diss = pd.read_csv(diss_gas_path, sep="\t", na_values="na", encoding="latin-1", low_memory=False).drop_duplicates("SampleID")
        df = pd.merge(df, df_diss, on="SampleID", how="left", suffixes=("", "_diss"))

    if tracer_output_path.exists():
        df_trc = pd.read_csv(tracer_output_path, sep="\t", na_values="na", encoding="latin-1", low_memory=False)
        tracer_columns = [
            col for col in (
                "SampleID",
                "RecordNum",
                "CFC11_s_pmolkg", "CFC11_s_Err_pmolkg", "TRC_CFC11_pptv", "TRC_CFC11_Err_pptv",
                "CFC12_s_pmolkg", "CFC12_s_Err_pmolkg", "TRC_CFC12_pptv", "TRC_CFC12_Err_pptv",
                "CFC113_s_pmolkg", "CFC113_s_Err_pmolkg", "TRC_CFC113_pptv", "TRC_CFC113_Err_pptv",
                "SF6_s_fmolkg", "SF6_s_Err_fmolkg", "TRC_SF6_pptv", "TRC_SF6_Err_pptv",
            )
            if col in df_trc.columns
        ]
        merge_keys = ["SampleID", "RecordNum"] if "RecordNum" in df.columns and "RecordNum" in tracer_columns else ["SampleID"]
        df = pd.merge(
            df,
            df_trc[tracer_columns].drop_duplicates(merge_keys),
            on=merge_keys,
            how="left",
        )

    mapping = {
        "SampleID": "site_id",
        "3H_TU": "tritium_TU", "3H_err_TU": "tritium_sigma_TU",
        "3He_trit_TU": "he3_trit_TU", "3He_trit_err_TU": "he3_trit_sigma_TU",
        "SF6_pptv": "sf6_pptv", "SF6_err_pptv": "sf6_sigma_pptv",
        "CFC-11_pptv": "cfc11_pptv", "CFC-11_err_pptv": "cfc11_sigma_pptv",
        "CFC-12_pptv": "cfc12_pptv", "CFC-12_err_pptv": "cfc12_sigma_pptv",
        "CFC-113_pptv": "cfc113_pptv", "CFC-113_err_pptv": "cfc113_sigma_pptv",
        "14C_pmC": "c14_pmc", "14C_err_pmC": "c14_sigma_pmc",
        "Corrected_14C_sample_pmC": "corrected_c14_pmc",
        "Corrected_Ao_pmC": "corrected_a0_pmc",
        "4He_ccpg": "he4_ccpg", "4He_err_ccpg": "he4_sigma_ccpg",
        "Rpt_TotAge_yrs": "reference_age_years",
        "LatDD83": "lat", "LongDD83": "lon", "Depth_m": "depth_m",
        "DG_Model": "dgm_name",
        "DG_GasesMod": "dgm_gases_mod",
        "DG_Meas_He_ccpg": "raw_he4_ccpg", "DG_Meas_He_Err_ccpg": "raw_he4_sigma_ccpg",
        "DG_Meas_Ne_ccpg": "ne_ccstp_g", "DG_Meas_Ne_Err_ccpg": "ne_sigma_ccstp_g",
        "DG_Meas_Ar_ccpg": "ar_ccstp_g", "DG_Meas_Ar_Err_ccpg": "ar_sigma_ccstp_g",
        "DG_Meas_Kr_ccpg": "kr_ccstp_g", "DG_Meas_Kr_Err_ccpg": "kr_sigma_ccstp_g",
        "DG_Meas_Xe_ccpg": "xe_ccstp_g", "DG_Meas_Xe_Err_ccpg": "xe_sigma_ccstp_g",
        "DG_I_Salinity": "salinity_ppt",
        "LPM_He_SolnRate_ccpgpyr": "he4_accumulation_rate_ccpg_per_year",
        "LPM_U_ppm": "he4_u_ppm", "LPM_Th_ppm": "he4_th_ppm",
        "LPM_Porosity": "he4_porosity", "LPM_BulkDensity": "he4_bulk_density",
        "CFC11_s_pmolkg": "raw_cfc11_dissolved", "CFC11_s_Err_pmolkg": "raw_cfc11_dissolved_sigma",
        "CFC12_s_pmolkg": "raw_cfc12_dissolved", "CFC12_s_Err_pmolkg": "raw_cfc12_dissolved_sigma",
        "CFC113_s_pmolkg": "raw_cfc113_dissolved", "CFC113_s_Err_pmolkg": "raw_cfc113_dissolved_sigma",
        "SF6_s_fmolkg": "raw_sf6_dissolved", "SF6_s_Err_fmolkg": "raw_sf6_dissolved_sigma",
        "TRC_CFC11_pptv": "dgm_cfc11_pptv", "TRC_CFC11_Err_pptv": "dgm_cfc11_sigma_pptv",
        "TRC_CFC12_pptv": "dgm_cfc12_pptv", "TRC_CFC12_Err_pptv": "dgm_cfc12_sigma_pptv",
        "TRC_CFC113_pptv": "dgm_cfc113_pptv", "TRC_CFC113_Err_pptv": "dgm_cfc113_sigma_pptv",
        "TRC_SF6_pptv": "dgm_sf6_pptv", "TRC_SF6_Err_pptv": "dgm_sf6_sigma_pptv",
    }
    df.rename(columns=mapping, inplace=True)
    df["sample_year"] = df["SampleDate"].map(_decimal_year)
    
    dgm_columns = [f"dgm_{gas}_pptv" for gas in ("sf6", "cfc11", "cfc12", "cfc113")]
    dgm_available = df[[col for col in dgm_columns if col in df.columns]].notna().any(axis=1)
    df["dissolved_gas_correction"] = np.where(
        dgm_available,
        "usgs_table6_dgm",
        "reported_atmospheric_equivalent",
    )

    he4_rates = pd.to_numeric(
        df.get("he4_accumulation_rate_ccpg_per_year", pd.Series(np.nan, index=df.index)),
        errors="coerce",
    )
    direct_he4_rate = he4_rates.notna() & (he4_rates > 0)
    regional_medians = df.assign(he4_rate=he4_rates).groupby("StudyUnit")["he4_rate"].transform("median")
    regional_he4_rate = ~direct_he4_rate & regional_medians.notna() & (regional_medians > 0)
    df["he4_accumulation_rate_ccpg_per_year"] = he4_rates.fillna(regional_medians).fillna(1e-11)
    df["he4_source"] = np.select(
        [direct_he4_rate, regional_he4_rate],
        ["usgs_lpm_site_solution_rate", "usgs_lpm_study_unit_median"],
        default="screening_default_rate",
    )

    # Raw dissolved concentrations are retained with their native units, but
    # they are not mislabeled as atmospheric-equivalent pptv.  A paired raw-gas
    # benchmark is supported only when a source supplies an independent raw
    # atmospheric-equivalent field explicitly.
    for gas in ("sf6", "cfc11", "cfc12", "cfc113"):
        canonical = f"{gas}_pptv"
        sigma = f"{gas}_sigma_pptv"
        if f"raw_{canonical}" not in df.columns:
            df[f"raw_{canonical}"] = np.nan
        if f"raw_{sigma}" not in df.columns:
            df[f"raw_{sigma}"] = np.nan
    df["raw_gas_atmospheric_equivalent_available"] = df[
        [f"raw_{gas}_pptv" for gas in ("sf6", "cfc11", "cfc12", "cfc113")]
    ].notna().any(axis=1)

    return df.dropna(subset=["sample_year"])



def _extract_fit_parameters(fit_result):
    """Extract mean_age_years and secondary parameter from a JointLpmFit."""
    params = getattr(fit_result, "parameters", {}) or {}
    tau = float(params.get("mean_age_years", float("nan")))
    secondary_name = ""
    secondary_value = float("nan")
    for key, value in params.items():
        if key == "mean_age_years":
            continue
        if key in ("dispersion", "shape", "piston_fraction", "capture_fraction", "binary_fraction"):
            secondary_name = key
            secondary_value = float(value)
            break
    return tau, secondary_name, secondary_value


def _extract_model_aiccs(fits_list):
    """Return a JSON dict of model -> AICc for Akaike weight computation."""
    aiccs = {}
    for fit in (fits_list or []):
        model = getattr(fit, "model", "")
        aic = getattr(fit, "aic", float("nan"))
        if model and math.isfinite(aic):
            aiccs[model] = float(aic)
    return json.dumps(aiccs) if aiccs else "{}"

def _fit_prepared_benchmark_row(
    row,
    age_steps=M3_DEFAULT_AGE_STEPS,
    factors=None,
    model_strategy: str | None = None,
    oldwater_priors: Mapping[str, OldGroundwaterPrior] | None = None,
):
    factors = dict(factors or {})
    sample_year = float(row["sample_year"])
    ref_age = _parse_age(row.get("reference_age_years"))
    obs = {
        "sample_year": sample_year,
        "tritium_TU": row.get("tritium_TU"),
        "tritium_sigma_TU": row.get("tritium_sigma_TU"),
        "he3_trit_TU": row.get("he3_trit_TU"),
        "he3_trit_sigma_TU": row.get("he3_trit_sigma_TU"),
        "sf6_pptv": row.get("sf6_pptv"),
        "sf6_sigma_pptv": row.get("sf6_sigma_pptv"),
        "cfc11_pptv": row.get("cfc11_pptv"),
        "cfc11_sigma_pptv": row.get("cfc11_sigma_pptv"),
        "cfc12_pptv": row.get("cfc12_pptv"),
        "cfc12_sigma_pptv": row.get("cfc12_sigma_pptv"),
        "cfc113_pptv": row.get("cfc113_pptv"),
        "cfc113_sigma_pptv": row.get("cfc113_sigma_pptv"),
        "c14_pmc": row.get("c14_pmc"),
        "c14_sigma_pmc": row.get("c14_sigma_pmc"),
        "corrected_c14_pmc": row.get("corrected_c14_pmc"),
        "corrected_a0_pmc": row.get("corrected_a0_pmc"),
        "he4_ccpg": row.get("he4_ccpg"),
        "he4_sigma_ccpg": row.get("he4_sigma_ccpg"),
        "he4_accumulation_rate_ccpg_per_year": row.get("he4_accumulation_rate_ccpg_per_year"),
        "he4_background_ccpg": row.get("he4_background_ccpg", 0.0),
    }

    obs, factors = _apply_design_factors(obs, row, factors)

    # Phase 5: hierarchical old-water priors
    c14_mode = factors.get("c14_correction_mode", "selected")
    he4_mode = factors.get("he4_mode", "calibrated")
    prior = None
    if (c14_mode == "hierarchical" or he4_mode == "hierarchical") and oldwater_priors is not None:
        prior, prior_scope = _lookup_oldwater_prior(
            oldwater_priors,
            str(row.get("StudyUnit", "")),
            str(row.get("AqGroup", "")),
        )

    obs, c14_initial_pmc, c14_diag = prepare_c14_observation(
        obs,
        mode=c14_mode,
        candidate_corrected_pmcs=row.get("c14_candidate_corrected_pmc_json"),
        candidate_initial_pmcs=row.get("c14_candidate_a0_pmc_json"),
        candidate_models=row.get("c14_candidate_models_json"),
        prior=prior,
    )
    obs, he4_diag = apply_he4_uncertainty_mode(
        obs,
        mode=he4_mode,
        prior=prior,
    )

    prediction_scale_factors: dict[str, float] = {}
    if factors.get("reported_input_scaling", False):
        prediction_scale_factors = _reported_prediction_scale_factors(row)
        # Table 4's 14C scale is the source model's A0 adjustment. Apply it
        # directly to A0 (rather than as a post-prediction multiplier) so the
        # optimizer sees the same physical parameterization as TracerLPM.
        c14_scale = prediction_scale_factors.pop("14C", None)
        if c14_scale is not None:
            c14_initial_pmc = 100.0 * c14_scale
            c14_diag["c14_effective_source"] = "usgs_table4_lpm_scale_factor"
            c14_diag["c14_effective_a0_pmc"] = c14_initial_pmc
        applied_scales = {
            **prediction_scale_factors,
            **({"14C": c14_scale} if c14_scale is not None else {}),
        }
        c14_diag["reported_input_scaling_applied"] = bool(applied_scales)
        c14_diag["reported_input_scale_factors_json"] = json.dumps(
            applied_scales, sort_keys=True
        )
    else:
        c14_diag["reported_input_scaling_applied"] = False
        c14_diag["reported_input_scale_factors_json"] = "{}"
    
    # Phase 1: TracerLPM-parity tracer set masking
    reported_tracer_keys = _reported_tracer_tokens(row.get("LPM_TracersMod"))
    if factors.get("tracer_set") == "reported":
        obs = _apply_reported_tracer_mask(obs, row.get("LPM_TracersMod"))
    else:
        obs = _apply_tracer_set(obs, factors.get("tracer_set", "reported"))

    weighted_obs, diag = calculate_tracer_reliability_weights(obs, sample_year)
    proxy_coherence, proxy_names, proxy_values = _proxy_age_coherence(weighted_obs, sample_year)
    diag["young_gas_proxy_coherence"] = proxy_coherence

    strategy = model_strategy or factors.get("lpm_strategy", "reported")
    reported_model_str = str(row.get("LPM_Name", "") or "").strip()
    model = _supported_reported_model(reported_model_str)
    
    if strategy == "reported":
        diag["reported_model_supported"] = (model == reported_model_str.upper())
        if not diag["reported_model_supported"]:
            diag["reported_model_fallback_reason"] = f"Unsupported model: {reported_model_str}"

    parameter_template: dict[str, float] | None = None
    free_parameters: tuple[str, ...] | None = None
    if strategy == "reported" and factors.get("reported_parameter_configuration", False):
        configured_template, configured_free = _reported_fit_configuration(row, model)
        if configured_template:
            parameter_template = configured_template
            free_parameters = configured_free
            diag["reported_parameter_configuration_applied"] = True
            diag["reported_free_parameters_json"] = json.dumps(list(configured_free))
        else:
            diag["reported_parameter_configuration_applied"] = False
            diag["reported_free_parameters_json"] = "[]"
    else:
        diag["reported_parameter_configuration_applied"] = False
        diag["reported_free_parameters_json"] = "[]"
            
    use_helium4 = factors.get("he4_mode", "calibrated") != "disabled"
    fit_note = ""
    fits = []

    # Phase 3: Site-specific input histories
    histories, input_meta = _get_site_histories(row)

    # Phase 6: age-fraction constraints
    age_fraction_obs = None
    if factors.get("use_age_fractions"):
        age_fraction_obs = {
            "anthropocene": _finite_float(row.get("FracAnthropocene")),
            "holocene": _finite_float(row.get("FracHolocene")),
            "pleistocene": _finite_float(row.get("FracPleistocene")),
            "sigma_fraction": 0.10,
        }

    try:
        if strategy == "selection":
            fits = fit_lpm_models(
                weighted_obs,
                sample_year=sample_year,
                models=("PFM", "EM", "DM", "GA", "EPM", "PEM"),
                age_steps=age_steps,
                use_helium4=use_helium4,
                initial_c14_pmc=c14_initial_pmc,
                refine=True,
                histories=histories,
                age_fraction_obs=age_fraction_obs,
                prediction_scale_factors=prediction_scale_factors,
                parameter_template=parameter_template,
                free_parameters=free_parameters,
            )
            fit_multi = fits[0]
            # Phase 7: gated AICc/LOO BMA for selection scenarios
            bma = compute_gated_bma_age(fits)
            est_saturated_age = bma["age_years"]
            diag["bma_used"] = bma.get("bma_used")
            diag["bma_skip_reason"] = bma.get("bma_skip_reason", "")
            diag["bma_n_models"] = bma.get("bma_n_models")
            diag["bma_log_age_span"] = bma.get("bma_log_age_span")
            diag["bma_information_criterion"] = bma.get("bma_information_criterion")
            diag["top_model_age_years"] = bma.get("top_model_age_years")
            diag["bma_age_years"] = bma.get("age_years")
        else:
            fit_multi = fit_lpm_model(
                weighted_obs,
                sample_year=sample_year,
                model=model,
                age_steps=age_steps,
                use_helium4=use_helium4,
                initial_c14_pmc=c14_initial_pmc,
                refine=True,
                histories=histories,
                age_fraction_obs=age_fraction_obs,
                prediction_scale_factors=prediction_scale_factors,
                parameter_template=parameter_template,
                free_parameters=free_parameters,
            )
            est_saturated_age = _fit_age(fit_multi)
            fits = [fit_multi]
    except Exception as exc:
        fit_note = str(exc)
        fit_multi = fit_lpm_model(
            {},
            sample_year=sample_year,
            model=model,
            age_steps=age_steps,
            use_helium4=False,
            initial_c14_pmc=c14_initial_pmc,
            prediction_scale_factors=prediction_scale_factors,
            parameter_template=parameter_template,
            free_parameters=free_parameters,
        )
        est_saturated_age = np.nan

    fit_tau, fit_secondary_name, fit_secondary_value = _extract_fit_parameters(fit_multi)
    fit_model_aiccs = _extract_model_aiccs(fits)
    fit_identifiable, fit_identifiability_reason = _fit_identifiability(
        fit_multi,
        reported_model_supported=(
            diag.get("reported_model_supported") if strategy == "reported" else True
        ),
    )

    # Phase 6: compute predicted age fractions from best fit
    pred_fracs = {"anthropocene": np.nan, "holocene": np.nan, "pleistocene": np.nan}
    age_fraction_loss = np.nan
    if getattr(fit_multi, "converged", False):
        try:
            pred_fracs = age_fraction_predictions(
                getattr(fit_multi, "model", model),
                fit_multi.parameters,
            )
            if age_fraction_obs is not None:
                loss = 0.0
                sigma = 0.10
                for key in ("anthropocene", "holocene", "pleistocene"):
                    obs_val = age_fraction_obs.get(key)
                    if obs_val is not None and math.isfinite(obs_val):
                        loss += ((obs_val - pred_fracs.get(key, 0.0)) / sigma) ** 2
                age_fraction_loss = float(loss)
        except Exception:
            pass

    # Phase 2: UZ travel-time and age-target handling
    est_age, age_diag = _apply_age_target_mode(est_saturated_age, row, factors)
    est_age_raw = est_age
    if not fit_identifiable:
        est_age = np.nan
    age_diag["est_age_multi_raw"] = est_age_raw
    age_diag["est_age_multi"] = est_age
    diag.update(age_diag)
    ref_age_eff = diag.get("reference_age_years", ref_age)

    # est_age is already computed above via BMA or single model plus UZ corrections
    est_age_3h = _estimate_3h_age(weighted_obs, sample_year)
    old_diag = diagnose_old_groundwater_constraints(weighted_obs, c14_initial_pmc=c14_initial_pmc)
    error_multi = abs(est_age - ref_age_eff) if math.isfinite(est_age) and math.isfinite(ref_age_eff) else np.nan
    he4_value = _finite_float(weighted_obs.get("he4_ccpg"))
    he4_rate = _finite_float(weighted_obs.get("he4_accumulation_rate_ccpg_per_year"))

    # Phase 4: gas likelihood diagnostics
    gas_counts = _gas_likelihood_counts(weighted_obs)

    res = {
        "site_id": row["site_id"],
        # Provenance travels with every row so results can be reported per
        # release, and so scenarios can be restricted to their true support
        # (Western publishes no coordinates; MRVA has no reported-LPM table).
        "source_release": row.get("source_release", "national_public_supply"),
        "has_coordinates": row.get("has_coordinates", True),
        "has_reported_lpm": row.get("has_reported_lpm", True),
        "StudyUnit": row.get("StudyUnit", row.get("StudyUnit_x", "")),
        "AqGroup": row.get("AqGroup", row.get("AqGroup_x", "")),
        "Depth_m": row.get("depth_m", row.get("Depth_m", np.nan)),
        "LPM_Name": row.get("LPM_Name", ""),
        "LPM_TracersMod": row.get("LPM_TracersMod", ""),
        "Rpt_TotAge_yrs": row.get("reference_age_years", row.get("Rpt_TotAge_yrs", np.nan)),
        "Rept_TotAge_Err_yrs": row.get("Rept_TotAge_Err_yrs", np.nan),
        "Rpt_UZtt_yrs": row.get("Rpt_UZtt_yrs", np.nan),
        "Rpt_UZtt_Err_yrs": row.get("Rpt_UZtt_Err_yrs", np.nan),
        "Rpt_ChiSquare": row.get("Rpt_ChiSquare", np.nan),
        "Rpt_Probability": row.get("Rpt_Probability", np.nan),
        "FracAnthropocene": row.get("FracAnthropocene", np.nan),
        "FracHolocene": row.get("FracHolocene", np.nan),
        "FracPleistocene": row.get("FracPleistocene", np.nan),
        "AgeCat": row.get("AgeCat", ""),
        "Comments": row.get("Comments", ""),
        "ref_age": ref_age_eff,
        "est_age_multi": est_age,
        "est_age_multi_raw": est_age_raw,
        "est_age_saturated_years": diag.get("est_age_saturated_years"),
        "est_age_total_years": diag.get("est_age_total_years"),
        "reported_uztt_years": diag.get("reported_uztt_years"),
        "age_target_mode": diag.get("age_target_mode"),
        "uztt_mode": diag.get("uztt_mode"),
        "reported_model_supported": diag.get("reported_model_supported"),
        "reported_model_fallback_reason": diag.get("reported_model_fallback_reason"),
        "input_history_mode": input_meta.get("input_history_mode"),
        "input_history_region": input_meta.get("input_history_region"),
        "input_history_source": input_meta.get("input_history_source"),
        "est_age_3h": est_age_3h,
        "fit_aic": getattr(fit_multi, "aic", np.nan),
        "fit_objective": getattr(fit_multi, "objective", np.nan),
        "fit_rmse_standardized": getattr(fit_multi, "rmse_standardized", np.nan),
        "fit_effective_n_params": getattr(fit_multi, "effective_n_params", 0),
        "fit_age_profile_log10_span": getattr(fit_multi, "age_profile_log10_span", np.nan),
        "fit_age_profile_n_near_optimal": getattr(fit_multi, "age_profile_n_near_optimal", 0),
        "fit_identifiable": fit_identifiable,
        "fit_identifiability_reason": fit_identifiability_reason,
        "fit_mean_age_years": fit_tau,
        "fit_secondary_param_name": fit_secondary_name,
        "fit_secondary_param_value": fit_secondary_value,
        "fit_model_aiccs_json": fit_model_aiccs,
        "young_gas_proxy_coherence": proxy_coherence,
        "age_class": _age_class(ref_age_eff),
        "factor_gas_correction_mode": factors.get("factor_gas_correction_mode", "corrected"),
        "gas_correction_input_source": factors.get("gas_correction_input_source", ""),
        "gas_correction_comparison_supported": factors.get("gas_correction_comparison_supported", True),
        "scenario_supported": factors.get("scenario_supported", True),
        "scenario_support_reason": factors.get("scenario_support_reason", ""),
        "modern_fraction": _modern_fraction_proxy(est_age),
        "modern_age": est_age_3h,
        "old_age": old_diag.get("old_groundwater_apparent_c14_age", np.nan),
        "reported_model_name": reported_model_str,
        "reported_model": reported_model_str,
        "multi_model": getattr(fit_multi, "model", model),
        "model_strategy": strategy,
        "reported_tracer_mask": _safe_join(sorted(reported_tracer_keys)),
        "tracer_mode": row.get("LPM_TracersMod", ""),
        "tracer_set_used": factors.get("tracer_set", "reported"),
        "n_tracers": getattr(fit_multi, "n_tracers", 0),
        "c14_initial_pmc": c14_initial_pmc,
        "c14_correction_mode": factors.get("c14_correction_mode", "selected"),
        "c14_correction_model": row.get("Model", ""),
        "c14_correction_strategy": c14_diag.get("c14_strategy", factors.get("c14_correction_mode", "selected")),
        "c14_candidate_count": c14_diag.get("c14_candidate_count", 0),
        "c14_candidate_models": c14_diag.get("c14_candidate_models", ""),
        "c14_candidate_spread_pmc": c14_diag.get("c14_candidate_spread_pmc", 0.0),
        "c14_effective_source": c14_diag.get("c14_effective_source", ""),
        "c14_effective_pmc": c14_diag.get("c14_effective_pmc", np.nan),
        "reported_input_scaling_applied": c14_diag.get("reported_input_scaling_applied", False),
        "reported_input_scale_factors_json": c14_diag.get("reported_input_scale_factors_json", "{}"),
        "reported_lpm_optimization_parameters": row.get("reported_lpm_optimization_parameters", ""),
        "reported_lpm_initial_age_c1_years": row.get("reported_lpm_initial_age_c1_years", np.nan),
        "reported_lpm_initial_model_param1_c1": row.get("reported_lpm_initial_model_param1_c1", np.nan),
        "reported_lpm_initial_fraction_c1": row.get("reported_lpm_initial_fraction_c1", np.nan),
        "reported_lpm_initial_age_c2_years": row.get("reported_lpm_initial_age_c2_years", np.nan),
        "reported_parameter_configuration_applied": diag.get("reported_parameter_configuration_applied", False),
        "reported_free_parameters_json": diag.get("reported_free_parameters_json", "[]"),
        "he4_calibrated": bool(use_helium4 and he4_value is not None),
        "he4_uncertainty_mode": he4_diag.get("he4_uncertainty_mode", factors.get("he4_mode", "calibrated")),
        "he4_accumulation_rate": he4_rate if he4_rate is not None else np.nan,
        "he4_rate_uncertainty_fraction": he4_diag.get("he4_rate_uncertainty_fraction", np.nan),
        "he4_sigma_effective_ccpg": he4_diag.get("he4_sigma_effective_ccpg", np.nan),
        "he4_source": row.get("he4_source", ""),
        "dissolved_gas_correction": row.get("dissolved_gas_correction", ""),
        "dgm_name": row.get("dgm_name", ""),
        "young_gas_masking_mode": factors.get("gas_masking_mode", "soft_reliability"),
        "gas_likelihood_mode": "mixed" if sum(gas_counts.values()) > 0 else "gaussian",
        "censored_gas_count": gas_counts.get("upper_censored", 0) + gas_counts.get("lower_censored", 0),
        "contaminated_gas_count": gas_counts.get("contaminated_mixture", 0),
        "pred_frac_anthropocene": pred_fracs.get("anthropocene", np.nan),
        "pred_frac_holocene": pred_fracs.get("holocene", np.nan),
        "pred_frac_pleistocene": pred_fracs.get("pleistocene", np.nan),
        "age_fraction_loss": age_fraction_loss,
        "bma_used": np.nan,
        "bma_skip_reason": "",
        "bma_n_models": np.nan,
        "bma_log_age_span": np.nan,
        "bma_information_criterion": "",
        "top_model_age_years": np.nan,
        "bma_age_years": np.nan,
        "refinement_attempted": getattr(fit_multi, "refinement_attempted", False),
        "refinement_success": getattr(fit_multi, "refinement_success", False),
        "refinement_message": getattr(fit_multi, "refinement_message", ""),
        "young_gas_proxy_names": proxy_names,
        "young_gas_proxy_count": len([name for name in proxy_names.split("|") if name]),
        "young_gas_proxy_values": proxy_values,
        "young_gas_selected_mode": factors.get("factor_gas_correction_mode", "corrected"),
        "young_gas_screen_reason": "",
        "young_gas_corrected_proxy_coherence": np.nan,
        "young_gas_raw_proxy_coherence": np.nan,
        "young_gas_corrected_fit_aic": np.nan,
        "young_gas_raw_fit_aic": np.nan,
        "young_gas_corrected_fit_objective": np.nan,
        "young_gas_raw_fit_objective": np.nan,
        "young_gas_delta_aic_raw_minus_corrected": np.nan,
        "young_gas_delta_objective_raw_minus_corrected": np.nan,
        "fit_converged": bool(getattr(fit_multi, "converged", False)),
        "error_multi": error_multi,
        "log10_error": _log10_abs_error(est_age, ref_age_eff),
        "within_factor_2": _ratio_within(est_age, ref_age_eff, 2.0),
        "within_factor_10": _ratio_within(est_age, ref_age_eff, 10.0),
        "failure_reason": fit_note or fit_identifiability_reason or getattr(fit_multi, "note", ""),
    }
    res.update(diag)
    res.update(c14_diag)
    res.update(he4_diag)
    res.update(old_diag)
    return res

def fit_benchmark_row(
    row,
    age_steps=M3_DEFAULT_AGE_STEPS,
    factors=None,
    model_strategy: str | None = None,
    oldwater_priors: Mapping[str, OldGroundwaterPrior] | None = None,
    **kwargs,
):
    factors = dict(factors or {})
    strategy = model_strategy or factors.get("lpm_strategy", "reported")
    if factors.get("gas_correction_mode") == "screened_dgm":
        if not _has_screenable_gas_difference(row):
            corrected_factors = {**factors, "gas_correction_mode": "usgs_dgm"}
            res = _fit_prepared_benchmark_row(
                row,
                age_steps=age_steps,
                factors=corrected_factors,
                model_strategy=strategy,
                oldwater_priors=oldwater_priors,
            )
            res["young_gas_screen_reason"] = "no_screenable_difference"
            res["young_gas_selected_mode"] = "usgs_dgm"
            return res

        corrected_factors = {**factors, "gas_correction_mode": "usgs_dgm"}
        raw_factors = {**factors, "gas_correction_mode": "raw"}
        res_corr = _fit_prepared_benchmark_row(row, age_steps=age_steps, factors=corrected_factors, model_strategy=strategy, oldwater_priors=oldwater_priors)
        res_raw = _fit_prepared_benchmark_row(row, age_steps=age_steps, factors=raw_factors, model_strategy=strategy, oldwater_priors=oldwater_priors)
        selected, reason = _choose_screened_gas_result(res_corr, res_raw)
        res = dict(res_raw if selected == "raw" else res_corr)
        res.update(
            {
                "young_gas_selected_mode": selected,
                "young_gas_screen_reason": reason,
                "young_gas_corrected_proxy_coherence": res_corr.get("young_gas_proxy_coherence", np.nan),
                "young_gas_raw_proxy_coherence": res_raw.get("young_gas_proxy_coherence", np.nan),
                "young_gas_corrected_fit_aic": res_corr.get("fit_aic", np.nan),
                "young_gas_raw_fit_aic": res_raw.get("fit_aic", np.nan),
                "young_gas_corrected_fit_objective": res_corr.get("fit_objective", np.nan),
                "young_gas_raw_fit_objective": res_raw.get("fit_objective", np.nan),
                "young_gas_delta_aic_raw_minus_corrected": res_raw.get("fit_aic", np.nan) - res_corr.get("fit_aic", np.nan),
                "young_gas_delta_objective_raw_minus_corrected": res_raw.get("fit_objective", np.nan) - res_corr.get("fit_objective", np.nan),
            }
        )
        return res
    return _fit_prepared_benchmark_row(row, age_steps=age_steps, factors=factors, model_strategy=strategy, oldwater_priors=oldwater_priors)

if __name__ == "__main__":
    df = load_usgs_national_dataset()
    print(f"Verified: {len(df)} rows loaded.")
