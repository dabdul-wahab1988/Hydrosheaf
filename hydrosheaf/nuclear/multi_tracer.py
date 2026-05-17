"""Multi-tracer groundwater age inference helpers.

This module extends the single-nuclide tools with common young-water tracers
used in public groundwater-age datasets: 3H/3He, SF6, CFCs, and optional 4He
accumulation screening. The atmospheric gas histories are intentionally compact
default histories for reproducible screening. Users can replace them with
station-specific histories when available.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional

import numpy as np

from .input_history import InputHistory
from .invert import infer_age_from_tracer
from .nuclides import CARBON14, TRITIUM
from .tracer_inputs import standardize_gas_observations


@dataclass(frozen=True)
class TracerAgeEstimate:
    """A single tracer-derived apparent-age estimate."""

    tracer: str
    age_years: float
    ci_low_years: float
    ci_high_years: float
    method: str
    observed_value: float
    predicted_value: Optional[float] = None
    multimodal: bool = False
    note: str = ""


ATMOSPHERIC_TRACER_HISTORIES: Dict[str, Dict[str, List[float]]] = {
    # Mixing ratios are approximate Northern Hemisphere histories in pptv.
    # They provide deterministic screening estimates, not a substitute for
    # local recharge-temperature and excess-air corrections.
    "SF6": {
        "years": [1960, 1970, 1980, 1990, 2000, 2010, 2020, 2025],
        "values": [0.0, 0.03, 0.55, 2.1, 4.2, 7.0, 10.2, 12.0],
    },
    "CFC11": {
        "years": [1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2025],
        "values": [0.0, 2.0, 25.0, 135.0, 225.0, 260.0, 256.0, 240.0, 226.0, 220.0],
    },
    "CFC12": {
        "years": [1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2025],
        "values": [0.0, 8.0, 125.0, 285.0, 430.0, 525.0, 540.0, 528.0, 505.0, 496.0],
    },
    "CFC113": {
        "years": [1955, 1965, 1975, 1985, 1995, 2005, 2015, 2025],
        "values": [0.0, 1.0, 12.0, 55.0, 85.0, 82.0, 74.0, 68.0],
    },
    "85KR": {
        "years": [1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020, 2025],
        "values": [0.0, 5.0, 15.0, 25.0, 35.0, 50.0, 70.0, 95.0, 110.0],
    },
}


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def build_atmospheric_tracer_input(tracer: str) -> InputHistory:
    """Return the default atmospheric input history for a gas tracer."""

    key = tracer.upper().replace("-", "")
    if key not in ATMOSPHERIC_TRACER_HISTORIES:
        raise ValueError(
            "Unsupported atmospheric tracer. Expected one of: "
            + ", ".join(sorted(ATMOSPHERIC_TRACER_HISTORIES))
        )
    spec = ATMOSPHERIC_TRACER_HISTORIES[key]
    return InputHistory(np.array(spec["years"], dtype=float), np.array(spec["values"], dtype=float))


def infer_3h_3he_age(
    tritium_tu: float,
    helium3_trit_tu: float,
    tritium_sigma_tu: Optional[float] = None,
    helium3_sigma_tu: Optional[float] = None,
) -> Optional[TracerAgeEstimate]:
    """Infer apparent age from the tritiogenic 3He / tritium ratio."""

    tritium = _finite_float(tritium_tu)
    helium3 = _finite_float(helium3_trit_tu)
    if tritium is None or helium3 is None or tritium <= 0 or helium3 < 0:
        return None

    lambda_y = math.log(2.0) / TRITIUM.half_life_years
    ratio = helium3 / tritium
    age = math.log1p(ratio) / lambda_y

    tritium_sigma = max(_finite_float(tritium_sigma_tu) or 0.0, 0.08 * tritium, 0.05)
    helium3_sigma = max(_finite_float(helium3_sigma_tu) or 0.0, 0.12 * max(helium3, 0.05), 0.05)
    # Delta-method uncertainty for age = ln(1 + H/T) / lambda.
    denom = tritium + helium3
    d_age_d_he = 1.0 / (lambda_y * denom)
    d_age_d_tritium = -helium3 / (lambda_y * tritium * denom)
    sigma_age = math.sqrt((d_age_d_he * helium3_sigma) ** 2 + (d_age_d_tritium * tritium_sigma) ** 2)
    ci_low = max(0.0, age - 1.96 * sigma_age)
    ci_high = max(ci_low + 0.1, age + 1.96 * sigma_age)

    return TracerAgeEstimate(
        tracer="3H/3He",
        age_years=float(age),
        ci_low_years=float(ci_low),
        ci_high_years=float(ci_high),
        method="closed_system_ingrowth",
        observed_value=float(helium3),
        predicted_value=float(helium3),
    )


def infer_atmospheric_tracer_age(
    observed_pptv: float,
    sample_year: float,
    tracer: str,
    observed_sigma_pptv: Optional[float] = None,
    history: Optional[InputHistory] = None,
    max_age_years: float = 85.0,
    grid_step_years: float = 0.1,
) -> Optional[TracerAgeEstimate]:
    """Infer apparent recharge age from a dissolved-gas atmospheric equivalent.

    The observation must already be corrected to an atmospheric mixing ratio
    equivalent, as the USGS E1 table fields are.
    """

    observed = _finite_float(observed_pptv)
    sample = _finite_float(sample_year)
    if observed is None or sample is None or observed < 0:
        return None
    tracer_key = tracer.upper().replace("-", "")
    hist = history or build_atmospheric_tracer_input(tracer_key)
    min_year = max(float(hist.years.min()), sample - max_age_years)
    max_year = min(float(hist.years.max()), sample)
    if min_year >= max_year:
        return None

    recharge_years = np.arange(min_year, max_year + grid_step_years, grid_step_years)
    predicted = np.interp(recharge_years, hist.years, hist.values)
    sigma = max(_finite_float(observed_sigma_pptv) or 0.0, 0.08 * max(observed, 0.1), 0.05)
    log_like = -0.5 * ((observed - predicted) / sigma) ** 2
    max_ll = float(np.nanmax(log_like))
    posterior = np.exp(log_like - max_ll)
    if not np.isfinite(posterior).all() or posterior.sum() <= 0:
        return None
    posterior = posterior / posterior.sum()
    ages = sample - recharge_years
    idx = int(np.argmax(posterior))
    cdf = np.cumsum(posterior)
    low_idx = int(np.searchsorted(cdf, 0.025))
    high_idx = int(np.searchsorted(cdf, 0.975))

    peak_count = 0
    for i in range(1, len(posterior) - 1):
        if posterior[i] > posterior[i - 1] and posterior[i] > posterior[i + 1]:
            if posterior[i] > 0.05 * posterior[idx]:
                peak_count += 1

    ci_values = sorted([float(ages[low_idx]), float(ages[min(high_idx, len(ages) - 1)])])
    return TracerAgeEstimate(
        tracer=tracer_key,
        age_years=float(max(0.0, ages[idx])),
        ci_low_years=max(0.0, ci_values[0]),
        ci_high_years=max(ci_values[0] + 0.1, ci_values[1]),
        method="atmospheric_history_pfm",
        observed_value=float(observed),
        predicted_value=float(predicted[idx]),
        multimodal=peak_count > 1,
        note="atmospheric-equivalent screening age; replace with local history when available",
    )


def infer_helium4_accumulation_age(
    helium4_ccpg: float,
    background_ccpg: float = 4.6e-8,
    accumulation_rate_ccpg_per_year: float = 2.0e-11,
) -> Optional[TracerAgeEstimate]:
    """Estimate a screening accumulation age from radiogenic 4He excess.

    4He accumulation depends strongly on lithology, porosity, U/Th content, and
    crustal flux. This helper is off by default in the multi-tracer combiner.
    """

    helium = _finite_float(helium4_ccpg)
    if helium is None or helium <= background_ccpg or accumulation_rate_ccpg_per_year <= 0:
        return None
    excess = helium - background_ccpg
    age = excess / accumulation_rate_ccpg_per_year
    return TracerAgeEstimate(
        tracer="4He",
        age_years=float(age),
        ci_low_years=float(max(0.0, age * 0.25)),
        ci_high_years=float(age * 4.0),
        method="radiogenic_accumulation_screening",
        observed_value=float(helium),
        predicted_value=float(background_ccpg + age * accumulation_rate_ccpg_per_year),
        note="highly uncertain unless accumulation rate is site-calibrated",
    )


def _from_single_nuclide_result(
    tracer: str,
    obs_value: float,
    result: Mapping[str, Any],
) -> TracerAgeEstimate:
    return TracerAgeEstimate(
        tracer=tracer,
        age_years=float(result["tau_map_years"]),
        ci_low_years=float(result["tau_ci_low_years"]),
        ci_high_years=float(result["tau_ci_high_years"]),
        method=str(result.get("model_used", "lpm")),
        observed_value=float(obs_value),
        predicted_value=_finite_float(result.get("predicted_conc_at_map")),
        multimodal=bool(result.get("multimodal", False)),
    )


def combine_tracer_age_estimates(estimates: Iterable[TracerAgeEstimate]) -> Optional[Dict[str, Any]]:
    """Combine tracer ages in log-age space with uncertainty-width weights."""

    usable = [estimate for estimate in estimates if estimate.age_years >= 0]
    if not usable:
        return None

    log_ages: List[float] = []
    weights: List[float] = []
    lows: List[float] = []
    highs: List[float] = []
    for estimate in usable:
        age = max(estimate.age_years, 0.1)
        low = max(estimate.ci_low_years, 0.1)
        high = max(estimate.ci_high_years, low + 0.1)
        width = max(math.log10(high) - math.log10(low), 0.05)
        log_ages.append(math.log10(age))
        weights.append(1.0 / width)
        lows.append(low)
        highs.append(high)

    combined_log_age = float(np.average(np.array(log_ages), weights=np.array(weights)))
    return {
        "age_years": float(10**combined_log_age),
        "ci_low_years": float(min(lows)),
        "ci_high_years": float(max(highs)),
        "tracers": [estimate.tracer for estimate in usable],
        "method": "+".join(estimate.tracer for estimate in usable),
        "n_tracers": len(usable),
        "multimodal": any(estimate.multimodal for estimate in usable),
        "estimates": [estimate.__dict__ for estimate in usable],
    }


def infer_multi_tracer_age(
    observations: Mapping[str, Any],
    sample_year: float,
    use_helium4: bool = False,
) -> Dict[str, Any]:
    """Infer a combined apparent age from all supported tracer observations."""

    observations = standardize_gas_observations(observations)
    estimates: List[TracerAgeEstimate] = []
    gas_estimates: List[TracerAgeEstimate] = []
    skipped_estimates: List[Dict[str, Any]] = []
    tritium_modern = False
    c14_age_years: Optional[float] = None

    tritium = _finite_float(observations.get("tritium_TU"))
    if tritium is not None and tritium >= 0:
        tritium_modern = tritium >= 0.5
        tritium_sigma = max(_finite_float(observations.get("tritium_sigma_TU")) or 0.0, 0.15 * max(tritium, 0.1), 0.1)
        result = infer_age_from_tracer(
            tritium,
            tritium_sigma,
            sample_year,
            TRITIUM,
            model="EM",
            tau_range=(0.0, 200.0),
            search_steps=400,
        )
        estimates.append(_from_single_nuclide_result("3H", tritium, result))

    he3_age = infer_3h_3he_age(
        observations.get("tritium_TU"),
        observations.get("he3_trit_TU"),
        observations.get("tritium_sigma_TU"),
        observations.get("he3_trit_sigma_TU"),
    )
    if he3_age is not None:
        tritium_modern = True
        estimates.append(he3_age)

    c14 = _finite_float(observations.get("c14_pmc"))
    if c14 is not None and 0 < c14 <= 130:
        c14_sigma = max(_finite_float(observations.get("c14_sigma_pmc")) or 0.0, 0.05 * c14, 1.0)
        result = infer_age_from_tracer(
            c14,
            c14_sigma,
            sample_year,
            CARBON14,
            model="PFM",
            tau_range=(0.0, 50000.0),
            search_steps=900,
            initial_activity=100.0,
        )
        c14_estimate = _from_single_nuclide_result("14C", c14, result)
        c14_age_years = c14_estimate.age_years
        estimates.append(c14_estimate)

    gas_specs = [
        ("SF6", "sf6_pptv", "sf6_sigma_pptv"),
        ("CFC11", "cfc11_pptv", "cfc11_sigma_pptv"),
        ("CFC12", "cfc12_pptv", "cfc12_sigma_pptv"),
        ("CFC113", "cfc113_pptv", "cfc113_sigma_pptv"),
    ]
    for tracer, value_key, sigma_key in gas_specs:
        estimate = infer_atmospheric_tracer_age(
            observations.get(value_key),
            sample_year,
            tracer,
            observations.get(sigma_key),
        )
        if estimate is not None:
            gas_estimates.append(estimate)

    for estimate in gas_estimates:
        # SF6/CFCs are powerful young-water tracers, but they are also prone to
        # contamination and can represent a very small young fraction in mixed
        # samples. Do not let an uncorroborated gas tracer collapse an old 14C
        # age estimate unless 3H/3He/tritium evidence also supports young water.
        if c14_age_years is not None and c14_age_years > 1000.0 and not tritium_modern:
            skipped_estimates.append({**estimate.__dict__, "skip_reason": "uncorroborated_young_gas_with_old_14c"})
            continue
        estimates.append(estimate)

    if use_helium4:
        he4_age = infer_helium4_accumulation_age(observations.get("he4_ccpg"))
        if he4_age is not None:
            estimates.append(he4_age)

    combined = combine_tracer_age_estimates(estimates)
    if combined is None:
        empty = {
            "age_years": np.nan,
            "ci_low_years": np.nan,
            "ci_high_years": np.nan,
            "method": "no_supported_tracer",
            "n_tracers": 0,
            "tracers": [],
            "multimodal": False,
            "estimates": [],
            "skipped_estimates": skipped_estimates,
        }
        from .diagnostics import diagnose_tracer_disagreement

        empty["disagreement"] = diagnose_tracer_disagreement(empty)
        empty["flags"] = empty["disagreement"]["flags"]
        return empty
    combined["skipped_estimates"] = skipped_estimates
    from .diagnostics import diagnose_tracer_disagreement

    combined["disagreement"] = diagnose_tracer_disagreement(combined)
    combined["flags"] = combined["disagreement"]["flags"]
    return combined

def historical_max_concentration(tracer: str, sample_year: float) -> float:
    """Get the maximum historical concentration of a gas tracer up to the sample year."""
    try:
        history = build_atmospheric_tracer_input(tracer)
    except ValueError:
        return float('inf')
    years, values = np.asarray(history.years, dtype=float), np.asarray(history.values, dtype=float)
    if math.isfinite(sample_year):
        mask = years <= sample_year + 1.0e-6
        if bool(np.any(mask)):
            return float(np.max(values[mask]))
    return float(np.max(values))


def calculate_tracer_reliability_weights(
    observations: Mapping[str, Any],
    sample_year: float,
    reference_age: float = float("nan"),
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Calculate soft reliability weights for gas tracers based on historical limits
    and proxy coherence.
    """
    out = dict(observations)
    masked: list[str] = []
    reasons: list[str] = []

    if not math.isfinite(reference_age):
        independent_proxies: list[float] = []
        age_3h3he = float("nan")
        tritium = _finite_float(out.get("tritium_TU"))
        he3 = _finite_float(out.get("he3_trit_TU"))
        if tritium is not None and he3 is not None and tritium > 0 and he3 >= 0:
            lambda_y = math.log(2.0) / TRITIUM.half_life_years
            age_3h3he = math.log1p(he3 / tritium) / lambda_y
        
        if math.isfinite(age_3h3he) and age_3h3he >= 0.5:
            independent_proxies.append(float(age_3h3he))
        if independent_proxies:
            reference_age = float(np.median(independent_proxies))

    gas_specs = [
        ("SF6", "sf6_pptv"),
        ("CFC11", "cfc11_pptv"),
        ("CFC12", "cfc12_pptv"),
        ("CFC113", "cfc113_pptv"),
    ]
    for tracer, value_key in gas_specs:
        value = _finite_float(out.get(value_key))
        if value is None or value < 0:
            continue
        tracer_reasons: list[str] = []
        max_hist = historical_max_concentration(tracer, sample_year)
        
        tracer_weight = 1.0
        if value > max_hist * 1.02:
            tracer_weight *= 0.05
            tracer_reasons.append("above_historical_max_downweighted")
            
        proxy_age = float("nan")
        try:
            est = infer_atmospheric_tracer_age(value, sample_year, tracer)
            if est is not None:
                proxy_age = est.age_years
        except Exception:
            pass

        if (
            math.isfinite(reference_age)
            and reference_age >= 5.0
            and math.isfinite(proxy_age)
            and proxy_age <= max(0.5, reference_age * 0.2)
        ):
            tracer_weight *= 0.1
            tracer_reasons.append("excessively_modern_vs_independent_proxy_downweighted")
            
        if tracer_reasons:
            out[f"{tracer.lower()}_weight"] = tracer_weight
            masked.append(tracer)
            reasons.append(f"{tracer}:{'+'.join(tracer_reasons)}")
            
    return out, {
        "young_gas_masked_tracers": "|".join(masked),
        "young_gas_masked_reason": "|".join(reasons),
        "young_gas_masked_count": len(masked),
    }

