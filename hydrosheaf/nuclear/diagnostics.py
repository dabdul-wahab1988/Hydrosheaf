"""Diagnostics for nuclear-tracer age interpretation.

The helpers in this module keep manuscript-facing cautions explicit: tracer
ages can disagree, several LPM families can fit equally well, and some samples
are not identifiable from the available tracer suite.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from .joint_lpm import JointLpmFit, fit_lpm_models
from .multi_tracer import TracerAgeEstimate, combine_tracer_age_estimates


YOUNG_TRACERS = {"3H", "3H/3HE", "SF6", "CFC11", "CFC12", "CFC113"}
OLD_TRACERS = {"14C", "4HE"}


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _tracer_key(value: Any) -> str:
    return str(value or "").upper().replace("-", "")


def _estimate_from_mapping(row: Mapping[str, Any]) -> Optional[TracerAgeEstimate]:
    age = _finite_float(row.get("age_years"))
    low = _finite_float(row.get("ci_low_years"))
    high = _finite_float(row.get("ci_high_years"))
    observed = _finite_float(row.get("observed_value"))
    if age is None or low is None or high is None or observed is None:
        return None
    return TracerAgeEstimate(
        tracer=str(row.get("tracer", "")),
        age_years=age,
        ci_low_years=low,
        ci_high_years=high,
        method=str(row.get("method", "unknown")),
        observed_value=observed,
        predicted_value=_finite_float(row.get("predicted_value")),
        multimodal=bool(row.get("multimodal", False)),
        note=str(row.get("note", "")),
    )


def _log10_age(age: float) -> float:
    return math.log10(max(float(age), 0.1))


def diagnose_tracer_disagreement(
    age_result: Mapping[str, Any],
    *,
    young_age_threshold_years: float = 80.0,
    old_age_threshold_years: float = 1000.0,
    conflict_log10_threshold: float = 1.0,
    influence_log10_threshold: float = 0.25,
) -> Dict[str, Any]:
    """Score tracer-age disagreement and contamination/mixing flags.

    ``age_result`` is the dictionary returned by ``infer_multi_tracer_age``.
    The function is intentionally conservative: it does not decide that a
    young tracer is wrong, only that the combination needs explicit discussion.
    """

    estimates = [
        estimate
        for row in age_result.get("estimates", [])
        for estimate in [_estimate_from_mapping(row)]
        if estimate is not None
    ]
    skipped = list(age_result.get("skipped_estimates", []))
    tracer_rows: List[Dict[str, Any]] = []
    for estimate in estimates:
        tracer_rows.append(
            {
                "tracer": estimate.tracer,
                "age_years": float(estimate.age_years),
                "ci_low_years": float(estimate.ci_low_years),
                "ci_high_years": float(estimate.ci_high_years),
                "log10_age": _log10_age(estimate.age_years),
                "multimodal": bool(estimate.multimodal),
                "used": True,
            }
        )
    for row in skipped:
        age = _finite_float(row.get("age_years"))
        if age is None:
            continue
        tracer_rows.append(
            {
                "tracer": str(row.get("tracer", "")),
                "age_years": age,
                "ci_low_years": _finite_float(row.get("ci_low_years")),
                "ci_high_years": _finite_float(row.get("ci_high_years")),
                "log10_age": _log10_age(age),
                "multimodal": bool(row.get("multimodal", False)),
                "used": False,
                "skip_reason": str(row.get("skip_reason", "")),
            }
        )

    used_log_ages = [row["log10_age"] for row in tracer_rows if row.get("used")]
    conflict_score = max(used_log_ages) - min(used_log_ages) if len(used_log_ages) >= 2 else 0.0

    young_rows = [
        row
        for row in tracer_rows
        if _tracer_key(row.get("tracer")) in YOUNG_TRACERS
        and row["age_years"] <= young_age_threshold_years
    ]
    old_rows = [
        row
        for row in tracer_rows
        if _tracer_key(row.get("tracer")) in OLD_TRACERS
        and row["age_years"] >= old_age_threshold_years
    ]
    skipped_reasons = {str(row.get("skip_reason", "")) for row in skipped}

    influence: List[Dict[str, Any]] = []
    combined_age = _finite_float(age_result.get("age_years"))
    if combined_age is not None and len(estimates) >= 2:
        combined_log = _log10_age(combined_age)
        for idx, estimate in enumerate(estimates):
            leave_one_out = combine_tracer_age_estimates(
                other for j, other in enumerate(estimates) if j != idx
            )
            loo_age = _finite_float((leave_one_out or {}).get("age_years"))
            if loo_age is None:
                continue
            delta_log = _log10_age(loo_age) - combined_log
            influence.append(
                {
                    "tracer": estimate.tracer,
                    "leave_one_out_age_years": loo_age,
                    "delta_log10_age": float(delta_log),
                    "high_influence": abs(delta_log) >= influence_log10_threshold,
                }
            )

    flags = {
        "tracer_conflict_high": conflict_score >= conflict_log10_threshold,
        "young_old_contradiction": bool(young_rows and old_rows),
        "possible_gas_contamination": any(
            _tracer_key(row.get("tracer")) in {"SF6", "CFC11", "CFC12", "CFC113"}
            for row in young_rows
        )
        and bool(old_rows),
        "mixed_age_ambiguity": bool(young_rows and old_rows)
        or any(row.get("multimodal") for row in tracer_rows),
        "low_tracer_count": int(age_result.get("n_tracers", 0) or 0) < 2,
        "multimodal_tracer": any(row.get("multimodal") for row in tracer_rows),
        "uncorroborated_young_gas_with_old_14c": (
            "uncorroborated_young_gas_with_old_14c" in skipped_reasons
        ),
        "high_influence_tracer": any(row.get("high_influence") for row in influence),
    }

    return {
        "conflict_score_log10_age": float(conflict_score),
        "young_tracers": young_rows,
        "old_tracers": old_rows,
        "per_tracer": tracer_rows,
        "per_tracer_influence": influence,
        "flags": flags,
        "interpretation": _interpret_disagreement(flags),
    }


def _interpret_disagreement(flags: Mapping[str, bool]) -> str:
    if flags.get("young_old_contradiction"):
        return (
            "Young and old tracer systems imply incompatible apparent ages; "
            "report as mixed, contaminated, or non-identifiable unless resolved "
            "by independent dissolved-gas and hydrostratigraphic evidence."
        )
    if flags.get("tracer_conflict_high"):
        return "Tracer ages span more than one order of magnitude; avoid a single unqualified age claim."
    if flags.get("low_tracer_count"):
        return "Only one usable tracer constrains age; identifiability is weak."
    return "No high-severity tracer disagreement detected by the configured thresholds."


def _fit_mean_age(fit: JointLpmFit) -> Optional[float]:
    params = fit.parameters
    if not fit.converged:
        return None
    if fit.model.upper().startswith("BMM-"):
        age1 = _finite_float(params.get("mean_age_1_years"))
        age2 = _finite_float(params.get("mean_age_2_years"))
        fraction = _finite_float(params.get("binary_fraction"))
        if age1 is None or age2 is None or fraction is None:
            return None
        fraction = min(max(fraction, 0.0), 1.0)
        return fraction * age1 + (1.0 - fraction) * age2
    return _finite_float(params.get("mean_age_years"))


def diagnose_lpm_identifiability(
    fits: Sequence[JointLpmFit],
    *,
    aic_delta_threshold: float = 4.0,
    broad_age_log10_threshold: float = 0.35,
    min_tracers: int = 2,
) -> Dict[str, Any]:
    """Diagnose whether joint LPM fitting is identifiable."""

    converged = [
        fit
        for fit in fits
        if fit.converged and math.isfinite(fit.aic) and _fit_mean_age(fit) is not None
    ]
    converged = sorted(converged, key=lambda fit: fit.aic)
    if not converged:
        return {
            "best_model": None,
            "acceptable_models": [],
            "age_span_log10": float("nan"),
            "posterior_width_log10": float("nan"),
            "flags": {
                "no_converged_lpm": True,
                "low_tracer_count": True,
                "multiple_acceptable_lpm_families": False,
                "broad_acceptable_age_range": False,
                "parameter_equifinality": False,
                "not_identifiable_from_available_tracers": True,
            },
        }

    best = converged[0]
    best_aic = best.aic
    acceptable = [fit for fit in converged if fit.aic - best_aic <= aic_delta_threshold]
    ages = [_fit_mean_age(fit) for fit in acceptable]
    finite_ages = [age for age in ages if age is not None and age > 0]
    age_span = (
        max(_log10_age(age) for age in finite_ages) - min(_log10_age(age) for age in finite_ages)
        if len(finite_ages) >= 2
        else 0.0
    )
    min_tracer_count = min((fit.n_tracers for fit in acceptable), default=0)

    acceptable_rows = []
    for fit in acceptable:
        mean_age = _fit_mean_age(fit)
        acceptable_rows.append(
            {
                "model": fit.model,
                "delta_aic": float(fit.aic - best_aic),
                "mean_age_years": float(mean_age) if mean_age is not None else float("nan"),
                "objective": float(fit.objective),
                "rmse_standardized": float(fit.rmse_standardized),
                "parameters": dict(fit.parameters),
            }
        )

    flags = {
        "no_converged_lpm": False,
        "low_tracer_count": min_tracer_count < min_tracers,
        "multiple_acceptable_lpm_families": len({fit.model for fit in acceptable}) > 1,
        "broad_acceptable_age_range": age_span >= broad_age_log10_threshold,
        "parameter_equifinality": len(acceptable) > 1 and age_span >= 0.15,
        "not_identifiable_from_available_tracers": False,
    }
    flags["not_identifiable_from_available_tracers"] = any(
        flags[key]
        for key in (
            "low_tracer_count",
            "multiple_acceptable_lpm_families",
            "broad_acceptable_age_range",
            "parameter_equifinality",
        )
    )

    return {
        "best_model": best.model,
        "best_mean_age_years": float(_fit_mean_age(best) or float("nan")),
        "acceptable_models": acceptable_rows,
        "age_span_log10": float(age_span),
        "posterior_width_log10": float(age_span),
        "flags": flags,
        "interpretation": _interpret_identifiability(flags),
    }


def _interpret_identifiability(flags: Mapping[str, bool]) -> str:
    if flags.get("not_identifiable_from_available_tracers"):
        return (
            "Joint LPM age is weakly identifiable; report model ambiguity, "
            "acceptable age range, and tracer-removal sensitivity instead of a "
            "single definitive age."
        )
    return "Joint LPM fit is identifiable under the configured AIC and age-span thresholds."


TRACER_OBSERVATION_KEYS: Dict[str, List[str]] = {
    "3H": ["tritium_TU", "tritium_sigma_TU"],
    "3H/3HE": ["he3_trit_TU", "he3_trit_sigma_TU"],
    "14C": ["c14_pmc", "c14_sigma_pmc"],
    "SF6": ["sf6_pptv", "sf6_sigma_pptv"],
    "CFC11": ["cfc11_pptv", "cfc11_sigma_pptv"],
    "CFC12": ["cfc12_pptv", "cfc12_sigma_pptv"],
    "CFC113": ["cfc113_pptv", "cfc113_sigma_pptv"],
    "4HE": ["he4_ccpg", "he4_sigma_ccpg"],
}


def tracer_removal_sensitivity(
    observations: Mapping[str, Any],
    *,
    sample_year: float,
    models: Optional[Sequence[str]] = None,
    max_age_years: Optional[float] = None,
    age_steps: int = 60,
    use_helium4: bool = False,
) -> Dict[str, Any]:
    """Refit joint LPMs after removing each tracer family."""

    full_fits = fit_lpm_models(
        observations,
        sample_year=sample_year,
        models=models,
        max_age_years=max_age_years,
        age_steps=age_steps,
        use_helium4=use_helium4,
    )
    full_diag = diagnose_lpm_identifiability(full_fits)
    full_age = _finite_float(full_diag.get("best_mean_age_years"))
    rows: List[Dict[str, Any]] = []

    for tracer_key, keys in TRACER_OBSERVATION_KEYS.items():
        if not any(key in observations for key in keys):
            continue
        reduced = dict(observations)
        for key in keys:
            reduced.pop(key, None)
        fits = fit_lpm_models(
            reduced,
            sample_year=sample_year,
            models=models,
            max_age_years=max_age_years,
            age_steps=age_steps,
            use_helium4=use_helium4,
        )
        diag = diagnose_lpm_identifiability(fits)
        reduced_age = _finite_float(diag.get("best_mean_age_years"))
        delta_log = (
            _log10_age(reduced_age) - _log10_age(full_age)
            if reduced_age is not None and full_age is not None
            else float("nan")
        )
        rows.append(
            {
                "removed_tracer": tracer_key,
                "best_model": diag.get("best_model"),
                "best_mean_age_years": reduced_age if reduced_age is not None else float("nan"),
                "delta_log10_age": float(delta_log),
                "identifiability_flags": diag.get("flags", {}),
            }
        )

    return {
        "full_fit": full_diag,
        "removal_sensitivity": rows,
        "high_sensitivity": any(
            math.isfinite(row["delta_log10_age"]) and abs(row["delta_log10_age"]) >= 0.25
            for row in rows
        ),
    }
