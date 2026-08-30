"""Joint multi-tracer lumped-parameter-model fitting.

This module fits one age-distribution model against all available tracers in a
sample. It is intentionally dependency-light: the optimizer is a deterministic
grid search so it can run in validation scripts without requiring SciPy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..log import get_logger
from .input_history import InputHistory, build_default_tritium_input

logger = get_logger(__name__)
from .multi_tracer import build_atmospheric_tracer_input
from .nuclides import CARBON14, TRITIUM, ARGON39, KRYPTON85
from .tracer_inputs import (
    load_tracer_histories,
    normalize_tracer_key,
    standardize_gas_observations,
)


SUPPORTED_LPM_MODELS = (
    "PFM",
    "EM",
    "DM",
    "GA",
    "EPM",
    "PEM",
    "EMM",
    "BMM-DM-DM",
    "BMM-PEM-PEM",
)


@dataclass(frozen=True)
class TracerFitObservation:
    tracer: str
    value: float
    sigma: float
    units: str
    note: str = ""
    weight: float = 1.0
    likelihood: str = "gaussian"


@dataclass(frozen=True)
class TracerFitResidual:
    tracer: str
    observed: float
    predicted: float
    sigma: float
    standardized_residual: float


@dataclass(frozen=True)
class JointLpmFit:
    model: str
    parameters: Dict[str, float]
    objective: float
    rmse_standardized: float
    aic: float
    bic: float
    aicc: float
    effective_n_params: int
    n_tracers: int
    residuals: List[TracerFitResidual]
    converged: bool
    note: str = ""
    refinement_attempted: bool = False
    refinement_success: bool = False
    refinement_message: str = ""
    age_profile_log10_span: float = float("nan")
    age_profile_n_near_optimal: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model": self.model,
            "parameters": dict(self.parameters),
            "objective": float(self.objective),
            "rmse_standardized": float(self.rmse_standardized),
            "aic": float(self.aic),
            "bic": float(self.bic),
            "aicc": float(self.aicc),
            "effective_n_params": int(self.effective_n_params),
            "n_tracers": int(self.n_tracers),
            "converged": bool(self.converged),
            "note": self.note,
            "refinement_attempted": bool(self.refinement_attempted),
            "refinement_success": bool(self.refinement_success),
            "refinement_message": self.refinement_message,
            "age_profile_log10_span": float(self.age_profile_log10_span),
            "age_profile_n_near_optimal": int(self.age_profile_n_near_optimal),
            "residuals": [r.__dict__ for r in self.residuals],
        }


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _sigma(
    observations: Mapping[str, Any],
    key: str,
    value: float,
    *,
    relative: float,
    floor: float,
) -> float:
    configured = _finite_float(observations.get(key))
    if configured is not None and configured > 0:
        return float(configured)
    return float(max(relative * max(abs(value), floor), floor))


def build_lpm_tracer_observations(
    observations: Mapping[str, Any],
    *,
    use_helium4: bool = False,
) -> List[TracerFitObservation]:
    """Normalize Hydrosheaf/E1 tracer fields into weighted fit observations."""
    observations = standardize_gas_observations(observations)
    out: List[TracerFitObservation] = []

    tritium = _finite_float(observations.get("tritium_TU"))
    if tritium is not None and tritium >= 0:
        out.append(
            TracerFitObservation(
                "3H",
                tritium,
                _sigma(observations, "tritium_sigma_TU", tritium, relative=0.15, floor=0.1),
                "TU",
            )
        )

    he3 = _finite_float(observations.get("he3_trit_TU"))
    if he3 is not None and he3 >= 0:
        out.append(
            TracerFitObservation(
                "3H/3He",
                he3,
                _sigma(observations, "he3_trit_sigma_TU", he3, relative=0.20, floor=0.1),
                "TU",
                note="tritiogenic helium in TU-equivalent units",
            )
        )

    c14 = _finite_float(observations.get("c14_pmc"))
    if c14 is not None and 0 < c14 <= 130:
        out.append(
            TracerFitObservation(
                "14C",
                c14,
                _sigma(observations, "c14_sigma_pmc", c14, relative=0.05, floor=1.0),
                "pmc",
            )
        )

    ar39 = _finite_float(observations.get("ar39_pmc"))
    if ar39 is not None and 0 < ar39 <= 130:
        out.append(
            TracerFitObservation(
                "39Ar",
                ar39,
                _sigma(observations, "ar39_sigma_pmc", ar39, relative=0.10, floor=5.0),
                "pmc",
                note="cosmogenic noble gas isotope for intermediate ages",
            )
        )

    kr85 = _finite_float(observations.get("kr85_pptv"))
    if kr85 is not None and kr85 >= 0:
        out.append(
            TracerFitObservation(
                "85Kr",
                kr85,
                _sigma(observations, "kr85_sigma_pptv", kr85, relative=0.10, floor=0.05),
                "dpm/cc",
                note="expects corrected atmospheric equivalent",
            )
        )

    gas_specs = [
        ("SF6", "sf6_pptv", "sf6_sigma_pptv"),
        ("CFC11", "cfc11_pptv", "cfc11_sigma_pptv"),
        ("CFC12", "cfc12_pptv", "cfc12_sigma_pptv"),
        ("CFC113", "cfc113_pptv", "cfc113_sigma_pptv"),
    ]
    
    from .multi_tracer import calculate_tracer_reliability_weights, historical_max_concentration
    sample_year = _finite_float(observations.get("sample_year")) or 2026.0
    weighted_obs, _ = calculate_tracer_reliability_weights(observations, sample_year)

    for tracer, value_key, sigma_key in gas_specs:
        value = _finite_float(observations.get(value_key))
        if value is None or value < 0:
            continue
        weight = _finite_float(weighted_obs.get(f"{tracer.lower()}_weight"))
        if weight is None or weight < 0:
            weight = 1.0

        # Phase 4: mark suspect gases with likelihoods rather than treating
        # them as ordinary Gaussian observations.
        max_hist = historical_max_concentration(tracer, sample_year)
        likelihood = str(weighted_obs.get(f"{tracer.lower()}_likelihood") or "gaussian")
        if likelihood not in {"gaussian", "upper_censored", "lower_censored", "contaminated_mixture"}:
            likelihood = "gaussian"
        if value > max_hist * 1.02:
            likelihood = "upper_censored"

        out.append(
            TracerFitObservation(
                tracer,
                value,
                _sigma(observations, sigma_key, value, relative=0.10, floor=0.05),
                "pptv",
                note="expects atmospheric-equivalent corrected gas input",
                weight=weight,
                likelihood=likelihood,
            )
        )

    if use_helium4:
        he4 = _finite_float(observations.get("he4_ccpg"))
        if he4 is not None and he4 > 0:
            out.append(
                TracerFitObservation(
                    "4He",
                    he4,
                    _sigma(observations, "he4_sigma_ccpg", he4, relative=0.50, floor=1.0e-9),
                    "cc/g",
                    note="screening accumulation tracer; requires site calibration",
                )
            )

    return out


def _age_grid(max_age_years: float, steps: int) -> np.ndarray:
    max_age = max(float(max_age_years), 1.0)
    young = np.linspace(0.0, min(10.0, max_age), max(6, min(steps // 5, 30)))
    older = np.geomspace(0.25, max_age, max(steps, 20))
    linear = np.linspace(0.0, max_age, max(steps, 20))
    annual = np.arange(0.0, max_age + 1.0, 1.0) if max_age <= 500.0 else np.array([])
    return np.unique(np.round(np.concatenate([young, older, linear, annual]), 6))


def _integration_grid(mean_age: float, max_age_years: float) -> Tuple[np.ndarray, np.ndarray]:
    """Return an age grid and trapezoidal quadrature weights.

    A single uniform step tied to the oldest supported tracer made a 50,000-year
    fit use 10-year cells everywhere.  That quantized all young-water
    predictions, including piston-flow ages below 10 years.  This grid preserves
    sub-annual resolution through the modern-tracer window while using coarser
    cells only in the old-water tail.
    """
    requested_upper = max(float(max_age_years), 1.0)
    upper = min(requested_upper, max(100.0, 8.0 * max(float(mean_age), 1.0)))
    young = np.arange(0.0, min(100.0, upper) + 0.125, 0.25)
    intermediate = (
        np.arange(100.0, min(2000.0, upper) + 0.5, 1.0)
        if upper > 100.0
        else np.array([], dtype=float)
    )
    old = (
        np.arange(2000.0, upper + 5.0, 10.0)
        if upper > 2000.0
        else np.array([], dtype=float)
    )
    landmarks = [upper]
    for cutoff in (70.0, 11700.0):
        if cutoff <= upper:
            landmarks.append(cutoff)
    taus = np.unique(np.concatenate([young, intermediate, old, np.asarray(landmarks, dtype=float)]))
    taus = taus[(taus >= 0.0) & (taus <= upper)]

    weights = np.empty_like(taus)
    weights[0] = 0.5 * (taus[1] - taus[0])
    weights[-1] = 0.5 * (taus[-1] - taus[-2])
    weights[1:-1] = 0.5 * (taus[2:] - taus[:-2])
    return taus, weights


def _normalize_pdf(g: np.ndarray, integration_weights: np.ndarray | float) -> np.ndarray:
    g = np.nan_to_num(g, nan=0.0, posinf=0.0, neginf=0.0)
    g[g < 0] = 0.0
    area = float(np.sum(g * integration_weights))
    if area <= 0:
        return g
    return g / area


def _component_pdf(model: str, mean_age: float, taus: np.ndarray, params: Mapping[str, float]) -> np.ndarray:
    model_key = model.upper()
    mean = max(float(mean_age), 1.0e-6)
    g = np.zeros_like(taus, dtype=float)

    if model_key == "PFM":
        idx = int(np.argmin(np.abs(taus - mean)))
        g[idx] = 1.0
        return g

    if model_key == "EM":
        g = (1.0 / mean) * np.exp(-taus / mean)
        return g

    if model_key == "DM":
        dispersion = max(float(params.get("dispersion", 0.1)), 1.0e-4)
        ts = np.maximum(taus, 1.0e-9)
        term1 = 1.0 / (mean * np.sqrt(4.0 * math.pi * dispersion * (ts / mean)))
        term2 = np.exp(-((1.0 - ts / mean) ** 2) / (4.0 * dispersion * (ts / mean)))
        g = term1 * term2
        g[taus == 0] = 0.0
        return g

    if model_key == "GA":
        # Gamma distribution: shape (a) and scale (theta = mean/a)
        # Using scipy.stats if possible, else manual.
        shape = max(float(params.get("shape", 1.0)), 0.1)
        scale = mean / shape
        from scipy.stats import gamma
        return gamma.pdf(taus, shape, scale=scale)

    if model_key == "EPM":
        piston_fraction = min(max(float(params.get("piston_fraction", 0.5)), 0.0), 0.95)
        piston_delay = piston_fraction * mean
        em_mean = max(mean - piston_delay, 1.0e-6)
        active = taus >= piston_delay
        g[active] = (1.0 / em_mean) * np.exp(-(taus[active] - piston_delay) / em_mean)
        return g

    if model_key == "PEM":
        capture_fraction = min(max(float(params.get("capture_fraction", 0.7)), 0.05), 1.0)
        em_mean = max(mean, 1.0e-6)
        cutoff = -em_mean * math.log(max(1.0 - capture_fraction, 1.0e-9))
        active = taus <= cutoff
        g[active] = (1.0 / em_mean) * np.exp(-taus[active] / em_mean)
        return g

    if model_key == "EMM":
        young_fraction = min(max(float(params.get("young_fraction", 0.7)), 0.0), 1.0)
        old_age_ratio = max(float(params.get("old_age_ratio", 5.0)), 1.1)
        old_mean = mean * old_age_ratio
        young_mean = max((mean - (1.0 - young_fraction) * old_mean) / max(young_fraction, 1.0e-6), 0.1)
        g = young_fraction * (1.0 / young_mean) * np.exp(-taus / young_mean)
        g += (1.0 - young_fraction) * (1.0 / old_mean) * np.exp(-taus / old_mean)
        return g

    raise ValueError(f"Unsupported LPM component model: {model}")


def transit_time_pdf(
    model: str,
    parameters: Mapping[str, float],
    taus: np.ndarray,
    integration_weights: np.ndarray | float,
) -> np.ndarray:
    """Return a normalized transit-time distribution for a model family."""
    model_key = model.upper()
    if model_key.startswith("BMM-"):
        parts = model_key.split("-")
        if len(parts) != 3:
            raise ValueError("Binary-mixture model names must look like BMM-DM-DM.")
        fraction = min(max(float(parameters.get("binary_fraction", 0.5)), 0.0), 1.0)
        g1 = _normalize_pdf(
            _component_pdf(parts[1], float(parameters["mean_age_1_years"]), taus, parameters),
            integration_weights,
        )
        params2 = {
            key[:-2] if key.endswith("_2") else key: value
            for key, value in parameters.items()
        }
        g2 = _normalize_pdf(
            _component_pdf(parts[2], float(parameters["mean_age_2_years"]), taus, params2),
            integration_weights,
        )
        return _normalize_pdf(fraction * g1 + (1.0 - fraction) * g2, integration_weights)

    g = _component_pdf(model_key, float(parameters["mean_age_years"]), taus, parameters)
    return _normalize_pdf(g, integration_weights)


def age_fraction_predictions(
    model: str,
    parameters: Mapping[str, float],
    max_age_years: float = 50000.0,
) -> dict[str, float]:
    """Return predicted Anthropocene/Holocene/Pleistocene fractions for a parameter set."""
    model_key = model.upper()
    if model_key.startswith("BMM-"):
        age_scale = max(
            float(parameters.get("mean_age_1_years", 1.0)),
            float(parameters.get("mean_age_2_years", 1.0)),
        )
    else:
        age_scale = float(parameters.get("mean_age_years", 1.0))
    taus, integration_weights = _integration_grid(age_scale, max_age_years)
    pdf = transit_time_pdf(model_key, parameters, taus, integration_weights)

    # Simple fixed-age cutoffs for first implementation
    anthropocene_mask = taus <= 70.0
    holocene_mask = (taus > 70.0) & (taus <= 11700.0)
    pleistocene_mask = taus > 11700.0

    frac_anthro = float(np.sum(pdf[anthropocene_mask] * integration_weights[anthropocene_mask])) if np.any(anthropocene_mask) else 0.0
    frac_holocene = float(np.sum(pdf[holocene_mask] * integration_weights[holocene_mask])) if np.any(holocene_mask) else 0.0
    frac_pleistocene = float(np.sum(pdf[pleistocene_mask] * integration_weights[pleistocene_mask])) if np.any(pleistocene_mask) else 0.0
    total = frac_anthro + frac_holocene + frac_pleistocene
    if total > 0:
        frac_anthro /= total
        frac_holocene /= total
        frac_pleistocene /= total
    return {
        "anthropocene": frac_anthro,
        "holocene": frac_holocene,
        "pleistocene": frac_pleistocene,
    }


def _predict_from_history(
    sample_year: float,
    history: InputHistory,
    pdf: np.ndarray,
    taus: np.ndarray,
    integration_weights: np.ndarray | float,
    decay_lambda_per_year: float = 0.0,
    *,
    daughter: bool = False,
    left_value: Optional[float] = None,
) -> float:
    recharge_years = float(sample_year) - taus
    left = float(history.values[0]) if left_value is None else float(left_value)
    input_values = np.interp(recharge_years, history.years, history.values, left=left)
    decay = np.exp(-decay_lambda_per_year * taus)
    if daughter:
        values = input_values * (1.0 - decay)
    else:
        values = input_values * decay
    return float(np.sum(values * pdf * integration_weights))


def _predict_c14(pdf: np.ndarray, taus: np.ndarray, integration_weights: np.ndarray | float, initial_pmc: float) -> float:
    lambda_y = math.log(2.0) / CARBON14.half_life_years
    return float(np.sum(initial_pmc * np.exp(-lambda_y * taus) * pdf * integration_weights))


def _predict_ar39(pdf: np.ndarray, taus: np.ndarray, integration_weights: np.ndarray | float, initial_pmc: float) -> float:
    lambda_y = math.log(2.0) / ARGON39.half_life_years
    return float(np.sum(initial_pmc * np.exp(-lambda_y * taus) * pdf * integration_weights))


def tracer_response_kernel(
    tracer: str,
    age_grid_years: Sequence[float],
    sample_year: float,
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    initial_c14_pmc: float = 100.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
) -> np.ndarray:
    """Return the predicted tracer response for unit mass in each age bin.

    The returned vector is the linear forward operator used by set-valued TTD
    inference: for an age-bin probability-mass vector ``w``, the predicted
    concentration is ``kernel @ w``.  It deliberately shares the tracer input
    histories and decay assumptions used by :func:`predict_lpm_tracers`.
    """
    ages = np.asarray(age_grid_years, dtype=float)
    if ages.ndim != 1 or ages.size == 0:
        raise ValueError("age_grid_years must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(ages)) or np.any(ages < 0.0):
        raise ValueError("age_grid_years must be finite and non-negative.")
    if not math.isfinite(float(sample_year)):
        raise ValueError("sample_year must be finite.")

    key = normalize_tracer_key(tracer)
    hists = {
        normalize_tracer_key(history_key): history
        for history_key, history in dict(histories or {}).items()
    }
    recharge_years = float(sample_year) - ages

    if key in {"3H", "3H/3He"}:
        history = hists.get("3H") or build_default_tritium_input()
        inputs = np.interp(
            recharge_years,
            history.years,
            history.values,
            left=float(history.values[0]),
        )
        decay = np.exp(-(math.log(2.0) / TRITIUM.half_life_years) * ages)
        kernel = inputs * (decay if key == "3H" else 1.0 - decay)
    elif key == "14C":
        kernel = float(initial_c14_pmc) * np.exp(
            -(math.log(2.0) / CARBON14.half_life_years) * ages
        )
    elif key == "39Ar":
        kernel = 100.0 * np.exp(
            -(math.log(2.0) / ARGON39.half_life_years) * ages
        )
    elif key in {"SF6", "CFC11", "CFC12", "CFC113", "85Kr"}:
        history = hists.get(key) or build_atmospheric_tracer_input(key)
        inputs = np.interp(
            recharge_years,
            history.years,
            history.values,
            left=float(history.values[0]),
        )
        if key == "85Kr":
            inputs = inputs * np.exp(
                -(math.log(2.0) / KRYPTON85.half_life_years) * ages
            )
        kernel = inputs
    elif key == "4He":
        kernel = (
            float(helium4_background_ccpg)
            + float(helium4_accumulation_rate_ccpg_per_year) * ages
        )
    else:
        raise ValueError(f"Unsupported tracer response kernel: {tracer!r}.")

    scales = {
        normalize_tracer_key(scale_key): _finite_float(scale_value)
        for scale_key, scale_value in dict(prediction_scale_factors or {}).items()
    }
    scale = scales.get(key)
    if scale is not None and scale > 0.0:
        kernel = kernel * scale
    return np.asarray(kernel, dtype=float)


def _mean_from_distribution(pdf: np.ndarray, taus: np.ndarray, integration_weights: np.ndarray | float) -> float:
    return float(np.sum(taus * pdf * integration_weights))


def predict_lpm_tracers(
    model: str,
    parameters: Mapping[str, float],
    sample_year: float,
    tracers: Iterable[str],
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    initial_c14_pmc: float = 100.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    max_age_years: float = 50000.0,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
) -> Dict[str, float]:
    """Predict tracer values from one LPM parameter set.

    ``prediction_scale_factors`` reproduces source-model recharge/input scaling
    when a benchmark publishes it explicitly.  Scaling model predictions (not
    measurements) preserves residuals in the observations' reported units.
    """
    model_key = model.upper()
    if model_key.startswith("BMM-"):
        age_scale = max(
            float(parameters.get("mean_age_1_years", 1.0)),
            float(parameters.get("mean_age_2_years", 1.0)),
        )
    else:
        age_scale = float(parameters.get("mean_age_years", 1.0))
    taus, integration_weights = _integration_grid(age_scale, max_age_years)
    pdf = transit_time_pdf(model_key, parameters, taus, integration_weights)
    hists = {normalize_tracer_key(key): value for key, value in dict(histories or {}).items()}
    out: Dict[str, float] = {}
    lambda_tritium = math.log(2.0) / TRITIUM.half_life_years

    requested = {tracer.upper().replace("-", "") for tracer in tracers}
    if "3H" in requested or "3H/3HE" in requested:
        tritium_hist = hists.get("3H") or hists.get("TRITIUM") or build_default_tritium_input()
        if "3H" in requested:
            out["3H"] = _predict_from_history(
                sample_year,
                tritium_hist,
                pdf,
                taus,
                integration_weights,
                decay_lambda_per_year=lambda_tritium,
                left_value=float(tritium_hist.values[0]),
            )
        if "3H/3HE" in requested:
            out["3H/3He"] = _predict_from_history(
                sample_year,
                tritium_hist,
                pdf,
                taus,
                integration_weights,
                decay_lambda_per_year=lambda_tritium,
                daughter=True,
                left_value=float(tritium_hist.values[0]),
            )

    if "14C" in requested:
        out["14C"] = _predict_c14(pdf, taus, integration_weights, initial_c14_pmc)

    if "39AR" in requested:
        # Standard initial Ar-39 activity is 100 pmc
        out["39Ar"] = _predict_ar39(pdf, taus, integration_weights, 100.0)

    for tracer in ("SF6", "CFC11", "CFC12", "CFC113", "85KR"):
        if tracer in requested:
            history = hists.get(tracer) or build_atmospheric_tracer_input(tracer)
            
            decay_lambda = 0.0
            if tracer == "85KR":
                decay_lambda = math.log(2.0) / KRYPTON85.half_life_years
                
            out[tracer] = _predict_from_history(
                sample_year, 
                history, 
                pdf, 
                taus, 
                integration_weights,
                decay_lambda_per_year=decay_lambda,
                left_value=float(history.values[0]),
            )
            if tracer == "85KR":
                # Rename back to 85Kr for consistency with observations
                out["85Kr"] = out.pop("85KR")

    if "4HE" in requested:
        mean_age = _mean_from_distribution(pdf, taus, integration_weights)
        out["4He"] = helium4_background_ccpg + helium4_accumulation_rate_ccpg_per_year * mean_age

    output_keys = {
        "3H": "3H",
        "3H/3HE": "3H/3He",
        "14C": "14C",
        "39AR": "39Ar",
        "SF6": "SF6",
        "CFC11": "CFC11",
        "CFC12": "CFC12",
        "CFC113": "CFC113",
        "85KR": "85Kr",
    }
    for tracer, raw_scale in dict(prediction_scale_factors or {}).items():
        normalized = str(tracer).upper().replace("-", "")
        output_key = output_keys.get(normalized)
        scale = _finite_float(raw_scale)
        if output_key in out and scale is not None and scale > 0:
            out[output_key] *= scale

    return out


def _parameter_candidates(model: str, age_grid: np.ndarray) -> Tuple[List[Dict[str, float]], int]:
    model_key = model.upper()
    candidates: List[Dict[str, float]] = []

    if model_key in {"PFM", "EM"}:
        return ([{"mean_age_years": float(age)} for age in age_grid], 1)

    if model_key == "DM":
        for age in age_grid:
            for dispersion in (0.01, 0.03, 0.1, 0.3, 1.0):
                candidates.append({"mean_age_years": float(age), "dispersion": dispersion})
        return candidates, 2

    if model_key == "GA":
        for age in age_grid:
            for shape in (0.5, 1.0, 2.0, 5.0):
                candidates.append({"mean_age_years": float(age), "shape": shape})
        return candidates, 2

    if model_key == "EPM":
        for age in age_grid:
            for piston_fraction in (0.1, 0.3, 0.5, 0.7, 0.9):
                candidates.append({"mean_age_years": float(age), "piston_fraction": piston_fraction})
        return candidates, 2

    if model_key == "PEM":
        for age in age_grid:
            for capture_fraction in (0.1, 0.25, 0.5, 0.75, 0.95):
                candidates.append({"mean_age_years": float(age), "capture_fraction": capture_fraction})
        return candidates, 2

    if model_key == "EMM":
        for age in age_grid:
            for young_fraction in (0.25, 0.5, 0.75, 0.9):
                for old_age_ratio in (2.0, 5.0, 10.0):
                    candidates.append(
                        {
                            "mean_age_years": float(age),
                            "young_fraction": young_fraction,
                            "old_age_ratio": old_age_ratio,
                        }
                    )
        return candidates, 3

    if model_key.startswith("BMM-"):
        base_ages = age_grid[:: max(1, len(age_grid) // 45)]
        for age1 in base_ages:
            for ratio in (2.0, 5.0, 10.0, 25.0):
                age2 = min(float(age1) * ratio + 1.0, float(age_grid[-1]))
                if age2 <= age1:
                    continue
                for fraction in (0.2, 0.5, 0.8):
                    row = {
                        "binary_fraction": fraction,
                        "mean_age_1_years": float(age1),
                        "mean_age_2_years": float(age2),
                    }
                    if "DM" in model_key:
                        row["dispersion"] = 0.1
                        row["dispersion_2"] = 0.1
                    if "PEM" in model_key:
                        row["capture_fraction"] = 0.75
                        row["capture_fraction_2"] = 0.75
                    candidates.append(row)
        return candidates, 5

    raise ValueError(
        "Unsupported LPM model. Expected one of: " + ", ".join(SUPPORTED_LPM_MODELS)
    )


def _pack_parameters(
    model: str,
    params: Mapping[str, float],
    max_age_years: float = 50000.0,
    free_parameters: Optional[Sequence[str]] = None,
) -> tuple[np.ndarray, list[str], list[tuple[float, float]]]:
    """Pack model parameters into a flat array for optimization with transforms."""
    m = model.upper()
    is_bmm = m.startswith("BMM-")
    keys: list[str] = []
    values: list[float] = []
    bounds: list[tuple[float, float]] = []

    if is_bmm:
        age1 = float(params.get("mean_age_1_years", 1.0))
        age2 = float(params.get("mean_age_2_years", 10.0))
        frac = float(params.get("binary_fraction", 0.5))

        keys.extend(["mean_age_1_years", "mean_age_2_years", "binary_fraction"])
        values.extend([math.log10(max(age1, 0.01)), math.log10(max(age2, 0.01)), math.log(max(frac, 1e-6) / max(1 - frac, 1e-6))])
        log_min = math.log10(0.01)
        log_max = math.log10(max(max_age_years, 1.0))
        bounds.extend([(log_min, log_max), (log_min, log_max), (-6.0, 6.0)])

        if "dispersion" in params:
            keys.append("dispersion")
            values.append(float(params["dispersion"]))
            bounds.append((1e-4, 2.0))
        if "dispersion_2" in params:
            keys.append("dispersion_2")
            values.append(float(params["dispersion_2"]))
            bounds.append((1e-4, 2.0))
        if "capture_fraction" in params:
            keys.append("capture_fraction")
            values.append(float(params["capture_fraction"]))
            bounds.append((0.05, 1.0))
        if "capture_fraction_2" in params:
            keys.append("capture_fraction_2")
            values.append(float(params["capture_fraction_2"]))
            bounds.append((0.05, 1.0))
    else:
        if "mean_age_years" in params:
            keys.append("mean_age_years")
            values.append(float(params["mean_age_years"]))
            bounds.append((0.01, max(max_age_years, 1.0)))

        for key in ("dispersion", "shape", "piston_fraction", "capture_fraction"):
            if key in params:
                keys.append(key)
                values.append(float(params[key]))
                if key == "dispersion":
                    bounds.append((1e-4, 2.0))
                elif key == "shape":
                    bounds.append((0.1, 10.0))
                elif key == "piston_fraction":
                    bounds.append((0.0, 0.95))
                elif key == "capture_fraction":
                    bounds.append((0.05, 1.0))

    if free_parameters is not None:
        allowed = set(free_parameters)
        selected = [
            (value, key, bound)
            for value, key, bound in zip(values, keys, bounds)
            if key in allowed
        ]
        values = [item[0] for item in selected]
        keys = [item[1] for item in selected]
        bounds = [item[2] for item in selected]

    return np.array(values, dtype=float), keys, bounds


def _unpack_parameters(
    model: str,
    x: np.ndarray,
    keys: list[str],
    bounds: list[tuple[float, float]],
) -> dict[str, float]:
    """Unpack array into parameter dict with model-specific constraints."""
    m = model.upper()
    is_bmm = m.startswith("BMM-")
    p: dict[str, float] = {}
    for key, v, (lb, ub) in zip(keys, x, bounds):
        clipped = float(min(max(v, lb), ub))
        if key == "mean_age_1_years":
            p[key] = 10.0 ** clipped
        elif key == "mean_age_2_years":
            p[key] = 10.0 ** clipped
        elif key == "binary_fraction":
            p[key] = 1.0 / (1.0 + math.exp(-clipped))
        else:
            p[key] = clipped

    if is_bmm:
        age1 = p.get("mean_age_1_years")
        age2 = p.get("mean_age_2_years")
        if age1 is not None and age2 is not None and age2 <= age1:
            p["mean_age_2_years"] = age1 + 0.1
        if "binary_fraction" in p:
            p["binary_fraction"] = min(max(p["binary_fraction"], 0.001), 0.999)
    else:
        if "mean_age_years" in p:
            p["mean_age_years"] = max(0.01, p["mean_age_years"])
        if "dispersion" in p:
            p["dispersion"] = min(max(p["dispersion"], 1e-4), 2.0)
        if "shape" in p:
            p["shape"] = min(max(p["shape"], 0.1), 10.0)
        if "piston_fraction" in p:
            p["piston_fraction"] = min(max(p["piston_fraction"], 0.0), 0.95)
        if "capture_fraction" in p:
            p["capture_fraction"] = min(max(p["capture_fraction"], 0.05), 1.0)

    return p


def _aicc(aic: float, n: int, k: int) -> float:
    """Corrected AIC for small sample sizes."""
    if n <= k + 1:
        return float("inf")
    return aic + (2 * k * (k + 1)) / (n - k - 1)


def compute_gated_bma_age(
    fits: list[JointLpmFit],
    *,
    max_delta_aicc: float = 4.0,
    max_log_age_span: float = 0.7,
) -> dict[str, Any]:
    """Compute a gated Bayesian model average age.

    Rules:
    - Ignore non-converged fits.
    - Use AICc when at least two models have finite AICc; otherwise use AIC.
    - Use only fits within the configured information-criterion delta.
    - If contributing model ages span more than max_log_age_span, use top model.
    - If fewer than two models pass gates, use top model.
    - Record whether BMA was used or skipped.
    """
    converged = [f for f in fits if f.converged]
    finite_aicc = [f for f in converged if math.isfinite(f.aicc)]
    if len(finite_aicc) >= 2:
        criterion_name = "aicc"
        eligible = finite_aicc
    else:
        criterion_name = "aic_fallback"
        eligible = [f for f in converged if math.isfinite(f.aic)]

    if not eligible:
        return {
            "age_years": float("nan"),
            "bma_used": False,
            "bma_skip_reason": "no_converged_fits",
            "bma_n_models": 0,
            "bma_log_age_span": float("nan"),
            "top_model_age_years": float("nan"),
            "bma_information_criterion": "none",
        }

    def criterion_value(fit: JointLpmFit) -> float:
        return float(fit.aicc if criterion_name == "aicc" else fit.aic)

    eligible = sorted(eligible, key=criterion_value)
    best_criterion = criterion_value(eligible[0])
    within_delta = [
        fit for fit in eligible
        if criterion_value(fit) - best_criterion <= max_delta_aicc
    ]

    if len(within_delta) < 2:
        top = eligible[0]
        top_age = _extract_scalar_age(top.parameters)
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "insufficient_models_within_delta_aicc",
            "bma_n_models": len(within_delta),
            "bma_log_age_span": float("nan"),
            "top_model_age_years": top_age,
            "bma_information_criterion": criterion_name,
        }

    age_fits = [
        (fit, _extract_scalar_age(fit.parameters))
        for fit in within_delta
    ]
    age_fits = [(fit, age) for fit, age in age_fits if math.isfinite(age) and age > 0]
    if len(age_fits) < 2:
        top_age = age_fits[0][1] if age_fits else float("nan")
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "insufficient_finite_ages",
            "bma_n_models": len(age_fits),
            "bma_log_age_span": float("nan"),
            "top_model_age_years": top_age,
            "bma_information_criterion": criterion_name,
        }

    ages = [age for _, age in age_fits]
    log_ages = [math.log10(a) for a in ages]
    log_span = max(log_ages) - min(log_ages)
    if log_span > max_log_age_span:
        top_age = ages[0]
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "log_age_span_exceeded",
            "bma_n_models": len(age_fits),
            "bma_log_age_span": log_span,
            "top_model_age_years": top_age,
            "bma_information_criterion": criterion_name,
        }

    deltas = [criterion_value(fit) - best_criterion for fit, _ in age_fits]
    weights = [math.exp(-0.5 * delta) for delta in deltas]
    total_w = sum(weights)
    if total_w <= 0 or not math.isfinite(total_w):
        top_age = ages[0]
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "weight_computation_failed",
            "bma_n_models": len(age_fits),
            "bma_log_age_span": log_span,
            "top_model_age_years": top_age,
            "bma_information_criterion": criterion_name,
        }

    bma_age = sum(w * a for w, a in zip(weights, ages)) / total_w
    return {
        "age_years": bma_age,
        "bma_used": True,
        "bma_skip_reason": "",
        "bma_n_models": len(age_fits),
        "bma_log_age_span": log_span,
        "top_model_age_years": ages[0],
        "bma_information_criterion": criterion_name,
    }


def _extract_scalar_age(parameters: Dict[str, float]) -> float:
    """Extract a scalar mean age from parameters."""
    if "mean_age_years" in parameters:
        return float(parameters["mean_age_years"])
    if "mean_age_1_years" in parameters and "mean_age_2_years" in parameters:
        fraction = float(parameters.get("binary_fraction", 0.5))
        return fraction * parameters["mean_age_1_years"] + (1.0 - fraction) * parameters["mean_age_2_years"]
    return float("nan")


def _standardized_observation_loss(obs: TracerFitObservation, predicted: float) -> float:
    """Return a loss contribution for a single observation-prediction pair."""
    standardized = (obs.value - predicted) / obs.sigma
    if obs.likelihood == "upper_censored":
        # No penalty when predicted <= observed; penalty only when predicted > observed
        if predicted <= obs.value:
            return 0.0
        return float(standardized * standardized)
    if obs.likelihood == "lower_censored":
        # No penalty when predicted >= observed; penalty only when predicted < observed
        if predicted >= obs.value:
            return 0.0
        return float(standardized * standardized)
    if obs.likelihood == "contaminated_mixture":
        # Student-t-like robust loss (Cauchy-like tails)
        return float(math.log1p(standardized * standardized))
    # Default gaussian
    return float(standardized * standardized)


def _evaluate_lpm_candidate(
    model: str,
    parameters: Mapping[str, float],
    sample_year: float,
    tracer_names: list[str],
    fit_observations: list[TracerFitObservation],
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    initial_c14_pmc: float = 100.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    max_age_years: float = 50000.0,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
) -> tuple[float, List[TracerFitResidual], bool]:
    """Score one candidate with the same loss used by grid and refinement."""
    predicted = predict_lpm_tracers(
        model,
        parameters,
        sample_year,
        tracer_names,
        histories=histories,
        initial_c14_pmc=initial_c14_pmc,
        helium4_background_ccpg=helium4_background_ccpg,
        helium4_accumulation_rate_ccpg_per_year=helium4_accumulation_rate_ccpg_per_year,
        max_age_years=max_age_years,
        prediction_scale_factors=prediction_scale_factors,
    )
    objective = 0.0
    residuals: List[TracerFitResidual] = []
    for obs in fit_observations:
        pred = predicted.get(obs.tracer)
        if pred is None or not math.isfinite(pred):
            return float("inf"), [], False
        standardized = (obs.value - pred) / obs.sigma
        objective += _standardized_observation_loss(obs, pred) * obs.weight
        residuals.append(
            TracerFitResidual(
                tracer=obs.tracer,
                observed=float(obs.value),
                predicted=float(pred),
                sigma=float(obs.sigma),
                standardized_residual=float(standardized),
            )
        )

    if age_fraction_obs is not None:
        frac_pred = age_fraction_predictions(model, parameters, max_age_years=max_age_years)
        sigma_fraction = _finite_float(age_fraction_obs.get("sigma_fraction")) or 0.10
        sigma_fraction = max(sigma_fraction, 1.0e-9)
        for key in ("anthropocene", "holocene", "pleistocene"):
            obs_val = _finite_float(age_fraction_obs.get(key))
            if obs_val is not None:
                pred_val = frac_pred.get(key, 0.0)
                residual = (obs_val - pred_val) / sigma_fraction
                objective += residual * residual
    return float(objective), residuals, True


def _refine_best_candidate(
    model: str,
    best_params: Dict[str, float],
    sample_year: float,
    tracer_names: list[str],
    fit_observations: list[TracerFitObservation],
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    initial_c14_pmc: float = 100.0,
    helium4_background_ccpg: float = 4.6e-8,
    helium4_accumulation_rate_ccpg_per_year: float = 2.0e-11,
    max_age_years: float = 50000.0,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
    parameter_template: Optional[Mapping[str, float]] = None,
    free_parameters: Optional[Sequence[str]] = None,
) -> tuple[Dict[str, float], bool, str]:
    """Refine a grid-best candidate with local optimization.

    Returns (refined_params, success, message).
    """
    try:
        from scipy.optimize import minimize
    except Exception:
        return best_params, False, "scipy not available"

    x0, keys, bounds = _pack_parameters(
        model,
        best_params,
        max_age_years=max_age_years,
        free_parameters=free_parameters,
    )
    if not keys:
        return best_params, False, "no free parameters to refine"

    def objective_fn(x: np.ndarray) -> float:
        params = dict(parameter_template or {})
        params.update(_unpack_parameters(model, x, keys, bounds))
        if str(model).upper().startswith("BMM-"):
            age1 = _finite_float(params.get("mean_age_1_years"))
            age2 = _finite_float(params.get("mean_age_2_years"))
            if age1 is not None and age2 is not None and age2 <= age1:
                return 1.0e12
        objective, _, valid = _evaluate_lpm_candidate(
            model,
            params,
            sample_year,
            tracer_names,
            fit_observations,
            histories=histories,
            initial_c14_pmc=initial_c14_pmc,
            helium4_background_ccpg=helium4_background_ccpg,
            helium4_accumulation_rate_ccpg_per_year=helium4_accumulation_rate_ccpg_per_year,
            max_age_years=max_age_years,
            age_fraction_obs=age_fraction_obs,
            prediction_scale_factors=prediction_scale_factors,
        )
        return objective if valid else 1.0e12

    # Nelder-Mead first
    try:
        res_nm = minimize(objective_fn, x0, method="Nelder-Mead", options={"maxiter": 300, "xatol": 1e-4, "fatol": 1e-4})
        if res_nm.success and res_nm.fun < objective_fn(x0):
            refined = dict(parameter_template or {})
            refined.update(_unpack_parameters(model, res_nm.x, keys, bounds))
            return refined, True, "Nelder-Mead refined"
        x0 = res_nm.x
    except Exception:
        pass

    # L-BFGS-B fallback on bounded transformed parameters
    try:
        res_bfgs = minimize(objective_fn, x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": 300})
        if res_bfgs.success and res_bfgs.fun < objective_fn(x0):
            refined = dict(parameter_template or {})
            refined.update(_unpack_parameters(model, res_bfgs.x, keys, bounds))
            return refined, True, "L-BFGS-B refined"
    except Exception as exc:
        return best_params, False, f"refinement failed: {exc}"

    return best_params, False, "refinement did not improve objective"


def fit_lpm_model(
    observations: Mapping[str, Any],
    *,
    sample_year: float,
    model: str = "DM",
    histories: Optional[Mapping[str, InputHistory]] = None,
    history_paths: Optional[Iterable[str]] = None,
    max_age_years: Optional[float] = None,
    age_steps: int = 90,
    age_profile_objective_window: float = 4.0,
    use_helium4: bool = False,
    initial_c14_pmc: float = 100.0,
    refine: bool = False,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
    parameter_template: Optional[Mapping[str, float]] = None,
    free_parameters: Optional[Sequence[str]] = None,
) -> JointLpmFit:
    """Fit one LPM family jointly to all supported tracer observations."""
    fit_observations = build_lpm_tracer_observations(observations, use_helium4=use_helium4)
    resolved_histories = {normalize_tracer_key(key): value for key, value in dict(histories or {}).items()}
    if history_paths:
        resolved_histories.update(load_tracer_histories(history_paths))
    if not fit_observations:
        logger.debug(f"fit_lpm_model: no fit_observations found for {observations}")
        return JointLpmFit(
            model=model.upper(),
            parameters={},
            objective=float("nan"),
            rmse_standardized=float("nan"),
            aic=float("nan"),
            bic=float("nan"),
            aicc=float("nan"),
            effective_n_params=0,
            n_tracers=0,
            residuals=[],
            converged=False,
            note="no supported tracer observations",
            refinement_attempted=False,
            refinement_success=False,
            refinement_message="",
        )

    has_old_tracer = any(obs.tracer in {"14C", "4He"} for obs in fit_observations)
    max_age = float(max_age_years if max_age_years is not None else (50000.0 if has_old_tracer else 250.0))
    helium4_background = _finite_float(observations.get("he4_background_ccpg"))
    if helium4_background is None:
        helium4_background = 4.6e-8
    helium4_rate = _finite_float(observations.get("he4_accumulation_rate_ccpg_per_year"))
    if helium4_rate is None or helium4_rate <= 0:
        helium4_rate = 2.0e-11
    ages = _age_grid(max_age, age_steps)
    candidates, n_params = _parameter_candidates(model, ages)
    template = {
        str(key): float(value)
        for key, value in dict(parameter_template or {}).items()
        if _finite_float(value) is not None
    }
    if free_parameters is not None:
        free_set = set(free_parameters)
        n_params = len(free_set)
        configured_candidates: list[dict[str, float]] = []
        seen: set[tuple[tuple[str, float], ...]] = set()
        for candidate in candidates:
            configured = dict(candidate)
            for key, value in template.items():
                if key not in free_set:
                    configured[key] = value
            signature = tuple(sorted((key, round(float(value), 12)) for key, value in configured.items()))
            if signature not in seen:
                seen.add(signature)
                configured_candidates.append(configured)
        candidates = configured_candidates
    tracer_names = [obs.tracer for obs in fit_observations]

    best: Optional[Tuple[float, Dict[str, float], List[TracerFitResidual]]] = None
    grid_profile: list[tuple[float, float]] = []
    for params in candidates:
        objective, residuals, valid = _evaluate_lpm_candidate(
            model,
            params,
            sample_year,
            tracer_names,
            fit_observations,
            histories=resolved_histories,
            initial_c14_pmc=initial_c14_pmc,
            helium4_background_ccpg=helium4_background,
            helium4_accumulation_rate_ccpg_per_year=helium4_rate,
            max_age_years=max_age,
            age_fraction_obs=age_fraction_obs,
            prediction_scale_factors=prediction_scale_factors,
        )
        if not valid:
            continue
        scalar_age = _extract_scalar_age(params)
        if math.isfinite(scalar_age) and scalar_age > 0:
            grid_profile.append((float(objective), float(scalar_age)))
        if best is None or objective < best[0]:
            best = (float(objective), dict(params), residuals)

    if best is None:
        return JointLpmFit(
            model=model.upper(),
            parameters={},
            objective=float("nan"),
            rmse_standardized=float("nan"),
            aic=float("nan"),
            bic=float("nan"),
            aicc=float("nan"),
            effective_n_params=n_params,
            n_tracers=len(fit_observations),
            residuals=[],
            converged=False,
            note="no finite model prediction",
            refinement_attempted=False,
            refinement_success=False,
            refinement_message="",
        )

    objective, params, residuals = best
    profile_threshold = float(objective) + max(float(age_profile_objective_window), 0.0)
    near_optimal_ages = [
        age for candidate_objective, age in grid_profile
        if candidate_objective <= profile_threshold
    ]
    age_profile_n_near_optimal = len(near_optimal_ages)
    if near_optimal_ages:
        age_profile_log10_span = float(
            max(math.log10(age) for age in near_optimal_ages)
            - min(math.log10(age) for age in near_optimal_ages)
        )
    else:
        age_profile_log10_span = float("nan")

    # Phase 8: local refinement
    refinement_attempted = False
    refinement_success = False
    refinement_message = ""
    if refine:
        refinement_attempted = True
        refined_params, refinement_success, refinement_message = _refine_best_candidate(
            model,
            params,
            sample_year,
            tracer_names,
            fit_observations,
            histories=resolved_histories,
            initial_c14_pmc=initial_c14_pmc,
            helium4_background_ccpg=helium4_background,
            helium4_accumulation_rate_ccpg_per_year=helium4_rate,
            max_age_years=max_age,
            age_fraction_obs=age_fraction_obs,
            prediction_scale_factors=prediction_scale_factors,
            parameter_template=template,
            free_parameters=free_parameters,
        )
        if refinement_success:
            new_objective, new_residuals, valid = _evaluate_lpm_candidate(
                model,
                refined_params,
                sample_year,
                tracer_names,
                fit_observations,
                histories=resolved_histories,
                initial_c14_pmc=initial_c14_pmc,
                helium4_background_ccpg=helium4_background,
                helium4_accumulation_rate_ccpg_per_year=helium4_rate,
                max_age_years=max_age,
                age_fraction_obs=age_fraction_obs,
                prediction_scale_factors=prediction_scale_factors,
            )
            if valid and new_objective < objective:
                params = refined_params
                objective = new_objective
                residuals = new_residuals
            else:
                refinement_success = False
                refinement_message = "refinement candidate did not improve the unified objective"

    n = len(residuals)
    rss_per_obs = max(objective / max(n, 1), 1.0e-12)
    aic = n * math.log(rss_per_obs) + 2.0 * n_params
    bic = n * math.log(rss_per_obs) + n_params * math.log(max(n, 2))
    aicc_val = _aicc(aic, n, n_params)
    return JointLpmFit(
        model=model.upper(),
        parameters=params,
        objective=objective,
        rmse_standardized=math.sqrt(rss_per_obs),
        aic=float(aic),
        bic=float(bic),
        aicc=float(aicc_val),
        effective_n_params=n_params,
        n_tracers=n,
        residuals=residuals,
        converged=True,
        refinement_attempted=refinement_attempted,
        refinement_success=refinement_success,
        refinement_message=refinement_message,
        age_profile_log10_span=age_profile_log10_span,
        age_profile_n_near_optimal=age_profile_n_near_optimal,
    )


def fit_lpm_models(
    observations: Mapping[str, Any],
    *,
    sample_year: float,
    models: Optional[Sequence[str]] = None,
    histories: Optional[Mapping[str, InputHistory]] = None,
    history_paths: Optional[Iterable[str]] = None,
    max_age_years: Optional[float] = None,
    age_steps: int = 90,
    age_profile_objective_window: float = 4.0,
    use_helium4: bool = False,
    initial_c14_pmc: float = 100.0,
    refine: bool = False,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
    prediction_scale_factors: Optional[Mapping[str, float]] = None,
    parameter_template: Optional[Mapping[str, float]] = None,
    free_parameters: Optional[Sequence[str]] = None,
) -> List[JointLpmFit]:
    """Fit and rank several LPM families by AIC."""
    model_list = list(models or SUPPORTED_LPM_MODELS[:6])
    fits = [
        fit_lpm_model(
            observations,
            sample_year=sample_year,
            model=model,
            histories=histories,
            history_paths=history_paths,
            max_age_years=max_age_years,
            age_steps=age_steps,
            age_profile_objective_window=age_profile_objective_window,
            use_helium4=use_helium4,
            initial_c14_pmc=initial_c14_pmc,
            refine=refine,
            age_fraction_obs=age_fraction_obs,
            prediction_scale_factors=prediction_scale_factors,
            parameter_template=parameter_template,
            free_parameters=free_parameters,
        )
        for model in model_list
    ]
    return sorted(
        fits,
        key=lambda fit: (
            not fit.converged,
            fit.aic if math.isfinite(fit.aic) else float("inf"),
            fit.objective if math.isfinite(fit.objective) else float("inf"),
        ),
    )
