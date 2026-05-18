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

        # Phase 4: mark supersaturated gas as upper-censored
        max_hist = historical_max_concentration(tracer, sample_year)
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


def _integration_grid(mean_age: float, max_age_years: float) -> Tuple[np.ndarray, float]:
    upper = max(100.0, 6.0 * max(mean_age, 1.0), max_age_years)
    upper = min(upper, max(max_age_years, mean_age * 8.0, 100.0))
    if upper <= 250.0:
        d_tau = 0.25
    elif upper <= 2000.0:
        d_tau = 1.0
    else:
        d_tau = 10.0
    taus = np.arange(0.0, upper + d_tau, d_tau)
    return taus, d_tau


def _normalize_pdf(g: np.ndarray, d_tau: float) -> np.ndarray:
    g = np.nan_to_num(g, nan=0.0, posinf=0.0, neginf=0.0)
    g[g < 0] = 0.0
    area = float(np.sum(g) * d_tau)
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
        return _normalize_pdf(g, float(np.diff(taus[:2])[0]) if len(taus) > 1 else 1.0)

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
    d_tau: float,
) -> np.ndarray:
    """Return a normalized transit-time distribution for a model family."""
    model_key = model.upper()
    if model_key.startswith("BMM-"):
        parts = model_key.split("-")
        if len(parts) != 3:
            raise ValueError("Binary-mixture model names must look like BMM-DM-DM.")
        fraction = min(max(float(parameters.get("binary_fraction", 0.5)), 0.0), 1.0)
        g1 = _component_pdf(parts[1], float(parameters["mean_age_1_years"]), taus, parameters)
        params2 = {
            key[:-2] if key.endswith("_2") else key: value
            for key, value in parameters.items()
        }
        g2 = _component_pdf(parts[2], float(parameters["mean_age_2_years"]), taus, params2)
        return _normalize_pdf(fraction * g1 + (1.0 - fraction) * g2, d_tau)

    g = _component_pdf(model_key, float(parameters["mean_age_years"]), taus, parameters)
    return _normalize_pdf(g, d_tau)


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
    taus, d_tau = _integration_grid(age_scale, max_age_years)
    pdf = transit_time_pdf(model_key, parameters, taus, d_tau)

    # Simple fixed-age cutoffs for first implementation
    anthropocene_mask = taus <= 70.0
    holocene_mask = (taus > 70.0) & (taus <= 11700.0)
    pleistocene_mask = taus > 11700.0

    frac_anthro = float(np.sum(pdf[anthropocene_mask]) * d_tau) if np.any(anthropocene_mask) else 0.0
    frac_holocene = float(np.sum(pdf[holocene_mask]) * d_tau) if np.any(holocene_mask) else 0.0
    frac_pleistocene = float(np.sum(pdf[pleistocene_mask]) * d_tau) if np.any(pleistocene_mask) else 0.0
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
    d_tau: float,
    decay_lambda_per_year: float = 0.0,
    *,
    daughter: bool = False,
) -> float:
    recharge_years = float(sample_year) - taus
    input_values = np.interp(recharge_years, history.years, history.values, left=0.0)
    decay = np.exp(-decay_lambda_per_year * taus)
    if daughter:
        values = input_values * (1.0 - decay)
    else:
        values = input_values * decay
    return float(np.sum(values * pdf) * d_tau)


def _predict_c14(pdf: np.ndarray, taus: np.ndarray, d_tau: float, initial_pmc: float) -> float:
    lambda_y = math.log(2.0) / CARBON14.half_life_years
    return float(np.sum(initial_pmc * np.exp(-lambda_y * taus) * pdf) * d_tau)


def _predict_ar39(pdf: np.ndarray, taus: np.ndarray, d_tau: float, initial_pmc: float) -> float:
    lambda_y = math.log(2.0) / ARGON39.half_life_years
    return float(np.sum(initial_pmc * np.exp(-lambda_y * taus) * pdf) * d_tau)


def _mean_from_distribution(pdf: np.ndarray, taus: np.ndarray, d_tau: float) -> float:
    return float(np.sum(taus * pdf) * d_tau)


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
) -> Dict[str, float]:
    """Predict tracer values from one LPM parameter set."""
    model_key = model.upper()
    if model_key.startswith("BMM-"):
        age_scale = max(
            float(parameters.get("mean_age_1_years", 1.0)),
            float(parameters.get("mean_age_2_years", 1.0)),
        )
    else:
        age_scale = float(parameters.get("mean_age_years", 1.0))
    taus, d_tau = _integration_grid(age_scale, max_age_years)
    pdf = transit_time_pdf(model_key, parameters, taus, d_tau)
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
                d_tau,
                decay_lambda_per_year=lambda_tritium,
            )
        if "3H/3HE" in requested:
            out["3H/3He"] = _predict_from_history(
                sample_year,
                tritium_hist,
                pdf,
                taus,
                d_tau,
                decay_lambda_per_year=lambda_tritium,
                daughter=True,
            )

    if "14C" in requested:
        out["14C"] = _predict_c14(pdf, taus, d_tau, initial_c14_pmc)

    if "39AR" in requested:
        # Standard initial Ar-39 activity is 100 pmc
        out["39Ar"] = _predict_ar39(pdf, taus, d_tau, 100.0)

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
                d_tau,
                decay_lambda_per_year=decay_lambda
            )
            if tracer == "85KR":
                # Rename back to 85Kr for consistency with observations
                out["85Kr"] = out.pop("85KR")

    if "4HE" in requested:
        mean_age = _mean_from_distribution(pdf, taus, d_tau)
        out["4He"] = helium4_background_ccpg + helium4_accumulation_rate_ccpg_per_year * mean_age

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
        age1 = p.get("mean_age_1_years", 0.01)
        age2 = p.get("mean_age_2_years", 0.01)
        if age2 <= age1:
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
    - Use only fits within delta_aicc <= max_delta_aicc.
    - If contributing model ages span more than max_log_age_span, use top model.
    - If fewer than two models pass gates, use top model.
    - Record whether BMA was used or skipped.
    """
    converged = [f for f in fits if f.converged and math.isfinite(f.aicc)]
    if not converged:
        return {
            "age_years": float("nan"),
            "bma_used": False,
            "bma_skip_reason": "no_converged_fits",
            "bma_n_models": 0,
            "bma_log_age_span": float("nan"),
            "top_model_age_years": float("nan"),
        }

    best_aicc = min(f.aicc for f in converged)
    within_delta = [f for f in converged if f.aicc - best_aicc <= max_delta_aicc]

    if len(within_delta) < 2:
        top = converged[0]
        top_age = _extract_scalar_age(top.parameters)
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "insufficient_models_within_delta_aicc",
            "bma_n_models": len(within_delta),
            "bma_log_age_span": float("nan"),
            "top_model_age_years": top_age,
        }

    ages = [_extract_scalar_age(f.parameters) for f in within_delta]
    ages = [a for a in ages if math.isfinite(a) and a > 0]
    if len(ages) < 2:
        top_age = ages[0] if ages else float("nan")
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "insufficient_finite_ages",
            "bma_n_models": len(within_delta),
            "bma_log_age_span": float("nan"),
            "top_model_age_years": top_age,
        }

    log_ages = [math.log10(a) for a in ages]
    log_span = max(log_ages) - min(log_ages)
    if log_span > max_log_age_span:
        top_age = ages[0]
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "log_age_span_exceeded",
            "bma_n_models": len(within_delta),
            "bma_log_age_span": log_span,
            "top_model_age_years": top_age,
        }

    # Compute AICc weights
    delta_aiccs = [f.aicc - best_aicc for f in within_delta]
    weights = [math.exp(-0.5 * d) for d in delta_aiccs]
    total_w = sum(weights)
    if total_w <= 0 or not math.isfinite(total_w):
        top_age = ages[0]
        return {
            "age_years": top_age,
            "bma_used": False,
            "bma_skip_reason": "weight_computation_failed",
            "bma_n_models": len(within_delta),
            "bma_log_age_span": log_span,
            "top_model_age_years": top_age,
        }

    bma_age = sum(w * a for w, a in zip(weights, ages)) / total_w
    return {
        "age_years": bma_age,
        "bma_used": True,
        "bma_skip_reason": "",
        "bma_n_models": len(within_delta),
        "bma_log_age_span": log_span,
        "top_model_age_years": ages[0],
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


def _refine_best_candidate(
    model: str,
    best_params: Dict[str, float],
    sample_year: float,
    tracer_names: list[str],
    fit_observations: list[TracerFitObservation],
    *,
    histories: Optional[Mapping[str, InputHistory]] = None,
    initial_c14_pmc: float = 100.0,
    max_age_years: float = 50000.0,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
) -> tuple[Dict[str, float], bool, str]:
    """Refine a grid-best candidate with local optimization.

    Returns (refined_params, success, message).
    """
    try:
        from scipy.optimize import minimize
    except Exception:
        return best_params, False, "scipy not available"

    x0, keys, bounds = _pack_parameters(model, best_params, max_age_years=max_age_years)

    def objective_fn(x: np.ndarray) -> float:
        params = _unpack_parameters(model, x, keys, bounds)
        predicted = predict_lpm_tracers(
            model, params, sample_year, tracer_names,
            histories=histories, initial_c14_pmc=initial_c14_pmc, max_age_years=max_age_years,
        )
        obj = 0.0
        for obs in fit_observations:
            pred = predicted.get(obs.tracer)
            if pred is None or not math.isfinite(pred):
                return 1e12
            obj += _standardized_observation_loss(obs, pred) * obs.weight
        if age_fraction_obs is not None:
            frac_pred = age_fraction_predictions(model, params, max_age_years=max_age_years)
            sigma_fraction = float(age_fraction_obs.get("sigma_fraction", 0.10))
            for key in ("anthropocene", "holocene", "pleistocene"):
                obs_val = age_fraction_obs.get(key)
                if obs_val is not None and math.isfinite(obs_val):
                    pred_val = frac_pred.get(key, 0.0)
                    resid = (obs_val - pred_val) / sigma_fraction
                    obj += resid * resid
        return obj

    # Nelder-Mead first
    try:
        res_nm = minimize(objective_fn, x0, method="Nelder-Mead", options={"maxiter": 300, "xatol": 1e-4, "fatol": 1e-4})
        if res_nm.success and res_nm.fun < objective_fn(x0):
            refined = _unpack_parameters(model, res_nm.x, keys, bounds)
            return refined, True, "Nelder-Mead refined"
        x0 = res_nm.x
    except Exception as exc:
        pass

    # L-BFGS-B fallback on bounded transformed parameters
    try:
        res_bfgs = minimize(objective_fn, x0, method="L-BFGS-B", bounds=bounds, options={"maxiter": 300})
        if res_bfgs.success and res_bfgs.fun < objective_fn(x0):
            refined = _unpack_parameters(model, res_bfgs.x, keys, bounds)
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
    use_helium4: bool = False,
    initial_c14_pmc: float = 100.0,
    refine: bool = False,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
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
    ages = _age_grid(max_age, age_steps)
    candidates, n_params = _parameter_candidates(model, ages)
    tracer_names = [obs.tracer for obs in fit_observations]

    best: Optional[Tuple[float, Dict[str, float], List[TracerFitResidual]]] = None
    for params in candidates:
        predicted = predict_lpm_tracers(
            model,
            params,
            sample_year,
            tracer_names,
            histories=resolved_histories,
            initial_c14_pmc=initial_c14_pmc,
            max_age_years=max_age,
        )
        
        residuals: List[TracerFitResidual] = []
        objective = 0.0
        missing_prediction = False
        for obs in fit_observations:
            pred = predicted.get(obs.tracer)
            if pred is None or not math.isfinite(pred):
                missing_prediction = True
                break
            loss = _standardized_observation_loss(obs, pred) * obs.weight
            objective += loss
            standardized = (obs.value - pred) / obs.sigma
            residuals.append(
                TracerFitResidual(
                    tracer=obs.tracer,
                    observed=float(obs.value),
                    predicted=float(pred),
                    sigma=float(obs.sigma),
                    standardized_residual=float(standardized),
                )
            )
        # Phase 6: age-fraction constraints
        if age_fraction_obs is not None:
            frac_pred = age_fraction_predictions(model, params, max_age_years=max_age)
            sigma_fraction = float(age_fraction_obs.get("sigma_fraction", 0.10))
            for key in ("anthropocene", "holocene", "pleistocene"):
                obs_val = age_fraction_obs.get(key)
                if obs_val is not None and math.isfinite(obs_val):
                    pred_val = frac_pred.get(key, 0.0)
                    resid = (obs_val - pred_val) / sigma_fraction
                    objective += resid * resid

        if missing_prediction:
            continue
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
            max_age_years=max_age,
            age_fraction_obs=age_fraction_obs,
        )
        if refinement_success:
            params = refined_params
            # Recompute objective with refined parameters
            predicted = predict_lpm_tracers(
                model,
                params,
                sample_year,
                tracer_names,
                histories=resolved_histories,
                initial_c14_pmc=initial_c14_pmc,
                max_age_years=max_age,
            )
            new_residuals: List[TracerFitResidual] = []
            new_objective = 0.0
            missing_prediction = False
            for obs in fit_observations:
                pred = predicted.get(obs.tracer)
                if pred is None or not math.isfinite(pred):
                    missing_prediction = True
                    break
                standardized = (obs.value - pred) / obs.sigma
                new_objective += (standardized * standardized) * obs.weight
                new_residuals.append(
                    TracerFitResidual(
                        tracer=obs.tracer,
                        observed=float(obs.value),
                        predicted=float(pred),
                        sigma=float(obs.sigma),
                        standardized_residual=float(standardized),
                    )
                )
            if not missing_prediction and new_objective < objective:
                objective = new_objective
                residuals = new_residuals

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
    use_helium4: bool = False,
    initial_c14_pmc: float = 100.0,
    refine: bool = False,
    age_fraction_obs: Optional[Mapping[str, float]] = None,
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
            use_helium4=use_helium4,
            initial_c14_pmc=initial_c14_pmc,
            age_fraction_obs=age_fraction_obs,
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
