"""Dissolved-gas recharge-condition and tracer-correction models.

The implementation follows the same modelling structure used in USGS dissolved
gas workflows: estimate recharge conditions from conservative dissolved gases,
then use those conditions to convert environmental tracers to atmospheric
equivalent inputs for age modelling. It is not a byte-for-byte DGMETA clone; the
default solubility constants are compact approximations and can be replaced with
study-specific calibration tables.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .tracer_inputs import GAS_TRACERS, standardize_gas_observations

MOLAR_VOLUME_CC_STP_PER_MOL = 22414.0
R_GAS_L_ATM_MOL_K = 0.082057366

MOLAR_MASS_G_MOL = {
    "HE": 4.002602,
    "NE": 20.1797,
    "AR": 39.948,
    "KR": 83.798,
    "XE": 131.293,
    "N2": 28.0134,
    "O2": 31.9988,
    "SF6": 146.055,
    "CFC11": 137.368,
    "CFC12": 120.913,
    "CFC113": 187.376,
}

AIR_MOLE_FRACTION = {
    "HE": 5.24e-6,
    "NE": 18.18e-6,
    "AR": 9.34e-3,
    "KR": 1.14e-6,
    "XE": 0.087e-6,
    "N2": 0.78084,
    "O2": 0.20946,
    "SF6": 1.0e-12,
    "CFC11": 1.0e-12,
    "CFC12": 1.0e-12,
    "CFC113": 1.0e-12,
}

# Approximate water-equilibrium concentrations at 20 C and 1 atm dry air in
# ccSTP/g water. Values are intentionally replaceable by user solubility tables.
DEFAULT_SOLUBILITY_CCSTP_G_AT_20C = {
    "HE": 4.6e-8,
    "NE": 1.9e-7,
    "AR": 3.2e-4,
    "KR": 7.3e-8,
    "XE": 1.05e-8,
    "N2": 1.50e-2,
    "O2": 7.7e-3,
    # Per pptv atmospheric mixing ratio, not ambient modern concentration.
    "SF6": 2.4e-15,
    "CFC11": 1.1e-12,
    "CFC12": 4.8e-13,
    "CFC113": 8.0e-13,
}

# Positive values increase solubility in colder water.
DEFAULT_TEMP_COEFF_PER_C = {
    "HE": 0.012,
    "NE": 0.018,
    "AR": 0.030,
    "KR": 0.040,
    "XE": 0.052,
    "N2": 0.025,
    "O2": 0.030,
    "SF6": 0.045,
    "CFC11": 0.050,
    "CFC12": 0.047,
    "CFC113": 0.055,
}


@dataclass(frozen=True)
class DissolvedGasObservation:
    gas: str
    value: float
    sigma: float
    unit: str = "ccSTP/g"
    include_in_fit: bool = True


@dataclass(frozen=True)
class DissolvedGasResidual:
    gas: str
    observed_ccstp_g: float
    predicted_ccstp_g: float
    sigma_ccstp_g: float
    standardized_residual: float


@dataclass(frozen=True)
class DissolvedGasFit:
    model: str
    recharge_temperature_c: float
    excess_air_ccstp_g: float
    recharge_pressure_atm: float
    salinity_ppt: float
    fractionation: float
    excess_n2_ccstp_g: float
    objective: float
    rmse_standardized: float
    aic: float
    bic: float
    n_gases: int
    residuals: List[DissolvedGasResidual]
    corrected_tracers: Dict[str, float]
    corrected_tracer_sigmas: Dict[str, float]
    flags: List[str]
    converged: bool
    note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "model": self.model,
            "recharge_temperature_c": float(self.recharge_temperature_c),
            "excess_air_ccstp_g": float(self.excess_air_ccstp_g),
            "recharge_pressure_atm": float(self.recharge_pressure_atm),
            "salinity_ppt": float(self.salinity_ppt),
            "fractionation": float(self.fractionation),
            "excess_n2_ccstp_g": float(self.excess_n2_ccstp_g),
            "objective": float(self.objective),
            "rmse_standardized": float(self.rmse_standardized),
            "aic": float(self.aic),
            "bic": float(self.bic),
            "n_gases": int(self.n_gases),
            "corrected_tracers": dict(self.corrected_tracers),
            "corrected_tracer_sigmas": dict(self.corrected_tracer_sigmas),
            "flags": list(self.flags),
            "converged": bool(self.converged),
            "note": self.note,
            "residuals": [r.__dict__ for r in self.residuals],
        }


def _gas_key(gas: str) -> str:
    return str(gas).strip().upper().replace("-", "").replace("_", "")


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def water_vapor_pressure_atm(temp_c: float) -> float:
    """Saturation water vapor pressure in atm using Tetens' approximation."""
    temp = float(temp_c)
    kpa = 0.61078 * math.exp((17.2694 * temp) / (temp + 237.29))
    return float(kpa / 101.325)


def pressure_from_elevation_atm(elevation_m: float) -> float:
    """Approximate atmospheric pressure from elevation."""
    elev = float(elevation_m)
    return float((1.0 - 2.25577e-5 * elev) ** 5.25588)


def ccstp_g_from_mg_l(value_mg_l: float, gas: str, water_density_kg_l: float = 1.0) -> float:
    """Convert mg/L dissolved gas to ccSTP/g water."""
    key = _gas_key(gas)
    if key not in MOLAR_MASS_G_MOL:
        raise ValueError(f"No molar mass available for gas {gas!r}.")
    mol_l = float(value_mg_l) * 1.0e-3 / MOLAR_MASS_G_MOL[key]
    cc_l = mol_l * MOLAR_VOLUME_CC_STP_PER_MOL
    return float(cc_l / (1000.0 * float(water_density_kg_l)))


def mg_l_from_ccstp_g(value_ccstp_g: float, gas: str, water_density_kg_l: float = 1.0) -> float:
    """Convert ccSTP/g water to mg/L dissolved gas."""
    key = _gas_key(gas)
    if key not in MOLAR_MASS_G_MOL:
        raise ValueError(f"No molar mass available for gas {gas!r}.")
    mol_g = float(value_ccstp_g) / MOLAR_VOLUME_CC_STP_PER_MOL
    g_g = mol_g * MOLAR_MASS_G_MOL[key]
    return float(g_g * 1000.0 * 1000.0 * float(water_density_kg_l))


def to_ccstp_g(value: float, gas: str, unit: str = "ccSTP/g") -> float:
    """Normalize a dissolved-gas concentration to ccSTP/g water."""
    unit_key = str(unit).strip().lower().replace(" ", "")
    if unit_key in {"ccstp/g", "cm3stp/g", "cm^3stp/g", "mlstp/g"}:
        return float(value)
    if unit_key in {"ccstp/kg", "cm3stp/kg", "mlstp/kg"}:
        return float(value) / 1000.0
    if unit_key in {"mg/l", "mgperliter", "mgl"}:
        return ccstp_g_from_mg_l(float(value), gas)
    if unit_key in {"umol/l", "µmol/l", "umolperliter"}:
        cc_l = float(value) * 1.0e-6 * MOLAR_VOLUME_CC_STP_PER_MOL
        return cc_l / 1000.0
    raise ValueError(f"Unsupported dissolved gas unit {unit!r}.")


def equilibrium_concentration_ccstp_g(
    gas: str,
    temp_c: float,
    *,
    salinity_ppt: float = 0.0,
    pressure_atm: float = 1.0,
    solubility_table: Optional[Mapping[str, float]] = None,
) -> float:
    """Approximate air-water equilibrium concentration in ccSTP/g."""
    key = _gas_key(gas)
    solubilities = dict(DEFAULT_SOLUBILITY_CCSTP_G_AT_20C)
    if solubility_table:
        solubilities.update({_gas_key(k): float(v) for k, v in solubility_table.items()})
    if key not in solubilities:
        raise ValueError(f"No solubility available for gas {gas!r}.")
    base = float(solubilities[key])
    coeff = float(DEFAULT_TEMP_COEFF_PER_C.get(key, 0.03))
    dry_pressure = max(float(pressure_atm) - water_vapor_pressure_atm(float(temp_c)), 0.05)
    pressure_factor = dry_pressure / max(1.0 - water_vapor_pressure_atm(20.0), 0.05)
    salinity_factor = math.exp(-0.010 * max(float(salinity_ppt), 0.0))
    return float(base * math.exp(coeff * (20.0 - float(temp_c))) * pressure_factor * salinity_factor)


def _ua_excess_air_component(gas: str, excess_air_ccstp_g: float) -> float:
    key = _gas_key(gas)
    x_air = AIR_MOLE_FRACTION.get(key)
    if x_air is None:
        raise ValueError(f"No atmospheric mole fraction available for gas {gas!r}.")
    return float(max(excess_air_ccstp_g, 0.0) * x_air)


def _ce_excess_air_component(
    gas: str,
    excess_air_ccstp_g: float,
    equilibrium_ccstp_g: float,
    fractionation: float,
) -> float:
    key = _gas_key(gas)
    x_air = AIR_MOLE_FRACTION.get(key)
    if x_air is None:
        raise ValueError(f"No atmospheric mole fraction available for gas {gas!r}.")
    a = max(float(excess_air_ccstp_g), 0.0)
    f = min(max(float(fractionation), 0.0), 0.999999)
    eq = max(float(equilibrium_ccstp_g), 1.0e-30)
    return float(((1.0 - f) * a * x_air) / (1.0 + (f * a * x_air / eq)))


def _excess_air_component(
    gas: str,
    excess_air_ccstp_g: float,
    model: str,
    fractionation: float,
    *,
    equilibrium_ccstp_g: Optional[float] = None,
) -> float:
    """Return UA or CE gas-specific excess-air contribution in ccSTP/g."""
    model_key = model.upper()
    if model_key == "UA":
        return _ua_excess_air_component(gas, excess_air_ccstp_g)
    if model_key == "CE":
        if equilibrium_ccstp_g is None:
            raise ValueError("CE excess-air calculation requires equilibrium_ccstp_g.")
        return _ce_excess_air_component(gas, excess_air_ccstp_g, equilibrium_ccstp_g, fractionation)
    raise ValueError("Dissolved-gas model must be `UA` or `CE`.")


def excess_air_concentration_ccstp_g(
    gas: str,
    excess_air_ccstp_g: float,
    *,
    model: str = "UA",
    fractionation: float = 0.0,
    equilibrium_ccstp_g: Optional[float] = None,
) -> float:
    """Return the gas-specific excess-air contribution in ccSTP/g."""
    return _excess_air_component(
        gas,
        excess_air_ccstp_g,
        model,
        fractionation,
        equilibrium_ccstp_g=equilibrium_ccstp_g,
    )


def predict_dissolved_gas_ccstp_g(
    gas: str,
    *,
    recharge_temperature_c: float,
    excess_air_ccstp_g: float = 0.0,
    model: str = "UA",
    fractionation: float = 0.0,
    salinity_ppt: float = 0.0,
    pressure_atm: float = 1.0,
    excess_n2_ccstp_g: float = 0.0,
    solubility_table: Optional[Mapping[str, float]] = None,
) -> float:
    """Predict dissolved gas from recharge conditions."""
    key = _gas_key(gas)
    eq = equilibrium_concentration_ccstp_g(
        key,
        recharge_temperature_c,
        salinity_ppt=salinity_ppt,
        pressure_atm=pressure_atm,
        solubility_table=solubility_table,
    )
    excess = _excess_air_component(key, excess_air_ccstp_g, model, fractionation, equilibrium_ccstp_g=eq)
    denit = float(excess_n2_ccstp_g) if key == "N2" else 0.0
    return float(eq + excess + denit)


def build_dissolved_gas_observations(
    observations: Mapping[str, Any],
    *,
    default_relative_sigma: float = 0.03,
) -> List[DissolvedGasObservation]:
    """Extract dissolved-gas observations from a mapping."""
    out: List[DissolvedGasObservation] = []
    gases = ("HE", "NE", "AR", "KR", "XE", "N2", "O2")
    for gas in gases:
        lower = gas.lower()
        value = None
        unit = str(observations.get(f"{lower}_unit") or "ccSTP/g")
        for key in (f"{lower}_ccstp_g", f"{lower}_dissolved", lower):
            value = _finite_float(observations.get(key))
            if value is not None:
                if key.endswith("_ccstp_g"):
                    unit = "ccSTP/g"
                break
        if value is None or value < 0:
            continue
        converted = to_ccstp_g(value, gas, unit)
        sigma_value = None
        for key in (f"{lower}_sigma_ccstp_g", f"{lower}_dissolved_sigma", f"{lower}_sigma"):
            sigma_value = _finite_float(observations.get(key))
            if sigma_value is not None:
                sigma_unit = "ccSTP/g" if key.endswith("_ccstp_g") else unit
                sigma_value = to_ccstp_g(abs(sigma_value), gas, sigma_unit)
                break
        sigma = float(max(sigma_value or 0.0, default_relative_sigma * max(converted, 1.0e-12), 1.0e-12))
        out.append(DissolvedGasObservation(gas=gas, value=converted, sigma=sigma))
    return out


def _grid(values: Optional[Sequence[float]], default: Sequence[float]) -> List[float]:
    return [float(v) for v in (values if values is not None else default)]


def _fit_one_model(
    observations: Sequence[DissolvedGasObservation],
    *,
    model: str,
    temp_grid_c: Sequence[float],
    excess_air_grid_ccstp_g: Sequence[float],
    fractionation_grid: Sequence[float],
    excess_n2_grid_ccstp_g: Sequence[float],
    salinity_ppt: float,
    pressure_atm: float,
    solubility_table: Optional[Mapping[str, float]],
) -> DissolvedGasFit:
    best: Optional[Tuple[float, float, float, float, float, List[DissolvedGasResidual]]] = None
    model_key = model.upper()
    frac_values = [0.0] if model_key == "UA" else list(fractionation_grid)
    n2_values = list(excess_n2_grid_ccstp_g) if any(obs.gas == "N2" for obs in observations) else [0.0]

    for temp in temp_grid_c:
        for excess_air in excess_air_grid_ccstp_g:
            for frac in frac_values:
                for excess_n2 in n2_values:
                    residuals: List[DissolvedGasResidual] = []
                    objective = 0.0
                    for obs in observations:
                        pred = predict_dissolved_gas_ccstp_g(
                            obs.gas,
                            recharge_temperature_c=temp,
                            excess_air_ccstp_g=excess_air,
                            model=model_key,
                            fractionation=frac,
                            salinity_ppt=salinity_ppt,
                            pressure_atm=pressure_atm,
                            excess_n2_ccstp_g=excess_n2,
                            solubility_table=solubility_table,
                        )
                        std = (obs.value - pred) / obs.sigma
                        objective += std * std
                        residuals.append(
                            DissolvedGasResidual(
                                gas=obs.gas,
                                observed_ccstp_g=float(obs.value),
                                predicted_ccstp_g=float(pred),
                                sigma_ccstp_g=float(obs.sigma),
                                standardized_residual=float(std),
                            )
                        )
                    if best is None or objective < best[0]:
                        best = (float(objective), float(temp), float(excess_air), float(frac), float(excess_n2), residuals)

    if best is None:
        return DissolvedGasFit(
            model=model_key,
            recharge_temperature_c=float("nan"),
            excess_air_ccstp_g=float("nan"),
            recharge_pressure_atm=float(pressure_atm),
            salinity_ppt=float(salinity_ppt),
            fractionation=float("nan"),
            excess_n2_ccstp_g=float("nan"),
            objective=float("nan"),
            rmse_standardized=float("nan"),
            aic=float("nan"),
            bic=float("nan"),
            n_gases=len(observations),
            residuals=[],
            corrected_tracers={},
            corrected_tracer_sigmas={},
            flags=["no_finite_solution"],
            converged=False,
        )

    objective, temp, excess_air, frac, excess_n2, residuals = best
    n = len(residuals)
    k = 2 + (1 if model_key == "CE" else 0) + (1 if any(obs.gas == "N2" for obs in observations) else 0)
    rss_per_obs = max(objective / max(n, 1), 1.0e-12)
    aic = n * math.log(rss_per_obs) + 2.0 * k
    bic = n * math.log(rss_per_obs) + k * math.log(max(n, 2))
    flags = _diagnostic_flags(residuals, temp, excess_air, frac, n)
    return DissolvedGasFit(
        model=model_key,
        recharge_temperature_c=temp,
        excess_air_ccstp_g=excess_air,
        recharge_pressure_atm=float(pressure_atm),
        salinity_ppt=float(salinity_ppt),
        fractionation=frac,
        excess_n2_ccstp_g=excess_n2,
        objective=objective,
        rmse_standardized=math.sqrt(rss_per_obs),
        aic=float(aic),
        bic=float(bic),
        n_gases=n,
        residuals=residuals,
        corrected_tracers={},
        corrected_tracer_sigmas={},
        flags=flags,
        converged=True,
        note="Dissolved-gas model with UA and closed-system-equilibration CE excess-air options.",
    )


def _diagnostic_flags(
    residuals: Sequence[DissolvedGasResidual],
    temp_c: float,
    excess_air_ccstp_g: float,
    fractionation: float,
    n_gases: int,
) -> List[str]:
    flags: List[str] = []
    if n_gases < 2:
        flags.append("underconstrained_less_than_two_gases")
    if n_gases < 4:
        flags.append("limited_gas_suite")
    if temp_c <= 0.0 or temp_c >= 35.0:
        flags.append("recharge_temperature_at_grid_boundary_or_extreme")
    if excess_air_ccstp_g <= 0.0:
        flags.append("zero_or_negative_excess_air_solution")
    if fractionation >= 0.95:
        flags.append("strong_fractionation_solution")
    if any(abs(r.standardized_residual) > 3.0 for r in residuals):
        flags.append("large_standardized_residual")
    if any(r.gas == "N2" and r.standardized_residual > 3.0 for r in residuals):
        flags.append("possible_excess_n2_or_denitrification")
    return flags


def fit_dissolved_gas_model(
    observations: Mapping[str, Any],
    *,
    models: Sequence[str] = ("UA", "CE"),
    temperature_range_c: Tuple[float, float] = (0.0, 35.0),
    temperature_step_c: float = 0.5,
    excess_air_range_ccstp_g: Tuple[float, float] = (0.0, 0.02),
    excess_air_steps: int = 61,
    fractionation_grid: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
    excess_n2_range_ccstp_g: Tuple[float, float] = (0.0, 0.02),
    excess_n2_steps: int = 21,
    salinity_ppt: Optional[float] = None,
    elevation_m: Optional[float] = None,
    pressure_atm: Optional[float] = None,
    solubility_table: Optional[Mapping[str, float]] = None,
) -> List[DissolvedGasFit]:
    """Fit recharge temperature/excess air from dissolved-gas observations."""
    obs = build_dissolved_gas_observations(observations)
    if not obs:
        return [
            DissolvedGasFit(
                model=str(models[0]).upper() if models else "UA",
                recharge_temperature_c=float("nan"),
                excess_air_ccstp_g=float("nan"),
                recharge_pressure_atm=float("nan"),
                salinity_ppt=float("nan"),
                fractionation=float("nan"),
                excess_n2_ccstp_g=float("nan"),
                objective=float("nan"),
                rmse_standardized=float("nan"),
                aic=float("nan"),
                bic=float("nan"),
                n_gases=0,
                residuals=[],
                corrected_tracers={},
                corrected_tracer_sigmas={},
                flags=["no_dissolved_gas_observations"],
                converged=False,
            )
        ]

    salinity = float(salinity_ppt if salinity_ppt is not None else (_finite_float(observations.get("salinity_ppt")) or 0.0))
    if pressure_atm is None:
        elevation = float(elevation_m if elevation_m is not None else (_finite_float(observations.get("elevation_m")) or 0.0))
        pressure = pressure_from_elevation_atm(elevation)
    else:
        pressure = float(pressure_atm)

    temp_grid = list(np.arange(temperature_range_c[0], temperature_range_c[1] + temperature_step_c, temperature_step_c))
    excess_air_grid = list(np.linspace(excess_air_range_ccstp_g[0], excess_air_range_ccstp_g[1], max(2, excess_air_steps)))
    excess_n2_grid = list(np.linspace(excess_n2_range_ccstp_g[0], excess_n2_range_ccstp_g[1], max(1, excess_n2_steps)))
    fits = [
        _fit_one_model(
            obs,
            model=model,
            temp_grid_c=temp_grid,
            excess_air_grid_ccstp_g=excess_air_grid,
            fractionation_grid=fractionation_grid,
            excess_n2_grid_ccstp_g=excess_n2_grid,
            salinity_ppt=salinity,
            pressure_atm=pressure,
            solubility_table=solubility_table,
        )
        for model in models
    ]
    return sorted(
        fits,
        key=lambda fit: (
            not fit.converged,
            fit.aic if math.isfinite(fit.aic) else float("inf"),
            fit.objective if math.isfinite(fit.objective) else float("inf"),
        ),
    )


def atmospheric_equivalent_from_dissolved(
    gas: str,
    dissolved_ccstp_g: float,
    fit: DissolvedGasFit,
    *,
    solubility_table: Optional[Mapping[str, float]] = None,
) -> float:
    """Convert dissolved environmental tracer concentration to dry-air pptv."""
    key = _gas_key(gas)
    eq_per_pptv = equilibrium_concentration_ccstp_g(
        key,
        fit.recharge_temperature_c,
        salinity_ppt=fit.salinity_ppt,
        pressure_atm=fit.recharge_pressure_atm,
        solubility_table=solubility_table,
    )
    if eq_per_pptv <= 0:
        raise ValueError(f"No positive environmental-tracer solubility for {gas!r}.")
    excess_per_pptv = _excess_air_component(
        key,
        fit.excess_air_ccstp_g,
        fit.model,
        fit.fractionation,
        equilibrium_ccstp_g=eq_per_pptv,
    )
    total_per_pptv = eq_per_pptv + excess_per_pptv
    if total_per_pptv <= 0:
        raise ValueError(f"No positive environmental-tracer dissolved response for {gas!r}.")
    return float(max(float(dissolved_ccstp_g), 0.0) / total_per_pptv)


def correct_environmental_tracers(
    observations: Mapping[str, Any],
    fit: DissolvedGasFit,
    *,
    solubility_table: Optional[Mapping[str, float]] = None,
) -> Dict[str, Any]:
    """Add corrected atmospheric-equivalent tracer fields to observations."""
    out = standardize_gas_observations(observations)
    corrected: Dict[str, float] = {}
    sigmas: Dict[str, float] = {}
    for tracer in GAS_TRACERS:
        lower = tracer.lower()
        raw = _finite_float(out.get(f"{lower}_dissolved"))
        if raw is None:
            continue
        unit = str(out.get(f"{lower}_unit") or "ccSTP/g")
        raw_cc = to_ccstp_g(raw, tracer, unit)
        pptv = atmospheric_equivalent_from_dissolved(tracer, raw_cc, fit, solubility_table=solubility_table)
        out[f"{lower}_pptv"] = pptv
        corrected[tracer] = pptv
        sigma_raw = _finite_float(out.get(f"{lower}_dissolved_sigma"))
        if sigma_raw is not None:
            sigma_cc = to_ccstp_g(abs(sigma_raw), tracer, unit)
            eq_per_pptv = equilibrium_concentration_ccstp_g(
                tracer,
                fit.recharge_temperature_c,
                salinity_ppt=fit.salinity_ppt,
                pressure_atm=fit.recharge_pressure_atm,
                solubility_table=solubility_table,
            )
            excess_per_pptv = _excess_air_component(
                tracer,
                fit.excess_air_ccstp_g,
                fit.model,
                fit.fractionation,
                equilibrium_ccstp_g=eq_per_pptv,
            )
            total_per_pptv = eq_per_pptv + excess_per_pptv
            sigmas[tracer] = float(sigma_cc / total_per_pptv)
            out[f"{lower}_sigma_pptv"] = sigmas[tracer]
    out["dissolved_gas_model"] = fit.to_dict()
    out["dissolved_gas_model"]["corrected_tracers"] = corrected
    out["dissolved_gas_model"]["corrected_tracer_sigmas"] = sigmas
    return out


def fit_and_correct_dissolved_gases(
    observations: Mapping[str, Any],
    **kwargs: Any,
) -> Tuple[DissolvedGasFit, Dict[str, Any]]:
    """Fit the best dissolved-gas model and return corrected observations."""
    fits = fit_dissolved_gas_model(observations, **kwargs)
    best = fits[0]
    corrected = correct_environmental_tracers(observations, best, solubility_table=kwargs.get("solubility_table"))
    best = DissolvedGasFit(
        **{
            **best.__dict__,
            "corrected_tracers": corrected.get("dissolved_gas_model", {}).get("corrected_tracers", {}),
            "corrected_tracer_sigmas": corrected.get("dissolved_gas_model", {}).get("corrected_tracer_sigmas", {}),
        }
    )
    return best, corrected


def monte_carlo_dissolved_gas_fit(
    observations: Mapping[str, Any],
    *,
    n: int = 100,
    seed: int = 42,
    **fit_kwargs: Any,
) -> Dict[str, Any]:
    """Propagate dissolved-gas analytical uncertainty through model fitting."""
    base_obs = build_dissolved_gas_observations(observations)
    if not base_obs:
        return {"n": 0, "fits": [], "summary": {}, "flags": ["no_dissolved_gas_observations"]}
    rng = np.random.default_rng(seed)
    fits: List[DissolvedGasFit] = []
    for _ in range(max(1, int(n))):
        sampled = dict(observations)
        for obs in base_obs:
            lower = obs.gas.lower()
            sampled[f"{lower}_ccstp_g"] = max(0.0, float(rng.normal(obs.value, obs.sigma)))
            sampled[f"{lower}_sigma_ccstp_g"] = obs.sigma
        fit = fit_dissolved_gas_model(sampled, **fit_kwargs)[0]
        if fit.converged:
            fits.append(fit)
    if not fits:
        return {"n": 0, "fits": [], "summary": {}, "flags": ["no_converged_monte_carlo_fits"]}

    def q(values: Sequence[float]) -> Dict[str, float]:
        arr = np.array(values, dtype=float)
        return {
            "p05": float(np.percentile(arr, 5)),
            "p50": float(np.percentile(arr, 50)),
            "p95": float(np.percentile(arr, 95)),
        }

    summary = {
        "recharge_temperature_c": q([fit.recharge_temperature_c for fit in fits]),
        "excess_air_ccstp_g": q([fit.excess_air_ccstp_g for fit in fits]),
        "fractionation": q([fit.fractionation for fit in fits]),
        "excess_n2_ccstp_g": q([fit.excess_n2_ccstp_g for fit in fits]),
    }
    return {
        "n": len(fits),
        "summary": summary,
        "models": {model: sum(1 for fit in fits if fit.model == model) for model in sorted({fit.model for fit in fits})},
        "flags": sorted({flag for fit in fits for flag in fit.flags}),
        "fits": [fit.to_dict() for fit in fits],
    }
