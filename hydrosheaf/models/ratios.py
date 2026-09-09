"""Geochemical ratio diagnostics for sparse, multi-panel groundwater data.

Ratios are derived from the same canonical ``mmol/L`` observations used by
the reaction model.  They are reported as diagnostics and optional null/prior
evidence; they are *not* appended to the reaction matrix, which prevents
double-counting the numerator and denominator ions.  Missing or non-positive
concentrations remain missing rather than being imputed.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Optional


_MG_PER_MMOL = {
    "Ca": 40.078,
    "Mg": 24.305,
    "Na": 22.989769,
    "K": 39.0983,
    "HCO3": 61.0168,
    "Cl": 35.453,
    "SO4": 96.06,
    "NO3": 62.0049,
    "F": 18.998403,
    "Fe": 55.845,
    "PO4": 94.9714,
    "SiO2": 60.0843,
    "Sr": 87.62,
}


RATIO_DEFINITIONS: dict[str, tuple[str, str, float]] = {
    # name: numerator, denominator, denominator-equivalent multiplier
    "Na_Cl": ("Na", "Cl", 1.0),
    "Ca_Mg": ("Ca", "Mg", 1.0),
    "Mg_Ca": ("Mg", "Ca", 1.0),
    "K_Na": ("K", "Na", 1.0),
    "HCO3_CaMg_equiv": ("HCO3", "CaMg_equiv", 1.0),
    "Ca_Sr": ("Ca", "Sr", 1.0),
    "Sr_Ca": ("Sr", "Ca", 1.0),
    "SiO2_HCO3": ("SiO2", "HCO3", 1.0),
    "F_Ca": ("F", "Ca", 1.0),
    "NO3_Cl": ("NO3", "Cl", 1.0),
    "SO4_Cl": ("SO4", "Cl", 1.0),
}


def _finite(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _value_mmol(sample: Mapping[str, Any], ion: str) -> Optional[float]:
    """Read a canonical mmol/L value, with explicit raw mg/L fallbacks."""

    direct = _finite(sample.get(ion))
    if direct is not None and direct >= 0.0:
        return direct
    for key in (f"{ion}_mmol_L", f"{ion}_mmolL"):
        direct = _finite(sample.get(key))
        if direct is not None and direct >= 0.0:
            return direct
    for key in (f"{ion}_mg_L", f"{ion}_mgL", f"{ion}_mg/L"):
        raw = _finite(sample.get(key))
        if raw is not None and raw >= 0.0:
            return raw / _MG_PER_MMOL[ion]
    # The two native Northern Ghana diagnostic names are retained as accepted
    # fallbacks so ratio diagnostics can be run before harmonisation.
    native = {"Sr": "Sr_mg_L", "SiO2": "SiO2_mg_L"}.get(ion)
    if native:
        raw = _finite(sample.get(native))
        if raw is not None and raw >= 0.0:
            return raw / _MG_PER_MMOL[ion]
    return None


def _denominator(sample: Mapping[str, Any], name: str) -> Optional[float]:
    if name == "CaMg_equiv":
        ca = _value_mmol(sample, "Ca")
        mg = _value_mmol(sample, "Mg")
        if ca is None or mg is None:
            return None
        # Ca and Mg each carry two charge equivalents per mmol.
        return 2.0 * (ca + mg)
    return _value_mmol(sample, name)


@dataclass(frozen=True)
class RatioDiagnostics:
    """Ratios and explicit missingness for one observation."""

    values: Mapping[str, float]
    missing: tuple[str, ...]
    log_values: Mapping[str, float]

    @property
    def n_available(self) -> int:
        return len(self.values)

    def as_dict(self, *, include_logs: bool = True) -> dict[str, Any]:
        result: dict[str, Any] = {
            "values": dict(self.values),
            "missing": list(self.missing),
            "n_available": self.n_available,
        }
        if include_logs:
            result["log_values"] = dict(self.log_values)
        return result


def compute_geochemical_ratios(sample: Mapping[str, Any]) -> RatioDiagnostics:
    """Compute non-log and natural-log ratios from a canonical sample."""

    values: dict[str, float] = {}
    log_values: dict[str, float] = {}
    missing: list[str] = []
    for name, (numerator, denominator, multiplier) in RATIO_DEFINITIONS.items():
        num = _value_mmol(sample, numerator)
        den = _denominator(sample, denominator)
        if num is None or den is None or den <= 0.0:
            missing.append(name)
            continue
        ratio = num / (multiplier * den)
        if not math.isfinite(ratio) or ratio <= 0.0:
            missing.append(name)
            continue
        values[name] = float(ratio)
        log_values[name] = float(math.log(ratio))
    return RatioDiagnostics(values=values, missing=tuple(missing), log_values=log_values)


def log_ratio(sample: Mapping[str, Any], numerator: str, denominator: str) -> Optional[float]:
    """Return ``ln(numerator/denominator)`` or ``None`` when unidentifiable."""

    num = _value_mmol(sample, numerator)
    den = _denominator(sample, denominator)
    if num is None or den is None or num <= 0.0 or den <= 0.0:
        return None
    return float(math.log(num / den))


def compare_ratio_diagnostics(
    upstream: Mapping[str, Any], downstream: Mapping[str, Any]
) -> dict[str, Any]:
    """Compare endpoint ratios in log space without imputing missing pairs."""

    first = compute_geochemical_ratios(upstream)
    second = compute_geochemical_ratios(downstream)
    differences = {
        name: float(second.log_values[name] - first.log_values[name])
        for name in first.log_values.keys() & second.log_values.keys()
    }
    if differences:
        rmse = math.sqrt(sum(value * value for value in differences.values()) / len(differences))
        similarity = math.exp(-0.5 * rmse * rmse)
    else:
        rmse = None
        similarity = None
    return {
        "differences": differences,
        "n_pairs": len(differences),
        "logratio_rmse": rmse,
        "similarity": similarity,
        "missing_upstream": list(first.missing),
        "missing_downstream": list(second.missing),
    }


def augment_with_ratio_diagnostics(sample: Mapping[str, Any]) -> dict[str, Any]:
    """Return a copy with auditable ratio fields; never alters the input."""

    result = dict(sample)
    diagnostics = compute_geochemical_ratios(result)
    for name, value in diagnostics.values.items():
        result[f"ratio_{name}"] = value
    result["ratio_available_count"] = diagnostics.n_available
    result["ratio_missing"] = list(diagnostics.missing)
    return result


__all__ = [
    "RATIO_DEFINITIONS",
    "RatioDiagnostics",
    "augment_with_ratio_diagnostics",
    "compare_ratio_diagnostics",
    "compute_geochemical_ratios",
    "log_ratio",
]
