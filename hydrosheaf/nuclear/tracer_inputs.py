"""Tracer-history loading and gas-observation normalization."""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Union

import numpy as np

from .input_history import InputHistory


GAS_TRACERS = ("SF6", "CFC11", "CFC12", "CFC113")
TRACER_ALIASES = {
    "3H": "3H",
    "H3": "3H",
    "TRITIUM": "3H",
    "SF6": "SF6",
    "CFC11": "CFC11",
    "CFC_11": "CFC11",
    "CFC-11": "CFC11",
    "CFC12": "CFC12",
    "CFC_12": "CFC12",
    "CFC-12": "CFC12",
    "CFC113": "CFC113",
    "CFC_113": "CFC113",
    "CFC-113": "CFC113",
    "14C": "14C",
    "C14": "14C",
    "CARBON14": "14C",
    "4HE": "4He",
    "HE4": "4He",
}


def normalize_tracer_key(tracer: str) -> str:
    """Return Hydrosheaf's canonical tracer key."""
    key = str(tracer).strip().upper().replace(" ", "").replace("/", "")
    return TRACER_ALIASES.get(key, TRACER_ALIASES.get(key.replace("-", "_"), str(tracer).strip()))


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _clean_column(name: str) -> str:
    return str(name).strip().lower().replace(" ", "_").replace("-", "_")


def _history_from_rows(rows: Sequence[Mapping[str, Any]]) -> Dict[str, InputHistory]:
    grouped: Dict[str, Dict[str, List[float]]] = {}
    for row in rows:
        tracer = normalize_tracer_key(str(row.get("tracer", "")))
        year = _finite_float(row.get("year"))
        value = _finite_float(row.get("value"))
        if not tracer or year is None or value is None:
            continue
        sigma = _finite_float(row.get("sigma"))
        bucket = grouped.setdefault(tracer, {"years": [], "values": [], "sigma": []})
        bucket["years"].append(float(year))
        bucket["values"].append(float(value))
        bucket["sigma"].append(float(sigma if sigma is not None else 0.0))
    return {
        tracer: InputHistory(
            np.array(data["years"], dtype=float),
            np.array(data["values"], dtype=float),
            np.array(data["sigma"], dtype=float),
        )
        for tracer, data in grouped.items()
        if data["years"]
    }


def _history_from_wide_table(rows: Sequence[Mapping[str, Any]]) -> Dict[str, InputHistory]:
    if not rows:
        return {}
    first = rows[0]
    tracer_columns = []
    for col in first:
        clean = _clean_column(col)
        if clean in {"year", "date"} or clean.endswith("_sigma"):
            continue
        tracer = normalize_tracer_key(col)
        if tracer:
            tracer_columns.append((col, tracer))

    histories: Dict[str, InputHistory] = {}
    for col, tracer in tracer_columns:
        years: List[float] = []
        values: List[float] = []
        sigmas: List[float] = []
        sigma_col = None
        for candidate in (f"{col}_sigma", f"{col}_sd", f"{col}_err"):
            if candidate in first:
                sigma_col = candidate
                break
        for row in rows:
            year = _finite_float(row.get("year") or row.get("date"))
            value = _finite_float(row.get(col))
            if year is None or value is None:
                continue
            sigma = _finite_float(row.get(sigma_col)) if sigma_col else None
            years.append(float(year))
            values.append(float(value))
            sigmas.append(float(sigma if sigma is not None else 0.0))
        if years:
            histories[tracer] = InputHistory(
                np.array(years, dtype=float),
                np.array(values, dtype=float),
                np.array(sigmas, dtype=float),
            )
    return histories


def load_tracer_histories_csv(path: Union[str, Path]) -> Dict[str, InputHistory]:
    """Load local tracer input histories from a CSV file.

    Supported formats:
    - long: `tracer,year,value,sigma`
    - wide: `year,3H,SF6,CFC12,...` with optional `<tracer>_sigma`
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = [{_clean_column(k): v for k, v in row.items()} for row in reader]
    if not rows:
        return {}
    columns = set(rows[0])
    if {"tracer", "year", "value"}.issubset(columns):
        return _history_from_rows(rows)
    if "year" in columns or "date" in columns:
        return _history_from_wide_table(rows)
    raise ValueError(
        "Tracer history CSV must be long (`tracer,year,value`) or wide (`year,<tracer>...`)."
    )


def load_tracer_histories(paths: Iterable[Union[str, Path]]) -> Dict[str, InputHistory]:
    """Load and merge local tracer histories from one or more CSV files."""
    histories: Dict[str, InputHistory] = {}
    for path in paths:
        histories.update(load_tracer_histories_csv(path))
    return histories


def dissolved_gas_to_atmospheric_equivalent(
    dissolved_value: float,
    solubility_value_per_pptv: float,
    *,
    excess_air_fraction: float = 0.0,
) -> float:
    """Convert dissolved gas to atmospheric-equivalent pptv with explicit solubility.

    Hydrosheaf does not infer gas solubility from temperature/salinity here. The
    caller must provide a tracer-specific solubility coefficient in the same
    dissolved-concentration units per pptv. Excess air is represented as a
    fractional additive term, so 0.1 means 10 percent excess-air contribution.
    """
    dissolved = _finite_float(dissolved_value)
    solubility = _finite_float(solubility_value_per_pptv)
    excess = _finite_float(excess_air_fraction) or 0.0
    if dissolved is None or solubility is None:
        raise ValueError("dissolved_value and solubility_value_per_pptv must be finite numbers.")
    if dissolved < 0:
        raise ValueError("dissolved_value must be non-negative.")
    if solubility <= 0:
        raise ValueError("solubility_value_per_pptv must be positive.")
    if excess <= -0.99:
        raise ValueError("excess_air_fraction must be greater than -0.99.")
    return float(dissolved / (solubility * (1.0 + excess)))


def standardize_gas_observations(observations: Mapping[str, Any]) -> Dict[str, Any]:
    """Return observations with canonical atmospheric-equivalent gas fields.

    For each gas tracer, Hydrosheaf first keeps existing `<tracer>_pptv` values,
    then accepts `<tracer>_atm_equiv_pptv`, and finally computes `<tracer>_pptv`
    from raw fields only when both `<tracer>_dissolved` and
    `<tracer>_solubility_per_pptv` are supplied.
    """
    out = dict(observations)
    notes: List[str] = list(out.get("gas_correction_notes") or [])
    for tracer in GAS_TRACERS:
        lower = tracer.lower()
        canonical = f"{lower}_pptv"
        if _finite_float(out.get(canonical)) is not None:
            continue

        alias = f"{lower}_atm_equiv_pptv"
        alias_value = _finite_float(out.get(alias))
        if alias_value is not None:
            out[canonical] = alias_value
            notes.append(f"{tracer}: used provided atmospheric-equivalent field `{alias}`.")
            continue

        dissolved = _finite_float(out.get(f"{lower}_dissolved"))
        solubility = _finite_float(out.get(f"{lower}_solubility_per_pptv"))
        if dissolved is None or solubility is None:
            continue
        excess = _finite_float(out.get(f"{lower}_excess_air_fraction")) or 0.0
        out[canonical] = dissolved_gas_to_atmospheric_equivalent(
            dissolved,
            solubility,
            excess_air_fraction=excess,
        )
        sigma = _finite_float(out.get(f"{lower}_dissolved_sigma"))
        if sigma is not None:
            out[f"{lower}_sigma_pptv"] = abs(
                dissolved_gas_to_atmospheric_equivalent(
                    sigma,
                    solubility,
                    excess_air_fraction=excess,
                )
            )
        notes.append(
            f"{tracer}: converted dissolved gas with user-supplied solubility/excess-air fields."
        )
    if notes:
        out["gas_correction_notes"] = notes
    return out
