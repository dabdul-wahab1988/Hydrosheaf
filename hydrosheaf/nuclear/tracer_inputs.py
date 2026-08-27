"""Tracer-history loading and gas-observation normalization."""
from __future__ import annotations

import csv
import math
import os
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Union

import numpy as np

from .input_history import InputHistory, build_default_tritium_input


@dataclass(frozen=True)
class SiteInputContext:
    site_id: str
    sample_year: float
    latitude: float | None = None
    longitude: float | None = None
    study_unit: str = ""
    aquifer_group: str = ""
    recharge_temperature_c: float | None = None
    elevation_m: float | None = None
    # The national USGS release records the tritium input as a 2-degree by
    # 5-degree quadrangle label in Table 4.  Keeping that label with the site
    # context allows the primary benchmark to use the same public input family
    # instead of silently substituting a nearest-station reconstruction.
    tritium_input_source: str | None = None


USGS_QUADRANGLE_HISTORY_PATH = (
    Path(__file__).resolve().parents[2]
    / "M3.1"
    / "m3_age_benchmark"
    / "external"
    / "usgs_tritium_deposition_2023"
    / "Table_2_TritiumInQuadrangles.csv"
)


def _normalise_quadrangle_label(value: Any) -> str:
    """Normalise the label spelling used by USGS Table 4 and Table 2."""
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return ""
    return re.sub(r"\s+", "", text.replace(",", "_"))


def _quadrangle_from_coordinates(latitude: float | None, longitude: float | None) -> str:
    """Return the USGS quadrangle label containing a CONUS coordinate."""
    if latitude is None or longitude is None:
        return ""
    try:
        lat = float(latitude)
        lon_w = abs(float(longitude))
    except (TypeError, ValueError):
        return ""
    if not (math.isfinite(lat) and math.isfinite(lon_w)):
        return ""

    lat_bins = [
        (49.0, 47.0), (47.0, 45.0), (45.0, 43.0), (43.0, 41.0),
        (41.0, 39.0), (39.0, 37.0), (37.0, 35.0), (35.0, 33.0),
        (33.0, 31.0), (31.0, 29.0), (29.0, 27.0),
    ]
    lon_bins = [
        (125.0, 120.0), (120.0, 115.0), (115.0, 110.0),
        (110.0, 105.0), (105.0, 100.0), (100.0, 95.0),
        (95.0, 90.0), (90.0, 85.0), (85.0, 80.0), (80.0, 75.0),
        (75.0, 70.0), (70.0, 65.0),
    ]
    lat_label = ""
    if lat >= 49.0:
        lat_label = "49-47"
    elif lat < 27.0:
        lat_label = "<27"
    else:
        for upper, lower in lat_bins:
            if lower <= lat < upper:
                lat_label = f"{int(upper)}-{int(lower)}"
                break
    lon_label = ""
    if lon_w >= 125.0:
        lon_label = "125-120"
    elif lon_w < 65.0:
        lon_label = "<65"
    else:
        for west, east in lon_bins:
            if east <= lon_w < west:
                lon_label = f"{int(west)}-{int(east)}"
                break
    return f"{lat_label}_{lon_label}" if lat_label and lon_label else ""


@lru_cache(maxsize=1)
def _load_usgs_quadrangle_histories() -> dict[str, InputHistory]:
    """Load the public USGS monthly quadrangle histories once per process."""
    histories: dict[str, dict[str, list[float]]] = {}
    if not USGS_QUADRANGLE_HISTORY_PATH.exists():
        return {}
    with USGS_QUADRANGLE_HISTORY_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        quadrangles = [
            name for name in (reader.fieldnames or [])
            if name not in {"Date", "Year", "Month"}
        ]
        for name in quadrangles:
            histories[name] = {"years": [], "values": []}
        for row in reader:
            try:
                year = float(row.get("Year", ""))
                month = float(row.get("Month", ""))
            except (TypeError, ValueError):
                continue
            target_year = year + (month - 0.5) / 12.0
            for name in quadrangles:
                try:
                    value = float(row.get(name, ""))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(value):
                    histories[name]["years"].append(target_year)
                    histories[name]["values"].append(value)

    out: dict[str, InputHistory] = {}
    for name, data in histories.items():
        if not data["years"]:
            continue
        years = np.asarray(data["years"], dtype=float)
        values = np.asarray(data["values"], dtype=float)
        # The USGS release begins in August 1953.  A fixed pre-bomb background
        # is preferable to endpoint clamping for older age-grid nodes; the
        # background is negligible after tritium decay and is documented in the
        # run manifest as part of the input-history contract.
        pre_years = np.arange(1850.0, years.min(), 1.0 / 12.0)
        pre_values = np.full(pre_years.shape, 4.0, dtype=float)
        out[name] = InputHistory(
            np.concatenate([pre_years, years]),
            np.concatenate([pre_values, values]),
            sigma=np.zeros(len(pre_years) + len(years), dtype=float),
        )
    return out


def _usgs_quadrangle_history(context: SiteInputContext) -> tuple[InputHistory | None, str]:
    """Resolve the exact USGS quadrangle history named by a site or its coordinates."""
    histories = _load_usgs_quadrangle_histories()
    if not histories:
        return None, ""
    source = _normalise_quadrangle_label(context.tritium_input_source)
    if source not in histories:
        source = _quadrangle_from_coordinates(context.latitude, context.longitude)
    history = histories.get(source)
    return history, source


def build_site_tracer_histories(context: SiteInputContext) -> dict[str, InputHistory]:
    """Return site-aware tracer input histories.

    The primary M3.1 benchmark uses the USGS 2-degree by 5-degree quadrangle
    history named in Table 4.  Setting ``HYDROSHEAF_TRITIUM_HISTORY_MODE=wiser``
    explicitly selects the prior nearest-WISER reconstruction for sensitivity
    analysis.  Gas histories remain compact atmospheric histories, but are
    adjusted by broad latitude band to avoid treating southern/tropical
    recharge as identical to the Northern Hemisphere reference curve.
    """
    lat = context.latitude
    region = "northern_continental"
    if lat is not None and lat < 0:
        region = "southern_hemisphere"
    elif lat is not None and abs(lat) <= 23.5:
        region = "tropical"

    tritium_history: InputHistory | None = None
    mode = os.environ.get("HYDROSHEAF_TRITIUM_HISTORY_MODE", "usgs_quadrangle").strip().lower()
    if mode != "wiser":
        tritium_history, _ = _usgs_quadrangle_history(context)
    if tritium_history is None:
        try:
            from .wiser_loader import WISER_NA

            if WISER_NA is not None and context.latitude is not None and context.longitude is not None:
                tritium_history = WISER_NA.get_nearest_input_history(
                    context.latitude,
                    context.longitude,
                    "Tritium",
                )
                tritium_history = extend_history_to_present(tritium_history, region)
                tritium_history = _compact_history(tritium_history)
            else:
                raise ValueError("nearest WISER lookup requires bundled data and coordinates")
        except Exception:
            tritium_history = build_default_tritium_input(region)

    histories: dict[str, InputHistory] = {
        "3H": tritium_history,
    }

    for tracer in ("SF6", "CFC11", "CFC12", "CFC113", "85Kr"):
        histories[tracer] = _localized_gas_history(tracer, context)

    return histories


def site_input_history_metadata(context: SiteInputContext) -> dict[str, str]:
    """Return metadata tags describing which histories were selected."""
    lat = context.latitude
    region = _history_region(context)
    mode = os.environ.get("HYDROSHEAF_TRITIUM_HISTORY_MODE", "usgs_quadrangle").strip().lower()
    _, source = _usgs_quadrangle_history(context) if mode != "wiser" else (None, "")
    if source:
        return {
            "input_history_mode": "usgs_quadrangle",
            "input_history_region": source,
            "input_history_source": "USGS Table 2 quadrangle histories (Michel et al. 2018; Jurgens 2023 extension)",
        }
    try:
        from .wiser_loader import WISER_NA

        has_wiser = WISER_NA is not None and context.latitude is not None and context.longitude is not None
    except Exception:
        has_wiser = False
    return {
        "input_history_mode": "wiser_nearest_plus_site_adjusted_gases" if has_wiser else "hemisphere_scaled",
        "input_history_region": region,
        "input_history_source": "wiser_north_america_nearest" if has_wiser else f"default_{region}_fallback",
    }


def _history_region(context: SiteInputContext) -> str:
    lat = context.latitude
    if lat is not None and lat < 0:
        return "southern_hemisphere"
    if lat is not None and abs(lat) <= 23.5:
        return "tropical"
    return "northern_hemisphere"


def _localized_gas_history(tracer: str, context: SiteInputContext) -> InputHistory:
    from .multi_tracer import build_atmospheric_tracer_input

    base = build_atmospheric_tracer_input(tracer)
    region = _history_region(context)
    scale = 1.0
    lag_years = 0.0
    if region == "southern_hemisphere":
        scale = 0.94
        lag_years = 1.25
    elif region == "tropical":
        scale = 0.97
        lag_years = 0.5

    years = np.asarray(base.years, dtype=float) + lag_years
    values = np.asarray(base.values, dtype=float) * scale
    sigma = np.maximum(np.asarray(base.sigma, dtype=float) * scale, np.abs(values) * 0.03)
    return InputHistory(years, values, sigma)


def extend_history_to_present(
    history: InputHistory,
    region: str,
    *,
    present_year: float = 2026.0,
    step_years: float = 0.5,
) -> InputHistory:
    """Continue a station record past its last observation.

    ``InputHistory.interpolate`` delegates to :func:`numpy.interp`, which clamps
    out-of-range targets to the nearest endpoint value.  For GNIP/WISER tritium
    stations this is a serious defect: most North American stations stopped
    reporting between 1965 and 1999, while the benchmark samples were collected
    2004-2020.  Without a continuation, every recharge year after the record end
    is assigned the station's *last observed* value -- for a station whose record
    ends in 1965 that is a bomb-era 245 TU standing in for a modern background of
    a few TU.  Any piston-flow fit with a near-zero mean age then predicts an
    impossibly high tritium concentration.

    The continuation follows the regional post-bomb decline curve, rescaled so
    that it joins the station record continuously at the splice year.  This keeps
    the extension consistent with the framework's own regional input model rather
    than inventing a separate decay law, and it converges to the modern
    background instead of freezing a bomb-era value.
    """
    years = np.asarray(history.years, dtype=float)
    values = np.asarray(history.values, dtype=float)
    sigma = np.asarray(history.sigma, dtype=float)
    if years.size == 0:
        return history

    splice_year = float(years.max())
    if splice_year >= present_year:
        return history

    reference = build_default_tritium_input(region)
    ref_years = np.asarray(reference.years, dtype=float)
    ref_values = np.asarray(reference.values, dtype=float)

    ref_at_splice = float(np.interp(splice_year, ref_years, ref_values))
    station_at_splice = float(values[np.argmax(years)])
    # Rescale the regional curve onto the station level at the splice point.
    # Guard against a degenerate reference value and against extreme rescaling
    # driven by a single noisy final observation.
    scale = station_at_splice / ref_at_splice if ref_at_splice > 0.0 else 1.0
    scale = float(min(max(scale, 0.1), 10.0))

    tail_years = np.arange(splice_year + step_years, present_year + step_years, step_years)
    if tail_years.size == 0:
        return history
    tail_values = np.maximum(np.interp(tail_years, ref_years, ref_values) * scale, 0.0)
    # Past the bomb peak the continuation must not rise above the level already
    # reached at the splice point.
    if splice_year >= 1964.0:
        tail_values = np.minimum.accumulate(tail_values)
    tail_sigma = np.maximum(np.abs(tail_values) * 0.25, 0.5)

    return InputHistory(
        np.concatenate([years, tail_years]),
        np.concatenate([values, tail_values]),
        np.concatenate([sigma, tail_sigma]),
    )


def _compact_history(history: InputHistory, *, step_years: float = 0.5, max_points: int = 220) -> InputHistory:
    """Aggregate dense station histories to a stable grid for repeated fitting."""
    years = np.asarray(history.years, dtype=float)
    values = np.asarray(history.values, dtype=float)
    sigma = np.asarray(history.sigma, dtype=float)
    if len(years) <= max_points:
        return history

    bins = np.round(years / step_years) * step_years
    unique_bins, inverse = np.unique(bins, return_inverse=True)
    out_values = np.zeros_like(unique_bins, dtype=float)
    out_sigma = np.zeros_like(unique_bins, dtype=float)
    for idx in range(len(unique_bins)):
        mask = inverse == idx
        out_values[idx] = float(np.mean(values[mask]))
        out_sigma[idx] = float(np.sqrt(np.mean(sigma[mask] ** 2))) if np.any(mask) else 0.0
    return InputHistory(unique_bins, out_values, out_sigma)


GAS_TRACERS = ("SF6", "CFC11", "CFC12", "CFC113", "85KR")
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
    "39AR": "39Ar",
    "AR39": "39Ar",
    "ARGON39": "39Ar",
    "85KR": "85Kr",
    "KR85": "85Kr",
    "KRYPTON85": "85Kr",
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
