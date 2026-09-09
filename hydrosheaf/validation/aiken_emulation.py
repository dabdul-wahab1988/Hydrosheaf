"""Secondary, model-conditioned age/transport emulation for the Aiken release.

This module deliberately sits beside the Aiken parser rather than changing
HydroSheaf's primary benchmark metrics.  It turns the release's validated
MODPATH 5 particle records and CFC apparent-age intervals into auditable
*hypotheses*:

* endpoint records provide direct release-cell to final-cell transport
  hypotheses for each particle;
* pathline records provide indirect segment/path summaries;
* CFC values remain broad apparent-age screening intervals; and
* all results retain model direction, time-unit conversion, censoring, and
  source-member provenance.

The release does not contain independent well-to-well connectivity labels.
Consequently this module never emits a direct well--well edge truth table and
never authorises an integrated field-accuracy score.  A well-to-model-cell
crosswalk is considered explicit only when the exact run/well identity is
present in the MODPATH ``.loc`` starting-location member.  Coordinates are not
used for nearest-neighbour matching and no missing crosswalk is guessed.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

import pandas as pd

from .aiken_reference import (
    REFERENCE_TYPE,
    decode_modpath5_ipcode,
    parse_apparent_age_label,
)


EMULATION_SCHEMA = "aiken-model-conditioned-age-transport-emulation-v1"
EMULATION_TYPE = "secondary_model_conditioned_age_transport"
MODEL_TIME_REFERENCE_YEAR = 365.25
CLAIM_BOUNDARY = (
    "Aiken outputs are a calibrated MODFLOW-NWT/MODPATH5 model-conditioned "
    "reference. Endpoint and pathline results support implementation, "
    "direction, travel-time, and CFC-interval concordance diagnostics; they "
    "do not establish independent field age, flow, reaction, or direct "
    "well-to-well adjacency truth. Integrated field scoring is prohibited."
)


_RUN_FROM_MEMBER_RE = re.compile(r"(?:^|/)(?:output\.)?([^/]+_MP)(?:/|$)", re.IGNORECASE)
_TIME_FACTORS_TO_DAYS: dict[str, float] = {
    "day": 1.0,
    "days": 1.0,
    "d": 1.0,
    "hour": 1.0 / 24.0,
    "hours": 1.0 / 24.0,
    "h": 1.0 / 24.0,
    "minute": 1.0 / 1440.0,
    "minutes": 1.0 / 1440.0,
    "min": 1.0 / 1440.0,
    "second": 1.0 / 86400.0,
    "seconds": 1.0 / 86400.0,
    "s": 1.0 / 86400.0,
    "year": MODEL_TIME_REFERENCE_YEAR,
    "years": MODEL_TIME_REFERENCE_YEAR,
    "yr": MODEL_TIME_REFERENCE_YEAR,
    "yrs": MODEL_TIME_REFERENCE_YEAR,
}


def _empty(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=columns)


def _finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _text(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    result = str(value).strip()
    return result or None


def _normalise_time_unit(value: Any) -> str | None:
    text = _text(value)
    if text is None:
        return None
    return re.sub(r"[^a-z]", "", text.lower())


def time_unit_conversion_to_years(units: Any) -> dict[str, Any]:
    """Return a declared model-time conversion, failing closed when unknown.

    The Aiken model reference declares ``days``.  The conversion uses the
    fixed Julian-year denominator 365.25 and is recorded in every derived
    table/manifest so that model time is never confused with calendar age.
    """

    normalised = _normalise_time_unit(units)
    if normalised is None or normalised not in _TIME_FACTORS_TO_DAYS:
        return {
            "source_units": _text(units),
            "normalised_units": normalised,
            "days_per_source_unit": None,
            "years_per_source_unit": None,
            "conversion_status": "ABSTAIN_UNKNOWN_TIME_UNIT",
            "conversion_rule": None,
        }
    days_per_unit = _TIME_FACTORS_TO_DAYS[normalised]
    years_per_unit = days_per_unit / MODEL_TIME_REFERENCE_YEAR
    if normalised in {"year", "years", "yr", "yrs"}:
        rule = "source years retained as years"
    elif normalised in {"day", "days", "d"}:
        rule = "source days / 365.25"
    else:
        rule = f"source {normalised} * {days_per_unit:g} days / 365.25"
    return {
        "source_units": _text(units),
        "normalised_units": normalised,
        "days_per_source_unit": days_per_unit,
        "years_per_source_unit": years_per_unit,
        "conversion_status": "VALIDATED",
        "conversion_rule": rule,
    }


def _run_id_from_member(value: Any) -> str | None:
    text = _text(value)
    if text is None:
        return None
    normalised = text.replace("\\", "/")
    match = _RUN_FROM_MEMBER_RE.search(normalised)
    if match:
        return match.group(1)
    parts = normalised.split("/")
    for part in reversed(parts[:-1]):
        if part.lower().endswith("_mp"):
            return part
    return parts[-2] if len(parts) > 1 else None


def _well_id_from_member(value: Any) -> str | None:
    text = _text(value)
    if text is None:
        return None
    match = re.search(r"\b(?:AK|LEX)-\d+\b", Path(text.replace("\\", "/")).stem, re.IGNORECASE)
    return match.group(0).upper() if match else None


def _frame(value: Any, columns: list[str] | None = None) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        result = value.copy()
    elif value is None:
        result = pd.DataFrame()
    elif isinstance(value, Mapping):
        result = pd.DataFrame([dict(value)])
    else:
        result = pd.DataFrame(value)
    if columns:
        for column in columns:
            if column not in result:
                result[column] = pd.NA
    return result


def _reference_frames(reference: Any) -> dict[str, pd.DataFrame]:
    """Accept an ``AikenReference`` or a mapping of canonical frames."""

    names = (
        "well_metadata",
        "model_locations",
        "model_references",
        "endpoints",
        "pathlines",
        "cfc_ages",
    )
    if isinstance(reference, Mapping):
        return {name: _frame(reference.get(name)) for name in names}
    return {name: _frame(getattr(reference, name, None)) for name in names}


CROSSWALK_COLUMNS = [
    "run_id",
    "well_id",
    "node_id",
    "model_row",
    "model_column",
    "model_layer",
    "local_x",
    "local_y",
    "local_z",
    "location_ordinal",
    "source_member",
    "reference_member",
    "model_time_units",
    "well_identity_status",
    "model_location_status",
    "crosswalk_status",
    "crosswalk_basis",
    "crosswalk_is_explicit",
    "independent_well_to_model_truth",
    "reference_type",
]


def build_well_model_crosswalk(
    well_metadata: pd.DataFrame,
    model_locations: pd.DataFrame,
    model_references: pd.DataFrame | None = None,
    endpoints: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build an exact-ID, model-starting-location crosswalk.

    Rows from ``.loc`` are retained one-for-one because a run can release
    particles in more than one layer.  A missing ``.loc`` row is represented
    as ``NOT_ESTABLISHED_ENDPOINT_ONLY`` when an endpoint run exists; this is
    useful for auditability but is never promoted to an explicit crosswalk.
    """

    wells = _frame(well_metadata)
    locations = _frame(model_locations)
    references = _frame(model_references)
    endpoint_frame = _frame(endpoints)
    well_counts = (
        wells.groupby("well_id", dropna=False).size().to_dict()
        if "well_id" in wells
        else {}
    )
    well_nodes = {}
    if "well_id" in wells:
        for row in wells.itertuples(index=False):
            well = _text(getattr(row, "well_id", None))
            if well and well not in well_nodes:
                well_nodes[well] = _text(getattr(row, "node_id", None))
    reference_units: dict[str, Any] = {}
    reference_members: dict[str, Any] = {}
    if not references.empty:
        for row in references.itertuples(index=False):
            run = _text(getattr(row, "simulation_id", None))
            if run:
                reference_units[run] = _text(getattr(row, "time_units", None))
                reference_members[run] = _text(getattr(row, "source_member", None))

    records: list[dict[str, Any]] = []
    seen_runs: set[str] = set()
    if not locations.empty:
        for ordinal, row in enumerate(locations.to_dict(orient="records"), start=1):
            member = _text(row.get("source_member"))
            run = _text(row.get("run_id")) or _run_id_from_member(member)
            well = _text(row.get("well_id")) or _well_id_from_member(member)
            if run:
                seen_runs.add(run)
            well_count = int(well_counts.get(well, 0)) if well else 0
            well_status = (
                "EXACT_UNIQUE_WELL_ID" if well_count == 1 else
                ("AMBIGUOUS_WELL_ID" if well_count > 1 else "WELL_ID_NOT_IN_METADATA")
            )
            cell_values = [row.get(name) for name in ("model_row", "model_column", "model_layer")]
            has_cell = all(_finite(value) is not None for value in cell_values)
            location_status = "EXPLICIT_LOC_STARTING_CELL" if has_cell else "LOC_CELL_INCOMPLETE"
            crosswalk_status = (
                "EXPLICIT_MODEL_STARTING_CELL" if has_cell and well_count == 1 and run else
                ("MODEL_STARTING_CELL_WITH_ID_AMBIGUITY" if has_cell else "NOT_ESTABLISHED")
            )
            records.append(
                {
                    "run_id": run,
                    "well_id": well,
                    "node_id": well_nodes.get(well),
                    "model_row": int(_finite(row.get("model_row"))) if _finite(row.get("model_row")) is not None else None,
                    "model_column": int(_finite(row.get("model_column"))) if _finite(row.get("model_column")) is not None else None,
                    "model_layer": int(_finite(row.get("model_layer"))) if _finite(row.get("model_layer")) is not None else None,
                    "local_x": _finite(row.get("local_x")),
                    "local_y": _finite(row.get("local_y")),
                    "local_z": _finite(row.get("local_z")),
                    "location_ordinal": ordinal,
                    "source_member": member,
                    "reference_member": reference_members.get(run),
                    "model_time_units": reference_units.get(run),
                    "well_identity_status": well_status,
                    "model_location_status": location_status,
                    "crosswalk_status": crosswalk_status,
                    "crosswalk_basis": "exact_run_well_id_plus_MODPATH5_loc_cell",
                    "crosswalk_is_explicit": crosswalk_status == "EXPLICIT_MODEL_STARTING_CELL",
                    "independent_well_to_model_truth": False,
                    "reference_type": REFERENCE_TYPE,
                }
            )

    endpoint_runs: list[tuple[str | None, str | None]] = []
    if not endpoint_frame.empty:
        for row in endpoint_frame.to_dict(orient="records"):
            run = _text(row.get("run_id")) or _run_id_from_member(row.get("source_member"))
            well = _text(row.get("well_id")) or _well_id_from_member(row.get("source_member"))
            if run and run not in seen_runs:
                endpoint_runs.append((run, well))
                seen_runs.add(run)
    for run, well in endpoint_runs:
        well_count = int(well_counts.get(well, 0)) if well else 0
        records.append(
            {
                "run_id": run,
                "well_id": well,
                "node_id": well_nodes.get(well),
                "model_row": None,
                "model_column": None,
                "model_layer": None,
                "local_x": None,
                "local_y": None,
                "local_z": None,
                "location_ordinal": None,
                "source_member": None,
                "reference_member": reference_members.get(run),
                "model_time_units": reference_units.get(run),
                "well_identity_status": (
                    "EXACT_UNIQUE_WELL_ID" if well_count == 1 else
                    ("AMBIGUOUS_WELL_ID" if well_count > 1 else "WELL_ID_NOT_IN_METADATA")
                ),
                "model_location_status": "NO_LOC_STARTING_MEMBER",
                "crosswalk_status": "NOT_ESTABLISHED_ENDPOINT_ONLY",
                "crosswalk_basis": "endpoint_run_without_explicit_LOC_crosswalk",
                "crosswalk_is_explicit": False,
                "independent_well_to_model_truth": False,
                "reference_type": REFERENCE_TYPE,
            }
        )
    if not records:
        return _empty(CROSSWALK_COLUMNS)
    return pd.DataFrame(records, columns=CROSSWALK_COLUMNS)


CFC_INTERVAL_COLUMNS = [
    "well_id",
    "node_id",
    "sample_date",
    "sample_year",
    "apparent_age_label",
    "recharge_year_min",
    "recharge_year_max",
    "age_low_years",
    "age_high_years",
    "age_midpoint_years",
    "age_interval_width_years",
    "age_interval_status",
    "age_bounds_basis",
    "independent_age_truth",
    "reference_type",
    "source_workbook",
    "source_sheet",
    "source_row",
]


def build_cfc_age_intervals(
    cfc_ages: pd.DataFrame,
    *,
    reference_year: int = 2015,
) -> pd.DataFrame:
    """Convert published qualitative CFC recharge labels to age intervals.

    If a sample date is absent, the release's declared 2015 comparison year
    is used only as a labelled reference-year fallback.  No exact age or
    independent age truth is manufactured.
    """

    ages = _frame(cfc_ages)
    if ages.empty:
        return _empty(CFC_INTERVAL_COLUMNS)
    records: list[dict[str, Any]] = []
    for row in ages.to_dict(orient="records"):
        label = row.get("apparent_age_label")
        parsed = {
            "recharge_year_min": row.get("recharge_year_min"),
            "recharge_year_max": row.get("recharge_year_max"),
        }
        if _finite(parsed["recharge_year_min"]) is None or _finite(parsed["recharge_year_max"]) is None:
            parsed_label = parse_apparent_age_label(label)
            parsed["recharge_year_min"] = parsed_label.get("recharge_year_min")
            parsed["recharge_year_max"] = parsed_label.get("recharge_year_max")
        lower_recharge = _finite(parsed["recharge_year_min"])
        upper_recharge = _finite(parsed["recharge_year_max"])
        date_value = pd.to_datetime(row.get("sample_date"), errors="coerce")
        sample_year = int(date_value.year) if not pd.isna(date_value) else int(reference_year)
        date_basis = "sample_date_year" if not pd.isna(date_value) else "release_reference_year_fallback_2015"
        if lower_recharge is None or upper_recharge is None:
            low_age = high_age = midpoint = width = None
            status = "ABSTAIN_AGE_INTERVAL_UNAVAILABLE"
            basis = "no_published_recharge_year_bounds"
        else:
            low_age = max(0.0, float(sample_year) - upper_recharge)
            high_age = max(low_age, float(sample_year) - lower_recharge)
            midpoint = (low_age + high_age) / 2.0
            width = high_age - low_age
            status = "CFC_APPARENT_AGE_SCREENING_INTERVAL"
            basis = f"{date_basis}_minus_published_recharge_year_interval"
        records.append(
            {
                "well_id": _text(row.get("well_id")),
                "node_id": _text(row.get("node_id")),
                "sample_date": row.get("sample_date"),
                "sample_year": sample_year,
                "apparent_age_label": _text(label),
                "recharge_year_min": lower_recharge,
                "recharge_year_max": upper_recharge,
                "age_low_years": low_age,
                "age_high_years": high_age,
                "age_midpoint_years": midpoint,
                "age_interval_width_years": width,
                "age_interval_status": status,
                "age_bounds_basis": basis,
                "independent_age_truth": False,
                "reference_type": REFERENCE_TYPE,
                "source_workbook": row.get("source_workbook"),
                "source_sheet": row.get("source_sheet"),
                "source_row": row.get("source_row"),
            }
        )
    return pd.DataFrame(records, columns=CFC_INTERVAL_COLUMNS)


DIRECT_COLUMNS = [
    "run_id",
    "well_id",
    "particle_ordinal",
    "particle_id",
    "particle_id_status",
    "release_node",
    "final_node",
    "release_x",
    "release_y",
    "release_z_local",
    "final_x",
    "final_y",
    "final_z_local",
    "total_tracking_time_raw",
    "release_time_raw",
    "travel_time_model_units",
    "model_time_units",
    "travel_time_days",
    "travel_time_years",
    "travel_time_observation_status",
    "event_observed",
    "censoring_status",
    "termination_status",
    "ipcode",
    "ipcode_idcode",
    "ipcode_nslast",
    "tracking_direction_code",
    "tracking_direction",
    "direction_status",
    "hypothesis_kind",
    "well_to_well_edge_status",
    "direct_adjacency_truth_status",
    "reference_type",
    "source_archive",
    "source_member",
    "reference_time",
    "crosswalk_status",
]

SEGMENT_COLUMNS = [
    "run_id",
    "well_id",
    "particle_ordinal",
    "particle_id",
    "segment_ordinal",
    "from_model_node",
    "to_model_node",
    "time_start_raw",
    "time_end_raw",
    "segment_time_model_units",
    "model_time_units",
    "segment_time_days",
    "segment_time_years",
    "segment_observation_status",
    "endpoint_event_observed",
    "endpoint_censoring_status",
    "endpoint_termination_status",
    "endpoint_ipcode",
    "tracking_direction_code",
    "tracking_direction",
    "direction_status",
    "hypothesis_kind",
    "direct_adjacency_truth_status",
    "reference_type",
    "source_archive",
    "source_member",
]

SUMMARY_COLUMNS = [
    "run_id",
    "well_id",
    "hypothesis_kind",
    "model_time_units",
    "time_conversion_status",
    "time_conversion_rule",
    "n_particle_records",
    "n_observed_termination",
    "n_right_censored",
    "n_non_target_termination",
    "n_usable_travel_times",
    "travel_time_years_median",
    "travel_time_years_p10",
    "travel_time_years_p90",
    "travel_time_years_min",
    "travel_time_years_max",
    "direction_values",
    "crosswalk_status",
    "cfc_age_match_status",
    "cfc_age_low_years",
    "cfc_age_high_years",
    "n_age_concordant_observed",
    "n_age_outside_observed",
    "integrated_field_scoring_allowed",
    "direct_adjacency_truth_status",
    "reference_type",
]


def _prepare_model_references(model_references: pd.DataFrame | None) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    frame = _frame(model_references)
    if frame.empty:
        return result
    for row in frame.to_dict(orient="records"):
        run = _text(row.get("simulation_id")) or _text(row.get("run_id"))
        if not run:
            continue
        conversion = time_unit_conversion_to_years(row.get("time_units"))
        result[run] = {
            **conversion,
            "model_time_units": _text(row.get("time_units")),
            "reference_member": _text(row.get("source_member")),
        }
    return result


def _endpoint_records(
    endpoints: pd.DataFrame,
    model_references: Mapping[str, Mapping[str, Any]],
    crosswalk: pd.DataFrame,
) -> pd.DataFrame:
    frame = _frame(endpoints)
    if frame.empty:
        return _empty(DIRECT_COLUMNS)
    crosswalk_by_run = {}
    if not crosswalk.empty:
        for row in crosswalk.to_dict(orient="records"):
            run = _text(row.get("run_id"))
            if run and run not in crosswalk_by_run:
                crosswalk_by_run[run] = _text(row.get("crosswalk_status"))
    records: list[dict[str, Any]] = []
    for ordinal, row in enumerate(frame.to_dict(orient="records"), start=1):
        run = _text(row.get("run_id")) or _run_id_from_member(row.get("source_member"))
        well = _text(row.get("well_id")) or _well_id_from_member(row.get("source_member"))
        particle_ordinal = _finite(row.get("particle_ordinal"))
        if particle_ordinal is None:
            particle_ordinal = _finite(row.get("record_ordinal")) or ordinal
        particle_ordinal = int(particle_ordinal)
        time_total = _finite(row.get("total_tracking_time"))
        if time_total is None:
            time_total = _finite(row.get("time"))
        release_time = _finite(row.get("release_time")) or 0.0
        duration = time_total - release_time if time_total is not None else None
        model_info = dict(model_references.get(run or "", {}))
        if not model_info:
            model_info = time_unit_conversion_to_years(row.get("model_time_units") or row.get("time_units"))
            model_info["model_time_units"] = _text(row.get("model_time_units") or row.get("time_units"))
        scale = _finite(model_info.get("years_per_source_unit"))
        travel_years = duration * scale if duration is not None and scale is not None and duration >= 0 else None
        travel_days = duration * _finite(model_info.get("days_per_source_unit")) if duration is not None and _finite(model_info.get("days_per_source_unit")) is not None and duration >= 0 else None
        ipcode = _finite(row.get("ipcode"))
        decoded_ipcode = decode_modpath5_ipcode(int(ipcode)) if ipcode is not None else {}
        censoring = _text(row.get("particle_censoring_status")) or decoded_ipcode.get("particle_censoring_status") or "unknown"
        termination = _text(row.get("termination_status")) or decoded_ipcode.get("termination_status") or "unknown_ipcode"
        if censoring == "observed_termination":
            event_observed: bool | None = True
            observation_status = "exact_model_termination_time" if travel_years is not None else "ABSTAIN_TIME_INVALID"
        elif censoring == "right_censored_active":
            event_observed = False
            observation_status = "right_censored_lower_bound" if travel_years is not None else "ABSTAIN_TIME_INVALID"
        elif censoring in {"censored_dry_cell", "unreleased"}:
            event_observed = False
            observation_status = "non_target_termination"
        else:
            event_observed = None
            observation_status = "ABSTAIN_CENSORING_UNKNOWN"
        ipcode_fields = decoded_ipcode if ipcode is not None else {
            "ipcode_idcode": row.get("ipcode_idcode"),
            "ipcode_nslast": row.get("ipcode_nslast"),
            "termination_status": termination,
            "particle_censoring_status": censoring,
        }
        direction = _text(row.get("tracking_direction")) or "unknown_tracking_direction"
        direction_status = "DECLARED_FROM_RUN_CONFIGURATION" if not direction.startswith("unknown") else "ABSTAIN_DIRECTION_UNKNOWN"
        records.append(
            {
                "run_id": run,
                "well_id": well,
                "particle_ordinal": particle_ordinal,
                "particle_id": row.get("particle_id"),
                "particle_id_status": _text(row.get("particle_id_status")),
                "release_node": row.get("release_node") if row.get("release_node") is not None else row.get("initial_cell"),
                "final_node": row.get("final_node") if row.get("final_node") is not None else row.get("final_cell"),
                "release_x": row.get("release_x") if row.get("release_x") is not None else row.get("x0"),
                "release_y": row.get("release_y") if row.get("release_y") is not None else row.get("y0"),
                "release_z_local": row.get("release_z_local") if row.get("release_z_local") is not None else row.get("zloc0"),
                "final_x": row.get("final_x") if row.get("final_x") is not None else row.get("x"),
                "final_y": row.get("final_y") if row.get("final_y") is not None else row.get("y"),
                "final_z_local": row.get("final_z_local") if row.get("final_z_local") is not None else row.get("zloc"),
                "total_tracking_time_raw": time_total,
                "release_time_raw": release_time,
                "travel_time_model_units": duration,
                "model_time_units": model_info.get("model_time_units"),
                "travel_time_days": travel_days,
                "travel_time_years": travel_years,
                "travel_time_observation_status": observation_status,
                "event_observed": event_observed,
                "censoring_status": censoring,
                "termination_status": ipcode_fields.get("termination_status", termination),
                "ipcode": int(ipcode) if ipcode is not None else row.get("ipcode"),
                "ipcode_idcode": ipcode_fields.get("ipcode_idcode"),
                "ipcode_nslast": ipcode_fields.get("ipcode_nslast"),
                "tracking_direction_code": row.get("tracking_direction_code"),
                "tracking_direction": direction,
                "direction_status": direction_status,
                "hypothesis_kind": "direct_endpoint_transport",
                "well_to_well_edge_status": "PROHIBITED_NO_WELL_TO_WELL_TRUTH",
                "direct_adjacency_truth_status": "ABSTAIN",
                "reference_type": REFERENCE_TYPE,
                "source_archive": row.get("source_archive"),
                "source_member": row.get("source_member"),
                "reference_time": row.get("reference_time"),
                "crosswalk_status": crosswalk_by_run.get(run, "NOT_ESTABLISHED"),
            }
        )
    return pd.DataFrame(records, columns=DIRECT_COLUMNS)


def _pathline_segments(
    pathlines: pd.DataFrame,
    model_references: Mapping[str, Mapping[str, Any]],
    crosswalk: pd.DataFrame,
) -> pd.DataFrame:
    frame = _frame(pathlines)
    if frame.empty:
        return _empty(SEGMENT_COLUMNS)
    rows: list[dict[str, Any]] = []
    crosswalk_by_run = {}
    if not crosswalk.empty:
        for record in crosswalk.to_dict(orient="records"):
            run = _text(record.get("run_id"))
            if run and run not in crosswalk_by_run:
                crosswalk_by_run[run] = _text(record.get("crosswalk_status"))
    work = frame.copy()
    if "run_id" not in work:
        work["run_id"] = pd.NA
    if "particle_ordinal" not in work:
        if "particle_id" in work:
            work["particle_ordinal"] = work.groupby("run_id", dropna=False)["particle_id"].transform(
                lambda values: pd.factorize(values, sort=False)[0] + 1
            )
        else:
            work["particle_ordinal"] = pd.NA
    if "record_ordinal" not in work:
        work["record_ordinal"] = range(1, len(work) + 1)
    for (run, particle_ordinal), group in work.groupby(["run_id", "particle_ordinal"], dropna=False, sort=False):
        group = group.sort_values("record_ordinal", kind="stable")
        records = group.to_dict(orient="records")
        if len(records) < 2:
            continue
        run_text = _text(run)
        model_info = dict(model_references.get(run_text or "", {}))
        if not model_info:
            units = records[0].get("model_time_units") or records[0].get("time_units")
            model_info = time_unit_conversion_to_years(units)
            model_info["model_time_units"] = _text(units)
        scale = _finite(model_info.get("years_per_source_unit"))
        days_scale = _finite(model_info.get("days_per_source_unit"))
        for index, (left, right) in enumerate(zip(records[:-1], records[1:]), start=1):
            start = _finite(left.get("tracking_time"))
            end = _finite(right.get("tracking_time"))
            if start is None:
                start = abs(_finite(left.get("time")) or 0.0)
            if end is None:
                end = abs(_finite(right.get("time")) or 0.0)
            duration = end - start
            valid = duration >= 0
            rows.append(
                {
                    "run_id": run_text,
                    "well_id": _text(left.get("well_id")) or _well_id_from_member(left.get("source_member")),
                    "particle_ordinal": int(_finite(particle_ordinal)) if _finite(particle_ordinal) is not None else None,
                    "particle_id": left.get("particle_id"),
                    "segment_ordinal": index,
                    "from_model_node": left.get("node") if left.get("node") is not None else left.get("global_node"),
                    "to_model_node": right.get("node") if right.get("node") is not None else right.get("global_node"),
                    "time_start_raw": left.get("time_raw") if left.get("time_raw") is not None else left.get("time"),
                    "time_end_raw": right.get("time_raw") if right.get("time_raw") is not None else right.get("time"),
                    "segment_time_model_units": duration if valid else None,
                    "model_time_units": model_info.get("model_time_units"),
                    "segment_time_days": duration * days_scale if valid and days_scale is not None else None,
                    "segment_time_years": duration * scale if valid and scale is not None else None,
                    "segment_observation_status": "validated_segment_increment" if valid and scale is not None else (
                        "ABSTAIN_UNKNOWN_TIME_UNIT" if scale is None else "ABSTAIN_NONMONOTONIC_TIME"
                    ),
                    "endpoint_event_observed": (
                        True if _text(left.get("particle_censoring_status")) == "observed_termination" else
                        (False if _text(left.get("particle_censoring_status")) in {
                            "right_censored_active", "censored_dry_cell", "unreleased"
                        } else None)
                    ),
                    "endpoint_censoring_status": _text(left.get("particle_censoring_status")) or "unknown",
                    "endpoint_termination_status": _text(left.get("termination_status")) or "unknown_ipcode",
                    "endpoint_ipcode": left.get("ipcode"),
                    "tracking_direction_code": left.get("tracking_direction_code"),
                    "tracking_direction": _text(left.get("tracking_direction")) or "unknown_tracking_direction",
                    "direction_status": "DECLARED_FROM_RUN_CONFIGURATION" if _text(left.get("tracking_direction")) and not str(left.get("tracking_direction")).startswith("unknown") else "ABSTAIN_DIRECTION_UNKNOWN",
                    "hypothesis_kind": "indirect_pathline_segment_transport",
                    "direct_adjacency_truth_status": "ABSTAIN",
                    "reference_type": REFERENCE_TYPE,
                    "source_archive": left.get("source_archive"),
                    "source_member": left.get("source_member"),
                }
            )
    return pd.DataFrame(rows, columns=SEGMENT_COLUMNS)


def _percentile(values: pd.Series, quantile: float) -> float | None:
    if values.empty:
        return None
    result = values.quantile(quantile)
    return _finite(result)


def _attach_age_concordance(
    direct: pd.DataFrame,
    cfc_intervals: pd.DataFrame,
) -> dict[tuple[str | None, str | None], dict[str, Any]]:
    lookup: dict[tuple[str | None, str | None], dict[str, Any]] = {}
    if cfc_intervals.empty:
        return lookup
    counts = cfc_intervals.groupby("well_id", dropna=False).size().to_dict()
    for row in cfc_intervals.to_dict(orient="records"):
        key = (_text(row.get("run_id")), _text(row.get("well_id")))
        # CFC ages are usually keyed by well rather than run.  The fallback
        # below allows direct hypotheses to use the same interval without a
        # fabricated run identifier.
        key_well = (None, _text(row.get("well_id")))
        if int(counts.get(row.get("well_id"), 0)) != 1:
            continue
        value = {
            "status": _text(row.get("age_interval_status")),
            "low": _finite(row.get("age_low_years")),
            "high": _finite(row.get("age_high_years")),
        }
        lookup[key] = value
        lookup[key_well] = value
    return lookup


def build_transport_hypotheses(
    endpoints: pd.DataFrame,
    pathlines: pd.DataFrame,
    model_references: pd.DataFrame | None,
    crosswalk: pd.DataFrame,
    cfc_intervals: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build direct endpoint, indirect segment, and run-level summaries."""

    model_info = _prepare_model_references(model_references)
    direct = _endpoint_records(endpoints, model_info, crosswalk)
    segments = _pathline_segments(pathlines, model_info, crosswalk)
    age_lookup = _attach_age_concordance(direct, cfc_intervals)
    if not direct.empty:
        age_statuses: list[str] = []
        age_lows: list[float | None] = []
        age_highs: list[float | None] = []
        for row in direct.to_dict(orient="records"):
            age = age_lookup.get((_text(row.get("run_id")), _text(row.get("well_id")))) or age_lookup.get((None, _text(row.get("well_id"))))
            travel = _finite(row.get("travel_time_years"))
            if age is None or age.get("low") is None or age.get("high") is None:
                status = "ABSTAIN_NO_UNIQUE_CFC_INTERVAL"
                low = high = None
            elif row.get("event_observed") is True and travel is not None:
                low, high = age["low"], age["high"]
                status = "consistent_with_cfc_interval" if low <= travel <= high else "outside_cfc_interval"
            elif row.get("censoring_status") == "right_censored_active" and travel is not None:
                low, high = age["low"], age["high"]
                status = "censored_compatible_with_cfc_interval" if high > travel else "censored_inconsistent_with_cfc_interval"
            else:
                low, high = age["low"], age["high"]
                status = "ABSTAIN_NON_TARGET_OR_UNKNOWN_TIME"
            age_statuses.append(status)
            age_lows.append(low)
            age_highs.append(high)
        direct["cfc_age_concordance_status"] = age_statuses
        direct["cfc_age_low_years"] = age_lows
        direct["cfc_age_high_years"] = age_highs
    summary_rows: list[dict[str, Any]] = []
    direct_groups = []
    if not direct.empty:
        direct_groups.append(("direct_endpoint_transport", direct))
    if not segments.empty:
        segments_for_summary = segments.rename(columns={"segment_time_years": "travel_time_years"}).copy()
        segments_for_summary["event_observed"] = segments_for_summary["endpoint_event_observed"]
        segments_for_summary["censoring_status"] = "not_applicable_segment"
        segments_for_summary["cfc_age_concordance_status"] = "ABSTAIN_SEGMENT_NOT_ENDPOINT_AGE"
        direct_groups.append(("indirect_pathline_segment_transport", segments_for_summary))
    for kind, group in direct_groups:
        group = group.copy()
        # Keep a row per run/well; run-level transport summaries are not
        # pooled across differently configured simulations.
        group_keys = [column for column in ("run_id", "well_id") if column in group]
        grouped = group.groupby(group_keys, dropna=False, sort=False) if group_keys else [((), group)]
        for key, sub in grouped:
            if group_keys:
                key_values = key if isinstance(key, tuple) else (key,)
                run = _text(key_values[group_keys.index("run_id")]) if "run_id" in group_keys else None
                well = _text(key_values[group_keys.index("well_id")]) if "well_id" in group_keys else None
            else:
                run = well = None
            units = _text(sub["model_time_units"].dropna().iloc[0]) if "model_time_units" in sub and sub["model_time_units"].notna().any() else None
            conversion = time_unit_conversion_to_years(units)
            travel_values = pd.to_numeric(sub.get("travel_time_years", pd.Series(dtype=float)), errors="coerce").dropna()
            if kind == "direct_endpoint_transport":
                observed_count = int(sub.get("event_observed", pd.Series(dtype=bool)).fillna(False).sum())
                right_count = int((sub.get("censoring_status") == "right_censored_active").sum())
                non_target_count = int(sub.get("travel_time_observation_status", pd.Series(dtype=str)).eq("non_target_termination").sum())
                concordant = int(sub.get("cfc_age_concordance_status", pd.Series(dtype=str)).eq("consistent_with_cfc_interval").sum())
                outside = int(sub.get("cfc_age_concordance_status", pd.Series(dtype=str)).eq("outside_cfc_interval").sum())
                cfc_status = next((str(value) for value in sub.get("cfc_age_concordance_status", pd.Series(dtype=str)).tolist() if str(value) not in {"ABSTAIN_NO_UNIQUE_CFC_INTERVAL", "nan"}), "ABSTAIN_NO_UNIQUE_CFC_INTERVAL")
                cfc_low = _finite(sub.get("cfc_age_low_years", pd.Series(dtype=float)).dropna().iloc[0]) if "cfc_age_low_years" in sub and sub["cfc_age_low_years"].notna().any() else None
                cfc_high = _finite(sub.get("cfc_age_high_years", pd.Series(dtype=float)).dropna().iloc[0]) if "cfc_age_high_years" in sub and sub["cfc_age_high_years"].notna().any() else None
                crosswalk_status = _text(sub.get("crosswalk_status", pd.Series(dtype=str)).dropna().iloc[0]) if "crosswalk_status" in sub and sub["crosswalk_status"].notna().any() else "NOT_ESTABLISHED"
                direction_values = "|".join(sorted({str(value) for value in sub.get("tracking_direction", pd.Series(dtype=str)).dropna()}))
            else:
                unique_particle_sub = sub.drop_duplicates(
                    [column for column in ("run_id", "particle_ordinal") if column in sub]
                )
                event_observed = pd.Series(
                    unique_particle_sub.get(
                        "event_observed", pd.Series(dtype=bool)
                    ),
                    copy=False,
                ).astype("boolean")
                observed_count = int(event_observed.fillna(False).sum())
                right_count = int(
                    unique_particle_sub.get(
                        "endpoint_censoring_status", pd.Series(dtype=str)
                    ).eq("right_censored_active").sum()
                )
                non_target_count = int(
                    unique_particle_sub.get(
                        "endpoint_censoring_status", pd.Series(dtype=str)
                    ).isin({"censored_dry_cell", "unreleased"}).sum()
                )
                concordant = outside = 0
                cfc_status = "ABSTAIN_SEGMENT_NOT_ENDPOINT_AGE"
                cfc_low = cfc_high = None
                crosswalk_status = "NOT_ESTABLISHED"
                direction_values = "|".join(
                    sorted(
                        {
                            str(value)
                            for value in sub.get(
                                "tracking_direction", pd.Series(dtype=str)
                            ).dropna()
                        }
                    )
                )
            summary_rows.append(
                {
                    "run_id": run,
                    "well_id": well,
                    "hypothesis_kind": kind,
                    "model_time_units": units,
                    "time_conversion_status": conversion["conversion_status"],
                    "time_conversion_rule": conversion["conversion_rule"],
                    "n_particle_records": int(len(sub)),
                    "n_observed_termination": observed_count,
                    "n_right_censored": right_count,
                    "n_non_target_termination": non_target_count,
                    "n_usable_travel_times": int(len(travel_values)),
                    "travel_time_years_median": _finite(travel_values.median()) if not travel_values.empty else None,
                    "travel_time_years_p10": _percentile(travel_values, 0.10),
                    "travel_time_years_p90": _percentile(travel_values, 0.90),
                    "travel_time_years_min": _finite(travel_values.min()) if not travel_values.empty else None,
                    "travel_time_years_max": _finite(travel_values.max()) if not travel_values.empty else None,
                    "direction_values": direction_values,
                    "crosswalk_status": crosswalk_status,
                    "cfc_age_match_status": cfc_status,
                    "cfc_age_low_years": cfc_low,
                    "cfc_age_high_years": cfc_high,
                    "n_age_concordant_observed": concordant,
                    "n_age_outside_observed": outside,
                    "integrated_field_scoring_allowed": False,
                    "direct_adjacency_truth_status": "ABSTAIN",
                    "reference_type": REFERENCE_TYPE,
                }
            )
    summary = pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS)
    return direct, segments, summary


@dataclass(frozen=True)
class AikenEmulationResult:
    """Output tables and immutable manifest for one emulation run."""

    crosswalk: pd.DataFrame
    cfc_age_intervals: pd.DataFrame
    direct_hypotheses: pd.DataFrame
    indirect_segments: pd.DataFrame
    transport_summary: pd.DataFrame
    manifest: dict[str, Any]


def _manifest(
    *,
    emulation_id: str,
    crosswalk: pd.DataFrame,
    cfc_age_intervals: pd.DataFrame,
    direct: pd.DataFrame,
    segments: pd.DataFrame,
    summary: pd.DataFrame,
    source: Any = None,
) -> dict[str, Any]:
    explicit = int(crosswalk.get("crosswalk_is_explicit", pd.Series(dtype=bool)).fillna(False).sum()) if not crosswalk.empty else 0
    crosswalk_status = (
        "EXPLICIT" if explicit and explicit == len(crosswalk) else
        ("PARTIAL" if explicit else "NOT_ESTABLISHED")
    )
    units = sorted({str(value) for value in summary.get("model_time_units", pd.Series(dtype=str)).dropna()}) if not summary.empty else []
    conversions = sorted({str(value) for value in summary.get("time_conversion_status", pd.Series(dtype=str)).dropna()}) if not summary.empty else []
    source_path = str(getattr(source, "source", "")) if source is not None else None
    return {
        "schema": EMULATION_SCHEMA,
        "emulation_id": str(emulation_id),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "emulation_type": EMULATION_TYPE,
        "reference_type": REFERENCE_TYPE,
        "source": source_path,
        "source_doi": "10.5066/P9U0GHLU",
        "publication_doi": "10.3133/sir20225036",
        "crosswalk": {
            "status": crosswalk_status,
            "n_rows": int(len(crosswalk)),
            "n_explicit_rows": explicit,
            "basis": "exact run/well identity plus MODPATH5 .loc starting cell only",
            "nearest_coordinate_matching_used": False,
            "independent_well_to_model_truth": False,
        },
        "model_time": {
            "source_units_observed": units,
            "conversion_statuses": conversions,
            "year_denominator_days": MODEL_TIME_REFERENCE_YEAR,
            "rule": "declared model units converted to years only when the model reference declares a supported unit",
        },
        "outputs": {
            "crosswalk_rows": int(len(crosswalk)),
            "cfc_age_interval_rows": int(len(cfc_age_intervals)),
            "direct_endpoint_hypothesis_rows": int(len(direct)),
            "indirect_pathline_segment_rows": int(len(segments)),
            "transport_summary_rows": int(len(summary)),
        },
        "independence_and_claims": {
            "independent_direct_adjacency_truth": False,
            "independent_age_truth": False,
            "independent_reaction_truth": False,
            "direct_well_to_well_truth_emitted": False,
            "integrated_field_scoring_allowed": False,
            "integrated_scoring_allowed": False,
            "field_benchmark_status": "PROHIBITED",
        },
        "claim_boundary": CLAIM_BOUNDARY,
    }


def _sha256_file(path: Path) -> str:
    """Hash a generated artifact without reading the source archives."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_aiken_model_conditioned_emulation(
    reference: Any,
    *,
    emulation_id: str = "AIKEN-MODEL-CONDITIONED-AGE-TRANSPORT",
    reference_year: int = 2015,
) -> AikenEmulationResult:
    """Run the secondary emulation from an ``AikenReference``-like object."""

    frames = _reference_frames(reference)
    crosswalk = build_well_model_crosswalk(
        frames["well_metadata"],
        frames["model_locations"],
        frames["model_references"],
        frames["endpoints"],
    )
    cfc_intervals = build_cfc_age_intervals(frames["cfc_ages"], reference_year=reference_year)
    direct, segments, summary = build_transport_hypotheses(
        frames["endpoints"],
        frames["pathlines"],
        frames["model_references"],
        crosswalk,
        cfc_intervals,
    )
    manifest = _manifest(
        emulation_id=emulation_id,
        crosswalk=crosswalk,
        cfc_age_intervals=cfc_intervals,
        direct=direct,
        segments=segments,
        summary=summary,
        source=reference,
    )
    return AikenEmulationResult(
        crosswalk=crosswalk,
        cfc_age_intervals=cfc_intervals,
        direct_hypotheses=direct,
        indirect_segments=segments,
        transport_summary=summary,
        manifest=manifest,
    )


def write_aiken_emulation_outputs(
    result: AikenEmulationResult,
    output: Path | str,
) -> dict[str, str]:
    """Write emulation tables and a manifest without hashing source archives.

    The generated tables are hashed after serialization.  The manifest itself
    is deliberately excluded from ``output_files`` because including its own
    digest would create a circular value; the run directory and the table
    hashes are the immutable artifact boundary.
    """

    root = Path(output).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    tables = {
        "aiken_model_crosswalk": result.crosswalk,
        "aiken_cfc_age_intervals": result.cfc_age_intervals,
        "aiken_direct_endpoint_hypotheses": result.direct_hypotheses,
        "aiken_indirect_pathline_segments": result.indirect_segments,
        "aiken_transport_summary": result.transport_summary,
    }
    generated: dict[str, str] = {}
    for name, frame in tables.items():
        filename = f"{name}.csv"
        frame.to_csv(root / filename, index=False)
        generated[name] = filename
    result.manifest["output_files"] = {
        name: {
            "path": filename,
            "sha256": _sha256_file(root / filename),
            "size_bytes": int((root / filename).stat().st_size),
        }
        for name, filename in generated.items()
    }
    result.manifest["output_hash_scope"] = (
        "Generated CSV tables only; source archives and this manifest are not "
        "hashed by this writer."
    )
    manifest_name = "aiken_model_conditioned_emulation_manifest.json"
    (root / manifest_name).write_text(
        json.dumps(result.manifest, indent=2, ensure_ascii=False, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    generated["manifest"] = manifest_name
    return generated


__all__ = [
    "AikenEmulationResult",
    "CLAIM_BOUNDARY",
    "CFC_INTERVAL_COLUMNS",
    "CROSSWALK_COLUMNS",
    "DIRECT_COLUMNS",
    "EMULATION_SCHEMA",
    "EMULATION_TYPE",
    "MODEL_TIME_REFERENCE_YEAR",
    "SEGMENT_COLUMNS",
    "SUMMARY_COLUMNS",
    "build_cfc_age_intervals",
    "build_transport_hypotheses",
    "build_well_model_crosswalk",
    "run_aiken_model_conditioned_emulation",
    "time_unit_conversion_to_years",
    "write_aiken_emulation_outputs",
]
