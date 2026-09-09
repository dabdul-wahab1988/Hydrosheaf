"""Read-only inventory and harmonisation diagnostic for the field datasets.

This module inventories the four field-data products used by the groundwater
workflows.  It does not alter source files, infer undocumented units, merge
observations, or create a modelling input table.  Harmonisation is limited to
reporting conservative aliases for a small canonical vocabulary and marking
where a field's unit or meaning is not declared.

Run from the repository root with::

    python M7/m7_nonuniqueness_benchmark/scripts/field_data_inventory.py

Use ``--format json`` for machine-readable output.  The default Markdown
report is deliberately printed to stdout so a caller can redirect it to a
review artifact without this diagnostic writing into a locked results tree.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping
import json
from pathlib import Path
import re
import unicodedata
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DATASETS: dict[str, Path] = {
    "LowerAnayari/manu": Path("data/FieldData/LowerAnayari/manu.csv"),
    "NorthenGhana/NorthernGhana": Path(
        "data/FieldData/NorthenGhana/NorthernGhana.xlsx"
    ),
    "NorthernGhanaNew/compiled UER data_new": Path(
        "data/FieldData/NorthernGhanaNew/compiled UER data_new.xlsx"
    ),
    "Talensi_MiningArea/talensi": Path(
        "data/FieldData/Talensi_MiningArea/talensi.csv"
    ),
}

REQUESTED_FIELDS = ("SiO2", "Sr", "Fe", "depth", "date")

# Aliases are intentionally conservative.  A name that resembles a field is
# not enough to establish a unit; unit inference is handled separately.
CANONICAL_ALIASES: dict[str, tuple[str, ...]] = {
    "sample_id": (
        "sample id",
        "sample_id",
        "well_id",
        "well id",
        "code",
        # The UER workbook uses a sequential ``No.`` column.  It is reported
        # as an identifier candidate, not upgraded to a globally stable ID.
        "no",
        "number",
    ),
    "site_group": ("station", "community", "community_code", "town"),
    "latitude": ("latitude", "lat", "latituted", "y coordinate"),
    "longitude": ("longitude", "lon", "long", "x coordinate"),
    "elevation": ("elevation", "elevation_m", "elev", "elevation m"),
    "depth": (
        "borehole_depth_m",
        "borehole depth m",
        "well_depth",
        "well depth",
        "screen_depth",
        "screen depth",
        "depth",
    ),
    "date": ("date", "sample_date", "sample date", "date sampled"),
    "hydraulic_head": (
        "head_meas",
        "hydraulic_head",
        "hydraulic head",
        "water_table_elevation",
        "piezometric_head",
    ),
    "water_level": (
        "static_water_level_m",
        "static water level m",
        "static water level",
        "water_level",
        "water level",
    ),
    "SiO2": (
        "sio2_mg_l",
        "sio2",
        "si o2",
        "dissolved sio2",
        "silica",
    ),
    "Sr": ("sr_mg_l", "sr", "sr2+", "strontium"),
    "Fe": ("fe_mg_l", "fe", "fe2+", "iron"),
    "distance_river": (
        "distance_river_km",
        "distance river km",
        "distance to river",
    ),
    "rainfall": ("amount(mm)", "rainfall", "rainfall mm", "amount mm"),
}

UNIT_TOKENS = {
    "m",
    "mg/l",
    "mg/litre",
    "us/cm",
    "µs/cm",
    "μs/cm",
    "°c",
    "c",
    "‰",
    "tu",
    "mm",
}


def _normalise_name(value: Any) -> str:
    text = unicodedata.normalize("NFKC", str(value)).strip().lower()
    text = text.replace("μ", "µ")
    text = re.sub(r"[\u2010-\u2015\u2212]", "-", text)
    text = re.sub(r"[^a-z0-9µ%+]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _normalise_unit(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = unicodedata.normalize("NFKC", str(value)).strip().lower()
    text = text.replace("μ", "µ")
    text = re.sub(r"\s+", "", text)
    return text or None


def _unit_from_header(column: Any) -> str | None:
    name = _normalise_name(column)
    if "mg l" in name or "mg/l" in name:
        return "mg/L"
    if "us cm" in name or "µs cm" in name or "µs/cm" in name:
        return "µS/cm"
    if "per mil" in name or "permil" in name or "‰" in str(column):
        return "‰"
    if "temperature c" in name or name.endswith(" temp c"):
        return "°C"
    if name.endswith(" km") or " km " in f" {name} ":
        return "km"
    if name.endswith(" mm") or "amount mm" in name:
        return "mm"
    # A suffix of _m is a useful declared unit, but avoid treating arbitrary
    # names such as "community" as metres.
    if name.endswith(" m") or name.endswith("_m"):
        return "m"
    return None


def _unit_from_value(value: Any) -> str | None:
    normalised = _normalise_unit(value)
    if normalised not in UNIT_TOKENS:
        return None
    return {
        "m": "m",
        "mg/l": "mg/L",
        "mg/litre": "mg/L",
        "us/cm": "µS/cm",
        "µs/cm": "µS/cm",
        "°c": "°C",
        "c": "°C",
        "‰": "‰",
        "tu": "TU",
        "mm": "mm",
    }[normalised]


def _detect_unit_row(frame: pd.DataFrame) -> int | None:
    """Return a likely unit-row index, never guessing from numeric values."""

    if frame.empty:
        return None
    for index in frame.index[:3]:
        values = [value for value in frame.loc[index].tolist() if not pd.isna(value)]
        if not values:
            continue
        known = sum(_unit_from_value(value) is not None for value in values)
        if known >= 3 and known / len(values) >= 0.35:
            return int(index)
    return None


def _read_csv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="cp1252")


def _read_tables(path: Path) -> list[tuple[str, pd.DataFrame]]:
    if path.suffix.lower() == ".csv":
        return [(path.stem, _read_csv(path))]
    workbook = pd.ExcelFile(path)
    return [(str(sheet), pd.read_excel(path, sheet_name=sheet)) for sheet in workbook.sheet_names]


def _source_columns(frame: pd.DataFrame, canonical: str) -> list[str]:
    aliases = {_normalise_name(alias) for alias in CANONICAL_ALIASES[canonical]}
    return [
        str(column)
        for column in frame.columns
        if _normalise_name(column) in aliases
    ]


def _column_unit(
    column: str,
    frame: pd.DataFrame,
    unit_row: int | None,
) -> str | None:
    header_unit = _unit_from_header(column)
    if header_unit is not None:
        return header_unit
    if unit_row is not None:
        return _unit_from_value(frame.loc[unit_row, column])
    return None


def _coverage_record(
    frame: pd.DataFrame,
    canonical: str,
    unit_row: int | None,
) -> dict[str, Any]:
    columns = _source_columns(frame, canonical)
    effective = frame.drop(index=unit_row) if unit_row is not None else frame
    source_units = {
        column: _column_unit(column, frame, unit_row) for column in columns
    }
    nonmissing = {
        column: int(effective[column].notna().sum()) for column in columns
    }
    declared = {unit for unit in source_units.values() if unit is not None}
    if not columns:
        unit_status = "not_applicable"
    elif not declared:
        unit_status = "not_declared"
    elif len(declared) == 1:
        unit_status = "declared"
    else:
        unit_status = "mixed"
    return {
        "available": bool(columns),
        "source_columns": columns,
        "source_units": source_units,
        "nonmissing_rows": nonmissing,
        "unit_status": unit_status,
    }


def _duplicate_summary(
    frame: pd.DataFrame,
    *,
    unit_row: int | None,
    sample_columns: Iterable[str],
) -> dict[str, Any]:
    effective = frame.drop(index=unit_row) if unit_row is not None else frame
    duplicated = effective.duplicated(keep="first")
    key_summary: dict[str, Any] = {}
    for column in sample_columns:
        values = effective[column].dropna().astype(str).str.strip()
        counts = values.value_counts()
        repeated = counts[counts > 1]
        key_summary[column] = {
            "nonmissing": int(values.size),
            "unique": int(values.nunique()),
            "duplicate_groups": int(repeated.size),
            "duplicate_extra_rows": int((repeated - 1).sum()),
            "repeated_values": {
                str(key): int(value) for key, value in repeated.items()
            },
        }
    return {
        "exact_duplicate_groups": int(
            effective.loc[duplicated].drop_duplicates().shape[0]
        ),
        "exact_duplicate_extra_rows": int(duplicated.sum()),
        "key_duplicates": key_summary,
    }


def _flow_proxies(
    sheet: str,
    frame: pd.DataFrame,
    coverage: Mapping[str, Mapping[str, Any]],
) -> list[dict[str, Any]]:
    source = {
        field: list(details.get("source_columns", []))
        for field, details in coverage.items()
    }
    proxies: list[dict[str, Any]] = []
    if source["hydraulic_head"]:
        proxies.append(
            {
                "name": "explicit_hydraulic_head",
                "source_columns": source["hydraulic_head"],
                "status": "available",
                "caveat": "head field is present; datum and measurement uncertainty require metadata",
            }
        )
    if source["elevation"] and source["water_level"]:
        proxies.append(
            {
                "name": "derived_water_table_head",
                "source_columns": source["elevation"] + source["water_level"],
                "status": "derivable",
                "caveat": "requires Static_Water_Level to be depth below the elevation datum; do not treat as measured hydraulic head without confirmation",
            }
        )
    if source["elevation"] and source["latitude"] and source["longitude"]:
        proxies.append(
            {
                "name": "elevation_gradient_context",
                "source_columns": source["elevation"]
                + source["latitude"]
                + source["longitude"],
                "status": "weak_context_only",
                "caveat": "topographic/elevation gradient is not a hydraulic-head measurement and cannot establish well-to-well flow direction by itself",
            }
        )
    if source["distance_river"]:
        proxies.append(
            {
                "name": "distance_to_river_context",
                "source_columns": source["distance_river"],
                "status": "context_only",
                "caveat": "distance to river is not a directed flow observation",
            }
        )
    if source["rainfall"]:
        proxies.append(
            {
                "name": "rainfall_recharge_context",
                "source_columns": source["rainfall"],
                "status": "context_only",
                "caveat": "rainfall can inform recharge timing but does not identify groundwater edge direction",
            }
        )
    # Monitoring-well sheets often contain one date column plus a set of well
    # columns whose semantics/units are not in the workbook.  Detect them
    # without assigning a physical interpretation.
    date_columns = set(source["date"])
    numeric_context = [
        str(column)
        for column in frame.columns
        if str(column) not in date_columns
        and pd.api.types.is_numeric_dtype(frame[column])
        and frame[column].notna().any()
    ]
    if (
        source["date"]
        and len(numeric_context) >= 3
        and "monitor" in sheet.lower()
    ):
        proxies.append(
            {
                "name": "monitoring_time_series",
                "source_columns": source["date"] + numeric_context,
                "status": "potential_but_unverified",
                "caveat": "time-series values may support hydraulic dynamics only after well identity, variable semantics, and units are documented",
            }
        )
    return proxies


def _sheet_inventory(sheet: str, frame: pd.DataFrame) -> dict[str, Any]:
    unit_row = _detect_unit_row(frame)
    effective = frame.drop(index=unit_row) if unit_row is not None else frame
    coverage = {
        field: _coverage_record(frame, field, unit_row)
        for field in (
            "sample_id",
            "site_group",
            "latitude",
            "longitude",
            "elevation",
            "depth",
            "date",
            "hydraulic_head",
            "water_level",
            "SiO2",
            "Sr",
            "Fe",
            "distance_river",
            "rainfall",
        )
    }
    sample_columns = coverage["sample_id"]["source_columns"]
    duplicates = _duplicate_summary(
        frame,
        unit_row=unit_row,
        sample_columns=sample_columns,
    )
    return {
        "sheet": sheet,
        "raw_rows": int(len(frame)),
        "effective_rows": int(len(effective)),
        "unit_row_index": unit_row,
        "unit_row_detected": unit_row is not None,
        "columns": [str(column) for column in frame.columns],
        "canonical_fields": coverage,
        "flow_proxies": _flow_proxies(sheet, frame, coverage),
        "duplicates": duplicates,
    }


def _cross_sheet_overlaps(
    sheets: list[Mapping[str, Any]],
    frames: Mapping[str, pd.DataFrame],
) -> list[dict[str, Any]]:
    overlaps: list[dict[str, Any]] = []
    for left_index, left in enumerate(sheets):
        left_name = str(left["sheet"])
        left_columns = left["canonical_fields"]["sample_id"]["source_columns"]
        if not left_columns:
            continue
        left_values = (
            frames[left_name][left_columns[0]]
            .dropna()
            .astype(str)
            .str.strip()
        )
        left_set = set(left_values)
        for right in sheets[left_index + 1 :]:
            right_name = str(right["sheet"])
            right_columns = right["canonical_fields"]["sample_id"]["source_columns"]
            if not right_columns:
                continue
            right_values = (
                frames[right_name][right_columns[0]]
                .dropna()
                .astype(str)
                .str.strip()
            )
            overlap = sorted(left_set & set(right_values))
            if overlap:
                overlaps.append(
                    {
                        "left_sheet": left_name,
                        "right_sheet": right_name,
                        "sample_id_overlap_count": len(overlap),
                        "sample_id_overlap_examples": overlap[:10],
                        "interpretation": "repeated identifiers across sheets; determine whether these are seasonal/repeated measurements before stacking",
                    }
                )
    return overlaps


def _dataset_summary(
    name: str,
    path: Path,
    tables: list[tuple[str, pd.DataFrame]],
) -> dict[str, Any]:
    frames = {sheet: frame for sheet, frame in tables}
    sheets = [_sheet_inventory(sheet, frame) for sheet, frame in tables]
    union: dict[str, dict[str, Any]] = {}
    for field in sheets[0]["canonical_fields"] if sheets else {}:
        available_sheets = [
            sheet["sheet"]
            for sheet in sheets
            if sheet["canonical_fields"][field]["available"]
        ]
        units = sorted(
            {
                unit
                for sheet in sheets
                for unit in sheet["canonical_fields"][field]["source_units"].values()
                if unit is not None
            }
        )
        coverage_by_sheet: dict[str, dict[str, Any]] = {}
        for sheet in sheets:
            item = sheet["canonical_fields"][field]
            source_columns = item["source_columns"]
            nonmissing = item["nonmissing_rows"]
            n_with_any = max(
                (int(nonmissing[column]) for column in source_columns),
                default=0,
            )
            denominator = int(sheet["effective_rows"])
            coverage_by_sheet[str(sheet["sheet"])] = {
                "rows_with_value": n_with_any,
                "effective_rows": denominator,
                "fraction": (
                    float(n_with_any / denominator) if denominator else None
                ),
            }
        union[field] = {
            "available_in_any_sheet": bool(available_sheets),
            "available_sheets": available_sheets,
            "source_columns": sorted(
                {
                    column
                    for sheet in sheets
                    for column in sheet["canonical_fields"][field]["source_columns"]
                }
            ),
            "declared_units": units,
            "coverage_by_sheet": coverage_by_sheet,
            "unit_status": (
                "not_applicable"
                if not available_sheets
                else "declared"
                if units
                else "not_declared"
            ),
        }

    direct_requested = {
        field: {
            "available_in_any_sheet": union.get(field, {}).get(
                "available_in_any_sheet", False
            ),
            "available_sheets": union.get(field, {}).get("available_sheets", []),
            "source_columns": union.get(field, {}).get("source_columns", []),
            "declared_units": union.get(field, {}).get("declared_units", []),
            "coverage_by_sheet": union.get(field, {}).get(
                "coverage_by_sheet", {}
            ),
            "unit_status": union.get(field, {}).get("unit_status", "not_applicable"),
        }
        for field in REQUESTED_FIELDS
    }
    flow_proxy_summaries = [
        {
            "sheet": sheet["sheet"],
            **proxy,
        }
        for sheet in sheets
        for proxy in sheet["flow_proxies"]
    ]
    return {
        "dataset": name,
        "path": str(path),
        "exists": bool(path.exists()),
        "file_type": path.suffix.lower().lstrip("."),
        "raw_rows": int(sum(sheet["raw_rows"] for sheet in sheets)),
        "effective_rows": int(sum(sheet["effective_rows"] for sheet in sheets)),
        "sheets": sheets,
        "canonical_fields": union,
        "requested_field_availability": direct_requested,
        "flow_proxies": flow_proxy_summaries,
        "cross_sheet_sample_id_overlaps": _cross_sheet_overlaps(sheets, frames),
    }


def inventory_field_data(
    repo_root: Path | str = REPO_ROOT,
    datasets: Mapping[str, Path | str] | None = None,
) -> dict[str, Any]:
    """Return a read-only inventory for the configured field-data products."""

    root = Path(repo_root).resolve()
    configured = datasets or DEFAULT_DATASETS
    report: dict[str, Any] = {
        "schema": "field-data-inventory-v1",
        "repo_root": str(root),
        "datasets": [],
    }
    for name, relative in configured.items():
        path = Path(relative)
        if not path.is_absolute():
            path = root / path
        if not path.exists():
            report["datasets"].append(
                {
                    "dataset": name,
                    "path": str(path),
                    "exists": False,
                    "error": "source file is missing",
                }
            )
            continue
        tables = _read_tables(path)
        report["datasets"].append(_dataset_summary(name, path, tables))
    return report


def _fmt_units(units: Iterable[str]) -> str:
    values = list(units)
    return ", ".join(values) if values else "not declared"


def render_markdown(report: Mapping[str, Any]) -> str:
    """Render a compact human-readable inventory without changing the report."""

    lines = [
        "# Field-data inventory and harmonisation diagnostic",
        "",
        "Read-only report; no source file was modified and undocumented units are not inferred.",
        "",
        "| Dataset | Sheet(s) | Raw rows | Effective rows | Exact duplicate extras |",
        "|---|---:|---:|---:|---:|",
    ]
    for dataset in report.get("datasets", []):
        if not dataset.get("exists"):
            lines.append(
                f"| {dataset['dataset']} | missing | — | — | — |"
            )
            continue
        sheet_names = ", ".join(str(sheet["sheet"]) for sheet in dataset["sheets"])
        duplicate_rows = sum(
            int(sheet["duplicates"]["exact_duplicate_extra_rows"])
            for sheet in dataset["sheets"]
        )
        lines.append(
            f"| {dataset['dataset']} | {sheet_names} | {dataset['raw_rows']} | {dataset['effective_rows']} | {duplicate_rows} |"
        )

    lines.extend(
        [
            "",
            "## Requested field coverage",
            "",
            "| Dataset | SiO2 | Sr | Fe | Depth | Date |",
            "|---|---|---|---|---|---|",
        ]
    )
    for dataset in report.get("datasets", []):
        availability = dataset.get("requested_field_availability", {})
        if not availability:
            lines.append(f"| {dataset['dataset']} | missing source | — | — | — | — |")
            continue
        values = []
        for field in REQUESTED_FIELDS:
            item = availability[field]
            if not item["available_in_any_sheet"]:
                values.append("no")
            else:
                sheet_text = "; ".join(
                    f"{sheet}: {details['rows_with_value']}/{details['effective_rows']}"
                    for sheet, details in item.get("coverage_by_sheet", {}).items()
                    if sheet in item["available_sheets"]
                )
                values.append(
                    f"yes ({sheet_text}; {_fmt_units(item['declared_units'])})"
                )
        lines.append(f"| {dataset['dataset']} | " + " | ".join(values) + " |")

    lines.extend(["", "## Flow-proxy availability", ""])
    for dataset in report.get("datasets", []):
        proxies = dataset.get("flow_proxies", [])
        lines.append(f"**{dataset['dataset']}**")
        if not proxies:
            lines.append("- No flow proxy field detected.")
            continue
        for proxy in proxies:
            columns = ", ".join(proxy.get("source_columns", []))
            lines.append(
                f"- `{proxy['name']}` ({proxy['status']}): {columns}. {proxy['caveat']}"
            )

    lines.extend(["", "## Sheet-level canonical fields and units", ""])
    for dataset in report.get("datasets", []):
        if not dataset.get("exists"):
            continue
        lines.append(f"**{dataset['dataset']}**")
        for sheet in dataset["sheets"]:
            present = []
            for field, item in sheet["canonical_fields"].items():
                if item["available"]:
                    units = _fmt_units(
                        unit for unit in item["source_units"].values() if unit
                    )
                    present.append(f"{field}={','.join(item['source_columns'])} [{units}]")
            unit_row = (
                f"unit row {sheet['unit_row_index']} removed"
                if sheet["unit_row_detected"]
                else "no unit row detected"
            )
            lines.append(
                f"- `{sheet['sheet']}` ({sheet['effective_rows']} effective rows; {unit_row}): "
                + ("; ".join(present) if present else "no canonical fields")
            )
        for overlap in dataset.get("cross_sheet_sample_id_overlaps", []):
            lines.append(
                f"- Cross-sheet sample-ID overlap {overlap['left_sheet']} vs {overlap['right_sheet']}: "
                f"{overlap['sample_id_overlap_count']} ({overlap['interpretation']})."
            )
    return "\n".join(lines) + "\n"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=REPO_ROOT,
        help="repository root containing data/FieldData (default: inferred)",
    )
    parser.add_argument(
        "--format",
        choices=("markdown", "json"),
        default="markdown",
        help="stdout format (default: markdown)",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    report = inventory_field_data(args.repo_root)
    if args.format == "json":
        print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    else:
        print(render_markdown(report), end="")
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised by CLI smoke tests
    raise SystemExit(main())


__all__ = [
    "CANONICAL_ALIASES",
    "DEFAULT_DATASETS",
    "REQUESTED_FIELDS",
    "inventory_field_data",
    "main",
    "render_markdown",
]
