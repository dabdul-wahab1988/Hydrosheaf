"""Auditable harmonisation of the available groundwater field datasets.

The repository contains four field-data packages with deliberately different
schemas and analytical coverage.  This module is the single, read-only
boundary between those native files and HydroSheaf's internal sample schema.
It converts concentrations to ``mmol/L``, preserves native values/metadata,
records source hashes, and leaves unavailable observations as ``None`` rather
than manufacturing zeros.  The loader does not infer flow truth or lithology
from chemistry; mapped geology is accepted only through an explicitly supplied
join table.

The public helpers are intentionally dependency-light (``pandas`` is used only
for CSV/XLSX decoding) so that field-data audits can run without the inference
stack.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import pandas as pd

from .units import mgL_to_mmolL


FIELD_DATA_RELATIVE_PATHS: dict[str, Path] = {
    "lower_anayari": Path("LowerAnayari") / "manu.csv",
    # The spelling is part of the on-disk provenance and is retained for
    # compatibility with the canonical repository package.
    "northen_ghana": Path("NorthenGhana") / "NorthernGhana.xlsx",
    "northern_ghana_new": Path("NorthernGhanaNew") / "compiled UER data_new.xlsx",
    "talensi_mining_area": Path("Talensi_MiningArea") / "talensi.csv",
}

CANONICAL_IONS: tuple[str, ...] = (
    "Ca",
    "Mg",
    "Na",
    "K",
    "HCO3",
    "Cl",
    "SO4",
    "NO3",
    "F",
    "Fe",
    "PO4",
    "SiO2",
    "Sr",
)

# Species not in data.units because they are optional diagnostic tracers in the
# core model.  Values are mg per mmol (numerically the molar mass in g/mol).
_OPTIONAL_MOLAR_MASS = {"SiO2": 60.0843, "Sr": 87.62}


def _finite(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _safe_mmol(value: Any, ion: str, *, unit: str = "mg/L") -> Optional[float]:
    number = _finite(value)
    if number is None or number < 0.0:
        return None
    if unit.lower().replace(" ", "") in {"mmol/l", "mmol_l", "mmoll"}:
        return number
    if ion in _OPTIONAL_MOLAR_MASS:
        return number / _OPTIONAL_MOLAR_MASS[ion]
    try:
        return mgL_to_mmolL(number, ion)
    except KeyError:
        return None


def _native_columns(frame: pd.DataFrame) -> tuple[str, ...]:
    return tuple(str(column) for column in frame.columns)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sample_id(prefix: str, value: Any) -> str:
    text = str(value).strip()
    return f"{prefix}_{text}" if text else f"{prefix}_unknown"


def _parse_coordinate(value: Any) -> Optional[float]:
    """Parse decimal or simple DMS coordinates without guessing a hemisphere."""

    direct = _finite(value)
    if direct is not None:
        return direct
    text = str(value or "").strip().replace("º", "°")
    if not text:
        return None
    import re

    match = re.match(
        r"^\s*(?P<deg>-?\d+(?:\.\d+)?)\s*(?:°|deg)?\s*"
        r"(?:(?P<min>\d+(?:\.\d+)?)\s*['′])?\s*"
        r"(?:(?P<sec>\d+(?:\.\d+)?)\s*[\"″])?\s*"
        r"(?P<hem>[NSEW])?\s*$",
        text,
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    degrees = float(match.group("deg"))
    minutes = float(match.group("min") or 0.0)
    seconds = float(match.group("sec") or 0.0)
    if minutes >= 60.0 or seconds >= 60.0:
        return None
    magnitude = abs(degrees) + minutes / 60.0 + seconds / 3600.0
    sign = -1.0 if degrees < 0.0 else 1.0
    hemisphere = (match.group("hem") or "").upper()
    if hemisphere in {"S", "W"}:
        sign = -1.0
    elif hemisphere in {"N", "E"}:
        sign = 1.0
    return sign * magnitude


def _base_record(
    *,
    dataset: str,
    sample_id: str,
    site_id: str,
    season: str = "unknown",
    native_row: Optional[int] = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "dataset": dataset,
        "sample_id": sample_id,
        # ``node_id`` is the graph-safe observation key.  It is distinct from
        # ``site_id`` because one physical site may have repeated seasonal or
        # campaign observations.
        "node_id": sample_id,
        "site_id": site_id,
        "season": season,
        "native_row": native_row,
        "sample_date": None,
        "sample_year": None,
        "lat": None,
        "lon": None,
        "elevation": None,
        "well_depth": None,
        "static_water_level": None,
        "screen_top": None,
        "screen_bottom": None,
        "pH": None,
        "EC": None,
        "TDS": None,
        "temp_c": None,
        "18O": None,
        "2H": None,
        "coordinate_quality": "unavailable",
        "geology_join_status": "unavailable",
        "geology_source_hash": None,
    }
    for ion in CANONICAL_IONS:
        record[ion] = None
    return record


def _fill_mgL(record: dict[str, Any], source: Mapping[str, Any], mapping: Mapping[str, str]) -> None:
    for ion, column in mapping.items():
        record[ion] = _safe_mmol(source.get(column), ion)


def _fill_common(
    record: dict[str, Any],
    source: Mapping[str, Any],
    *,
    pH: str = "pH",
    ec: str = "EC",
    tds: str = "TDS",
    temp: str = "Temp",
    d18: str = "d18O",
    d2h: str = "d2H",
) -> None:
    record["pH"] = _finite(source.get(pH))
    record["EC"] = _finite(source.get(ec))
    record["TDS"] = _finite(source.get(tds))
    record["temp_c"] = _finite(source.get(temp))
    record["18O"] = _finite(source.get(d18))
    record["2H"] = _finite(source.get(d2h))


def _load_lower_anayari(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    frame = pd.read_csv(path)
    records: list[dict[str, Any]] = []
    mapping = {
        "Na": "Na",
        "K": "K",
        "Ca": "Ca",
        "Mg": "Mg",
        "F": "F",
        "Cl": "Cl",
        "HCO3": "HCO3",
        "NO3": "NO3",
        "SO4": "SO4",
        "Fe": "Fe",
    }
    for index, row in frame.iterrows():
        native = row.to_dict()
        sample_id = _sample_id("LA", native.get("Sample ID"))
        record = _base_record(
            dataset="lower_anayari",
            sample_id=sample_id,
            site_id=sample_id,
            native_row=int(index) + 2,
        )
        record["station"] = native.get("Station")
        record["native_sample_id"] = str(native.get("Sample ID") or "").strip()
        record["lat"] = _parse_coordinate(native.get("Y coordinate"))
        record["lon"] = _parse_coordinate(native.get("X coordinate"))
        record["elevation"] = _finite(native.get("Elevation"))
        record["coordinate_quality"] = "decimal_source_labels"
        _fill_common(record, native, ec="EC", tds="TDS", temp="Temp")
        _fill_mgL(record, native, mapping)
        records.append(record)
    return records, {"native_columns": _native_columns(frame), "sheet_names": []}


def _load_talensi(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    frame = pd.read_csv(path)
    records: list[dict[str, Any]] = []
    mapping = {
        "HCO3": "HCO3",
        "Na": "Na",
        "K": "K",
        "Ca": "Ca",
        "Mg": "Mg",
        "Cl": "Cl",
        "SO4": "SO4",
        "NO3": "NO3",
        "Fe": "Fe",
    }
    for index, row in frame.iterrows():
        native = row.to_dict()
        code = str(native.get("Code") or f"row-{index + 1}").strip()
        record = _base_record(
            dataset="talensi_mining_area",
            sample_id=_sample_id("TA", code),
            site_id=_sample_id("TA", code),
            native_row=int(index) + 2,
        )
        record["town"] = native.get("Town")
        record["native_sample_id"] = code
        record["lat"] = _parse_coordinate(native.get("Latitude"))
        # The native table stores a positive west-longitude magnitude.  Keep
        # the signed value used for spatial work and document the correction.
        longitude = _parse_coordinate(native.get("Longitude"))
        record["lon"] = -abs(longitude) if longitude is not None else None
        record["coordinate_quality"] = "decimal_west_magnitude_corrected"
        record["elevation"] = _finite(native.get("Elevation"))
        record["Eh"] = _finite(native.get("Eh"))
        record["Sal"] = _finite(native.get("Sal"))
        _fill_common(record, native, ec="EC", tds="TDS", temp="Temp")
        _fill_mgL(record, native, mapping)
        records.append(record)
    return records, {"native_columns": _native_columns(frame), "sheet_names": []}


def _load_northern_ghana(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    sheet_columns: dict[str, tuple[str, ...]] = {}
    mapping = {
        "Ca": "Ca_mg_L",
        "Mg": "Mg_mg_L",
        "Na": "Na_mg_L",
        "K": "K_mg_L",
        "HCO3": "HCO3_mg_L",
        "Cl": "Cl_mg_L",
        "SO4": "SO4_mg_L",
        "NO3": "NO3_mg_L",
        "F": "F_mg_L",
        "Sr": "Sr_mg_L",
        "SiO2": "SiO2_mg_L",
    }
    for season in ("Dry", "Wet"):
        frame = pd.read_excel(path, sheet_name=season)
        sheet_columns[season] = _native_columns(frame)
        for index, row in frame.iterrows():
            native = row.to_dict()
            well = str(native.get("Well_ID") or f"row-{index + 1}").strip()
            sample_id = _sample_id("NG", f"{well}_{season.lower()}")
            record = _base_record(
                dataset="northen_ghana",
                sample_id=sample_id,
                site_id=_sample_id("NG", well),
                season=season.lower(),
                native_row=int(index) + 2,
            )
            record["well_id"] = well
            record["native_sample_id"] = well
            record["region"] = native.get("Region")
            record["district"] = native.get("District")
            record["community_code"] = native.get("Community_Code")
            record["lat"] = _parse_coordinate(native.get("Latitude"))
            record["lon"] = _parse_coordinate(native.get("Longitude"))
            record["elevation"] = _finite(native.get("Elevation_m"))
            record["well_depth"] = _finite(native.get("Borehole_Depth_m"))
            record["static_water_level"] = _finite(native.get("Static_Water_Level_m"))
            record["distance_river_km"] = _finite(native.get("Distance_River_km"))
            record["distance_farm_km"] = _finite(native.get("Distance_Farm_km"))
            record["distance_settlement_km"] = _finite(native.get("Distance_Settlement_km"))
            record["coordinate_quality"] = "decimal_native"
            _fill_common(
                record,
                native,
                ec="EC_uS_cm",
                tds="TDS_mg_L",
                temp="Temperature_C",
                d18="d18O_permil",
                d2h="d2H_permil",
            )
            _fill_mgL(record, native, mapping)
            records.append(record)
    return records, {"native_columns": sheet_columns, "sheet_names": ["Dry", "Wet"]}


def _load_northern_ghana_new(
    path: Path,
    *,
    geology_join_path: Optional[Path] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    frame = pd.read_excel(path, sheet_name="GW", header=0)
    # Row 2 in the native workbook is a units legend rather than a sample.
    frame = frame[pd.to_numeric(frame.get("No."), errors="coerce").notna()].copy()
    mapping = {
        "Ca": "Ca2+",
        "Mg": "Mg 2+",
        "Na": "Na+",
        "K": "K+ ",
        "HCO3": "HCO3-",
        "Cl": "Cl-",
        "SO4": "SO42-",
        "NO3": "NO3-",
        "F": "F-",
    }
    records: list[dict[str, Any]] = []
    for index, row in frame.iterrows():
        native = row.to_dict()
        number = int(float(native["No."]))
        sample_id = _sample_id("NGN", f"{number:03d}")
        record = _base_record(
            dataset="northern_ghana_new",
            sample_id=sample_id,
            site_id=sample_id,
            native_row=int(index) + 2,
        )
        record["sample_no"] = number
        record["community"] = native.get("Community")
        record["sample_type"] = native.get("Type")
        lon_raw = _parse_coordinate(native.get("Longitude"))
        lat_raw = _parse_coordinate(native.get("Latituted"))
        # The first block of the workbook has latitude/longitude labels
        # reversed (positive ~10 first, negative ~1 second).  DMS rows near
        # the end are parsed by hemisphere and are already ordered.  The
        # range rule is explicit and auditable, not a free-form spatial guess.
        if lon_raw is not None and lat_raw is not None and lon_raw > 0.0 and lat_raw < 0.0:
            record["lat"], record["lon"] = lon_raw, lat_raw
            record["coordinate_quality"] = "range_reordered_source_labels"
        elif lon_raw is not None and lat_raw is not None:
            record["lon"], record["lat"] = lon_raw, lat_raw
            record["coordinate_quality"] = "decimal_or_dms_native"
        record["elevation"] = _finite(native.get("Elev"))
        _fill_common(
            record,
            native,
            ec="EC",
            tds="TDS",
            temp="Temp",
            d18="δ18O",
            d2h="δ2H",
        )
        _fill_mgL(record, native, mapping)
        record["3H"] = _finite(native.get("    3H"))
        record["d15N"] = _finite(native.get("d15N-NO3 (Air)"))
        record["d18O_NO3"] = _finite(native.get("d18O-NO3 (VSMOW)"))
        records.append(record)

    metadata: dict[str, Any] = {
        "native_columns": _native_columns(frame),
        "sheet_names": ["GW", "rain", "monitoring wells"],
    }
    if geology_join_path is not None and geology_join_path.exists():
        _attach_geology_join(records, geology_join_path, key_mode="sample_no")
        metadata["geology_join_path"] = str(geology_join_path)
        metadata["geology_join_sha256"] = _sha256(geology_join_path)
    return records, metadata


def _attach_geology_join(
    records: list[dict[str, Any]], join_path: Path, *, key_mode: str = "sample_id"
) -> None:
    """Attach a precomputed spatial join by ``sample_no`` only.

    Geometry is intentionally not recomputed here.  The join table must carry
    its own source CRS/predicate/hash metadata; unmatched records remain
    unknown and therefore receive neutral priors downstream.
    """

    join = pd.read_csv(join_path)
    if key_mode == "sample_no":
        lookup = {
            ("", int(float(row["sample_no"]))): row.to_dict()
            for _, row in join.iterrows()
            if _finite(row.get("sample_no")) is not None
        }
    elif key_mode == "well_season":
        lookup = {
            (
                str(row.get("Well_ID") or "").strip(),
                str(row.get("source_sheet") or "").strip().lower(),
            ): row.to_dict()
            for _, row in join.iterrows()
        }
    else:
        source_column = "Sample ID" if "Sample ID" in join.columns else "Code"
        lookup = {
            (str(row.get(source_column) or "").strip(), ""): row.to_dict()
            for _, row in join.iterrows()
        }
    copied = (
        "geology_code_1000",
        "geology_symbol",
        "geology_stratigraphic_unit",
        "geology_age_upper_ma",
        "geology_age_lower_ma",
        "geology_stratigraphic_formation",
        "geology_tectonic_domain",
        "geology_sub_domain",
        "geology_metamorphic_grade",
        "geology_legend_text",
        "geology_source_feature_index",
        "geology_boundary_distance_m",
        "geology_nearest_distance_m",
        "geology_nearest_symbol",
        "geology_join_status",
        "geology_match_count",
        "geology_source_layer",
        "geology_source_crs",
        "geology_point_crs",
        "geology_join_predicate",
        # The UER sidecar records the audited coordinate repair separately
        # from the chemistry workbook.  Preserve it and use the candidate
        # coordinate only when both values are finite; a failed DMS parse must
        # remain missing rather than being guessed from a map polygon.
        "latitude_candidate_dd",
        "longitude_candidate_dd",
        "coordinate_assignment_method",
        "coordinate_label_flag",
        "longitude_parse_kind",
        "latitude_parse_kind",
        "longitude_parse_flag",
        "latitude_parse_flag",
        "geology_coordinate_assumption",
    )
    join_hash = _sha256(join_path)
    for record in records:
        if key_mode == "sample_no":
            row = lookup.get(("", int(record.get("sample_no", -1))))
        elif key_mode == "well_season":
            row = lookup.get(
                (
                    str(record.get("well_id") or "").strip(),
                    str(record.get("season") or "").strip().lower(),
                )
            )
        else:
            row = lookup.get((str(record.get("native_sample_id") or "").strip(), ""))
        if row is None:
            continue
        for key in copied:
            if key in row:
                value = row[key]
                if pd.isna(value):
                    value = None
                record[key] = value
        record["geology_source_hash"] = join_hash
        record["geology_join_status"] = str(row.get("geology_join_status") or "UNKNOWN")
        candidate_lat = _finite(row.get("latitude_candidate_dd"))
        candidate_lon = _finite(row.get("longitude_candidate_dd"))
        coordinate_flag = " ".join(
            str(row.get(key) or "")
            for key in (
                "coordinate_label_flag",
                "longitude_parse_flag",
                "latitude_parse_flag",
            )
        )
        if candidate_lat is not None and candidate_lon is not None:
            record["lat"] = candidate_lat
            record["lon"] = candidate_lon
            record["coordinate_quality"] = (
                "audited_sidecar_candidate_review"
                if "OUT_OF_RANGE" in coordinate_flag.upper()
                else "audited_sidecar_candidate"
            )
        if row.get("coordinate_label_flag") is not None:
            record["coordinate_label_flag"] = row.get("coordinate_label_flag")


@dataclass(frozen=True)
class FieldDataset:
    """Harmonised records plus immutable source provenance."""

    name: str
    source_path: str
    source_sha256: str
    records: tuple[Mapping[str, Any], ...]
    metadata: Mapping[str, Any] = field(default_factory=dict)
    auxiliary_tables: Mapping[str, tuple[Mapping[str, Any], ...]] = field(default_factory=dict)

    @property
    def n_records(self) -> int:
        return len(self.records)

    def coverage(self) -> dict[str, int]:
        keys = tuple(
            dict.fromkeys(
                CANONICAL_IONS
                + ("pH", "18O", "2H", "elevation", "well_depth", "static_water_level")
            )
        )
        return {
            key: sum(_finite(row.get(key)) is not None for row in self.records)
            for key in keys
        }

    def as_records(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.records]

    def manifest_record(self) -> dict[str, Any]:
        return {
            "dataset": self.name,
            "source_path": self.source_path,
            "source_sha256": self.source_sha256,
            "n_records": self.n_records,
            "coverage": self.coverage(),
            "metadata": dict(self.metadata),
            "auxiliary_tables": {
                key: len(value) for key, value in self.auxiliary_tables.items()
            },
        }


def default_field_data_root(repo_root: Optional[Path] = None) -> Path:
    if repo_root is not None:
        return Path(repo_root)
    return Path(__file__).resolve().parents[2] / "data" / "FieldData"


def _default_geology_join_path(field_root: Path) -> Optional[Path]:
    # This helper is retained for the alternate workbook; canonical packages
    # use the package-specific sidecars selected in ``load_field_dataset``.
    candidate = field_root.parent.parent / "outputs" / "2026-09-05_northern_ghana_geology_join" / "NorthernGhanaNew_geology_join.csv"
    return candidate if candidate.exists() else None


def _load_auxiliary_tables(path: Path, sheet_names: Sequence[str]) -> dict[str, tuple[Mapping[str, Any], ...]]:
    tables: dict[str, tuple[Mapping[str, Any], ...]] = {}
    for sheet in sheet_names:
        if sheet == "GW":
            continue
        try:
            frame = pd.read_excel(path, sheet_name=sheet)
        except Exception:
            continue
        # Keep auxiliary series (rainfall/monitoring) available to temporal
        # workflows, but never promote them to sample chemistry implicitly.
        tables[sheet] = tuple(row.to_dict() for _, row in frame.iterrows())
    return tables


def load_field_dataset(
    name: str,
    *,
    field_root: Optional[Path] = None,
    geology_join_path: Optional[Path] = None,
) -> FieldDataset:
    """Load one of the four canonical field-data packages."""

    key = str(name).strip().lower()
    aliases = {
        "loweranayari": "lower_anayari",
        "manu": "lower_anayari",
        "northenghana": "northen_ghana",
        "northern_ghana": "northen_ghana",
        "northernghananew": "northern_ghana_new",
        "northern_ghana_new": "northern_ghana_new",
        "talensi": "talensi_mining_area",
        "talensi_miningarea": "talensi_mining_area",
    }
    key = aliases.get(key, key)
    if key not in FIELD_DATA_RELATIVE_PATHS:
        raise ValueError(f"Unknown field dataset {name!r}; choose {sorted(FIELD_DATA_RELATIVE_PATHS)}")
    root = default_field_data_root(field_root)
    source = root / FIELD_DATA_RELATIVE_PATHS[key]
    if not source.exists():
        raise FileNotFoundError(source)

    join_root = root.parent.parent / "outputs" / "2026-09-05_multi_geology_join"
    if key == "lower_anayari":
        records, metadata = _load_lower_anayari(source)
        join = geology_join_path or (join_root / "LowerAnayari_geology_join.csv")
        if join.exists():
            _attach_geology_join(records, join, key_mode="sample_id")
            metadata["geology_join_path"] = str(join)
            metadata["geology_join_sha256"] = _sha256(join)
    elif key == "talensi_mining_area":
        records, metadata = _load_talensi(source)
        join = geology_join_path or (join_root / "Talensi_geology_join.csv")
        if join.exists():
            _attach_geology_join(records, join, key_mode="sample_id")
            metadata["geology_join_path"] = str(join)
            metadata["geology_join_sha256"] = _sha256(join)
    elif key == "northen_ghana":
        records, metadata = _load_northern_ghana(source)
        join = geology_join_path or (join_root / "NorthernGhanaData_geology_join.csv")
        if join.exists():
            _attach_geology_join(records, join, key_mode="well_season")
            metadata["geology_join_path"] = str(join)
            metadata["geology_join_sha256"] = _sha256(join)
    else:
        join_path = geology_join_path
        if join_path is None:
            join_path = _default_geology_join_path(root)
        records, metadata = _load_northern_ghana_new(source, geology_join_path=join_path)

    auxiliary = {}
    if key == "northern_ghana_new":
        auxiliary = _load_auxiliary_tables(source, metadata.get("sheet_names", ()))
    metadata.setdefault("canonical_ion_unit", "mmol/L")
    metadata.setdefault("source_ion_unit", "mg/L")
    metadata.setdefault(
        "source_ion_unit_status",
        "declared_by_header_or_units_row"
        if key in {"northen_ghana", "northern_ghana_new"}
        else "not_declared_in_native_header_loader_contract_only",
    )
    metadata.setdefault("missing_value_policy", "preserve_none_no_zero_imputation")
    return FieldDataset(
        name=key,
        source_path=str(source),
        source_sha256=_sha256(source),
        records=tuple(records),
        metadata=metadata,
        auxiliary_tables=auxiliary,
    )


def load_all_field_datasets(
    *,
    field_root: Optional[Path] = None,
    geology_join_path: Optional[Path] = None,
) -> dict[str, FieldDataset]:
    """Load all four packages without mutating any source file."""

    return {
        key: load_field_dataset(
            key,
            field_root=field_root,
            geology_join_path=geology_join_path,
        )
        for key in FIELD_DATA_RELATIVE_PATHS
    }


def flatten_field_datasets(
    datasets: Mapping[str, FieldDataset],
) -> list[dict[str, Any]]:
    """Return all harmonised records with collision-safe graph identifiers.

    The source packages remain separate in provenance and in the ``dataset``
    field.  Flattening is only a convenience for a downstream graph builder;
    it does not pool units or fill missing ions.  Repeated physical sites (for
    example Dry/Wet observations) retain distinct ``node_id``/``sample_id``
    values.  Duplicate graph identifiers are rejected rather than overwritten.
    """

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for dataset_name, dataset in datasets.items():
        if not isinstance(dataset, FieldDataset):
            raise TypeError(f"{dataset_name!r} is not a FieldDataset")
        for source_row in dataset.records:
            row = dict(source_row)
            node_id = str(row.get("node_id") or row.get("sample_id") or "").strip()
            if not node_id:
                raise ValueError(
                    f"dataset {dataset_name!r} contains a record without node_id/sample_id"
                )
            if node_id in seen:
                raise ValueError(f"duplicate graph node_id across field datasets: {node_id!r}")
            row["node_id"] = node_id
            seen.add(node_id)
            rows.append(row)
    return rows


def field_data_manifest(
    *,
    field_root: Optional[Path] = None,
    geology_join_path: Optional[Path] = None,
) -> dict[str, Any]:
    """Return JSON-serialisable provenance and coverage for all packages."""

    datasets = load_all_field_datasets(
        field_root=field_root,
        geology_join_path=geology_join_path,
    )
    return {
        "schema_version": "hydrosheaf.field-data.v1",
        "datasets": [datasets[key].manifest_record() for key in FIELD_DATA_RELATIVE_PATHS],
    }


__all__ = [
    "CANONICAL_IONS",
    "FIELD_DATA_RELATIVE_PATHS",
    "FieldDataset",
    "default_field_data_root",
    "field_data_manifest",
    "flatten_field_datasets",
    "load_all_field_datasets",
    "load_field_dataset",
]
