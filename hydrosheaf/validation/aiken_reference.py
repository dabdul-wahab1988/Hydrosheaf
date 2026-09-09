"""Read-only adapter for the USGS Aiken County model/data release.

The Aiken release is a useful *calibrated-model reference* for HydroSheaf:
it contains MODFLOW-NWT/MODPATH5 files alongside groundwater chemistry and
CFC apparent-age tables.  This module deliberately keeps those evidence
types separate.  In particular, a MODPATH recharge-to-well path is not
converted into a direct well--well edge and no reaction or adjacency truth is
manufactured from the chemistry tables.

The adapter accepts either an extracted release directory or the directory
containing the release's ``ancillary.zip``/``model.zip`` files.  Workbook
members are read through :mod:`zipfile` into memory (the workbooks are small).
The small fixed-width MODPATH 5 ``.ept``/``.pth`` particle files are validated
and decoded directly from the archive; multi-gigabyte MODFLOW heads, budgets,
and binary outputs are not extracted or decoded.  ``inventory_aiken_source``
can hash every ZIP member by streaming its contents, and therefore provides a
reproducible provenance inventory without creating a second copy of the archive.

The public functions intentionally return pandas data frames because the
rest of HydroSheaf's benchmark code uses that representation.  Every
canonical table contains ``reference_type`` and explicit ``ABSTAIN``/unknown
fields where the source cannot support a truth claim.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date, datetime, time
from io import BytesIO
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from typing import Any, Iterable, Mapping, Sequence
import unicodedata
import zipfile

import openpyxl
import pandas as pd


REFERENCE_TYPE = "calibrated_model_reference"
SOURCE_DOI = "10.5066/P9U0GHLU"
PUBLICATION_DOI = "10.3133/sir20225036"
SOURCE_URL = (
    "https://www.usgs.gov/data/modflow-nwt-and-modpath5-used-evaluate-"
    "groundwater-availability-geochemistry-and-flow-pathways"
)

CORE_WORKBOOK_NAME = "Tables 2 and 8 thru 15 v 8 18 2020.xlsx"
VOC_WORKBOOK_NAME = "Table 10 Split into 7 tables revised 1 12 21.xlsx"

# MODPATH 5 binary particle-coordinate files are Fortran stream records in
# the Aiken release.  They have an 80-byte text header followed by a single
# precision reference time (TREF), then fixed-width little-endian records.
# The release is known to contain the 32-byte pathline and 56-byte endpoint
# variants documented by USGS.  Keep the layouts explicit here rather than
# relying on a text-oriented or version-inference reader: a silent byte-shift
# in a particle file would create plausible-looking but invalid coordinates.
MODPATH5_HEADER_BYTES = 84
MODPATH5_PATHLINE_RECORD_BYTES = 32
MODPATH5_ENDPOINT_RECORD_BYTES = 56
MODPATH5_BINARY_ENDIAN = "little"
MODPATH5_BINARY_PARSER = "hydrosheaf.validation.aiken_reference.modpath5"
MODPATH5_BINARY_PARSER_VERSION = "1.0"
MODPATH5_IPCODE_ENCODING = "nslast_times_10_plus_idcode"
MODPATH5_TRACKING_DIRECTIONS = {
    1: "forward_in_direction_of_flow",
    2: "backward_toward_recharge",
}
_MODPATH5_HEADER_TEXT_BYTES = 80
_MODPATH5_PATHLINE_STRUCT = struct.Struct("<iffff fii".replace(" ", ""))
_MODPATH5_ENDPOINT_STRUCT = struct.Struct("<ii7fiiiif")

_MISSING_MARKERS = frozenset(
    {
        "",
        "-",
        "--",
        "—",
        "–",
        "‒",
        "−",
        "�",
        "na",
        "n/a",
        "not available",
        "no data",
    }
)
_WELL_RE = re.compile(r"\b(?:AK|LEX)-\d+\b", flags=re.IGNORECASE)
_DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", flags=re.IGNORECASE)
_NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")
_AGE_MODIFIER_RE = re.compile(
    r"(?P<modifier>early|mid|late)\s*(?P<decade>\d{4})\s*(?:['’]?[sS])?"
)
_AGE_TOKEN_RE = re.compile(
    r"(?P<modifier>early|mid|late)\s*(?P<decade>\d{4})\s*(?:['’]?[sS])?",
    flags=re.IGNORECASE,
)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_stream(handle: Any, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    for block in iter(lambda: handle.read(chunk_size), b""):
        digest.update(block)
    return digest.hexdigest()


def _sha256_file(path: Path) -> str:
    with path.open("rb") as handle:
        return _sha256_stream(handle)


def _decode_text(data: bytes) -> str:
    """Decode release text files, including the UTF-16 readme variant."""

    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16", errors="replace")
    for encoding in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("latin-1", errors="replace")


def _normalise_text(value: Any) -> str:
    if value is None:
        return ""
    text = unicodedata.normalize("NFKC", str(value)).strip()
    return text.replace("\u00a0", " ")


def _normalise_header(value: Any) -> str:
    text = _normalise_text(value).lower()
    text = text.replace("°", " deg ").replace("µ", "u").replace("μ", "u")
    text = re.sub(r"[†*+]+", "", text)
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def _clean_scalar(value: Any) -> Any:
    """Convert spreadsheet scalars to JSON/data-frame friendly values."""

    if value is None:
        return None
    if isinstance(value, (datetime, date, time)):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        text = _normalise_text(value)
        return None if text in _MISSING_MARKERS else text
    return value


def _identifier(value: Any) -> str | None:
    text = _normalise_text(value)
    if not text or text.lower() in _MISSING_MARKERS:
        return None
    match = _WELL_RE.search(text)
    if match:
        return match.group(0).upper()
    return None


def _usgs_site_id(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, float) and value.is_integer():
        text = str(int(value))
    else:
        text = _normalise_text(value)
        if text.endswith(".0") and text[:-2].isdigit():
            text = text[:-2]
    text = re.sub(r"\s+", "", text)
    return text if text and text.lower() not in _MISSING_MARKERS else None


def parse_usgs_site_coordinates(value: Any) -> dict[str, Any]:
    """Decode the usual USGS site-number DMS prefix when it is unambiguous.

    USGS site identifiers in these tables encode latitude as DDMMSS and
    longitude as DDDMMSS.  The parser returns ``None`` rather than guessing
    when a suffix, non-numeric identifier, or invalid DMS component is found.
    """

    site_id = _usgs_site_id(value)
    result: dict[str, Any] = {
        "latitude_dd": None,
        "longitude_dd": None,
        "coordinate_parse_status": "unknown",
    }
    if site_id is None:
        return result
    digits = re.fullmatch(r"(\d{6})(\d{7})(?:\d+)?", site_id)
    if not digits:
        result["coordinate_parse_status"] = "unparsed_site_id"
        return result
    lat_dms, lon_dms = digits.groups()
    lat_deg, lat_min, lat_sec = int(lat_dms[:2]), int(lat_dms[2:4]), int(lat_dms[4:6])
    lon_deg, lon_min, lon_sec = int(lon_dms[:3]), int(lon_dms[3:5]), int(lon_dms[5:7])
    if lat_min >= 60 or lat_sec >= 60 or lon_min >= 60 or lon_sec >= 60:
        result["coordinate_parse_status"] = "invalid_dms"
        return result
    result.update(
        {
            "latitude_dd": lat_deg + lat_min / 60.0 + lat_sec / 3600.0,
            "longitude_dd": -(lon_deg + lon_min / 60.0 + lon_sec / 3600.0),
            "coordinate_parse_status": "parsed_from_site_id",
        }
    )
    return result


def _parse_date(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (datetime, date)):
        return value.date().isoformat() if isinstance(value, datetime) else value.isoformat()
    text = _normalise_text(value)
    if not text or text.lower() in _MISSING_MARKERS:
        return None
    parsed = pd.to_datetime(text, errors="coerce")
    if pd.isna(parsed):
        return None
    return parsed.date().isoformat()


def _parse_time(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        return value.time().replace(microsecond=0).isoformat()
    if isinstance(value, time):
        return value.replace(microsecond=0).isoformat()
    text = _normalise_text(value)
    if not text or text.lower() in _MISSING_MARKERS:
        return None
    parsed = pd.to_datetime(text, errors="coerce")
    if pd.isna(parsed):
        return None
    return parsed.time().replace(microsecond=0).isoformat()


def _combine_datetime(date_value: Any, time_value: Any) -> str | None:
    parsed_date = _parse_date(date_value)
    if parsed_date is None:
        return None
    parsed_time = _parse_time(time_value)
    return f"{parsed_date}T{parsed_time}" if parsed_time else parsed_date


def _number(value: Any) -> float | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    match = _NUMBER_RE.search(_normalise_text(value).replace(",", ""))
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def parse_qualified_value(value: Any) -> dict[str, Any]:
    """Parse a numeric, censored, interval, ``C``, or ``NP`` value.

    ``value_numeric`` is populated only for an exact numeric observation.  A
    reporting-limit value such as ``<0.04`` is represented by ``value_upper``
    and is never silently treated as zero.  ``C`` means CFC above
    air--water equilibrium; it is a source-quality flag, not a numeric
    concentration.  ``NP`` means that the source says dating was not
    possible.
    """

    raw = None if value is None else _normalise_text(value)
    result: dict[str, Any] = {
        "raw_value": raw,
        "value_numeric": None,
        "value_lower": None,
        "value_upper": None,
        "qualifier": None,
        "status": "missing",
        "flags": [],
    }
    if raw is None or raw.lower() in _MISSING_MARKERS:
        result["flags"] = ["unknown_or_missing"]
        return result
    upper = raw.upper().strip()
    if upper == "NP" or "NOT POSSIBLE" in upper:
        result.update(status="not_possible", flags=["NP", "age_not_identifiable"])
        return result
    if upper == "C":
        result.update(status="above_equilibrium", flags=["C", "non_atmospheric_source_possible"])
        return result
    interval = re.fullmatch(
        r"\s*(?:[<>=~]\s*)?([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*[-–—]\s*"
        r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s*",
        raw,
    )
    if interval:
        lower, higher = float(interval.group(1)), float(interval.group(2))
        if lower > higher:
            lower, higher = higher, lower
        result.update(
            value_lower=lower,
            value_upper=higher,
            value_numeric=(lower + higher) / 2.0,
            status="interval",
            flags=["interval_value"],
        )
        return result
    comparator = re.match(r"^\s*([<>]=?|~=?)\s*", raw)
    qualifier = comparator.group(1) if comparator else None
    numeric = _number(raw)
    if numeric is None:
        result["status"] = "unparsed"
        result["flags"] = ["unknown_value_format"]
        return result
    result["qualifier"] = qualifier
    if qualifier and qualifier.startswith("<"):
        result.update(value_upper=numeric, status="below_reporting_limit", flags=["left_censored"])
    elif qualifier and qualifier.startswith(">"):
        result.update(value_lower=numeric, status="above_reporting_limit", flags=["right_censored"])
    elif qualifier and qualifier.startswith("~"):
        result.update(value_numeric=numeric, value_lower=numeric, value_upper=numeric, status="approximate", flags=["approximate_value"])
    else:
        result.update(value_numeric=numeric, value_lower=numeric, value_upper=numeric, status="observed")
    return result


def parse_apparent_age_label(label: Any) -> dict[str, Any]:
    """Turn a CFC recharge-date label into an explicit interval.

    The source uses qualitative labels such as ``Early 1980's`` and
    ``Late 1960's to early 1970's``.  We retain the original label and expose
    a broad calendar interval.  These bounds are a representation of the
    published label, not a newly inferred age uncertainty model.
    """

    raw = None if label is None else _normalise_text(label)
    result: dict[str, Any] = {
        "apparent_age_label": raw,
        "apparent_age_kind": "unknown",
        "recharge_year_min": None,
        "recharge_year_max": None,
        "recharge_year_midpoint": None,
        "age_parse_status": "missing" if not raw else "unknown",
        "age_flags": [],
    }
    if not raw or raw.lower() in _MISSING_MARKERS or raw.upper() in {"NP", "C"}:
        result["age_flags"] = ["age_unknown"]
        return result
    tokens = list(_AGE_TOKEN_RE.finditer(raw))
    if not tokens:
        years = re.findall(r"\b(\d{4})\b", raw)
        if len(years) == 1:
            year = int(years[0])
            result.update(
                apparent_age_kind="point",
                recharge_year_min=year,
                recharge_year_max=year,
                recharge_year_midpoint=float(year),
                age_parse_status="point",
            )
            return result
        result["age_flags"] = ["age_label_unparsed"]
        return result

    bounds: list[tuple[int, int]] = []
    for token in tokens:
        decade = int(token.group("decade"))
        modifier = token.group("modifier").lower()
        if modifier == "early":
            bounds.append((decade, decade + 3))
        elif modifier == "mid":
            bounds.append((decade + 3, decade + 7))
        else:
            bounds.append((decade + 7, decade + 10))
    low = min(item[0] for item in bounds)
    high = max(item[1] for item in bounds)
    result.update(
        apparent_age_kind="interval" if len(tokens) > 1 or bounds[0][0] != bounds[0][1] else "point",
        recharge_year_min=low,
        recharge_year_max=high,
        recharge_year_midpoint=(low + high) / 2.0,
        age_parse_status="interval",
        age_flags=["qualitative_recharge_date_interval"],
    )
    return result


def parse_modelgeoref(text: str) -> dict[str, Any]:
    """Parse the release's four-corner geographic reference file."""

    datum_match = re.search(r"^#\s*Datum\s+(.+?)\s*$", text, flags=re.IGNORECASE | re.MULTILINE)
    corners: dict[str, dict[str, float]] = {}
    pattern = re.compile(
        r"^\s*(upper_left|upper_right|lower_right|lower_left)\s+"
        r"([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s*$",
        flags=re.IGNORECASE | re.MULTILINE,
    )
    for match in pattern.finditer(text):
        name = match.group(1).lower()
        corners[name] = {"longitude": float(match.group(2)), "latitude": float(match.group(3))}
    result: dict[str, Any] = {
        "datum": datum_match.group(1).strip() if datum_match else None,
        "coordinate_reference_system": None,
        "corners": corners,
        "status": "COMPLETE" if len(corners) == 4 else ("PARTIAL" if corners else "MISSING"),
        "source_kind": "model_domain_corner_coordinates",
    }
    if result["datum"]:
        result["coordinate_reference_system"] = f"{result['datum']} geographic; EPSG not declared"
    return result


_MODPATH5_PATHLINE_COLUMNS = [
    "run_id",
    "simulation_id",
    "well_id",
    "source_archive",
    "source_member",
    "file_type",
    "header_text",
    "reference_time",
    "binary_endian",
    "binary_parser",
    "binary_parser_version",
    "record_bytes",
    "ipcode_encoding",
    "record_ordinal",
    "particle_ordinal",
    "particle_ordinal_basis",
    "tracking_direction_code",
    "tracking_direction",
    "particle_id",
    "particle_id_status",
    "particle_id_derived",
    "x",
    "y",
    "zloc",
    "z",
    "time",
    "node",
    "cumulative_timestep",
    "x_global",
    "y_global",
    "z_local",
    "z_global",
    "global_node",
    "time_raw",
    "tracking_time",
    "time_is_intermediate_requested_point",
    "particle_censoring_status",
    "ipcode",
    "ipcode_idcode",
    "ipcode_nslast",
    "termination_status",
    "reference_type",
    "binary_parse_status",
]

_MODPATH5_ENDPOINT_COLUMNS = [
    "run_id",
    "simulation_id",
    "well_id",
    "source_archive",
    "source_member",
    "file_type",
    "header_text",
    "reference_time",
    "binary_endian",
    "binary_parser",
    "binary_parser_version",
    "record_bytes",
    "ipcode_encoding",
    "record_ordinal",
    "particle_ordinal",
    "particle_ordinal_basis",
    "tracking_direction_code",
    "tracking_direction",
    "particle_id",
    "particle_id_status",
    "particle_id_derived",
    "final_zone",
    "final_node",
    "final_x",
    "final_y",
    "final_z_local",
    "total_tracking_time",
    "release_x",
    "release_y",
    "release_z_local",
    "release_node",
    "release_zone",
    "release_cumulative_timestep",
    "ipcode",
    "release_time",
    # Common MODPATH field aliases are retained for callers that already use
    # the text/Flopy readers.  z/z0 are deliberately documented as local
    # coordinates below; binary MODPATH 5 does not contain global endpoint z.
    "x0",
    "y0",
    "z0",
    "zloc0",
    "x",
    "y",
    "z",
    "zloc",
    "time",
    "initial_cell",
    "final_cell",
    "status",
    "ipcode_idcode",
    "ipcode_nslast",
    "termination_status",
    "particle_censoring_status",
    "reference_type",
    "binary_parse_status",
]


def _coerce_binary_payload(payload: bytes | bytearray | memoryview | Path | str) -> bytes:
    """Return bytes for a binary fixture or file path without decoding text."""

    if isinstance(payload, (bytes, bytearray, memoryview)):
        return bytes(payload)
    return Path(payload).expanduser().read_bytes()


def _validate_modpath5_binary_header(
    payload: bytes,
    *,
    record_bytes: int,
    source_member: str | None = None,
) -> tuple[str, float]:
    """Validate the fixed MODPATH 5 header and return text plus ``TREF``.

    MODPATH 5 writes an 80-byte character field followed by one little-endian
    IEEE-754 single precision reference time.  This check is intentionally
    strict about the version and reference-time representation.  In
    particular, a file with a valid-looking prefix but a shifted payload is
    rejected before any particle record is interpreted.
    """

    label = f" for {source_member}" if source_member else ""
    if len(payload) < MODPATH5_HEADER_BYTES:
        raise ValueError(
            f"MODPATH 5 binary file{label} is truncated: expected at least "
            f"{MODPATH5_HEADER_BYTES} bytes, got {len(payload)}"
        )
    header_bytes = payload[:_MODPATH5_HEADER_TEXT_BYTES]
    header_text = header_bytes.decode("ascii", errors="replace").rstrip("\x00 ")
    if not re.match(r"^MODPATH\s+5\.0(?:\s|$)", header_text, flags=re.IGNORECASE):
        raise ValueError(
            f"MODPATH 5 binary file{label} has an unsupported header: "
            f"{header_text!r}"
        )
    reference_time = struct.unpack_from("<f", payload, _MODPATH5_HEADER_TEXT_BYTES)[0]
    if not math.isfinite(reference_time):
        raise ValueError(
            f"MODPATH 5 binary file{label} has a non-finite reference time "
            f"({reference_time!r})"
        )
    data_bytes = len(payload) - MODPATH5_HEADER_BYTES
    if data_bytes % record_bytes:
        raise ValueError(
            f"MODPATH 5 binary file{label} has an invalid payload length: "
            f"{data_bytes} bytes after the {MODPATH5_HEADER_BYTES}-byte header "
            f"is not divisible by {record_bytes}"
        )
    return header_text, float(reference_time)


def decode_modpath5_ipcode(ipcode: int | None) -> dict[str, Any]:
    """Decode an MODPATH 5 ``IPCODE`` without discarding its raw integer.

    ``IPCODE`` in the Aiken MODPATH 5.0 output uses the encoded value
    ``10 * NSLAST + IDCODE`` for positive values.  This interpretation is
    checked against the release's ``.sum`` summaries (for example, 10 and 540
    are active particles, whereas 51/41/31/21 are normal terminations), rather
    than inferred from endpoint row order.  The returned censoring labels are
    deliberately conservative: they describe whether an observed termination
    is available, not an independently measured field travel-time truth.
    """

    result: dict[str, Any] = {
        "ipcode_idcode": None,
        "ipcode_nslast": None,
        "termination_status": "unknown_ipcode",
        "particle_censoring_status": "unknown",
    }
    if ipcode is None:
        return result
    value = int(ipcode)
    if value == -2:
        result.update(
            ipcode_idcode=-2,
            termination_status="unreleased",
            particle_censoring_status="unreleased",
        )
        return result
    if value == -1:
        result.update(
            ipcode_idcode=-1,
            termination_status="stranded_inactive_dry_cell",
            particle_censoring_status="censored_dry_cell",
        )
        return result
    if value == 0:
        result.update(
            ipcode_idcode=0,
            termination_status="active_at_stop_time",
            particle_censoring_status="right_censored_active",
        )
        return result
    if value > 0:
        idcode = value % 10
        nslast = value // 10
        status = {
            0: "active_at_stop_time",
            1: "discharged_normally",
            2: "stopped_in_specified_zone",
        }.get(idcode, "unknown_positive_idcode")
        result.update(
            ipcode_idcode=idcode,
            ipcode_nslast=nslast,
            termination_status=status,
            particle_censoring_status=(
                "observed_termination" if idcode in {1, 2} else "unknown"
            ),
        )
        if idcode == 0:
            result["particle_censoring_status"] = "right_censored_active"
    return result


def _binary_run_id(source_member: str | None, fallback: str | None = None) -> str | None:
    """Derive the archive run directory while preserving the member path."""

    if source_member:
        parts = source_member.replace("\\", "/").split("/")
        for part in reversed(parts[:-1]):
            if part.lower().startswith("output."):
                return part[len("output.") :]
            if part.lower().endswith("_mp"):
                return part
        if len(parts) > 1:
            return parts[-2]
    return fallback


def _binary_well_id(source_member: str | None) -> str | None:
    if not source_member:
        return None
    match = _WELL_RE.search(Path(source_member.replace("\\", "/")).stem)
    return match.group(0).upper() if match else None


def _parse_modpath5_tracking_direction(text: str) -> tuple[int | None, str]:
    """Read the direction choice from a MODPATH response file when present.

    Aiken contains one forward run (AK-831) among backward-to-recharge runs.
    The response file is the authoritative run configuration; particle
    coordinates alone cannot distinguish the two directions.  Keep unknown
    direction explicit rather than inferring it from endpoint geometry.
    """

    lines = text.splitlines()
    for index, line in enumerate(lines):
        if "IN WHICH DIRECTION SHOULD PARTICLES BE TRACKED" not in line.upper():
            continue
        for candidate in lines[index + 1 : index + 14]:
            value = candidate.strip()
            if not value or value.startswith(("@", "*")):
                continue
            if value in {"1", "2"}:
                code = int(value)
                return code, MODPATH5_TRACKING_DIRECTIONS[code]
        break
    return None, "unknown_tracking_direction"


def _model_tracking_directions(source: "_AikenReader") -> dict[str, tuple[int | None, str]]:
    """Return response-file tracking directions keyed by MODPATH run ID."""

    directions: dict[str, tuple[int | None, str]] = {}
    for member in source.member_names("model"):
        if not member.lower().endswith(".rsp"):
            continue
        try:
            payload = source.read_member(member)
        except (OSError, KeyError):
            continue
        run_id = _binary_run_id(member)
        if run_id is not None:
            directions[run_id] = _parse_modpath5_tracking_direction(
                payload.decode("latin-1", errors="replace")
            )
    return directions


def parse_modpath5_pathline_binary(
    payload: bytes | bytearray | memoryview | Path | str,
    *,
    run_id: str | None = None,
    source_archive: str | None = None,
    source_member: str | None = None,
    tracking_direction_code: int | None = None,
    tracking_direction: str | None = None,
) -> pd.DataFrame:
    """Parse a MODPATH 5 binary pathline (``.pth``) payload.

    The eight fields are ``particle, x, y, zloc, z, time, node,
    cumulative_timestep``.  The raw signed time is preserved because MODPATH
    uses a negative time flag for requested intermediate points; the absolute
    physical tracking time is exposed separately as ``tracking_time``.
    """

    data = _coerce_binary_payload(payload)
    header_text, reference_time = _validate_modpath5_binary_header(
        data,
        record_bytes=MODPATH5_PATHLINE_RECORD_BYTES,
        source_member=source_member,
    )
    run = run_id or _binary_run_id(source_member)
    well_id = _binary_well_id(source_member)
    particle_ordinals: dict[int, int] = {}
    rows: list[dict[str, Any]] = []
    for ordinal, offset in enumerate(
        range(MODPATH5_HEADER_BYTES, len(data), MODPATH5_PATHLINE_RECORD_BYTES),
        start=1,
    ):
        particle_id, x, y, zloc, z, time_value, node, cumulative_timestep = (
            _MODPATH5_PATHLINE_STRUCT.unpack_from(data, offset)
        )
        if particle_id not in particle_ordinals:
            particle_ordinals[particle_id] = len(particle_ordinals) + 1
        time_raw = float(time_value)
        row = {
            "run_id": run,
            "simulation_id": run,
            "well_id": well_id,
            "source_archive": source_archive,
            "source_member": source_member,
            "file_type": "MODPATH5_BINARY_PATHLINE",
            "header_text": header_text,
            "reference_time": reference_time,
            "binary_endian": MODPATH5_BINARY_ENDIAN,
            "binary_parser": MODPATH5_BINARY_PARSER,
            "binary_parser_version": MODPATH5_BINARY_PARSER_VERSION,
            "record_bytes": MODPATH5_PATHLINE_RECORD_BYTES,
            "ipcode_encoding": MODPATH5_IPCODE_ENCODING,
            "record_ordinal": ordinal,
            "particle_ordinal": particle_ordinals[particle_id],
            "particle_ordinal_basis": "first_seen_native_pathline_particle_id",
            "tracking_direction_code": tracking_direction_code,
            "tracking_direction": tracking_direction or "unknown_tracking_direction",
            "particle_id": int(particle_id),
            "particle_id_status": "native_pathline_particle_index",
            "particle_id_derived": False,
            "x": float(x),
            "y": float(y),
            "zloc": float(zloc),
            "z": float(z),
            "time": time_raw,
            "node": int(node),
            "cumulative_timestep": int(cumulative_timestep),
            "x_global": float(x),
            "y_global": float(y),
            "z_local": float(zloc),
            "z_global": float(z),
            "global_node": int(node),
            "time_raw": time_raw,
            "tracking_time": abs(time_raw),
            "time_is_intermediate_requested_point": time_raw < 0.0,
            "particle_censoring_status": "unknown_no_endpoint_record",
            "ipcode": None,
            "ipcode_idcode": None,
            "ipcode_nslast": None,
            "termination_status": "unknown_no_endpoint_record",
            "reference_type": REFERENCE_TYPE,
            "binary_parse_status": "validated",
        }
        rows.append(row)
    return pd.DataFrame(rows, columns=_MODPATH5_PATHLINE_COLUMNS)


def parse_modpath5_endpoint_binary(
    payload: bytes | bytearray | memoryview | Path | str,
    *,
    run_id: str | None = None,
    source_archive: str | None = None,
    source_member: str | None = None,
    tracking_direction_code: int | None = None,
    tracking_direction: str | None = None,
) -> pd.DataFrame:
    """Parse a MODPATH 5 binary endpoint (``.ept``) payload.

    MODPATH 5 endpoint records omit the global final and release ``z``
    coordinates.  The local coordinates are therefore retained under both
    their release names (``release_z_local``/``final_z_local``) and the
    conventional ``zloc0``/``zloc`` aliases; ``z`` and the global release
    coordinate remain explicitly missing rather than reconstructed.
    """

    data = _coerce_binary_payload(payload)
    header_text, reference_time = _validate_modpath5_binary_header(
        data,
        record_bytes=MODPATH5_ENDPOINT_RECORD_BYTES,
        source_member=source_member,
    )
    run = run_id or _binary_run_id(source_member)
    well_id = _binary_well_id(source_member)
    rows: list[dict[str, Any]] = []
    for ordinal, offset in enumerate(
        range(MODPATH5_HEADER_BYTES, len(data), MODPATH5_ENDPOINT_RECORD_BYTES),
        start=1,
    ):
        (
            final_zone,
            final_node,
            final_x,
            final_y,
            final_z_local,
            total_tracking_time,
            release_x,
            release_y,
            release_z_local,
            release_node,
            release_zone,
            release_cumulative_timestep,
            ipcode,
            release_time,
        ) = _MODPATH5_ENDPOINT_STRUCT.unpack_from(data, offset)
        ipcode_fields = decode_modpath5_ipcode(int(ipcode))
        rows.append(
            {
                "run_id": run,
                "simulation_id": run,
                "well_id": well_id,
                "source_archive": source_archive,
                "source_member": source_member,
                "file_type": "MODPATH5_BINARY_ENDPOINT",
                "header_text": header_text,
                "reference_time": reference_time,
                "binary_endian": MODPATH5_BINARY_ENDIAN,
                "binary_parser": MODPATH5_BINARY_PARSER,
                "binary_parser_version": MODPATH5_BINARY_PARSER_VERSION,
                "record_bytes": MODPATH5_ENDPOINT_RECORD_BYTES,
                "ipcode_encoding": MODPATH5_IPCODE_ENCODING,
                "record_ordinal": ordinal,
                "particle_ordinal": ordinal,
                "particle_ordinal_basis": "endpoint_record_row_order_crosswalk_only",
                "tracking_direction_code": tracking_direction_code,
                "tracking_direction": tracking_direction or "unknown_tracking_direction",
                # MODPATH 5 binary endpoint records contain no particle ID.
                # Keep the row ordinal as the explicit join key and leave the
                # native ID missing rather than presenting a derived ordinal
                # as a source identifier.
                "particle_id": None,
                "particle_id_status": "not_present_in_modpath5_endpoint",
                "particle_id_derived": False,
                "final_zone": int(final_zone),
                "final_node": int(final_node),
                "final_x": float(final_x),
                "final_y": float(final_y),
                "final_z_local": float(final_z_local),
                "total_tracking_time": float(total_tracking_time),
                "release_x": float(release_x),
                "release_y": float(release_y),
                "release_z_local": float(release_z_local),
                "release_node": int(release_node),
                "release_zone": int(release_zone),
                "release_cumulative_timestep": int(release_cumulative_timestep),
                "ipcode": int(ipcode),
                "release_time": float(release_time),
                "x0": float(release_x),
                "y0": float(release_y),
                "z0": float(release_z_local),
                "zloc0": float(release_z_local),
                "x": float(final_x),
                "y": float(final_y),
                "z": None,
                "zloc": float(final_z_local),
                "time": float(total_tracking_time),
                "initial_cell": int(release_node),
                "final_cell": int(final_node),
                "status": int(ipcode),
                **ipcode_fields,
                "reference_type": REFERENCE_TYPE,
                "binary_parse_status": "validated",
            }
        )
    return pd.DataFrame(rows, columns=_MODPATH5_ENDPOINT_COLUMNS)


def parse_readme(text: str) -> dict[str, Any]:
    """Extract non-interpretive release metadata from ``readme.txt``."""

    dois = []
    for match in _DOI_RE.findall(text):
        clean = match.rstrip(".,;)")
        if clean.lower() not in {item.lower() for item in dois}:
            dois.append(clean)
    def _archive_date(label: str) -> str | None:
        match = re.search(rf"{label}\s*:\s*(\d{{4}}-\d{{2}}-\d{{2}})", text, flags=re.IGNORECASE)
        return match.group(1) if match else None

    number_words = {
        "one": 1,
        "two": 2,
        "three": 3,
        "four": 4,
        "five": 5,
        "six": 6,
        "seven": 7,
        "eight": 8,
        "nine": 9,
        "ten": 10,
        "eleven": 11,
        "twelve": 12,
        "thirteen": 13,
        "fourteen": 14,
        "fifteen": 15,
        "sixteen": 16,
        "seventeen": 17,
        "eighteen": 18,
        "nineteen": 19,
        "twenty": 20,
    }

    def _count(pattern: str) -> int | None:
        # The release uses both numerals and words (for example, "six
        # MODFLOW-NWT simulations"), so preserve either representation.
        word_pattern = pattern.replace(r"(\d+)", r"(\d+|" + "|".join(number_words) + ")")
        match = re.search(word_pattern, text, flags=re.IGNORECASE)
        if not match:
            return None
        value = match.group(1).lower()
        return int(value) if value.isdigit() else number_words.get(value)

    return {
        "publication_doi": next((item for item in dois if item.lower() == PUBLICATION_DOI.lower()), None),
        "data_release_doi": next((item for item in dois if item.lower() == SOURCE_DOI.lower()), None),
        "dois": dois,
        "archive_created": _archive_date("Archive created"),
        "archive_updated": _archive_date("Archive updated"),
        "archive_released": _archive_date("Archive released"),
        "n_modflow_simulations": _count(r"(\d+)\s+MODFLOW-NWT\s+simulations"),
        "n_modpath_simulations": _count(r"(\d+)\s+MODPATH\s+simulations"),
        "n_public_supply_wells": _count(r"(\d+)\s+PSWs\b"),
        "n_monitoring_wells": _count(r"(\d+)\s+(?:existing\s+)?monitoring wells"),
        "model_families": [
            name for name in ("MODFLOW-NWT", "MODPATH5", "PEST") if name.lower() in text.lower()
        ],
        "age_and_pathway_statement_present": bool(
            re.search(r"CFC.*age.*particle-tracking|particle-tracking.*CFC", text, re.I | re.S)
        ),
        "raw_text_sha256": _sha256_bytes(text.encode("utf-8")),
    }


def _find_header_start(rows: Sequence[Sequence[Any]]) -> int | None:
    for index, row in enumerate(rows):
        labels = {_normalise_header(value) for value in row if _normalise_text(value)}
        if "well_id" in labels and (
            "sample_date" in labels
            or "county_number_ak_n" in labels
            or "county_number_for_well_ak_n_or_lex_n" in labels
        ):
            return index
    return None


def _header_end(sheet_name: str, rows: Sequence[Sequence[Any]], start: int) -> int:
    if _normalise_text(sheet_name).lower() == "table 13":
        # Table 13 has a three-row grouped header.  Include the row containing
        # CFC-11/12/113 labels but do not absorb the section labels below it.
        for index in range(start, min(len(rows), start + 4)):
            text = " ".join(_normalise_header(value) for value in rows[index])
            if "cfc_11" in text and "cfc_12" in text:
                return index + 1
    return start + 1


def _unique_headers(headers: Sequence[Any]) -> list[str]:
    output: list[str] = []
    counts: dict[str, int] = {}
    for index, value in enumerate(headers):
        name = _normalise_header(value) or f"unnamed_{index + 1}"
        counts[name] = counts.get(name, 0) + 1
        output.append(name if counts[name] == 1 else f"{name}_{counts[name]}")
    return output


@dataclass
class _ParsedSheet:
    workbook_name: str
    sheet_name: str
    header_start: int | None
    header_end: int | None
    headers: list[str]
    raw_headers: list[Any]
    rows: list[dict[str, Any]]


def _read_workbook_sheets(data: bytes, workbook_name: str) -> dict[str, _ParsedSheet]:
    workbook = openpyxl.load_workbook(BytesIO(data), read_only=True, data_only=True)
    output: dict[str, _ParsedSheet] = {}
    try:
        for sheet_name in workbook.sheetnames:
            worksheet = workbook[sheet_name]
            matrix = [list(row) for row in worksheet.iter_rows(values_only=True)]
            start = _find_header_start(matrix)
            if start is None:
                output[sheet_name] = _ParsedSheet(workbook_name, sheet_name, None, None, [], [], [])
                continue
            end = _header_end(sheet_name, matrix, start)
            width = max((len(row) for row in matrix[start:end]), default=0)
            raw_headers: list[Any] = []
            for column in range(width):
                value = None
                for row in matrix[start:end]:
                    if column < len(row) and _normalise_text(row[column]):
                        value = row[column]
                raw_headers.append(value)
            headers = _unique_headers(raw_headers)
            rows: list[dict[str, Any]] = []
            role = "unknown"
            organization: str | None = None
            for row_number, row in enumerate(matrix[end:], start=end + 1):
                values = list(row) + [None] * max(0, width - len(row))
                row_text = " ".join(_normalise_text(value) for value in values if _normalise_text(value))
                lower = row_text.lower()
                if "monitoring wells" in lower:
                    role = "monitoring_well"
                elif "public-supply wells" in lower or "public supply wells" in lower:
                    role = "public_supply_well"
                elif lower.strip() == "surface-water sample" or "surface-water" in lower:
                    role = "surface_water"
                candidate_ids = [_identifier(value) for value in values]
                well_id = next((item for item in candidate_ids if item), None)
                sample_date_index = next(
                    (idx for idx, header in enumerate(headers) if header == "sample_date"), None
                )
                has_sample_date = sample_date_index is not None and _parse_date(values[sample_date_index]) is not None
                is_surface = (
                    has_sample_date
                    and (
                        ("surface" in lower and "water" in lower)
                        or any(_normalise_text(value).upper() == "SW" for value in values)
                    )
                )
                if not well_id and not is_surface:
                    # Section headings, units, MRL/MCL rows and empty spacer
                    # rows are intentionally ignored.
                    if row_text and not lower.startswith(("mrl", "mcl")) and not lower.startswith("public-supply wells"):
                        organization = row_text
                    continue
                if not well_id and is_surface:
                    role = "surface_water"
                if not well_id and is_surface:
                    # Preserve surface samples in chemistry but keep their
                    # node identity explicitly unknown.
                    well_id = None
                if well_id is not None and not any(_identifier(value) for value in values):
                    continue
                # An organization is a non-data line between a role heading
                # and the first data row.  Keep it only when no identifier is
                # present; this avoids treating well descriptions as utilities.
                if well_id is None and row_text:
                    organization = row_text
                rows.append(
                    {
                        "row_number": row_number,
                        "values": values,
                        "well_id": well_id,
                        "site_role": role,
                        "organization": organization,
                    }
                )
            output[sheet_name] = _ParsedSheet(
                workbook_name,
                sheet_name,
                start,
                end,
                headers,
                raw_headers,
                rows,
            )
    finally:
        workbook.close()
    return output


def _col_index(sheet: _ParsedSheet, *names: str) -> int | None:
    normalised = {_normalise_header(name) for name in names}
    for index, header in enumerate(sheet.headers):
        if header in normalised:
            return index
    for index, header in enumerate(sheet.headers):
        if any(name in header for name in normalised):
            return index
    return None


def _row_value(sheet: _ParsedSheet, row: Mapping[str, Any], *names: str) -> Any:
    index = _col_index(sheet, *names)
    if index is None:
        return None
    values = row["values"]
    return values[index] if index < len(values) else None


def _base_row(sheet: _ParsedSheet, row: Mapping[str, Any]) -> dict[str, Any]:
    well_id = row.get("well_id")
    date_value = _row_value(sheet, row, "sample_date")
    time_value = _row_value(sheet, row, "sample_time")
    site_id = _usgs_site_id(
        _row_value(sheet, row, "usgs_site_id_number", "usgs_site_id")
    )
    result: dict[str, Any] = {
        "node_id": f"AIKEN:{well_id}" if well_id else None,
        "well_id": well_id,
        "county_number": _normalise_text(
            _row_value(
                sheet,
                row,
                "county_number_ak_n",
                "county_number_for_well_ak_n_or_lex_n",
            )
        )
        or None,
        "usgs_site_id": site_id,
        "sample_date": _parse_date(date_value),
        "sample_time": _parse_time(time_value),
        "sample_datetime": _combine_datetime(date_value, time_value),
        "site_role": row.get("site_role", "unknown"),
        "organization": row.get("organization"),
        "reference_type": REFERENCE_TYPE,
        "source_workbook": sheet.workbook_name,
        "source_sheet": sheet.sheet_name,
        "source_row": row.get("row_number"),
        "adjacency_truth_status": "ABSTAIN",
        "reaction_truth_status": "ABSTAIN",
    }
    result.update(parse_usgs_site_coordinates(site_id))
    return result


def _empty_frame(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame({column: pd.Series(dtype="object") for column in columns})


def _make_well_metadata(sheet: _ParsedSheet | None) -> pd.DataFrame:
    columns = [
        "node_id",
        "well_id",
        "well_description",
        "county_number",
        "usgs_site_id",
        "year_installed",
        "land_surface_altitude_ft_ngvd29",
        "total_depth_ft_bls",
        "pump_type",
        "aquifer_screened",
        "site_role",
        "organization",
        "latitude_dd",
        "longitude_dd",
        "coordinate_parse_status",
        "screen_interval_status",
        "reference_type",
        "adjacency_truth_status",
        "reaction_truth_status",
        "source_workbook",
        "source_sheet",
        "source_row",
    ]
    if sheet is None:
        return _empty_frame(columns)
    records = []
    for row in sheet.rows:
        base = _base_row(sheet, row)
        if not base["well_id"]:
            continue
        record = dict(base)
        record.update(
            {
                "well_description": _clean_scalar(_row_value(sheet, row, "well_id")),
                "year_installed": _number(_row_value(sheet, row, "year_installed")),
                "land_surface_altitude_ft_ngvd29": _number(
                    _row_value(sheet, row, "altitude_of_land_surface_feet_ngvd_29", "altitude_of_land_surface_feet_ngvd_29")
                ),
                "total_depth_ft_bls": _number(
                    _row_value(sheet, row, "total_depth_of_completed_well_feet_below_land_surface_altitude")
                ),
                "pump_type": _clean_scalar(_row_value(sheet, row, "pump_type")),
                "aquifer_screened": _clean_scalar(
                    _row_value(sheet, row, "aquifer_screened_by_well_open_hole_if_bedrock")
                ),
                "screen_interval_status": "unknown_total_depth_only",
            }
        )
        records.append(record)
    frame = pd.DataFrame(records)
    for column in columns:
        if column not in frame:
            frame[column] = pd.NA
    return frame[columns].drop_duplicates(subset=["well_id"], keep="first").reset_index(drop=True)


def _make_field_samples(sheet: _ParsedSheet | None) -> pd.DataFrame:
    columns = [
        "node_id", "well_id", "county_number", "usgs_site_id", "sample_date", "sample_time",
        "sample_datetime", "sample_collection_method", "water_temperature_c",
        "specific_conductance_us_cm", "ph", "dissolved_oxygen_mg_l",
        "dissolved_oxygen_pct_saturation", "comments", "site_role", "reference_type",
        "source_workbook", "source_sheet", "source_row", "adjacency_truth_status", "reaction_truth_status",
    ]
    if sheet is None:
        return _empty_frame(columns)
    records = []
    for row in sheet.rows:
        base = _base_row(sheet, row)
        if not base["well_id"]:
            continue
        record = dict(base)
        record.update(
            {
                "sample_collection_method": _clean_scalar(_row_value(sheet, row, "sample_collection_method")),
                "water_temperature_c": _number(
                    _row_value(sheet, row, "water_temperature_c", "water_temperature_deg_c")
                ),
                "specific_conductance_us_cm": _number(_row_value(sheet, row, "specific_conductance_us_cm")),
                "ph": _number(_row_value(sheet, row, "ph")),
                "dissolved_oxygen_mg_l": _number(_row_value(sheet, row, "dissolved_oxygen_mg_l")),
                "dissolved_oxygen_pct_saturation": _number(
                    _row_value(
                        sheet,
                        row,
                        "dissolved_oxygen_pct_saturation_at_sample_temperature",
                        "dissolved_oxygen_saturation_at_sample_temperature",
                    )
                ),
                "comments": _clean_scalar(_row_value(sheet, row, "comments")),
            }
        )
        records.append(record)
    frame = pd.DataFrame(records)
    for column in columns:
        if column not in frame:
            frame[column] = pd.NA
    return frame[columns].reset_index(drop=True)


def _make_cfc_ages(sheet: _ParsedSheet | None) -> pd.DataFrame:
    columns = [
        "node_id", "well_id", "county_number", "usgs_site_id", "sample_date", "sample_time", "sample_datetime",
        "cfc_11_raw", "cfc_11_value_numeric", "cfc_11_status", "cfc_11_flags",
        "cfc_12_raw", "cfc_12_value_numeric", "cfc_12_status", "cfc_12_flags",
        "cfc_113_raw", "cfc_113_value_numeric", "cfc_113_status", "cfc_113_flags",
        "piston_cfc_11_raw", "piston_cfc_12_raw", "piston_cfc_113_raw", "cfc_used_for_ages",
        "apparent_age_label", "apparent_age_kind", "recharge_year_min", "recharge_year_max",
        "recharge_year_midpoint", "age_parse_status", "age_flags", "age_likelihood_status",
        "independent_age_truth", "reference_type", "adjacency_truth_status", "reaction_truth_status",
        "source_workbook", "source_sheet", "source_row",
    ]
    if sheet is None:
        return _empty_frame(columns)
    # Table 13's grouped headers are disambiguated by positional occurrence.
    # The first occurrence of each CFC header is concentration; the second is
    # piston elapsed time.  Unique-header suffixes are stable because the
    # grouped header is read left-to-right.
    def _index_with_suffix(name: str, suffix: int) -> int | None:
        target = name if suffix == 1 else f"{name}_{suffix}"
        try:
            return sheet.headers.index(target)
        except ValueError:
            return None

    records = []
    for row in sheet.rows:
        base = _base_row(sheet, row)
        if not base["well_id"]:
            continue
        values = row["values"]
        def _at(index: int | None) -> Any:
            return values[index] if index is not None and index < len(values) else None

        parsed: dict[str, dict[str, Any]] = {}
        for name in ("cfc_11", "cfc_12", "cfc_113"):
            parsed[name] = parse_qualified_value(_at(_index_with_suffix(name, 1)))
        age = parse_apparent_age_label(_row_value(sheet, row, "assigned_cfc_apparent_groundwater_age_date"))
        flags = set(age.get("age_flags", []))
        for item in parsed.values():
            flags.update(item.get("flags", []))
        usable_concentration = any(
            item["status"] in {"observed", "interval", "approximate"} for item in parsed.values()
        )
        record = dict(base)
        record.update(
            {
                "cfc_11_raw": parsed["cfc_11"]["raw_value"],
                "cfc_11_value_numeric": parsed["cfc_11"]["value_numeric"],
                "cfc_11_status": parsed["cfc_11"]["status"],
                "cfc_11_flags": "|".join(parsed["cfc_11"]["flags"]),
                "cfc_12_raw": parsed["cfc_12"]["raw_value"],
                "cfc_12_value_numeric": parsed["cfc_12"]["value_numeric"],
                "cfc_12_status": parsed["cfc_12"]["status"],
                "cfc_12_flags": "|".join(parsed["cfc_12"]["flags"]),
                "cfc_113_raw": parsed["cfc_113"]["raw_value"],
                "cfc_113_value_numeric": parsed["cfc_113"]["value_numeric"],
                "cfc_113_status": parsed["cfc_113"]["status"],
                "cfc_113_flags": "|".join(parsed["cfc_113"]["flags"]),
                "piston_cfc_11_raw": _clean_scalar(_at(_index_with_suffix("cfc_11", 2))),
                "piston_cfc_12_raw": _clean_scalar(_at(_index_with_suffix("cfc_12", 2))),
                "piston_cfc_113_raw": _clean_scalar(_at(_index_with_suffix("cfc_113", 2))),
                "cfc_used_for_ages": _clean_scalar(_row_value(sheet, row, "cfcs_used_for_ages")),
                **age,
                "age_flags": "|".join(sorted(flags)),
                "age_likelihood_status": "screening_interval_available" if usable_concentration and age["age_parse_status"] in {"interval", "point"} else "ABSTAIN",
                "independent_age_truth": False,
            }
        )
        records.append(record)
    frame = pd.DataFrame(records)
    for column in columns:
        if column not in frame:
            frame[column] = pd.NA
    return frame[columns].reset_index(drop=True)


def _make_pathways(sheet: _ParsedSheet | None) -> pd.DataFrame:
    columns = [
        "node_id", "well_id", "county_number", "usgs_site_id", "groundwater_age_years_before_2015",
        "recharge_date_label", "pathway_extent_ft", "flow_velocity_ft_per_year", "comment",
        "pathway_status", "pathway_reference_kind", "independent_flow_truth", "direct_adjacency_truth_status",
        "reference_type", "source_workbook", "source_sheet", "source_row",
    ]
    if sheet is None:
        return _empty_frame(columns)
    records = []
    for row in sheet.rows:
        base = _base_row(sheet, row)
        if not base["well_id"]:
            continue
        age_raw = _clean_scalar(_row_value(sheet, row, "groundwater_age_elapsed_time_since_recharge_yrs_before_2015"))
        extent_raw = _clean_scalar(_row_value(sheet, row, "groundwater_flow_pathway_extent_from_most_distal_recharge_area_to_the_well_ft"))
        velocity_raw = _clean_scalar(_row_value(sheet, row, "groundwater_flow_velocity_estimated_ft_yr"))
        extent = parse_qualified_value(extent_raw)
        velocity = parse_qualified_value(velocity_raw)
        if extent["status"] == "not_possible" or velocity["status"] == "not_possible":
            status = "not_possible"
        elif extent["value_numeric"] is not None and velocity["value_numeric"] is not None:
            status = "reported"
        else:
            status = "missing_or_partial"
        record = dict(base)
        record.update(
            {
                "groundwater_age_years_before_2015": _number(age_raw),
                "recharge_date_label": _clean_scalar(_row_value(sheet, row, "recharge_date")),
                "pathway_extent_ft": extent["value_numeric"],
                "flow_velocity_ft_per_year": velocity["value_numeric"],
                "comment": _clean_scalar(_row_value(sheet, row, "comment")),
                "pathway_status": status,
                "pathway_reference_kind": "model_derived_recharge_to_well_pathway_summary",
                "independent_flow_truth": False,
                "direct_adjacency_truth_status": "ABSTAIN",
            }
        )
        records.append(record)
    frame = pd.DataFrame(records)
    for column in columns:
        if column not in frame:
            frame[column] = pd.NA
    return frame[columns].reset_index(drop=True)


def _chemical_family(sheet_name: str, parameter: str) -> str:
    """Classify a source column without letting analyte names leak across families.

    The release uses several wide tables.  In particular, bromomethane and
    dibromomethane are VOC analytes in Table 10, not dissolved gases, and the
    isotope headers include their units (``dH, permil``/``dO, permil``).
    Table context therefore takes precedence over substring heuristics.
    """

    sheet_lower = sheet_name.lower()
    lower = f"{sheet_name} {parameter}".lower()
    if sheet_lower.startswith("table 10"):
        return "voc"
    if sheet_lower == "table 8":
        # These are field parameters, even when a label contains "oxygen";
        # they are not the dissolved-gas panel in Table 14.
        return "field_physical_or_other"
    if "radium" in lower or "ra-" in lower:
        return "radium"
    if (
        "isotope" in lower
        or "permil" in lower
        or parameter.lower().strip() in {"dh", "do", "δh", "δo"}
        or parameter.lower().startswith(("dh,", "do,"))
    ):
        return "stable_isotope"
    if "nitrate" in lower:
        return "nitrate"
    if "cfc" in lower or "chlorofluorocarbon" in lower:
        return "cfc"
    if any(
        token in lower
        for token in ("methane", "carbon dioxide", "nitrogen", "oxygen", "argon", "ch4", "co2")
    ):
        return "dissolved_gas"
    return "field_physical_or_other"


def _chemistry_role(sheet_name: str, parameter_key: str, parameter: str) -> str:
    """Return the measurement/derivation role of a canonical chemistry row."""

    sheet_lower = sheet_name.lower()
    key = parameter_key.lower()
    parameter_lower = parameter.lower()
    if sheet_lower.startswith("table 10"):
        return "voc_concentration"
    if sheet_lower == "table 9":
        if "ratio" in key or "/" in parameter:
            return "radium_ratio"
        return "radionuclide_concentration"
    if sheet_lower == "table 11":
        return "stable_isotope"
    if sheet_lower == "table 12":
        return "nitrate_concentration"
    if sheet_lower == "table 13":
        if key in {"cfc_11", "cfc_12", "cfc_113"}:
            return "cfc_concentration"
        if key in {"cfc_11_2", "cfc_12_2", "cfc_113_2"}:
            return "piston_elapsed_years"
        if key.startswith("from_cfc"):
            return "young_water_fraction"
        return "cfc_derived_metric"
    if sheet_lower == "table 14":
        if "temperature" in key or "temperature" in parameter_lower:
            return "field_measurement"
        return "dissolved_gas_concentration"
    if sheet_lower == "table 8":
        return "field_measurement"
    return "analyte_concentration"


def _unit_for(sheet_name: str, parameter: str, role: str | None = None) -> str | None:
    """Resolve units using the source-table role before text heuristics."""

    if role == "radium_ratio":
        return "dimensionless"
    if role == "radionuclide_concentration":
        return "pCi/L"
    if role == "stable_isotope":
        return "permil"
    if role == "cfc_concentration":
        return "pg/kg"
    if role == "piston_elapsed_years":
        return "years"
    if role == "young_water_fraction":
        return "%"
    if role == "voc_concentration":
        return "ug/L"
    if role == "nitrate_concentration" or role == "dissolved_gas_concentration":
        return "mg/L"

    lower = f"{sheet_name} {parameter}".lower()
    if "dissolved oxygen" in lower and "saturation" in lower:
        return "%"
    if "dissolved oxygen" in lower or "dissolved_oxygen" in lower:
        return "mg/L"
    if "temperature" in lower:
        return "deg C"
    if "conductance" in lower:
        return "uS/cm"
    if parameter.lower().strip() == "ph":
        return "pH"
    if "permil" in lower or "isotope" in lower:
        return "permil"
    if "ratio" in lower or "/" in parameter:
        return "dimensionless"
    return None


def _make_chemistry(sheets: Iterable[_ParsedSheet]) -> pd.DataFrame:
    columns = [
        "node_id", "well_id", "county_number", "usgs_site_id", "sample_date", "sample_time", "sample_datetime",
        "parameter", "parameter_key", "parameter_role", "chemical_family", "unit", "raw_value", "value_numeric", "value_lower", "value_upper",
        "qualifier", "measurement_status", "measurement_flags", "site_role", "reference_type",
        "reaction_truth_status", "adjacency_truth_status", "source_workbook", "source_sheet", "source_row",
        "preferred_source_for_duplicate_voc",
    ]
    metadata_headers = {
        "well_id", "usgs_site_id", "usgs_site_id_number", "sample_date", "sample_time",
        "sample_collection_method", "county_number_ak_n", "county_number_ak_n_or_lex_n",
        "county_number_for_well_ak_n_or_lex_n",
        "pump_status", "comments", "aquifer_screened_by_well_open_hole_if_bedrock", "sample_collection_location", "cfc_used_for_ages",
        "cfcs_used_for_ages", "assigned_cfc_apparent_groundwater_age_date", "recharge_date",
        "comment", "groundwater_age_elapsed_time_since_recharge_yrs_before_2015",
        "groundwater_flow_pathway_extent_from_most_distal_recharge_area_to_the_well_ft",
        "groundwater_flow_velocity_estimated_ft_yr", "recharge",
    }
    records = []
    for sheet in sheets:
        if not sheet.rows:
            continue
        for row in sheet.rows:
            base = _base_row(sheet, row)
            values = row["values"]
            for index, header in enumerate(sheet.headers):
                if header in metadata_headers or header.startswith("unnamed_"):
                    continue
                if index >= len(values):
                    continue
                raw = values[index]
                parsed = parse_qualified_value(raw)
                if parsed["status"] == "missing" and raw is None:
                    continue
                parameter = _normalise_text(sheet.raw_headers[index])
                if not parameter:
                    continue
                parameter_key = header
                role = _chemistry_role(sheet.sheet_name, parameter_key, parameter)
                family = _chemical_family(sheet.sheet_name, parameter)
                is_split_voc = "split into 7 tables" in sheet.workbook_name.lower()
                records.append(
                    {
                        **base,
                        "parameter": parameter,
                        "parameter_key": parameter_key,
                        "parameter_role": role,
                        "chemical_family": family,
                        "unit": _unit_for(sheet.sheet_name, parameter, role),
                        "raw_value": parsed["raw_value"],
                        "value_numeric": parsed["value_numeric"],
                        "value_lower": parsed["value_lower"],
                        "value_upper": parsed["value_upper"],
                        "qualifier": parsed["qualifier"],
                        "measurement_status": parsed["status"],
                        "measurement_flags": "|".join(parsed["flags"]),
                        "reaction_truth_status": "ABSTAIN",
                        "adjacency_truth_status": "ABSTAIN",
                        "preferred_source_for_duplicate_voc": bool(is_split_voc and family == "voc"),
                    }
                )
    frame = pd.DataFrame(records)
    for column in columns:
        if column not in frame:
            frame[column] = pd.NA
    return frame[columns].reset_index(drop=True)


def _parse_model_reference_text(text: str, source_member: str) -> dict[str, Any]:
    record: dict[str, Any] = {
        "source_member": source_member,
        "reference_type": REFERENCE_TYPE,
        "model_reference_status": "parsed",
    }
    for line in text.splitlines():
        clean = line.split("#", 1)[0].strip()
        if not clean:
            continue
        if clean.upper().startswith("ESRI:"):
            record["coordinate_reference_system"] = clean
            continue
        parts = clean.split(None, 1)
        if len(parts) == 2:
            record[parts[0].strip().lower()] = parts[1].strip()
    for key in ("xul", "yul", "rotation"):
        if key in record:
            record[key] = _number(record[key])
    return record


def _parse_model_text(source: "_AikenReader") -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reference_rows: list[dict[str, Any]] = []
    location_rows: list[dict[str, Any]] = []
    name_rows: list[dict[str, Any]] = []
    for member in source.member_names("model"):
        lower = member.lower()
        if lower.endswith("usgs.model.reference"):
            text = _decode_text(source.read_member(member))
            record = _parse_model_reference_text(text, member)
            match = _WELL_RE.search(member)
            record["well_id"] = match.group(0).upper() if match else None
            record["simulation_id"] = member.split("/")[-2] if "/" in member else None
            reference_rows.append(record)
        elif lower.endswith(".loc"):
            match = _WELL_RE.search(member)
            well_id = match.group(0).upper() if match else None
            text = _decode_text(source.read_member(member))
            for line_number, line in enumerate(text.splitlines(), start=1):
                values = line.split()
                if len(values) < 3 or not all(_number(item) is not None for item in values[:3]):
                    continue
                numeric = [_number(item) for item in values]
                location_rows.append(
                    {
                        "node_id": f"AIKEN:{well_id}" if well_id else None,
                        "well_id": well_id,
                        "model_row": int(numeric[0]) if numeric[0] is not None else None,
                        "model_column": int(numeric[1]) if numeric[1] is not None else None,
                        "model_layer": int(numeric[2]) if numeric[2] is not None else None,
                        "local_x": numeric[3] if len(numeric) > 3 else None,
                        "local_y": numeric[4] if len(numeric) > 4 else None,
                        "local_z": numeric[5] if len(numeric) > 5 else None,
                        "source_member": member,
                        "source_row": line_number,
                        "reference_type": REFERENCE_TYPE,
                        "location_reference_kind": "MODPATH5_starting_location",
                        "direct_adjacency_truth_status": "ABSTAIN",
                    }
                )
        elif lower.endswith(".nam"):
            text = source.read_member(member).decode("latin-1", errors="replace")
            simulation_id = member.split("/")[-2] if "/" in member else None
            for line_number, line in enumerate(text.splitlines(), start=1):
                clean = line.strip()
                if not clean or clean.startswith("#"):
                    continue
                parts = clean.split()
                if len(parts) < 2:
                    continue
                name_rows.append(
                    {
                        "simulation_id": simulation_id,
                        "package": parts[0],
                        "unit": _number(parts[1]),
                        "file": " ".join(parts[2:]) if len(parts) > 2 else None,
                        "source_member": member,
                        "source_row": line_number,
                        "reference_type": REFERENCE_TYPE,
                    }
                )
    location_columns = [
        "node_id", "well_id", "model_row", "model_column", "model_layer", "local_x", "local_y", "local_z",
        "source_member", "source_row", "reference_type", "location_reference_kind", "direct_adjacency_truth_status",
    ]
    reference_columns = [
        "source_member", "reference_type", "model_reference_status", "well_id", "simulation_id", "xul", "yul",
        "rotation", "length_units", "time_units", "start_date", "start_time", "model", "coordinate_reference_system",
    ]
    name_columns = ["simulation_id", "package", "unit", "file", "source_member", "source_row", "reference_type"]
    return (
        pd.DataFrame(location_rows, columns=location_columns),
        pd.DataFrame(reference_rows, columns=reference_columns),
        pd.DataFrame(name_rows, columns=name_columns),
    )


def _parse_model_binary(
    source: "_AikenReader",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read the small MODPATH 5 ``.pth``/``.ept`` members in ``output.zip``.

    The Aiken output archive also contains multi-gigabyte MODFLOW heads and
    budgets.  Only the fixed-width particle-coordinate members are selected
    here, and each is read directly from the ZIP stream.  A malformed member
    is recorded in the diagnostics table rather than silently omitted; valid
    members remain available for model-conditioned pathway diagnostics.
    """

    pathline_rows: list[pd.DataFrame] = []
    endpoint_rows: list[pd.DataFrame] = []
    diagnostics: list[dict[str, Any]] = []
    tracking_directions = _model_tracking_directions(source)
    members = [
        member
        for member in source.member_names("output")
        if member.lower().endswith((".pth", ".ept"))
    ]
    for member in members:
        suffix = Path(member.replace("\\", "/")).suffix.lower()
        file_type = "pathline" if suffix == ".pth" else "endpoint"
        metadata: dict[str, Any] = {}
        try:
            payload, metadata = source.read_member_with_metadata(member)
            run_id = _binary_run_id(member)
            direction_code, direction_label = tracking_directions.get(
                run_id,
                (None, "unknown_tracking_direction"),
            )
            if suffix == ".pth":
                frame = parse_modpath5_pathline_binary(
                    payload,
                    run_id=run_id,
                    source_archive=metadata.get("source_archive"),
                    source_member=metadata.get("source_member", member),
                    tracking_direction_code=direction_code,
                    tracking_direction=direction_label,
                )
                pathline_rows.append(frame)
            else:
                frame = parse_modpath5_endpoint_binary(
                    payload,
                    run_id=run_id,
                    source_archive=metadata.get("source_archive"),
                    source_member=metadata.get("source_member", member),
                    tracking_direction_code=direction_code,
                    tracking_direction=direction_label,
                )
                endpoint_rows.append(frame)
            diagnostics.append(
                {
                    "run_id": run_id,
                    "file_type": f"MODPATH5_BINARY_{file_type.upper()}",
                    "source_archive": metadata.get("source_archive"),
                    "source_member": metadata.get("source_member", member),
                    "binary_endian": MODPATH5_BINARY_ENDIAN,
                    "binary_parser": MODPATH5_BINARY_PARSER,
                    "binary_parser_version": MODPATH5_BINARY_PARSER_VERSION,
                    "record_bytes": (
                        MODPATH5_PATHLINE_RECORD_BYTES
                        if file_type == "pathline"
                        else MODPATH5_ENDPOINT_RECORD_BYTES
                    ),
                    "ipcode_encoding": MODPATH5_IPCODE_ENCODING,
                    "tracking_direction_code": direction_code,
                    "tracking_direction": direction_label,
                    "status": "validated",
                    "record_count": int(len(frame)),
                    "error": None,
                    "reference_type": REFERENCE_TYPE,
                }
            )
        except (OSError, KeyError, ValueError, struct.error) as error:
            diagnostics.append(
                {
                    "run_id": _binary_run_id(member),
                    "file_type": f"MODPATH5_BINARY_{file_type.upper()}",
                    "source_archive": metadata.get("source_archive"),
                    "source_member": metadata.get("source_member", member),
                    "binary_endian": MODPATH5_BINARY_ENDIAN,
                    "binary_parser": MODPATH5_BINARY_PARSER,
                    "binary_parser_version": MODPATH5_BINARY_PARSER_VERSION,
                    "record_bytes": (
                        MODPATH5_PATHLINE_RECORD_BYTES
                        if file_type == "pathline"
                        else MODPATH5_ENDPOINT_RECORD_BYTES
                    ),
                    "ipcode_encoding": MODPATH5_IPCODE_ENCODING,
                    "tracking_direction_code": tracking_directions.get(
                        _binary_run_id(member),
                        (None, "unknown_tracking_direction"),
                    )[0],
                    "tracking_direction": tracking_directions.get(
                        _binary_run_id(member),
                        (None, "unknown_tracking_direction"),
                    )[1],
                    "status": "rejected",
                    "record_count": 0,
                    "error": str(error),
                    "reference_type": REFERENCE_TYPE,
                }
            )

    endpoints = (
        pd.concat(endpoint_rows, ignore_index=True)
        if endpoint_rows
        else pd.DataFrame(columns=_MODPATH5_ENDPOINT_COLUMNS)
    )
    pathlines = (
        pd.concat(pathline_rows, ignore_index=True)
        if pathline_rows
        else pd.DataFrame(columns=_MODPATH5_PATHLINE_COLUMNS)
    )
    if not endpoints.empty and not pathlines.empty:
        # Endpoint files have no native particle ID.  The only defensible
        # cross-file key is the documented one-record-per-particle endpoint
        # row order matched to the first-seen particle order in the pathline
        # file.  Do not join on the derived endpoint ``particle_id`` field.
        endpoint_by_ordinal = {
            (row.run_id, int(row.particle_ordinal)): row
            for row in endpoints.itertuples(index=False)
        }
        for index, row in pathlines.iterrows():
            endpoint = endpoint_by_ordinal.get(
                (row["run_id"], int(row["particle_ordinal"]))
            )
            if endpoint is None:
                continue
            pathlines.at[index, "ipcode"] = endpoint.ipcode
            pathlines.at[index, "ipcode_idcode"] = endpoint.ipcode_idcode
            pathlines.at[index, "ipcode_nslast"] = endpoint.ipcode_nslast
            pathlines.at[index, "termination_status"] = endpoint.termination_status
            pathlines.at[index, "particle_censoring_status"] = endpoint.particle_censoring_status
    diagnostic_frame = pd.DataFrame(
        diagnostics,
        columns=[
            "run_id",
            "file_type",
            "source_archive",
            "source_member",
            "binary_endian",
            "binary_parser",
            "binary_parser_version",
            "record_bytes",
            "ipcode_encoding",
            "tracking_direction_code",
            "tracking_direction",
            "status",
            "record_count",
            "error",
            "reference_type",
        ],
    )
    return endpoints, pathlines, diagnostic_frame


class _AikenReader:
    """Resolve members from extracted or ZIP-backed Aiken sources."""

    def __init__(self, source: Path | str):
        self.source = Path(source).expanduser().resolve()
        if not self.source.exists():
            raise FileNotFoundError(self.source)

    def _archives(self) -> list[Path]:
        if self.source.is_file() and self.source.suffix.lower() == ".zip":
            return [self.source]
        if not self.source.is_dir():
            return []
        return sorted(self.source.rglob("*.zip"))

    def _direct_candidates(self, expected: str) -> list[Path]:
        if not self.source.is_dir():
            return []
        expected_path = Path(expected)
        candidates = [self.source / expected_path]
        # An extracted release may have retained SIR2022-5036/ as an extra
        # root; use a basename fallback only when it is unique.
        candidates.extend(
            path for path in self.source.rglob(expected_path.name) if path not in candidates
        )
        return [path for path in candidates if path.is_file()]

    @staticmethod
    def _normalised_member(name: str) -> str:
        return name.replace("\\", "/").lstrip("./").lower()

    def _find_member(self, expected: str, archive_hint: str | None = None) -> tuple[Path, str] | None:
        expected_norm = self._normalised_member(expected)
        for archive in self._archives():
            if archive_hint and archive.name.lower() != archive_hint.lower():
                continue
            with zipfile.ZipFile(archive) as handle:
                names = [info.filename for info in handle.infolist() if not info.is_dir()]
            exact = [name for name in names if self._normalised_member(name) == expected_norm]
            suffix = [name for name in names if self._normalised_member(name).endswith("/" + expected_norm)]
            candidates = exact or suffix
            if not candidates:
                candidates = [name for name in names if Path(name).name.lower() == Path(expected).name.lower()]
            if candidates:
                return archive, sorted(candidates, key=lambda item: (len(item), item))[0]
        return None

    def read_named(self, expected: str, *, archive_hint: str | None = None) -> tuple[bytes, dict[str, Any]] | None:
        direct = self._direct_candidates(expected)
        if direct:
            path = direct[0]
            return path.read_bytes(), {
                "kind": "extracted_file",
                "path": str(path),
                "member": None,
                "sha256": _sha256_file(path),
            }
        found = self._find_member(expected, archive_hint)
        if not found:
            return None
        archive, member = found
        with zipfile.ZipFile(archive) as handle:
            data = handle.read(member)
        return data, {
            "kind": "zip_member",
            "path": str(archive),
            "member": member,
            "sha256": _sha256_bytes(data),
        }

    def read_member(self, member: str) -> bytes:
        if self.source.is_dir():
            direct = self.source / Path(member)
            if direct.is_file():
                return direct.read_bytes()
        wanted = self._normalised_member(member)
        for archive in self._archives():
            with zipfile.ZipFile(archive) as handle:
                names = {self._normalised_member(name): name for name in handle.namelist()}
                actual = names.get(wanted)
                if actual is not None:
                    return handle.read(actual)
        raise KeyError(f"Aiken archive member not found: {member}")

    def read_member_with_metadata(self, member: str) -> tuple[bytes, dict[str, Any]]:
        """Read one member and retain its archive/member provenance.

        Binary transport parsing uses this method so a row can be traced back
        to the exact release member.  It mirrors :meth:`read_member` and does
        not extract a ZIP member to disk.
        """

        if self.source.is_dir():
            direct = self.source / Path(member)
            if direct.is_file():
                return direct.read_bytes(), {
                    "source_archive": None,
                    "source_member": str(Path(member).as_posix()),
                    "kind": "extracted_file",
                }
        wanted = self._normalised_member(member)
        for archive in self._archives():
            with zipfile.ZipFile(archive) as handle:
                names = {self._normalised_member(name): name for name in handle.namelist()}
                actual = names.get(wanted)
                if actual is not None:
                    return handle.read(actual), {
                        "source_archive": str(archive),
                        "source_member": actual,
                        "kind": "zip_member",
                    }
        raise KeyError(f"Aiken archive member not found: {member}")

    def member_names(self, prefix: str | None = None) -> list[str]:
        names: list[str] = []
        for archive in self._archives():
            with zipfile.ZipFile(archive) as handle:
                for info in handle.infolist():
                    if info.is_dir():
                        continue
                    name = info.filename
                    if prefix is None or self._normalised_member(name).startswith(self._normalised_member(prefix)):
                        names.append(name)
        if self.source.is_dir():
            root = self.source / prefix if prefix else self.source
            if root.exists():
                names.extend(str(path.relative_to(self.source)).replace("\\", "/") for path in root.rglob("*") if path.is_file())
        return sorted(set(names))


def inventory_zip(
    path: Path | str,
    *,
    hash_members: bool = True,
    max_member_hash_bytes: int | None = None,
    hash_archive: bool | None = None,
) -> dict[str, Any]:
    """Inventory one ZIP without extraction.

    Member hashes are SHA-256 hashes of the uncompressed member bytes.  A
    size limit can be supplied for a quick audit; skipped hashes are explicit
    in the returned record rather than represented as false values.  Set
    ``hash_archive`` explicitly when the archive-content hash should be
    collected independently of member hashes.
    """

    archive = Path(path).expanduser().resolve()
    if hash_archive is None:
        # Avoid an implicit multi-gigabyte read when a caller explicitly
        # requests metadata/member inventory without content hashes.
        hash_archive = hash_members
    if not archive.is_file():
        raise FileNotFoundError(archive)
    members: list[dict[str, Any]] = []
    with zipfile.ZipFile(archive) as handle:
        for info in handle.infolist():
            member: dict[str, Any] = {
                "name": info.filename,
                "is_dir": info.is_dir(),
                "compressed_size_bytes": int(info.compress_size),
                "size_bytes": int(info.file_size),
                "crc32": f"{info.CRC:08x}",
                "sha256": None,
                "hash_status": "not_requested" if not hash_members else "pending",
            }
            if not info.is_dir() and hash_members:
                if max_member_hash_bytes is not None and info.file_size > max_member_hash_bytes:
                    member["hash_status"] = "skipped_size_limit"
                else:
                    with handle.open(info, "r") as stream:
                        member["sha256"] = _sha256_stream(stream)
                        member["hash_status"] = "complete"
            members.append(member)
    return {
        "path": str(archive),
        "archive_name": archive.name,
        "size_bytes": archive.stat().st_size,
        "sha256": _sha256_file(archive) if hash_archive else None,
        "member_count": len(members),
        "uncompressed_member_bytes": sum(item["size_bytes"] for item in members),
        "members": members,
    }


def inventory_aiken_source(
    source: Path | str,
    *,
    hash_members: bool = True,
    max_member_hash_bytes: int | None = None,
    hash_archives: bool | None = None,
) -> dict[str, Any]:
    """Return a complete, provenance-oriented inventory for an Aiken source.

    ``hash_archives`` defaults to ``hash_members`` so that requesting a
    metadata-only inventory does not silently stream the multi-gigabyte Aiken
    archives.  Both archive and member hashing remain opt-in at the call site
    when a full content fingerprint is required.
    """

    reader = _AikenReader(source)
    if hash_archives is None:
        hash_archives = hash_members
    archives = [
        inventory_zip(
            path,
            hash_members=hash_members,
            max_member_hash_bytes=max_member_hash_bytes,
            hash_archive=hash_archives,
        )
        for path in reader._archives()
    ]
    files: list[dict[str, Any]] = []
    if reader.source.is_dir():
        for name in ("readme.txt", "modelgeoref.txt"):
            matches = reader._direct_candidates(name)
            if matches:
                path = matches[0]
                files.append(
                    {
                        "path": str(path),
                        "name": name,
                        "size_bytes": path.stat().st_size,
                        "sha256": _sha256_file(path),
                    }
                )
    readme = reader.read_named("readme.txt")
    modelgeoref = reader.read_named("modelgeoref.txt")
    return {
        "schema": "aiken-source-inventory-v1",
        "source": str(reader.source),
        "reference_type": REFERENCE_TYPE,
        "source_doi": SOURCE_DOI,
        "publication_doi": PUBLICATION_DOI,
        "hash_policy": {
            "archive_sha256": "complete" if hash_archives else "not_requested",
            "member_sha256": "complete" if hash_members and max_member_hash_bytes is None else ("requested_with_size_limit" if hash_members else "not_requested"),
            "member_hash_bytes_limit": max_member_hash_bytes,
            "member_hash_is_uncompressed_content": True,
        },
        "archives": archives,
        "standalone_files": files,
        "resolved_metadata": {
            "readme": readme[1] if readme else None,
            "modelgeoref": modelgeoref[1] if modelgeoref else None,
        },
    }


@dataclass
class AikenReference:
    """Canonical, panel-separated Aiken reference data."""

    source: Path
    reference_type: str
    readme: dict[str, Any]
    modelgeoref: dict[str, Any]
    inventory: dict[str, Any] | None
    well_metadata: pd.DataFrame
    field_samples: pd.DataFrame
    cfc_ages: pd.DataFrame
    pathways: pd.DataFrame
    chemistry: pd.DataFrame
    model_locations: pd.DataFrame
    model_references: pd.DataFrame
    model_name_files: pd.DataFrame
    native_tables: dict[str, pd.DataFrame] = field(default_factory=dict, repr=False)
    notes: tuple[str, ...] = ()
    # Particle-coordinate outputs are a separate model-conditioned panel.
    # They are intentionally not folded into ``pathways`` or any independent
    # well-to-well truth table.
    endpoints: pd.DataFrame = field(default_factory=pd.DataFrame, repr=False)
    pathlines: pd.DataFrame = field(default_factory=pd.DataFrame, repr=False)
    binary_parse_diagnostics: pd.DataFrame = field(default_factory=pd.DataFrame, repr=False)

    @property
    def modpath_endpoints(self) -> pd.DataFrame:
        """Alias with an explicit MODPATH name for downstream callers."""

        return self.endpoints

    @property
    def modpath_pathlines(self) -> pd.DataFrame:
        """Alias with an explicit MODPATH name for downstream callers."""

        return self.pathlines

    @property
    def audit(self) -> dict[str, Any]:
        """Return the evidence boundary and simple coverage counts."""

        return {
            "reference_type": self.reference_type,
            "status": "COMPLETE" if len(self.well_metadata) else "PARTIAL",
            "n_wells": int(len(self.well_metadata)),
            "n_field_samples": int(len(self.field_samples)),
            "n_cfc_age_records": int(len(self.cfc_ages)),
            "n_pathway_records": int(len(self.pathways)),
            "n_chemistry_records": int(len(self.chemistry)),
            "n_model_locations": int(len(self.model_locations)),
            "n_model_reference_files": int(len(self.model_references)),
            "n_binary_endpoint_records": int(len(self.endpoints)),
            "n_binary_pathline_records": int(len(self.pathlines)),
            "n_binary_parse_diagnostics": int(len(self.binary_parse_diagnostics)),
            "binary_tracking_direction_counts": (
                self.endpoints["tracking_direction"].value_counts(dropna=False).to_dict()
                if "tracking_direction" in self.endpoints
                else {}
            ),
            "binary_censoring_status_counts": (
                self.endpoints["particle_censoring_status"].value_counts(dropna=False).to_dict()
                if "particle_censoring_status" in self.endpoints
                else {}
            ),
            "capabilities": {
                "well_identity": bool(len(self.well_metadata)),
                "sample_dates": bool(len(self.field_samples) and self.field_samples["sample_date"].notna().any()),
                "cfc_apparent_age_intervals": bool(len(self.cfc_ages) and self.cfc_ages["age_parse_status"].isin(["interval", "point"]).any()),
                "model_transport_reference": bool(len(self.pathways) or len(self.model_locations)),
                "binary_modpath_particle_records": bool(len(self.endpoints) or len(self.pathlines)),
                "binary_modpath_outputs_validated": bool(
                    len(self.binary_parse_diagnostics)
                    and self.binary_parse_diagnostics["status"].eq("validated").all()
                ),
                "independent_direct_adjacency_truth": False,
                "independent_age_truth": False,
                "independent_reaction_truth": False,
                "complete_carbonate_major_ion_panel": False,
            },
            "integrated_scoring_allowed": False,
            "missing_or_prohibited": [
                "independent_direct_adjacency_labels",
                "independent_flow_truth",
                "independent_reaction_truth",
                "complete_Ca_Mg_alkalinity_DIC_SiO2_Sr_major_ion_panel",
            ],
            "claim_boundary": (
                "Aiken is a calibrated MODFLOW-NWT/MODPATH5 reference. It supports "
                "model-conditioned pathway, travel-time, age, chemistry, and "
                "scenario-consistency diagnostics. It does not establish direct "
                "well-to-well adjacency or independent field accuracy."
            ),
            "notes": list(self.notes),
        }

    def to_manifest(self) -> dict[str, Any]:
        return {
            "schema": "aiken-reference-manifest-v1",
            "source": str(self.source),
            "reference_type": self.reference_type,
            "source_doi": SOURCE_DOI,
            "publication_doi": PUBLICATION_DOI,
            "readme": self.readme,
            "modelgeoref": self.modelgeoref,
            "audit": self.audit,
            "binary_transport": {
                "endpoint_records": int(len(self.endpoints)),
                "pathline_records": int(len(self.pathlines)),
                "diagnostics": self.binary_parse_diagnostics.to_dict(orient="records"),
            },
            "inventory": self.inventory,
        }


def _select_sheet(sheets: Mapping[str, _ParsedSheet], expected: str) -> _ParsedSheet | None:
    expected_norm = _normalise_text(expected).lower()
    for name, sheet in sheets.items():
        if _normalise_text(name).lower() == expected_norm:
            return sheet
    return next((sheet for name, sheet in sheets.items() if expected_norm in _normalise_text(name).lower()), None)


def load_aiken_reference(
    source: Path | str,
    *,
    include_inventory: bool = False,
    hash_members: bool = False,
    max_member_hash_bytes: int | None = None,
    hash_archives: bool | None = None,
) -> AikenReference:
    """Load Aiken panels without extracting the large model/output archives.

    ``include_inventory`` is opt-in because hashing every uncompressed member
    of the supplied release can require reading many gigabytes.  When it is
    requested, member hashes are streamed and never extracted; set
    ``hash_members=True`` for complete member SHA-256 provenance.  Set
    ``hash_archives=True`` independently when archive-content hashes are also
    required; otherwise it follows ``hash_members``.
    """

    reader = _AikenReader(source)
    readme_payload = reader.read_named("readme.txt")
    georef_payload = reader.read_named("modelgeoref.txt")
    readme_text = _decode_text(readme_payload[0]) if readme_payload else ""
    georef_text = _decode_text(georef_payload[0]) if georef_payload else ""
    readme = parse_readme(readme_text)
    modelgeoref = parse_modelgeoref(georef_text)

    core_payload = reader.read_named(CORE_WORKBOOK_NAME, archive_hint="ancillary.zip")
    voc_payload = reader.read_named(VOC_WORKBOOK_NAME, archive_hint="ancillary.zip")
    parsed_core = _read_workbook_sheets(core_payload[0], CORE_WORKBOOK_NAME) if core_payload else {}
    parsed_voc = _read_workbook_sheets(voc_payload[0], VOC_WORKBOOK_NAME) if voc_payload else {}
    all_sheets: dict[str, _ParsedSheet] = {}
    all_sheets.update({f"core::{name}": sheet for name, sheet in parsed_core.items()})
    all_sheets.update({f"voc::{name}": sheet for name, sheet in parsed_voc.items()})

    table2 = _select_sheet(parsed_core, "Table 2")
    table8 = _select_sheet(parsed_core, "Table 8")
    table13 = _select_sheet(parsed_core, "Table 13")
    table15 = _select_sheet(parsed_core, "Table 15")
    well_metadata = _make_well_metadata(table2)
    field_samples = _make_field_samples(table8)
    cfc_ages = _make_cfc_ages(table13)
    pathways = _make_pathways(table15)

    chemistry_sheets: list[_ParsedSheet] = []
    for name, sheet in parsed_core.items():
        if name in {"Table 8", "Table 9", "Table 10 pt 1", "Table 10 pt 2", "Table 11", "Table 12", "Table 13", "Table 14"}:
            chemistry_sheets.append(sheet)
    chemistry_sheets.extend(parsed_voc.values())
    chemistry = _make_chemistry(chemistry_sheets)

    model_source = _AikenReader(source)
    model_locations, model_references, model_name_files = _parse_model_text(model_source)
    endpoints, pathlines, binary_parse_diagnostics = _parse_model_binary(model_source)

    native_tables: dict[str, pd.DataFrame] = {}
    for key, sheet in all_sheets.items():
        rows = []
        for row in sheet.rows:
            values = row["values"]
            record = {header: _clean_scalar(values[index] if index < len(values) else None) for index, header in enumerate(sheet.headers)}
            record.update(
                {
                    "source_workbook": sheet.workbook_name,
                    "source_sheet": sheet.sheet_name,
                    "source_row": row["row_number"],
                    "reference_type": REFERENCE_TYPE,
                }
            )
            rows.append(record)
        native_tables[key] = pd.DataFrame(rows)

    notes_list = [
        "MODFLOW-NWT binary head and budget outputs were not extracted or decoded.",
        "MODPATH5 binary endpoint and pathline outputs were decoded with strict header and payload-length validation.",
        "Table 15 pathway extents and velocities are model-derived recharge-to-well summaries.",
        "Aiken chemistry and CFC ages do not provide independent direct-adjacency or reaction truth.",
        "Table 2 provides total depth and aquifer labels; screened top/bottom intervals remain unknown.",
    ]
    rejected_binary = int(
        binary_parse_diagnostics["status"].eq("rejected").sum()
        if not binary_parse_diagnostics.empty
        else 0
    )
    if rejected_binary:
        notes_list.append(
            f"{rejected_binary} MODPATH5 binary particle file(s) were rejected; see binary_parse_diagnostics."
        )
    notes = tuple(notes_list)
    inventory = (
        inventory_aiken_source(
            source,
            hash_members=hash_members,
            max_member_hash_bytes=max_member_hash_bytes,
            hash_archives=hash_archives,
        )
        if include_inventory
        else None
    )
    return AikenReference(
        source=reader.source,
        reference_type=REFERENCE_TYPE,
        readme=readme,
        modelgeoref=modelgeoref,
        inventory=inventory,
        well_metadata=well_metadata,
        field_samples=field_samples,
        cfc_ages=cfc_ages,
        pathways=pathways,
        chemistry=chemistry,
        model_locations=model_locations,
        model_references=model_references,
        model_name_files=model_name_files,
        native_tables=native_tables,
        notes=notes,
        endpoints=endpoints,
        pathlines=pathlines,
        binary_parse_diagnostics=binary_parse_diagnostics,
    )


def write_aiken_manifest(reference: AikenReference, path: Path | str) -> None:
    """Write only the manifest; canonical tables remain caller-controlled."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(reference.to_manifest(), indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


__all__ = [
    "AikenReference",
    "CORE_WORKBOOK_NAME",
    "MODPATH5_ENDPOINT_RECORD_BYTES",
    "MODPATH5_BINARY_ENDIAN",
    "MODPATH5_BINARY_PARSER",
    "MODPATH5_BINARY_PARSER_VERSION",
    "MODPATH5_HEADER_BYTES",
    "MODPATH5_IPCODE_ENCODING",
    "MODPATH5_PATHLINE_RECORD_BYTES",
    "MODPATH5_TRACKING_DIRECTIONS",
    "PUBLICATION_DOI",
    "REFERENCE_TYPE",
    "SOURCE_DOI",
    "SOURCE_URL",
    "VOC_WORKBOOK_NAME",
    "decode_modpath5_ipcode",
    "inventory_aiken_source",
    "inventory_zip",
    "load_aiken_reference",
    "parse_apparent_age_label",
    "parse_modelgeoref",
    "parse_modpath5_endpoint_binary",
    "parse_modpath5_pathline_binary",
    "parse_qualified_value",
    "parse_readme",
    "parse_usgs_site_coordinates",
    "write_aiken_manifest",
]
