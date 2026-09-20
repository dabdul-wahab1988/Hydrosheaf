"""Opt-in adapter for the external WATRES temporal-TTD benchmark.

WATRES is the catchment-scale benchmark accompanying Duchemin et al. (2026),
``Data-Driven Estimation of Time-Variable Catchment Transit Time
Distributions``.  It is useful as an *external temporal-TTD component*
benchmark.  It is not a groundwater-network truth set and must not be used to
validate HydroSheaf graph topology, edge-wise groundwater TTDs, or field
transfer.

The public Zenodo catalogue lists one ``data.zip`` artifact (displayed as
1.9 GB) for record 15658651.  The archive is intentionally neither packaged
nor fetched at import time.  Callers must explicitly opt into both metadata
fetching and archive downloading, and downloaded archives are checked against
the catalogue MD5 before they are made visible at their requested destination.

Because the public catalogue exposes the opaque ZIP artifact but not a stable
machine-readable internal file schema, this module deliberately performs no
automatic column or member guessing.  ``inspect_watres_archive`` reports a
ZIP layout without extraction.  ``load_watres_csv`` only reads a member after
the caller supplies an explicit, auditable column mapping.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.request import Request, urlopen
import zipfile

WATRES_BENCHMARK_SCOPE = (
    "External catchment-level temporal-TTD component benchmark only; it is "
    "not groundwater-graph validation."
)


class WATRESIntegrationError(RuntimeError):
    """Base exception for the deliberately narrow WATRES integration."""


class WATRESArchiveUnavailable(WATRESIntegrationError):
    """Raised when an explicitly requested local archive is absent."""


class WATRESDownloadRefused(WATRESIntegrationError):
    """Raised when a network transfer was not explicitly authorized."""


class WATRESIntegrityError(WATRESIntegrationError):
    """Raised when an archive violates the pinned checksum or size policy."""


class WATRESMetadataError(WATRESIntegrationError):
    """Raised when fetched metadata is not a JSON object."""


class WATRESSchemaUnsupported(WATRESIntegrationError):
    """Raised rather than guessing an unverified WATRES archive schema."""


@dataclass(frozen=True, slots=True)
class WATRESManifest:
    """Pinned identity and transfer policy for one external WATRES release.

    ``archive_max_bytes`` is a *ceiling*, not a claimed byte-exact catalogue
    size.  Zenodo's public record displayed a rounded ``1.9 GB`` size when
    this manifest was authored, so the exact MD5 is the immutable identity
    check and the byte ceiling protects against an accidental multi-gigabyte
    transfer or unexpectedly changed artifact.
    """

    manifest_version: str
    record_id: int
    doi: str
    release_version: str
    archive_name: str
    archive_md5: str
    archive_display_size: str
    archive_max_bytes: int
    record_url: str
    metadata_api_url: str
    archive_url: str
    benchmark_scope: str
    excluded_claims: tuple[str, ...]


# Source facts, checked against the public Zenodo catalogue on 2026-09-17:
# DOI 10.5281/zenodo.15658651; version v1; one artifact, data.zip; MD5
# 0154de3424de65f112ee9a409f2bf26f; displayed size 1.9 GB.  The archive is
# intentionally not represented by an exact byte-size assertion because the
# catalogue display is rounded and no record API response was bundled here.
WATRES_MANIFEST_V1 = WATRESManifest(
    manifest_version="hydrosheaf.watres-manifest/v1",
    record_id=15658651,
    doi="10.5281/zenodo.15658651",
    release_version="v1",
    archive_name="data.zip",
    archive_md5="0154de3424de65f112ee9a409f2bf26f",
    archive_display_size="1.9 GB (Zenodo catalogue display; rounded)",
    # 2.2 GB accommodates a decimal/GiB display ambiguity while still being a
    # meaningful guard for the pinned 1.9 GB artifact.
    archive_max_bytes=2_200_000_000,
    record_url="https://zenodo.org/records/15658651",
    metadata_api_url="https://zenodo.org/api/records/15658651",
    archive_url=("https://zenodo.org/records/15658651/files/data.zip?download=1"),
    benchmark_scope=WATRES_BENCHMARK_SCOPE,
    excluded_claims=(
        "groundwater graph-topology validation",
        "edge-wise groundwater TTD recovery validation",
        "groundwater field-transfer validation",
        "independent groundwater-age truth",
    ),
)


@dataclass(frozen=True, slots=True)
class WATRESMetadataValidation:
    """Result of comparing a supplied Zenodo JSON object to a pinned manifest."""

    manifest_version: str
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    observed_record_id: str | None
    observed_doi: str | None
    observed_release_version: str | None
    observed_archive_size_bytes: int | None

    @property
    def valid(self) -> bool:
        """Whether identity-critical record, DOI, archive, and checksum match."""

        return not self.errors


@dataclass(frozen=True, slots=True)
class WATRESArchiveIntegrity:
    """Checksum and size evidence for a local WATRES archive."""

    archive_path: Path
    size_bytes: int
    md5: str
    manifest_version: str


@dataclass(frozen=True, slots=True)
class WATRESDownloadReceipt:
    """Evidence emitted after a caller-authorized, verified download."""

    archive_path: Path
    bytes_downloaded: int
    md5: str
    manifest_version: str
    source_url: str


@dataclass(frozen=True, slots=True)
class WATRESArchiveEntry:
    """A ZIP central-directory entry; no member content is extracted."""

    name: str
    compressed_size_bytes: int
    uncompressed_size_bytes: int


@dataclass(frozen=True, slots=True)
class WATRESArchiveInspection:
    """Safe central-directory inspection of a local external archive."""

    archive_path: Path
    archive_size_bytes: int
    entries: tuple[WATRESArchiveEntry, ...]
    checksum_verified: bool
    automatic_schema_supported: bool
    schema_note: str
    benchmark_scope: str


@dataclass(frozen=True, slots=True)
class WATRESCSVAdapterSpec:
    """Explicit CSV member and column mapping for a caller-audited archive layout.

    The names are intentionally not inferred.  A research record should retain
    this mapping alongside its benchmark manifest so another run can reproduce
    exactly which WATRES member and columns were used.
    """

    member: str
    timestamp_column: str
    input_tracer_column: str
    output_tracer_column: str
    precipitation_column: str | None = None
    streamflow_column: str | None = None
    potential_evapotranspiration_column: str | None = None

    def __post_init__(self) -> None:
        required_text = {
            "member": self.member,
            "timestamp_column": self.timestamp_column,
            "input_tracer_column": self.input_tracer_column,
            "output_tracer_column": self.output_tracer_column,
        }
        blank = [name for name, value in required_text.items() if not value.strip()]
        if blank:
            raise ValueError(f"WATRES adapter fields must be non-empty: {blank}.")
        hydrologic_columns = (
            self.precipitation_column,
            self.streamflow_column,
            self.potential_evapotranspiration_column,
        )
        if any(value is None for value in hydrologic_columns) and any(
            value is not None for value in hydrologic_columns
        ):
            raise ValueError(
                "Map precipitation, streamflow, and potential evapotranspiration "
                "together, or omit all three for a temporal-TTD component-only read."
            )
        if any(value is not None and not value.strip() for value in hydrologic_columns):
            raise ValueError(
                "Optional WATRES hydrologic column names must be non-empty."
            )


@dataclass(frozen=True, slots=True)
class WATRESTemporalSeries:
    """Explicitly mapped input/output series from an external catchment archive."""

    timestamps: tuple[datetime, ...]
    input_values: tuple[float, ...]
    output_values: tuple[float, ...]
    archive_member: str
    timestamp_column: str
    input_tracer_column: str
    output_tracer_column: str
    precipitation_values: tuple[float, ...] | None = None
    streamflow_values: tuple[float, ...] | None = None
    potential_evapotranspiration_values: tuple[float, ...] | None = None
    benchmark_scope: str = WATRES_BENCHMARK_SCOPE

    def __post_init__(self) -> None:
        n = len(self.timestamps)
        if n == 0:
            raise ValueError("WATRES temporal series must contain at least one row.")
        if len(self.input_values) != n or len(self.output_values) != n:
            raise ValueError("WATRES temporal series columns must have equal lengths.")
        if not all(math.isfinite(value) for value in self.input_values):
            raise ValueError("WATRES input tracer values must be finite.")
        if not all(math.isfinite(value) for value in self.output_values):
            raise ValueError("WATRES output tracer values must be finite.")
        hydrologic_values = (
            self.precipitation_values,
            self.streamflow_values,
            self.potential_evapotranspiration_values,
        )
        if any(values is None for values in hydrologic_values) and any(
            values is not None for values in hydrologic_values
        ):
            raise ValueError(
                "WATRES hydrologic values must be supplied together or omitted together."
            )
        for name, values in (
            ("precipitation", self.precipitation_values),
            ("streamflow", self.streamflow_values),
            ("potential evapotranspiration", self.potential_evapotranspiration_values),
        ):
            if values is not None and len(values) != n:
                raise ValueError(f"WATRES {name} values must match timestamp length.")
            if values is not None and not all(math.isfinite(value) for value in values):
                raise ValueError(f"WATRES {name} values must be finite.")


def _as_nonempty_string(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _normalise_md5(value: Any) -> str | None:
    text = _as_nonempty_string(value)
    if text is None:
        return None
    if text.lower().startswith("md5:"):
        text = text.split(":", 1)[1]
    text = text.lower()
    if len(text) != 32 or any(
        character not in "0123456789abcdef" for character in text
    ):
        return None
    return text


def _metadata_doi(
    payload: Mapping[str, Any], metadata: Mapping[str, Any]
) -> str | None:
    direct = _as_nonempty_string(metadata.get("doi"))
    if direct:
        return direct
    pids = payload.get("pids")
    if isinstance(pids, Mapping):
        doi_entry = pids.get("doi")
        if isinstance(doi_entry, Mapping):
            return _as_nonempty_string(doi_entry.get("identifier"))
        return _as_nonempty_string(doi_entry)
    return None


def _metadata_files(payload: Mapping[str, Any]) -> tuple[Mapping[str, Any], ...]:
    """Accommodate common Zenodo API file list shapes without guessing content."""

    files_value = payload.get("files", ())
    if isinstance(files_value, Mapping):
        files_value = files_value.get("entries", ())
    if isinstance(files_value, Mapping):
        files_value = tuple(files_value.values())
    if not isinstance(files_value, Sequence) or isinstance(files_value, (str, bytes)):
        return ()
    return tuple(item for item in files_value if isinstance(item, Mapping))


def validate_watres_metadata(
    payload: Mapping[str, Any],
    *,
    manifest: WATRESManifest = WATRES_MANIFEST_V1,
) -> WATRESMetadataValidation:
    """Validate provided Zenodo metadata against the immutable WATRES manifest.

    This is a pure function.  It does not contact Zenodo, mutate the payload,
    or silently update the pinned manifest.  Absence of a byte-exact file size
    is a warning because the public catalogue exposed a rounded display size;
    a size exceeding the transfer ceiling is an error.
    """

    if not isinstance(payload, Mapping):
        raise WATRESMetadataError("Zenodo metadata must be a mapping/object.")

    errors: list[str] = []
    warnings: list[str] = []
    metadata_value = payload.get("metadata", {})
    metadata = metadata_value if isinstance(metadata_value, Mapping) else {}
    if not isinstance(metadata_value, Mapping):
        errors.append("Metadata payload has no object-valued 'metadata' field.")

    observed_record_id = _as_nonempty_string(payload.get("id"))
    if observed_record_id != str(manifest.record_id):
        errors.append(
            f"Record id mismatch: expected {manifest.record_id}, got "
            f"{observed_record_id!r}."
        )

    observed_doi = _metadata_doi(payload, metadata)
    if (observed_doi or "").lower() != manifest.doi.lower():
        errors.append(f"DOI mismatch: expected {manifest.doi}, got {observed_doi!r}.")

    observed_release_version = _as_nonempty_string(metadata.get("version"))
    if observed_release_version != manifest.release_version:
        errors.append(
            "Release version mismatch: expected "
            f"{manifest.release_version!r}, got {observed_release_version!r}."
        )

    resource_type = metadata.get("resource_type")
    if isinstance(resource_type, Mapping):
        resource_type = resource_type.get("id") or resource_type.get("title")
    resource_type_text = _as_nonempty_string(resource_type)
    if resource_type_text and resource_type_text.lower() != "dataset":
        warnings.append(
            "The record identifies a resource type other than dataset: "
            f"{resource_type_text!r}."
        )

    archive_entries = [
        entry
        for entry in _metadata_files(payload)
        if _as_nonempty_string(entry.get("key") or entry.get("name"))
        == manifest.archive_name
    ]
    if len(archive_entries) != 1:
        errors.append(
            "Expected exactly one pinned archive "
            f"{manifest.archive_name!r}; found {len(archive_entries)}."
        )
        observed_archive_size: int | None = None
    else:
        archive = archive_entries[0]
        observed_md5 = _normalise_md5(archive.get("checksum") or archive.get("md5"))
        if observed_md5 != manifest.archive_md5:
            errors.append(
                "Archive MD5 mismatch: expected "
                f"{manifest.archive_md5}, got {observed_md5!r}."
            )
        raw_size = archive.get("size")
        try:
            observed_archive_size = int(raw_size) if raw_size is not None else None
        except (TypeError, ValueError):
            observed_archive_size = None
            warnings.append(f"Archive size is not an integer: {raw_size!r}.")
        if observed_archive_size is None:
            warnings.append(
                "Zenodo metadata omitted an exact archive byte size; only the "
                "pinned MD5 and local transfer ceiling can be checked."
            )
        elif observed_archive_size <= 0:
            errors.append(
                f"Archive size must be positive, got {observed_archive_size}."
            )
        elif observed_archive_size > manifest.archive_max_bytes:
            errors.append(
                "Archive size exceeds the pinned transfer ceiling: "
                f"{observed_archive_size} > {manifest.archive_max_bytes}."
            )

    return WATRESMetadataValidation(
        manifest_version=manifest.manifest_version,
        errors=tuple(errors),
        warnings=tuple(warnings),
        observed_record_id=observed_record_id,
        observed_doi=observed_doi,
        observed_release_version=observed_release_version,
        observed_archive_size_bytes=observed_archive_size,
    )


def fetch_watres_metadata(
    *,
    allow_network: bool = False,
    manifest: WATRESManifest = WATRES_MANIFEST_V1,
    timeout_seconds: float = 30.0,
    maximum_response_bytes: int = 4_000_000,
) -> Mapping[str, Any]:
    """Fetch the small Zenodo JSON metadata document only after explicit opt-in.

    It never downloads ``data.zip``.  Call :func:`validate_watres_metadata`
    separately so callers can retain both the raw metadata and validation
    evidence in their own provenance record.
    """

    if not allow_network:
        raise WATRESDownloadRefused(
            "WATRES metadata access is opt-in. Pass allow_network=True to "
            "fetch the Zenodo JSON document explicitly."
        )
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive.")
    if maximum_response_bytes <= 0:
        raise ValueError("maximum_response_bytes must be positive.")

    request = Request(
        manifest.metadata_api_url,
        headers={"User-Agent": "HydroSheaf-WATRES-adapter/1"},
    )
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            body = response.read(maximum_response_bytes + 1)
    except OSError as exc:
        raise WATRESMetadataError(
            f"Unable to fetch WATRES metadata from {manifest.metadata_api_url}: {exc}"
        ) from exc
    if len(body) > maximum_response_bytes:
        raise WATRESMetadataError(
            "WATRES metadata response exceeded the configured safety limit of "
            f"{maximum_response_bytes} bytes."
        )
    try:
        parsed = json.loads(body.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise WATRESMetadataError("Zenodo response was not valid UTF-8 JSON.") from exc
    if not isinstance(parsed, Mapping):
        raise WATRESMetadataError("Zenodo metadata response was not a JSON object.")
    return parsed


def _md5_file(path: Path, *, chunk_size: int = 1_048_576) -> str:
    digest = hashlib.md5()  # nosec B324 - catalogued artifact integrity, not security
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _assert_archive_size(path: Path, manifest: WATRESManifest) -> int:
    if not path.is_file():
        raise WATRESArchiveUnavailable(
            f"WATRES archive is unavailable at {path}. It is not bundled with "
            "HydroSheaf; download it explicitly or supply a verified local copy."
        )
    size_bytes = path.stat().st_size
    if size_bytes <= 0:
        raise WATRESIntegrityError(f"WATRES archive {path} is empty.")
    if size_bytes > manifest.archive_max_bytes:
        raise WATRESIntegrityError(
            "WATRES archive exceeds the pinned transfer ceiling: "
            f"{size_bytes} > {manifest.archive_max_bytes} bytes."
        )
    return size_bytes


def verify_watres_archive(
    archive_path: str | Path,
    *,
    manifest: WATRESManifest = WATRES_MANIFEST_V1,
) -> WATRESArchiveIntegrity:
    """Check a local archive's size and exact pinned MD5 without extracting it."""

    path = Path(archive_path)
    size_bytes = _assert_archive_size(path, manifest)
    digest = _md5_file(path)
    if digest != manifest.archive_md5:
        raise WATRESIntegrityError(
            "WATRES archive MD5 mismatch: expected "
            f"{manifest.archive_md5}, got {digest}."
        )
    return WATRESArchiveIntegrity(
        archive_path=path,
        size_bytes=size_bytes,
        md5=digest,
        manifest_version=manifest.manifest_version,
    )


def download_watres_archive(
    destination: str | Path,
    *,
    allow_download: bool = False,
    manifest: WATRESManifest = WATRES_MANIFEST_V1,
    timeout_seconds: float = 60.0,
    chunk_size: int = 1_048_576,
    overwrite: bool = False,
) -> WATRESDownloadReceipt:
    """Download WATRES only after an explicit, caller-authorized request.

    The transfer is streamed into a sibling ``.part`` file.  The archive is
    promoted atomically only after its byte ceiling and pinned MD5 both pass.
    The function never runs during import or ordinary benchmark execution.
    """

    if not allow_download:
        raise WATRESDownloadRefused(
            "WATRES data.zip is approximately 1.9 GB and is never downloaded "
            "implicitly. Pass allow_download=True to authorize this transfer."
        )
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive.")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive.")

    path = Path(destination)
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite existing destination {path}; pass "
            "overwrite=True after checking the target."
        )
    part_path = path.with_name(path.name + ".part")
    if part_path.exists():
        raise FileExistsError(
            f"Refusing to replace existing partial download {part_path}. "
            "Inspect or remove it explicitly before retrying."
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    request = Request(
        manifest.archive_url,
        headers={"User-Agent": "HydroSheaf-WATRES-adapter/1"},
    )
    created_partial = False
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            content_length = response.headers.get("Content-Length")
            if content_length is not None:
                try:
                    declared_bytes = int(content_length)
                except ValueError as exc:
                    raise WATRESIntegrityError(
                        f"Invalid Content-Length from WATRES source: {content_length!r}."
                    ) from exc
                if declared_bytes <= 0 or declared_bytes > manifest.archive_max_bytes:
                    raise WATRESIntegrityError(
                        "Refusing WATRES download with Content-Length "
                        f"{declared_bytes}; allowed maximum is "
                        f"{manifest.archive_max_bytes}."
                    )

            digest = hashlib.md5()  # nosec B324 - catalogue identity verification
            total_bytes = 0
            with part_path.open("xb") as handle:
                created_partial = True
                while chunk := response.read(chunk_size):
                    total_bytes += len(chunk)
                    if total_bytes > manifest.archive_max_bytes:
                        raise WATRESIntegrityError(
                            "WATRES download exceeded the pinned transfer ceiling of "
                            f"{manifest.archive_max_bytes} bytes."
                        )
                    digest.update(chunk)
                    handle.write(chunk)
    except Exception:
        if created_partial:
            part_path.unlink(missing_ok=True)
        raise

    actual_md5 = digest.hexdigest()
    if total_bytes <= 0:
        part_path.unlink(missing_ok=True)
        raise WATRESIntegrityError("WATRES download completed with zero bytes.")
    if actual_md5 != manifest.archive_md5:
        part_path.unlink(missing_ok=True)
        raise WATRESIntegrityError(
            "Downloaded WATRES archive MD5 mismatch: expected "
            f"{manifest.archive_md5}, got {actual_md5}."
        )

    os.replace(part_path, path)
    return WATRESDownloadReceipt(
        archive_path=path,
        bytes_downloaded=total_bytes,
        md5=actual_md5,
        manifest_version=manifest.manifest_version,
        source_url=manifest.archive_url,
    )


def _check_safe_member_name(name: str) -> None:
    normalised = name.replace("\\", "/")
    pure_path = PurePosixPath(normalised)
    if pure_path.is_absolute() or ".." in pure_path.parts:
        raise WATRESSchemaUnsupported(
            f"Archive member has an unsafe path and will not be handled: {name!r}."
        )


def inspect_watres_archive(
    archive_path: str | Path,
    *,
    manifest: WATRESManifest = WATRES_MANIFEST_V1,
    verify_checksum: bool = False,
    max_entries: int = 100_000,
    max_total_uncompressed_bytes: int = 20_000_000_000,
    max_member_uncompressed_bytes: int = 20_000_000_000,
    max_compression_ratio: float = 1_000.0,
) -> WATRESArchiveInspection:
    """Inspect ZIP metadata without extraction or automatic schema inference.

    The public WATRES record does not publish an internal CSV/NetCDF schema,
    so this result always marks automatic schema support as false.  It is a
    discovery aid only; use an explicit :class:`WATRESCSVAdapterSpec` after
    auditing the resulting member list.
    """

    if max_entries <= 0:
        raise ValueError("max_entries must be positive.")
    if max_total_uncompressed_bytes <= 0 or max_member_uncompressed_bytes <= 0:
        raise ValueError("Uncompressed-size limits must be positive.")
    if max_compression_ratio <= 0:
        raise ValueError("max_compression_ratio must be positive.")

    path = Path(archive_path)
    archive_size = _assert_archive_size(path, manifest)
    if verify_checksum:
        verify_watres_archive(path, manifest=manifest)

    try:
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
    except zipfile.BadZipFile as exc:
        raise WATRESSchemaUnsupported(
            f"WATRES archive {path} is not a readable ZIP file."
        ) from exc
    if len(infos) > max_entries:
        raise WATRESSchemaUnsupported(
            f"Archive has {len(infos)} members, exceeding max_entries={max_entries}."
        )

    total_uncompressed = 0
    entries: list[WATRESArchiveEntry] = []
    for info in infos:
        _check_safe_member_name(info.filename)
        if info.flag_bits & 0x1:
            raise WATRESSchemaUnsupported(
                f"Encrypted archive member is unsupported: {info.filename!r}."
            )
        if info.file_size > max_member_uncompressed_bytes:
            raise WATRESSchemaUnsupported(
                "Archive member exceeds the configured uncompressed-size limit: "
                f"{info.filename!r}."
            )
        total_uncompressed += info.file_size
        if total_uncompressed > max_total_uncompressed_bytes:
            raise WATRESSchemaUnsupported(
                "Archive total uncompressed size exceeds the configured limit of "
                f"{max_total_uncompressed_bytes} bytes."
            )
        if info.file_size and info.compress_size:
            ratio = info.file_size / info.compress_size
            if ratio > max_compression_ratio:
                raise WATRESSchemaUnsupported(
                    "Archive member exceeds the configured compression-ratio limit: "
                    f"{info.filename!r} ({ratio:.1f})."
                )
        entries.append(
            WATRESArchiveEntry(
                name=info.filename,
                compressed_size_bytes=info.compress_size,
                uncompressed_size_bytes=info.file_size,
            )
        )

    return WATRESArchiveInspection(
        archive_path=path,
        archive_size_bytes=archive_size,
        entries=tuple(entries),
        checksum_verified=verify_checksum,
        automatic_schema_supported=False,
        schema_note=(
            "No WATRES internal archive schema is pinned by the public Zenodo "
            "catalogue. Automatic member/column inference is deliberately disabled; "
            "audit the layout and provide an explicit adapter specification."
        ),
        benchmark_scope=manifest.benchmark_scope,
    )


def _parse_iso8601(value: str, *, row_number: int, column: str) -> datetime:
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        return datetime.fromisoformat(text)
    except ValueError as exc:
        raise WATRESSchemaUnsupported(
            f"Row {row_number} column {column!r} is not ISO-8601: {value!r}."
        ) from exc


def _parse_finite_float(value: str, *, row_number: int, column: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise WATRESSchemaUnsupported(
            f"Row {row_number} column {column!r} is not numeric: {value!r}."
        ) from exc
    if not math.isfinite(parsed):
        raise WATRESSchemaUnsupported(
            f"Row {row_number} column {column!r} must be finite, got {value!r}."
        )
    return parsed


def load_watres_csv(
    archive_path: str | Path,
    spec: WATRESCSVAdapterSpec,
    *,
    manifest: WATRESManifest = WATRES_MANIFEST_V1,
    verify_checksum: bool = True,
    max_rows: int = 5_000_000,
    max_member_uncompressed_bytes: int = 20_000_000_000,
    max_compression_ratio: float = 1_000.0,
) -> WATRESTemporalSeries:
    """Read one explicitly mapped CSV member without extracting the archive.

    This is an adapter, not automatic WATRES schema support.  It rejects absent
    members, duplicate member names, missing mapped columns, malformed values,
    mixed naive/aware timestamps, and non-increasing timestamps rather than
    silently interpolating, filling, sorting, or guessing. The optional
    precipitation/streamflow/PET mapping is all-or-nothing because those are
    the hydrologic inputs described for WATRES; leaving all three unset is
    permitted only for HydroSheaf's bounded temporal-TTD component read.
    """

    if max_rows <= 0 or max_member_uncompressed_bytes <= 0:
        raise ValueError("CSV row and member-size limits must be positive.")
    if max_compression_ratio <= 0:
        raise ValueError("max_compression_ratio must be positive.")
    _check_safe_member_name(spec.member)
    path = Path(archive_path)
    _assert_archive_size(path, manifest)
    if verify_checksum:
        verify_watres_archive(path, manifest=manifest)

    try:
        with zipfile.ZipFile(path) as archive:
            matches = [
                info for info in archive.infolist() if info.filename == spec.member
            ]
            if len(matches) != 1:
                raise WATRESSchemaUnsupported(
                    "Explicit WATRES CSV member must occur exactly once; "
                    f"{spec.member!r} occurred {len(matches)} times."
                )
            info = matches[0]
            if info.flag_bits & 0x1:
                raise WATRESSchemaUnsupported(
                    f"Encrypted CSV member is unsupported: {spec.member!r}."
                )
            if info.file_size > max_member_uncompressed_bytes:
                raise WATRESSchemaUnsupported(
                    "CSV member exceeds the configured uncompressed-size limit: "
                    f"{info.file_size} > {max_member_uncompressed_bytes}."
                )
            if info.file_size and info.compress_size:
                ratio = info.file_size / info.compress_size
                if ratio > max_compression_ratio:
                    raise WATRESSchemaUnsupported(
                        "CSV member exceeds the configured compression-ratio limit: "
                        f"{ratio:.1f} > {max_compression_ratio}."
                    )

            with archive.open(info, "r") as raw:
                with io.TextIOWrapper(raw, encoding="utf-8-sig", newline="") as text:
                    reader = csv.DictReader(text)
                    headers = set(reader.fieldnames or ())
                    required = {
                        spec.timestamp_column,
                        spec.input_tracer_column,
                        spec.output_tracer_column,
                    }
                    hydrologic_columns = (
                        spec.precipitation_column,
                        spec.streamflow_column,
                        spec.potential_evapotranspiration_column,
                    )
                    if all(column is not None for column in hydrologic_columns):
                        required.update(
                            column
                            for column in hydrologic_columns
                            if column is not None
                        )
                    missing = sorted(required - headers)
                    if missing:
                        raise WATRESSchemaUnsupported(
                            "CSV member does not satisfy the explicit adapter mapping; "
                            f"missing columns: {missing}."
                        )

                    timestamps: list[datetime] = []
                    input_values: list[float] = []
                    output_values: list[float] = []
                    precipitation_values: list[float] | None = (
                        [] if spec.precipitation_column is not None else None
                    )
                    streamflow_values: list[float] | None = (
                        [] if spec.streamflow_column is not None else None
                    )
                    potential_evapotranspiration_values: list[float] | None = (
                        []
                        if spec.potential_evapotranspiration_column is not None
                        else None
                    )
                    aware: bool | None = None
                    previous: datetime | None = None
                    for row_number, row in enumerate(reader, start=2):
                        if row_number - 1 > max_rows:
                            raise WATRESSchemaUnsupported(
                                f"CSV member exceeds max_rows={max_rows}."
                            )
                        timestamp_text = row.get(spec.timestamp_column)
                        input_text = row.get(spec.input_tracer_column)
                        output_text = row.get(spec.output_tracer_column)
                        if (
                            timestamp_text is None
                            or input_text is None
                            or output_text is None
                        ):
                            raise WATRESSchemaUnsupported(
                                f"Row {row_number} is missing an explicitly mapped value."
                            )
                        timestamp = _parse_iso8601(
                            timestamp_text,
                            row_number=row_number,
                            column=spec.timestamp_column,
                        )
                        timestamp_is_aware = (
                            timestamp.tzinfo is not None
                            and timestamp.utcoffset() is not None
                        )
                        if aware is None:
                            aware = timestamp_is_aware
                        elif aware != timestamp_is_aware:
                            raise WATRESSchemaUnsupported(
                                "CSV timestamp column mixes timezone-aware and naive "
                                f"values at row {row_number}."
                            )
                        if previous is not None and timestamp <= previous:
                            raise WATRESSchemaUnsupported(
                                "CSV timestamp column must be strictly increasing; "
                                f"row {row_number} is not later than the preceding row."
                            )
                        timestamps.append(timestamp)
                        input_values.append(
                            _parse_finite_float(
                                input_text,
                                row_number=row_number,
                                column=spec.input_tracer_column,
                            )
                        )
                        output_values.append(
                            _parse_finite_float(
                                output_text,
                                row_number=row_number,
                                column=spec.output_tracer_column,
                            )
                        )
                        if precipitation_values is not None:
                            assert spec.precipitation_column is not None
                            precipitation_text = row.get(spec.precipitation_column)
                            if precipitation_text is None:
                                raise WATRESSchemaUnsupported(
                                    f"Row {row_number} is missing precipitation input."
                                )
                            precipitation_values.append(
                                _parse_finite_float(
                                    precipitation_text,
                                    row_number=row_number,
                                    column=spec.precipitation_column,
                                )
                            )
                        if streamflow_values is not None:
                            assert spec.streamflow_column is not None
                            streamflow_text = row.get(spec.streamflow_column)
                            if streamflow_text is None:
                                raise WATRESSchemaUnsupported(
                                    f"Row {row_number} is missing streamflow input."
                                )
                            streamflow_values.append(
                                _parse_finite_float(
                                    streamflow_text,
                                    row_number=row_number,
                                    column=spec.streamflow_column,
                                )
                            )
                        if potential_evapotranspiration_values is not None:
                            assert spec.potential_evapotranspiration_column is not None
                            potential_evapotranspiration_text = row.get(
                                spec.potential_evapotranspiration_column
                            )
                            if potential_evapotranspiration_text is None:
                                raise WATRESSchemaUnsupported(
                                    f"Row {row_number} is missing potential evapotranspiration input."
                                )
                            potential_evapotranspiration_values.append(
                                _parse_finite_float(
                                    potential_evapotranspiration_text,
                                    row_number=row_number,
                                    column=(spec.potential_evapotranspiration_column),
                                )
                            )
                        previous = timestamp
    except zipfile.BadZipFile as exc:
        raise WATRESSchemaUnsupported(
            f"WATRES archive {path} is not a readable ZIP file."
        ) from exc

    return WATRESTemporalSeries(
        timestamps=tuple(timestamps),
        input_values=tuple(input_values),
        output_values=tuple(output_values),
        archive_member=spec.member,
        timestamp_column=spec.timestamp_column,
        input_tracer_column=spec.input_tracer_column,
        output_tracer_column=spec.output_tracer_column,
        precipitation_values=(
            tuple(precipitation_values) if precipitation_values is not None else None
        ),
        streamflow_values=(
            tuple(streamflow_values) if streamflow_values is not None else None
        ),
        potential_evapotranspiration_values=(
            tuple(potential_evapotranspiration_values)
            if potential_evapotranspiration_values is not None
            else None
        ),
        benchmark_scope=manifest.benchmark_scope,
    )
