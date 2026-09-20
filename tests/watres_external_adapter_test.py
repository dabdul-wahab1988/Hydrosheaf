"""Fixture-only tests for the opt-in external WATRES adapter.

These tests deliberately never contact Zenodo and never request the 1.9 GB
external archive.  They use small ZIP fixtures with a test-specific manifest
so integrity checks and explicit schema mappings are still exercised.
"""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
import hashlib
from pathlib import Path
import zipfile

import pytest

from hydrosheaf.benchmarks.watres import (
    WATRES_BENCHMARK_SCOPE,
    WATRES_MANIFEST_V1,
    WATRESArchiveUnavailable,
    WATRESCSVAdapterSpec,
    WATRESDownloadRefused,
    WATRESSchemaUnsupported,
    download_watres_archive,
    fetch_watres_metadata,
    inspect_watres_archive,
    load_watres_csv,
    validate_watres_metadata,
    verify_watres_archive,
)


def _write_archive(
    path: Path, text: str, *, member: str = "virtual/traces.csv"
) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(member, text)


def _fixture_manifest(path: Path):
    return replace(
        WATRES_MANIFEST_V1,
        archive_md5=hashlib.md5(
            path.read_bytes()
        ).hexdigest(),  # nosec B324 - fixture checksum
        archive_max_bytes=max(path.stat().st_size + 1024, 4096),
    )


def _valid_metadata(manifest, *, archive_size: int = 1024):
    return {
        "id": manifest.record_id,
        "metadata": {
            "doi": manifest.doi,
            "version": manifest.release_version,
            "resource_type": {"id": "dataset"},
        },
        "files": [
            {
                "key": manifest.archive_name,
                "checksum": f"md5:{manifest.archive_md5}",
                "size": archive_size,
            }
        ],
    }


def test_manifest_is_frozen_and_explicitly_bounded_to_component_validation():
    assert WATRES_MANIFEST_V1.manifest_version == "hydrosheaf.watres-manifest/v1"
    assert WATRES_MANIFEST_V1.archive_name == "data.zip"
    assert WATRES_MANIFEST_V1.archive_md5 == "0154de3424de65f112ee9a409f2bf26f"
    assert "catchment-level temporal-TTD component" in WATRES_BENCHMARK_SCOPE
    assert "groundwater-graph validation" in WATRES_BENCHMARK_SCOPE
    assert "groundwater graph-topology validation" in WATRES_MANIFEST_V1.excluded_claims
    with pytest.raises(FrozenInstanceError):
        WATRES_MANIFEST_V1.archive_name = "other.zip"  # type: ignore[misc]


def test_metadata_validator_accepts_a_pinned_fixture_and_rejects_identity_drift():
    valid = validate_watres_metadata(_valid_metadata(WATRES_MANIFEST_V1))
    assert valid.valid
    assert valid.observed_archive_size_bytes == 1024
    assert valid.errors == ()

    drifted_payload = _valid_metadata(WATRES_MANIFEST_V1)
    drifted_payload["metadata"]["doi"] = "10.5281/zenodo.other"
    drifted_payload["files"][0]["checksum"] = "md5:00000000000000000000000000000000"
    drifted = validate_watres_metadata(drifted_payload)
    assert not drifted.valid
    assert any("DOI mismatch" in error for error in drifted.errors)
    assert any("MD5 mismatch" in error for error in drifted.errors)


def test_missing_archive_fails_honestly_without_network(tmp_path: Path):
    with pytest.raises(WATRESArchiveUnavailable, match="not bundled"):
        inspect_watres_archive(tmp_path / "not-present.zip")


def test_downloader_requires_explicit_authorization_before_any_network_or_file_write(
    tmp_path: Path,
):
    destination = tmp_path / "data.zip"
    with pytest.raises(WATRESDownloadRefused, match="never downloaded implicitly"):
        download_watres_archive(destination)
    assert not destination.exists()
    assert not destination.with_name("data.zip.part").exists()


def test_metadata_fetch_is_opt_in_too():
    with pytest.raises(WATRESDownloadRefused, match="metadata access is opt-in"):
        fetch_watres_metadata()


def test_explicit_csv_mapping_reads_a_verified_fixture_without_auto_schema_guessing(
    tmp_path: Path,
):
    archive_path = tmp_path / "fixture.zip"
    _write_archive(
        archive_path,
        "timestamp,input_d18O,output_d18O\n"
        "2024-01-01T00:00:00+00:00,-8.1,-8.0\n"
        "2024-01-02T00:00:00+00:00,-7.9,-8.05\n",
    )
    manifest = _fixture_manifest(archive_path)

    integrity = verify_watres_archive(archive_path, manifest=manifest)
    assert integrity.md5 == manifest.archive_md5

    inspection = inspect_watres_archive(
        archive_path,
        manifest=manifest,
        verify_checksum=True,
    )
    assert inspection.checksum_verified
    assert not inspection.automatic_schema_supported
    assert inspection.entries[0].name == "virtual/traces.csv"
    assert "Automatic member/column inference" in inspection.schema_note

    series = load_watres_csv(
        archive_path,
        WATRESCSVAdapterSpec(
            member="virtual/traces.csv",
            timestamp_column="timestamp",
            input_tracer_column="input_d18O",
            output_tracer_column="output_d18O",
        ),
        manifest=manifest,
    )
    assert series.input_values == (-8.1, -7.9)
    assert series.output_values == (-8.0, -8.05)
    assert series.benchmark_scope == WATRES_BENCHMARK_SCOPE


def test_explicit_adapter_rejects_unmapped_columns_instead_of_guessing(tmp_path: Path):
    archive_path = tmp_path / "fixture.zip"
    _write_archive(
        archive_path,
        "time,inflow\n2024-01-01T00:00:00,1.0\n",
    )
    manifest = _fixture_manifest(archive_path)
    with pytest.raises(WATRESSchemaUnsupported, match="missing columns"):
        load_watres_csv(
            archive_path,
            WATRESCSVAdapterSpec(
                member="virtual/traces.csv",
                timestamp_column="time",
                input_tracer_column="inflow",
                output_tracer_column="outflow",
            ),
            manifest=manifest,
        )


def test_explicit_mapping_can_retain_all_watres_hydrologic_input_series(
    tmp_path: Path,
):
    archive_path = tmp_path / "fixture.zip"
    _write_archive(
        archive_path,
        "time,precipitation,streamflow,pet,input,output\n"
        "2024-01-01T00:00:00,2.5,1.2,0.8,-8.1,-8.0\n"
        "2024-01-02T00:00:00,0.0,1.1,0.9,-7.9,-8.05\n",
    )
    manifest = _fixture_manifest(archive_path)
    series = load_watres_csv(
        archive_path,
        WATRESCSVAdapterSpec(
            member="virtual/traces.csv",
            timestamp_column="time",
            input_tracer_column="input",
            output_tracer_column="output",
            precipitation_column="precipitation",
            streamflow_column="streamflow",
            potential_evapotranspiration_column="pet",
        ),
        manifest=manifest,
    )
    assert series.precipitation_values == (2.5, 0.0)
    assert series.streamflow_values == (1.2, 1.1)
    assert series.potential_evapotranspiration_values == (0.8, 0.9)


def test_partial_hydrologic_mapping_is_rejected_before_reading_the_archive():
    with pytest.raises(ValueError, match="together"):
        WATRESCSVAdapterSpec(
            member="virtual/traces.csv",
            timestamp_column="time",
            input_tracer_column="input",
            output_tracer_column="output",
            precipitation_column="precipitation",
        )
