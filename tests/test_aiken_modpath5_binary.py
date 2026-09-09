"""Focused tests for the bounded MODPATH 5 binary Aiken reader."""

from __future__ import annotations

import struct
from pathlib import Path
import zipfile

import pytest

from hydrosheaf.validation.aiken_reference import (
    MODPATH5_ENDPOINT_RECORD_BYTES,
    MODPATH5_HEADER_BYTES,
    MODPATH5_PATHLINE_RECORD_BYTES,
    decode_modpath5_ipcode,
    load_aiken_reference,
    parse_modpath5_endpoint_binary,
    parse_modpath5_pathline_binary,
)


def _header(reference_time: float = 42003.0) -> bytes:
    return b"MODPATH 5.0".ljust(80, b" ") + struct.pack("<f", reference_time)


def test_pathline_reader_preserves_raw_local_global_coordinates_and_time_flag() -> None:
    payload = _header(17.5) + b"".join(
        [
            struct.pack("<ifffffii", 7, 11.5, 12.5, 0.25, 100.0, 0.0, 321, 4),
            struct.pack("<ifffffii", 7, 11.75, 12.75, 0.30, 101.0, -2.5, 322, 4),
            struct.pack("<ifffffii", 9, 13.0, 14.0, 0.40, 102.0, 5.0, 323, 5),
        ]
    )

    frame = parse_modpath5_pathline_binary(
        payload,
        source_archive="output.zip",
        source_member="output/output.AK-455_MP/AK-455.pth",
    )

    assert len(payload) == MODPATH5_HEADER_BYTES + 3 * MODPATH5_PATHLINE_RECORD_BYTES
    assert list(frame["particle_ordinal"]) == [1, 1, 2]
    assert list(frame["particle_id"]) == [7, 7, 9]
    assert set(frame["particle_ordinal_basis"]) == {"first_seen_native_pathline_particle_id"}
    assert frame.loc[1, "time_raw"] == pytest.approx(-2.5)
    assert frame.loc[1, "tracking_time"] == pytest.approx(2.5)
    assert bool(frame.loc[1, "time_is_intermediate_requested_point"])
    assert frame.loc[0, "z_local"] == pytest.approx(0.25)
    assert frame.loc[0, "z_global"] == pytest.approx(100.0)
    assert frame.loc[0, "global_node"] == 321
    assert frame.loc[0, "run_id"] == "AK-455_MP"
    assert frame.loc[0, "source_member"].endswith("AK-455.pth")
    assert set(frame["binary_parse_status"]) == {"validated"}


def test_endpoint_reader_preserves_ipcode_censoring_and_release_z_local() -> None:
    payload = _header(42003.0) + struct.pack(
        "<ii7fiiiif",
        2,  # final zone
        999,  # final global node
        91.0,
        92.0,
        0.75,  # final local z; global final z is not in MODPATH 5 binary
        365.0,
        81.0,
        82.0,
        0.125,  # release local z; retain without reconstructing global z
        111,  # release global node
        3,  # release zone
        6,  # release cumulative timestep
        212,  # IDCODE 2, NSLAST 21 (10*NSLAST + IDCODE)
        10.0,  # release time
    )

    frame = parse_modpath5_endpoint_binary(
        payload,
        source_archive="output.zip",
        source_member="output/output.AK-455_MP/AK-455.ept",
    )

    assert len(payload) == MODPATH5_HEADER_BYTES + MODPATH5_ENDPOINT_RECORD_BYTES
    row = frame.iloc[0]
    assert row["particle_ordinal"] == 1
    assert row["particle_ordinal_basis"] == "endpoint_record_row_order_crosswalk_only"
    assert row["particle_id"] is None
    assert row["particle_id_status"] == "not_present_in_modpath5_endpoint"
    assert not bool(row["particle_id_derived"])
    assert row["final_node"] == 999
    assert row["release_node"] == 111
    assert row["release_z_local"] == pytest.approx(0.125)
    assert row["zloc0"] == pytest.approx(0.125)
    assert row["final_z_local"] == pytest.approx(0.75)
    assert row["z"] is None
    assert row["ipcode"] == 212
    assert row["ipcode_idcode"] == 2
    assert row["ipcode_nslast"] == 21
    assert row["termination_status"] == "stopped_in_specified_zone"
    assert row["particle_censoring_status"] == "observed_termination"
    assert row["binary_endian"] == "little"
    assert row["binary_parser_version"] == "1.0"
    assert row["record_bytes"] == MODPATH5_ENDPOINT_RECORD_BYTES
    assert row["ipcode_encoding"] == "nslast_times_10_plus_idcode"
    assert row["run_id"] == "AK-455_MP"


def test_aiken_loader_joins_endpoint_status_by_particle_ordinal(tmp_path: Path) -> None:
    endpoint_payload = _header() + b"".join(
        [
            struct.pack(
                "<ii7fiiiif",
                1,
                100,
                1.0,
                2.0,
                0.1,
                5.0,
                3.0,
                4.0,
                0.2,
                10,
                1,
                2,
                51,
                0.0,
            ),
            struct.pack(
                "<ii7fiiiif",
                1,
                101,
                1.0,
                2.0,
                0.1,
                6.0,
                3.0,
                4.0,
                0.3,
                11,
                1,
                2,
                0,
                0.0,
            ),
        ]
    )
    pathline_payload = _header() + b"".join(
        [
            struct.pack("<ifffffii", 9, 3.0, 4.0, 0.2, 10.0, 0.0, 10, 2),
            struct.pack("<ifffffii", 9, 3.5, 4.5, 0.25, 10.5, 1.0, 11, 2),
            struct.pack("<ifffffii", 4, 3.0, 4.0, 0.3, 20.0, 0.0, 11, 2),
        ]
    )
    root = tmp_path / "aiken"
    root.mkdir()
    with zipfile.ZipFile(root / "output.zip", "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("output/output.AK-455_MP/AK-455.ept", endpoint_payload)
        archive.writestr("output/output.AK-455_MP/AK-455.pth", pathline_payload)

    reference = load_aiken_reference(root)

    assert reference.endpoints["particle_id"].isna().all()
    assert list(reference.pathlines["particle_ordinal"]) == [1, 1, 2]
    assert list(reference.pathlines["ipcode"]) == [51, 51, 0]
    assert list(reference.pathlines["particle_censoring_status"]) == [
        "observed_termination",
        "observed_termination",
        "right_censored_active",
    ]
    assert set(reference.binary_parse_diagnostics["status"]) == {"validated"}


@pytest.mark.parametrize(
    ("ipcode", "idcode", "termination", "censoring"),
    [
        (-2, -2, "unreleased", "unreleased"),
        (-1, -1, "stranded_inactive_dry_cell", "censored_dry_cell"),
        (0, 0, "active_at_stop_time", "right_censored_active"),
        (10, 0, "active_at_stop_time", "right_censored_active"),
        (51, 1, "discharged_normally", "observed_termination"),
    ],
)
def test_ipcode_decoding_retains_explicit_censoring_boundary(
    ipcode: int,
    idcode: int,
    termination: str,
    censoring: str,
) -> None:
    parsed = decode_modpath5_ipcode(ipcode)
    assert parsed["ipcode_idcode"] == idcode
    assert parsed["termination_status"] == termination
    assert parsed["particle_censoring_status"] == censoring


def test_aiken_ipcode_uses_units_digit_for_idcode() -> None:
    parsed = decode_modpath5_ipcode(540)
    assert parsed["ipcode_idcode"] == 0
    assert parsed["ipcode_nslast"] == 54
    assert parsed["termination_status"] == "active_at_stop_time"


def test_binary_header_and_payload_validation_rejects_shifted_or_wrong_files() -> None:
    valid_pathline = _header() + struct.pack("<ifffffii", 1, 1.0, 2.0, 0.5, 3.0, 4.0, 5, 6)

    with pytest.raises(ValueError, match="unsupported header"):
        parse_modpath5_pathline_binary(b"MODPATH 6.0".ljust(80, b" ") + struct.pack("<f", 0.0) + valid_pathline[84:])

    with pytest.raises(ValueError, match="invalid payload length"):
        parse_modpath5_pathline_binary(valid_pathline + b"x")

    with pytest.raises(ValueError, match="non-finite reference time"):
        parse_modpath5_pathline_binary(
            b"MODPATH 5.0".ljust(80, b" ") + struct.pack("<f", float("inf")) + valid_pathline[84:]
        )
