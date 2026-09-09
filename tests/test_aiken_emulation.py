"""Tests for the isolated Aiken model-conditioned emulation."""

from __future__ import annotations

from pathlib import Path
import hashlib

import pandas as pd
import pytest

from hydrosheaf.validation.aiken_emulation import (
    build_cfc_age_intervals,
    build_transport_hypotheses,
    build_well_model_crosswalk,
    run_aiken_model_conditioned_emulation,
    time_unit_conversion_to_years,
    write_aiken_emulation_outputs,
)


def _frames() -> dict[str, pd.DataFrame]:
    wells = pd.DataFrame(
        [
            {"well_id": "AK-1", "node_id": "AIKEN:AK-1"},
            {"well_id": "AK-2", "node_id": "AIKEN:AK-2"},
        ]
    )
    locations = pd.DataFrame(
        [
            {
                "well_id": "AK-1",
                "model_row": 10,
                "model_column": 20,
                "model_layer": 3,
                "local_x": 0.25,
                "local_y": 0.50,
                "local_z": 0.75,
                "source_member": "model/AK-1_MP/AK-1.loc",
            }
        ]
    )
    references = pd.DataFrame(
        [
            {
                "simulation_id": "AK-1_MP",
                "time_units": "days",
                "source_member": "model/AK-1_MP/usgs.model.reference",
            },
            {
                "simulation_id": "AK-2_MP",
                "time_units": "days",
                "source_member": "model/AK-2_MP/usgs.model.reference",
            },
        ]
    )
    endpoints = pd.DataFrame(
        [
            {
                "run_id": "AK-1_MP",
                "well_id": "AK-1",
                "particle_ordinal": 1,
                "particle_id": None,
                "particle_id_status": "not_present_in_modpath5_endpoint",
                "release_node": 100,
                "final_node": 200,
                "total_tracking_time": 365.25,
                "release_time": 0.0,
                "ipcode": 51,
                "particle_censoring_status": "observed_termination",
                "termination_status": "discharged_normally",
                "tracking_direction_code": 2,
                "tracking_direction": "backward_toward_recharge",
                "reference_time": 0.0,
                "source_member": "output/output.AK-1_MP/AK-1.ept",
                "source_archive": "output.zip",
            },
            {
                "run_id": "AK-1_MP",
                "well_id": "AK-1",
                "particle_ordinal": 2,
                "particle_id": None,
                "particle_id_status": "not_present_in_modpath5_endpoint",
                "release_node": 100,
                "final_node": 201,
                "total_tracking_time": 730.5,
                "release_time": 0.0,
                "ipcode": 10,
                "particle_censoring_status": "right_censored_active",
                "termination_status": "active_at_stop_time",
                "tracking_direction_code": 2,
                "tracking_direction": "backward_toward_recharge",
                "reference_time": 0.0,
                "source_member": "output/output.AK-1_MP/AK-1.ept",
                "source_archive": "output.zip",
            },
        ]
    )
    pathlines = pd.DataFrame(
        [
            {
                "run_id": "AK-1_MP",
                "well_id": "AK-1",
                "particle_ordinal": 1,
                "particle_id": 7,
                "record_ordinal": 1,
                "node": 100,
                "time": 0.0,
                "tracking_time": 0.0,
                "tracking_direction": "backward_toward_recharge",
                "tracking_direction_code": 2,
                "source_member": "output/output.AK-1_MP/AK-1.pth",
                "source_archive": "output.zip",
            },
            {
                "run_id": "AK-1_MP",
                "well_id": "AK-1",
                "particle_ordinal": 1,
                "particle_id": 7,
                "record_ordinal": 2,
                "node": 150,
                "time": 365.25,
                "tracking_time": 365.25,
                "tracking_direction": "backward_toward_recharge",
                "tracking_direction_code": 2,
                "source_member": "output/output.AK-1_MP/AK-1.pth",
                "source_archive": "output.zip",
            },
            {
                "run_id": "AK-1_MP",
                "well_id": "AK-1",
                "particle_ordinal": 1,
                "particle_id": 7,
                "record_ordinal": 3,
                "node": 200,
                "time": 730.5,
                "tracking_time": 730.5,
                "tracking_direction": "backward_toward_recharge",
                "tracking_direction_code": 2,
                "source_member": "output/output.AK-1_MP/AK-1.pth",
                "source_archive": "output.zip",
            },
        ]
    )
    cfc = pd.DataFrame(
        [
            {
                "well_id": "AK-1",
                "node_id": "AIKEN:AK-1",
                "sample_date": "2015-07-01",
                "apparent_age_label": "Early 2010's",
                "recharge_year_min": 2010,
                "recharge_year_max": 2013,
                "source_sheet": "Table 13",
            }
        ]
    )
    return {
        "well_metadata": wells,
        "model_locations": locations,
        "model_references": references,
        "endpoints": endpoints,
        "pathlines": pathlines,
        "cfc_ages": cfc,
    }


def test_time_conversion_is_explicit_and_fails_closed() -> None:
    conversion = time_unit_conversion_to_years("days")
    assert conversion["conversion_status"] == "VALIDATED"
    assert conversion["years_per_source_unit"] == pytest.approx(1.0 / 365.25)
    assert conversion["conversion_rule"] == "source days / 365.25"
    assert time_unit_conversion_to_years("fortnights")["conversion_status"] == "ABSTAIN_UNKNOWN_TIME_UNIT"


def test_crosswalk_requires_exact_loc_and_does_not_guess_endpoint_cell() -> None:
    frames = _frames()
    endpoint_only = frames["endpoints"].iloc[[0]].copy()
    endpoint_only["run_id"] = "AK-2_MP"
    endpoint_only["well_id"] = "AK-2"
    endpoint_only["source_member"] = "output/output.AK-2_MP/AK-2.ept"
    endpoint_only["particle_ordinal"] = 1
    endpoints = pd.concat([frames["endpoints"], endpoint_only], ignore_index=True)
    crosswalk = build_well_model_crosswalk(
        frames["well_metadata"],
        frames["model_locations"],
        frames["model_references"],
        endpoints,
    )
    explicit = crosswalk[crosswalk["run_id"] == "AK-1_MP"].iloc[0]
    endpoint_only = crosswalk[crosswalk["run_id"] == "AK-2_MP"].iloc[0]
    assert explicit["crosswalk_status"] == "EXPLICIT_MODEL_STARTING_CELL"
    assert bool(explicit["crosswalk_is_explicit"])
    assert endpoint_only["crosswalk_status"] == "NOT_ESTABLISHED_ENDPOINT_ONLY"
    assert not bool(endpoint_only["crosswalk_is_explicit"])


def test_cfc_label_becomes_screening_interval_not_independent_truth() -> None:
    interval = build_cfc_age_intervals(_frames()["cfc_ages"])
    row = interval.iloc[0]
    assert row["age_low_years"] == pytest.approx(2.0)
    assert row["age_high_years"] == pytest.approx(5.0)
    assert row["age_interval_status"] == "CFC_APPARENT_AGE_SCREENING_INTERVAL"
    assert row["independent_age_truth"] is False or not bool(row["independent_age_truth"])


def test_emulation_separates_direct_endpoint_and_indirect_segments_with_censoring() -> None:
    frames = _frames()
    crosswalk = build_well_model_crosswalk(
        frames["well_metadata"],
        frames["model_locations"],
        frames["model_references"],
        frames["endpoints"],
    )
    cfc = build_cfc_age_intervals(frames["cfc_ages"])
    direct, segments, summary = build_transport_hypotheses(
        frames["endpoints"],
        frames["pathlines"],
        frames["model_references"],
        crosswalk,
        cfc,
    )
    assert set(direct["hypothesis_kind"]) == {"direct_endpoint_transport"}
    assert list(direct["travel_time_years"]) == pytest.approx([1.0, 2.0])
    assert list(direct["event_observed"]) == [True, False]
    assert direct.loc[1, "travel_time_observation_status"] == "right_censored_lower_bound"
    assert len(segments) == 2
    assert set(segments["hypothesis_kind"]) == {"indirect_pathline_segment_transport"}
    assert list(segments["segment_time_years"]) == pytest.approx([1.0, 1.0])
    direct_summary = summary[summary["hypothesis_kind"] == "direct_endpoint_transport"].iloc[0]
    assert direct_summary["n_observed_termination"] == 1
    assert direct_summary["n_right_censored"] == 1
    assert direct_summary["time_conversion_status"] == "VALIDATED"
    assert direct_summary["integrated_field_scoring_allowed"] is False or not bool(direct_summary["integrated_field_scoring_allowed"])
    assert direct_summary["direct_adjacency_truth_status"] == "ABSTAIN"


def test_emulation_manifest_and_writer_mark_field_scoring_prohibited(tmp_path: Path) -> None:
    result = run_aiken_model_conditioned_emulation(_frames())
    assert result.manifest["reference_type"] == "calibrated_model_reference"
    assert result.manifest["independence_and_claims"]["integrated_field_scoring_allowed"] is False
    assert result.manifest["independence_and_claims"]["direct_well_to_well_truth_emitted"] is False
    assert result.manifest["crosswalk"]["status"] == "EXPLICIT"
    generated = write_aiken_emulation_outputs(result, tmp_path / "output")
    assert (tmp_path / "output" / generated["manifest"]).exists()
    assert (tmp_path / "output" / generated["aiken_direct_endpoint_hypotheses"]).exists()
    output_record = result.manifest["output_files"]["aiken_direct_endpoint_hypotheses"]
    output_path = tmp_path / "output" / output_record["path"]
    assert len(output_record["sha256"]) == 64
    assert output_record["size_bytes"] == output_path.stat().st_size
    digest = hashlib.sha256(output_path.read_bytes()).hexdigest()
    assert output_record["sha256"] == digest
