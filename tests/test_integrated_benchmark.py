from __future__ import annotations

from pathlib import Path

from hydrosheaf.validation.integrated_benchmark import (
    audit_integrated_panels,
    audit_reference_panel,
    source_manifest_entry,
)


def test_model_reference_can_be_complete_without_being_field_truth() -> None:
    audit = audit_reference_panel(
        "modpath",
        "calibrated_model_reference",
        nodes=[{"node_id": "a"}, {"node_id": "b"}],
        edges=[{"u": "a", "v": "b"}],
        metadata={"model_name": "MODPATH", "source_doi": "10.0000/example"},
    )
    assert audit.status == "COMPLETE"
    assert audit.capabilities["directed_edge_table"] is True
    assert audit.capabilities["independent_direct_adjacency_truth"] is False


def test_separate_complete_panels_are_not_silently_integrated() -> None:
    first = audit_reference_panel(
        "age",
        "calibrated_model_reference",
        nodes=[{"node_id": "age:a"}],
        observations=[{"node_id": "age:a", "lat": 1.0, "lon": 2.0, "sample_date": "2020-01-01"}],
        metadata={
            "required_components": ("nodes", "observations"),
            "model_name": "LPM",
            "source_doi": "10.0000/age",
        },
    )
    second = audit_reference_panel(
        "flow",
        "calibrated_model_reference",
        nodes=[{"node_id": "flow:a"}],
        edges=[{"u": "flow:a", "v": "flow:b"}],
        metadata={"model_name": "MODPATH", "source_doi": "10.0000/flow"},
    )
    summary = audit_integrated_panels([first, second])
    assert summary["status"] == "PARTIAL"
    assert summary["integrated_scoring_allowed"] is False
    assert summary["crosswalk_declared"] is False


def test_missing_source_is_explicit(tmp_path: Path) -> None:
    entry = source_manifest_entry(
        tmp_path / "not-present.txt",
        source_id="missing",
        role="test",
    )
    assert entry["exists"] is False
    assert entry["sha256"] is None
