from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.models.evidence_lifted import (
    evidence_lifted_resolution,
    reaction_identifiability_diagnostics,
    reaction_panel_diagnostics,
)


def _carbonate_matrix() -> tuple[np.ndarray, list[str]]:
    # Calcite and dolomite are exactly equivalent in this deliberately reduced
    # observed panel; albite supplies a distinct non-carbonate signature.
    return (
        np.asarray(
            [
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
            ]
        ),
        ["calcite", "dolomite", "albite"],
    )


def test_panel_reports_rank_nullity_coherence_and_abstains_on_carbonate_class() -> None:
    matrix, labels = _carbonate_matrix()
    report = reaction_panel_diagnostics(
        matrix,
        labels,
        ion_order=["Ca", "HCO3", "SiO2"],
        observed_ions=["Ca", "HCO3"],
    )

    assert report["status"] == "ABSTAIN"
    assert report["carbonate_status"] == "ABSTAIN"
    assert report["observed_ion_panel"] == ["Ca", "HCO3"]
    assert report["n_observed_ions"] == 2
    assert report["rank"] == 2
    assert report["nullity"] == 1
    assert report["rank_deficient"] is False
    assert report["ambiguous_equivalence_class_count"] == 1
    assert report["unresolved_carbonate_classes"][0]["carbonate_members"] == [
        "calcite",
        "dolomite",
    ]


def test_independent_evidence_can_lift_a_carbonate_class_without_changing_structure() -> None:
    matrix, labels = _carbonate_matrix()
    report = reaction_panel_diagnostics(
        matrix,
        labels,
        ion_order=["Ca", "HCO3", "SiO2"],
        observed_ions=["Ca", "HCO3"],
        evidence_scores={"calcite": 0.9, "dolomite": 0.1, "albite": 0.5},
        evidence_sources={
            "calcite": ["independent_tracer"],
            "dolomite": ["independent_tracer"],
        },
    )

    assert report["status"] == "RUN"
    assert report["carbonate_status"] == "NO_UNRESOLVED_CLASS"
    assert report["structural_status"] == "equivalence_classes_present"
    resolution = next(
        row
        for row in report["evidence_lifted_resolution"]
        if row["n_members"] == 2
    )
    assert resolution["resolution_status"] == "evidence_lifted_resolved"
    assert resolution["independent_evidence_available"] is True


def test_relative_scores_without_independent_source_do_not_clear_abstention() -> None:
    matrix, labels = _carbonate_matrix()
    report = reaction_panel_diagnostics(
        matrix,
        labels,
        ion_order=["Ca", "HCO3", "SiO2"],
        evidence_scores={"calcite": 0.9, "dolomite": 0.1, "albite": 0.5},
    )

    assert report["status"] == "ABSTAIN"
    assert report["evidence_status"] == "relative_scores_only"
    assert report["unresolved_carbonate_classes"]


def test_missing_evidence_score_is_neutral_and_recorded() -> None:
    result = evidence_lifted_resolution(
        ["calcite", "dolomite"],
        {"calcite": None},
    )

    assert result.evidence_scores["calcite"] == pytest.approx(0.5)
    assert result.evidence_scores["dolomite"] == pytest.approx(0.5)
    assert result.evidence_lifted_resolution_index == pytest.approx(0.0)
    assert result.missing_evidence_members == ("calcite", "dolomite")


def test_identifiability_diagnostics_reject_nonfinite_matrix() -> None:
    with pytest.raises(ValueError, match="finite"):
        reaction_identifiability_diagnostics(
            [[1.0, np.nan]],
            ["calcite"],
        )


def test_panel_validates_explicit_observed_ion_schema() -> None:
    matrix, labels = _carbonate_matrix()
    with pytest.raises(ValueError, match="ion_order length"):
        reaction_panel_diagnostics(matrix, labels, ion_order=["Ca"])
    with pytest.raises(ValueError, match="absent from ion_order"):
        reaction_panel_diagnostics(
            matrix,
            labels,
            ion_order=["Ca", "HCO3", "SiO2"],
            observed_ions=["Cl"],
        )
