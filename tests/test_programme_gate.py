from __future__ import annotations

import pytest

from hydrosheaf.validation import (
    prepare_truth_blind_rows,
    programme_stages_from_status,
    required_stages_completed,
    score_age_posteriors,
    score_reaction_families,
    score_topology,
)


def test_prepare_rows_remaps_and_permutates_observed_tracer_only() -> None:
    rows = prepare_truth_blind_rows(
        [{"site_id": "A", "tritium_TU": 1.0}, {"site_id": "B", "tritium_TU": 2.0}],
        tracer_source="tritium_TU",
        tracer_target="3H",
        permute_tracer_seed=12,
    )

    assert {row["3H"] for row in rows} == {1.0, 2.0}
    assert all("true_" not in key for row in rows for key in row)

    with pytest.raises(ValueError, match="truth fields"):
        prepare_truth_blind_rows([{"site_id": "A", "true_age": 4.0}])


def test_stage_helper_requires_all_requested_stages() -> None:
    status = {
        "network_fit": {"status": "completed", "requested": True},
        "nuclear_age": {"status": "completed", "requested": True},
    }
    stages = programme_stages_from_status(status, name_prefix="case:")

    assert [stage.name for stage in stages] == ["case:network_fit", "case:nuclear_age"]
    assert required_stages_completed(status, ("network_fit", "nuclear_age"))
    assert not required_stages_completed(status, ("network_fit", "sheaf_refinement"))


def test_missing_truth_blind_instrumentation_fails_closed() -> None:
    stages = programme_stages_from_status(
        {"network_fit": {"status": "completed"}},
    )
    assert stages[0].truth_blind is False


def test_topology_and_age_scores_are_separate_and_auditable() -> None:
    topology = score_topology(
        [("A", "B"), ("B", "C")],
        [{"u": "A", "v": "B"}, {"u": "A", "v": "C"}],
        [("A", "B")],
    )
    age = score_age_posteriors(
        {"A": 5.0, "B": 10.0},
        {
            "A": {"mean_age_years": 5.5, "age_95_low": 2.0, "age_95_high": 9.0},
            "B": {"mean_age_years": 9.0, "age_95_low": 5.0, "age_95_high": 13.0},
        },
    )

    assert topology["candidate_recall"] == 0.5
    assert topology["selected_recall"] == 0.5
    assert age["status"] == "scored"
    assert age["n"] == 2
    assert age["all_finite"] is True


def test_reaction_family_score_uses_edge_results_without_truth_inference() -> None:
    class Result:
        def __init__(self, edge_id: str, label: str, extent: float) -> None:
            self.edge_id = edge_id
            self.z_labels = [label]
            self.z_extents = [extent]

    report = score_reaction_families(
        {"A->B": "carbonate_weathering", "B->C": "iron_reduction"},
        [Result("A->B", "calcite", 0.4), Result("B->C", "iron_reduction", 0.2)],
    )

    assert report["status"] == "scored"
    assert report["n"] == 2
    assert report["metrics"]["accuracy"] == 1.0
