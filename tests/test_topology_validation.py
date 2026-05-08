from __future__ import annotations

from hydrosheaf.graph.types import Edge
from hydrosheaf.validation import (
    apply_modpath_informed_graph_priors,
    build_modpath_informed_graph_priors,
    edge_confusion,
    normalize_directed_edges,
    validate_independent_graph_against_modpath,
)


def test_normalize_directed_edges_accepts_common_edge_formats():
    edges = [
        ("A", "B"),
        {"u": "B", "v": "C"},
        {"edge_id": "C->D"},
        Edge(edge_id="D->E", u="D", v="E"),
    ]

    assert normalize_directed_edges(edges) == {
        ("A", "B"),
        ("B", "C"),
        ("C", "D"),
        ("D", "E"),
    }


def test_edge_confusion_reports_false_positives_and_false_negatives():
    report = edge_confusion(
        reference_edges=[("A", "B"), ("B", "C"), ("C", "D")],
        inferred_edges=[("A", "B"), ("B", "D"), ("C", "D")],
        candidate_edges=[
            ("A", "B"),
            ("B", "C"),
            ("C", "D"),
            ("B", "D"),
            ("A", "C"),
        ],
    )

    assert report["tp"] == 2.0
    assert report["fp"] == 1.0
    assert report["fn"] == 1.0
    assert report["false_positives"] == [("B", "D")]
    assert report["false_negatives"] == [("B", "C")]
    assert report["false_positive_rate"] > 0.0


def test_validate_independent_graph_against_modpath_keeps_mode_explicit():
    report = validate_independent_graph_against_modpath(
        inferred_edges=[("A", "B"), ("A", "C"), ("C", "D")],
        modpath_reference_edges=[("A", "B"), ("B", "C"), ("C", "D")],
        candidate_edges=[("A", "B"), ("A", "C"), ("B", "C"), ("C", "D")],
        edge_lengths={
            ("A", "B"): 100.0,
            ("B", "C"): 120.0,
            ("C", "D"): 130.0,
            ("A", "C"): 900.0,
        },
    )

    assert report["validation_mode"] == "independent_graph_inference"
    assert report["metrics"]["fp"] == 1.0
    assert report["metrics"]["fn"] == 1.0
    assert report["failure_modes"]["false_positive_edges"] == [("A", "C")]
    assert report["scale_mismatch"]["scale_mismatch"]


def test_modpath_informed_priors_are_not_independent_validation():
    priors = build_modpath_informed_graph_priors(
        [("A", "B"), ("B", "C")],
        travel_time_days={"A->B": 10.0, "B->C": 20.0},
    )

    assert [prior.edge_id() for prior in priors] == ["A->B", "B->C"]
    assert priors[0].tt_mean_days == 10.0

    report = apply_modpath_informed_graph_priors(
        hydrosheaf_edges=[("A", "B")],
        modpath_edges=[("A", "B"), ("B", "C")],
        mode="merge",
    )

    assert report["validation_mode"] == "modpath_informed_graph_prior"
    assert report["not_independent_validation"] is True
    assert set(report["output_edges"]) == {"A->B", "B->C"}
