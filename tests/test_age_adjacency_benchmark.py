from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "M7" / "m7_nonuniqueness_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from age_adjacency_benchmark import (  # noqa: E402
    classify_candidate_pair,
    classify_candidate_relation,
    evaluate_age_adjacency,
    label_candidate_relations,
)


TRUTH_EDGES = (("A", "B"), ("B", "C"), ("C", "D"))
PATH_METADATA = (
    {"node_id": "A", "particle": 0, "milestone": 0},
    {"node_id": "B", "particle": 0, "milestone": 1},
    {"node_id": "C", "particle": 0, "milestone": 2},
    {"node_id": "D", "particle": 0, "milestone": 3},
    {"node_id": "X", "particle": 1, "milestone": 0},
    {"node_id": "Y", "particle": 1, "milestone": 1},
)


def test_relation_labels_separate_direct_reachability_and_direction() -> None:
    assert classify_candidate_relation("A", "B", TRUTH_EDGES) == "direct_adjacent"
    assert (
        classify_candidate_relation("A", "C", TRUTH_EDGES)
        == "transitive_reachable"
    )
    evidence = classify_candidate_pair("A", "D", TRUTH_EDGES)
    assert evidence.label == "transitive_reachable"
    assert evidence.truth_path_length == 3
    assert classify_candidate_relation("C", "A", TRUTH_EDGES) == "reverse_incompatible"


def test_path_metadata_identifies_cross_path_and_same_path_unrelated_pairs() -> None:
    cross_path = classify_candidate_pair(
        "A",
        "X",
        TRUTH_EDGES,
        path_metadata=PATH_METADATA,
    )
    assert cross_path.label == "cross_path"
    assert cross_path.status == "known"
    assert cross_path.same_truth_path is False

    same_path_unrelated = classify_candidate_pair(
        "X",
        "Y",
        TRUTH_EDGES,
        path_metadata=PATH_METADATA,
    )
    assert same_path_unrelated.label == "unrelated"
    assert same_path_unrelated.status == "known"
    assert same_path_unrelated.same_truth_path is True


def test_missing_path_metadata_is_reported_as_ambiguous() -> None:
    evidence = classify_candidate_pair("A", "X", TRUTH_EDGES)
    assert evidence.label == "unknown"
    assert evidence.status == "unknown"

    disconnected = classify_candidate_pair(
        "A",
        "E",
        (*TRUTH_EDGES, ("E", "F")),
    )
    assert disconnected.label == "unrelated"
    assert disconnected.status == "ambiguous"

    # Both endpoints occur in the supplied graph, so non-reachability can be
    # used as a negative for reachability, but cross-path versus unrelated is
    # not identifiable without path IDs.
    evidence = classify_candidate_pair("A", "D", TRUTH_EDGES)
    assert evidence.label == "transitive_reachable"
    evidence = classify_candidate_pair("A", "B", TRUTH_EDGES)
    assert evidence.status == "known"

    unrelated = classify_candidate_pair("B", "D", TRUTH_EDGES)
    assert unrelated.label == "transitive_reachable"
    unrelated = classify_candidate_pair("D", "B", TRUTH_EDGES)
    assert unrelated.label == "reverse_incompatible"


def test_labeling_uses_truth_only_and_parses_edge_ids() -> None:
    rows = [
        {"edge_id": "A->B", "direct_probability": 0.9, "age_cost": 0.0},
        {"edge_id": "A->C", "direct_probability": 0.8, "age_cost": 0.0},
    ]
    labelled = label_candidate_relations(rows, TRUTH_EDGES)
    assert [row["evaluation_relation"] for row in labelled] == [
        "direct_adjacent",
        "transitive_reachable",
    ]
    assert labelled[0]["edge_id"] == "A->B"
    assert labelled[0]["age_cost"] == 0.0
    assert labelled[1]["evaluation_truth_path_length"] == 2
    assert labelled[1]["evaluation_direct_target"] == 0
    assert labelled[1]["evaluation_reachable_target"] == 1


def test_metrics_expose_the_multistep_skip_failure_mode() -> None:
    # A -> C and A -> D are reachable but not direct.  Their high direct
    # probabilities intentionally reproduce the observed age-gate failure:
    # temporal compatibility accepts the skips instead of rejecting them.
    rows = label_candidate_relations(
        [
            {
                "edge_id": "A->B",
                "direct_probability": 0.90,
                "reachability_probability": 0.95,
                "direction_probability": 0.95,
            },
            {
                "edge_id": "A->C",
                "direct_probability": 0.80,
                "reachability_probability": 0.95,
                "direction_probability": 0.95,
            },
            {
                "edge_id": "A->D",
                "direct_probability": 0.70,
                "reachability_probability": 0.90,
                "direction_probability": 0.95,
            },
            {
                "edge_id": "C->A",
                "direct_probability": 0.10,
                "reachability_probability": 0.05,
                "direction_probability": 0.05,
            },
            {
                "edge_id": "A->X",
                "direct_probability": 0.20,
                "reachability_probability": 0.05,
                "direction_probability": None,
            },
            {
                "edge_id": "X->Y",
                "direct_probability": 0.20,
                "reachability_probability": 0.05,
                "direction_probability": None,
            },
        ],
        TRUTH_EDGES,
        path_metadata=PATH_METADATA,
    )
    metrics = evaluate_age_adjacency(
        rows,
        true_edges=TRUTH_EDGES,
        threshold=0.5,
    )

    assert metrics["n_relation_direct_adjacent"] == 1
    assert metrics["n_relation_transitive_reachable"] == 2
    assert metrics["n_relation_reverse_incompatible"] == 1
    assert metrics["n_relation_cross_path"] == 1
    assert metrics["n_relation_unrelated"] == 1
    assert metrics["truth_direct_count"] == 3
    assert metrics["candidate_direct_count"] == 1
    assert metrics["candidate_direct_recall"] == pytest.approx(1.0 / 3.0)
    assert metrics["truth_reachable_count"] == 6
    assert metrics["candidate_reachable_count"] == 3
    assert metrics["candidate_reachable_contained_count"] == 3
    assert metrics["candidate_reachable_recall"] == pytest.approx(0.5)

    # The same signal is directionally correct but cannot reject either skip.
    assert metrics["direction_n"] == 4
    assert metrics["direction_consistency"] == pytest.approx(1.0)
    assert metrics["false_skip_n"] == 2
    assert metrics["false_skip_rejection_rate"] == pytest.approx(0.0)
    assert metrics["false_skip_mean_direct_probability"] == pytest.approx(0.75)

    # Reachability treats both direct and skips as positive; direct adjacency
    # treats skips as negative.  These must be separate estimands.
    assert metrics["reachability_pr_auc"] is not None
    assert metrics["direct_adjacency_pr_auc"] is not None
    assert metrics["direct_adjacency_n"] == 6
    assert metrics["reachability_n"] == 6


def test_metrics_are_explicit_when_scores_are_missing_or_single_class() -> None:
    rows = label_candidate_relations(
        [{"edge_id": "A->B"}, {"edge_id": "A->C"}],
        TRUTH_EDGES,
    )
    metrics = evaluate_age_adjacency(rows)
    assert metrics["direct_adjacency_n"] == 0
    assert metrics["direct_adjacency_roc_auc"] is None
    assert metrics["direct_adjacency_pr_auc"] is None
    assert metrics["false_skip_rejection_rate"] is None
    assert metrics["direction_consistency"] is None

    # Brier/ECE are still defined with a single class if a probability is
    # available, but rank metrics are deliberately returned as undefined.
    only_direct = label_candidate_relations(
        [{"edge_id": "A->B", "direct_probability": 0.75}],
        TRUTH_EDGES,
    )
    single = evaluate_age_adjacency(only_direct)
    assert single["direct_adjacency_n"] == 1
    assert single["direct_adjacency_brier"] == pytest.approx(0.0625)
    assert single["direct_adjacency_roc_auc"] is None
    assert single["direct_adjacency_pr_auc"] is None
