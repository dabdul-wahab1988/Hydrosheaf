from hydrosheaf.validation.metrics import topology_metrics
from hydrosheaf.validation.programme_gate import score_topology


def test_topology_metrics_reports_false_positive_rate_for_full_universe() -> None:
    metrics = topology_metrics(
        [("A", "B")],
        [("A", "B"), ("A", "C")],
        candidate_edges=[("A", "B"), ("A", "C"), ("B", "C")],
    )

    assert metrics["tp"] == 1.0
    assert metrics["fp"] == 1.0
    assert metrics["tn"] == 1.0
    assert metrics["false_positive_rate"] == 0.5
    assert metrics["false_discovery_rate"] == 0.5


def test_score_topology_uses_same_candidate_universe_for_selection_metrics() -> None:
    metrics = score_topology(
        reference_edges=[("A", "B")],
        candidate_edges=[("A", "B"), ("A", "C"), ("B", "C")],
        selected_edges=[("A", "B"), ("A", "C")],
    )

    assert metrics["candidate_recall"] == 1.0
    assert metrics["selected_recall"] == 1.0
    assert metrics["selected_false_positive_rate"] == 0.5
    assert metrics["selected_false_discovery_rate"] == 0.5

