from __future__ import annotations

import pytest
from unittest.mock import patch

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.network_fit import infer_null_aware_edges
from hydrosheaf.inference.network_fit import infer_edges
from hydrosheaf.inference.null_aware import (
    DECISION_ABSENT,
    DECISION_PRESENT,
    NullAwareLogisticCalibrator,
    build_feature_rows,
)


def _edge(edge_id: str, head: float, null_chemistry: float) -> Edge:
    return Edge(
        edge_id=edge_id,
        u="A",
        v=edge_id,
        attrs={
            "flow_features": {
                "head_drop_m": head,
                "age_direction_support": 0.9 if head > 0 else 0.1,
                "transport_support": 0.8 if head > 0 else 0.2,
            },
            "null_features": {
                "null_chemistry_similarity": null_chemistry,
                "common_lithology": 1.0,
                "null_score": null_chemistry,
            },
        },
    )


def _calibrator() -> NullAwareLogisticCalibrator:
    candidates = [
        _edge("F1", 8.0, 0.05),
        _edge("F2", 7.0, 0.10),
        _edge("F3", 6.0, 0.08),
        _edge("N1", -1.0, 0.90),
        _edge("N2", 0.0, 0.85),
        _edge("N3", -2.0, 0.95),
    ]
    rows = build_feature_rows(candidates)
    return NullAwareLogisticCalibrator(
        l2=0.2,
        fit_scope="held_out_calibration",
        calibration_provenance={
            "independent": True,
            "generator_id": "unit-test-generator",
            "split_id": "held-out-calibration",
            "dataset_hash": "unit-test-hash",
        },
    ).fit(
        rows,
        [1, 1, 1, 0, 0, 0],
    )


def test_null_aware_pipeline_scores_every_candidate_and_selects_only_present() -> None:
    candidates = [_edge("F", 8.0, 0.05), _edge("N", -1.0, 0.95)]
    selected, report = infer_null_aware_edges(
        {},
        candidates,
        Config(
            null_aware_present_threshold=0.60,
            null_aware_absent_threshold=0.40,
            null_aware_null_reject_threshold=0.80,
            edge_max_neighbors=1,
        ),
        topology_calibrator=_calibrator(),
        max_neighbors=1,
    )

    assert report["candidate_count"] == 2
    assert len(report["records"]) == 2
    assert {row["edge_id"] for row in report["records"]} == {"F", "N"}
    assert selected
    assert selected[0].edge_id == "F"
    assert selected[0].attrs["null_aware_decision"] == DECISION_PRESENT
    assert report["records"][1]["decision"] == DECISION_ABSENT


def test_null_aware_pipeline_fails_closed_without_calibration() -> None:
    with pytest.raises(ValueError, match="requires an explicit fitted calibrator"):
        infer_null_aware_edges(
            {},
            [_edge("F", 8.0, 0.05)],
            Config(null_aware_require_calibration=True),
            topology_calibrator=None,
        )


def test_optional_diagnostic_mode_reports_abstention_without_calibration() -> None:
    selected, report = infer_null_aware_edges(
        {},
        [_edge("F", 8.0, 0.05)],
        Config(
            null_aware_require_calibration=False,
            null_aware_require_deployable_calibration=False,
        ),
        topology_calibrator=None,
    )

    assert selected == []
    assert report["abstained_count"] == 1
    assert report["records"][0]["null_explanation_score"] is None


def test_null_aware_pipeline_rejects_development_only_calibration() -> None:
    training_edges = [
        _edge("F1", 8.0, 0.05),
        _edge("F2", 7.0, 0.10),
        _edge("F3", 6.0, 0.08),
        _edge("N1", -1.0, 0.90),
        _edge("N2", 0.0, 0.85),
        _edge("N3", -2.0, 0.95),
    ]
    development = NullAwareLogisticCalibrator(l2=0.2).fit(
        build_feature_rows(training_edges),
        [1, 1, 1, 0, 0, 0],
    )
    with pytest.raises(ValueError, match="held-out/deployment"):
        infer_null_aware_edges(
            {},
            [_edge("F", 8.0, 0.05)],
            Config(),
            topology_calibrator=development,
        )


def test_null_aware_pipeline_does_not_mutate_input_edges() -> None:
    candidate = _edge("F", 8.0, 0.05)
    original_attrs = dict(candidate.attrs)
    infer_null_aware_edges(
        {},
        [candidate],
        Config(
            null_aware_present_threshold=0.60,
            null_aware_absent_threshold=0.40,
        ),
        topology_calibrator=_calibrator(),
    )
    assert candidate.attrs == original_attrs


def test_public_infer_edges_exposes_null_aware_method() -> None:
    candidates = [_edge("F", 8.0, 0.05), _edge("N", -1.0, 0.95)]
    config = Config(
        phreeqc_enabled=False,
        edge_max_neighbors=1,
        edge_map_candidate_multiplier=1,
        null_aware_present_threshold=0.60,
        null_aware_absent_threshold=0.40,
    )
    with patch(
        "hydrosheaf.inference.network_fit.infer_edges_probabilistic",
        return_value=candidates,
    ):
        selected = infer_edges(
            [{"site_id": "A"}],
            method="null_aware_sheaf",
            config=config,
            topology_calibrator=_calibrator(),
        )

    assert [edge.edge_id for edge in selected] == ["F"]
    assert selected[0].attrs["null_aware_calibration_status"] == "FITTED"


def test_public_infer_edges_can_return_the_complete_candidate_audit() -> None:
    candidates = [_edge("F", 8.0, 0.05), _edge("N", -1.0, 0.95)]
    config = Config(
        phreeqc_enabled=False,
        edge_max_neighbors=1,
        edge_map_candidate_multiplier=1,
        null_aware_present_threshold=0.60,
        null_aware_absent_threshold=0.40,
    )
    with patch(
        "hydrosheaf.inference.network_fit.infer_edges_probabilistic",
        return_value=candidates,
    ):
        selected, report = infer_edges(
            [{"site_id": "A"}],
            method="null_aware_sheaf",
            config=config,
            topology_calibrator=_calibrator(),
            return_report=True,
        )

    assert [edge.edge_id for edge in selected] == ["F"]
    assert report["candidate_count"] == 2
    assert {record["edge_id"] for record in report["records"]} == {"F", "N"}
