"""Tests for the dependency-light age-adjacency result-table audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "scripts"
    / "age_adjacency_protocol.py"
)
SPEC = importlib.util.spec_from_file_location("age_adjacency_protocol", SCRIPT)
assert SPEC and SPEC.loader
age_protocol = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(age_protocol)


FIELDS = [
    "edge_id",
    "u",
    "v",
    "seed",
    "split",
    "is_true_edge",
    "probability_age_adjacency",
    "relation_label",
]


def _rows() -> list[dict[str, str]]:
    return [
        {
            "edge_id": "case1-a",
            "u": "n0",
            "v": "n1",
            "seed": "1",
            "split": "locked_test",
            "is_true_edge": "1",
            "probability_age_adjacency": "0.90",
            "relation_label": "direct_adjacent",
        },
        {
            "edge_id": "case1-b",
            "u": "n0",
            "v": "n2",
            "seed": "1",
            "split": "locked_test",
            "is_true_edge": "0",
            "probability_age_adjacency": "0.80",
            "relation_label": "transitive_reachable",
        },
        {
            "edge_id": "case1-c",
            "u": "n2",
            "v": "n1",
            "seed": "1",
            "split": "locked_test",
            "is_true_edge": "0",
            "probability_age_adjacency": "0.10",
            "relation_label": "unrelated",
        },
        {
            "edge_id": "case2-a",
            "u": "n0",
            "v": "n1",
            "seed": "2",
            "split": "locked_test",
            "is_true_edge": "1",
            "probability_age_adjacency": "0.40",
            "relation_label": "direct_adjacent",
        },
    ]


def test_audit_separates_candidate_and_all_truth_denominators() -> None:
    result = age_protocol.audit_rows(
        _rows(),
        fieldnames=FIELDS,
        score_column="probability_age_adjacency",
        all_truth_count=3,
        threshold=0.5,
    )

    assert result["n_rows"] == 4
    assert result["n_cases"] == 2
    assert result["n_positive_candidates"] == 2
    assert result["all_truth_count"] == 3
    assert result["candidate_recall_against_all_truth"] == pytest.approx(2 / 3)
    assert result["splits"] == {"locked_test": 4}
    assert result["relation_counts"] == {
        "direct_adjacent": 2,
        "transitive_reachable": 1,
        "unrelated": 1,
    }
    assert result["classification"]["f1"] == pytest.approx(0.5)
    assert result["calibration"]["brier"] == pytest.approx(
        ((0.9 - 1) ** 2 + (0.8 - 0) ** 2 + (0.1 - 0) ** 2 + (0.4 - 1) ** 2) / 4
    )


def test_audit_infers_fieldnames_for_in_memory_rows() -> None:
    result = age_protocol.audit_rows(
        _rows(),
        score_column="probability_age_adjacency",
    )
    assert result["n_rows"] == 4


def test_audit_rejects_truth_bearing_prediction_name() -> None:
    with pytest.raises(age_protocol.ProtocolAuditError, match="truth-bearing"):
        age_protocol.audit_rows(
            _rows(),
            fieldnames=FIELDS,
            score_column="direct_adjacent_probability",
        )


def test_audit_rejects_out_of_range_probability_and_duplicate_edge() -> None:
    bad_score = _rows()
    bad_score[0] = dict(bad_score[0], probability_age_adjacency="1.1")
    with pytest.raises(age_protocol.ProtocolAuditError, match=r"in \[0, 1\]"):
        age_protocol.audit_rows(
            bad_score,
            fieldnames=FIELDS,
            score_column="probability_age_adjacency",
        )

    duplicate = _rows() + [dict(_rows()[0], edge_id="case1-a-copy")]
    with pytest.raises(age_protocol.ProtocolAuditError, match="duplicate candidate edge"):
        age_protocol.audit_rows(
            duplicate,
            fieldnames=FIELDS,
            score_column="probability_age_adjacency",
        )


def test_cli_writes_audit_json(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    csv_path = tmp_path / "edges.csv"
    csv_path.write_text(
        ",".join(FIELDS)
        + "\n"
        + "case1-a,n0,n1,1,locked_test,1,0.9,direct_adjacent\n"
        + "case1-b,n0,n2,1,locked_test,0,0.2,transitive_reachable\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "audit.json"

    assert (
        age_protocol.main(
            [
                "--input",
                str(csv_path),
                "--score-column",
                "probability_age_adjacency",
                "--truth-count",
                "1",
                "--output",
                str(output_path),
            ]
        )
        == 0
    )
    assert "\"valid\": true" in capsys.readouterr().out
    written = output_path.read_text(encoding="utf-8")
    assert '"candidate_recall_against_all_truth": 1.0' in written
