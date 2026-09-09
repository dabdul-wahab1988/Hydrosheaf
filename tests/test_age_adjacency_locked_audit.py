from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "scripts"
    / "age_adjacency_locked_audit.py"
)
SPEC = importlib.util.spec_from_file_location("age_adjacency_locked_audit", SCRIPT)
assert SPEC and SPEC.loader
audit_module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit_module)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_locked_audit_joins_truth_after_inference_and_isolates_output(tmp_path: Path) -> None:
    source = tmp_path / "locked"
    case = source / "cases" / "locked_test_1"
    _write_csv(
        source / "locked_test_edge_results.csv",
        ["edge_id", "u", "v", "seed", "split", "age_cost"],
        [
            {
                "edge_id": "A->B",
                "u": "A",
                "v": "B",
                "seed": "1",
                "split": "locked_test",
                "age_cost": 0.0,
            }
        ],
    )
    _write_csv(
        case / "heldout_truth.csv",
        ["edge_id", "u", "v", "true_process"],
        [{"edge_id": "A->B", "u": "A", "v": "B", "true_process": "none"}],
    )
    _write_csv(
        case / "modpath_pathline_truth.csv",
        ["node_id", "particle", "milestone"],
        [
            {"node_id": "A", "particle": 0, "milestone": 0},
            {"node_id": "B", "particle": 0, "milestone": 1},
        ],
    )
    (source / "manifest.json").write_text("{}\n", encoding="utf-8")

    destination = tmp_path / "audit"
    result = audit_module.audit_locked_results(source, destination)

    assert result["locked_results_modified"] is False
    assert result["metrics"]["n_relation_direct_adjacent"] == 1
    assert result["metrics"]["candidate_direct_recall"] == 1.0
    assert result["metrics"]["score_interpretation"].startswith(
        "historical age compatibility"
    )
    labelled = (destination / "relation_labeled_locked_test_edge_results.csv").read_text(
        encoding="utf-8"
    )
    assert "evaluation_relation" in labelled
    assert "direct_adjacent" in labelled
    assert (source / "locked_test_edge_results.csv").exists()


def test_locked_audit_refuses_nonempty_output_directory(tmp_path: Path) -> None:
    source = tmp_path / "locked"
    destination = tmp_path / "audit"
    destination.mkdir()
    (destination / "existing.txt").write_text("keep", encoding="utf-8")
    with pytest.raises(audit_module.LockedAuditError, match="new or empty"):
        audit_module.audit_locked_results(source, destination)
