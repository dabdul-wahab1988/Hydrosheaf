"""Run the relation-level age audit on an existing locked M7 result.

This command is intentionally read-only with respect to the locked result
directory.  It joins independent case truth *after* inference, labels direct
edges versus transitive skips, and evaluates the historical age-compatibility
score against those separate estimands.  The resulting score is called an
``age_compatibility_probability`` because the direction-gate output is not a
calibrated direct-adjacency probability.  A new output directory is required
for every audit.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from age_adjacency_benchmark import (  # noqa: E402
    evaluate_age_adjacency,
    label_candidate_relations,
)


PROTOCOL_NAME = "age-adjacency-locked-audit-v1"


class LockedAuditError(ValueError):
    """Raised when a locked result/case artifact is incomplete or inconsistent."""


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise LockedAuditError(f"{path}: CSV has no header")
        rows = [dict(row) for row in reader]
    if not rows:
        raise LockedAuditError(f"{path}: CSV contains no rows")
    return rows


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise LockedAuditError("cannot write an empty labelled result table")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(str(key))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows({key: row.get(key) for key in fieldnames} for row in rows)


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _compatibility_probability(value: Any) -> float | None:
    try:
        cost = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(cost):
        return None
    # The locked direction-gate score is a non-negative cost.  Clipping only
    # prevents underflow; it does not turn the score into directness evidence.
    return math.exp(-min(40.0, max(0.0, cost)))


def _truth_artifacts(case_dir: Path) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    truth_path = case_dir / "heldout_truth.csv"
    pathline_path = case_dir / "modpath_pathline_truth.csv"
    if not truth_path.exists() or not pathline_path.exists():
        raise LockedAuditError(
            f"{case_dir}: heldout_truth.csv and modpath_pathline_truth.csv are required"
        )
    truth_rows = _read_csv(truth_path)
    pathline_rows = _read_csv(pathline_path)
    edges: list[tuple[str, str]] = []
    for row in truth_rows:
        upstream = str(row.get("u", "")).strip()
        downstream = str(row.get("v", "")).strip()
        if not upstream or not downstream or upstream == downstream:
            raise LockedAuditError(f"{truth_path}: malformed truth edge {row!r}")
        edges.append((upstream, downstream))
    return edges, pathline_rows


def audit_locked_results(
    results_dir: str | Path,
    output_dir: str | Path,
    *,
    require_split: str = "locked_test",
) -> dict[str, object]:
    """Label and audit a locked result table into a new output directory."""

    source_dir = Path(results_dir).resolve()
    destination = Path(output_dir).resolve()
    if source_dir == destination:
        raise LockedAuditError("output_dir must be separate from the locked results directory")
    if destination.exists() and any(destination.iterdir()):
        raise LockedAuditError(
            f"output_dir must be new or empty to preserve audit provenance: {destination}"
        )
    input_path = source_dir / "locked_test_edge_results.csv"
    if not input_path.exists():
        raise LockedAuditError(f"missing locked result table: {input_path}")
    rows = _read_csv(input_path)
    if any(str(row.get("split", "")).strip() != require_split for row in rows):
        raise LockedAuditError(f"all rows must have split={require_split!r}")
    if "age_cost" not in rows[0]:
        raise LockedAuditError("locked result table must contain age_cost")

    case_root = source_dir / "cases"
    labelled: list[dict[str, Any]] = []
    all_truth_edges: list[tuple[str, str]] = []
    source_case_hashes: dict[str, dict[str, str]] = {}
    seeds = sorted({str(row.get("seed", "")).strip() for row in rows})
    if "" in seeds:
        raise LockedAuditError("every locked row must contain a seed")
    for seed in seeds:
        case_dir = case_root / f"{require_split}_{seed}"
        truth_edges, pathline_rows = _truth_artifacts(case_dir)
        all_truth_edges.extend(truth_edges)
        case_rows = [row for row in rows if str(row.get("seed")) == seed]
        case_labelled = label_candidate_relations(
            case_rows,
            truth_edges,
            path_metadata=pathline_rows,
        )
        for row in case_labelled:
            score = _compatibility_probability(row.get("age_cost"))
            row["age_compatibility_probability"] = score
            row["direction_probability"] = score
            row["prediction_status"] = "scored" if score is not None else "invalid"
            row["score_interpretation"] = "direction_compatibility_not_direct_probability"
        labelled.extend(case_labelled)
        source_case_hashes[seed] = {
            "heldout_truth.csv": _sha256(case_dir / "heldout_truth.csv"),
            "modpath_pathline_truth.csv": _sha256(
                case_dir / "modpath_pathline_truth.csv"
            ),
        }

    metrics = evaluate_age_adjacency(
        labelled,
        direct_probability_key="age_compatibility_probability",
        direction_probability_key="direction_probability",
        true_edges=all_truth_edges,
    )
    metrics["score_interpretation"] = (
        "historical age compatibility score evaluated against direct labels; "
        "not a calibrated direct-adjacency probability"
    )
    manifest_path = source_dir / "manifest.json"
    audit_manifest: dict[str, object] = {
        "protocol": PROTOCOL_NAME,
        "source_results_dir": str(source_dir),
        "source_result_sha256": _sha256(input_path),
        "source_manifest_sha256": (
            _sha256(manifest_path) if manifest_path.exists() else None
        ),
        "case_truth_hashes": source_case_hashes,
        "split": require_split,
        "seeds": seeds,
        "n_rows": len(labelled),
        "n_complete_truth_direct_edges": len(set(all_truth_edges)),
        "truth_join": "post_inference_evaluation_only",
        "locked_results_modified": False,
        "metrics": metrics,
    }
    destination.mkdir(parents=True, exist_ok=True)
    _write_csv(destination / "relation_labeled_locked_test_edge_results.csv", labelled)
    (destination / "relation_metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (destination / "audit_manifest.json").write_text(
        json.dumps(audit_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return audit_manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--split", default="locked_test")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    manifest = audit_locked_results(
        args.results_dir,
        args.output_dir,
        require_split=args.split,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
