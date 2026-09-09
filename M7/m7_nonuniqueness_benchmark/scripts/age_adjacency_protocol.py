"""Small, deterministic audit for age-adjacency benchmark result tables.

This module is deliberately an evaluation-side utility.  It does not fit a
model, select a threshold, infer an edge, or read truth-bearing columns as
features.  It checks the schema and computes descriptive metrics for a result
table whose predictions have already been produced.  Keeping this surface
dependency-free makes it useful before the full benchmark environment is
available and reduces the risk that an audit silently changes the locked
analysis.

Example
-------
    python age_adjacency_protocol.py \
        --input results/age_adjacency/edge_results.csv \
        --score-column probability_age_adjacency \
        --truth-count 103 \
        --output results/age_adjacency/audit.json

The input table may contain truth columns for *final evaluation*, but the
selected score column must be prediction-only.  ``is_true_edge`` is consumed
only as the evaluation label.  A direct-adjacency relation column, when
present, is summarized but never used to calculate the score.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


PROTOCOL_NAME = "age-adjacency-audit-v1"
DEFAULT_LABEL_COLUMN = "is_true_edge"
DEFAULT_RELATION_COLUMN = "relation_label"
DEFAULT_THRESHOLD = 0.5
ECE_BINS = 10

# These names are safe as evaluation labels/annotations, but must not be used
# as a prediction score.  The check is intentionally conservative: a new
# truth-bearing feature should require an explicit code review rather than be
# silently accepted by this audit.
TRUTH_ONLY_TOKENS = (
    "is_true",
    "true_edge",
    "true_process",
    "truth",
    "direct_adjacent",
    "transitive_reachable",
    "relation_label",
)


class ProtocolAuditError(ValueError):
    """Raised for malformed benchmark input or a truth-leaking score field."""


def _normalise_bool(value: Any, *, field: str, row_number: int) -> int:
    """Parse the small set of boolean spellings allowed in CSV results."""

    if isinstance(value, bool):
        return int(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return 1
    if text in {"0", "false", "no", "n"}:
        return 0
    raise ProtocolAuditError(
        f"row {row_number}: {field!r} must be boolean/0/1, got {value!r}"
    )


def _finite_float(value: Any, *, field: str, row_number: int) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ProtocolAuditError(
            f"row {row_number}: {field!r} must be numeric, got {value!r}"
        ) from exc
    if not math.isfinite(number):
        raise ProtocolAuditError(
            f"row {row_number}: {field!r} must be finite, got {value!r}"
        )
    return number


def _check_prediction_column(score_column: str) -> None:
    """Reject an obviously truth-bearing prediction column.

    This is a guardrail, not a proof of causal independence.  A score can
    still leak truth under an innocuous name, so the protocol requires code
    review and truth-blind generation in addition to this check.
    """

    lowered = score_column.strip().lower()
    if not lowered:
        raise ProtocolAuditError("score column must not be empty")
    if any(token in lowered for token in TRUTH_ONLY_TOKENS):
        raise ProtocolAuditError(
            f"score column {score_column!r} is truth-bearing by name; "
            "select a prediction-only column"
        )


def _require_columns(
    fieldnames: Sequence[str] | None,
    *,
    score_column: str,
    label_column: str,
    relation_column: str | None,
) -> list[str]:
    names = list(fieldnames or [])
    missing = {
        "edge_id",
        "u",
        "v",
        "seed",
        "split",
        score_column,
        label_column,
    } - set(names)
    if missing:
        raise ProtocolAuditError(
            "missing required columns: " + ", ".join(sorted(missing))
        )
    if relation_column and relation_column not in names:
        # Relation labels are recommended for the adjacency estimand but are
        # not required for auditing older direction-gate result tables.
        return names
    return names


def _read_csv(path: str | Path) -> tuple[list[dict[str, str]], list[str]]:
    input_path = Path(path)
    with input_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ProtocolAuditError(f"{input_path}: CSV has no header")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ProtocolAuditError(f"{input_path}: CSV contains no data rows")
    return rows, list(reader.fieldnames)


def _safe_case_key(row: Mapping[str, Any]) -> str:
    """Return a stable case key without assuming a particular seed type."""

    return f"{row.get('split', '')}:{row.get('seed', '')}"


def _ranking_metrics(labels: Sequence[int], scores: Sequence[float]) -> dict[str, float | None]:
    """Compute ROC-AUC and average precision without scikit-learn.

    Ties are handled as ties for ROC-AUC and by stable score order for average
    precision.  The latter is the standard stepwise precision-recall area for
    a ranked list and is deterministic for tied scores after input order.
    """

    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return {"roc_auc": None, "pr_auc": None}

    # Mann-Whitney formulation with tie averaging.
    indexed = sorted(enumerate(scores), key=lambda item: item[1])
    rank_sum = 0.0
    position = 0
    while position < len(indexed):
        end = position + 1
        score = indexed[position][1]
        while end < len(indexed) and indexed[end][1] == score:
            end += 1
        average_rank = (position + 1 + end) / 2.0
        rank_sum += sum(
            average_rank for index, _ in indexed[position:end] if labels[index] == 1
        )
        position = end
    roc_auc = (rank_sum - positives * (positives + 1) / 2.0) / (
        positives * negatives
    )

    order = sorted(range(len(scores)), key=lambda index: (-scores[index], index))
    seen_positive = 0
    pr_area = 0.0
    for rank, index in enumerate(order, start=1):
        if labels[index]:
            seen_positive += 1
            pr_area += seen_positive / rank
    pr_auc = pr_area / positives
    return {"roc_auc": float(roc_auc), "pr_auc": float(pr_auc)}


def _classification_metrics(
    labels: Sequence[int], scores: Sequence[float], threshold: float
) -> dict[str, float | int | None]:
    predictions = [int(score >= threshold) for score in scores]
    tp = sum(pred == 1 and label == 1 for pred, label in zip(predictions, labels))
    fp = sum(pred == 1 and label == 0 for pred, label in zip(predictions, labels))
    fn = sum(pred == 0 and label == 1 for pred, label in zip(predictions, labels))
    tn = sum(pred == 0 and label == 0 for pred, label in zip(predictions, labels))
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    f1 = (
        2.0 * precision * recall / (precision + recall)
        if precision is not None and recall is not None and precision + recall
        else None
    )
    return {
        "threshold": float(threshold),
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "tn": int(tn),
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def _calibration_metrics(labels: Sequence[int], scores: Sequence[float]) -> dict[str, float]:
    brier = statistics.fmean((score - label) ** 2 for label, score in zip(labels, scores))
    log_loss = statistics.fmean(
        -(
            label * math.log(min(max(score, 1e-15), 1.0))
            + (1 - label) * math.log(min(max(1.0 - score, 1e-15), 1.0))
        )
        for label, score in zip(labels, scores)
    )

    bin_labels: list[list[int]] = [[] for _ in range(ECE_BINS)]
    bin_scores: list[list[float]] = [[] for _ in range(ECE_BINS)]
    for label, score in zip(labels, scores):
        index = min(int(score * ECE_BINS), ECE_BINS - 1)
        bin_labels[index].append(label)
        bin_scores[index].append(score)
    ece = 0.0
    n = len(labels)
    for observed, predicted in zip(bin_labels, bin_scores):
        if observed:
            ece += len(observed) / n * abs(statistics.fmean(observed) - statistics.fmean(predicted))
    return {"brier": float(brier), "log_loss": float(log_loss), "ece_10": float(ece)}


def _relation_counts(rows: Iterable[Mapping[str, Any]], relation_column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        relation = str(row.get(relation_column, "")).strip() or "missing"
        counts[relation] = counts.get(relation, 0) + 1
    return dict(sorted(counts.items()))


def audit_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    fieldnames: Sequence[str] | None = None,
    score_column: str,
    label_column: str = DEFAULT_LABEL_COLUMN,
    relation_column: str | None = DEFAULT_RELATION_COLUMN,
    threshold: float = DEFAULT_THRESHOLD,
    all_truth_count: int | None = None,
) -> dict[str, Any]:
    """Validate and summarize already-generated result rows.

    Parameters
    ----------
    all_truth_count:
        Number of true edges in the complete case graph, including edges not
        present in the candidate table.  If omitted, ``candidate_recall`` is
        left null because candidate-contained truth is not the all-truth
        denominator.
    """

    _check_prediction_column(score_column)
    resolved_fieldnames = list(fieldnames or (rows[0].keys() if rows else []))
    _require_columns(
        resolved_fieldnames,
        score_column=score_column,
        label_column=label_column,
        relation_column=relation_column,
    )
    if not rows:
        raise ProtocolAuditError("result table contains no rows")
    if not 0.0 <= threshold <= 1.0:
        raise ProtocolAuditError("threshold must be in [0, 1]")
    if all_truth_count is not None and all_truth_count < 0:
        raise ProtocolAuditError("all_truth_count must be non-negative")

    labels: list[int] = []
    scores: list[float] = []
    case_keys: set[str] = set()
    seen_edges: set[tuple[str, str, str, str]] = set()
    for row_number, row in enumerate(rows, start=2):
        for name in ("edge_id", "u", "v", "seed", "split"):
            if not str(row.get(name, "")).strip():
                raise ProtocolAuditError(f"row {row_number}: {name!r} must not be empty")
        edge_key = (
            str(row["split"]),
            str(row["seed"]),
            str(row["u"]),
            str(row["v"]),
        )
        if edge_key in seen_edges:
            raise ProtocolAuditError(f"row {row_number}: duplicate candidate edge {edge_key!r}")
        seen_edges.add(edge_key)
        labels.append(_normalise_bool(row[label_column], field=label_column, row_number=row_number))
        score = _finite_float(row[score_column], field=score_column, row_number=row_number)
        if not 0.0 <= score <= 1.0:
            raise ProtocolAuditError(
                f"row {row_number}: {score_column!r} must be in [0, 1], got {score}"
            )
        scores.append(score)
        case_keys.add(_safe_case_key(row))

    candidate_truth_count = sum(labels)
    if all_truth_count is not None and candidate_truth_count > all_truth_count:
        raise ProtocolAuditError(
            "candidate truth count exceeds all-truth count; check the denominator"
        )
    candidate_recall = (
        candidate_truth_count / all_truth_count
        if all_truth_count not in (None, 0)
        else (1.0 if all_truth_count == 0 and candidate_truth_count == 0 else None)
    )
    summary: dict[str, Any] = {
        "n_rows": len(rows),
        "n_cases": len(case_keys),
        "n_positive_candidates": candidate_truth_count,
        "n_negative_candidates": len(rows) - candidate_truth_count,
        "candidate_prevalence": candidate_truth_count / len(rows),
        "all_truth_count": all_truth_count,
        "candidate_recall_against_all_truth": candidate_recall,
        "splits": {},
    }
    # Construct deterministic split counts explicitly for the JSON contract.
    split_counts: dict[str, int] = {}
    for row in rows:
        split = str(row["split"])
        split_counts[split] = split_counts.get(split, 0) + 1
    summary["splits"] = dict(sorted(split_counts.items()))
    summary["score_column"] = score_column
    summary["score_min"] = min(scores)
    summary["score_max"] = max(scores)
    summary["ranking"] = _ranking_metrics(labels, scores)
    summary["calibration"] = _calibration_metrics(labels, scores)
    summary["classification"] = _classification_metrics(labels, scores, threshold)
    if relation_column and relation_column in resolved_fieldnames:
        summary["relation_column"] = relation_column
        summary["relation_counts"] = _relation_counts(rows, relation_column)
    else:
        summary["relation_column"] = None
        summary["relation_counts"] = None
    return summary


def audit_csv(
    path: str | Path,
    *,
    score_column: str,
    label_column: str = DEFAULT_LABEL_COLUMN,
    relation_column: str | None = DEFAULT_RELATION_COLUMN,
    threshold: float = DEFAULT_THRESHOLD,
    all_truth_count: int | None = None,
) -> dict[str, Any]:
    """Read and audit a CSV result table."""

    rows, fieldnames = _read_csv(path)
    summary = audit_rows(
        rows,
        fieldnames=fieldnames,
        score_column=score_column,
        label_column=label_column,
        relation_column=relation_column,
        threshold=threshold,
        all_truth_count=all_truth_count,
    )
    return {
        "protocol": PROTOCOL_NAME,
        "input": str(Path(path)),
        "valid": True,
        "errors": [],
        "warnings": [
            "This audit is descriptive; uncertainty intervals and threshold "
            "selection must be computed in the preregistered grouped benchmark.",
            "A null all-truth recall means that the complete-graph truth count "
            "was not supplied.",
        ],
        "summary": summary,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="CSV result table")
    parser.add_argument("--score-column", required=True, help="prediction-only probability column")
    parser.add_argument("--label-column", default=DEFAULT_LABEL_COLUMN)
    parser.add_argument("--relation-column", default=DEFAULT_RELATION_COLUMN)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument(
        "--truth-count",
        type=int,
        default=None,
        help="complete-graph true-edge count for candidate recall",
    )
    parser.add_argument("--output", type=Path, default=None, help="optional JSON output path")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = audit_csv(
            args.input,
            score_column=args.score_column,
            label_column=args.label_column,
            relation_column=args.relation_column or None,
            threshold=args.threshold,
            all_truth_count=args.truth_count,
        )
    except (OSError, ProtocolAuditError) as exc:
        result = {
            "protocol": PROTOCOL_NAME,
            "input": str(args.input),
            "valid": False,
            "errors": [str(exc)],
            "warnings": [],
        }
        exit_code = 2
    else:
        exit_code = 0

    rendered = json.dumps(result, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - exercised through CLI smoke test
    raise SystemExit(main())
