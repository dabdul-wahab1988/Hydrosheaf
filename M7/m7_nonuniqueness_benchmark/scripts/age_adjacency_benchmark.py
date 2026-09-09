"""Evaluation-only labels and metrics for age/topology benchmarks.

The production inference path deliberately does not import this module.  The
functions here consume an already generated candidate edge table and an
independent truth graph supplied by a benchmark.  They are therefore safe to
use after inference to distinguish three estimands that are often conflated:

* direct adjacency (an edge in the generating graph),
* directed reachability (a path of one or more generating edges), and
* temporal direction (forward versus reverse pairs).

In particular, a pair such as ``A -> C`` in a truth graph containing
``A -> B -> C`` is labelled ``transitive_reachable`` rather than
``direct_adjacent``.  This is the negative control that exposes a node-age
ordering signal which is useful for direction but cannot, by itself, identify
the graph's cover relation.

No feature or score from a candidate row is used to create a truth label.
Truth edges and optional path metadata must be passed explicitly by the
caller.  Missing path metadata is reported as an ambiguity for the
cross-path/unrelated distinction rather than being silently treated as a
known relation.
"""

from __future__ import annotations

from collections import Counter, deque
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Any, Literal

import numpy as np

RelationLabel = Literal[
    "direct_adjacent",
    "transitive_reachable",
    "reverse_incompatible",
    "cross_path",
    "unrelated",
    "unknown",
    "invalid",
]

KNOWN_RELATIONS = frozenset(
    {
        "direct_adjacent",
        "transitive_reachable",
        "reverse_incompatible",
        "cross_path",
        "unrelated",
    }
)
FORWARD_RELATIONS = frozenset({"direct_adjacent", "transitive_reachable"})
SKIP_RELATION = "transitive_reachable"
UNKNOWN_RELATIONS = frozenset({"unknown", "invalid", ""})


@dataclass(frozen=True)
class RelationEvidence:
    """Truth-derived relation information for one candidate pair.

    ``status`` is ``"known"`` for a fully interpretable label,
    ``"ambiguous"`` when the graph establishes non-reachability but missing
    path metadata prevents distinguishing a cross-path pair, and
    ``"unknown"`` when the supplied truth metadata cannot support a label.
    """

    label: RelationLabel
    truth_path_length: int | None = None
    reverse_path_length: int | None = None
    same_truth_path: bool | None = None
    status: Literal["known", "ambiguous", "unknown"] = "unknown"
    reason: str = ""


def _coerce_node(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _pair_from_edge(edge: Any) -> tuple[str, str] | None:
    """Read an explicit truth edge from common tuple/mapping/string forms."""

    if isinstance(edge, Mapping):
        upstream = edge.get("u", edge.get("source", edge.get("from")))
        downstream = edge.get("v", edge.get("target", edge.get("to")))
    elif isinstance(edge, str):
        if "->" not in edge:
            return None
        upstream, downstream = edge.split("->", 1)
    else:
        try:
            values = list(edge)
        except TypeError:
            return None
        if len(values) != 2:
            return None
        upstream, downstream = values

    upstream = _coerce_node(upstream)
    downstream = _coerce_node(downstream)
    if upstream is None or downstream is None:
        return None
    return upstream, downstream


def _normalise_truth_edges(
    true_edges: Iterable[Any],
) -> tuple[frozenset[tuple[str, str]], dict[str, set[str]], set[str]]:
    edge_set: set[tuple[str, str]] = set()
    adjacency: dict[str, set[str]] = {}
    nodes: set[str] = set()
    for edge in true_edges:
        pair = _pair_from_edge(edge)
        if pair is None:
            continue
        upstream, downstream = pair
        if upstream == downstream:
            # Self loops cannot identify a forward/reverse relation and are
            # excluded from the truth graph used for path searches.
            nodes.update((upstream, downstream))
            continue
        edge_set.add(pair)
        adjacency.setdefault(upstream, set()).add(downstream)
        nodes.update((upstream, downstream))
    return frozenset(edge_set), adjacency, nodes


def _shortest_path_length(
    adjacency: Mapping[str, set[str]],
    source: str,
    target: str,
) -> int | None:
    if source == target:
        return 0
    queue: deque[tuple[str, int]] = deque([(source, 0)])
    visited = {source}
    while queue:
        node, distance = queue.popleft()
        for child in adjacency.get(node, ()):
            if child in visited:
                continue
            next_distance = distance + 1
            if child == target:
                return next_distance
            visited.add(child)
            queue.append((child, next_distance))
    return None


def _transitive_closure(
    adjacency: Mapping[str, set[str]],
    nodes: Iterable[str],
) -> set[tuple[str, str]]:
    """Return all ordered node pairs connected by a non-empty truth path."""

    reachable: set[tuple[str, str]] = set()
    for source in nodes:
        queue: deque[str] = deque(adjacency.get(source, ()))
        visited: set[str] = set()
        while queue:
            target = queue.popleft()
            if target in visited:
                continue
            visited.add(target)
            if target != source:
                reachable.add((source, target))
            queue.extend(adjacency.get(target, ()))
    return reachable


def _normalise_path_metadata(
    path_metadata: Mapping[Any, Any] | Iterable[Mapping[str, Any]] | None,
) -> dict[str, dict[str, Any]]:
    """Normalise pathline metadata indexed by ``node_id``.

    The independent MODFLOW generator returns a sequence of dictionaries,
    while callers may find a mapping indexed by node ID more convenient.  A
    mapping value can also be a compact path identifier (for example
    ``{"A": 0}``); it is interpreted as ``particle`` metadata.
    """

    if path_metadata is None:
        return {}

    normalised: dict[str, dict[str, Any]] = {}
    if isinstance(path_metadata, Mapping):
        entries = path_metadata.items()
        for key, value in entries:
            node_id = _coerce_node(key)
            if node_id is None:
                continue
            if isinstance(value, Mapping):
                metadata = dict(value)
                metadata.setdefault("node_id", node_id)
            else:
                metadata = {"node_id": node_id, "particle": value}
            normalised[node_id] = metadata
        return normalised

    for value in path_metadata:
        if not isinstance(value, Mapping):
            continue
        node_id = _coerce_node(
            value.get("node_id", value.get("site_id", value.get("sample_id")))
        )
        if node_id is not None:
            normalised[node_id] = dict(value)
    return normalised


def _path_identifier(metadata: Mapping[str, Any]) -> Any:
    for key in (
        "path_id",
        "trajectory_id",
        "particle",
        "particle_id",
        "particle_index",
        "path",
    ):
        if key in metadata and metadata[key] is not None:
            return metadata[key]
    return None


def classify_candidate_pair(
    upstream: Any,
    downstream: Any,
    true_edges: Iterable[Any],
    *,
    path_metadata: Mapping[Any, Any] | Iterable[Mapping[str, Any]] | None = None,
) -> RelationEvidence:
    """Classify a candidate pair using benchmark truth only.

    Parameters
    ----------
    upstream, downstream:
        Candidate node identifiers.
    true_edges:
        Generating directed edges.  Accepted forms are ``(u, v)`` pairs,
        ``{"u": ..., "v": ...}`` mappings, or ``"u->v"`` strings.
    path_metadata:
        Optional node-to-path metadata or a sequence of generator pathline
        rows.  It is used only to identify confirmed cross-path pairs.

    Notes
    -----
    A pair that is neither reachable nor reverse-reachable in a supplied
    graph is not automatically called cross-path.  Without two path IDs this
    remains an ``unrelated`` label with ``status="ambiguous"``.  This keeps
    the useful non-reachability negative while making the missing
    cross-path-vs-unrelated distinction visible to downstream metrics.
    """

    source = _coerce_node(upstream)
    target = _coerce_node(downstream)
    if source is None or target is None or source == target:
        return RelationEvidence(
            label="invalid",
            status="unknown",
            reason="candidate endpoints are missing or identical",
        )

    edge_set, adjacency, truth_nodes = _normalise_truth_edges(true_edges)
    if (source, target) in edge_set:
        return RelationEvidence(
            label="direct_adjacent",
            truth_path_length=1,
            same_truth_path=True,
            status="known",
            reason="candidate pair is an explicit generating edge",
        )

    forward_length = _shortest_path_length(adjacency, source, target)
    if forward_length is not None:
        return RelationEvidence(
            label="transitive_reachable",
            truth_path_length=forward_length,
            same_truth_path=True,
            status="known",
            reason="candidate pair is reachable through two or more generating edges",
        )

    reverse_length = _shortest_path_length(adjacency, target, source)
    if reverse_length is not None:
        return RelationEvidence(
            label="reverse_incompatible",
            reverse_path_length=reverse_length,
            same_truth_path=True,
            status="known",
            reason="generating graph supports the opposite direction",
        )

    metadata = _normalise_path_metadata(path_metadata)
    source_metadata = metadata.get(source)
    target_metadata = metadata.get(target)
    source_path = _path_identifier(source_metadata or {})
    target_path = _path_identifier(target_metadata or {})
    if source_path is not None and target_path is not None:
        same_path = source_path == target_path
        if not same_path:
            return RelationEvidence(
                label="cross_path",
                same_truth_path=False,
                status="known",
                reason="path metadata identifies different generating trajectories",
            )
        return RelationEvidence(
            label="unrelated",
            same_truth_path=True,
            status="known",
            reason="same generating trajectory but no directed path in supplied graph",
        )

    # If both endpoints are explicit truth nodes, non-reachability is useful
    # as a negative for the reachability estimand.  It is nevertheless
    # ambiguous for the finer cross-path/unrelated distinction.
    if source in truth_nodes and target in truth_nodes:
        return RelationEvidence(
            label="unrelated",
            same_truth_path=None,
            status="ambiguous",
            reason="truth graph establishes non-reachability but path IDs are absent",
        )

    return RelationEvidence(
        label="unknown",
        status="unknown",
        reason="one or both candidate endpoints are absent from supplied truth metadata",
    )


def classify_candidate_relation(
    upstream: Any,
    downstream: Any,
    true_edges: Iterable[Any],
    *,
    path_metadata: Mapping[Any, Any] | Iterable[Mapping[str, Any]] | None = None,
) -> RelationLabel:
    """Return only the truth-derived relation label.

    Use :func:`classify_candidate_pair` when status, path length, or reasons
    are needed.  This short wrapper is convenient in table-oriented code.
    """

    return classify_candidate_pair(
        upstream,
        downstream,
        true_edges,
        path_metadata=path_metadata,
    ).label


def _candidate_pair(row: Mapping[str, Any]) -> tuple[str | None, str | None]:
    upstream = row.get("u", row.get("source", row.get("from")))
    downstream = row.get("v", row.get("target", row.get("to")))
    if upstream is None or downstream is None:
        edge_id = row.get("edge_id")
        if isinstance(edge_id, str) and "->" in edge_id:
            upstream, downstream = edge_id.split("->", 1)
    return _coerce_node(upstream), _coerce_node(downstream)


def label_candidate_relations(
    rows: Iterable[Mapping[str, Any]],
    true_edges: Iterable[Any],
    *,
    path_metadata: Mapping[Any, Any] | Iterable[Mapping[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Attach evaluation-only relation columns to candidate feature rows.

    The input rows are copied.  The derived columns use an explicit
    ``evaluation_`` prefix so that they cannot be mistaken for features sent
    to an inference model.  The original feature values, including any score,
    are preserved unchanged.
    """

    # Materialise once because callers commonly pass a generator and truth
    # graph normalisation must be identical for every row.
    truth_edges = tuple(true_edges)
    # Materialise non-mapping metadata once.  In particular, generator
    # pathline rows must not be consumed by the first candidate and silently
    # disappear for all subsequent labels.
    if path_metadata is not None and not isinstance(path_metadata, Mapping):
        path_metadata = tuple(path_metadata)
    labelled: list[dict[str, Any]] = []
    for original in rows:
        row = dict(original)
        upstream, downstream = _candidate_pair(row)
        evidence = classify_candidate_pair(
            upstream,
            downstream,
            truth_edges,
            path_metadata=path_metadata,
        )
        row.update(
            {
                "evaluation_relation": evidence.label,
                "evaluation_relation_status": evidence.status,
                "evaluation_relation_reason": evidence.reason,
                "evaluation_truth_path_length": evidence.truth_path_length,
                "evaluation_reverse_path_length": evidence.reverse_path_length,
                "evaluation_same_truth_path": evidence.same_truth_path,
                "evaluation_direct_target": int(
                    evidence.label == "direct_adjacent"
                ),
                "evaluation_reachable_target": int(
                    evidence.label in FORWARD_RELATIONS
                ),
            }
        )
        labelled.append(row)
    return labelled


def _safe_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _pair_key(row: Mapping[str, Any]) -> tuple[str, str] | None:
    upstream, downstream = _candidate_pair(row)
    if upstream is None or downstream is None:
        return None
    return upstream, downstream


def _score_metrics(
    rows: Sequence[Mapping[str, Any]],
    *,
    positive_relations: frozenset[str],
    negative_relations: frozenset[str],
    probability_key: str | None,
    prefix: str,
) -> dict[str, float | int | None]:
    selected: list[tuple[int, float]] = []
    if probability_key is None:
        return {
            f"{prefix}_n": 0,
            f"{prefix}_roc_auc": None,
            f"{prefix}_pr_auc": None,
            f"{prefix}_brier": None,
            f"{prefix}_calibration_ece": None,
            f"{prefix}_positive_rate": None,
        }
    for row in rows:
        relation = str(row.get("evaluation_relation", "unknown"))
        if relation not in positive_relations | negative_relations:
            continue
        probability = _safe_float(row.get(probability_key))
        if probability is None or not 0.0 <= probability <= 1.0:
            continue
        selected.append((int(relation in positive_relations), probability))

    y = np.asarray([item[0] for item in selected], dtype=int)
    p = np.asarray([item[1] for item in selected], dtype=float)
    n = int(y.size)
    result: dict[str, float | int | None] = {
        f"{prefix}_n": n,
        f"{prefix}_roc_auc": None,
        f"{prefix}_pr_auc": None,
        f"{prefix}_brier": None,
        f"{prefix}_calibration_ece": None,
        f"{prefix}_positive_rate": None,
    }
    if n == 0:
        return result

    result[f"{prefix}_brier"] = float(np.mean((p - y) ** 2))
    result[f"{prefix}_positive_rate"] = float(np.mean(y))
    if np.unique(y).size < 2:
        # ROC-AUC and ranking AP are not identifiable with one class.  Brier
        # and calibration remain defined and are retained above.
        result[f"{prefix}_calibration_ece"] = _expected_calibration_error(y, p)
        return result

    # scikit-learn is a development dependency of the benchmark.  Keep the
    # import local so production HydroSheaf imports do not require it.
    from sklearn.metrics import average_precision_score, roc_auc_score

    result[f"{prefix}_roc_auc"] = float(roc_auc_score(y, p))
    result[f"{prefix}_pr_auc"] = float(average_precision_score(y, p))
    result[f"{prefix}_calibration_ece"] = _expected_calibration_error(y, p)
    return result


def _expected_calibration_error(
    truth: np.ndarray,
    probability: np.ndarray,
    *,
    n_bins: int = 10,
) -> float:
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    # Include p=1.0 in the last bin.  Values were range-checked before this
    # helper is called.
    assignments = np.minimum(np.digitize(probability, edges[1:-1]), n_bins - 1)
    ece = 0.0
    for index in range(n_bins):
        mask = assignments == index
        if not np.any(mask):
            continue
        ece += float(np.mean(mask)) * abs(
            float(np.mean(probability[mask])) - float(np.mean(truth[mask]))
        )
    return float(ece)


def evaluate_age_adjacency(
    rows: Iterable[Mapping[str, Any]],
    *,
    direct_probability_key: str = "direct_probability",
    reachability_probability_key: str | None = "reachability_probability",
    direction_probability_key: str = "direction_probability",
    threshold: float = 0.5,
    true_edges: Iterable[Any] | None = None,
) -> dict[str, float | int | None]:
    """Evaluate direction, reachability, and direct-adjacency separately.

    ``direct_probability_key`` must contain a probability that a candidate is
    a *direct* edge.  ``reachability_probability_key`` is a separate optional
    probability that the candidate is reachable by one or more generating
    edges.  It defaults to ``"reachability_probability"``; when that column is
    absent, reachability metrics are explicitly undefined rather than treating
    a direct-edge score as a reachability score.  A multi-step skip is a
    positive for reachability but a negative for direct adjacency.
    ``direction_probability_key`` should be the probability that the
    candidate's listed direction is forward; direction is evaluated only on
    truth-supported forward/reverse pairs.

    Unknown/invalid labels and missing/non-probability scores are excluded
    from the corresponding score metric and counted explicitly.  If
    ``true_edges`` is supplied, direct candidate containment and recall are
    computed against that complete explicit edge list.
    """

    materialised = [dict(row) for row in rows]
    relation_values = [str(row.get("evaluation_relation", "unknown")) for row in materialised]
    counts = Counter(relation_values)
    known_rows = [
        row
        for row in materialised
        if row.get("evaluation_relation") in KNOWN_RELATIONS
        and row.get("evaluation_relation_status", "known") == "known"
    ]

    result: dict[str, float | int | None] = {
        "n_rows": len(materialised),
        "n_labelled": len(known_rows),
        "n_ambiguous": int(
            sum(
                row.get("evaluation_relation_status") == "ambiguous"
                for row in materialised
            )
        ),
        "n_unknown": int(
            sum(
                row.get("evaluation_relation_status", "unknown") == "unknown"
                or row.get("evaluation_relation") in UNKNOWN_RELATIONS
                for row in materialised
            )
        ),
        "threshold": float(threshold),
        "candidate_direct_count": int(
            sum(row.get("evaluation_relation") == "direct_adjacent" for row in materialised)
        ),
        "candidate_direct_recall": None,
        "direction_n": 0,
        "direction_consistency": None,
        "direction_brier": None,
        "false_skip_n": int(
            sum(row.get("evaluation_relation") == SKIP_RELATION for row in materialised)
        ),
        "false_skip_rejection_rate": None,
        "false_skip_mean_direct_probability": None,
    }
    for relation in (
        "direct_adjacent",
        "transitive_reachable",
        "reverse_incompatible",
        "cross_path",
        "unrelated",
        "unknown",
        "invalid",
    ):
        result[f"n_relation_{relation}"] = int(counts.get(relation, 0))

    result.update(
        _score_metrics(
            materialised,
            positive_relations=frozenset({"direct_adjacent"}),
            negative_relations=frozenset(
                {
                    "transitive_reachable",
                    "reverse_incompatible",
                    "cross_path",
                    "unrelated",
                }
            ),
            probability_key=direct_probability_key,
            prefix="direct_adjacency",
        )
    )
    result.update(
        _score_metrics(
            materialised,
            positive_relations=FORWARD_RELATIONS,
            negative_relations=frozenset(
                {"reverse_incompatible", "cross_path", "unrelated"}
            ),
            probability_key=reachability_probability_key,
            prefix="reachability",
        )
    )

    direction_truth: list[int] = []
    direction_probability: list[float] = []
    for row in materialised:
        relation = str(row.get("evaluation_relation", "unknown"))
        if relation in FORWARD_RELATIONS:
            expected = 1
        elif relation == "reverse_incompatible":
            expected = 0
        else:
            continue
        probability = _safe_float(row.get(direction_probability_key))
        if probability is None or not 0.0 <= probability <= 1.0:
            continue
        direction_truth.append(expected)
        direction_probability.append(probability)
    if direction_truth:
        direction_y = np.asarray(direction_truth, dtype=int)
        direction_p = np.asarray(direction_probability, dtype=float)
        result["direction_n"] = int(direction_y.size)
        result["direction_consistency"] = float(
            np.mean((direction_p >= float(threshold)) == direction_y)
        )
        result["direction_brier"] = float(
            np.mean((direction_p - direction_y) ** 2)
        )

    skip_probabilities = [
        probability
        for row in materialised
        if row.get("evaluation_relation") == SKIP_RELATION
        for probability in [_safe_float(row.get(direct_probability_key))]
        if probability is not None and 0.0 <= probability <= 1.0
    ]
    if skip_probabilities:
        result["false_skip_rejection_rate"] = float(
            np.mean(np.asarray(skip_probabilities) < float(threshold))
        )
        result["false_skip_mean_direct_probability"] = float(
            np.mean(skip_probabilities)
        )

    if true_edges is not None:
        truth_edges_materialised = tuple(true_edges)
        truth_pairs, truth_adjacency, truth_nodes = _normalise_truth_edges(
            truth_edges_materialised
        )
        reachable_pairs = _transitive_closure(truth_adjacency, truth_nodes)
        candidate_pairs = {
            pair
            for row in materialised
            if row.get("evaluation_relation") == "direct_adjacent"
            and (pair := _pair_key(row)) is not None
        }
        contained = len(truth_pairs & candidate_pairs)
        result["truth_direct_count"] = int(len(truth_pairs))
        result["candidate_direct_count"] = int(len(candidate_pairs))
        result["candidate_direct_contained_count"] = int(contained)
        result["candidate_direct_recall"] = (
            float(contained / len(truth_pairs)) if truth_pairs else None
        )
        candidate_reachable_pairs = {
            pair
            for row in materialised
            if row.get("evaluation_relation") in FORWARD_RELATIONS
            and (pair := _pair_key(row)) is not None
        }
        reachable_contained = len(reachable_pairs & candidate_reachable_pairs)
        result["truth_reachable_count"] = int(len(reachable_pairs))
        result["candidate_reachable_count"] = int(len(candidate_reachable_pairs))
        result["candidate_reachable_contained_count"] = int(reachable_contained)
        result["candidate_reachable_recall"] = (
            float(reachable_contained / len(reachable_pairs))
            if reachable_pairs
            else None
        )
    return result


# Explicit aliases make the intended evaluation boundary discoverable to
# benchmark scripts without duplicating implementation.
label_evaluation_relations = label_candidate_relations
evaluate_relation_metrics = evaluate_age_adjacency


__all__ = [
    "FORWARD_RELATIONS",
    "KNOWN_RELATIONS",
    "RelationEvidence",
    "classify_candidate_pair",
    "classify_candidate_relation",
    "evaluate_age_adjacency",
    "evaluate_relation_metrics",
    "label_candidate_relations",
    "label_evaluation_relations",
]
