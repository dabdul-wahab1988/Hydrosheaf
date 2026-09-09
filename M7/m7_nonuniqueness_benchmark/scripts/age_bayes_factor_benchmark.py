"""Controlled synthetic benchmark for age-supported directness inference.

This benchmark isolates a narrow scientific question: can an observed age
increment distinguish a one-segment (direct) candidate from a multi-segment
(indirect) candidate when the candidate also supplies independent,
edge-specific transport hypotheses?

The generator creates monotone groundwater-like paths whose latent age
increments equal positive segment travel times.  Candidate rows expose only
the noisy endpoint ages, their uncertainties, and paired direct/indirect
transport hypotheses.  The direct/indirect relation is kept in a separate
truth ledger and is joined only by :func:`evaluate_predictions` after scores
have been produced.  Consequently the generated inference table contains no
``is_direct``/``truth`` field.

Four prespecified controls are evaluated on a held-out split:

``no_age``
    A null 0.5 prior with age evidence omitted.  It is marked as abstaining;
    its Brier/log-loss values are retained as a transparent baseline.
``order_only``
    The one-sided forward age-ordering probability.  It tests the known
    failure mode in which both direct and multi-step forward pairs receive
    similar support.
``full_bf``
    The existing HydroSheaf age-adjacency likelihood evaluator, comparing the
    observed increment with both supplied hypotheses.  Candidates whose
    hypotheses are not separated by the declared uncertainty threshold
    abstain rather than being forced into a direct/indirect decision.
``permuted_full_bf``
    The same full evaluator after a deterministic within-case permutation of
    transport hypotheses.  Ages and endpoint uncertainty are unchanged; the
    control tests whether edge-specific transport pairing carries information.

The output is a run-scoped manifest plus separate inference, truth, held-out
prediction, and metric files.  It is controlled-synthetic evidence only and
must not be presented as field validation or as proof that age alone identifies
direct adjacency.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Literal, Mapping, Sequence

import numpy as np

from hydrosheaf.validation.age_adjacency import (
    compute_age_adjacency_evidence,
    compute_direction_evidence,
)

# The relation evaluator is intentionally kept in the M7 script tier because
# it is an evaluation-only module and is not part of the production inference
# path.  The fallback keeps direct imports useful when this file is loaded from
# a test or notebook with a different working directory.
try:  # pragma: no cover - import branch depends on invocation context
    from age_adjacency_benchmark import (
        evaluate_age_adjacency,
        label_candidate_relations,
    )
except ImportError:  # pragma: no cover
    _SCRIPT_DIR = Path(__file__).resolve().parent
    if str(_SCRIPT_DIR) not in sys.path:
        sys.path.insert(0, str(_SCRIPT_DIR))
    from age_adjacency_benchmark import (  # type: ignore[no-redef]
        evaluate_age_adjacency,
        label_candidate_relations,
    )


PROTOCOL_NAME = "age-direct-versus-indirect-bayes-factor-v1"
METHODS = ("no_age", "order_only", "full_bf", "permuted_full_bf")
CASE_STRATA = ("separable_transport", "overlapping_transport", "mixed_transport")
IDENTIFIABILITY_STRATA = ("separable", "intermediate", "overlapping")
RELATIONS = ("direct_adjacent", "indirect_reachable")
DEFAULT_OUTPUT = Path(".codex_work/runs/RUN-AGE-BF-CONTROLLED-20260908-01")


def _clip_probability(value: float) -> float:
    return float(min(1.0 - 1.0e-12, max(1.0e-12, value)))


def _sigmoid(value: float) -> float:
    # Clipping makes JSON/metric output stable for very decisive likelihoods.
    bounded = min(40.0, max(-40.0, float(value)))
    if bounded >= 0.0:
        z = math.exp(-bounded)
        return float(1.0 / (1.0 + z))
    z = math.exp(bounded)
    return float(z / (1.0 + z))


def _finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


@dataclass(frozen=True)
class AgeCandidate:
    """Inference-visible candidate fields; no truth label is stored here."""

    candidate_id: str
    u: str
    v: str
    path_id: str
    hop_count: int
    upstream_age_years: float
    downstream_age_years: float
    upstream_sigma_years: float
    downstream_sigma_years: float
    process_sigma_years: float
    direct_travel_years: float
    indirect_travel_years: float
    direct_travel_sigma_years: float
    indirect_travel_sigma_years: float
    transport_separation_z: float
    identifiability_stratum: str
    # Endpoint-age error covariance is explicit even when the controlled
    # generator sets it to zero.  A non-zero value is required whenever the
    # two endpoint estimates share a laboratory, campaign, or model error.
    age_covariance_years2: float = 0.0

    def to_inference_row(self, *, case_id: str, case_stratum: str) -> dict[str, Any]:
        result = asdict(self)
        result.update(
            {
                "case_id": case_id,
                "case_stratum": case_stratum,
                "reference_type": "independent_synthetic_truth",
            }
        )
        # Defensive assertion: the inference object must never silently grow a
        # truth-bearing column through a future dataclass edit.
        forbidden = {"is_direct", "truth_label", "relation", "direct_target", "indirect_target"}
        if forbidden.intersection(result):
            raise AssertionError("AgeCandidate inference row contains a truth-bearing field")
        return result


@dataclass(frozen=True)
class AgeTruth:
    """Evaluation-only truth ledger kept separate from inference rows."""

    relation_by_candidate: Mapping[str, Literal["direct_adjacent", "indirect_reachable"]]
    direct_edges: tuple[tuple[str, str], ...]
    path_by_node: Mapping[str, str]
    latent_age_by_node: Mapping[str, float]

    def relation(self, candidate_id: str) -> str:
        return str(self.relation_by_candidate[candidate_id])


@dataclass(frozen=True)
class AgeCase:
    case_id: str
    seed: int
    split: Literal["development", "locked_test"]
    case_stratum: str
    candidates: tuple[AgeCandidate, ...]
    truth: AgeTruth

    def inference_rows(self) -> list[dict[str, Any]]:
        return [
            candidate.to_inference_row(
                case_id=self.case_id,
                case_stratum=self.case_stratum,
            )
            for candidate in self.candidates
        ]


def _identifiability(separation_z: float) -> str:
    if separation_z >= 3.0:
        return "separable"
    if separation_z <= 1.0:
        return "overlapping"
    return "intermediate"


def _positive_segment_times(rng: np.random.Generator, n: int) -> np.ndarray:
    # A modest, positive travel-time distribution in years.  It is a synthetic
    # reference scale, not a claim about any particular aquifer.
    return rng.lognormal(mean=math.log(7.0), sigma=0.20, size=n)


def _transport_gap(
    rng: np.random.Generator,
    actual_increment: float,
    case_stratum: str,
    *,
    relation: str,
) -> tuple[float, float, float, float]:
    """Create edge-specific direct/indirect hypotheses without labels in rows."""

    if case_stratum == "separable_transport":
        gap_fraction = float(rng.uniform(0.35, 0.55))
        hypothesis_sigma = 0.30
        endpoint_sigma = 0.28
        process_sigma = 0.18
    elif case_stratum == "overlapping_transport":
        gap_fraction = float(rng.uniform(0.01, 0.06))
        hypothesis_sigma = 1.15
        endpoint_sigma = 1.35
        process_sigma = 0.65
    else:
        gap_fraction = float(rng.uniform(0.12, 0.25))
        hypothesis_sigma = 0.65
        endpoint_sigma = 0.70
        process_sigma = 0.35
    gap = max(0.08, gap_fraction * actual_increment)
    direct_noise = float(rng.normal(0.0, hypothesis_sigma))
    indirect_noise = float(rng.normal(0.0, hypothesis_sigma))
    if relation == "direct_adjacent":
        direct = max(0.02, actual_increment + direct_noise)
        indirect = max(0.02, actual_increment + gap + indirect_noise)
    else:
        direct = max(0.02, actual_increment - gap + direct_noise)
        indirect = max(0.02, actual_increment + indirect_noise)
    return direct, indirect, endpoint_sigma, process_sigma


def generate_age_case(
    seed: int,
    *,
    split: Literal["development", "locked_test"] = "locked_test",
    case_stratum: str | None = None,
    n_paths: int = 2,
    nodes_per_path: int = 9,
) -> AgeCase:
    """Generate one sealed direct/indirect case.

    Direct candidates are adjacent nodes on a generating path.  Indirect
    candidates skip one node (two generating segments).  The candidate's
    ``hop_count`` and transport hypotheses describe the inference problem, but
    the relation label remains only in :class:`AgeTruth`.
    """

    if n_paths < 1:
        raise ValueError("n_paths must be positive")
    if nodes_per_path < 4:
        raise ValueError("nodes_per_path must be at least 4")
    if case_stratum is None:
        case_stratum = CASE_STRATA[int(seed) % len(CASE_STRATA)]
    if case_stratum not in CASE_STRATA:
        raise ValueError(f"unknown case_stratum: {case_stratum}")
    rng = np.random.default_rng(int(seed))
    case_id = f"AGEBF-{int(seed):06d}"
    latent_age_by_node: dict[str, float] = {}
    path_by_node: dict[str, str] = {}
    direct_edges: list[tuple[str, str]] = []
    path_nodes: dict[str, list[str]] = {}
    node_observations: dict[str, tuple[float, float]] = {}

    for path_index in range(n_paths):
        path_id = f"P{path_index:02d}"
        nodes = [f"{path_id}-N{index:02d}" for index in range(nodes_per_path)]
        path_nodes[path_id] = nodes
        age = float(rng.uniform(2.0, 8.0))
        segment_times = _positive_segment_times(rng, nodes_per_path - 1)
        for node_index, node_id in enumerate(nodes):
            if node_index:
                age += float(segment_times[node_index - 1])
            latent_age_by_node[node_id] = age
            path_by_node[node_id] = path_id
            if case_stratum == "separable_transport":
                endpoint_sigma = 0.28
            elif case_stratum == "overlapping_transport":
                endpoint_sigma = 1.35
            else:
                endpoint_sigma = 0.70
            observed = age + float(rng.normal(0.0, endpoint_sigma))
            node_observations[node_id] = (observed, endpoint_sigma)

        for index in range(nodes_per_path - 1):
            direct_edges.append((nodes[index], nodes[index + 1]))

    candidates: list[AgeCandidate] = []
    relation_by_candidate: dict[str, Literal["direct_adjacent", "indirect_reachable"]] = {}
    # Include one-hop direct pairs and two-hop skips.  A fixed order keeps the
    # candidate composition and labels reproducible across Python versions.
    for path_id, nodes in path_nodes.items():
        for hop_count, relation in ((1, "direct_adjacent"), (2, "indirect_reachable")):
            for index in range(nodes_per_path - hop_count):
                u, v = nodes[index], nodes[index + hop_count]
                candidate_id = f"{case_id}:{u}->{v}"
                latent_increment = latent_age_by_node[v] - latent_age_by_node[u]
                direct, indirect, hypothesis_sigma, process_sigma = _transport_gap(
                    rng,
                    latent_increment,
                    case_stratum,
                    relation=relation,
                )
                separation = abs(direct - indirect) / math.sqrt(
                    2.0 * hypothesis_sigma**2
                )
                observed_u, sigma_u = node_observations[u]
                observed_v, sigma_v = node_observations[v]
                candidate = AgeCandidate(
                    candidate_id=candidate_id,
                    u=u,
                    v=v,
                    path_id=path_id,
                    hop_count=hop_count,
                    upstream_age_years=float(observed_u),
                    downstream_age_years=float(observed_v),
                    upstream_sigma_years=float(sigma_u),
                    downstream_sigma_years=float(sigma_v),
                    process_sigma_years=float(process_sigma),
                    direct_travel_years=float(direct),
                    indirect_travel_years=float(indirect),
                    direct_travel_sigma_years=float(hypothesis_sigma),
                    indirect_travel_sigma_years=float(hypothesis_sigma),
                    transport_separation_z=float(separation),
                    identifiability_stratum=_identifiability(separation),
                )
                candidates.append(candidate)
                relation_by_candidate[candidate_id] = relation  # type: ignore[assignment]

    return AgeCase(
        case_id=case_id,
        seed=int(seed),
        split=split,
        case_stratum=case_stratum,
        candidates=tuple(candidates),
        truth=AgeTruth(
            relation_by_candidate=relation_by_candidate,
            direct_edges=tuple(direct_edges),
            path_by_node=path_by_node,
            latent_age_by_node=latent_age_by_node,
        ),
    )


def generate_cases(
    *,
    first_seed: int = 2026090801,
    n_cases: int = 24,
    n_development: int | None = None,
) -> tuple[AgeCase, ...]:
    """Generate a deterministic development/held-out case schedule."""

    if n_cases < 2:
        raise ValueError("n_cases must be at least 2")
    development = n_cases // 2 if n_development is None else int(n_development)
    if not 1 <= development < n_cases:
        raise ValueError("n_development must leave at least one locked-test case")
    cases = []
    for index in range(n_cases):
        split: Literal["development", "locked_test"] = (
            "development" if index < development else "locked_test"
        )
        cases.append(generate_age_case(first_seed + index, split=split))
    return tuple(cases)


def _transport_permutation(
    rows: Sequence[Mapping[str, Any]],
    *,
    seed: int,
) -> list[dict[str, Any]]:
    """Permute only edge-specific transport fields within one case."""

    if not rows:
        return []
    rng = np.random.default_rng(int(seed))
    permutation = rng.permutation(len(rows))
    fields = (
        "direct_travel_years",
        "indirect_travel_years",
        "direct_travel_sigma_years",
        "indirect_travel_sigma_years",
        "transport_separation_z",
        "identifiability_stratum",
    )
    output: list[dict[str, Any]] = []
    for target_index, source_index in enumerate(permutation):
        result = dict(rows[target_index])
        result["transport_source_candidate_id"] = rows[int(source_index)]["candidate_id"]
        for field_name in fields:
            result[field_name] = rows[int(source_index)][field_name]
        output.append(result)
    return output


def _score_candidate(
    row: Mapping[str, Any],
    *,
    method: str,
    min_separation_z: float,
) -> dict[str, Any]:
    """Produce a prediction without accessing a truth ledger."""

    if method not in METHODS:
        raise ValueError(f"unknown method: {method}")
    result = dict(row)
    result.update(
        {
            "method": method,
            # A probability is not emitted for an abstention.  Evaluation
            # metrics may use the declared 0.5 prior as a transparent null,
            # but persisted prediction rows must distinguish that prior from
            # a scored directness probability.
            "direct_probability": None,
            "log_bayes_factor_direct_vs_indirect": None,
            "score_status": "abstain",
            "abstain_reason": None,
            "evidence_class": method,
        }
    )
    if method == "no_age":
        result["abstain_reason"] = "age_observation_omitted"
        return result

    common = {
        "upstream_age_years": row["upstream_age_years"],
        "downstream_age_years": row["downstream_age_years"],
        "upstream_sigma_years": row["upstream_sigma_years"],
        "downstream_sigma_years": row["downstream_sigma_years"],
        "age_covariance_years2": row.get("age_covariance_years2", 0.0),
    }
    if method == "order_only":
        evidence = compute_direction_evidence(**common)
        result.update(
            {
                "direct_probability": evidence.forward_compatibility_probability,
                "direction_probability": evidence.forward_compatibility_probability,
                "score_status": "scored",
                "abstain_reason": None,
                "evidence_class": "temporal_order_only",
            }
        )
        return result

    evidence = compute_age_adjacency_evidence(
        **common,
        process_sigma_years=row["process_sigma_years"],
        direct_travel_years=row["direct_travel_years"],
        indirect_travel_years=row["indirect_travel_years"],
        direct_travel_sigma_years=row["direct_travel_sigma_years"],
        indirect_travel_sigma_years=row["indirect_travel_sigma_years"],
    )
    result["direction_probability"] = evidence.direction.forward_compatibility_probability
    separation = abs(
        float(row["direct_travel_years"]) - float(row["indirect_travel_years"])
    ) / math.sqrt(
        float(row["direct_travel_sigma_years"]) ** 2
        + float(row["indirect_travel_sigma_years"]) ** 2
    )
    if (
        evidence.log_bayes_factor_direct_vs_indirect is None
        or not evidence.comparison_discriminating
    ):
        result["abstain_reason"] = "likelihood_comparison_unavailable"
        return result
    if separation < float(min_separation_z):
        result["abstain_reason"] = "transport_hypotheses_overlap"
        return result
    log_bf = float(evidence.log_bayes_factor_direct_vs_indirect)
    result.update(
        {
            "direct_probability": _sigmoid(log_bf),
            "log_bayes_factor_direct_vs_indirect": log_bf,
            "score_status": "scored",
            "abstain_reason": None,
            "evidence_class": "age_plus_transport_bayes_factor",
            "bf_comparison_discriminating": True,
        }
    )
    return result


def score_case(
    case: AgeCase,
    *,
    method: str,
    min_separation_z: float = 1.0,
) -> list[dict[str, Any]]:
    """Score one case from inference rows only."""

    rows = case.inference_rows()
    if method == "permuted_full_bf":
        rows = _transport_permutation(rows, seed=case.seed + 811_733)
        score_method = "full_bf"
    else:
        score_method = method
    scored = [
        _score_candidate(
            row,
            method=score_method,
            min_separation_z=min_separation_z,
        )
        for row in rows
    ]
    for row in scored:
        row["method"] = method
        if method == "permuted_full_bf":
            row["evidence_class"] = "permuted_age_plus_transport_bayes_factor"
    return scored


def _ece(truth: np.ndarray, probability: np.ndarray, *, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    assignments = np.minimum(np.digitize(probability, edges[1:-1]), bins - 1)
    total = 0.0
    for index in range(bins):
        mask = assignments == index
        if np.any(mask):
            total += float(np.mean(mask)) * abs(
                float(np.mean(probability[mask])) - float(np.mean(truth[mask]))
            )
    return float(total)


def _basic_directness_metrics(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    labels = np.asarray(
        [int(row["truth_label"] == "direct_adjacent") for row in rows],
        dtype=int,
    )
    # For inclusive diagnostic metrics, an abstention is represented by the
    # predeclared equal-prior null (0.5).  The row itself remains null-valued
    # so downstream QA cannot mistake an abstention for a prediction.
    probabilities = np.asarray(
        [
            _clip_probability(
                float(row["direct_probability"])
                if _finite_float(row.get("direct_probability")) is not None
                else 0.5
            )
            for row in rows
        ],
        dtype=float,
    )
    status = [str(row.get("score_status", "abstain")) for row in rows]
    scored_mask = np.asarray([item == "scored" for item in status], dtype=bool)
    result: dict[str, Any] = {
        "n": int(labels.size),
        "n_scored": int(np.sum(scored_mask)),
        "n_abstain": int(np.sum(~scored_mask)),
        "abstain_rate": float(np.mean(~scored_mask)) if labels.size else None,
        "direct_prevalence": float(np.mean(labels)) if labels.size else None,
        "directness_pr_auc": None,
        "directness_roc_auc": None,
        "brier": None,
        "log_loss": None,
        "calibration_ece": None,
        "mean_direct_probability": float(np.mean(probabilities)) if labels.size else None,
        "bf_sign_n": 0,
        "bf_sign_accuracy": None,
        "bf_sign_coverage": 0.0,
    }
    if not labels.size:
        return result
    result["brier"] = float(np.mean((probabilities - labels) ** 2))
    result["log_loss"] = float(
        -np.mean(labels * np.log(probabilities) + (1 - labels) * np.log(1.0 - probabilities))
    )
    result["calibration_ece"] = _ece(labels, probabilities)
    if np.unique(labels).size > 1:
        from sklearn.metrics import average_precision_score, roc_auc_score

        result["directness_pr_auc"] = float(average_precision_score(labels, probabilities))
        result["directness_roc_auc"] = float(roc_auc_score(labels, probabilities))
    signs: list[bool] = []
    for row, label in zip(rows, labels, strict=True):
        log_bf = _finite_float(row.get("log_bayes_factor_direct_vs_indirect"))
        if log_bf is None or str(row.get("score_status")) != "scored":
            continue
        signs.append((log_bf > 0.0) == bool(label))
    result["bf_sign_n"] = len(signs)
    result["bf_sign_accuracy"] = float(np.mean(signs)) if signs else None
    result["bf_sign_coverage"] = float(len(signs) / labels.size)
    return result


def _truth_join(
    case_by_id: Mapping[str, AgeCase],
    prediction_rows: Iterable[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    joined: list[dict[str, Any]] = []
    for original in prediction_rows:
        row = dict(original)
        case = case_by_id.get(str(row.get("case_id")))
        if case is None:
            raise ValueError(f"prediction row references unknown case {row.get('case_id')!r}")
        candidate_id = str(row.get("candidate_id"))
        if candidate_id not in case.truth.relation_by_candidate:
            raise ValueError(f"prediction row references unknown candidate {candidate_id!r}")
        row["truth_label"] = case.truth.relation(candidate_id)
        joined.append(row)
    return joined


def _legacy_evaluator_metrics(
    case_by_id: Mapping[str, AgeCase],
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Run the existing relation evaluator as a compatibility diagnostic."""

    labelled: list[dict[str, Any]] = []
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row["case_id"]), []).append(row)
    for case_id, case_rows in grouped.items():
        case = case_by_id[case_id]
        transformed = []
        for row in case_rows:
            # Node IDs are reused across synthetic cases (P00-N01, etc.).
            # Prefix them before applying the legacy evaluator, otherwise
            # truth edges from different cases are silently merged and the
            # complete-graph denominators are wrong.
            scoped_u = f"{case_id}::{row['u']}"
            scoped_v = f"{case_id}::{row['v']}"
            transformed.append(
                {
                    **row,
                    "u": scoped_u,
                    "v": scoped_v,
                    "edge_id": f"{scoped_u}->{scoped_v}",
                    "direct_probability": row["direct_probability"],
                    "reachability_probability": row["direct_probability"],
                    "direction_probability": row.get("direction_probability"),
                }
            )
        path_metadata = [
            {"node_id": f"{case_id}::{node}", "path_id": path}
            for node, path in case.truth.path_by_node.items()
        ]
        scoped_truth_edges = [
            (f"{case_id}::{u}", f"{case_id}::{v}")
            for u, v in case.truth.direct_edges
        ]
        labelled.extend(
            label_candidate_relations(
                transformed,
                scoped_truth_edges,
                path_metadata=path_metadata,
            )
        )
    scoped_truth = [
        (f"{case_id}::{u}", f"{case_id}::{v}")
        for case_id, case in case_by_id.items()
        for u, v in case.truth.direct_edges
    ]
    return evaluate_age_adjacency(
        labelled,
        direct_probability_key="direct_probability",
        reachability_probability_key="reachability_probability",
        direction_probability_key="direction_probability",
        true_edges=scoped_truth,
    )


def evaluate_predictions(
    cases: Sequence[AgeCase],
    prediction_rows: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    """Join held-out truth after scoring and compute directness metrics."""

    case_by_id = {case.case_id: case for case in cases}
    joined = _truth_join(case_by_id, prediction_rows)
    by_method: dict[str, list[dict[str, Any]]] = {}
    for row in joined:
        by_method.setdefault(str(row["method"]), []).append(row)
    output: dict[str, Any] = {
        "protocol": PROTOCOL_NAME,
        "n_cases": len(cases),
        "n_rows": len(joined),
        "methods": {},
    }
    for method in METHODS:
        method_rows = by_method.get(method, [])
        metrics = _basic_directness_metrics(method_rows)
        strata: dict[str, Any] = {}
        for stratum in IDENTIFIABILITY_STRATA:
            stratum_rows = [
                row for row in method_rows if row.get("identifiability_stratum") == stratum
            ]
            strata[stratum] = _basic_directness_metrics(stratum_rows)
        case_strata: dict[str, Any] = {}
        for stratum in CASE_STRATA:
            stratum_rows = [row for row in method_rows if row.get("case_stratum") == stratum]
            case_strata[stratum] = _basic_directness_metrics(stratum_rows)
        output["methods"][method] = {
            **metrics,
            "identifiability_strata": strata,
            "case_strata": case_strata,
            "legacy_age_adjacency": _legacy_evaluator_metrics(case_by_id, method_rows)
            if method_rows
            else {},
        }
    output["truth_join"] = "post_inference_evaluation_only"
    output["claim_boundary"] = (
        "Controlled-synthetic directness evidence only; age ordering alone is "
        "not a direct-adjacency label, and a calibrated-model or field claim "
        "requires independent segment-level transport information."
    )
    return output


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(list(rows)).to_csv(path, index=False)


def _directed_closure(case: AgeCase) -> dict[tuple[str, str], int]:
    """Return shortest directed path lengths for the complete truth graph."""

    adjacency: dict[str, list[str]] = {}
    nodes: set[str] = set(case.truth.path_by_node)
    for upstream, downstream in case.truth.direct_edges:
        adjacency.setdefault(upstream, []).append(downstream)
        nodes.update((upstream, downstream))
    closure: dict[tuple[str, str], int] = {}
    for source in sorted(nodes):
        queue: list[tuple[str, int]] = [(source, 0)]
        visited = {source}
        while queue:
            current, distance = queue.pop(0)
            for target in adjacency.get(current, []):
                if target in visited:
                    continue
                visited.add(target)
                next_distance = distance + 1
                closure[(source, target)] = next_distance
                queue.append((target, next_distance))
    return closure


def _complete_truth_rows(cases: Sequence[AgeCase]) -> list[dict[str, Any]]:
    """Serialise all direct and reachable truth pairs, including absent candidates."""

    rows: list[dict[str, Any]] = []
    for case in cases:
        direct = set(case.truth.direct_edges)
        closure = _directed_closure(case)
        for (upstream, downstream), path_length in sorted(closure.items()):
            rows.append(
                {
                    "case_id": case.case_id,
                    "seed": case.seed,
                    "split": case.split,
                    "case_stratum": case.case_stratum,
                    "u": upstream,
                    "v": downstream,
                    "relation": "direct_adjacent" if (upstream, downstream) in direct else "transitive_reachable",
                    "truth_path_length": path_length,
                }
            )
    return rows


def _paired_case_bootstrap(
    cases: Sequence[AgeCase],
    joined_rows: Sequence[Mapping[str, Any]],
    *,
    n_resamples: int = 1000,
    seed: int = 2026090817,
) -> dict[str, Any]:
    """Estimate paired arm differences by resampling complete cases.

    Edges within one case share endpoint observations and transport draws, so
    resampling individual rows would understate uncertainty.  The locked
    cases, not their candidate edges, are the bootstrap units.
    """

    case_ids = [case.case_id for case in cases]
    by_method_case: dict[str, dict[str, list[Mapping[str, Any]]]] = {
        method: {case_id: [] for case_id in case_ids} for method in METHODS
    }
    for row in joined_rows:
        method = str(row.get("method"))
        case_id = str(row.get("case_id"))
        if method in by_method_case and case_id in by_method_case[method]:
            by_method_case[method][case_id].append(row)

    rng = np.random.default_rng(int(seed))
    contrasts = {
        "full_bf_minus_order_only": ("full_bf", "order_only"),
        "full_bf_minus_permuted_full_bf": ("full_bf", "permuted_full_bf"),
    }
    distributions: dict[str, list[float]] = {name: [] for name in contrasts}
    if not case_ids:
        return {
            "n_resamples": 0,
            "seed": int(seed),
            "unit": "complete_case",
            "contrasts": {},
        }
    for _ in range(int(n_resamples)):
        sampled_indices = rng.integers(0, len(case_ids), size=len(case_ids))
        sampled_ids = [case_ids[int(index)] for index in sampled_indices]
        aggregate: dict[str, list[Mapping[str, Any]]] = {method: [] for method in METHODS}
        for method in METHODS:
            for case_id in sampled_ids:
                aggregate[method].extend(by_method_case[method][case_id])
        values: dict[str, float | None] = {}
        for method in METHODS:
            metric = _basic_directness_metrics(aggregate[method])
            values[method] = _finite_float(metric.get("directness_pr_auc"))
        for name, (left, right) in contrasts.items():
            left_value = values[left]
            right_value = values[right]
            if left_value is not None and right_value is not None:
                distributions[name].append(float(left_value - right_value))
    result: dict[str, Any] = {
        "n_resamples": int(n_resamples),
        "seed": int(seed),
        "unit": "complete_case",
        "interval": "percentile_95",
        "contrasts": {},
    }
    for name, values in distributions.items():
        if not values:
            result["contrasts"][name] = {
                "n_valid_resamples": 0,
                "estimate": None,
                "ci95_low": None,
                "ci95_high": None,
            }
            continue
        result["contrasts"][name] = {
            "n_valid_resamples": len(values),
            "estimate": float(np.mean(values)),
            "ci95_low": float(np.quantile(values, 0.025)),
            "ci95_high": float(np.quantile(values, 0.975)),
        }
    return result


def _package_truth_counts(cases: Sequence[AgeCase]) -> tuple[int, int]:
    direct_count = sum(len(case.truth.direct_edges) for case in cases)
    reachable_count = sum(len(_directed_closure(case)) for case in cases)
    return int(direct_count), int(reachable_count)


def _write_tier_package(
    output_path: Path,
    cases: Sequence[AgeCase],
    predictions: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Write and audit a T2 package for the full Bayes-factor arm.

    The package is deliberately separate from the CSV metrics table.  The
    truth sidecar contains complete graph denominators, while prediction rows
    contain only inference outputs and explicit ``ABSTAIN`` statuses.
    """

    repo_root = Path(__file__).resolve().parents[3]
    protocol_path = repo_root / "docs" / "AGE_ADJACENCY_TWO_TIER_QA.md"
    inference_path = repo_root / "hydrosheaf" / "validation" / "age_adjacency.py"
    generator_path = Path(__file__).resolve()
    lock_path = repo_root / "uv.lock"
    run_id = output_path.name
    full_predictions = [
        dict(row) for row in predictions if str(row.get("method")) == "full_bf"
    ]
    age_records: dict[tuple[str, str], dict[str, Any]] = {}
    transport_records: list[dict[str, Any]] = []
    for case in cases:
        for candidate in case.candidates:
            for node, age, sigma in (
                (candidate.u, candidate.upstream_age_years, candidate.upstream_sigma_years),
                (candidate.v, candidate.downstream_age_years, candidate.downstream_sigma_years),
            ):
                age_records.setdefault(
                    (case.case_id, node),
                    {
                        "case_id": case.case_id,
                        "node_id": node,
                        "age_years": float(age),
                        "age_sigma_years": float(sigma),
                        "age_covariance_years2": float(candidate.age_covariance_years2),
                        "age_status": "observed",
                        "source_id": "synthetic_age_observations",
                    },
                )
            edge_id = f"{candidate.u}->{candidate.v}"
            if candidate.hop_count >= 2:
                intermediate = [f"{candidate.path_id}-N{int(candidate.u.rsplit('N', 1)[1]) + 1:02d}"]
            else:
                intermediate = [f"counterfactual_intermediate:{candidate.candidate_id}"]
            for hypothesis, travel, sigma, nodes in (
                ("direct", candidate.direct_travel_years, candidate.direct_travel_sigma_years, []),
                ("indirect", candidate.indirect_travel_years, candidate.indirect_travel_sigma_years, intermediate),
            ):
                transport_records.append(
                    {
                        "case_id": case.case_id,
                        "edge_id": edge_id,
                        "u": candidate.u,
                        "v": candidate.v,
                        "hypothesis": hypothesis,
                        "travel_time_years": float(travel),
                        "travel_sigma_years": float(sigma),
                        "evidence_source_id": "synthetic_transport_hypotheses",
                        "independent_of_endpoint_age": True,
                        "path_basis": "one_segment_rtd" if hypothesis == "direct" else "two_segment_rtd_convolution",
                        "intermediate_nodes": nodes,
                    }
                )
    package_predictions: list[dict[str, Any]] = []
    for row in full_predictions:
        scored = str(row.get("score_status")) == "scored"
        original_stratum = str(row.get("identifiability_stratum"))
        package_row: dict[str, Any] = {
            "case_id": str(row["case_id"]),
            "edge_id": f"{row['u']}->{row['v']}",
            "u": str(row["u"]),
            "v": str(row["v"]),
            "direction_probability": float(row.get("direction_probability", 0.5)),
            "prediction_status": "scored" if scored else "ABSTAIN",
            "adjacency_status": "scored" if scored else "insufficient_information",
            "identifiability_stratum": "edge_transport_comparison" if scored else (
                "unidentifiable" if original_stratum == "overlapping" else "transport_censored"
            ),
        }
        if scored:
            package_row["direct_probability"] = float(row["direct_probability"])
            package_row["log_bayes_factor_direct_vs_indirect"] = float(row["log_bayes_factor_direct_vs_indirect"])
        else:
            package_row["flags"] = [str(row.get("abstain_reason") or "directness_unavailable")]
        package_predictions.append(package_row)

    source_entries = []
    for source_id, role, path in (
        ("synthetic_age_observations", "endpoint_age_observations", output_path / "heldout_inference_candidates.csv"),
        ("synthetic_transport_hypotheses", "edge_transport_hypotheses", output_path / "heldout_inference_candidates.csv"),
    ):
        source_entries.append(
            {
                "source_id": source_id,
                "role": role,
                "path": path.name,
                "sha256": _sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )
    direct_count, reachable_count = _package_truth_counts(cases)
    complete_truth_path = output_path / "complete_truth.csv"
    _write_csv(complete_truth_path, _complete_truth_rows(cases))
    truth = {
        "schema": "age-adjacency-truth-v1",
        "run_id": run_id,
        "sealed": True,
        "direct_edge_count": direct_count,
        "reachable_ordered_pair_count": reachable_count,
        "truth_scope": "all held-out generating paths, including pairs absent from candidate set",
        "complete_truth_file": complete_truth_path.name,
        "complete_truth_sha256": _sha256_file(complete_truth_path),
    }
    truth_path = output_path / "truth_artifact.json"
    truth_path.write_text(json.dumps(truth, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    package: dict[str, Any] = {
        "schema": "age-adjacency-two-tier-v1",
        "protocol_id": "age-adjacency-two-tier-v1",
        "run_id": run_id,
        "tier": "T2_segment_transport",
        "claim_tier": "controlled_synthetic_component",
        "metadata": {
            "source_kind": "independent_synthetic_truth",
            "truth_sealed": True,
            "truth_access_mode": "evaluation_only",
            "candidate_set_frozen": True,
            "integrated_scoring_allowed": True,
            "aiken_emulation": False,
            "prediction_input_columns": [
                "upstream_age_years", "downstream_age_years", "upstream_sigma_years",
                "downstream_sigma_years", "age_covariance_years2", "process_sigma_years",
                "direct_travel_years", "indirect_travel_years", "direct_travel_sigma_years",
                "indirect_travel_sigma_years",
            ],
            "evaluation_only_columns": ["truth_label", "relation_label", "is_true_edge", "graph_closure"],
        },
        "units": {
            "age_years": "years",
            "age_sigma_years": "years",
            "age_covariance_years2": "years^2",
            "travel_time_years": "years",
            "travel_sigma_years": "years",
            "probability": "1",
        },
        "provenance": {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "protocol_hash": _sha256_file(protocol_path),
            "generator": {"name": "age_bayes_factor_benchmark", "version": PROTOCOL_NAME, "revision": _sha256_file(generator_path)},
            "inference": {"name": "hydrosheaf.age_adjacency", "version": "explicit_normal_likelihood_v1", "revision": _sha256_file(inference_path)},
            "environment": {"runtime": f"python-{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}", "dependency_lock_hash": _sha256_file(lock_path)},
            "inputs": source_entries,
        },
        "ages": list(age_records.values()),
        "transport_hypotheses": transport_records,
        "predictions": package_predictions,
        "truth_artifact": {
            "available": True,
            "sealed": True,
            "role": "evaluation_only",
            "path": truth_path.name,
            "sha256": _sha256_file(truth_path),
            "complete_truth_file": complete_truth_path.name,
            "complete_truth_sha256": _sha256_file(complete_truth_path),
        },
    }
    package_path = output_path / "package.json"
    package_path.write_text(json.dumps(package, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    script_dir = Path(__file__).resolve().parent
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    from audit_age_adjacency_tiers import audit_package

    qa = audit_package(package_path)
    (output_path / "qa.json").write_text(json.dumps(qa, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "package": package_path.name,
        "truth_artifact": truth_path.name,
        "complete_truth": complete_truth_path.name,
        "qa": qa,
    }


def run_benchmark(
    *,
    output: Path | str = DEFAULT_OUTPUT,
    first_seed: int = 2026090801,
    n_cases: int = 24,
    n_development: int | None = None,
    min_separation_z: float = 1.0,
) -> dict[str, Any]:
    """Run the benchmark and write a reproducible held-out package."""

    output_path = Path(output).expanduser().resolve()
    if output_path.exists() and any(output_path.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output_path}")
    cases = generate_cases(
        first_seed=first_seed,
        n_cases=n_cases,
        n_development=n_development,
    )
    heldout = tuple(case for case in cases if case.split == "locked_test")
    predictions: list[dict[str, Any]] = []
    for case in heldout:
        for method in METHODS:
            predictions.extend(
                score_case(
                    case,
                    method=method,
                    min_separation_z=min_separation_z,
                )
            )
    metrics = evaluate_predictions(heldout, predictions)
    # The case schedule is frozen before test scoring.  Bootstrap contrasts
    # use whole held-out cases as units and are therefore descriptive
    # uncertainty intervals, not a post-hoc tuning device.
    joined_for_bootstrap = _truth_join(
        {case.case_id: case for case in heldout}, predictions
    )
    metrics["paired_case_bootstrap"] = _paired_case_bootstrap(
        heldout,
        joined_for_bootstrap,
    )
    output_path.mkdir(parents=True, exist_ok=True)
    inference_rows = [row for case in heldout for row in case.inference_rows()]
    truth_rows = [
        {
            "case_id": case.case_id,
            "seed": case.seed,
            "split": case.split,
            "candidate_id": candidate.candidate_id,
            "truth_label": case.truth.relation(candidate.candidate_id),
            "identifiability_stratum": candidate.identifiability_stratum,
        }
        for case in heldout
        for candidate in case.candidates
    ]
    _write_csv(output_path / "heldout_inference_candidates.csv", inference_rows)
    _write_csv(output_path / "heldout_truth.csv", truth_rows)
    _write_csv(output_path / "heldout_predictions.csv", predictions)
    _write_csv(output_path / "complete_truth.csv", _complete_truth_rows(heldout))
    (output_path / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    manifest: dict[str, Any] = {
        "schema": "age-bayes-factor-run-manifest-v1",
        "protocol": PROTOCOL_NAME,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {
            "module": str(Path(__file__).resolve()),
            "first_seed": int(first_seed),
            "n_cases": int(n_cases),
            "n_development": int(n_development if n_development is not None else n_cases // 2),
            "n_locked_test_cases": len(heldout),
            "case_strata": list(CASE_STRATA),
            "identifiability_rule": "|direct_travel-indirect_travel| / sqrt(direct_sigma^2 + indirect_sigma^2)",
            "minimum_scorable_separation_z": float(min_separation_z),
        },
        "methods": list(METHODS),
        "truth": {
            "truth_file": "heldout_truth.csv",
            "complete_truth_file": "complete_truth.csv",
            "truth_join": "post_inference_evaluation_only",
            "inference_file_excludes_truth_labels": True,
            "direct_relation": "adjacent generating segment",
            "indirect_relation": "two-segment generating path skip",
        },
        "outputs": {
            "heldout_inference_candidates": "heldout_inference_candidates.csv",
            "heldout_truth": "heldout_truth.csv",
            "heldout_predictions": "heldout_predictions.csv",
            "complete_truth": "complete_truth.csv",
            "metrics": "metrics.json",
        },
        "evidence_boundary": (
            "Independent synthetic truth supports controlled direct-versus-indirect "
            "benchmarking only. No Aiken or Ghana field claim is made."
        ),
        "metrics": metrics,
    }
    tier_package = _write_tier_package(output_path, heldout, predictions)
    manifest["tier_package"] = tier_package
    manifest["complete_truth_counts"] = {
        "direct_edges": _package_truth_counts(heldout)[0],
        "reachable_ordered_pairs": _package_truth_counts(heldout)[1],
    }
    # Metrics include the bootstrap block added above; rewrite it after the
    # package build so the manifest remains the final provenance record.
    (output_path / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    (output_path / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    # Record output hashes after all files except the manifest itself exist.
    output_hashes = {
        path.name: _sha256_file(path)
        for path in sorted(output_path.iterdir())
        if path.is_file() and path.name != "manifest.json"
    }
    manifest["output_sha256"] = output_hashes
    (output_path / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--first-seed", type=int, default=2026090801)
    parser.add_argument("--n-cases", type=int, default=24)
    parser.add_argument("--n-development", type=int, default=None)
    parser.add_argument("--min-separation-z", type=float, default=1.0)
    parser.add_argument("--quick", action="store_true", help="run four cases with two development cases")
    args = parser.parse_args(argv)
    n_cases = 4 if args.quick else args.n_cases
    n_development = 2 if args.quick else args.n_development
    manifest = run_benchmark(
        output=args.output,
        first_seed=args.first_seed,
        n_cases=n_cases,
        n_development=n_development,
        min_separation_z=args.min_separation_z,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


__all__ = [
    "AgeCandidate",
    "AgeCase",
    "AgeTruth",
    "CASE_STRATA",
    "DEFAULT_OUTPUT",
    "IDENTIFIABILITY_STRATA",
    "METHODS",
    "PROTOCOL_NAME",
    "evaluate_predictions",
    "generate_age_case",
    "generate_cases",
    "run_benchmark",
    "score_case",
]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
