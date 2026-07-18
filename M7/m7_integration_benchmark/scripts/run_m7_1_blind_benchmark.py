"""Run the blind, replicated M7.1 sequential integration benchmark.

The development seeds are used only to freeze decision thresholds.  All
reported performance is computed on disjoint test aquifers.  A smaller heavy
tier audits actual topology-posterior, PHREEQC, and Bayesian-aging execution;
it is not used for superiority claims.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import logging
import platform
import subprocess
import time
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import networkx as nx
import numpy as np
from sklearn.metrics import (
    average_precision_score,
    matthews_corrcoef,
    precision_recall_fscore_support,
)
from sklearn.linear_model import LogisticRegression

from m7_1_inference import (
    HELD_OUT_IONS,
    BlindInferenceResult,
    benchmark_config,
    run_blind_inference,
)
from m7_1_truth_model import BlindAquifer, generate_blind_aquifer, poison_truth


ROOT = Path(__file__).resolve().parents[3]
M7_ROOT = Path(__file__).resolve().parents[1]
RESULTS = M7_ROOT / "results" / "m7_1_blind"
METHOD_KEYS = {
    "hydraulic": "hydraulic_probability",
    "age_only": "age_only_probability",
    "sheaf_multievidence": "sheaf_multievidence_probability",
    "chemistry": "chemistry_probability",
    "joint_equal_weight": "joint_probability",
    "joint_logistic": "joint_logistic_probability",
}
FUSION_FEATURES = (
    "hydraulic_probability",
    "age_only_probability",
    "chemistry_probability",
)


def _edge_truth(case: BlindAquifer, result: BlindInferenceResult) -> Tuple[List[int], Dict[str, List[float]]]:
    truth_ids = {f"{u}->{v}" for u, v in case.true_edges}
    labels = [int(str(row["edge_id"]) in truth_ids) for row in result.edge_scores]
    scores = {
        method: [float(row[key]) for row in result.edge_scores]
        for method, key in METHOD_KEYS.items()
    }
    return labels, scores


def _f1(labels: Sequence[int], scores: Sequence[float], threshold: float) -> float:
    predictions = [int(value >= threshold) for value in scores]
    return float(
        precision_recall_fscore_support(
            labels, predictions, average="binary", zero_division=0
        )[2]
    )


def freeze_thresholds(
    development: Sequence[Tuple[BlindAquifer, BlindInferenceResult]]
) -> Dict[str, float]:
    """Choose thresholds on development aquifers only using macro aquifer F1."""

    frozen: Dict[str, float] = {}
    for method in METHOD_KEYS:
        observed = sorted(
            {
                value
                for case, result in development
                for value in _edge_truth(case, result)[1][method]
            }
        )
        if not observed:
            frozen[method] = 0.5
            continue
        candidates = sorted({0.0, 1.0, *observed})
        scored: List[Tuple[float, float]] = []
        for threshold in candidates:
            macro_f1 = float(
                np.mean(
                    [
                        _f1(
                            _edge_truth(case, result)[0],
                            _edge_truth(case, result)[1][method],
                            threshold,
                        )
                        for case, result in development
                    ]
                )
            )
            scored.append((macro_f1, threshold))
        # Conservative deterministic tie-break: select the largest threshold.
        frozen[method] = float(max(scored, key=lambda item: (item[0], item[1]))[1])
    return frozen


def fit_fusion_on_development(
    development: Sequence[Tuple[BlindAquifer, BlindInferenceResult]]
) -> Dict[str, object]:
    """Fit and freeze a regularized evidence stacker on development truth."""

    x: List[List[float]] = []
    y: List[int] = []
    for case, result in development:
        truth_ids = {f"{u}->{v}" for u, v in case.true_edges}
        for row in result.edge_scores:
            x.append([float(row[name]) for name in FUSION_FEATURES])
            y.append(int(str(row["edge_id"]) in truth_ids))
    model = LogisticRegression(
        C=1.0,
        max_iter=1000,
        random_state=719,
    )
    model.fit(np.asarray(x, dtype=float), np.asarray(y, dtype=int))
    return {
        "intercept": float(model.intercept_[0]),
        "coefficients": {
            name: float(value)
            for name, value in zip(FUSION_FEATURES, model.coef_[0])
        },
        "regularization": "L2",
        "C": 1.0,
        "class_weight": None,
        "random_state": 719,
    }


def attach_fusion_scores(
    result: BlindInferenceResult, fusion_model: Mapping[str, object]
) -> None:
    coefficients = fusion_model["coefficients"]
    assert isinstance(coefficients, Mapping)
    for row in result.edge_scores:
        linear = float(fusion_model["intercept"]) + sum(
            float(coefficients[name]) * float(row[name])
            for name in FUSION_FEATURES
        )
        row["joint_logistic_probability"] = 1.0 / (
            1.0 + np.exp(-np.clip(linear, -40.0, 40.0))
        )


def _ece(labels: Sequence[int], scores: Sequence[float], bins: int = 10) -> float:
    labels_arr = np.asarray(labels, dtype=float)
    scores_arr = np.asarray(scores, dtype=float)
    total = max(1, len(labels_arr))
    value = 0.0
    for low, high in zip(np.linspace(0, 1, bins + 1)[:-1], np.linspace(0, 1, bins + 1)[1:]):
        mask = (scores_arr >= low) & (
            (scores_arr <= high) if high == 1.0 else (scores_arr < high)
        )
        if np.any(mask):
            value += float(np.sum(mask) / total) * abs(
                float(np.mean(scores_arr[mask])) - float(np.mean(labels_arr[mask]))
            )
    return value


def _reachability_jaccard(true_edges: Iterable[Tuple[str, str]], predicted_ids: Iterable[str]) -> float:
    true_graph = nx.DiGraph()
    true_graph.add_edges_from(true_edges)
    predicted_graph = nx.DiGraph()
    predicted_graph.add_nodes_from(true_graph.nodes())
    for edge_id in predicted_ids:
        u, v = edge_id.split("->", 1)
        predicted_graph.add_edge(u, v)

    def reachable(graph: nx.DiGraph) -> set[Tuple[str, str]]:
        return {
            (u, v)
            for u in graph.nodes()
            for v in nx.descendants(graph, u)
        }

    truth = reachable(true_graph)
    prediction = reachable(predicted_graph)
    union = truth | prediction
    return float(len(truth & prediction) / len(union)) if union else 1.0


def evaluate_case(
    case: BlindAquifer,
    result: BlindInferenceResult,
    thresholds: Mapping[str, float],
) -> List[Dict[str, object]]:
    labels, scores_by_method = _edge_truth(case, result)
    truth_ids = {f"{u}->{v}" for u, v in case.true_edges}
    candidate_ids = set(result.candidate_edges)
    candidate_recall = len(truth_ids & candidate_ids) / len(truth_ids)
    observations = {
        str(row["site_id"]): row for row in case.observations
    }
    score_by_id = {
        str(row["edge_id"]): row for row in result.edge_scores
    }
    held_out_errors: List[float] = []
    supported_correct: List[bool] = []
    supported_count = 0
    for edge_id in truth_ids & candidate_ids:
        score_row = score_by_id[edge_id]
        u, v = edge_id.split("->", 1)
        predicted_delta = score_row["held_out_reaction_delta"]
        assert isinstance(predicted_delta, Mapping)
        for ion in HELD_OUT_IONS:
            source_value = observations[u].get(ion)
            target_value = observations[v].get(ion)
            if source_value is None or target_value is None:
                continue
            observed_delta = float(target_value) - float(source_value)
            held_out_errors.append(
                float(predicted_delta[ion]) - observed_delta
            )
        truth_process = case.true_processes[edge_id]
        expected_label = {
            "carbonate_weathering": "calcite",
            "silicate_weathering": "albite",
        }.get(truth_process)
        if expected_label is not None:
            supported_count += 1
            extents = score_row["reaction_extents"]
            assert isinstance(extents, Mapping)
            diagnostic = {
                label: abs(float(extents.get(label, 0.0)))
                for label in ("calcite", "albite")
            }
            inferred_label = max(diagnostic, key=diagnostic.get)
            supported_correct.append(inferred_label == expected_label)
    held_out_rmse = (
        float(np.sqrt(np.mean(np.square(held_out_errors))))
        if held_out_errors
        else float("nan")
    )
    process_support_fraction = supported_count / max(1, len(truth_ids & candidate_ids))
    process_accuracy_supported = (
        float(np.mean(supported_correct)) if supported_correct else float("nan")
    )
    rows: List[Dict[str, object]] = []
    for method, scores in scores_by_method.items():
        threshold = float(thresholds[method])
        predictions = [int(value >= threshold) for value in scores]
        precision, recall, f1, _ = precision_recall_fscore_support(
            labels, predictions, average="binary", zero_division=0
        )
        selected_ids = {
            str(row["edge_id"])
            for row, prediction in zip(result.edge_scores, predictions)
            if prediction
        }
        rows.append(
            {
                "seed": case.seed,
                "archetype": case.archetype,
                "age_regime": case.age_regime,
                "method": method,
                "threshold": threshold,
                "candidate_recall": float(candidate_recall),
                "n_candidates": len(labels),
                "n_true_edges": len(truth_ids),
                "n_selected": int(sum(predictions)),
                "precision": float(precision),
                "recall": float(recall),
                "f1": float(f1),
                "mcc": float(matthews_corrcoef(labels, predictions)),
                "pr_auc": float(average_precision_score(labels, scores)),
                "brier": float(np.mean((np.asarray(scores) - np.asarray(labels)) ** 2)),
                "ece10": _ece(labels, scores),
                "structural_hamming_distance": len(truth_ids ^ selected_ids),
                "reachability_jaccard": _reachability_jaccard(case.true_edges, selected_ids),
                "held_out_reaction_rmse": held_out_rmse,
                "reaction_family_support_fraction": process_support_fraction,
                "reaction_family_accuracy_supported": process_accuracy_supported,
            }
        )
    return rows


def _interval(values: Sequence[float], rng: np.random.Generator) -> Dict[str, float]:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return {"mean": float("nan"), "ci95_low": float("nan"), "ci95_high": float("nan")}
    indices = rng.integers(0, arr.size, size=(2000, arr.size))
    boot = np.mean(arr[indices], axis=1)
    return {
        "mean": float(np.mean(arr)),
        "ci95_low": float(np.percentile(boot, 2.5)),
        "ci95_high": float(np.percentile(boot, 97.5)),
    }


def summarise(rows: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    rng = np.random.default_rng(731_991)
    metrics = (
        "candidate_recall",
        "precision",
        "recall",
        "f1",
        "mcc",
        "pr_auc",
        "brier",
        "ece10",
        "structural_hamming_distance",
        "reachability_jaccard",
        "held_out_reaction_rmse",
        "reaction_family_support_fraction",
        "reaction_family_accuracy_supported",
    )
    summary: Dict[str, object] = {}
    for method in METHOD_KEYS:
        method_rows = [row for row in rows if row["method"] == method]
        summary[method] = {
            metric: _interval([float(row[metric]) for row in method_rows], rng)
            for metric in metrics
        }
    hydraulic_by_seed = {
        int(row["seed"]): row for row in rows if row["method"] == "hydraulic"
    }
    paired: Dict[str, object] = {}
    for method in METHOD_KEYS:
        if method == "hydraulic":
            continue
        method_by_seed = {
            int(row["seed"]): row for row in rows if row["method"] == method
        }
        shared = sorted(set(hydraulic_by_seed) & set(method_by_seed))
        paired[method] = {
            metric: _interval(
                [
                    float(method_by_seed[seed][metric])
                    - float(hydraulic_by_seed[seed][metric])
                    for seed in shared
                ],
                rng,
            )
            for metric in ("f1", "mcc", "pr_auc", "brier", "ece10")
        }
    summary["paired_difference_vs_hydraulic"] = paired
    age_regime_summary: Dict[str, object] = {}
    for regime in sorted({str(row["age_regime"]) for row in rows}):
        regime_rows = [row for row in rows if row["age_regime"] == regime]
        age_regime_summary[regime] = {
            method: {
                metric: _interval(
                    [
                        float(row[metric])
                        for row in regime_rows
                        if row["method"] == method
                    ],
                    rng,
                )
                for metric in ("f1", "mcc", "pr_auc")
            }
            for method in METHOD_KEYS
        }
    summary["subgroups_by_age_regime"] = age_regime_summary
    return summary


def _hash_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_value(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unavailable"


def _source_worktree_status() -> List[str]:
    """Return pre-run changes, excluding only this benchmark's artifacts."""

    status = _git_value("status", "--porcelain", "--untracked-files=all")
    if status == "unavailable":
        return ["unavailable"]
    artifact_prefix = str(RESULTS.relative_to(ROOT)).replace("\\", "/") + "/"
    source_changes: List[str] = []
    for line in status.splitlines():
        path = line[3:].strip().strip('"').replace("\\", "/")
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        if not path.startswith(artifact_prefix):
            source_changes.append(line)
    return source_changes


def _runtime_source_files() -> List[Path]:
    """Enumerate all repository inputs capable of changing benchmark results."""

    files = [
        path
        for path in (ROOT / "hydrosheaf").rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    ]
    files.extend(Path(__file__).parent.glob("*.py"))
    files.extend([ROOT / "pyproject.toml", ROOT / "uv.lock"])
    return sorted({path.resolve() for path in files if path.is_file()})


def _package_versions() -> Dict[str, str]:
    names = ("hydrosheaf", "numpy", "scipy", "networkx", "scikit-learn", "pymc", "nutpie", "phreeqpython")
    versions: Dict[str, str] = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _run_case(
    seed: int,
    threshold: float,
    posterior: bool = False,
    heavy: bool = False,
    draws: int = 120,
    fusion_model: Mapping[str, object] | None = None,
) -> Tuple[BlindAquifer, BlindInferenceResult]:
    case = generate_blind_aquifer(seed)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Recharge date .* exceeds input history start")
        result = run_blind_inference(
            case.observations,
            threshold=threshold,
            seed=seed,
            run_posterior=posterior,
            run_heavy_modules=heavy,
            aging_draws=draws,
            fusion_model=fusion_model,
        )
    return case, result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development-replicates", type=int, default=20)
    parser.add_argument("--test-replicates", type=int, default=100)
    parser.add_argument("--heavy-replicates", type=int, default=8)
    parser.add_argument("--aging-draws", type=int, default=1000)
    parser.add_argument("--skip-heavy", action="store_true")
    args = parser.parse_args()
    if args.development_replicates < 1 or args.test_replicates < 1:
        raise ValueError("Development and test replicate counts must be positive.")
    if args.heavy_replicates < 0:
        raise ValueError("Heavy replicate count must be non-negative.")
    logging.getLogger("hydrosheaf.sheaf.topology_refine").setLevel(logging.ERROR)

    started = time.time()
    git_commit_before_run = _git_value("rev-parse", "HEAD")
    pre_run_source_status = _source_worktree_status()
    source_files = _runtime_source_files()
    source_sha256 = {
        str(path.relative_to(ROOT)): _hash_file(path) for path in source_files
    }
    RESULTS.mkdir(parents=True, exist_ok=True)
    development_seeds = list(range(10_000, 10_000 + args.development_replicates))
    test_seeds = list(range(20_000, 20_000 + args.test_replicates))
    heavy_seeds = list(range(30_000, 30_000 + args.heavy_replicates))
    if set(development_seeds) & set(test_seeds):
        raise AssertionError("Development and test seeds must be disjoint.")

    print(f"M7.1 development tier: {len(development_seeds)} aquifers")
    development = [_run_case(seed, threshold=0.5) for seed in development_seeds]
    fusion_model = fit_fusion_on_development(development)
    for _, result in development:
        attach_fusion_scores(result, fusion_model)
    thresholds = freeze_thresholds(development)
    print(f"Frozen thresholds: {thresholds}")

    # Truth-poisoning leakage test: hidden labels must not affect inference.
    leakage_case = generate_blind_aquifer(development_seeds[0])
    poisoned = poison_truth(leakage_case)
    original_result = run_blind_inference(
        leakage_case.observations,
        thresholds["joint_logistic"],
        leakage_case.seed,
        fusion_model=fusion_model,
    )
    poisoned_result = run_blind_inference(
        poisoned.observations,
        thresholds["joint_logistic"],
        poisoned.seed,
        fusion_model=fusion_model,
    )
    leakage_passed = (
        json.dumps(original_result.serializable(), sort_keys=True)
        == json.dumps(poisoned_result.serializable(), sort_keys=True)
    )
    if not leakage_passed:
        raise AssertionError("Truth-poisoning test failed: hidden truth leaked into inference.")

    print(f"M7.1 test tier: {len(test_seeds)} blind aquifers")
    test_pairs = [
        _run_case(
            seed,
            threshold=thresholds["joint_logistic"],
            fusion_model=fusion_model,
        )
        for seed in test_seeds
    ]
    replicate_rows = [
        row
        for case, result in test_pairs
        for row in evaluate_case(case, result, thresholds)
    ]
    summary = summarise(replicate_rows)

    module_totals = {
        name: sum(result.module_calls[name] for _, result in development + test_pairs)
        for name in original_result.module_calls
    }
    heavy_rows: List[Dict[str, object]] = []
    heavy_failures: List[Dict[str, object]] = []
    if not args.skip_heavy:
        print(f"M7.1 heavy module-audit tier: {len(heavy_seeds)} aquifers")
        for seed in heavy_seeds:
            try:
                case, result = _run_case(
                    seed,
                    threshold=thresholds["joint_logistic"],
                    posterior=True,
                    heavy=True,
                    draws=args.aging_draws,
                    fusion_model=fusion_model,
                )
                for name, count in result.module_calls.items():
                    module_totals[name] += count
                node_age_rows = {
                    node: values
                    for node, values in result.bayesian_ages.items()
                    if node != "_diagnostics"
                }
                errors = [
                    abs(float(values["mean_age_years"]) - case.true_ages_years[node])
                    for node, values in node_age_rows.items()
                ]
                coverage = [
                    float(values["age_95_low"])
                    <= case.true_ages_years[node]
                    <= float(values["age_95_high"])
                    for node, values in node_age_rows.items()
                ]
                age_diagnostics = result.bayesian_ages.get("_diagnostics", {})
                age_r_hat = age_diagnostics.get("age_r_hat_max")
                age_ess_bulk = age_diagnostics.get("age_ess_bulk_min")
                age_ess_tail = age_diagnostics.get("age_ess_tail_min")
                age_divergences = age_diagnostics.get("divergences")
                edge_r_hat_raw = result.posterior_diagnostics.get("edge_r_hat", {})
                edge_ess_raw = result.posterior_diagnostics.get("edge_ess", {})
                edge_r_hats = [
                    float(value)
                    for value in edge_r_hat_raw.values()
                    if value is not None and np.isfinite(float(value))
                ] if isinstance(edge_r_hat_raw, Mapping) else []
                edge_ess_values = [
                    float(value)
                    for value in edge_ess_raw.values()
                    if np.isfinite(float(value))
                ] if isinstance(edge_ess_raw, Mapping) else []
                edge_diagnostics_complete = bool(
                    isinstance(edge_r_hat_raw, Mapping)
                    and isinstance(edge_ess_raw, Mapping)
                    and edge_r_hat_raw
                    and edge_ess_raw
                    and len(edge_r_hats) == len(edge_r_hat_raw)
                    and len(edge_ess_values) == len(edge_ess_raw)
                    and set(edge_r_hat_raw) == set(edge_ess_raw)
                )
                posterior_map = result.posterior_diagnostics.get("map_edges", [])
                posterior_graph = nx.DiGraph()
                posterior_graph.add_nodes_from(case.true_ages_years)
                for edge_id in posterior_map:
                    u, v = str(edge_id).split("->", 1)
                    posterior_graph.add_edge(u, v)
                posterior_root = str(
                    max(
                        case.observations,
                        key=lambda row: float(row["head_meas"]),
                    )["site_id"]
                )
                posterior_reachable = {posterior_root} | nx.descendants(
                    posterior_graph, posterior_root
                )
                posterior_root_reaches_all = (
                    posterior_reachable == set(posterior_graph.nodes)
                )
                posterior_is_dag = nx.is_directed_acyclic_graph(
                    posterior_graph
                )
                posterior_is_connected = nx.is_weakly_connected(
                    posterior_graph
                )
                posterior_max_out_degree = max(
                    dict(posterior_graph.out_degree()).values(),
                    default=0,
                )
                posterior_structure_valid = bool(
                    posterior_graph.number_of_edges() >= 7
                    and posterior_is_dag
                    and posterior_is_connected
                    and posterior_max_out_degree <= 2
                    and posterior_root_reaches_all
                )
                posterior_r_hat = result.posterior_diagnostics.get(
                    "n_edges_r_hat"
                )
                posterior_ess = result.posterior_diagnostics.get("n_edges_ess")
                posterior_valid = bool(
                    posterior_r_hat is not None
                    and np.isfinite(float(posterior_r_hat))
                    and float(posterior_r_hat) <= 1.01
                    and posterior_ess is not None
                    and float(posterior_ess) >= 400.0
                    and edge_diagnostics_complete
                    and edge_r_hats
                    and max(edge_r_hats) <= 1.01
                    and edge_ess_values
                    and min(edge_ess_values) >= 400.0
                    and posterior_structure_valid
                )
                age_valid = bool(
                    len(errors) == len(case.true_ages_years)
                    and age_r_hat is not None
                    and np.isfinite(float(age_r_hat))
                    and float(age_r_hat) <= 1.01
                    and age_ess_bulk is not None
                    and float(age_ess_bulk) >= 400.0
                    and age_ess_tail is not None
                    and float(age_ess_tail) >= 400.0
                    and age_divergences == 0
                )
                phreeqc_valid = bool(
                    result.phreeqc_summary.get("success_fraction") == 1.0
                    and int(
                        result.phreeqc_summary.get(
                            "n_edge_fits_constrained", 0
                        )
                    )
                    > 0
                )
                coupling_valid = bool(
                    result.phreeqc_summary.get("aging_graph_source")
                    == "topology_posterior_map"
                    and "age_feedback_changed_edge_count"
                    in result.phreeqc_summary
                    and phreeqc_valid
                )
                heavy_rows.append(
                    {
                        "seed": seed,
                        "posterior_acceptance_rate": result.posterior_diagnostics.get(
                            "acceptance_rate"
                        ),
                        "posterior_n_edges_rhat": result.posterior_diagnostics.get(
                            "n_edges_r_hat"
                        ),
                        "posterior_n_edges_ess": result.posterior_diagnostics.get(
                            "n_edges_ess"
                        ),
                        "posterior_edge_rhat_max": max(edge_r_hats)
                        if edge_r_hats
                        else None,
                        "posterior_edge_rhat_fraction_le_1_01": float(
                            np.mean(np.asarray(edge_r_hats) <= 1.01)
                        )
                        if edge_r_hats
                        else None,
                        "posterior_edge_ess_min": min(edge_ess_values)
                        if edge_ess_values
                        else None,
                        "posterior_edge_ess_fraction_ge_400": float(
                            np.mean(np.asarray(edge_ess_values) >= 400.0)
                        )
                        if edge_ess_values
                        else None,
                        "posterior_edge_diagnostics_complete": edge_diagnostics_complete,
                        "posterior_map_n_edges": posterior_graph.number_of_edges(),
                        "posterior_map_is_dag": posterior_is_dag,
                        "posterior_map_is_weakly_connected": posterior_is_connected,
                        "posterior_map_max_out_degree": posterior_max_out_degree,
                        "posterior_map_root_node": posterior_root,
                        "posterior_map_root_reaches_all_nodes": posterior_root_reaches_all,
                        "posterior_structure_valid": posterior_structure_valid,
                        "posterior_diagnostics_valid": posterior_valid,
                        "marginal_edge_entropy": result.posterior_diagnostics.get(
                            "marginal_edge_entropy"
                        ),
                        "phreeqc_success_fraction": result.phreeqc_summary.get(
                            "success_fraction"
                        ),
                        "phreeqc_constrained_edge_fits": result.phreeqc_summary.get(
                            "n_edge_fits_constrained"
                        ),
                        "phreeqc_mean_absolute_objective_change": result.phreeqc_summary.get(
                            "mean_absolute_objective_change"
                        ),
                        "phreeqc_diagnostics_valid": phreeqc_valid,
                        "aging_graph_source": result.phreeqc_summary.get(
                            "aging_graph_source"
                        ),
                        "age_feedback_changed_edge_count": result.phreeqc_summary.get(
                            "age_feedback_changed_edge_count"
                        ),
                        "age_feedback_edge_jaccard": result.phreeqc_summary.get(
                            "age_feedback_edge_jaccard"
                        ),
                        "sequential_coupling_executed": coupling_valid,
                        "bayesian_age_rhat_max": age_r_hat,
                        "bayesian_age_ess_bulk_min": age_ess_bulk,
                        "bayesian_age_ess_tail_min": age_ess_tail,
                        "bayesian_age_divergences": age_divergences,
                        "bayesian_age_diagnostics_valid": age_valid,
                        "bayesian_age_mae_years": (
                            float(np.mean(errors)) if errors and age_valid else None
                        ),
                        "bayesian_age_interval_coverage": (
                            float(np.mean(coverage))
                            if coverage and age_valid
                            else None
                        ),
                        "exploratory_bayesian_age_mae_years": (
                            float(np.mean(errors)) if errors else None
                        ),
                        "exploratory_bayesian_age_interval_coverage": (
                            float(np.mean(coverage)) if coverage else None
                        ),
                        "bayesian_age_n_nodes": len(errors),
                    }
                )
            except Exception as exc:
                heavy_failures.append(
                    {"seed": seed, "type": type(exc).__name__, "message": str(exc)}
                )

    _write_csv(RESULTS / "replicate_metrics.csv", replicate_rows)
    _write_csv(RESULTS / "heavy_module_audit.csv", heavy_rows)
    (RESULTS / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    config_payload = vars(benchmark_config()).copy()
    database_path = Path(str(config_payload["phreeqc_database"]))
    try:
        config_payload["phreeqc_database"] = str(database_path.relative_to(ROOT))
    except ValueError:
        config_payload["phreeqc_database"] = database_path.name
    runtime_input_sha256_after_run = {
        str(path.relative_to(ROOT)): _hash_file(path)
        for path in _runtime_source_files()
    }
    runtime_inputs_unchanged = source_sha256 == runtime_input_sha256_after_run
    git_commit_after_run = _git_value("rev-parse", "HEAD")
    manifest = {
        "benchmark": "M7.1 blind replicated integration benchmark",
        "claim_scope": {
            "tier_a": "replicated held-out synthetic capability test",
            "tier_b": "sequential module-execution and convergence audit; not a superiority test",
            "field_validation": False,
        },
        "development_seeds": development_seeds,
        "test_seeds": test_seeds,
        "heavy_seeds": [] if args.skip_heavy else heavy_seeds,
        "thresholds_frozen_before_test_evaluation": thresholds,
        "fusion_model_frozen_before_test_evaluation": fusion_model,
        "truth_poisoning_test_passed": leakage_passed,
        "module_call_totals": module_totals,
        "heavy_failures": heavy_failures,
        "heavy_execution_successes": len(heavy_rows),
        "heavy_valid_posterior_cases": sum(
            bool(row["posterior_diagnostics_valid"]) for row in heavy_rows
        ),
        "heavy_valid_phreeqc_cases": sum(
            bool(row["phreeqc_diagnostics_valid"]) for row in heavy_rows
        ),
        "heavy_valid_bayesian_age_cases": sum(
            bool(row["bayesian_age_diagnostics_valid"]) for row in heavy_rows
        ),
        "heavy_sequential_coupling_cases": sum(
            bool(row["sequential_coupling_executed"]) for row in heavy_rows
        ),
        "configuration": config_payload,
        "git_commit": git_commit_before_run,
        "git_commit_before_run": git_commit_before_run,
        "git_commit_after_run": git_commit_after_run,
        "git_commit_unchanged_during_run": (
            git_commit_before_run == git_commit_after_run
        ),
        "source_worktree_clean_before_run": not pre_run_source_status,
        "source_worktree_status_before_run": pre_run_source_status,
        "runtime_inputs_unchanged_during_run": runtime_inputs_unchanged,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "package_versions": _package_versions(),
        "runtime_input_file_count": len(source_sha256),
        "runtime_input_sha256": source_sha256,
        "elapsed_seconds": time.time() - started,
    }
    (RESULTS / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"summary": summary, "manifest": manifest}, indent=2, default=str))
    if not args.skip_heavy and heavy_failures:
        raise RuntimeError(
            f"{len(heavy_failures)} heavy-tier aquifers failed; see {RESULTS / 'manifest.json'}"
        )
    if not runtime_inputs_unchanged or git_commit_before_run != git_commit_after_run:
        raise RuntimeError(
            "Runtime inputs or Git commit changed during execution; benchmark artifact is invalid."
        )


if __name__ == "__main__":
    main()
