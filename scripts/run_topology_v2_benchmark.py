"""Run the isolated topology-v2 calibration and held-out audit.

The script reads existing M4/M7 artifacts but writes only to a new run
directory under ``.codex_work/topology-v2`` by default.  It never modifies
legacy benchmark or manuscript outputs.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import subprocess
from typing import Iterable, Mapping, Sequence

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.graph.build import infer_edges_from_coordinates
from hydrosheaf.validation.topology_v2 import (
    MODEL_FEATURE_SETS,
    TopologyV2Config,
    TopologyV2LogisticCalibrator,
    bootstrap_case_metric_ci,
    build_topology_v2_feature_rows,
    generate_topology_v2_candidate_universe,
    topology_v2_metrics,
    tune_topology_thresholds,
)


DEFAULT_M7_ROOT = (
    PROJECT_ROOT
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "results"
    / "m7_3_locked"
)
DEFAULT_SAVAGE_ROOT = (
    PROJECT_ROOT
    / "M4"
    / "m4_topology_benchmark"
    / "results"
    / "public_archives"
    / "savage"
)
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / ".codex_work" / "topology-v2"

SPLITS: dict[str, tuple[str, ...]] = {
    "development": ("development_5201", "development_5202", "development_5203"),
    "validation": ("development_5204", "development_5205", "development_5206"),
    "locked_test": tuple(f"locked_test_{seed}" for seed in range(5301, 5313)),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _write_json(path: Path, payload: object) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _hash_payload(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _git_value(*args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _git_manifest() -> dict[str, object]:
    return {
        "head": _git_value("rev-parse", "HEAD"),
        "branch": _git_value("branch", "--show-current"),
        "status_porcelain": _git_value("status", "--short"),
        "worktree_dirty": bool(_git_value("status", "--porcelain")),
    }


def _normalise_truth(rows: Iterable[Mapping[str, object]]) -> set[str]:
    return {
        f"{str(row['u'])}->{str(row['v'])}"
        for row in rows
        if row.get("u") is not None and row.get("v") is not None
    }


def _case_paths(m7_root: Path, case_id: str) -> tuple[Path, Path]:
    case_root = m7_root / "cases" / case_id
    return case_root / "blind_observations.csv", case_root / "heldout_truth.csv"


def _load_blind_observations(m7_root: Path, case_id: str) -> list[dict[str, str]]:
    observations_path, truth_path = _case_paths(m7_root, case_id)
    if not observations_path.exists() or not truth_path.exists():
        raise FileNotFoundError(f"M7 case is incomplete: {case_id}")
    return _read_csv(observations_path)


def _load_scoring_truth(m7_root: Path, case_id: str) -> set[str]:
    _observations_path, truth_path = _case_paths(m7_root, case_id)
    return _normalise_truth(_read_csv(truth_path))


def _confusion(reference: set[str], inferred: set[str]) -> dict[str, object]:
    tp = len(reference & inferred)
    fp = len(inferred - reference)
    fn = len(reference - inferred)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "n_reference_edges": len(reference),
        "n_inferred_edges": len(inferred),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def _summary(values: Sequence[float]) -> dict[str, object]:
    if not values:
        return {"n": 0, "min": None, "median": None, "mean": None, "max": None}
    array = np.asarray(values, dtype=float)
    return {
        "n": int(len(array)),
        "min": float(np.min(array)),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "max": float(np.max(array)),
    }


def _edge_lengths(
    edge_ids: Iterable[str], coordinates: Mapping[str, tuple[float, float]]
) -> list[float]:
    values: list[float] = []
    for edge_id in edge_ids:
        if "->" not in edge_id:
            continue
        u, v = edge_id.split("->", 1)
        if u not in coordinates or v not in coordinates:
            continue
        values.append(float(np.hypot(
            coordinates[u][0] - coordinates[v][0],
            coordinates[u][1] - coordinates[v][1],
        )))
    return values


def reproduce_savage_baseline(savage_root: Path) -> dict[str, object]:
    """Reproduce v1 in memory and return an audit; do not write legacy paths."""

    node_path = savage_root / "modpath_node_mapping.csv"
    reference_path = savage_root / "modpath_reference_edges.csv"
    node_rows = _read_csv(node_path)
    reference_rows = _read_csv(reference_path)
    reference = _normalise_truth(reference_rows)
    samples: list[dict[str, object]] = []
    coordinates: dict[str, tuple[float, float]] = {}
    z_values: dict[str, float] = {}
    for row in node_rows:
        node_id = str(row["node_id"])
        x = float(row["x"])
        y = float(row["y"])
        z = _float(row.get("z"))
        z_proxy = 0.0 if z is None else z
        samples.append(
            {
                "site_id": node_id,
                "lat": y,
                "lon": x,
                "elevation": z_proxy,
            }
        )
        coordinates[node_id] = (x, y)
        z_values[node_id] = z_proxy
    inferred = {
        f"{edge.u}->{edge.v}"
        for edge in infer_edges_from_coordinates(
            samples,
            max_neighbors=2,
            allow_uphill=False,
        )
    }
    metrics = _confusion(reference, inferred)

    false_negatives = sorted(reference - inferred)
    false_positives = sorted(inferred - reference)
    fn_drops = {
        edge_id: z_values[edge_id.split("->", 1)[0]] - z_values[edge_id.split("->", 1)[1]]
        for edge_id in false_negatives
    }
    downhill_fn = [edge_id for edge_id, drop in fn_drops.items() if drop > 0.0]
    uphill_proxy_fn = [edge_id for edge_id, drop in fn_drops.items() if drop <= 0.0]
    reference_sinks = sorted(edge_id.split("->", 1)[1] for edge_id in reference)
    sink_set = set(reference_sinks)
    fp_to_sink = [edge_id for edge_id in false_positives if edge_id.split("->", 1)[1] in sink_set]
    fp_to_non_sink = [edge_id for edge_id in false_positives if edge_id.split("->", 1)[1] not in sink_set]

    all_nodes = sorted(coordinates)
    downhill_ranks: list[int] = []
    for edge_id in downhill_fn:
        u, v = edge_id.split("->", 1)
        ranked = sorted(
            (
                float(np.hypot(
                    coordinates[u][0] - coordinates[target][0],
                    coordinates[u][1] - coordinates[target][1],
                )),
                target,
            )
            for target in all_nodes
            if target != u and z_values[target] < z_values[u]
        )
        for rank, (_distance, target) in enumerate(ranked, start=1):
            if target == v:
                downhill_ranks.append(rank)
                break

    v2_universe = generate_topology_v2_candidate_universe(samples)
    return {
        "baseline_id": "topology_v1_baseline",
        "generator": "hydrosheaf.graph.build.infer_edges_from_coordinates",
        "parameters": {"max_neighbors": 2, "allow_uphill": False, "elevation_as_head_proxy": True},
        "source_hashes": {
            str(node_path): _sha256(node_path),
            str(reference_path): _sha256(reference_path),
        },
        "metrics": metrics,
        "expected_frozen_metrics": {
            "tp": 147,
            "fp": 155,
            "fn": 27,
            "precision": 147 / 302,
            "recall": 147 / 174,
            "f1": 2 * (147 / 302) * (147 / 174) / ((147 / 302) + (147 / 174)),
        },
        "matches_frozen_metrics": (
            metrics["tp"] == 147
            and metrics["fp"] == 155
            and metrics["fn"] == 27
        ),
        "candidate_universe_audit": {
            "v2_algorithm": v2_universe.algorithm,
            "v2_candidate_count": len(v2_universe.edges),
            "v2_candidate_graph_recall": v2_universe.candidate_graph_recall(reference),
            "v2_rejection_counts": dict(v2_universe.rejection_counts),
        },
        "failure_analysis": {
            "false_negative_count": len(false_negatives),
            "false_negatives_uphill_under_elevation_proxy": len(uphill_proxy_fn),
            "false_negatives_downhill_under_elevation_proxy": len(downhill_fn),
            "false_negative_proxy_head_drop_summary": _summary(list(fn_drops.values())),
            "downhill_false_negative_nearest_rank_summary": _summary(
                [float(rank) for rank in downhill_ranks]
            ),
            "reference_sink_count": len(sink_set),
            "false_positives_to_reference_sinks": len(fp_to_sink),
            "false_positives_to_non_sink_targets": len(fp_to_non_sink),
            "interpretation": (
                "The legacy downhill gate removes uphill-under-elevation reference edges "
                "before scoring; remaining false positives are dominated by fan-in/receptor "
                "ambiguity rather than a candidate-universe proof of connectivity."
            ),
        },
        "distance_statistics_coordinate_units": {
            "reference": _summary(_edge_lengths(reference, coordinates)),
            "inferred": _summary(_edge_lengths(inferred, coordinates)),
        },
        "protection": {
            "legacy_outputs_modified": False,
            "legacy_output_root": str(savage_root),
        },
    }


def _load_m7_cache(m7_root: Path) -> dict[str, dict[str, object]]:
    cache: dict[str, dict[str, object]] = {}
    for split, case_ids in SPLITS.items():
        for case_id in case_ids:
            observations = _load_blind_observations(m7_root, case_id)
            universe = generate_topology_v2_candidate_universe(observations)
            features = build_topology_v2_feature_rows(
                universe,
                observations,
                case_id=case_id,
            )
            # Release the sealed reference only after blind candidate and feature
            # construction has completed.  It is used solely to create labels.
            truth = _load_scoring_truth(m7_root, case_id)
            labels = np.asarray(
                [int(row.edge_id in truth) for row in features],
                dtype=float,
            )
            cache[case_id] = {
                "split": split,
                "observations": observations,
                "truth": truth,
                "universe": universe,
                "features": features,
                "labels": labels,
                "candidate_recall": universe.candidate_graph_recall(truth),
            }
    return cache


def _concatenate(cache: Mapping[str, Mapping[str, object]], case_ids: Sequence[str]):
    rows = []
    labels = []
    cases = []
    for case_id in case_ids:
        case = cache[case_id]
        case_rows = list(case["features"])
        case_labels = np.asarray(case["labels"], dtype=float)
        rows.extend(case_rows)
        labels.extend(case_labels.tolist())
        cases.extend([case_id] * len(case_rows))
    return tuple(rows), np.asarray(labels, dtype=float), tuple(cases)


def _feature_coverage(rows, feature_names: Sequence[str]) -> dict[str, float]:
    coverage: dict[str, float] = {}
    for name in feature_names:
        count = sum(
            1
            for row in rows
            if row.features.get(name) is not None and np.isfinite(float(row.features[name]))
        )
        coverage[name] = float(count / len(rows)) if rows else 0.0
    return coverage


def _probability_curve_rows(
    labels: Sequence[int | float],
    probabilities: Sequence[float],
) -> list[dict[str, object]]:
    """Return diagnostic PR/ROC operating curves.

    These curves are diagnostics, not selected inference thresholds.  The
    locked operating point is still chosen only by ``tune_topology_thresholds``
    on the validation split.
    """

    y = np.asarray(labels, dtype=float)
    p = np.asarray(probabilities, dtype=float)
    thresholds = np.unique(np.concatenate((np.linspace(0.01, 0.99, 99), p)))
    rows: list[dict[str, object]] = []
    positive_count = int(np.sum(y == 1.0))
    negative_count = int(np.sum(y == 0.0))
    for threshold in sorted(float(value) for value in thresholds if 0.0 < value < 1.0):
        selected = p >= threshold
        tp = int(np.sum(selected & (y == 1.0)))
        fp = int(np.sum(selected & (y == 0.0)))
        precision = tp / (tp + fp) if tp + fp else 1.0
        recall = tp / positive_count if positive_count else 0.0
        fpr = fp / negative_count if negative_count else 0.0
        rows.append(
            {
                "threshold": threshold,
                "precision": precision,
                "recall": recall,
                "fdr": fp / (tp + fp) if tp + fp else 0.0,
                "true_positive_rate": recall,
                "false_positive_rate": fpr,
                "tp": tp,
                "fp": fp,
                "n_selected": int(np.sum(selected)),
            }
        )
    return rows


def _per_case_metrics(
    cache: Mapping[str, Mapping[str, object]],
    case_ids: Sequence[str],
    feature_names: Sequence[str],
    calibrator: TopologyV2LogisticCalibrator,
    policy,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_id in case_ids:
        case = cache[case_id]
        case_rows = tuple(case["features"])
        probabilities = calibrator.predict_proba(case_rows)
        metrics = topology_v2_metrics(
            case["labels"],
            probabilities,
            policy=policy,
            candidate_recall=case["candidate_recall"],
        )
        rows.append(
            {
                "case_id": case_id,
                "split": case["split"],
                "candidate_recall": case["candidate_recall"],
                "pr_auc": metrics["pr_auc"],
                "roc_auc": metrics["roc_auc"],
                "brier": metrics["brier"],
                "ece": metrics["ece"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "fdr": metrics["fdr"],
                "abstention_rate": metrics["abstention_rate"],
                "coverage": metrics["coverage"],
                "tp": metrics["tp"],
                "fp": metrics["fp"],
                "fn": metrics["fn"],
                "tn": metrics["tn"],
            }
        )
    return rows


def run_model(
    model_name: str,
    cache: Mapping[str, Mapping[str, object]],
    *,
    bootstrap_replicates: int,
) -> tuple[dict[str, object], list[dict[str, object]], dict[str, object]]:
    feature_names = MODEL_FEATURE_SETS[model_name]
    development_rows, development_labels, _ = _concatenate(cache, SPLITS["development"])
    validation_rows, validation_labels, _ = _concatenate(cache, SPLITS["validation"])
    test_rows, test_labels, test_case_ids = _concatenate(cache, SPLITS["locked_test"])
    fit_hash = _hash_payload(
        [
            {
                "case_id": row.case_id,
                "edge_id": row.edge_id,
                "label": int(label),
            }
            for row, label in zip(development_rows, development_labels)
        ]
    )
    calibrator = TopologyV2LogisticCalibrator(feature_names=feature_names, l2=0.25)
    calibrator.fit(
        development_rows,
        development_labels,
        scope="held_out_calibration",
        independent=True,
        generator_id="M7.3_external_MODFLOW_MODPATH_blind_observations",
        split_id="development_5201_5202_5203",
        dataset_hash=fit_hash,
    )
    validation_probabilities = calibrator.predict_proba(validation_rows)
    policy, threshold_audit = tune_topology_thresholds(
        validation_labels,
        validation_probabilities,
        target_fdr=0.25,
        minimum_gap=0.10,
    )
    test_probabilities = calibrator.predict_proba(test_rows)
    validation_metrics = topology_v2_metrics(
        validation_labels,
        validation_probabilities,
        policy=policy,
        candidate_recall=float(
            np.mean([cache[case_id]["candidate_recall"] for case_id in SPLITS["validation"]])
        ),
    )
    test_metrics = topology_v2_metrics(
        test_labels,
        test_probabilities,
        policy=policy,
        candidate_recall=float(
            np.mean([cache[case_id]["candidate_recall"] for case_id in SPLITS["locked_test"]])
        ),
    )
    per_case = _per_case_metrics(
        cache,
        SPLITS["locked_test"],
        feature_names,
        calibrator,
        policy,
    )
    confidence_intervals = {
        metric: bootstrap_case_metric_ci(
            test_labels,
            test_probabilities,
            test_case_ids,
            metric=metric,
            policy=policy,
            n_bootstrap=bootstrap_replicates,
            seed=20260922 + list(MODEL_FEATURE_SETS).index(model_name),
        )
        for metric in ("pr_auc", "brier", "fdr", "recall")
    }
    test_coverage = _feature_coverage(test_rows, feature_names)
    unsupported_features = [name for name, coverage in test_coverage.items() if coverage == 0.0]
    record = {
        "model": model_name,
        "feature_names": list(feature_names),
        "calibrator": calibrator.to_dict(),
        "validation": {
            "metrics": validation_metrics,
            "threshold_policy": policy.to_dict(),
            "threshold_grid_evaluated": len(threshold_audit),
            "target_fdr_met": bool(float(validation_metrics["fdr"]) <= 0.25),
            "selection_mode": (
                "target_fdr_constrained"
                if float(validation_metrics["fdr"]) <= 0.25
                else "nonempty_lowest_fdr_fallback_target_infeasible"
            ),
        },
        "locked_test": {
            "metrics": test_metrics,
            "whole_case_bootstrap": confidence_intervals,
            "n_cases": len(SPLITS["locked_test"]),
            "probability_curve": _probability_curve_rows(test_labels, test_probabilities),
        },
        "feature_coverage_locked_test": test_coverage,
        "unsupported_feature_names_locked_test": unsupported_features,
        "candidate_graph_recall_by_locked_case": {
            case_id: cache[case_id]["candidate_recall"]
            for case_id in SPLITS["locked_test"]
        },
        "interpretation_guardrail": (
            "A model with zero locked-test coverage for an added feature is not evidence "
            "that the feature has no effect; its contribution is unidentifiable in this panel."
            if unsupported_features
            else "All selected feature channels have non-zero locked-test coverage."
        ),
    }
    return record, per_case, {"policy": policy.to_dict(), "audit": threshold_audit}


def _write_comparison(path: Path, records: Mapping[str, Mapping[str, object]]) -> None:
    rows = []
    for model_name, record in records.items():
        metrics = record["locked_test"]["metrics"]
        rows.append(
            {
                "model": model_name,
                "pr_auc": metrics["pr_auc"],
                "roc_auc": metrics["roc_auc"],
                "brier": metrics["brier"],
                "ece": metrics["ece"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "fdr": metrics["fdr"],
                "mcc_resolved": metrics["mcc_resolved"],
                "abstention_rate": metrics["abstention_rate"],
                "coverage": metrics["coverage"],
                "candidate_graph_recall": metrics["candidate_graph_recall"],
            }
        )
    _write_csv(path, rows, tuple(rows[0]) if rows else ("model",))


def _interpretation(records: Mapping[str, Mapping[str, object]]) -> str:
    lines = [
        "# Topology-v2 locked-test interpretation",
        "",
        "This report is a controlled M7.3 MODFLOW/MODPATH synthetic comparison. "
        "It is not field connectivity validation and does not establish a universal "
        "topology advantage.",
        "",
        "## What separates a true edge from a downhill-plausible edge?",
        "",
        "The primary comparison is the change in PR-AUC, Brier/ECE, FDR, abstention, "
        "and whole-case uncertainty as feature groups are added. Head direction and "
        "gradient are evidence, not a candidate-universe gate. Screen, capacity, "
        "geology, chemistry, and tracer channels are credited only when their locked "
        "case coverage is non-zero.",
        "",
        "| model | PR-AUC | Brier | ECE | precision | recall | FDR | abstention |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for model_name, record in records.items():
        metrics = record["locked_test"]["metrics"]
        lines.append(
            f"| {model_name} | {metrics['pr_auc']:.4f} | {metrics['brier']:.4f} | "
            f"{metrics['ece']:.4f} | {metrics['precision']:.4f} | {metrics['recall']:.4f} | "
            f"{metrics['fdr']:.4f} | {metrics['abstention_rate']:.4f} |"
        )
    lines.extend(
        [
            "",
            "Model E and later must be read together with their feature coverage. "
            "A zero-capacity channel is a missing-data result, not a negative hydraulic "
            "finding. Fan-in ambiguity remains a structural source of false positives "
            "unless independent receptor/sink evidence is available, and such evidence "
            "must not be introduced into the independent path.",
            "",
            "The locked-test result is therefore reported as conditional evidence: "
            "whether the declared observations improve ranking and calibrated selective "
            "decisions on these generated cases. It is not a claim of exact field flow "
            "direction, field age truth, or general superiority over all aquifers.",
            "",
        ]
    )
    return "\n".join(lines)


def run_audit(
    *,
    output_dir: Path,
    m7_root: Path,
    savage_root: Path,
    bootstrap_replicates: int,
    skip_savage: bool,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=False)
    design_path = PROJECT_ROOT / "docs" / "topology_v2_design.md"
    (output_dir / "design_doc.md").write_text(design_path.read_text(encoding="utf-8"), encoding="utf-8")

    provenance = {
        "run_id": output_dir.name,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(PROJECT_ROOT),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git": _git_manifest(),
        "source_files": {
            "topology_v2_module": str(PROJECT_ROOT / "hydrosheaf" / "validation" / "topology_v2.py"),
            "topology_v2_module_sha256": _sha256(PROJECT_ROOT / "hydrosheaf" / "validation" / "topology_v2.py"),
            "runner": str(Path(__file__).resolve()),
            "runner_sha256": _sha256(Path(__file__).resolve()),
            "design_doc_sha256": _sha256(design_path),
        },
        "inputs": {
            "m7_manifest": str(m7_root / "manifest.json"),
            "m7_manifest_sha256": _sha256(m7_root / "manifest.json"),
            "savage_root": str(savage_root),
        },
        "output_protection": {
            "isolated_output_directory": True,
            "legacy_outputs_overwritten": False,
            "manuscript_outputs_overwritten": False,
        },
    }
    _write_json(output_dir / "provenance.json", provenance)

    if skip_savage:
        baseline = {"skipped": True, "reason": "--skip-savage"}
    else:
        baseline = reproduce_savage_baseline(savage_root)
    _write_json(output_dir / "baseline_reproduction.json", baseline)
    if not skip_savage:
        _write_json(
            output_dir / "candidate_generation_audit.json",
            {
                "savage_v1": baseline["candidate_universe_audit"],
                "savage_failure_analysis": baseline["failure_analysis"],
                "m7_default_policy": {
                    "candidate_generator": "truth_blind_all_pairs_soft_admissibility_v2",
                    "max_neighbors": None,
                    "max_distance_m": None,
                    "hard_direction": False,
                    "candidate_recall_expected": 1.0,
                },
            },
        )

    split_manifest = {
        "generator": "M7.3_external_MODFLOW_MODPATH_blind_observations",
        "truth_blind_generation": True,
        "development_cases": list(SPLITS["development"]),
        "validation_cases": list(SPLITS["validation"]),
        "locked_test_cases": list(SPLITS["locked_test"]),
        "split_unit": "whole_case_aquifer",
        "no_edge_level_split": True,
        "m7_manifest_claimed_development_cases": 6,
        "m7_manifest_claimed_locked_test_cases": 12,
        "truth_files_used_only_after_generation": ["heldout_truth.csv", "modpath_pathline_truth.csv"],
    }
    _write_json(output_dir / "split_manifest.json", split_manifest)

    leakage_audit = {
        "independent_mode": True,
        "modpath_edges_used_in_candidate_generation": False,
        "modpath_edges_used_in_feature_construction": False,
        "modpath_edges_used_in_calibration": False,
        "prior_assisted_edges_used_in_independent_summary": False,
        "truth_loaded_for_scoring_after_blind_generation": True,
        "truth_fields_rejected_by_candidate_and_feature_boundary": True,
        "note": "M7 truth is loaded by the runner only to create scoring labels after blind rows/features exist.",
    }
    _write_json(output_dir / "leakage_audit.json", leakage_audit)

    cache = _load_m7_cache(m7_root)
    candidate_recall_by_case = {
        case_id: cache[case_id]["candidate_recall"]
        for case_id in (*SPLITS["development"], *SPLITS["validation"], *SPLITS["locked_test"])
    }
    candidate_audit = {
        "savage_v1": None if skip_savage else baseline["candidate_universe_audit"],
        "m7_cases": {
            case_id: {
                "split": cache[case_id]["split"],
                "n_nodes": len(cache[case_id]["universe"].nodes),
                "n_candidate_edges": len(cache[case_id]["universe"].edges),
                "candidate_graph_recall": cache[case_id]["candidate_recall"],
                "rejection_counts": dict(cache[case_id]["universe"].rejection_counts),
                "truth_blind_generation": True,
            }
            for case_id in candidate_recall_by_case
        },
        "m7_summary": {
            "candidate_count_min": min(
                len(cache[case_id]["universe"].edges) for case_id in candidate_recall_by_case
            ),
            "candidate_count_max": max(
                len(cache[case_id]["universe"].edges) for case_id in candidate_recall_by_case
            ),
            "candidate_recall_min": float(min(candidate_recall_by_case.values())),
            "candidate_recall_mean": float(np.mean(list(candidate_recall_by_case.values()))),
            "physical_gates_default": "soft_direction_soft_screen_soft_aquifer",
        },
    }
    _write_json(output_dir / "candidate_generation_audit.json", candidate_audit)
    _write_json(
        output_dir / "m7_candidate_recall.json",
        {
            "by_case": candidate_recall_by_case,
            "minimum": float(min(candidate_recall_by_case.values())),
            "mean": float(np.mean(list(candidate_recall_by_case.values()))),
        },
    )

    model_records: dict[str, dict[str, object]] = {}
    per_case_rows: list[dict[str, object]] = []
    threshold_records: dict[str, object] = {}
    for model_name in MODEL_FEATURE_SETS:
        record, per_case, threshold_record = run_model(
            model_name,
            cache,
            bootstrap_replicates=bootstrap_replicates,
        )
        model_records[model_name] = record
        threshold_records[model_name] = threshold_record
        for row in per_case:
            per_case_rows.append({"model": model_name, **row})
    _write_json(output_dir / "ablation_metrics.json", model_records)
    _write_json(
        output_dir / "threshold_protocol.json",
        {
            "tuning_split": "validation",
            "locked_test_used_for_threshold_selection": False,
            "rule": "maximum recall subject to FDR <= 0.25, with low-FDR fallback",
            "thresholds_strictly_interior": True,
            "models": threshold_records,
        },
    )
    _write_comparison(output_dir / "comparison_table.csv", model_records)
    curve_rows: list[dict[str, object]] = []
    reliability_rows: list[dict[str, object]] = []
    for model_name, record in model_records.items():
        for row in record["locked_test"]["probability_curve"]:
            curve_rows.append({"model": model_name, "split": "locked_test", **row})
        for row in record["locked_test"]["metrics"]["reliability"]:
            reliability_rows.append({"model": model_name, "split": "locked_test", **row})
    _write_csv(
        output_dir / "probability_curves.csv",
        curve_rows,
        tuple(curve_rows[0]) if curve_rows else ("model", "split", "threshold"),
    )
    _write_csv(
        output_dir / "reliability_bins.csv",
        reliability_rows,
        tuple(reliability_rows[0]) if reliability_rows else ("model", "split", "bin_lower"),
    )
    _write_csv(
        output_dir / "per_case_metrics.csv",
        per_case_rows,
        tuple(per_case_rows[0]) if per_case_rows else ("model", "case_id"),
    )
    (output_dir / "interpretation.md").write_text(
        _interpretation(model_records),
        encoding="utf-8",
    )
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m7-root", type=Path, default=DEFAULT_M7_ROOT)
    parser.add_argument("--savage-root", type=Path, default=DEFAULT_SAVAGE_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--bootstrap-replicates", type=int, default=500)
    parser.add_argument("--skip-savage", action="store_true")
    args = parser.parse_args()
    if args.bootstrap_replicates < 10:
        raise SystemExit("--bootstrap-replicates must be at least 10")
    if not args.m7_root.exists():
        raise SystemExit(f"M7 root not found: {args.m7_root}")
    if not args.skip_savage and not args.savage_root.exists():
        raise SystemExit(f"Savage root not found: {args.savage_root}")
    if args.output_dir is None:
        run_id = datetime.now(timezone.utc).strftime("run-%Y%m%dT%H%M%SZ")
        args.output_dir = DEFAULT_OUTPUT_ROOT / run_id
    output = run_audit(
        output_dir=args.output_dir,
        m7_root=args.m7_root,
        savage_root=args.savage_root,
        bootstrap_replicates=args.bootstrap_replicates,
        skip_savage=args.skip_savage,
    )
    print(f"topology-v2 audit written to {output}")


if __name__ == "__main__":
    main()
