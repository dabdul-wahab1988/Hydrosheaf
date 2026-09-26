"""Run the redesigned, archive-informed M4 topology benchmark.

The runner has an intentionally narrow contract:

* M4-A uses projected geometry and the sparse endpoint-z proxy only.
* M4-B uses the same all-pairs candidate universe plus actual FHD heads and
  public MODFLOW CBC source/sink records.
* Synthetic development cases fit the probability calibrator and freeze the
  tri-state thresholds before Savage is scored.
* MODPATH reference edges are loaded only by the evaluator after inference.

The legacy Phase 2b outputs are not modified.  Results are written to a new
``fair_redesign_*`` directory under the benchmark results tree.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable, Mapping, Sequence
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.physics.modflow_head import parse_fhd
from hydrosheaf.validation.m4_fair import (
    M4_A,
    M4_B,
    SavageProjectedFrame,
    aggregate_savage_cbc_context,
    build_savage_observations,
)
from hydrosheaf.validation.topology_v2 import (
    MODEL_FEATURE_SETS,
    TopologyV2Config,
    TopologyV2LogisticCalibrator,
    TopologyV2Scorer,
    build_topology_v2_feature_rows,
    generate_topology_v2_candidate_universe,
    topology_v2_metrics,
    tune_topology_thresholds,
)


RUN_VERSION = "m4_fair_redesign_v1"
RANDOM_SEED = 20260922
SAVAGE_RESULTS = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "public_archives" / "savage"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "fair_redesign_20260922"
NODE_PATH = SAVAGE_RESULTS / "modpath_node_mapping.csv"
REFERENCE_PATH = SAVAGE_RESULTS / "modpath_reference_edges.csv"
FHD_PATH = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "public_archives" / "savage" / "base.fhd"
CBC_PATH = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "public_archives" / "savage" / "output" / "output.D2010_base_calibration" / "base.cbc"
RAW_ARCHIVE_PATH = PROJECT_ROOT / "M4" / "data" / "Teir_1" / "output.zip"


def _jsonable(value: object) -> object:
    return json.loads(json.dumps(value, default=str, sort_keys=True))


def _sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_payload(value: object) -> str:
    payload = json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _edge_set(edges: Iterable[object]) -> set[str]:
    result: set[str] = set()
    for edge in edges:
        if isinstance(edge, str) and "->" in edge:
            left, right = edge.split("->", 1)
        elif isinstance(edge, Mapping):
            left = edge.get("u", edge.get("source"))
            right = edge.get("v", edge.get("target"))
            if left is None or right is None:
                continue
        elif isinstance(edge, Sequence) and len(edge) >= 2:
            left, right = edge[0], edge[1]
        else:
            left = getattr(edge, "u", None)
            right = getattr(edge, "v", None)
        if left is not None and right is not None and str(left) != str(right):
            result.add(f"{str(left)}->{str(right)}")
    return result


def _direct_confusion(predicted: set[str], reference: set[str]) -> dict[str, object]:
    true_positive = len(predicted & reference)
    false_positive = len(predicted - reference)
    false_negative = len(reference - predicted)
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else 0.0
    recall = true_positive / len(reference) if reference else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "tp": true_positive,
        "fp": false_positive,
        "fn": false_negative,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def _labels(rows: Sequence[Any], reference: set[str]) -> np.ndarray:
    return np.asarray([float(row.edge_id in reference) for row in rows], dtype=float)


def _synthetic_case(mode: str, case_name: str, *, missing_elevation: bool = False) -> tuple[tuple[dict[str, object], ...], set[str]]:
    """Build a declared development case, independent of Savage labels."""

    case_prefix = f"{mode.replace('-', '').lower()}_{case_name}"
    rows: list[dict[str, object]] = []
    n_nodes = 12
    for index in range(n_nodes):
        column = index % 4
        row_index = index // 4
        x = float(column * 100.0 + (index % 2) * 3.0)
        y = float(row_index * 120.0 + ((index // 2) % 2) * 2.0)
        potential = 200.0 - 0.75 * x - 0.55 * y
        observation: dict[str, object] = {
            "site_id": f"{case_prefix}_n{index:02d}",
            "x_m": x,
            "y_m": y,
        }
        if not (missing_elevation and index == 2):
            observation["elevation"] = potential
        if mode == M4_B:
            observation.update(
                {
                    "hydraulic_head": potential,
                    "head_sigma_m": 0.05,
                    "well_rate": -500.0 if index == n_nodes - 1 else 0.0,
                    "river_leakage": 0.0,
                    "recharge": 100.0 if index == 0 else 0.0,
                    "head_boundary_flux": 0.0,
                }
            )
        rows.append(observation)

    ids = [str(row["site_id"]) for row in rows]
    reference: set[str] = set()
    # A simple directed chain plus additional fan-in edges.  The labels are
    # kept in this function only for development calibration and are never
    # passed to candidate generation or feature construction.
    for index in range(n_nodes - 1):
        reference.add(f"{ids[index]}->{ids[index + 1]}")
    for index in (7, 8, 9):
        reference.add(f"{ids[index]}->{ids[-1]}")
    return tuple(rows), reference


def _fit_development_calibration(mode: str) -> dict[str, object]:
    if mode == M4_A:
        feature_names = MODEL_FEATURE_SETS["M4_A_sparse"]
    elif mode == M4_B:
        feature_names = MODEL_FEATURE_SETS["M4_B_archive_informed"]
    else:
        raise ValueError(mode)

    fit_rows: list[Any] = []
    fit_labels: list[float] = []
    fit_payload: list[object] = []
    for case_name in ("fit_1", "fit_2"):
        observations, reference = _synthetic_case(
            mode,
            case_name,
            missing_elevation=(mode == M4_A and case_name == "fit_2"),
        )
        universe = generate_topology_v2_candidate_universe(
            observations,
            config=TopologyV2Config(default_head_sigma_m=0.05),
        )
        rows = build_topology_v2_feature_rows(universe, observations, case_id=case_name)
        fit_rows.extend(rows)
        fit_labels.extend(_labels(rows, reference).tolist())
        fit_payload.append({"case_id": case_name, "observations": observations, "reference": sorted(reference)})

    dataset_hash = _sha256_payload(fit_payload)
    calibrator = TopologyV2LogisticCalibrator(feature_names=feature_names, l2=0.25)
    calibrator.fit(
        fit_rows,
        fit_labels,
        scope="held_out_calibration",
        independent=True,
        generator_id=f"{RUN_VERSION}:{mode}",
        split_id=f"{mode}:synthetic-fit-cases-v1",
        dataset_hash=dataset_hash,
    )

    threshold_observations, threshold_reference = _synthetic_case(
        mode,
        "threshold_1",
        missing_elevation=(mode == M4_A),
    )
    threshold_universe = generate_topology_v2_candidate_universe(
        threshold_observations,
        config=TopologyV2Config(default_head_sigma_m=0.05),
    )
    threshold_rows = build_topology_v2_feature_rows(
        threshold_universe,
        threshold_observations,
        case_id="threshold_1",
    )
    threshold_labels = _labels(threshold_rows, threshold_reference)
    threshold_probabilities = calibrator.predict_proba(threshold_rows)
    policy, threshold_audit = tune_topology_thresholds(
        threshold_labels,
        threshold_probabilities,
        target_fdr=0.25,
        minimum_gap=0.10,
    )
    return {
        "mode": mode,
        "feature_names": tuple(feature_names),
        "calibrator": calibrator,
        "policy": policy,
        "threshold_audit": threshold_audit,
        "calibration_metadata": {
            "calibration_source": "synthetic_development_cases_v1",
            "fit_case_ids": ["fit_1", "fit_2"],
            "threshold_case_id": "threshold_1",
            "fit_labels_are_not_savage_reference": True,
            "dataset_hash": dataset_hash,
            "split_id": f"{mode}:synthetic-fit-cases-v1",
            "threshold_target_fdr": 0.25,
            "threshold_minimum_gap": 0.10,
        },
    }


def _infer_without_reference(
    observations: Sequence[Mapping[str, object]],
    *,
    mode: str,
    calibration: Mapping[str, object],
) -> dict[str, object]:
    """Run candidate generation and scoring without accepting truth input."""

    universe = generate_topology_v2_candidate_universe(
        observations,
        config=TopologyV2Config(default_head_sigma_m=0.05),
    )
    rows = build_topology_v2_feature_rows(
        universe,
        observations,
        case_id="savage",
        config=TopologyV2Config(default_head_sigma_m=0.05),
    )
    calibrator = calibration["calibrator"]
    policy = calibration["policy"]
    assert isinstance(calibrator, TopologyV2LogisticCalibrator)
    scorer = TopologyV2Scorer(calibrator, policy=policy)
    probabilities = calibrator.predict_proba(rows)
    decisions = scorer.score_rows(rows)
    if len(decisions) != len(rows):
        raise RuntimeError("Topology-v2 scorer returned a misaligned result.")
    return {
        "mode": mode,
        "universe": universe,
        "rows": rows,
        "probabilities": probabilities,
        "decisions": decisions,
        "calibrator": calibrator,
        "policy": policy,
        "truth_blind_inference": True,
        "reference_edges_used": False,
    }


def _evaluate_after_inference(
    inference: Mapping[str, object],
    reference_edges: set[str],
) -> tuple[dict[str, object], pd.DataFrame]:
    universe = inference["universe"]
    rows = inference["rows"]
    probabilities = np.asarray(inference["probabilities"], dtype=float)
    decisions = inference["decisions"]
    assert hasattr(universe, "candidate_graph_recall")
    labels = _labels(rows, reference_edges)
    policy = inference["policy"]
    metrics = topology_v2_metrics(
        labels,
        probabilities,
        policy=policy,
        candidate_recall=universe.candidate_graph_recall(reference_edges),
    )
    predicted = {
        str(record["edge_id"])
        for record in decisions
        if record.get("decision") == "PRESENT"
    }
    confusion = _direct_confusion(predicted, reference_edges)
    summary = {
        "mode": inference["mode"],
        "n_nodes": len(universe.nodes),
        "n_candidate_edges": len(universe.edges),
        "n_reference_edges": len(reference_edges),
        "candidate_graph_recall": universe.candidate_graph_recall(reference_edges),
        "n_present": int(sum(record.get("decision") == "PRESENT" for record in decisions)),
        "n_absent": int(sum(record.get("decision") == "ABSENT" for record in decisions)),
        "n_abstain": int(sum(record.get("decision") == "ABSTAIN" for record in decisions)),
        "probability_pr_auc": metrics["pr_auc"],
        "probability_roc_auc": metrics["roc_auc"],
        "brier": metrics["brier"],
        "probability_log_loss": metrics["log_loss"],
        "probability_ece": metrics["ece"],
        "selective_precision": metrics["selective_precision"],
        "selective_recall": metrics["selective_recall"],
        "selective_f1": metrics["f1"],
        **confusion,
        "truth_blind_inference": True,
        "reference_edges_used_only_after_inference": True,
        "calibration_transfer_established": False,
        "result_status": "locked_rule_exploratory",
        "allowed_claim": (
            "Exploratory archive-conditioned reduced-order topology recovery under a "
            "synthetic-development locked rule; calibration transfer is not established, "
            "and MODPATH is a model-conditioned reference, not field truth."
        ),
    }

    score_rows: list[dict[str, object]] = []
    for row, probability, decision in zip(rows, probabilities, decisions):
        record = row.to_dict()
        record.pop("features", None)
        record.update(
            {
                "probability": float(probability),
                "decision": decision.get("decision"),
                "missing_features": json.dumps(decision.get("missing_features", [])),
                "in_reference_for_evaluation_only": row.edge_id in reference_edges,
            }
        )
        for feature_name, value in row.features.items():
            record[f"feature_{feature_name}"] = value
        score_rows.append(record)
    return summary, pd.DataFrame(score_rows)


def _load_savage_inputs() -> tuple[pd.DataFrame, tuple[dict[str, object], ...], dict[str, object]]:
    if not NODE_PATH.exists():
        raise FileNotFoundError(NODE_PATH)
    if not FHD_PATH.exists():
        raise FileNotFoundError(FHD_PATH)
    if not CBC_PATH.exists():
        raise FileNotFoundError(CBC_PATH)

    nodes = pd.read_csv(NODE_PATH)
    # parse_fhd returns integer one-based MODFLOW cell IDs; the public node
    # mapping names the same cells as ``cell_<id>``.
    head_map = {
        f"cell_{int(key)}": float(value)
        for key, value in parse_fhd(FHD_PATH).items()
    }
    node_ids = [str(value) for value in nodes["node_id"]]
    budget_context, budget_metadata = aggregate_savage_cbc_context(CBC_PATH, node_ids)
    frame = SavageProjectedFrame()
    # The return value is deliberately only observations and source metadata;
    # the reference set is held by the caller for post-inference evaluation.
    observations = build_savage_observations(
        nodes,
        mode=M4_B,
        heads=head_map,
        budget_context=budget_context,
        frame=frame,
    )
    metadata = {
        "frame": frame.to_dict(),
        "budget": budget_metadata,
        "head_count": len(head_map),
        "raw_archive_zip_exists": RAW_ARCHIVE_PATH.exists(),
    }
    return nodes, observations, metadata


def _load_reference_edges_for_evaluation() -> set[str]:
    """Load MODPATH edges only after both truth-blind inference runs finish."""

    if not REFERENCE_PATH.exists():
        raise FileNotFoundError(REFERENCE_PATH)
    reference_frame = pd.read_csv(REFERENCE_PATH)
    return {
        f"{str(row['u'])}->{str(row['v'])}"
        for _, row in reference_frame.iterrows()
    }


def _write_outputs(
    output_dir: Path,
    *,
    summaries: Sequence[Mapping[str, object]],
    score_tables: Mapping[str, pd.DataFrame],
    calibrations: Mapping[str, Mapping[str, object]],
    candidate_audits: Sequence[Mapping[str, object]],
    node_context: pd.DataFrame,
    manifest: Mapping[str, object],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(list(summaries)).to_csv(output_dir / "fair_benchmark_summary.csv", index=False)
    pd.DataFrame(list(candidate_audits)).to_csv(output_dir / "fair_candidate_audit.csv", index=False)
    node_context.to_csv(output_dir / "fair_archive_context.csv", index=False)
    for mode, table in score_tables.items():
        safe = "m4_a" if mode == M4_A else "m4_b"
        table.to_csv(output_dir / f"fair_edge_scores_{safe}.csv", index=False)

    calibration_rows: list[dict[str, object]] = []
    for mode, payload in calibrations.items():
        calibrator = payload["calibrator"]
        policy = payload["policy"]
        safe = "m4_a" if mode == M4_A else "m4_b"
        with (output_dir / f"fair_calibrator_{safe}.json").open("w", encoding="utf-8") as handle:
            json.dump(calibrator.to_dict(), handle, indent=2, default=str)
        audit = pd.DataFrame(list(payload["threshold_audit"]))
        audit.to_csv(output_dir / f"fair_threshold_audit_{safe}.csv", index=False)
        calibration_rows.append(
            {
                "mode": mode,
                "feature_names": json.dumps(list(payload["feature_names"])),
                "present_threshold": policy.present_threshold,
                "absent_threshold": policy.absent_threshold,
                "tuning_rule": policy.tuning_rule,
                "target_fdr": policy.target_fdr,
                "minimum_gap": policy.minimum_gap,
                "calibrator_deployment_eligible": calibrator.deployment_eligible,
                "calibration_transfer_established": False,
                **dict(payload["calibration_metadata"]),
            }
        )
    pd.DataFrame(calibration_rows).to_csv(output_dir / "fair_calibration.csv", index=False)
    with (output_dir / "fair_run_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(_jsonable(manifest), handle, indent=2)
    readme = """# M4 fair redesign run

This directory is separate from the legacy `results/public_archives/savage`
outputs.  M4-A uses projected geometry and the endpoint-z sparse proxy. M4-B
adds actual `base.fhd` heads and public MODFLOW CBC source/sink context.

Candidate generation is the deterministic all-directed-pairs topology-v2
universe. Direction, distance, and source/sink evidence are graded features;
they are not pre-inference truth gates. Calibration and threshold selection use
synthetic development cases, and the Savage MODPATH edge list is read only by
the post-inference evaluator. Therefore the reported MODPATH comparison is
still model-conditioned and must not be described as independent field truth.

Because this checkout contains only the Savage evaluation artifacts and not a
separate public development archive, calibration transfer to Savage is not
established. The scores are an auditable locked-rule experiment, not evidence
of improvement over the legacy benchmark.

The raw public `output.zip` is recorded in the manifest. If it is absent, the
run is based on the locally available derived FHD/CBC and mapping artifacts;
that is a provenance limitation, not a claim of full archive reproduction.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")


def run(output_dir: Path = DEFAULT_OUTPUT_DIR) -> Path:
    started = datetime.now(timezone.utc)
    nodes, archive_b_observations, archive_metadata = _load_savage_inputs()
    frame = SavageProjectedFrame()
    archive_a_observations = build_savage_observations(
        nodes,
        mode=M4_A,
        frame=frame,
    )

    calibrations = {
        M4_A: _fit_development_calibration(M4_A),
        M4_B: _fit_development_calibration(M4_B),
    }
    inferences = {
        M4_A: _infer_without_reference(
            archive_a_observations,
            mode=M4_A,
            calibration=calibrations[M4_A],
        ),
        M4_B: _infer_without_reference(
            archive_b_observations,
            mode=M4_B,
            calibration=calibrations[M4_B],
        ),
    }

    # This is intentionally after inference construction and scoring.  The
    # evaluator is the first component allowed to read MODPATH edge labels.
    reference_edges = _load_reference_edges_for_evaluation()

    summaries: list[dict[str, object]] = []
    score_tables: dict[str, pd.DataFrame] = {}
    candidate_audits: list[dict[str, object]] = []
    for mode, inference in inferences.items():
        summary, scores = _evaluate_after_inference(inference, reference_edges)
        summary["run_version"] = RUN_VERSION
        summary["calibration_source"] = calibrations[mode]["calibration_metadata"]["calibration_source"]
        summary["coordinate_system"] = frame.coordinate_system
        summaries.append(summary)
        score_tables[mode] = scores
        universe = inference["universe"]
        candidate_audits.append(
            {
                "mode": mode,
                "candidate_algorithm": universe.algorithm,
                "candidate_version": universe.version,
                "truth_blind": universe.truth_blind,
                "candidate_input_hash": universe.input_hash,
                "candidate_hash": universe.candidate_hash,
                "n_nodes": len(universe.nodes),
                "n_candidate_edges": len(universe.edges),
                "n_reference_edges_evaluation_only": len(reference_edges),
                "candidate_graph_recall_evaluation_only": universe.candidate_graph_recall(reference_edges),
                "rejection_counts": json.dumps(dict(universe.rejection_counts)),
                "reference_edges_used_during_inference": False,
            }
        )

    context_rows: list[dict[str, object]] = []
    for observation in archive_b_observations:
        context_rows.append(
            {
                "site_id": observation["site_id"],
                "well_rate": observation.get("well_rate", 0.0),
                "river_leakage": observation.get("river_leakage", 0.0),
                "recharge": observation.get("recharge", 0.0),
                "head_boundary_flux": observation.get("head_boundary_flux", 0.0),
                "hydraulic_head_m": observation.get("hydraulic_head"),
            }
        )

    source_paths = {
        "node_mapping": str(NODE_PATH),
        "reference_edges_evaluator_only": str(REFERENCE_PATH),
        "base_fhd": str(FHD_PATH),
        "base_cbc": str(CBC_PATH),
        "raw_output_zip": str(RAW_ARCHIVE_PATH),
    }
    manifest = {
        "run_version": RUN_VERSION,
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "random_seed": RANDOM_SEED,
        "modes": [M4_A, M4_B],
        "source_paths": source_paths,
        "source_sha256": {key: _sha256_file(Path(value)) for key, value in source_paths.items()},
        "raw_archive_zip_available": RAW_ARCHIVE_PATH.exists(),
        "archive_metadata": archive_metadata,
        "coordinate_frame": frame.to_dict(),
        "candidate_contract": {
            "algorithm": "truth_blind_all_pairs_soft_admissibility_v2",
            "all_directed_pairs_default": True,
            "distance_method": "Euclidean projected distance in metres",
            "uphill_or_downhill_hard_gate": False,
            "reference_edges_in_candidate_generation": False,
            "reference_edges_in_feature_construction": False,
            "reference_edges_in_threshold_selection": False,
        },
        "calibration_contract": {
            "fit_cases": "synthetic development cases",
            "threshold_case": "separate synthetic development case",
            "savage_evaluation_after_threshold_freeze": True,
            "calibration_transfer_to_savage_established": False,
        },
        "legacy_outputs_modified": False,
        "limitations": [
            "Raw public output.zip is absent in this checkout." if not RAW_ARCHIVE_PATH.exists() else "Raw public output.zip was present.",
            "No public screen/layer geometry or hydraulic conductivity fields were supplied to this run.",
            "Calibration transfer from synthetic development cases to Savage is not established.",
            "MODPATH remains a model-conditioned reference graph, not independent field truth.",
        ],
    }
    _write_outputs(
        output_dir,
        summaries=summaries,
        score_tables=score_tables,
        calibrations=calibrations,
        candidate_audits=candidate_audits,
        node_context=pd.DataFrame(context_rows),
        manifest=manifest,
    )
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Separate output directory; legacy Savage results are never overwritten.",
    )
    args = parser.parse_args()
    output_dir = run(args.output_dir)
    summary = pd.read_csv(output_dir / "fair_benchmark_summary.csv")
    print(f"Wrote fair M4 redesign outputs to {output_dir}")
    for _, row in summary.iterrows():
        print(
            f"{row['mode']}: F1={row['f1']:.3f} "
            f"P={row['precision']:.3f} R={row['recall']:.3f} "
            f"candidate_recall={row['candidate_graph_recall']:.3f}"
        )


if __name__ == "__main__":
    main()
