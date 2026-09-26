"""Run the truth-blind M4-C path-aware Savage topology experiment.

M4-C is intentionally separate from the frozen M4-A/M4-B directories.  It
uses public physical model inputs only during inference: the archived FHD
head field, the structured MODFLOW grid, CBC face-flow magnitudes, and CBC
well extraction context.  MODPATH reference edges are loaded only after all
path features and scores have been constructed.

The M4-C score is a predeclared physical path-support score, not a calibrated
probability. It combines source-relative beam-path support, mean local
transition efficiency, and a public-CBC terminal-sink indicator. The fixed
endpoint thresholds are recorded for diagnostic tri-state reporting and are
never selected from Savage labels.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from typing import Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.physics.modflow_head import (  # noqa: E402
    build_grid_geometry_from_params,
    parse_fhd,
)
from hydrosheaf.validation.m4_path_aware import (  # noqa: E402
    PathAwareConfig,
    PathAwareResult,
    build_path_aware_features,
)
from hydrosheaf.validation.topology_v2 import (  # noqa: E402
    TopologyV2CandidateEdge,
    TopologyV2CandidateUniverse,
    TopologyV2Config,
    build_topology_v2_feature_rows,
    generate_topology_v2_candidate_universe,
)
from run_m4_fair_topology_benchmark import (  # noqa: E402
    CBC_PATH,
    FHD_PATH,
    NODE_PATH,
    REFERENCE_PATH,
    _load_savage_inputs,
)


RUN_VERSION = "m4_c_path_aware_v2_sink_aware"
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT
    / "M4"
    / "m4_topology_benchmark"
    / "results"
    / "m4_c_path_aware_20260922"
)
M4_C_PRESENT_THRESHOLD = 0.50
M4_C_ABSENT_THRESHOLD = 0.10
M4_C_RELATIVE_WEIGHT = 0.45
M4_C_EFFICIENCY_WEIGHT = 0.25
M4_C_ENDPOINT_SINK_WEIGHT = 0.30


def _sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: object) -> object:
    return json.loads(json.dumps(value, default=str, sort_keys=True))


def _hash_payload(value: object) -> str:
    payload = json.dumps(_jsonable(value), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _labels(edge_ids: list[str], reference_edges: set[str]) -> np.ndarray:
    return np.asarray([float(edge_id in reference_edges) for edge_id in edge_ids], dtype=float)


def _path_score(path_attrs: Mapping[str, object]) -> float:
    if not bool(path_attrs.get("path_reached", False)):
        return 0.0
    relative = float(path_attrs.get("path_relative_support", 0.0) or 0.0)
    efficiency = float(path_attrs.get("path_transition_efficiency", 0.0) or 0.0)
    endpoint_sink = float(path_attrs.get("path_endpoint_sink_support", 0.0) or 0.0)
    score = (
        M4_C_RELATIVE_WEIGHT * relative
        + M4_C_EFFICIENCY_WEIGHT * efficiency
        + M4_C_ENDPOINT_SINK_WEIGHT * endpoint_sink
    )
    return float(max(0.0, min(1.0, score)))


def _attach_path_features(
    universe: TopologyV2CandidateUniverse,
    path_result: PathAwareResult,
) -> TopologyV2CandidateUniverse:
    edges: list[TopologyV2CandidateEdge] = []
    for edge in universe.edges:
        attrs = dict(edge.attrs)
        attrs.update(dict(path_result.edge_features.get(edge.edge_id, {})))
        attrs["m4_c_path_support_score"] = _path_score(attrs)
        edges.append(TopologyV2CandidateEdge(edge.u, edge.v, attrs))
    edges.sort(key=lambda edge: edge.edge_id)
    return replace(
        universe,
        edges=tuple(edges),
        algorithm="m4_c_truth_blind_all_pairs_path_aware",
        parameters={
            **dict(universe.parameters),
            "path_aware_features_attached": True,
            "path_reference_edges_used": False,
            "path_score_formula": (
                "0.45 * path_relative_support + 0.25 * path_transition_efficiency + "
                "0.30 * path_endpoint_sink_support"
            ),
        },
        candidate_hash=_hash_payload([edge.to_dict() for edge in edges]),
    )


def _fixed_decisions(scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    present = scores >= M4_C_PRESENT_THRESHOLD
    absent = scores <= M4_C_ABSENT_THRESHOLD
    abstain = ~(present | absent)
    return present, absent, abstain


def _confusion(scores: np.ndarray, labels: np.ndarray) -> dict[str, object]:
    present, absent, abstain = _fixed_decisions(scores)
    tp = int(np.sum(present & (labels == 1.0)))
    fp = int(np.sum(present & (labels == 0.0)))
    fn = int(np.sum(absent & (labels == 1.0)))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "n_present": int(np.sum(present)),
        "n_absent": int(np.sum(absent)),
        "n_abstain": int(np.sum(abstain)),
    }


def _load_reference_edges_for_evaluation() -> set[str]:
    """Read MODPATH labels only at the post-inference evaluation boundary."""

    reference = pd.read_csv(REFERENCE_PATH)
    return {
        f"{str(row['u'])}->{str(row['v'])}"
        for _, row in reference.iterrows()
    }


def _edge_table(
    universe: TopologyV2CandidateUniverse,
    feature_rows: tuple[object, ...],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for edge, feature_row in zip(universe.edges, feature_rows):
        attrs = dict(edge.attrs)
        record: dict[str, object] = {
            "edge_id": edge.edge_id,
            "u": edge.u,
            "v": edge.v,
            "m4_c_path_support_score": float(attrs.get("m4_c_path_support_score", 0.0)),
            "path_reached": bool(attrs.get("path_reached", False)),
            "path_reachability": attrs.get("path_reachability"),
            "path_relative_support": attrs.get("path_relative_support"),
            "path_log_support": attrs.get("path_log_support"),
            "path_steps": attrs.get("path_steps"),
            "path_transition_efficiency": attrs.get("path_transition_efficiency"),
            "path_activity_support": attrs.get("path_activity_support"),
            "path_endpoint_sink_support": attrs.get("path_endpoint_sink_support"),
            "distance_m": attrs.get("distance_m"),
            "head_delta_m": attrs.get("head_delta_m"),
            "gradient_m_per_km": attrs.get("gradient_m_per_km"),
            "source_pumping_strength": feature_row.features.get("source_pumping_strength"),
            "target_pumping_sink_strength": feature_row.features.get("target_pumping_sink_strength"),
            "in_reference_for_evaluation_only": False,
        }
        rows.append(record)
    return pd.DataFrame(rows)


def _feature_table(feature_rows: tuple[object, ...]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for feature_row in feature_rows:
        record: dict[str, object] = {
            "edge_id": feature_row.edge_id,
            "case_id": feature_row.case_id,
            "u": feature_row.u,
            "v": feature_row.v,
            "missing_features": json.dumps(list(feature_row.missing_features)),
            "candidate_prior_probability": feature_row.candidate_prior_probability,
        }
        for feature_name, value in feature_row.features.items():
            record[f"feature_{feature_name}"] = value
        rows.append(record)
    return pd.DataFrame(rows)


def _write_plot(table: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for label, color in ((0.0, "tab:gray"), (1.0, "tab:red")):
        values = table.loc[table["label_posthoc"] == label, "m4_c_path_support_score"]
        axes[0].hist(
            values,
            bins=np.linspace(0.0, 1.0, 31),
            alpha=0.62,
            density=True,
            color=color,
            label="MODPATH-present" if label else "MODPATH-absent",
        )
    axes[0].axvline(M4_C_PRESENT_THRESHOLD, color="black", linestyle="--", linewidth=0.9)
    axes[0].axvline(M4_C_ABSENT_THRESHOLD, color="black", linestyle=":", linewidth=0.9)
    axes[0].set_xlabel("M4-C physical path-support score")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Path-score separation")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.2)

    reached = table[table["path_reached"]]
    axes[1].scatter(
        reached["distance_m"],
        reached["m4_c_path_support_score"],
        s=8,
        alpha=0.32,
        c=reached["label_posthoc"],
        cmap="coolwarm",
    )
    axes[1].set_xlabel("Candidate distance (m)")
    axes[1].set_ylabel("M4-C path-support score")
    axes[1].set_title("Path score versus endpoint distance")
    axes[1].grid(alpha=0.2)
    fig.suptitle("M4-C diagnostic; MODPATH labels shown only post hoc")
    fig.tight_layout()
    fig.savefig(output_dir / "m4_c_path_score_separation.png", dpi=220)
    plt.close(fig)


def run(output_dir: Path = DEFAULT_OUTPUT_DIR) -> Path:
    started = datetime.now(timezone.utc)
    output_dir.mkdir(parents=True, exist_ok=True)
    nodes, observations, archive_metadata = _load_savage_inputs()
    head_map = {
        int(cell): float(head)
        for cell, head in parse_fhd(FHD_PATH).items()
    }
    grid = build_grid_geometry_from_params(
        ncol=183,
        nrow=202,
        nlay=8,
        dx=110.0,
        dy=73.333,
        rotation_deg=-12.0,
        origin_x=961030.4,
        origin_y=112955.0,
    )
    path_config = PathAwareConfig()

    # Everything through this point is truth-blind.  In particular, the path
    # result and all endpoint scores are constructed before the reference file
    # is read below.
    path_result = build_path_aware_features(
        observations,
        head_map=head_map,
        cbc_path=CBC_PATH,
        grid=grid,
        config=path_config,
    )
    base_universe = generate_topology_v2_candidate_universe(
        observations,
        config=TopologyV2Config(default_head_sigma_m=0.05),
    )
    universe = _attach_path_features(base_universe, path_result)
    feature_rows = build_topology_v2_feature_rows(
        universe,
        observations,
        case_id="savage",
        config=TopologyV2Config(default_head_sigma_m=0.05),
    )
    edge_ids = list(universe.edge_ids)
    scores = np.asarray(
        [float(edge.attrs.get("m4_c_path_support_score", 0.0)) for edge in universe.edges],
        dtype=float,
    )

    # Evaluation boundary: no reference edge is available to any construction
    # or scoring function above this line.
    reference_edges = _load_reference_edges_for_evaluation()
    labels = _labels(edge_ids, reference_edges)
    table = _edge_table(universe, feature_rows)
    table["label_posthoc"] = labels
    table["in_reference_for_evaluation_only"] = labels.astype(bool)
    table["present_at_predeclared_threshold"] = scores >= M4_C_PRESENT_THRESHOLD
    table["absent_at_predeclared_threshold"] = scores <= M4_C_ABSENT_THRESHOLD
    table["abstain_at_predeclared_threshold"] = ~(
        table["present_at_predeclared_threshold"]
        | table["absent_at_predeclared_threshold"]
    )

    confusion = _confusion(scores, labels)
    prevalence = float(np.mean(labels)) if len(labels) else float("nan")
    summary = {
        "run_version": RUN_VERSION,
        "n_nodes": len(universe.nodes),
        "n_candidate_edges": len(universe.edges),
        "n_reference_edges": len(reference_edges),
        "candidate_graph_recall_evaluation_only": universe.candidate_graph_recall(reference_edges),
        "path_reached_fraction": float(np.mean(table["path_reached"].to_numpy(dtype=bool))),
        "prevalence": prevalence,
        "path_score_pr_auc": float(average_precision_score(labels, scores)),
        "path_score_roc_auc": float(roc_auc_score(labels, scores)),
        "path_score_brier_as_bounded_score": float(brier_score_loss(labels, scores)),
        "mean_path_score": float(np.mean(scores)),
        "positive_mean_path_score": float(np.mean(scores[labels == 1.0])),
        "negative_mean_path_score": float(np.mean(scores[labels == 0.0])),
        "score_separation_mean_difference": float(
            np.mean(scores[labels == 1.0]) - np.mean(scores[labels == 0.0])
        ),
        "present_threshold_predeclared": M4_C_PRESENT_THRESHOLD,
        "absent_threshold_predeclared": M4_C_ABSENT_THRESHOLD,
        **confusion,
        "truth_blind_inference": True,
        "reference_edges_used_only_after_inference": True,
        "thresholds_retuned_on_savage": False,
        "calibration_transfer_established": False,
        "result_status": "locked_physical_path_support_experiment",
        "allowed_claim": (
            "Exploratory model-conditioned topology recovery using public FHD/CBC physical evidence; "
            "the MODPATH graph is a withheld model-conditioned reference, not field truth."
        ),
    }

    table.to_csv(output_dir / "m4_c_edge_scores.csv", index=False)
    _feature_table(feature_rows).to_csv(output_dir / "m4_c_feature_rows.csv", index=False)
    pd.DataFrame(list(path_result.source_diagnostics)).to_csv(
        output_dir / "m4_c_source_path_diagnostics.csv", index=False
    )
    pd.DataFrame([summary]).to_csv(output_dir / "m4_c_summary.csv", index=False)
    _write_plot(table, output_dir)

    source_paths = {
        "node_mapping": str(NODE_PATH),
        "base_fhd": str(FHD_PATH),
        "base_cbc": str(CBC_PATH),
        "reference_edges_evaluator_only": str(REFERENCE_PATH),
    }
    manifest = {
        "run_version": RUN_VERSION,
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_paths": source_paths,
        "source_sha256": {key: _sha256_file(Path(value)) for key, value in source_paths.items()},
        "code_paths": {
            "runner": str(Path(__file__).resolve()),
            "path_aware_module": str(
                PROJECT_ROOT / "hydrosheaf" / "validation" / "m4_path_aware.py"
            ),
            "topology_v2_module": str(
                PROJECT_ROOT / "hydrosheaf" / "validation" / "topology_v2.py"
            ),
        },
        "code_sha256": {
            "runner": _sha256_file(Path(__file__).resolve()),
            "path_aware_module": _sha256_file(
                PROJECT_ROOT / "hydrosheaf" / "validation" / "m4_path_aware.py"
            ),
            "topology_v2_module": _sha256_file(
                PROJECT_ROOT / "hydrosheaf" / "validation" / "topology_v2.py"
            ),
        },
        "artifacts": [
            "m4_c_summary.csv",
            "m4_c_edge_scores.csv",
            "m4_c_feature_rows.csv",
            "m4_c_source_path_diagnostics.csv",
            "m4_c_path_score_separation.png",
            "m4_c_manifest.json",
            "README.md",
        ],
        "archive_metadata": archive_metadata,
        "grid": {
            "ncol": grid.ncol,
            "nrow": grid.nrow,
            "nlay": grid.nlay,
            "dx": grid.dx,
            "dy": grid.dy,
            "rotation_deg": grid.rotation_deg,
            "origin_x": grid.origin_x,
            "origin_y": grid.origin_y,
        },
        "path_aware_config": path_config.to_dict(),
        "path_result": path_result.to_dict(),
        "score_definition": {
            "name": "predeclared_physical_path_support",
            "formula": (
                "0.45 * path_relative_support + 0.25 * path_transition_efficiency + "
                "0.30 * path_endpoint_sink_support"
            ),
            "unreached_path_score": 0.0,
            "present_threshold": M4_C_PRESENT_THRESHOLD,
            "absent_threshold": M4_C_ABSENT_THRESHOLD,
            "threshold_source": "predeclared before Savage evaluation; not tuned on reference labels",
        },
        "candidate_contract": {
            "candidate_generation": "all directed observation pairs",
            "path_features_attach_after_candidate_generation": True,
            "path_features_remove_candidates": False,
            "reference_edges_in_candidate_generation": False,
            "reference_edges_in_path_tracing": False,
            "reference_edges_in_score_formula": False,
            "reference_edges_in_threshold_selection": False,
            "reference_edges_used_only_for_posthoc_metrics": True,
        },
        "limitations": [
            "CBC face-flow magnitudes are used as non-directional activity; their signs are not treated as particle direction.",
            "Vertical transitions use a unit layer-index proxy because no public physical layer thickness field was supplied to this run.",
            "The path score is bounded support, not a calibrated posterior probability.",
            "Calibration transfer to Savage remains unestablished.",
            "MODPATH is a model-conditioned reference graph, not independent field truth.",
        ],
        "truth_blind_inference": True,
        "legacy_outputs_modified": False,
    }
    (output_dir / "m4_c_manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2), encoding="utf-8"
    )
    readme = """# M4-C path-aware HydroSheaf

M4-C is a separate, truth-blind diagnostic experiment. It uses the public
MODFLOW head field, head-gradient direction, structured-grid adjacency, CBC
face-flow magnitudes as non-directional activity, and CBC well extraction
context. All directed endpoint pairs remain in the candidate universe; path
features do not pre-filter candidates.

The score is a predeclared bounded physical path-support score, not a
calibrated probability. It includes a soft terminal preference for public
CBC sink cells. The thresholds in `m4_c_manifest.json` were fixed before the
Savage reference file was read. MODPATH labels are present only in the edge
table as post-hoc evaluation annotations.

The ranking metrics in `m4_c_summary.csv` are the main M4-C diagnostic. The
tri-state counts use the conservative predeclared score thresholds and are
not threshold-tuned on Savage; calibration transfer is therefore still an
open issue even when ranking improves.

This remains an exploratory model-conditioned benchmark. A high score means
that the public physical inputs support a grid path under the declared beam
and transition model; it does not establish field-truth recovery.
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    output = run(args.output_dir)
    summary = pd.read_csv(output / "m4_c_summary.csv").iloc[0]
    print(f"Wrote M4-C path-aware outputs to {output}")
    print(
        f"ROC-AUC={summary['path_score_roc_auc']:.3f} "
        f"PR-AUC={summary['path_score_pr_auc']:.3f} "
        f"F1={summary['f1']:.3f} "
        f"candidate_recall={summary['candidate_graph_recall_evaluation_only']:.3f}"
    )


if __name__ == "__main__":
    main()
