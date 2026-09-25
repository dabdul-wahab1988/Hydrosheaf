"""Freeze the current fair-redesign M4-A/M4-B diagnostic experiments.

This script copies the already-generated fair-redesign artifacts into a
versioned diagnostic snapshot and records hashes for the snapshot, source
inputs, and generating code.  It intentionally refuses to alter the legacy
Savage result directory and does not recompute or tune any model parameter.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SOURCE_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "fair_redesign_20260922"
DEFAULT_TARGET_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "frozen_m4_diagnostics_20260922"
LEGACY_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "public_archives" / "savage"

SOURCE_FILES = (
    "fair_benchmark_summary.csv",
    "fair_calibration.csv",
    "fair_candidate_audit.csv",
    "fair_run_manifest.json",
    "fair_archive_context.csv",
    "fair_edge_scores_m4_a.csv",
    "fair_edge_scores_m4_b.csv",
    "fair_calibrator_m4_a.json",
    "fair_calibrator_m4_b.json",
    "fair_threshold_audit_m4_a.csv",
    "fair_threshold_audit_m4_b.csv",
    "README.md",
)

CODE_FILES = (
    PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "scripts" / "run_m4_fair_topology_benchmark.py",
    PROJECT_ROOT / "hydrosheaf" / "validation" / "topology_v2.py",
    PROJECT_ROOT / "hydrosheaf" / "validation" / "m4_fair.py",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def _validate_source() -> tuple[dict[str, Any], dict[str, Any]]:
    if not SOURCE_DIR.exists():
        raise FileNotFoundError(f"Fair-redesign directory is missing: {SOURCE_DIR}")
    missing = [name for name in SOURCE_FILES if not (SOURCE_DIR / name).exists()]
    if missing:
        raise FileNotFoundError(f"Fair-redesign artifacts are missing: {missing}")
    manifest = _load_json(SOURCE_DIR / "fair_run_manifest.json")
    summary = (SOURCE_DIR / "fair_benchmark_summary.csv").read_text(encoding="utf-8")
    required_markers = (
        "M4-A_sparse_observation",
        "M4-B_archive_informed",
        "candidate_graph_recall",
        "calibration_transfer_established",
    )
    if not all(marker in summary for marker in required_markers):
        raise ValueError("Current fair summary does not contain the expected M4-A/M4-B contract.")
    if manifest.get("legacy_outputs_modified") is not False:
        raise ValueError("Refusing to freeze a run that claims legacy outputs were modified.")
    contract = manifest.get("candidate_contract", {})
    if contract.get("reference_edges_in_candidate_generation") is not False:
        raise ValueError("Fair candidate-generation contract is not truth-blind.")
    if contract.get("reference_edges_in_feature_construction") is not False:
        raise ValueError("Fair feature-construction contract is not truth-blind.")
    if contract.get("reference_edges_in_threshold_selection") is not False:
        raise ValueError("Fair threshold-selection contract is not truth-blind.")
    return manifest, {"summary_sha256": _sha256(SOURCE_DIR / "fair_benchmark_summary.csv")}


def freeze(target_dir: Path = DEFAULT_TARGET_DIR, *, force: bool = False) -> Path:
    manifest, checks = _validate_source()
    if target_dir.exists():
        existing = target_dir / "freeze_manifest.json"
        if not force and existing.exists():
            previous = _load_json(existing)
            if previous.get("source_artifact_hashes") == {
                name: _sha256(SOURCE_DIR / name) for name in SOURCE_FILES
            }:
                return target_dir
            raise FileExistsError(
                f"Frozen directory already exists with different hashes: {target_dir}. "
                "Use --force only to replace this diagnostic snapshot."
            )
        if not force:
            raise FileExistsError(f"Frozen directory already exists: {target_dir}")
        shutil.rmtree(target_dir)

    target_dir.mkdir(parents=True, exist_ok=False)
    for mode, names in {
        "M4-A": (
            "fair_edge_scores_m4_a.csv",
            "fair_calibrator_m4_a.json",
            "fair_threshold_audit_m4_a.csv",
        ),
        "M4-B": (
            "fair_edge_scores_m4_b.csv",
            "fair_calibrator_m4_b.json",
            "fair_threshold_audit_m4_b.csv",
        ),
    }.items():
        mode_dir = target_dir / mode
        mode_dir.mkdir()
        for name in names:
            shutil.copy2(SOURCE_DIR / name, mode_dir / name)
    shared_dir = target_dir / "shared"
    shared_dir.mkdir()
    for name in (
        "fair_benchmark_summary.csv",
        "fair_calibration.csv",
        "fair_candidate_audit.csv",
        "fair_archive_context.csv",
        "fair_run_manifest.json",
        "README.md",
    ):
        shutil.copy2(SOURCE_DIR / name, shared_dir / name)

    source_artifact_hashes = {name: _sha256(SOURCE_DIR / name) for name in SOURCE_FILES}
    frozen_artifact_hashes = {
        str(path.relative_to(target_dir)): _sha256(path)
        for path in target_dir.rglob("*")
        if path.is_file()
    }
    code_hashes = {
        str(path.relative_to(PROJECT_ROOT)): _sha256(path)
        for path in CODE_FILES
        if path.exists()
    }
    freeze_manifest = {
        "freeze_version": "m4_fair_diagnostics_freeze_v1",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "source_results_directory": str(SOURCE_DIR),
        "frozen_directory": str(target_dir),
        "modes": ["M4-A", "M4-B"],
        "source_artifact_hashes": source_artifact_hashes,
        "frozen_artifact_hashes": frozen_artifact_hashes,
        "generating_code_sha256": code_hashes,
        "legacy_results_directory": str(LEGACY_DIR),
        "legacy_outputs_modified_by_freeze": False,
        "reference_edges_used_during_inference": False,
        "thresholds_changed_during_freeze": False,
        "truth_labels_used_only_for_posthoc_evaluation": True,
        "source_manifest_contract": manifest.get("candidate_contract", {}),
        "freeze_checks": checks,
        "scientific_status": (
            "M4-A and M4-B are frozen diagnostic experiments. Their probabilities, "
            "thresholds, and outputs must not be altered after Savage diagnostics."
        ),
    }
    with (target_dir / "freeze_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(freeze_manifest, handle, indent=2, sort_keys=True)
    return target_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-dir", type=Path, default=DEFAULT_TARGET_DIR)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    target = freeze(args.target_dir, force=args.force)
    print(f"Frozen M4-A/M4-B diagnostics: {target}")


if __name__ == "__main__":
    main()
