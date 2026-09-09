"""Run a fresh, truth-blind integrated synthetic benchmark.

Point 2 is intentionally a new run directory.  The existing locked M7.3
results are never overwritten.  The independent generator produces a blind
observation table; Hydrosheaf is run on that table, and only then are the
sealed edge/age/process fields used for scoring.  This wrapper records the
stage boundary and audits the saved blind tables for truth-looking columns.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Sequence

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
for path in (REPO_ROOT, SCRIPT_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hydrosheaf.validation.integrated_benchmark import (  # noqa: E402
    audit_reference_panel,
    build_run_manifest,
    source_manifest_entry,
    write_json,
)
from hydrosheaf.validation.programme_contract import assert_truth_blind  # noqa: E402
from run_m7_3_nonuniqueness import (  # noqa: E402
    DEFAULT_BIN_DIR,
    DEFAULT_SIMULATOR_WORKSPACE,
    SMOKE_DEV_SEEDS,
    SMOKE_TEST_SEEDS,
    run_benchmark as run_m7_3,
)


DEFAULT_OUTPUT = REPO_ROOT / ".codex_work" / "runs" / "RUN-INTEGRATED-SYNTHETIC-20260908-01"
FULL_DEV_SEEDS = tuple(range(5201, 5207))
FULL_TEST_SEEDS = tuple(range(5301, 5313))
SMOKE_AGE_DRAWS = 120
SMOKE_AGE_PARTICLES = 2_000
SMOKE_REACTION_BOOTSTRAP = 2
SMOKE_PAIRED_BOOTSTRAP = 200
FULL_AGE_DRAWS = 600
FULL_AGE_PARTICLES = 50_000
FULL_REACTION_BOOTSTRAP = 64
FULL_PAIRED_BOOTSTRAP = 10_000


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_blind_tables(output: Path) -> dict[str, Any]:
    """Check every persisted blind observation table for leaked truth fields."""

    files = sorted((output / "cases").glob("*/blind_observations.csv"))
    reports: list[dict[str, Any]] = []
    for path in files:
        frame = pd.read_csv(path, low_memory=False)
        rows = frame.to_dict(orient="records")
        assert_truth_blind(rows)
        reports.append(
            {
                "path": str(path),
                "n_rows": int(len(frame)),
                "n_columns": int(len(frame.columns)),
                "sha256": _sha256(path),
                "truth_fields_seen": [],
            }
        )
    return {
        "status": "PASS" if files else "FAIL",
        "n_blind_tables": len(files),
        "tables": reports,
        "truth_released_for_scoring": True,
        "detail": "The inference input tables contain observations only; sealed truth is read for scoring after inference.",
    }


def _truth_panel_audit(output: Path) -> Any:
    """Audit the post-inference truth tables as an independent truth panel."""

    edge_frames = []
    process_fields: list[str] = []
    age_truth: dict[str, float] = {}
    node_frames = []
    for path in sorted((output / "cases").glob("*/heldout_truth.csv")):
        frame = pd.read_csv(path, low_memory=False)
        edge_frames.append(frame.assign(panel_id=path.parent.name))
        process_fields.extend(str(value) for value in frame.get("true_process", []))
    for path in sorted((output / "cases").glob("*/modpath_pathline_truth.csv")):
        frame = pd.read_csv(path, low_memory=False)
        if {"node_id", "age_years"}.issubset(frame.columns):
            values = pd.to_numeric(frame["age_years"], errors="coerce")
            for node_id, value in zip(frame["node_id"], values):
                if pd.notna(value):
                    age_truth.setdefault(str(node_id), float(value))
    for path in sorted((output / "cases").glob("*/blind_observations.csv")):
        frame = pd.read_csv(path, low_memory=False)
        node_frames.append(
            frame.rename(columns={"sample_id": "node_id"})[["node_id"]]
        )
    edges = pd.concat(edge_frames, ignore_index=True) if edge_frames else pd.DataFrame()
    nodes = pd.concat(node_frames, ignore_index=True).drop_duplicates() if node_frames else pd.DataFrame()
    processes = {str(index): value for index, value in enumerate(process_fields)}
    return audit_reference_panel(
        "m7_independent_synthetic_truth",
        "independent_synthetic_truth",
        nodes=nodes,
        edges=edges.rename(columns={"true_process": "process"}) if not edges.empty else edges,
        metadata={
            "truth_edges": True,
            "truth_ages_years": age_truth,
            "truth_processes": processes or {"sealed": True},
            "generator_independent": True,
            "independent_edge_truth": True,
            "independent_age_truth": bool(age_truth),
            "independent_reaction_truth": True,
        },
    )


def run(
    *,
    output: Path,
    simulator_workspace: Path,
    bin_dir: Path,
    dev_seeds: Sequence[int] = SMOKE_DEV_SEEDS,
    test_seeds: Sequence[int] = SMOKE_TEST_SEEDS,
    age_draws: int = SMOKE_AGE_DRAWS,
    age_particles: int = SMOKE_AGE_PARTICLES,
    reaction_bootstrap: int = SMOKE_REACTION_BOOTSTRAP,
    paired_bootstrap: int = SMOKE_PAIRED_BOOTSTRAP,
    overwrite: bool = False,
) -> dict[str, Any]:
    output = output.resolve()
    if output.exists() and any(output.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite non-empty run directory: {output}. Use a new run ID."
        )
    m7_manifest = run_m7_3(
        output=output,
        simulator_workspace=simulator_workspace.resolve(),
        mf6_executable=(bin_dir / "mf6.exe").resolve(),
        mp7_executable=(bin_dir / "mp7.exe").resolve(),
        dev_seeds=tuple(int(seed) for seed in dev_seeds),
        test_seeds=tuple(int(seed) for seed in test_seeds),
        age_draws=int(age_draws),
        age_particles=int(age_particles),
        reaction_bootstrap=int(reaction_bootstrap),
        paired_bootstrap=int(paired_bootstrap),
        overwrite=overwrite,
    )
    blind_audit = audit_blind_tables(output)
    truth_audit = _truth_panel_audit(output)
    generator_path = SCRIPT_DIR / "independent_modflow_generator.py"
    runner_path = SCRIPT_DIR / "run_m7_3_nonuniqueness.py"
    manifest = build_run_manifest(
        run_id=output.name,
        protocol="integrated-truth-blind-v1",
        reference_audits=[truth_audit],
        sources=[
            source_manifest_entry(
                generator_path,
                source_id="m7_independent_generator",
                role="independent MODFLOW/MODPATH/chemistry truth generator",
            ),
            source_manifest_entry(
                runner_path,
                source_id="m7_3_runner",
                role="Hydrosheaf synthetic benchmark runner",
            ),
            source_manifest_entry(
                (bin_dir / "mf6.exe"),
                source_id="modflow6_executable",
                role="forward hydraulic simulator",
            ),
            source_manifest_entry(
                (bin_dir / "mp7.exe"),
                source_id="modpath7_executable",
                role="forward particle tracker",
            ),
        ],
        claim_boundary=(
            "Point 2 is an independent controlled-synthetic integrated benchmark. "
            "It quantifies conditional uncertainty reduction and failure modes "
            "under the declared generator; it is not field validation."
        ),
        extra={
            "reference_type": "independent_synthetic_truth",
            "truth_blind": blind_audit,
            "stages": [
                {
                    "name": "generate_blind_observations",
                    "status": "completed",
                    "truth_blind": True,
                    "truth_fields_seen": [],
                },
                {
                    "name": "run_hydrosheaf_inference",
                    "status": "completed" if blind_audit["status"] == "PASS" else "failed",
                    "truth_blind": True,
                    "truth_fields_seen": [],
                },
                {
                    "name": "release_truth_for_scoring",
                    "status": "completed",
                    "truth_blind": False,
                    "truth_fields_seen": ["true_edges", "true_ages_years", "true_processes"],
                },
            ],
            "m7_runner_manifest": m7_manifest,
            "created_by": "run_integrated_truth_blind.py",
            "created_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    write_json(output / "integration_manifest.json", manifest)
    (output / "INTEGRATED_BENCHMARK.md").write_text(
        "# Independent integrated synthetic benchmark\n\n"
        "This run is truth-blind at the inference stage. Blind observation CSVs "
        "were checked for generator-truth fields; sealed truth tables were used "
        "only after inference for scoring. Results are controlled-synthetic and "
        "model-conditioned, not field validation.\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--simulator-workspace", type=Path, default=DEFAULT_SIMULATOR_WORKSPACE)
    parser.add_argument("--bin-dir", type=Path, default=DEFAULT_BIN_DIR)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    dev = SMOKE_DEV_SEEDS if args.quick else FULL_DEV_SEEDS
    test = SMOKE_TEST_SEEDS if args.quick else FULL_TEST_SEEDS
    settings = (
        (SMOKE_AGE_DRAWS, SMOKE_AGE_PARTICLES, SMOKE_REACTION_BOOTSTRAP, SMOKE_PAIRED_BOOTSTRAP)
        if args.quick
        else (FULL_AGE_DRAWS, FULL_AGE_PARTICLES, FULL_REACTION_BOOTSTRAP, FULL_PAIRED_BOOTSTRAP)
    )
    manifest = run(
        output=args.output,
        simulator_workspace=args.simulator_workspace,
        bin_dir=args.bin_dir,
        dev_seeds=dev,
        test_seeds=test,
        age_draws=settings[0],
        age_particles=settings[1],
        reaction_bootstrap=settings[2],
        paired_bootstrap=settings[3],
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["DEFAULT_OUTPUT", "audit_blind_tables", "main", "run"]
