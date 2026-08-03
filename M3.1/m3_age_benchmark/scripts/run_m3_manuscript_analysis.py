"""Run the M3 manuscript analysis pipeline."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = BENCHMARK_ROOT / "scripts"
RESULT_DIR = BENCHMARK_ROOT / "results"
DOCS_DIR = BENCHMARK_ROOT / "docs"
FIGURE_DIR = BENCHMARK_ROOT / "figures" / "Manuscript_Ready"
TABLE_DIR = BENCHMARK_ROOT / "tables" / "Manuscript_Ready"


def _run(cmd: list[str]) -> None:
    print("Running: " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def _artifact(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
    return {
        "path": str(path),
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else 0,
        "sha256": digest,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run M3 manuscript analysis.")
    parser.add_argument("--max-rows", type=int, default=80, help="Stratified row limit for the design matrix.")
    parser.add_argument("--full", action="store_true", help="Run all eligible USGS rows.")
    parser.add_argument(
        "--age-steps",
        type=int,
        default=90,
        help="Age-grid steps; defaults to the canonical M3 validation resolution.",
    )
    parser.add_argument("--skip-selection", action="store_true")
    parser.add_argument(
        "--skip-cross-validation",
        action="store_true",
        help="Reuse existing tracer-withholding outputs; their manifests remain part of the result lock.",
    )
    parser.add_argument("--cv-max-workers", type=int, default=None)
    parser.add_argument(
        "--regenerate-only",
        action="store_true",
        help="Reuse locked design results and regenerate downstream QA, tables, figures, and hashes.",
    )
    parser.add_argument(
        "--allow-coarse-full-grid",
        action="store_true",
        help="Allow a full run to use age_steps below the canonical M3 default.",
    )
    args = parser.parse_args(argv)

    for directory in (RESULT_DIR, DOCS_DIR, FIGURE_DIR, TABLE_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    design_cmd = [
        sys.executable,
        str(SCRIPT_DIR / "run_m3_design_matrix.py"),
        "--age-steps",
        str(args.age_steps),
    ]
    if args.full:
        design_cmd.append("--full")
    else:
        design_cmd.extend(["--max-rows", str(args.max_rows)])
    if args.allow_coarse_full_grid:
        design_cmd.append("--allow-coarse-full-grid")
    if args.skip_selection:
        design_cmd.append("--skip-selection")

    if not args.regenerate_only:
        _run(design_cmd)
    elif not (RESULT_DIR / "m3_design_matrix_results.csv").exists():
        parser.error("--regenerate-only requires m3_design_matrix_results.csv")
    _run([sys.executable, str(SCRIPT_DIR / "audit_m3_gas_corrections.py")])
    _run([
        sys.executable,
        str(SCRIPT_DIR / "run_m3_real_usgs_graph_benchmark.py"),
        "--pointwise-results",
        str(RESULT_DIR / "m3_design_matrix_results.csv"),
        "--scenario",
        "tracerlpm_strict_parity",
    ])
    if not args.skip_cross_validation:
        for tracer in ("3H", "SF6", "14C", "CFC11", "CFC12"):
            cv_cmd = [
                sys.executable,
                str(SCRIPT_DIR / "run_m3_cross_validation_benchmark.py"),
                "--withhold-tracer",
                tracer,
                "--age-steps",
                str(args.age_steps),
            ]
            if args.cv_max_workers is not None:
                cv_cmd.extend(["--max-workers", str(args.cv_max_workers)])
            _run(cv_cmd)
    _run([sys.executable, str(SCRIPT_DIR / "make_publication_tables.py")])
    _run([sys.executable, str(SCRIPT_DIR / "make_publication_figures.py")])

    artifacts = [
        RESULT_DIR / "m3_design_matrix_results.csv",
        RESULT_DIR / "m3_design_matrix_summary.csv",
        RESULT_DIR / "m3_design_matrix_pairwise_deltas.csv",
        RESULT_DIR / "m3_usgs_benchmark_results.csv",
        RESULT_DIR / "m3_gas_correction_audit.csv",
        RESULT_DIR / "m3_real_usgs_graph_benchmark.csv",
        RESULT_DIR / "m3_real_usgs_graph_edges.csv",
        RESULT_DIR / "m3_real_usgs_graph_benchmark_manifest.json",
        *[
            RESULT_DIR / f"m3_cv_benchmark_{tracer}{suffix}"
            for tracer in ("3H", "SF6", "14C", "CFC11", "CFC12")
            for suffix in (".csv", "_manifest.json")
        ],
        DOCS_DIR / "m3_design_matrix_qa.md",
        DOCS_DIR / "m3_gas_correction_audit.md",
        DOCS_DIR / "m3_graph_benchmark_qa.md",
        BENCHMARK_ROOT / "DECISIONS.md",
        BENCHMARK_ROOT / "configs" / "design_matrix.yaml",
        *sorted(TABLE_DIR.glob("*.csv")),
        *sorted(
            path
            for path in FIGURE_DIR.iterdir()
            if path.is_file() and path.suffix.lower() in {".png", ".pdf", ".tif", ".tiff"}
        ),
    ]
    manifest = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "max_rows": None if args.full else args.max_rows,
        "age_steps": args.age_steps,
        "skip_selection": args.skip_selection,
        "cross_validation_reused": args.skip_cross_validation,
        "fit_results_reused": args.regenerate_only,
        "artifacts": [_artifact(path) for path in artifacts],
    }
    manifest_path = RESULT_DIR / "m3_manuscript_analysis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Wrote M3 manuscript analysis manifest to {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
