"""Run the M4 analysis contracts needed by manuscript-ready figures and tables."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]


def _run(args: list[str]) -> None:
    print("$ " + " ".join(args), flush=True)
    subprocess.run(args, cwd=PROJECT_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run M4 topology validation outputs for tables, figures, and evidence registers."
    )
    parser.add_argument("--skip-controlled", action="store_true", help="Do not rebuild controlled topology outputs.")
    parser.add_argument(
        "--run-modpath-validation",
        action="store_true",
        help="Refresh endpoint-derived MODPATH validation using the selected M2 MODPATH output files.",
    )
    parser.add_argument(
        "--skip-external-archive",
        action="store_true",
        help="Do not ingest the M2 external USGS MODPATH archive validation into M4.",
    )
    parser.add_argument("--skip-tables", action="store_true", help="Do not rebuild manuscript-ready tables.")
    parser.add_argument("--skip-figures", action="store_true", help="Do not rebuild manuscript-ready figures.")
    args = parser.parse_args()

    if not args.skip_controlled:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "run_m4_topology_benchmark.py")])

    if args.run_modpath_validation:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "run_m4_modpath_validation.py")])

    if not args.skip_external_archive:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "run_m4_external_modpath_archive_validation.py")])

    if not args.skip_tables:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "make_publication_tables.py")])

    if not args.skip_figures:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "make_publication_figures.py")])

    print("M4 manuscript analysis complete.", flush=True)


if __name__ == "__main__":
    main()
