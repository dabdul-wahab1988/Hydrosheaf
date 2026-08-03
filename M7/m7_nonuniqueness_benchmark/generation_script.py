"""Manuscript-artifact generation entrypoint for the maintained M7 package.

This module documents and re-invokes the already-locked M7.3 pipeline; it
performs no new analysis of its own. Replay requires ``mf6.exe``/``mp7.exe``
on the path documented in ``README.md`` and is expected to reproduce the
locked outputs already committed under ``results/m7_3_locked``.

Usage:
    python generation_script.py                    # full locked M7.3 benchmark
    python generation_script.py --assets           # all publication assets only
    python generation_script.py --sheaf-vs-graph   # locked M7.4 comparator
    python generation_script.py --robust-hybrid-assets  # immutable M7.5 assets
    python generation_script.py --all              # analyses plus all assets
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"


def run_benchmark() -> int:
    return subprocess.run(
        [sys.executable, str(SCRIPTS / "run_m7_3_nonuniqueness.py")],
        cwd=ROOT,
        check=False,
    ).returncode


def run_publication_assets() -> int:
    m7_3_status = subprocess.run(
        [sys.executable, str(SCRIPTS / "make_m7_3_publication_assets.py")],
        cwd=ROOT,
        check=False,
    ).returncode
    if m7_3_status != 0:
        return m7_3_status
    m7_4_status = subprocess.run(
        [sys.executable, str(SCRIPTS / "make_m7_sheaf_vs_graph_assets.py")],
        cwd=ROOT,
        check=False,
    ).returncode
    if m7_4_status != 0:
        return m7_4_status
    m7_5_status = subprocess.run(
        [sys.executable, str(SCRIPTS / "make_m7_robust_hybrid_assets.py")],
        cwd=ROOT,
        check=False,
    ).returncode
    if m7_5_status != 0:
        return m7_5_status
    return subprocess.run(
        [sys.executable, str(SCRIPTS / "assemble_m7_supplement.py")],
        cwd=ROOT,
        check=False,
    ).returncode


def run_sheaf_vs_graph() -> int:
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "run_m7_sheaf_vs_graph.py"),
            "--overwrite",
        ],
        cwd=ROOT,
        check=False,
    ).returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--assets",
        action="store_true",
        help="Regenerate publication figures and tables from locked results only.",
    )
    parser.add_argument(
        "--sheaf-vs-graph",
        action="store_true",
        help="Replay only the prospectively locked M7.4 comparator.",
    )
    parser.add_argument(
        "--robust-hybrid-assets",
        action="store_true",
        help=(
            "Regenerate M7.5 assets from the immutable one-time locked-test "
            "outputs; never reruns the confirmatory test."
        ),
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Replay M7.3 and M7.4, then regenerate all publication assets.",
    )
    args = parser.parse_args()
    if args.assets:
        return run_publication_assets()
    if args.sheaf_vs_graph:
        return run_sheaf_vs_graph()
    if args.robust_hybrid_assets:
        return subprocess.run(
            [sys.executable, str(SCRIPTS / "make_m7_robust_hybrid_assets.py")],
            cwd=ROOT,
            check=False,
        ).returncode
    benchmark_status = run_benchmark()
    if benchmark_status != 0:
        return benchmark_status
    if args.all:
        sheaf_status = run_sheaf_vs_graph()
        if sheaf_status != 0:
            return sheaf_status
    return run_publication_assets()


if __name__ == "__main__":
    raise SystemExit(main())
