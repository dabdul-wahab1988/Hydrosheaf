"""Manuscript-artifact generation entrypoint for the M7.3 benchmark.

This module documents and re-invokes the already-locked M7.3 pipeline; it
performs no new analysis of its own. Replay requires ``mf6.exe``/``mp7.exe``
on the path documented in ``README.md`` and is expected to reproduce the
locked outputs already committed under ``results/m7_3_locked``.

Usage:
    python generation_script.py            # full locked benchmark
    python generation_script.py --assets    # regenerate figures/tables only
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
    return subprocess.run(
        [sys.executable, str(SCRIPTS / "make_m7_3_publication_assets.py")],
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
    args = parser.parse_args()
    if args.assets:
        return run_publication_assets()
    benchmark_status = run_benchmark()
    if benchmark_status != 0:
        return benchmark_status
    return run_publication_assets()


if __name__ == "__main__":
    raise SystemExit(main())
