"""Run the complete M5 analysis, table, and figure workflow."""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent


def run(script: str, extra: list[str] | None = None) -> None:
    command = [sys.executable, str(SCRIPT_DIR / script)]
    if extra:
        command.extend(extra)
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reuse-synthetic",
        action="store_true",
        help="Reuse existing synthetic CSV outputs and refresh field/displays.",
    )
    parser.add_argument(
        "--reuse-phreeqc",
        action="store_true",
        help="Reuse PHREEQC ground truth and regenerate model outputs/displays.",
    )
    args = parser.parse_args()
    if args.reuse_synthetic and args.reuse_phreeqc:
        parser.error("--reuse-synthetic and --reuse-phreeqc are mutually exclusive.")
    analysis_args = []
    if args.reuse_synthetic:
        analysis_args = ["--reuse-synthetic"]
    elif args.reuse_phreeqc:
        analysis_args = ["--reuse-phreeqc"]
    run("run_m5_inverse_reaction_benchmark.py", analysis_args)
    run("make_m5_publication_tables.py")
    run("make_m5_publication_figures.py")
    print("M5 workflow complete.")


if __name__ == "__main__":
    main()
