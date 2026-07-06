"""Run the complete M5 analysis, table, and figure workflow."""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
BENCHMARK_DIR = SCRIPT_DIR.parent


R_SCRIPT_CANDIDATES = [
    "Rscript",
    r"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.3.2\bin\Rscript.exe",
    r"C:\Program Files\R\R-4.4.0\bin\x64\Rscript.exe",
    r"C:\Program Files\R\R-4.4.0\bin\Rscript.exe",
]


def run(script: str, extra: list[str] | None = None) -> None:
    command = [sys.executable, str(SCRIPT_DIR / script)]
    if extra:
        command.extend(extra)
    subprocess.run(command, check=True)


def find_rscript() -> str | None:
    for candidate in R_SCRIPT_CANDIDATES:
        found = shutil.which(candidate)
        if found:
            return found
        path = Path(candidate)
        if path.exists():
            return str(path)
    return None


def run_r_figures() -> None:
    rscript = find_rscript()
    if rscript is None:
        print("Rscript not found; skipped Nature-style R figure layer.")
        return
    subprocess.run(
        [rscript, str(BENCHMARK_DIR / "r_figures" / "plot_m5_publication_figures.R")],
        cwd=BENCHMARK_DIR,
        check=True,
    )


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
    parser.add_argument(
        "--skip-r-figures",
        action="store_true",
        help="Skip the optional Nature-style R publication figure layer.",
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
    run("export_m5_results_database.py")
    run("make_m5_database_figures.py")
    run("make_m5_publication_figures.py")
    if not args.skip_r_figures:
        run_r_figures()
    run("export_m5_artifact_manifest.py")
    print("M5 workflow complete.")


if __name__ == "__main__":
    main()
