"""Run the complete M6 Q1 analysis, tables, and publication figures."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCHMARK = HERE.parent
R_FIGURES = BENCHMARK / "r_figures"
PY = sys.executable

ANALYSIS_STEPS = (
    "run_m6_field_transfer.py",
    "run_m6_synthetic_validation.py",
    "run_m6_robustness_diagnostics.py",
    "run_m6_null_sensitivity.py",
    "export_m6_figure_data.py",
    "make_m6_tables.py",
)
R_STEPS = (
    "plot_m6_publication_figures.R",
    "plot_m6_validation_figures.R",
    "plot_m6_supplementary_figures.R",
)


def _run(command: list[str], label: str) -> None:
    print(f"\n=== {label} ===", flush=True)
    subprocess.run(command, cwd=BENCHMARK, check=True)


def _rscript() -> str:
    found = shutil.which("Rscript")
    if found:
        return found
    windows = Path(r"C:\Program Files\R\R-4.6.1\bin\Rscript.exe")
    if windows.exists():
        return str(windows)
    raise RuntimeError("Rscript was not found; install R 4.3+ or add it to PATH.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--analysis-only",
        action="store_true",
        help="Run analysis and tables, but do not regenerate figures.",
    )
    mode.add_argument(
        "--figures-only",
        action="store_true",
        help="Regenerate figures from existing locked results.",
    )
    args = parser.parse_args()

    if not args.figures_only:
        for script in ANALYSIS_STEPS:
            _run([PY, str(HERE / script)], script)

    if not args.analysis_only:
        _run(
            [PY, str(HERE / "make_objective6_prequential_figure.py")],
            "main Figure 5: field prequential benchmark",
        )
        rscript = _rscript()
        for script in R_STEPS:
            _run([rscript, str(R_FIGURES / script)], script)

    print("\nM6 Q1 workflow completed.", flush=True)


if __name__ == "__main__":
    main()
