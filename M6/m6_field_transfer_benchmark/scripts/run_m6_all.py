"""One-command M6 field-transfer benchmark: analysis -> aux data -> tables.

R figures are generated separately:
  Rscript r_figures/plot_m6_publication_figures.R
  Rscript r_figures/plot_m6_supplementary_figures.R
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable


def run(script):
    print(f"\n=== {script} ===")
    subprocess.run([PY, str(HERE / script)], check=True)


def main():
    run("run_m6_field_transfer.py")
    run("run_m6_synthetic_validation.py")   # ground-truth validation (non-circular)
    run("run_m6_reviewer_response.py")       # gate on/off, transport, CBE, CIs
    run("export_m6_figure_data.py")
    run("make_m6_tables.py")
    print("\nM6 analysis + tables complete. Now run the R figure scripts:")
    print("  Rscript r_figures/plot_m6_publication_figures.R")
    print("  Rscript r_figures/plot_m6_validation_figures.R")
    print("  Rscript r_figures/plot_m6_supplementary_figures.R")


if __name__ == "__main__":
    main()
