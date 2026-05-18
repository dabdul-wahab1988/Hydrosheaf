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
        description="Run the M2 analysis contracts needed by publication tables and figures."
    )
    parser.add_argument(
        "--realisations",
        type=int,
        default=None,
        help="Override the locked synthetic benchmark realisation count.",
    )
    parser.add_argument(
        "--psi-trials",
        type=int,
        default=30,
        help="Monte Carlo trials per field edge for PSI.",
    )
    parser.add_argument(
        "--max-psi-edges-per-site",
        type=int,
        default=None,
        help="Optional smoke-test edge limit per field site.",
    )
    parser.add_argument("--skip-synthetic", action="store_true", help="Do not rebuild synthetic M2 outputs.")
    parser.add_argument("--skip-field", action="store_true", help="Do not rebuild Ghana field edge outputs.")
    parser.add_argument("--skip-psi", action="store_true", help="Do not rebuild field edge PSI outputs.")
    parser.add_argument("--skip-tables", action="store_true", help="Do not rebuild publication Markdown tables.")
    parser.add_argument("--skip-figures", action="store_true", help="Do not rebuild manuscript-ready figures.")
    parser.add_argument(
        "--run-m3-public-age",
        action="store_true",
        help="Refresh the canonical full M3 screened USGS public-age benchmark used by M2 Fig. 5.",
    )
    args = parser.parse_args()

    if not args.skip_synthetic:
        cmd = [sys.executable, str(BENCHMARK_ROOT / "scripts" / "run_m2_benchmark.py")]
        if args.realisations is not None:
            cmd.extend(["--realisations", str(args.realisations)])
        _run(cmd)

    if not args.skip_field:
        _run([sys.executable, str(PROJECT_ROOT / "scripts" / "analysis" / "run_m2_field_benchmarks.py")])

    if not args.skip_psi:
        cmd = [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "analysis" / "run_edge_psi.py"),
            "--trials",
            str(args.psi_trials),
        ]
        if args.max_psi_edges_per_site is not None:
            cmd.extend(["--max-edges-per-site", str(args.max_psi_edges_per_site)])
        _run(cmd)

    if args.run_m3_public_age:
        _run(
            [
                sys.executable,
                str(PROJECT_ROOT / "M3" / "m3_age_benchmark" / "scripts" / "run_m3_design_matrix.py"),
                "--full",
                "--age-steps",
                "35",
                "--scenario",
                "screened_dgm_gases",
                "--output",
                str(PROJECT_ROOT / "M3" / "m3_age_benchmark" / "results" / "m3_phase4_screened_full_results.csv"),
            ]
        )

    if not args.skip_tables:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "make_publication_tables.py")])

    if not args.skip_figures:
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "make_publication_figures.py")])
        _run([sys.executable, str(BENCHMARK_ROOT / "scripts" / "make_supplementary_figures.py")])

    print("M2 manuscript analysis complete.", flush=True)


if __name__ == "__main__":
    main()
