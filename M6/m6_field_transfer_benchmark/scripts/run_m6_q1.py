"""Run the approved refined-cohort M6 input-QA workflow only.

Legacy field-transfer tables and figures depend on a seasonal/older-cohort
design and are intentionally not regenerated. The supported route audits the
completed Central Region and Upper East Region workbooks without inferring
temporal transfer, field flow direction, or reaction truth.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCHMARK = HERE.parent
ANALYSIS_STEPS = ("run_m6_refined_cross_sectional.py",)
R_STEPS: tuple[str, ...] = ()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-only",
        action="store_true",
        help="Accepted for compatibility; the supported M6 workflow is input QA only.",
    )
    parser.add_argument(
        "--figures-only",
        action="store_true",
        help="Refuse legacy field figures, which were generated from excluded cohorts.",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    if args.figures_only:
        parser.error(
            "Legacy M6 figures are not valid for the refined cohorts. "
            "Regenerate displays only after a separately reviewed figure workflow is added."
        )
    command = [sys.executable, str(HERE / ANALYSIS_STEPS[0])]
    if args.output is not None:
        command.extend(["--output", str(args.output)])
    subprocess.run(command, cwd=BENCHMARK, check=True)
    print(
        "M6 refined field-input QA completed. No temporal, flow-direction, "
        "or field reaction-truth results were produced."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
