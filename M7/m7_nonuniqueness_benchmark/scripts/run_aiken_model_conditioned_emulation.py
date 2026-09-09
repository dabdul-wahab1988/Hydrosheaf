"""Run the secondary Aiken age/transport emulation.

This is an isolated supplementary panel.  It reads the local Aiken release,
decodes the small validated MODPATH 5 particle files, and writes explicit
crosswalk, endpoint, pathline-segment, CFC-interval, summary, and manifest
artifacts.  It does not hash multi-gigabyte archives, create well-to-well
truth labels, or modify any manuscript/result locked by another milestone.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.validation.aiken_emulation import (  # noqa: E402
    run_aiken_model_conditioned_emulation,
    write_aiken_emulation_outputs,
)
from hydrosheaf.validation.aiken_reference import load_aiken_reference  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the isolated Aiken calibrated-model age/transport emulation."
    )
    parser.add_argument(
        "--aiken-root",
        type=Path,
        required=True,
        help="Aiken release directory or ZIP-backed root containing output.zip/model.zip.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Empty or new directory for supplementary emulation artifacts.",
    )
    parser.add_argument(
        "--emulation-id",
        default="AIKEN-MODEL-CONDITIONED-AGE-TRANSPORT",
        help="Stable logical ID recorded in the manifest.",
    )
    parser.add_argument(
        "--reference-year",
        type=int,
        default=2015,
        help="Fallback sample/reference year used only when a CFC sample date is absent.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reference = load_aiken_reference(args.aiken_root.resolve())
    result = run_aiken_model_conditioned_emulation(
        reference,
        emulation_id=args.emulation_id,
        reference_year=args.reference_year,
    )
    generated = write_aiken_emulation_outputs(result, args.output)
    payload = {
        "manifest": result.manifest,
        "generated": generated,
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
