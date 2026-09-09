"""Orchestrate the controlled and Aiken model-conditioned age tests.

The two panels are intentionally run into separate child directories and are
linked only by this manifest.  The synthetic panel has sealed independent
truth and may score the direct-versus-indirect Bayes-factor estimand.  The
Aiken panel is a calibrated-model reference with no independent well-to-well
truth and is therefore restricted to transport/concordance diagnostics.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hydrosheaf.validation.aiken_emulation import (  # noqa: E402
    run_aiken_model_conditioned_emulation,
    write_aiken_emulation_outputs,
)
from hydrosheaf.validation.aiken_reference import load_aiken_reference  # noqa: E402

from age_bayes_factor_benchmark import run_benchmark  # noqa: E402


def run_two_tier(
    *,
    aiken_root: Path,
    output: Path,
    quick: bool = False,
    first_seed: int = 2026090801,
    n_cases: int = 24,
    n_development: int | None = None,
) -> dict[str, object]:
    """Run both panels without overwriting an existing orchestration root."""

    output = output.expanduser().resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if quick:
        n_cases = 4
        n_development = 2
    elif n_development is None:
        n_development = n_cases // 2

    synthetic_dir = output / "controlled_synthetic"
    aiken_dir = output / "aiken_model_conditioned"
    synthetic_manifest = run_benchmark(
        output=synthetic_dir,
        first_seed=first_seed,
        n_cases=n_cases,
        n_development=n_development,
    )
    reference = load_aiken_reference(aiken_root.expanduser().resolve())
    aiken_result = run_aiken_model_conditioned_emulation(reference)
    aiken_outputs = write_aiken_emulation_outputs(aiken_result, aiken_dir)
    orchestration = {
        "schema": "age-adjacency-two-tier-orchestration-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": "age-adjacency-two-tier-v1",
        "panels": {
            "controlled_synthetic": {
                "path": "controlled_synthetic",
                "run_id": synthetic_dir.name,
                "claim_tier": "controlled_synthetic_component",
                "integrated_scoring_allowed": True,
                "qa_decision": synthetic_manifest["tier_package"]["qa"]["decision"],
                "complete_truth_counts": synthetic_manifest["complete_truth_counts"],
            },
            "aiken_model_conditioned": {
                "path": "aiken_model_conditioned",
                "run_id": aiken_result.manifest["emulation_id"],
                "claim_tier": "calibrated_model_reference",
                "integrated_scoring_allowed": False,
                "direct_well_to_well_truth_emitted": False,
                "outputs": aiken_outputs,
                "counts": aiken_result.manifest["outputs"],
            },
        },
        "separation_rule": (
            "Never pool Aiken model-conditioned rows with synthetic truth rows; "
            "only the controlled-synthetic panel supplies an independent "
            "direct-adjacency denominator."
        ),
        "claim_boundary": (
            "The controlled panel tests the declared transport Bayes factor under "
            "synthetic truth. Aiken supplies calibrated-model transport and age "
            "concordance diagnostics, not independent field validation."
        ),
    }
    manifest_path = output / "orchestration_manifest.json"
    manifest_path.write_text(
        json.dumps(orchestration, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return orchestration


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aiken-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--first-seed", type=int, default=2026090801)
    parser.add_argument("--n-cases", type=int, default=24)
    parser.add_argument("--n-development", type=int, default=None)
    args = parser.parse_args(argv)
    result = run_two_tier(
        aiken_root=args.aiken_root,
        output=args.output,
        quick=args.quick,
        first_seed=args.first_seed,
        n_cases=args.n_cases,
        n_development=args.n_development,
    )
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())


__all__ = ["run_two_tier", "main"]
