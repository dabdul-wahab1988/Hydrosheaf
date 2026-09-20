"""Run the bounded static/dynamic graph-TTD virtual benchmark package.

This command does not download WATRES and does not treat Ghana field data as
validation.  Its output records controlled-synthetic execution, the claim
readiness gate, and the still-deferred field-validation boundary.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from hydrosheaf.benchmarks.ttd_graph_virtual_runner import (
    run_ttd_graph_virtual_benchmark,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs") / "ttd_graph_virtual_benchmark_v1",
        help="New output directory for immutable benchmark artifacts.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs") / "ttd_graph_virtual_benchmark_v1.json",
        help="Frozen JSON protocol configuration.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly allow replacement of artifacts at --output.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    result = run_ttd_graph_virtual_benchmark(
        args.output,
        config_path=args.config,
        overwrite=bool(args.overwrite),
    )
    print(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True, default=str))
    # A programme claim may correctly remain ABSTAIN while all execution
    # artifacts are valid.  Process completion, not scientific superiority,
    # determines this command's exit code.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
