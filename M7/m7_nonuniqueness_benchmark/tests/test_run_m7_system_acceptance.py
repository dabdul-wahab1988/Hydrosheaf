from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_m7_system_acceptance import (  # noqa: E402
    _bootstrap_contrast,
    _prepare_rows,
)


class _Case:
    seed = 9
    observations = (
        {"site_id": "A", "tritium_TU": 1.0},
        {"site_id": "B", "tritium_TU": 2.0},
    )


def test_prepare_rows_withholds_truth_and_does_not_mutate_case():
    rows = _prepare_rows(_Case(), permute_age=False)
    assert [row["3H"] for row in rows] == [1.0, 2.0]
    assert "3H" not in _Case.observations[0]


def test_bootstrap_contrast_preserves_paired_direction():
    frame = pd.DataFrame(
        [
            {"seed": seed, "condition": condition, "pr_auc": value}
            for seed, left, right in ((1, 0.8, 0.4), (2, 0.7, 0.5))
            for condition, value in (("left", left), ("right", right))
        ]
    )
    result = _bootstrap_contrast(
        frame, "left", "right", "pr_auc", samples=500, seed=1
    )
    assert result["mean_difference"] > 0.0
    assert result["ci95_low"] > 0.0
