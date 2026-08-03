from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from post_review_statistical_audit import (  # noqa: E402
    case_metrics_from_predictions,
    simultaneous_family,
)


def test_case_metrics_from_predictions_are_case_blocked() -> None:
    frame = pd.DataFrame(
        {
            "seed": [1, 1, 1, 1],
            "scenario": ["x"] * 4,
            "model": ["a"] * 2 + ["b"] * 2,
            "is_true_edge": [1, 0, 1, 0],
            "probability": [0.8, 0.2, 0.7, 0.3],
        }
    )
    result = case_metrics_from_predictions(frame, group_columns=("seed", "scenario"))
    assert set(result["model"]) == {"a", "b"}
    assert np.isfinite(result[["pr_auc", "brier", "log_loss"]]).all().all()


def test_simultaneous_interval_is_finite() -> None:
    rows = []
    for seed in range(1, 9):
        for model, shift in (("a", 0.05), ("b", 0.0)):
            rows.append(
                {
                    "seed": seed,
                    "scenario": "x" if seed <= 4 else "y",
                    "model": model,
                    "pr_auc": 0.5 + shift + 0.005 * seed,
                }
            )
    metrics = pd.DataFrame(rows)
    contrasts = pd.DataFrame(
        [
            {
                "scenario": "all",
                "left": "a",
                "right": "b",
                "metric": "pr_auc",
                "n_cases": 8,
                "mean_difference": 0.05,
                "ci95_low": 0.04,
                "ci95_high": 0.06,
            }
        ]
    )
    adjusted, critical = simultaneous_family(
        metrics, contrasts, family_name="test", rng_seed=7
    )
    assert critical >= 0
    assert np.isfinite(adjusted.loc[0, "simultaneous_ci95_low"])
    assert np.isfinite(adjusted.loc[0, "simultaneous_ci95_high"])
