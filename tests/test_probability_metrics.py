from __future__ import annotations

import math

import pytest

from hydrosheaf.validation.metrics import probability_metrics


def test_probability_metrics_are_finite_and_include_reliability_bins() -> None:
    result = probability_metrics([0, 1, 1, 0], [0.1, 0.9, 0.8, 0.2], n_bins=4)
    assert result["n"] == 4.0
    assert result["invalid"] == 0.0
    assert 0.0 <= float(result["brier"]) <= 1.0
    assert math.isfinite(float(result["log_loss"]))
    assert 0.0 <= float(result["ece"]) <= 1.0
    assert result["reliability_bins"]


def test_invalid_probabilities_are_counted_not_clipped() -> None:
    result = probability_metrics([0, 1, 1], [0.2, 1.2, "bad"], n_bins=5)
    assert result["n"] == 1.0
    assert result["invalid"] == 2.0


def test_probability_metrics_reject_nonpositive_bin_count() -> None:
    with pytest.raises(ValueError, match="n_bins"):
        probability_metrics([0], [0.5], n_bins=0)
