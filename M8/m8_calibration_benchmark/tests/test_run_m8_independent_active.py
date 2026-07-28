from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_m8_independent_active import (  # noqa: E402
    heuristic_portability_test,
    numerical_transport_truth,
)


def test_numerical_truth_is_bounded_and_front_increases():
    values = numerical_transport_truth((20.0, 60.0, 120.0, 240.0), cells=80)
    assert np.all(np.isfinite(values))
    assert np.all((values >= 0.0) & (values <= 1.0 + 1e-8))
    assert np.all(np.diff(values) >= -1e-8)


def test_existing_heuristic_is_not_silently_mapped_to_transport_time():
    _, audit = heuristic_portability_test()
    assert audit["status"] == "NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION"
    assert audit["explicit_transport_time_action"] is False
