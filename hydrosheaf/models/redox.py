"""Redox classification and constraint logic."""

import math
from typing import Dict, Iterable, Mapping, Tuple


def _finite(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def classify_redox(sample: Mapping[str, object]) -> str:
    """
    Classify the redox state of a sample based on Nitrate and Iron.
    Returns: 'oxic', 'reducing', or 'ambiguous'.
    """
    no3 = _finite(sample.get("NO3"))
    fe = _finite(sample.get("Fe"))

    # Thresholds in mmol/L
    if no3 is not None and no3 > 0.05:
        return "oxic"
    elif fe is not None and fe > 0.01:
        return "reducing"
    else:
        return "ambiguous"


def get_redox_constraints(
    sample_v: Mapping[str, object], labels: Iterable[str]
) -> Dict[str, Tuple[float, float]]:
    """
    Determine mineral bounds overrides based on redox state of the downstream sample.
    """
    state = classify_redox(sample_v)
    overrides = {}

    if state == "reducing":
        # Prevent aerobic pyrite oxidation if oxygen is likely absent
        for i, label in enumerate(labels):
            if "pyrite_oxidation_aerobic" in label:
                overrides[label] = (0.0, 0.0)  # Forced 0

    return overrides
