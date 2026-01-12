"""Redox classification and constraint logic."""

from typing import Dict, Iterable, Mapping, Optional, Tuple

def classify_redox(sample: Mapping[str, float]) -> str:
    """
    Classify the redox state of a sample based on Nitrate and Iron.
    Returns: 'oxic', 'reducing', or 'ambiguous'.
    """
    no3 = sample.get("NO3", 0.0)
    fe = sample.get("Fe", 0.0)
    
    # Thresholds in mmol/L
    if no3 > 0.05:
        return "oxic"
    elif fe > 0.01:
        return "reducing"
    else:
        return "ambiguous"

def get_redox_constraints(sample_v: Mapping[str, float], labels: Iterable[str]) -> Dict[str, Tuple[float, float]]:
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
