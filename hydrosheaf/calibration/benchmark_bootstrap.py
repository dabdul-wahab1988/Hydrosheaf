"""
Bootstrap uncertainty estimation for assumption-calibration benchmark metrics.

Resamples validation observations with replacement and recomputes
precision / recall / F1 / accuracy using pre-computed selected edge IDs
per variant. Sheaf refinement is **not** re-run; only the metric
aggregation is bootstrapped.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .adapters import TopologyCalibrationObservation

# Small-set threshold for bootstrap warning
_SMALL_SET_THRESHOLD = 20


@dataclass
class BootstrapVariantCI:
    """95% bootstrap percentile CI for a single variant's four metrics."""

    precision: Tuple[float, float]
    recall: Tuple[float, float]
    f1: Tuple[float, float]
    accuracy: Tuple[float, float]


@dataclass
class BootstrapDeltaCI:
    """Paired bootstrap delta CIs (assumption_calibrated - baseline)."""

    delta_precision: Tuple[float, float]
    delta_recall: Tuple[float, float]
    delta_f1: Tuple[float, float]
    delta_accuracy: Tuple[float, float]
    probability_delta_f1_gt_0: float
    probability_delta_precision_gt_0: float


# ── metric computation (no dataclass overhead in the bootstrap loop) ──

def _compute_metrics_from_ids(
    selected_ids: set,
    edge_ids: List[str],
    observed_present: List[float],
) -> Dict[str, float]:
    """Compute precision / recall / F1 / accuracy from pre-computed
    selected edge IDs against a (possibly resampled) set of observations.

    Only edges whose ``edge_id`` appears in *both* ``selected_ids`` and
    the observation list contribute to the confusion matrix — edges that
    are selected but never observed are excluded, consistent with the
    manuscript-safe evaluation rules.
    """
    tp = fp = tn = fn = 0
    for eid, obs_val in zip(edge_ids, observed_present):
        is_selected = eid in selected_ids
        if obs_val >= 0.5:
            if is_selected:
                tp += 1
            else:
                fn += 1
        else:
            if is_selected:
                fp += 1
            else:
                tn += 1

    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = (2.0 * precision * recall / (precision + recall)
          if precision + recall else 0.0)
    accuracy = (tp + tn) / (tp + fp + tn + fn) if (tp + fp + tn + fn) > 0 else float("nan")

    return {"precision": precision, "recall": recall, "f1": f1, "accuracy": accuracy}


# ── main bootstrap entry point ────────────────────────────────────────

def bootstrap_benchmark_metrics(
    selected_ids_by_variant: Dict[str, set],
    val_observations: List[TopologyCalibrationObservation],
    n_boot: int = 1000,
    seed: int = 123,
) -> Dict[str, Any]:
    """Bootstrap 95% percentile CIs for benchmark variant metrics and
    paired improvement deltas.

    Parameters
    ----------
    selected_ids_by_variant : dict
        Maps variant name to a set of selected edge IDs (pre-computed
        from sheaf refinement / thresholding).
    val_observations : list of TopologyCalibrationObservation
        Held-out validation labels.
    n_boot : int
        Number of bootstrap resamples.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict
        ``n_boot``, ``seed``, ``variant_cis``, ``delta_cis``, and an
        optional ``bootstrap_warning`` when the validation set is small.
    """
    rng = np.random.default_rng(seed)
    n_val = len(val_observations)

    warning: Optional[str] = None
    if n_val < _SMALL_SET_THRESHOLD:
        warning = (
            f"Validation set has {n_val} labeled edges "
            f"(< {_SMALL_SET_THRESHOLD}); bootstrap CIs may be unreliable."
        )

    # Pre-extract arrays for fast indexing inside the bootstrap loop
    obs_edge_ids = [obs.edge_id for obs in val_observations]
    obs_present = [obs.observed_present for obs in val_observations]

    variant_names = sorted(selected_ids_by_variant.keys())
    metric_names = ["precision", "recall", "f1", "accuracy"]

    # Accumulators: variant -> metric -> numpy array of bootstrap samples
    samples: Dict[str, Dict[str, np.ndarray]] = {
        vname: {m: np.empty(n_boot) for m in metric_names}
        for vname in variant_names
    }

    for i in range(n_boot):
        indices = rng.integers(0, n_val, size=n_val)
        resampled_ids = [obs_edge_ids[idx] for idx in indices]
        resampled_present = [obs_present[idx] for idx in indices]

        for vname in variant_names:
            metrics = _compute_metrics_from_ids(
                selected_ids_by_variant[vname],
                resampled_ids,
                resampled_present,
            )
            for m in metric_names:
                samples[vname][m][i] = metrics[m]

    # ── Per-variant percentile CIs ─────────────────────────────────
    variant_cis: Dict[str, BootstrapVariantCI] = {}
    for vname in variant_names:
        cis = {}
        for m in metric_names:
            lo, hi = np.percentile(samples[vname][m], [2.5, 97.5])
            cis[m] = (float(np.clip(lo, 0.0, 1.0)), float(np.clip(hi, 0.0, 1.0)))
        variant_cis[vname] = BootstrapVariantCI(**cis)

    # ── Paired delta CIs (calibrated - baseline) ───────────────────
    baseline = "baseline"
    calibrated = "assumption_calibrated"

    if baseline in samples and calibrated in samples:
        delta_samples: Dict[str, np.ndarray] = {}
        for m in metric_names:
            delta_samples[m] = samples[calibrated][m] - samples[baseline][m]

        delta_cis = {}
        for m in metric_names:
            lo, hi = np.percentile(delta_samples[m], [2.5, 97.5])
            delta_cis[m] = (float(np.clip(lo, -1.0, 1.0)), float(np.clip(hi, -1.0, 1.0)))

        prob_f1 = float(np.mean(delta_samples["f1"] > 0))
        prob_prec = float(np.mean(delta_samples["precision"] > 0))

        delta_cis_obj = BootstrapDeltaCI(
            delta_precision=delta_cis["precision"],
            delta_recall=delta_cis["recall"],
            delta_f1=delta_cis["f1"],
            delta_accuracy=delta_cis["accuracy"],
            probability_delta_f1_gt_0=prob_f1,
            probability_delta_precision_gt_0=prob_prec,
        )
    else:
        delta_cis_obj = BootstrapDeltaCI(
            delta_precision=(0.0, 0.0),
            delta_recall=(0.0, 0.0),
            delta_f1=(0.0, 0.0),
            delta_accuracy=(0.0, 0.0),
            probability_delta_f1_gt_0=0.0,
            probability_delta_precision_gt_0=0.0,
        )

    result: Dict[str, Any] = {
        "n_boot": n_boot,
        "seed": seed,
        "variant_cis": variant_cis,
        "delta_cis": delta_cis_obj,
    }
    if warning:
        result["bootstrap_warning"] = warning

    return result
