"""Reusable validation metrics for Hydrosheaf benchmark scripts."""

from __future__ import annotations

import math
from typing import Iterable, Mapping, Optional, Sequence, Tuple


def _finite_pairs(
    observed: Iterable[object],
    predicted: Iterable[object],
) -> Tuple[list, list]:
    obs = []
    pred = []
    for o_raw, p_raw in zip(observed, predicted):
        try:
            o_val = float(o_raw)
            p_val = float(p_raw)
        except (TypeError, ValueError):
            continue
        if math.isfinite(o_val) and math.isfinite(p_val):
            obs.append(o_val)
            pred.append(p_val)
    return obs, pred


def regression_metrics(
    observed: Iterable[object],
    predicted: Iterable[object],
    *,
    log10: bool = False,
) -> Mapping[str, float]:
    """Compute manuscript-friendly regression diagnostics.

    Non-finite values are skipped pairwise. When ``log10`` is true, only
    positive observed and predicted pairs are retained.
    """

    obs, pred = _finite_pairs(observed, predicted)
    if log10:
        positive = [(o, p) for o, p in zip(obs, pred) if o > 0 and p > 0]
        obs = [math.log10(o) for o, _ in positive]
        pred = [math.log10(p) for _, p in positive]
    n = len(obs)
    if n == 0:
        return {
            "n": 0.0,
            "mae": float("nan"),
            "rmse": float("nan"),
            "bias": float("nan"),
            "nse": float("nan"),
        }

    residuals = [p - o for o, p in zip(obs, pred)]
    mae = sum(abs(res) for res in residuals) / n
    rmse = math.sqrt(sum(res * res for res in residuals) / n)
    bias = sum(residuals) / n
    mean_obs = sum(obs) / n
    denom = sum((o - mean_obs) ** 2 for o in obs)
    nse = 1.0 - (sum(res * res for res in residuals) / denom) if denom > 0 else float("nan")
    return {"n": float(n), "mae": mae, "rmse": rmse, "bias": bias, "nse": nse}


def interval_coverage(
    observed: Iterable[object],
    lower: Iterable[object],
    upper: Iterable[object],
) -> Mapping[str, float]:
    """Compute empirical coverage for prediction or credible intervals."""

    total = 0
    covered = 0
    widths = []
    for obs_raw, lo_raw, hi_raw in zip(observed, lower, upper):
        try:
            obs = float(obs_raw)
            lo = float(lo_raw)
            hi = float(hi_raw)
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(obs) and math.isfinite(lo) and math.isfinite(hi)):
            continue
        if hi < lo:
            lo, hi = hi, lo
        total += 1
        covered += int(lo <= obs <= hi)
        widths.append(hi - lo)
    if total == 0:
        return {"n": 0.0, "coverage": float("nan"), "mean_width": float("nan")}
    return {
        "n": float(total),
        "coverage": covered / total,
        "mean_width": sum(widths) / total,
    }


def classification_metrics(
    expected: Iterable[object],
    predicted: Iterable[object],
    *,
    positive_label: Optional[object] = None,
) -> Mapping[str, float]:
    """Compute accuracy and binary precision/recall/F1 when possible."""

    pairs = [(e, p) for e, p in zip(expected, predicted) if e is not None and p is not None]
    n = len(pairs)
    if n == 0:
        return {"n": 0.0, "accuracy": float("nan")}
    accuracy = sum(1 for e, p in pairs if e == p) / n
    out = {"n": float(n), "accuracy": accuracy}
    if positive_label is None:
        return out
    tp = sum(1 for e, p in pairs if e == positive_label and p == positive_label)
    fp = sum(1 for e, p in pairs if e != positive_label and p == positive_label)
    fn = sum(1 for e, p in pairs if e == positive_label and p != positive_label)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    out.update({"precision": precision, "recall": recall, "f1": f1})
    return out


def topology_metrics(
    reference_edges: Iterable[Sequence[object]],
    inferred_edges: Iterable[Sequence[object]],
) -> Mapping[str, float]:
    """Compute directed-edge precision, recall, and F1 for graph benchmarks."""

    ref = {(str(edge[0]), str(edge[1])) for edge in reference_edges if len(edge) >= 2}
    inf = {(str(edge[0]), str(edge[1])) for edge in inferred_edges if len(edge) >= 2}
    tp = len(ref & inf)
    fp = len(inf - ref)
    fn = len(ref - inf)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "reference_edges": float(len(ref)),
        "inferred_edges": float(len(inf)),
        "tp": float(tp),
        "fp": float(fp),
        "fn": float(fn),
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }
