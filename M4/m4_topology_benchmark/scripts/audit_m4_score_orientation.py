"""M4 score-orientation diagnosis and physically declared re-scoring.

This script does two things and nothing else.

1.  Label-free monotonicity audit.  A physically valid edge score must be
    non-decreasing in `feature_direction_probability` (the probability that the
    source head exceeds the target head) and in
    `feature_target_pumping_sink_strength`.  This is checked on the complete
    candidate universe without using the reference labels at all.

2.  Scoring of physically declared criteria.  `direction_probability` is the
    direct implementation of the Darcy direction criterion and is therefore
    declared a priori as the primary physical score.  Its discrimination is
    reported together with that of the other declared channels.  No weight is
    fitted on Savage.

Outputs are written to the directory given on the command line.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
FROZEN = REPO / "M4" / "m4_topology_benchmark" / "results" / "frozen_m4_diagnostics_20260922"
LABEL = "in_reference_for_evaluation_only"


def roc_auc(labels: np.ndarray, scores: np.ndarray) -> float | None:
    from sklearn.metrics import roc_auc_score

    if len(np.unique(labels)) < 2:
        return None
    return float(roc_auc_score(labels, scores))


def pr_auc(labels: np.ndarray, scores: np.ndarray) -> float | None:
    from sklearn.metrics import average_precision_score

    if len(np.unique(labels)) < 2:
        return None
    return float(average_precision_score(labels, scores))


def threshold_metrics(labels: np.ndarray, scores: np.ndarray, thr: float) -> dict[str, float]:
    pred = scores >= thr
    tp = int(np.sum(pred & (labels == 1)))
    fp = int(np.sum(pred & (labels == 0)))
    fn = int(np.sum(~pred & (labels == 1)))
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "threshold": float(thr),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "n_present": int(pred.sum()),
    }


def fdr_threshold(labels: np.ndarray, scores: np.ndarray, target_fdr: float) -> float:
    """Smallest threshold whose achieved false-discovery rate is <= target.

    Selected on the candidate universe using a pre-declared rule.  Returned so
    the caller can record the rule and the resulting threshold explicitly.
    """
    order = np.argsort(-scores)
    s_sorted = scores[order]
    l_sorted = labels[order]
    tp = np.cumsum(l_sorted == 1)
    fp = np.cumsum(l_sorted == 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fdr = fp / np.maximum(tp + fp, 1)
    ok = np.where(fdr <= target_fdr)[0]
    if len(ok) == 0:
        return float(s_sorted[0] + 1.0)
    idx = ok[-1]
    if idx + 1 >= len(s_sorted):
        return float(s_sorted[idx])
    return float(s_sorted[idx + 1])


def monotonicity(table: pd.DataFrame) -> dict[str, object]:
    """Label-free check of score monotonicity in each physical feature."""
    out: dict[str, object] = {}
    score = pd.to_numeric(table["probability"], errors="coerce").to_numpy(float)
    for col in (
        "feature_direction_probability",
        "feature_target_pumping_sink_strength",
        "feature_head_z_score",
        "feature_elevation_direction_probability",
    ):
        if col not in table.columns:
            continue
        feat = pd.to_numeric(table[col], errors="coerce").to_numpy(float)
        mask = np.isfinite(feat) & np.isfinite(score)
        if mask.sum() < 50 or len(np.unique(feat[mask])) < 3:
            out[col] = {"usable": False, "n": int(mask.sum())}
            continue
        # Spearman-style: rank correlation between feature and frozen score.
        rf = pd.Series(feat[mask]).rank().to_numpy()
        rs = pd.Series(score[mask]).rank().to_numpy()
        rho = float(np.corrcoef(rf, rs)[0, 1])
        # Decile means to expose non-monotonicity directly.
        q = pd.qcut(pd.Series(feat[mask]), 5, labels=False, duplicates="drop")
        decile_means = pd.Series(score[mask]).groupby(q).mean().round(6).to_dict()
        out[col] = {
            "usable": True,
            "n": int(mask.sum()),
            "spearman_feature_vs_frozen_score": rho,
            "quintile_mean_frozen_score": {str(k): float(v) for k, v in decile_means.items()},
            "monotone_increasing": bool(rho > 0.0),
        }
    return out


def analyse(mode: str, path: Path, out_dir: Path) -> dict[str, object]:
    table = pd.read_csv(path)
    if LABEL not in table.columns:
        raise SystemExit(f"{path} lacks {LABEL}")
    labels = table[LABEL].astype(float).to_numpy()
    result: dict[str, object] = {
        "mode": mode,
        "source": str(path),
        "n_candidates": int(len(table)),
        "n_positives": int((labels == 1).sum()),
        "prevalence": float((labels == 1).mean()),
        "monotonicity_label_free": monotonicity(table),
        "declared_scores": {},
    }

    candidates = {
        "frozen_probability": pd.to_numeric(table["probability"], errors="coerce").to_numpy(float),
        "diagnostic_complement_1_minus_frozen": 1.0
        - pd.to_numeric(table["probability"], errors="coerce").to_numpy(float),
    }
    for col in (
        "candidate_prior_probability",
        "feature_direction_probability",
        "feature_head_z_score",
        "feature_gradient_m_per_km",
        "feature_target_pumping_sink_strength",
        "feature_elevation_direction_probability",
    ):
        if col in table.columns:
            candidates[col] = pd.to_numeric(table[col], errors="coerce").to_numpy(float)

    for name, score in candidates.items():
        m = np.isfinite(score)
        if m.sum() < 50:
            result["declared_scores"][name] = {"usable": False, "n": int(m.sum())}
            continue
        lab = labels[m]
        s = score[m]
        entry: dict[str, object] = {
            "usable": True,
            "n": int(m.sum()),
            "roc_auc": roc_auc(lab, s),
            "pr_auc": pr_auc(lab, s),
            "prevalence": float(lab.mean()),
            "pr_auc_lift_over_prevalence": (pr_auc(lab, s) or 0.0) / max(float(lab.mean()), 1e-12),
        }
        if name.startswith("feature_") or name.startswith("candidate_") or name == "frozen_probability":
            # Natural 0.5 decision point for probability-like channels.
            if s.min() <= 0.5 <= s.max():
                entry["at_0.5"] = threshold_metrics(lab, s, 0.5)
        result["declared_scores"][name] = entry

    return result


def main() -> int:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "M4" / "m4_topology_benchmark" / "results" / "m4_orientation_audit"
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = []
    for mode, fname in (("M4-A", "M4-A/fair_edge_scores_m4_a.csv"), ("M4-B", "M4-B/fair_edge_scores_m4_b.csv")):
        payload.append(analyse(mode, FROZEN / fname, out_dir))
    (out_dir / "m4_orientation_audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = ["# M4 score-orientation audit", ""]
    for entry in payload:
        lines.append(f"## {entry['mode']}")
        lines.append("")
        lines.append(f"Candidates {entry['n_candidates']}, positives {entry['n_positives']}, prevalence {entry['prevalence']:.5f}.")
        lines.append("")
        lines.append("| Score channel | n | ROC-AUC | PR-AUC | PR-AUC / prevalence |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        for name, m in entry["declared_scores"].items():
            if not m.get("usable"):
                continue
            roc = m["roc_auc"]
            pr = m["pr_auc"]
            lines.append(
                f"| {name} | {m['n']} | "
                f"{roc:.4f} | " if roc is not None else f"| {name} | {m['n']} | n/a | "
            )
            lines[-1] = (
                f"| {name} | {m['n']} | {('%.4f' % roc) if roc is not None else 'n/a'} | "
                f"{('%.4f' % pr) if pr is not None else 'n/a'} | {m['pr_auc_lift_over_prevalence']:.2f} |"
            )
        lines.append("")
        lines.append("Label-free monotonicity of the frozen score in each physical feature:")
        lines.append("")
        lines.append("| Feature | n | Spearman vs frozen score | Monotone increasing |")
        lines.append("| --- | ---: | ---: | --- |")
        for col, m in entry["monotonicity_label_free"].items():
            if not m.get("usable"):
                continue
            lines.append(
                f"| {col} | {m['n']} | {m['spearman_feature_vs_frozen_score']:+.4f} | "
                f"{'yes' if m['monotone_increasing'] else 'NO'} |"
            )
        lines.append("")
    (out_dir / "m4_orientation_audit.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nwritten -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
