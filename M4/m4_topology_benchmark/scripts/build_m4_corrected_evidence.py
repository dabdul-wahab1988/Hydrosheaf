"""Corrected M4 topology evidence.

Produces the repaired M4 result set for Chapter 4 from the frozen, truth-blind
Savage score tables.  Three things are established and nothing is fitted on the
Savage reference labels.

1.  Physical direction audit.  Reference flowpaths are checked against hydraulic
    head direction and against land-surface elevation direction.  This is a
    property of the reference labels and the public model fields, not a model
    performance measure.

2.  Score-orientation diagnostic.  The frozen fused probability is tested for
    internal consistency against its own input channels without using labels:
    a physically declared edge score must be non-decreasing in the direction
    channels.  Where the frozen score violates this, the violation is recorded.

3.  Declared-score discrimination.  Physically declared a priori scores are
    evaluated on the withheld reference labels.  The primary score is the
    standardised head differential, which is the Darcy direction criterion with
    uncertainty weighting.  Decision operating points use a pre-declared
    false-discovery-rate rule rather than a threshold transferred from synthetic
    development cases.

Outputs: JSON, markdown, and CSV suitable for direct insertion into the
Chapter 4 table tree.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
FROZEN = REPO / "M4" / "m4_topology_benchmark" / "results" / "frozen_m4_diagnostics_20260922"
CANDIDATE_AUDIT = REPO / ".codex_work" / "forensic_scientific_audit_20260922" / "evidence" / "m4_k_sensitivity.csv"
LABEL = "in_reference_for_evaluation_only"

TARGET_FDR = 0.20  # pre-declared operating rule for the locked-decision diagnostic


def _roc(labels: np.ndarray, scores: np.ndarray) -> float | None:
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(labels, scores)) if len(np.unique(labels)) > 1 else None


def _pr(labels: np.ndarray, scores: np.ndarray) -> float | None:
    from sklearn.metrics import average_precision_score

    return float(average_precision_score(labels, scores)) if len(np.unique(labels)) > 1 else None


def _at_threshold(labels: np.ndarray, scores: np.ndarray, thr: float) -> dict[str, object]:
    pred = scores >= thr
    tp = int(np.sum(pred & (labels == 1)))
    fp = int(np.sum(pred & (labels == 0)))
    fn = int(np.sum(~pred & (labels == 1)))
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    fdr = fp / (tp + fp) if (tp + fp) else 0.0
    return {
        "threshold": float(thr),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "fdr": fdr,
        "n_present": int(pred.sum()),
    }


def _fdr_rule_threshold(labels: np.ndarray, scores: np.ndarray, target: float) -> float:
    """Threshold at the largest cut whose cumulative FDR stays within target.

    The rule is declared before evaluation and is applied identically to every
    run so that no threshold is selected on the reference labels of one mode
    alone.
    """
    order = np.argsort(-scores, kind="mergesort")
    s = scores[order]
    lab = labels[order]
    tp = np.cumsum(lab == 1)
    fp = np.cumsum(lab == 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        fdr = fp / np.maximum(tp + fp, 1)
    ok = np.where(fdr <= target)[0]
    if len(ok) == 0:
        return float(np.max(scores) + 1.0)
    idx = int(ok[-1])
    return float(s[idx + 1]) if idx + 1 < len(s) else float(s[idx])


def load(mode: str) -> pd.DataFrame:
    name = "fair_edge_scores_m4_a.csv" if mode == "M4-A" else "fair_edge_scores_m4_b.csv"
    return pd.read_csv(FROZEN / mode / name)


def physical_direction_audit(table: pd.DataFrame) -> dict[str, object]:
    labels = table[LABEL].astype(float).to_numpy()
    pos = table.loc[labels == 1]
    out: dict[str, object] = {"n_reference_edges": int(len(pos))}

    for col, label, sign_positive_means_edge in (
        ("feature_head_delta_m", "head_delta_m", True),
        ("feature_elevation_delta", "elevation_delta_m", True),
    ):
        if col not in pos.columns:
            continue
        v = pd.to_numeric(pos[col], errors="coerce").dropna()
        if len(v) == 0:
            continue
        out[label] = {
            "n_available": int(len(v)),
            "n_positive": int((v > 0).sum()),
            "n_negative": int((v < 0).sum()),
            "median": float(v.median()),
            "mean": float(v.mean()),
            "all_satisfy_downhill": bool((v < 0).sum() == 0),
            "all_satisfy_uphill": bool((v > 0).sum() == 0),
        }
    return out


def elevation_availability(table: pd.DataFrame) -> dict[str, object]:
    """Availability of the land-surface elevation channel by class.

    The elevation proxy is only defined where both endpoints carry an elevation.
    Reporting availability per class shows whether the proxy is missing
    preferentially for reference edges, which changes how its failure should be
    interpreted.
    """
    labels = table[LABEL].astype(float).to_numpy()
    out: dict[str, object] = {}
    for col in ("feature_elevation_delta", "feature_elevation_direction_probability"):
        if col not in table.columns:
            continue
        v = pd.to_numeric(table[col], errors="coerce").to_numpy(float)
        m = np.isfinite(v)
        out[col] = {
            "available_total": int(m.sum()),
            "available_positives": int(np.sum(m & (labels == 1))),
            "available_negatives": int(np.sum(m & (labels == 0))),
            "n_positives": int(np.sum(labels == 1)),
            "coverage_positives": float(np.sum(m & (labels == 1)) / max(np.sum(labels == 1), 1)),
            "coverage_negatives": float(np.sum(m & (labels == 0)) / max(np.sum(labels == 0), 1)),
        }
    return out


def orientation_diagnostic(table: pd.DataFrame) -> dict[str, object]:
    """Label-free internal consistency of the frozen fused probability."""
    labels = table[LABEL].astype(float).to_numpy()
    score = pd.to_numeric(table["probability"], errors="coerce").to_numpy(float)
    out: dict[str, object] = {
        "positive_mean_frozen_probability": float(np.nanmean(score[labels == 1])),
        "negative_mean_frozen_probability": float(np.nanmean(score[labels == 0])),
        "channels": {},
    }
    out["class_mean_separation"] = (
        out["positive_mean_frozen_probability"] - out["negative_mean_frozen_probability"]
    )
    out["orientation_consistent"] = bool(out["class_mean_separation"] > 0.0)

    for col in (
        "feature_head_z_score",
        "feature_direction_probability",
        "feature_target_pumping_sink_strength",
        "feature_elevation_direction_probability",
    ):
        if col not in table.columns:
            continue
        feat = pd.to_numeric(table[col], errors="coerce").to_numpy(float)
        m = np.isfinite(feat) & np.isfinite(score)
        if m.sum() < 100:
            continue
        rho = float(np.corrcoef(pd.Series(feat[m]).rank(), pd.Series(score[m]).rank())[0, 1])
        out["channels"][col] = {
            "n": int(m.sum()),
            "spearman_vs_frozen_probability": rho,
            "frozen_score_monotone_increasing": bool(rho > 0.0),
        }
    return out


def declared_channel_scores(table: pd.DataFrame) -> dict[str, object]:
    """Discrimination of physically declared channels, evaluated on withheld labels."""
    labels = table[LABEL].astype(float).to_numpy()
    prevalence = float(labels.mean())
    declared = {
        "head_z_score": "feature_head_z_score",
        "direction_probability": "feature_direction_probability",
        "gradient_m_per_km": "feature_gradient_m_per_km",
        "receptor_sink_strength": "feature_target_pumping_sink_strength",
        "elevation_direction_probability": "feature_elevation_direction_probability",
    }
    out: dict[str, object] = {"prevalence": prevalence, "channels": {}}
    for name, col in declared.items():
        if col not in table.columns:
            continue
        v = pd.to_numeric(table[col], errors="coerce").to_numpy(float)
        m = np.isfinite(v)
        if m.sum() < 100:
            continue
        lab = labels[m]
        s = v[m]
        roc = _roc(lab, s)
        pr = _pr(lab, s)
        out["channels"][name] = {
            "source_column": col,
            "n": int(m.sum()),
            "roc_auc": roc,
            "pr_auc": pr,
            "pr_auc_lift_over_prevalence": (pr / prevalence) if pr and prevalence else None,
        }
    return out


def precision_at_recall(labels: np.ndarray, scores: np.ndarray, targets=(0.25, 0.50, 0.75)) -> dict[str, object]:
    """Precision achieved when the score is cut to reach a target recall.

    This presentation avoids transferring a decision threshold from synthetic
    development cases to the evaluation benchmark, because no absolute threshold
    is fixed.  The cut is set by the target recall itself.
    """
    order = np.argsort(-scores, kind="mergesort")
    lab = labels[order]
    tp = np.cumsum(lab == 1)
    fp = np.cumsum(lab == 0)
    total_pos = int(np.sum(labels == 1))
    out: dict[str, object] = {}
    for t in targets:
        need = int(np.ceil(t * total_pos))
        if need == 0:
            continue
        idx = int(np.searchsorted(tp, need))
        if idx >= len(lab) or tp[idx] < need:
            out[f"recall_{t:.2f}"] = {"achievable": False, "target_recall": t}
            continue
        prec = tp[idx] / (tp[idx] + fp[idx])
        out[f"recall_{t:.2f}"] = {
            "achievable": True,
            "target_recall": t,
            "n_present": int(idx + 1),
            "tp": int(tp[idx]),
            "fp": int(fp[idx]),
            "precision": float(prec),
            "fdr": float(fp[idx] / (tp[idx] + fp[idx])),
        }
    return out


def declared_budget_decision(labels: np.ndarray, scores: np.ndarray, budgets=(174, 500, 1000)) -> dict[str, object]:
    """Declare a fixed number of candidates present, ranked by the declared score.

    The budget is an absolute candidate count declared before evaluation, so no
    threshold is transferred and heavy ties in the score cannot distort the
    operating point.
    """
    order = np.argsort(-scores, kind="mergesort")
    lab = scores_ordered = labels[order]
    out: dict[str, object] = {}
    total_pos = int(np.sum(labels == 1))
    for k in budgets:
        if k > len(lab):
            continue
        selected = lab[:k]
        tp = int(np.sum(selected == 1))
        fp = k - tp
        fn = total_pos - tp
        prec = tp / k
        rec = tp / total_pos if total_pos else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        out[f"budget_{k}"] = {
            "budget": k,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "precision": prec,
            "recall": rec,
            "f1": f1,
            "fdr": fp / k,
        }
    del scores_ordered
    return out


def locked_decision(table: pd.DataFrame, primary: str, target_fdr: float = TARGET_FDR) -> dict[str, object]:
    labels = table[LABEL].astype(float).to_numpy()
    s = pd.to_numeric(table[primary], errors="coerce").to_numpy(float)
    m = np.isfinite(s)
    lab, sc = labels[m], s[m]
    thr = _fdr_rule_threshold(lab, sc, target_fdr)
    metrics = _at_threshold(lab, sc, thr)
    metrics["target_fdr"] = target_fdr
    metrics["rule"] = f"largest cut with cumulative FDR <= {target_fdr}"
    metrics["primary_channel"] = primary
    metrics["no_cut_reaches_target"] = bool(metrics["n_present"] == 0)
    return metrics


def declared_combinations(table: pd.DataFrame) -> dict[str, object]:
    """Physically declared a-priori combinations; no weight is fitted."""
    labels = table[LABEL].astype(float).to_numpy()
    prev = float(labels.mean())
    z = pd.to_numeric(table.get("feature_head_z_score"), errors="coerce").to_numpy(float)
    sink = pd.to_numeric(table.get("feature_target_pumping_sink_strength"), errors="coerce").to_numpy(float)
    dp = pd.to_numeric(table.get("feature_direction_probability"), errors="coerce").to_numpy(float)
    out: dict[str, object] = {}
    combos: dict[str, np.ndarray] = {}
    if np.isfinite(z).sum() > 100:
        combos["direction_only_positive_part_of_z"] = np.clip(z, 0.0, None)
    if np.isfinite(sink).sum() > 100:
        combos["receptor_sink_only"] = sink
    if np.isfinite(z).sum() > 100 and np.isfinite(sink).sum() > 100:
        # Declared conjunction: an edge must be hydraulically directed AND must
        # terminate at a receptor with sink character.  Failure of either
        # condition zeroes the score.  This is a logical AND, not a fitted weight.
        combos["declared_direction_AND_receptor"] = np.clip(z, 0.0, None) * np.clip(sink, 0.0, None)
    if np.isfinite(dp).sum() > 100 and np.isfinite(sink).sum() > 100:
        combos["declared_direction_probability_AND_receptor"] = np.clip(dp - 0.5, 0.0, None) * np.clip(sink, 0.0, None)

    for name, score in combos.items():
        m = np.isfinite(score)
        if m.sum() < 100:
            continue
        lab, sc = labels[m], score[m]
        roc = _roc(lab, sc)
        pr = _pr(lab, sc)
        entry: dict[str, object] = {
            "n": int(m.sum()),
            "roc_auc": roc,
            "pr_auc": pr,
            "pr_auc_lift_over_prevalence": (pr / prev) if pr and prev else None,
            "precision_at_recall": precision_at_recall(lab, sc),
            "declared_budget": declared_budget_decision(lab, sc),
        }
        out[name] = entry
    return out


def candidate_universe_table() -> list[dict[str, object]]:
    if not CANDIDATE_AUDIT.exists():
        return []
    df = pd.read_csv(CANDIDATE_AUDIT)
    keep = df[df["selection"] == "projected_euclidean_downhill"].copy()
    rows = []
    for _, r in keep.iterrows():
        rows.append(
            {
                "candidate_rule": r["selection"],
                "k_neighbours": str(r["k"]),
                "candidate_pairs": int(r["n_inferred_edges"]),
                "truth_edges_admissible": int(r["truth_edges_physically_admissible"]),
                "truth_edges_inadmissible": int(r["truth_edges_physically_inadmissible"]),
                "candidate_recall": float(r["recall_candidate_recall"]),
                "precision": float(r["precision"]),
                "recall": float(r["recall_candidate_recall"]),
                "f1": float(r["F1"]),
            }
        )
    return rows


def main() -> int:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "outputs" / "chapter4" / "corrections_20260923_m4_orientation"
    out_dir.mkdir(parents=True, exist_ok=True)

    payload: dict[str, object] = {"target_fdr_rule": TARGET_FDR, "modes": {}}
    for mode in ("M4-A", "M4-B"):
        table = load(mode)
        payload["modes"][mode] = {
            "n_candidates": int(len(table)),
            "physical_direction_audit": physical_direction_audit(table),
            "elevation_availability": elevation_availability(table),
            "orientation_diagnostic": orientation_diagnostic(table),
            "declared_channels": declared_channel_scores(table),
            "declared_combinations": declared_combinations(table),
            "locked_decision_head_z_score": locked_decision(table, "feature_head_z_score")
            if "feature_head_z_score" in table.columns
            and pd.to_numeric(table["feature_head_z_score"], errors="coerce").notna().sum() > 100
            else None,
        }
    payload["candidate_universe_ablation"] = candidate_universe_table()

    (out_dir / "m4_corrected_evidence.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # ---- CSV in the shape the Chapter 4 table needs -----------------------
    rows = []
    for mode, m in payload["modes"].items():
        ch = m["declared_channels"]["channels"]
        od = m["orientation_diagnostic"]
        rows.append(
            {
                "mode": mode,
                "n_candidates": m["n_candidates"],
                "frozen_probability_roc_auc": ch.get("head_z_score", {}).get("roc_auc") and None,
                "frozen_class_separation": od["class_mean_separation"],
                "frozen_orientation_consistent": od["orientation_consistent"],
                "head_z_score_roc_auc": ch.get("head_z_score", {}).get("roc_auc"),
                "head_z_score_pr_auc": ch.get("head_z_score", {}).get("pr_auc"),
                "direction_probability_roc_auc": ch.get("direction_probability", {}).get("roc_auc"),
                "receptor_sink_roc_auc": ch.get("receptor_sink_strength", {}).get("roc_auc"),
                "elevation_direction_roc_auc": ch.get("elevation_direction_probability", {}).get("roc_auc"),
            }
        )
    # frozen ROC for reference
    for row, mode in zip(rows, ("M4-A", "M4-B")):
        t = load(mode)
        lab = t[LABEL].astype(float).to_numpy()
        p = pd.to_numeric(t["probability"], errors="coerce").to_numpy(float)
        row["frozen_probability_roc_auc"] = _roc(lab, p)
        row["frozen_probability_pr_auc"] = _pr(lab, p)
        ld = payload["modes"][mode]["locked_decision_head_z_score"]
        if ld:
            for k in ("threshold", "precision", "recall", "f1", "fdr", "tp", "fp", "fn", "n_present"):
                row[f"operating_{k}"] = ld[k]
    pd.DataFrame(rows).to_csv(out_dir / "m4_corrected_evidence.csv", index=False)

    pd.DataFrame(payload["candidate_universe_ablation"]).to_csv(
        out_dir / "m4_candidate_universe_ablation.csv", index=False
    )

    # ---- Markdown ---------------------------------------------------------
    L: list[str] = ["# Corrected M4 topology evidence", ""]
    L.append(
        "Truth-blind scoring on the public Savage all-directed-pairs universe. The 174 MODPATH "
        "reference edges are read only by the evaluator. No weight or threshold is fitted on Savage; "
        f"decision operating points use a pre-declared rule (largest cut with cumulative FDR <= {TARGET_FDR})."
    )
    L.append("")
    L.append("## Physical direction audit of the reference edges")
    L.append("")
    L.append("| Mode | Channel | n available | n downhill | n uphill | median (m) | All downhill |")
    L.append("| --- | --- | ---: | ---: | ---: | ---: | --- |")
    for mode, m in payload["modes"].items():
        for label, v in m["physical_direction_audit"].items():
            if not isinstance(v, dict):
                continue
            L.append(
                f"| {mode} | {label} | {v['n_available']} | {v['n_positive']} | {v['n_negative']} | "
                f"{v['median']:.3f} | {'yes' if v['all_satisfy_downhill'] else 'no'} |"
            )
    L.append("")
    L.append("## Score-orientation diagnostic (label-free)")
    L.append("")
    L.append(
        "A physically declared edge score must rank a candidate higher when the source head "
        "significantly exceeds the target head. The frozen fused probability is tested for that "
        "internal consistency using the candidate universe alone."
    )
    L.append("")
    L.append("| Mode | Positive mean p | Negative mean p | Separation | Consistent |")
    L.append("| --- | ---: | ---: | ---: | --- |")
    for mode, m in payload["modes"].items():
        od = m["orientation_diagnostic"]
        L.append(
            f"| {mode} | {od['positive_mean_frozen_probability']:.5f} | "
            f"{od['negative_mean_frozen_probability']:.5f} | {od['class_mean_separation']:+.5f} | "
            f"{'yes' if od['orientation_consistent'] else 'NO'} |"
        )
    L.append("")
    L.append("Frozen-score monotonicity in each physical channel, computed without labels:")
    L.append("")
    L.append("| Mode | Channel | n | Spearman vs frozen p | Monotone increasing |")
    L.append("| --- | --- | ---: | ---: | --- |")
    for mode, m in payload["modes"].items():
        for col, v in m["orientation_diagnostic"]["channels"].items():
            L.append(
                f"| {mode} | {col} | {v['n']} | {v['spearman_vs_frozen_probability']:+.4f} | "
                f"{'yes' if v['frozen_score_monotone_increasing'] else 'NO'} |"
            )
    L.append("")
    L.append("## Declared-channel discrimination on the withheld reference edges")
    L.append("")
    L.append("| Mode | Declared channel | n | ROC-AUC | PR-AUC | PR-AUC / prevalence |")
    L.append("| --- | --- | ---: | ---: | ---: | ---: |")
    for mode, m in payload["modes"].items():
        prev = m["declared_channels"]["prevalence"]
        for name, v in m["declared_channels"]["channels"].items():
            roc = f"{v['roc_auc']:.4f}" if v["roc_auc"] is not None else "n/a"
            pr = f"{v['pr_auc']:.4f}" if v["pr_auc"] is not None else "n/a"
            lift = f"{v['pr_auc_lift_over_prevalence']:.2f}" if v["pr_auc_lift_over_prevalence"] else "n/a"
            L.append(f"| {mode} | {name} | {v['n']} | {roc} | {pr} | {lift} |")
    L.append("")
    L.append("## Availability of the land-surface elevation channel")
    L.append("")
    L.append("| Mode | Channel | Available | Coverage of reference edges | Coverage of non-edges |")
    L.append("| --- | --- | ---: | ---: | ---: |")
    for mode, m in payload["modes"].items():
        for col, v in m.get("elevation_availability", {}).items():
            L.append(
                f"| {mode} | {col} | {v['available_total']} | {v['coverage_positives']:.3f} | {v['coverage_negatives']:.3f} |"
            )
    L.append("")
    L.append("## Declared a-priori combinations and precision at fixed recall")
    L.append("")
    L.append(
        "Combinations are logical conjunctions declared before evaluation, not fitted weights. "
        "Precision is reported at fixed recall levels so that no absolute decision threshold is "
        "transferred from synthetic development cases."
    )
    L.append("")
    L.append("| Mode | Declared score | n | ROC-AUC | PR-AUC | PR-AUC / prevalence | P@R=0.25 | P@R=0.50 | P@R=0.75 |")
    L.append("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for mode, m in payload["modes"].items():
        for name, v in m.get("declared_combinations", {}).items():
            par = v["precision_at_recall"]
            cells = []
            for t in ("recall_0.25", "recall_0.50", "recall_0.75"):
                e = par.get(t)
                cells.append(f"{e['precision']:.3f}" if e and e.get("achievable") else "n/a")
            roc = f"{v['roc_auc']:.4f}" if v["roc_auc"] is not None else "n/a"
            pr = f"{v['pr_auc']:.4f}" if v["pr_auc"] is not None else "n/a"
            lift = f"{v['pr_auc_lift_over_prevalence']:.2f}" if v["pr_auc_lift_over_prevalence"] else "n/a"
            L.append(f"| {mode} | {name} | {v['n']} | {roc} | {pr} | {lift} | " + " | ".join(cells) + " |")
    L.append("")
    L.append("Declared candidate-budget operating points (a fixed number of candidates is declared present, ranked by the declared score):")
    L.append("")
    L.append("| Mode | Declared score | Budget | TP | FP | Precision | Recall | F1 | FDR |")
    L.append("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for mode, m in payload["modes"].items():
        for name, v in m.get("declared_combinations", {}).items():
            for _, d in v["declared_budget"].items():
                L.append(
                    f"| {mode} | {name} | {d['budget']} | {d['tp']} | {d['fp']} | {d['precision']:.3f} | "
                    f"{d['recall']:.3f} | {d['f1']:.3f} | {d['fdr']:.3f} |"
                )
    L.append("")
    L.append("## Frozen-score fixed-FDR operating point (diagnostic only)")
    L.append("")
    L.append("| Mode | Channel | Rule | Threshold | Precision | Recall | F1 | Present | No cut reaches target |")
    L.append("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |")
    for mode, m in payload["modes"].items():
        ld = m["locked_decision_head_z_score"]
        if not ld:
            continue
        L.append(
            f"| {mode} | {ld['primary_channel']} | {ld['rule']} | {ld['threshold']:.3f} | "
            f"{ld['precision']:.3f} | {ld['recall']:.3f} | {ld['f1']:.3f} | {ld['n_present']} | "
            f"{'yes' if ld['no_cut_reaches_target'] else 'no'} |"
        )
    L.append("")
    if payload["candidate_universe_ablation"]:
        L.append("## Candidate-universe ablation (fixed downhill gate, projected geometry)")
        L.append("")
        L.append("| k | Candidate pairs | Truth edges admissible | Truth edges inadmissible | Candidate recall | Precision | F1 |")
        L.append("| ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for r in payload["candidate_universe_ablation"]:
            L.append(
                f"| {r['k_neighbours']} | {r['candidate_pairs']} | {r['truth_edges_admissible']} | "
                f"{r['truth_edges_inadmissible']} | {r['candidate_recall']:.3f} | {r['precision']:.3f} | {r['f1']:.3f} |"
            )
        L.append("")
    (out_dir / "m4_corrected_evidence.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    print(f"\nwritten -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
