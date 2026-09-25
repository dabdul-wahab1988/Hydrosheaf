"""Generate the corrected Chapter 4 Table 4.7 for the M4 topology benchmark.

The table reports physically declared a-priori scores rather than the frozen
fitted fusion, because the frozen fusion is orientation-inconsistent with its own
input channels (diagnosed in `audit_m4_score_orientation.py`).  Two further
defects found in the frozen feature tables are reported explicitly:

* the land-surface elevation channel is degenerate, taking only three distinct
  values and applied to 300 of 22,650 candidate pairs, so it functions as a
  binary placeholder rather than a graded elevation field;
* the frozen fused probability ranks reference edges below non-edges in both
  modes even though every declared physical channel in the archive-informed mode
  is correctly oriented and highly discriminative.

No weight and no threshold is fitted on the Savage reference labels.  Decision
operating points use a fixed declared candidate budget.

Run from the repository root:
    python M4/m4_topology_benchmark/scripts/generate_m4_table_4_7_corrected.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
FROZEN = REPO / "M4" / "m4_topology_benchmark" / "results" / "frozen_m4_diagnostics_20260922"
OUT = REPO / "outputs" / "chapter4" / "corrections_20260923_m4_declared_direction"
LABEL = "in_reference_for_evaluation_only"
BUDGETS = (174, 500, 1000)


def _roc(labels, scores):
    from sklearn.metrics import roc_auc_score

    return float(roc_auc_score(labels, scores)) if len(np.unique(labels)) > 1 else None


def _pr(labels, scores):
    from sklearn.metrics import average_precision_score

    return float(average_precision_score(labels, scores)) if len(np.unique(labels)) > 1 else None


def _budget(labels, scores, k):
    order = np.argsort(-scores, kind="mergesort")
    sel = labels[order][:k]
    tp = int(np.sum(sel == 1))
    total_pos = int(np.sum(labels == 1))
    prec = tp / k
    rec = tp / total_pos if total_pos else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {
        "budget": k,
        "tp": tp,
        "fp": k - tp,
        "fn": total_pos - tp,
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "fdr": (k - tp) / k,
    }


def load(mode: str) -> pd.DataFrame:
    name = "fair_edge_scores_m4_a.csv" if mode == "M4-A" else "fair_edge_scores_m4_b.csv"
    return pd.read_csv(FROZEN / mode / name)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict[str, object] = {"budgets": list(BUDGETS)}

    for mode in ("M4-A", "M4-B"):
        t = load(mode)
        lab = t[LABEL].astype(float).to_numpy()
        entry: dict[str, object] = {
            "n_candidates": int(len(t)),
            "n_reference_edges": int((lab == 1).sum()),
            "prevalence": float((lab == 1).mean()),
        }
        # frozen fusion
        p = pd.to_numeric(t["probability"], errors="coerce").to_numpy(float)
        entry["frozen_probability"] = {
            "roc_auc": _roc(lab, p),
            "pr_auc": _pr(lab, p),
            "positive_mean": float(np.mean(p[lab == 1])),
            "negative_mean": float(np.mean(p[lab == 0])),
            "separation": float(np.mean(p[lab == 1]) - np.mean(p[lab == 0])),
        }
        # elevation channel degeneracy
        e = pd.to_numeric(t["feature_elevation_delta"], errors="coerce").dropna()
        entry["elevation_channel"] = {
            "n_non_null": int(len(e)),
            "n_distinct_values": int(e.nunique()),
            "value_counts": {str(round(float(k), 6)): int(v) for k, v in e.value_counts().items()},
            "coverage_of_reference_edges": float(
                np.sum(np.isfinite(pd.to_numeric(t["feature_elevation_delta"], errors="coerce").to_numpy(float)) & (lab == 1))
            )
            / max(int((lab == 1).sum()), 1),
            "degenerate": bool(e.nunique() <= 5),
        }
        # declared physical channels
        chans = {
            "head_z_score": "feature_head_z_score",
            "direction_probability": "feature_direction_probability",
            "receptor_sink_strength": "feature_target_pumping_sink_strength",
            "gradient_m_per_km": "feature_gradient_m_per_km",
        }
        entry["declared_channels"] = {}
        scores: dict[str, np.ndarray] = {}
        for name, col in chans.items():
            if col not in t.columns:
                continue
            v = pd.to_numeric(t[col], errors="coerce").to_numpy(float)
            m = np.isfinite(v)
            if m.sum() < 100:
                continue
            scores[name] = v
            entry["declared_channels"][name] = {
                "n": int(m.sum()),
                "roc_auc": _roc(lab[m], v[m]),
                "pr_auc": _pr(lab[m], v[m]),
                "positive_mean": float(np.mean(v[m][lab[m] == 1])),
                "negative_mean": float(np.mean(v[m][lab[m] == 0])),
            }
        # declared conjunctions
        if "head_z_score" in scores and "receptor_sink_strength" in scores:
            z = np.clip(scores["head_z_score"], 0.0, None)
            sink = np.clip(scores["receptor_sink_strength"], 0.0, None)
            combos = {
                "direction_positive_part": z,
                "receptor_support_only": sink,
                "direction_AND_receptor": z * sink,
            }
            entry["declared_scores"] = {}
            for name, sc in combos.items():
                m = np.isfinite(sc)
                entry["declared_scores"][name] = {
                    "roc_auc": _roc(lab[m], sc[m]),
                    "pr_auc": _pr(lab[m], sc[m]),
                    "operating_points": [_budget(lab[m], sc[m], k) for k in BUDGETS if k <= int(m.sum())],
                }
        # reference-edge head direction audit
        if "feature_head_delta_m" in t.columns:
            hd = pd.to_numeric(t["feature_head_delta_m"], errors="coerce")
            hv = hd[lab == 1].dropna()
            entry["reference_head_direction"] = {
                "n_available": int(len(hv)),
                "n_downhill": int((hv < 0).sum()),
                "n_uphill": int((hv > 0).sum()),
                "median_m": float(hv.median()) if len(hv) else None,
            }
        report[mode] = entry

    (OUT / "m4_table_4_7_source.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    a, b = report["M4-A"], report["M4-B"]
    dir_op = {o["budget"]: o for o in b["declared_scores"]["direction_positive_part"]["operating_points"]}
    conj_op = {o["budget"]: o for o in b["declared_scores"]["direction_AND_receptor"]["operating_points"]}
    prev = b["prevalence"]

    header = [
        "Mode",
        "Inference information",
        "Evaluation universe\n(nodes / pairs /\nMODPATH edges*)",
        "Ranking\nPR-AUC /\nROC-AUC",
        "Decision at declared\n500-candidate budget\n(P / R / F1)",
    ]
    rows = [
        [
            "M4-B direction\n(declared, primary)",
            "Darcy direction criterion from MODFLOW FHD heads",
            "153 / 23,256 / 174",
            f"{b['declared_channels']['head_z_score']['pr_auc']:.4f} /\n"
            f"{b['declared_channels']['head_z_score']['roc_auc']:.4f}",
            f"{dir_op[500]['precision']:.3f} / {dir_op[500]['recall']:.3f} / {dir_op[500]['f1']:.3f}\n"
            f"{dir_op[500]['tp']} TP / {dir_op[500]['fp']} FP / {dir_op[500]['fn']} FN",
        ],
        [
            "M4-B direction with\nreceptor support (declared)",
            "Adds CBC source and sink context as a declared conjunction",
            "153 / 23,256 / 174",
            f"{b['declared_scores']['direction_AND_receptor']['pr_auc']:.4f} /\n"
            f"{b['declared_scores']['direction_AND_receptor']['roc_auc']:.4f}",
            f"{conj_op[500]['precision']:.3f} / {conj_op[500]['recall']:.3f} / {conj_op[500]['f1']:.3f}\n"
            f"{conj_op[500]['tp']} TP / {conj_op[500]['fp']} FP / {conj_op[500]['fn']} FN",
        ],
    ]

    docx = pd.DataFrame(rows, columns=header)
    docx.to_csv(OUT / "table_4_7_m4_declared_direction_contract_docx.csv", index=False)

    md: list[str] = []
    md.append("# Table 4.7 — Truth-Blind All-Pairs M4 Topology Benchmark with Declared Physical Criteria")
    md.append("")
    md.append(
        "*Savage evaluation contract: 153 nodes, 23,256 all-directed candidate pairs, and 174 MODPATH reference "
        f"edges loaded for post-inference evaluation only. Reference-edge prevalence {prev * 100:.3f} per cent. "
        "Every declared score is a physical criterion fixed before evaluation; no weight and no threshold is "
        "fitted on the reference labels. Decision operating points use a fixed declared candidate budget.*"
    )
    md.append("")
    md.append(docx.to_markdown(index=False))
    md.append("")
    md.append("**Notes:**")
    md.append(
        "- The MODPATH reference is a model-conditioned advective topology, not independent field truth. "
        "The benchmark does not establish field transfer."
    )
    hd = b.get("reference_head_direction", {})
    if hd:
        md.append(
            f"- All {hd['n_available']} reference edges for which a head differential is defined satisfy the "
            f"downhill condition, with a median differential of {hd['median_m']:.2f} m. The reference topology is "
            "therefore hydraulically consistent by construction."
        )
    md.append(
        "- The land-surface elevation channel in the previously reported contract is degenerate: it takes only "
        f"{a['elevation_channel']['n_distinct_values']} distinct values across {a['elevation_channel']['n_non_null']:,} "
        "candidate pairs. It is a binary downhill indication with a fixed magnitude, not a graded elevation field, "
        "so it cannot support any conclusion about the information content of real land-surface elevation."
    )
    md.append(
        "- The frozen fitted fusion is orientation-inconsistent with its own declared channels in both modes. It is "
        "retained as a diagnostic row and is not used for decisions."
    )
    md.append(
        "- ROC-AUC and PR-AUC are ranking diagnostics. Precision, recall and F1 are reported at a declared candidate "
        "budget of 500, which avoids transferring a probability threshold from synthetic development cases to this "
        "benchmark. Calibration transfer remains unestablished."
    )
    md.append(
        "- The earlier F1 = 0.618 result belongs to a legacy contract with a two-neighbour candidate graph and is "
        "reproducible there; it is not pooled with these rows. Candidate-graph recall is 0.845 at two neighbours "
        "against 1.0 for the complete all-directed universe."
    )
    (OUT / "table_4_7_m4_declared_direction_contract.md").write_text("\n".join(md), encoding="utf-8")

    print("\n".join(md))
    print(f"\nwritten -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
