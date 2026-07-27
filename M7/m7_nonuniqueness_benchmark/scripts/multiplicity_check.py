"""Exact permutation test + Benjamini-Hochberg correction for the M7.3
case-block contrasts, run in response to the manuscript review's
Methodology issues 4-5 and Results/Discussion issue 1.

Uses only the already-locked, already-reported per-case metrics
(evidence_case_metrics.csv). No locked-test value is altered; this
performs an independent statistical re-check of the existing results.
"""

import itertools
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
metrics = pd.read_csv(ROOT / "results/m7_3_locked/evidence_case_metrics.csv")

CONTRASTS = [
    ("native", "HAC", "HC", "native_incremental_age"),
    ("native", "HAC", "HA", "native_incremental_chemistry"),
    ("native", "HAC", "AC", "native_incremental_hydraulics"),
    ("age_permuted", "HAC", "HC", "permuted_age_increment"),
    ("hydraulic_permuted", "HAC", "AC", "permuted_hydraulic_increment"),
    ("joint_misspecified", "HAC", "C", "joint_misspecification"),
]
METRIC_COLS = ["pr_auc", "brier", "log_loss", "mean_edge_entropy"]

def exact_sign_flip_p(diffs: np.ndarray) -> float:
    """Two-sided exact permutation (sign-flip) test on the mean paired difference."""
    n = len(diffs)
    observed = np.mean(diffs)
    count = 0
    total = 0
    for signs in itertools.product([1, -1], repeat=n):
        total += 1
        permuted_mean = np.mean(np.array(signs) * diffs)
        if abs(permuted_mean) >= abs(observed) - 1e-12:
            count += 1
    return count / total

rows = []
for condition, full_panel, base_panel, label in CONTRASTS:
    sub = metrics[metrics["condition"] == condition]
    full = sub[sub["panel"] == full_panel].set_index("seed")
    base = sub[sub["panel"] == base_panel].set_index("seed")
    common = sorted(set(full.index) & set(base.index))
    for metric in METRIC_COLS:
        diffs = (full.loc[common, metric] - base.loc[common, metric]).to_numpy(float)
        p = exact_sign_flip_p(diffs)
        rows.append({
            "contrast": label,
            "metric": metric,
            "n_cases": len(diffs),
            "mean_diff": float(np.mean(diffs)),
            "exact_perm_p_two_sided": p,
        })

result = pd.DataFrame(rows)

# Benjamini-Hochberg FDR correction across this predeclared family of 24 tests
p_values = result["exact_perm_p_two_sided"].to_numpy(float)
m = len(p_values)
order = np.argsort(p_values)
ranked_p = p_values[order]
bh_thresh = (np.arange(1, m + 1) / m) * 0.05
below = ranked_p <= bh_thresh
if np.any(below):
    largest_ok = np.max(np.where(below)[0])
    q_star = ranked_p[largest_ok]
else:
    q_star = -1.0
adjusted = np.empty(m)
running_min = 1.0
for i in range(m - 1, -1, -1):
    running_min = min(running_min, ranked_p[i] * m / (i + 1))
    adjusted[i] = running_min
adjusted_full = np.empty(m)
adjusted_full[order] = adjusted
result["bh_adjusted_p"] = adjusted_full
result["survives_bh_0.05"] = result["bh_adjusted_p"] <= 0.05

pd.set_option("display.width", 160)
pd.set_option("display.max_rows", 40)
print(result.round(5).to_string(index=False))
result.to_csv(ROOT / "manuscript" / "supplementary" / "multiplicity_check_results.csv", index=False)
print("\nSaved:", ROOT / "manuscript" / "supplementary" / "multiplicity_check_results.csv")
