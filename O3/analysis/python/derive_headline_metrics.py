"""TAB-2 / FIG-3: independent, uncalibrated capture-type vs correctness-type
metrics for the three inference layers, on each component's own native scale.

The three components do not share one statistic. Topology is a literal edge
classification (precision/recall). Reaction reports macro-averaged Phase F1
and false-discovery rate per scenario (Method 1) and, separately, per-reaction
confusion counts that can be pooled into a micro-averaged precision/recall
(Method 2); both are reported, and the difference between them is disclosed
rather than hidden, following the M2.3 D5 precedent for discrepancy
reporting. Age has no classification truth at all; its two axes are the
fraction of attempted fits that clear the reference-free reportability guard
("capture", i.e. how much of the attempted inference space yields an
actionable claim), and the within-factor-2 agreement of reported estimates on
an independent held-out cross-validation split ("correctness given
reported"). These are stated as a disclosed analogy in Methods, not a claim
that the three numbers are the same statistic.

Run:  .venv/Scripts/python.exe O3/analysis/python/derive_headline_metrics.py
"""

from __future__ import annotations

import pandas as pd

from _common import M3, M4, M5, OUT, wilson_ci, write

OUT_DISC = OUT / "evidence_discrepancies.csv"


def topology_metrics() -> pd.DataFrame:
    df = pd.read_csv(M4 / "tables/Manuscript_Ready"
                      / "Main_Table2_Topology_Performance_Failure_Summary.csv")
    row = df.loc[df["scenario"] == "Head gradient"].iloc[0]
    floor = df.loc[df["scenario"] == "Sink-aware baseline"].iloc[0]
    # Recover TP from the reported precision/recall/FP/FN (consistent to
    # rounding): precision = TP/(TP+FP), recall = TP/(TP+FN).
    fp, fn = float(row["false_positives"]), float(row["false_negatives"])
    tp = round(float(row["precision"]) * fp / (1 - float(row["precision"])))
    n_precision, n_recall = tp + fp, tp + fn
    prec_ci = wilson_ci(int(tp), int(n_precision))
    rec_ci = wilson_ci(int(tp), int(n_recall))
    return pd.DataFrame([
        dict(component="Topology", axis="capture", metric="recall",
             value=float(row["recall"]), ci_low=rec_ci[0], ci_high=rec_ci[1],
             n=int(n_recall), source="Main_Table2, Head gradient row",
             note="Fraction of true MODPATH edges recovered; 95% Wilson CI "
                  f"from {int(tp)}/{int(n_recall)} true edges"),
        dict(component="Topology", axis="correctness", metric="precision",
             value=float(row["precision"]), ci_low=prec_ci[0], ci_high=prec_ci[1],
             n=int(n_precision), source="Main_Table2, Head gradient row",
             note="Fraction of inferred edges that are true MODPATH edges; "
                  f"95% Wilson CI from {int(tp)}/{int(n_precision)} inferred edges"),
        dict(component="Topology", axis="informed_floor", metric="F1",
             value=float(floor["f1"]), ci_low=float("nan"), ci_high=float("nan"),
             n=pd.NA, source="Main_Table2, Sink-aware baseline row",
             note="Receptor-set-only floor, zero hydraulic information, not "
                  "independent evidence"),
    ])


def age_metrics() -> pd.DataFrame:
    design = pd.read_csv(M3 / "results/m3_design_matrix_summary.csv")
    strict = design.loc[design["scenario_id"] == "tracerlpm_strict_parity"].iloc[0]
    n_reportable, n_total = int(strict["identifiable_rows"]), int(strict["total_rows"])
    coverage = n_reportable / n_total
    coverage_ci = wilson_ci(n_reportable, n_total)
    held = pd.read_csv(M3 / "results/m3_usgs_calibrated_parity_summary.csv")
    unc = held.loc[held["mode"] == "uncalibrated_strict_parity_on_heldout_folds"].iloc[0]
    n_held = int(unc["N"])
    k_held = round(float(unc["within_factor_2"]) * n_held)
    wf2_ci = wilson_ci(k_held, n_held)
    return pd.DataFrame([
        dict(component="Age / residence time", axis="capture", metric="reportability_rate",
             value=coverage, ci_low=coverage_ci[0], ci_high=coverage_ci[1], n=n_total,
             source="m3_design_matrix_summary.csv, tracerlpm_strict_parity",
             note=f"{n_reportable} of {n_total} "
                  "fitted rows cleared the reference-free reportability guard; "
                  "95% Wilson CI"),
        dict(component="Age / residence time", axis="correctness", metric="within_factor_2",
             value=float(unc["within_factor_2"]), ci_low=wf2_ci[0], ci_high=wf2_ci[1],
             n=n_held,
             source="m3_usgs_calibrated_parity_summary.csv, held-out uncalibrated",
             note=f"n = {n_held} held-out cross-validated pairs, "
                  "independent (D4 primary result); 95% Wilson CI"),
        dict(component="Age / residence time", axis="correctness_r2", metric="log10_r2",
             value=float(unc["log10_r2"]), ci_low=float("nan"), ci_high=float("nan"),
             n=n_held, source="m3_usgs_calibrated_parity_summary.csv",
             note="Held-out uncalibrated log10 R2, independent; no closed-form "
                  "CI computed for R2 itself"),
    ])


def reaction_metrics() -> pd.DataFrame:
    t1 = pd.read_csv(M5 / "tables/table1_comparative_inverse_performance.csv")
    t1.columns = [c.strip() for c in t1.columns]
    guarded = t1.loc[t1["Method"] == "Hydrosheaf Guarded"].iloc[0]

    def bounds(cell: str) -> tuple[float, float, float]:
        # "0.361 [0.319, 0.402]" -> (mean, low, high)
        mean_s, rest = cell.split("[")
        lo_s, hi_s = rest.rstrip("]").split(",")
        return float(mean_s.strip()), float(lo_s.strip()), float(hi_s.strip())

    fdr_mean, fdr_lo, fdr_hi = bounds(guarded["False-discovery rate, mean [95% CI]"])
    macro_precision = 1.0 - fdr_mean
    # Precision = 1 - FDR, so the CI bounds invert and swap.
    prec_ci = (round(1 - fdr_hi, 3), round(1 - fdr_lo, 3))

    s6 = pd.read_csv(M5 / "tables/tableS6_reaction_confusion_matrices.csv")
    pooled = s6.groupby("Method")[["TP", "FP", "FN"]].sum()
    p = pooled.loc["hydrosheaf_guarded"]
    micro_precision = p["TP"] / (p["TP"] + p["FP"])
    micro_recall = p["TP"] / (p["TP"] + p["FN"])
    recall_ci = wilson_ci(int(p["TP"]), int(p["TP"] + p["FN"]))
    micro_precision_ci = wilson_ci(int(p["TP"]), int(p["TP"] + p["FP"]))

    rows = [
        dict(component="Reaction", axis="capture", metric="recall_pooled",
             value=float(micro_recall), ci_low=recall_ci[0], ci_high=recall_ci[1],
             n=int(p["TP"] + p["FN"]),
             source="tableS6_reaction_confusion_matrices.csv, hydrosheaf_guarded, "
                    "pooled TP/FP/FN across 16 reactions",
             note="Micro-averaged recall; not identical to any single "
                  "already-published headline number, reported for the "
                  "capture axis because no macro-averaged recall is published; "
                  "95% Wilson CI on the pooled count"),
        dict(component="Reaction", axis="correctness", metric="precision_macro",
             value=float(macro_precision), ci_low=prec_ci[0], ci_high=prec_ci[1],
             n=pd.NA,
             source="table1_comparative_inverse_performance.csv, Hydrosheaf Guarded, "
                    "1 - mean false-discovery rate",
             note="Macro-averaged across 240 scenarios; this is M5's own "
                  "published headline correctness axis; CI inverted from the "
                  "published false-discovery-rate 95% CI"),
        dict(component="Reaction", axis="correctness_pooled_crosscheck",
             metric="precision_pooled", value=float(micro_precision),
             ci_low=micro_precision_ci[0], ci_high=micro_precision_ci[1],
             n=int(p["TP"] + p["FP"]),
             source="tableS6, pooled TP/FP/FN",
             note="Cross-check only; differs from precision_macro because "
                  "macro (per-scenario mean) and micro (pooled count) "
                  "averaging are different, disclosed conventions"),
    ]
    return pd.DataFrame(rows)


def main() -> None:
    df = pd.concat([topology_metrics(), age_metrics(), reaction_metrics()],
                    ignore_index=True)
    write(df, "headline_metrics.csv")

    # Discrepancy disclosure between the two M5 precision conventions (D2/D5-style).
    disc = pd.DataFrame([dict(
        quantity="Reaction benchmark Hydrosheaf Guarded precision (reaction phases)",
        value_a="0.639 (macro mean across 240 scenarios; the reaction benchmark's own table 1)",
        value_b="0.586 (micro, pooled TP/FP/FN across 16 reactions; tableS6)",
        resolution="Both are correct under their own stated averaging "
                   "convention; O3 reports the macro value as the primary "
                   "correctness axis and the pooled value as a cross-check.")])
    disc.to_csv(OUT_DISC, index=False)
    print(f"wrote {OUT_DISC}")


if __name__ == "__main__":
    main()
