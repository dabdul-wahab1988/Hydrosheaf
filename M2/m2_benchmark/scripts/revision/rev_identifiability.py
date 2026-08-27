"""Revision analysis 1 (CAGEO-D-26-00847): inverse-model identifiability diagnostics.

Answers Reviewer 2 Major 3 (rank, condition number, reaction correlation,
sensitivity to the reaction dictionary) and supports the honest reframing of
reaction-extent recovery (Reviewer 1 Major 10/15, Reviewer 2 Minor 1).

Outputs (M2/m2_benchmark/results/revision/):
  - identifiability_dictionary.csv     rank, condition number, collinearity pairs
  - reaction_recovery_summary.csv      extent-level recovery metrics (locked canonical)
  - reaction_correlations.csv          empirical Pearson correlations of recovered extents
  - family_recovery.csv                process-family level recovery metrics
  - dictionary_sensitivity.csv         leave-one-out dictionary sensitivity (sample refits)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
BENCHMARK_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(BENCHMARK_ROOT / "scripts"))

from run_m2_benchmark import (  # noqa: E402
    add_reactions,
    build_noiseless_nodes,
    df_to_samples,
    edge_inputs_for_variant,
    load_truth,
    make_config,
    make_realisations,
    node_vectors,
    reaction_maps,
)
from hydrosheaf import fit_network  # noqa: E402

OUT_DIR = BENCHMARK_ROOT / "results" / "revision"

FAMILY = {
    "calcite": "carbonate", "dolomite": "carbonate",
    "gypsum": "evaporite", "halite": "evaporite",
    "albite": "silicate", "anorthite": "silicate",
    "pyrite_net": "redox", "pyrite_oxidation_aerobic": "redox",
    "sulfate_reduction": "redox", "iron_reduction": "redox", "denit": "redox",
    "NO3src": "nitrate",
    "CaNa_exch": "exchange", "NaCa_exch": "exchange",
    "MgNa_exch": "exchange", "NaMg_exch": "exchange",
}


def corr_r2(x: pd.Series, y: pd.Series) -> float:
    p = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")})
    p = p.replace([np.inf, -np.inf], np.nan).dropna()
    if len(p) < 2 or p["x"].std(ddof=0) <= 1e-12 or p["y"].std(ddof=0) <= 1e-12:
        return float("nan")
    return float(np.corrcoef(p["x"], p["y"])[0, 1] ** 2)


def dictionary_diagnostics(config, labels, matrix) -> pd.DataFrame:
    rows = []
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if j <= i:
                continue
            cos = abs(float(np.dot(matrix[i], matrix[j]))) / (
                float(np.linalg.norm(matrix[i]) * np.linalg.norm(matrix[j])) + 1e-15
            )
            if cos > 0.7:
                rows.append({"reaction_a": a, "reaction_b": b, "abs_cosine": round(cos, 3)})
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    truth = load_truth(BENCHMARK_ROOT / "config" / "ground_truth.yaml")
    config = make_config(truth)
    nodes = build_noiseless_nodes(truth, config)
    config = make_config(truth, node_vectors(nodes, truth["ion_order"]))
    rmap, labels = reaction_maps(config)
    matrix = np.array(list(rmap.values()), dtype=float)

    # ---- 1. dictionary diagnostics -------------------------------------
    rank = int(np.linalg.matrix_rank(matrix))
    cond = float(np.linalg.cond(matrix))
    dict_rows = [
        {"metric": "n_reactions", "value": matrix.shape[0]},
        {"metric": "n_ions", "value": matrix.shape[1]},
        {"metric": "rank", "value": rank},
        {"metric": "rank_deficiency", "value": matrix.shape[1] - rank},
        {"metric": "condition_number", "value": cond},
        {"metric": "n_collinear_pairs_abs_cos_gt_0.7", "value": len(dictionary_diagnostics(config, labels, matrix))},
    ]
    pd.DataFrame(dict_rows).to_csv(OUT_DIR / "identifiability_dictionary.csv", index=False)
    dictionary_diagnostics(config, labels, matrix).to_csv(
        OUT_DIR / "identifiability_dictionary.csv".replace(".csv", "_pairs.csv"), index=False
    )

    # ---- 2. extent-level recovery summary (canonical, from fresh run) ---
    rr = pd.read_csv(BENCHMARK_ROOT / "results" / "reaction_recovery.csv")
    active = rr[rr["true_extent_mmolL"].abs() > 0.01]
    inactive = rr[rr["true_extent_mmolL"].abs() <= 0.01]
    fp_rate = float((inactive["recovered_extent_mmolL"].abs() > 0.05).mean())
    summary = pd.DataFrame(
        [
            {"metric": "n_active_rows", "value": len(active)},
            {"metric": "n_inactive_rows", "value": len(inactive)},
            {"metric": "active_r2_corr", "value": corr_r2(active["true_extent_mmolL"], active["recovered_extent_mmolL"])},
            {"metric": "active_mae_mmolL", "value": float(np.mean(np.abs(active["true_extent_mmolL"] - active["recovered_extent_mmolL"])))},
            {"metric": "active_rmse_mmolL", "value": float(np.sqrt(np.mean((active["true_extent_mmolL"] - active["recovered_extent_mmolL"]) ** 2)))},
            {"metric": "false_positive_rate_0.05", "value": fp_rate},
        ]
    )
    summary.to_csv(OUT_DIR / "reaction_recovery_summary.csv", index=False)

    # ---- 3. empirical reaction correlation ------------------------------
    wide = rr.pivot_table(index=["realisation", "edge_id"], columns="reaction_label", values="recovered_extent_mmolL")
    corr = wide.corr()
    corr.to_csv(OUT_DIR / "reaction_correlations.csv")

    # ---- 4. family-level recovery ---------------------------------------
    rr["family"] = rr["reaction_label"].map(FAMILY)
    fam = rr.groupby(["realisation", "edge_id", "family"], as_index=False).agg(
        true=("true_extent_mmolL", "sum"), inf=("recovered_extent_mmolL", "sum")
    )
    fam_active = fam[fam["true"].abs() > 0.01]
    fam["abs_true"] = fam["true"].abs()
    dominant = fam.loc[fam.groupby(["realisation", "edge_id"])["abs_true"].idxmax()]
    dominant_hit = float(((np.sign(dominant["true"]) == np.sign(dominant["inf"])) & (dominant["abs_true"] > 0.01)).mean())
    fam_rows = pd.DataFrame(
        [
            {"metric": "n_active_family_rows", "value": len(fam_active)},
            {"metric": "family_r2_corr", "value": corr_r2(fam_active["true"], fam_active["inf"])},
            {"metric": "family_mae_mmolL", "value": float(np.mean(np.abs(fam_active["true"] - fam_active["inf"])))},
            {"metric": "family_sign_match_rate", "value": float((np.sign(fam_active["true"]) == np.sign(fam_active["inf"])).mean())},
            {"metric": "dominant_family_hit_rate", "value": dominant_hit},
        ]
    )
    fam_rows.to_csv(OUT_DIR / "family_recovery.csv", index=False)

    # ---- 5. dictionary sensitivity (leave-one-out, sample of realisations)
    n_sample = 5
    mineral_labels = [lbl for lbl in labels if lbl in truth["active_minerals"]]
    variants = {"full_dictionary": set(truth["active_minerals"])}
    for m in mineral_labels:
        variants[f"minus_{m}"] = {x for x in truth["active_minerals"] if x != m}
    sens_rows = []
    realisations_df = make_realisations(truth, nodes, n_sample)
    for variant, minerals in variants.items():
        cfg = make_config(truth, node_vectors(nodes, truth["ion_order"]))
        cfg.active_minerals = list(minerals)
        recs = []
        for realisation in range(n_sample):
            realised = realisations_df[realisations_df["realisation"] == realisation].copy()
            samples = df_to_samples(realised, truth["ion_order"])
            edges = edge_inputs_for_variant(truth, "full", realisation, realised)
            results = fit_network(samples, edges, cfg)
            for res in results:
                zmap = dict(zip(res.z_labels, res.z_extents))
                for lbl, val in zmap.items():
                    recs.append({"realisation": realisation, "edge_id": res.edge_id, "reaction": lbl, "recovered": float(val)})
        rec_df = pd.DataFrame(recs)
        truth_edges = pd.read_csv(BENCHMARK_ROOT / "data" / "ground_truth_edges.csv")
        true_map = {}
        for _, e in truth_edges.iterrows():
            for k, v in e.items():
                if isinstance(v, str) and v.startswith("reaction_") and v.endswith("_true"):
                    pass
        # build true extents from truth spec
        true_rows = []
        specs = {}
        for edge in truth["generation_edges"]:
            specs[edge["edge_id"]] = edge
        for _, r in rec_df.iterrows():
            spec = specs.get(r["edge_id"], {})
            true_rows.append({"edge_id": r["edge_id"], "reaction": r["reaction"], "true": float(spec.get("reactions", {}).get(r["reaction"], 0.0))})
        merged = rec_df.merge(pd.DataFrame(true_rows), on=["edge_id", "reaction"])
        act = merged[merged["true"].abs() > 0.01]
        if len(act):
            sens_rows.append(
                {
                    "dictionary_variant": variant,
                    "n_active_rows": len(act),
                    "active_r2_corr": corr_r2(act["true"], act["recovered"]),
                    "active_mae_mmolL": float(np.mean(np.abs(act["true"] - act["recovered"]))),
                }
            )
        else:
            sens_rows.append({"dictionary_variant": variant, "n_active_rows": 0, "active_r2_corr": np.nan, "active_mae_mmolL": np.nan})
    pd.DataFrame(sens_rows).to_csv(OUT_DIR / "dictionary_sensitivity.csv", index=False)
    print("Wrote revision identifiability diagnostics to", OUT_DIR)


if __name__ == "__main__":
    main()
