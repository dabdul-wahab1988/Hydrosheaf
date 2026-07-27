"""Run the complete M6 Hydrosheaf field-transfer & robustness benchmark.

Executes six experiments on the real Ghanaian datasets and writes all results
CSVs consumed by the R figure/table scripts. Deterministic (seed=1234).
Northern Ghana chemistry comes only from the canonical raw workbook
(data/FieldData/NorthenGhana/NorthernGhana.xlsx); no aquifer/geology/
lithology metadata, prior process labels, or externally provided graph
edges are used (see DECISIONS.md).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as m6  # noqa: E402

RESULTS = m6.RESULTS_DIR
RESULTS.mkdir(parents=True, exist_ok=True)


# ============================================================ E1: readiness ===
def experiment1_readiness(data):
    print("E1 dataset readiness / missingness ...")
    tier_cap = {"northern_ghana": "tier4_full_metadata", "manu": "tier2_fluoride",
                "talensi": "tier1_isotopes"}
    all_vars = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe",
                "d18O", "d2H", "Sr_mgL", "SiO2_mgL", "Calcite_SI", "Region"]
    read_rows, avail_rows, cbe_rows = [], [], []
    for name, df in data.items():
        cbe = df.apply(m6.charge_balance_error, axis=1)
        cls = cbe.map(m6.cbe_class)
        for v, val in cbe.items():
            cbe_rows.append({"dataset": name, "cbe_pct": val, "cbe_class": cls[v]})
        read_rows.append({
            "dataset": name, "n_samples": len(df), "n_sites": df.site_id.nunique(),
            "seasonal": df.season.nunique() > 1 and "unknown" not in df.season.unique(),
            "native_tier": tier_cap[name],
            "n_quantitative": int((cls == "quantitative").sum()),
            "n_screening": int((cls == "screening").sum()),
            "n_exploratory": int((cls == "exploratory").sum()),
            "cbe_median_abs": float(np.nanmedian(np.abs(cbe))),
        })
        for v in all_vars:
            col = df[v] if v in df else pd.Series([np.nan] * len(df))
            frac = float(col.notna().mean())
            avail_rows.append({"dataset": name, "variable": v,
                               "frac_present": frac,
                               "status": "present" if frac > 0.95 else
                                         "partial" if frac > 0.05 else "absent"})
    pd.DataFrame(read_rows).to_csv(RESULTS / "m6_dataset_readiness.csv", index=False)
    pd.DataFrame(avail_rows).to_csv(RESULTS / "m6_variable_availability.csv", index=False)
    pd.DataFrame(cbe_rows).to_csv(RESULTS / "m6_cbe_distribution.csv", index=False)


# ================================================ E2: NG seasonal transfer ===
def experiment2_ng_transfer(ng, clf, cmap, rng):
    print("E2 Northern Ghana seasonal x region transfer ...")
    pairs = m6.seasonal_well_pairs(ng)
    cfg = m6.TIERS["tier4_full_metadata"]
    rows = []
    t = time.time()
    for i, (w, (wet, dry)) in enumerate(pairs.items()):
        f = m6.conservative_factor(wet, dry)
        r = m6.fit_unit(wet, dry, cfg["ions"], m6._upstream_si(dry), clf, cmap, rng,
                        n_boot=16, lifters=cfg["lifters"], evidence_row=dry)
        if r is None:
            continue
        rows.append({
            "well_id": w, "region": wet["Region"], "district": wet.get("District"),
            "latitude": wet.get("Latitude"), "longitude": wet.get("Longitude"),
            "evapo_factor": f, "dominant_family": r.dominant_family,
            "dominant_process": r.dominant_process, "resolution_class": r.resolution_class,
            "base_class": r.base_resolution_class, "mrs": r.mrs,
            "support_stability": r.support_stability, "heldout_rmse": r.heldout_rmse,
            "reactive_rmse": r.rmse, "n_support": len(r.support),
            "support": ";".join(sorted(r.support)),
        })
        if (i + 1) % 40 == 0:
            print(f"    {i+1}/{len(pairs)} ({time.time()-t:.0f}s)")
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m6_ng_field_pairs.csv", index=False)

    # per region x class support
    sup = (df.groupby(["region", "resolution_class"]).size()
           .rename("n").reset_index())
    sup.to_csv(RESULTS / "m6_ng_class_support.csv", index=False)
    agg = (df.groupby("region")
           .agg(n_wells=("well_id", "count"), mean_mrs=("mrs", "mean"),
                mean_stability=("support_stability", "mean"),
                frac_partial=("resolution_class",
                              lambda s: float((s == "partially_identifiable").mean())),
                top_family=("dominant_family", lambda s: s.mode().iat[0] if len(s) else ""))
           .reset_index())
    agg.to_csv(RESULTS / "m6_ng_aquifer_season_summary.csv", index=False)
    # family x region for figure 3
    fam = (df.groupby(["region", "dominant_family"]).size().rename("n").reset_index())
    fam.to_csv(RESULTS / "m6_ng_family_by_aquifer.csv", index=False)
    return df


# ==================================================== E3: tier ablation ======
def experiment3_tier_ablation(ng, clf, cmap, rng):
    print("E3 tier ablation (Tier 4 -> Tier 0) ...")
    pairs = m6.seasonal_well_pairs(ng)
    rows = []
    for w, (wet, dry) in pairs.items():
        per_tier = {}
        for tier in m6.TIER_ORDER:
            cfg = m6.TIERS[tier]
            si = m6._upstream_si(dry) if "si" in cfg["lifters"] else None
            r = m6.fit_unit(wet, dry, cfg["ions"], si, clf, cmap, rng, n_boot=12,
                            lifters=cfg["lifters"], evidence_row=dry, heldout=True)
            if r is None:
                continue
            per_tier[tier] = r
            rows.append({"well_id": w, "region": wet["Region"], "tier": tier,
                         "resolution_class": r.resolution_class,
                         "dominant_family": r.dominant_family, "mrs": r.mrs,
                         "support_stability": r.support_stability,
                         "n_support": len(r.support),
                         "evidence_corroborated": r.evidence_corroborated})
    abl = pd.DataFrame(rows)
    abl.to_csv(RESULTS / "m6_tier_ablation.csv", index=False)

    # transitions vs Tier 4 reference
    ref = abl[abl.tier == "tier4_full_metadata"].set_index("well_id")
    trans = []
    for tier in m6.TIER_ORDER:
        sub = abl[abl.tier == tier].set_index("well_id")
        common = ref.index.intersection(sub.index)
        class_flip = float((ref.loc[common, "resolution_class"].values
                            != sub.loc[common, "resolution_class"].values).mean())
        fam_flip = float((ref.loc[common, "dominant_family"].values
                          != sub.loc[common, "dominant_family"].values).mean())
        cls_counts = sub["resolution_class"].value_counts().to_dict()
        trans.append({
            "tier": tier, "n": len(sub),
            "frac_class_changed_vs_tier4": class_flip,
            "frac_family_changed_vs_tier4": fam_flip,
            "frac_non_identifiable": float((sub.resolution_class == "non_identifiable").mean()),
            "frac_partial": float((sub.resolution_class == "partially_identifiable").mean()),
            "mean_mrs": float(sub.mrs.mean()),
            "mean_stability": float(sub.support_stability.mean()),
            **{f"n_{k}": v for k, v in cls_counts.items()},
        })
    pd.DataFrame(trans).to_csv(RESULTS / "m6_tier_ablation_transitions.csv", index=False)
    return abl


# ================================================= E4: edge uncertainty ======
def _node_lookup(samples):
    return {r["sample_id"]: r for r in samples.to_dict("records")}


def experiment4_edge_uncertainty(ng, clf, cmap, rng, cap=140):
    print("E4 edge & flow-path uncertainty ...")
    dry = ng[ng.season == "Dry"].reset_index(drop=True)
    nodes = _node_lookup(dry)

    # Three Hydrosheaf-generated edge sets, compared against each other.
    # An earlier revision also compared against a fourth "provided_graph"
    # edge set imported from a retired antecedent-study workbook (see
    # DECISIONS.md); that borrowed topology has been removed rather than
    # substituted, so this experiment now asks how much Hydrosheaf's own
    # edge-construction choice matters, not how it compares to an external
    # graph.
    chemistry_edges = m6.chemistry_knn_edges(dry, k=3)
    edge_sets = {
        "chemistry_knn": chemistry_edges,
        "geographic_nearest": m6.geographic_edges(dry, k=3),
        "random_perturbed": m6.random_edges(dry, n=len(chemistry_edges), rng=rng),
    }
    cfg = m6.TIERS["tier4_full_metadata"]
    rows = []
    for setname, edges in edge_sets.items():
        edges = list(edges)
        rng.shuffle(edges)
        edges = edges[:cap]
        for (s, t) in edges:
            src, tgt = nodes.get(s), nodes.get(t)
            if src is None or tgt is None:
                continue
            si = m6._upstream_si(tgt)
            r = m6.fit_unit(src, tgt, cfg["ions"], si, clf, cmap, rng, n_boot=8,
                            lifters=cfg["lifters"], evidence_row=tgt, heldout=False)
            if r is None:
                continue
            rows.append({"edge_set": setname, "source": s, "target": t,
                         "dominant_family": r.dominant_family,
                         "resolution_class": r.resolution_class, "mrs": r.mrs,
                         "support_stability": r.support_stability})
    ed = pd.DataFrame(rows)
    ed.to_csv(RESULTS / "m6_edge_sensitivity.csv", index=False)

    # network summary: family distribution + divergence vs the chemistry-knn
    # edge set (the most evidence-based of the three Hydrosheaf-generated
    # sets, used as the reference now that the borrowed provided-graph
    # comparator has been removed).
    fam = (ed.groupby(["edge_set", "dominant_family"]).size().rename("n").reset_index())
    fam["frac"] = fam.groupby("edge_set")["n"].transform(lambda s: s / s.sum())
    fam.to_csv(RESULTS / "m6_edge_family_distribution.csv", index=False)

    ref = fam[fam.edge_set == "chemistry_knn"].set_index("dominant_family")["frac"]
    summ = []
    for setname in edge_sets:
        cur = fam[fam.edge_set == setname].set_index("dominant_family")["frac"]
        idx = ref.index.union(cur.index)
        tvd = 0.5 * float(np.abs(ref.reindex(idx, fill_value=0)
                                 - cur.reindex(idx, fill_value=0)).sum())
        sub = ed[ed.edge_set == setname]
        summ.append({"edge_set": setname, "n_edges": len(sub),
                     "tvd_vs_chemistry_knn": tvd,
                     "mean_mrs": float(sub.mrs.mean()),
                     "mean_stability": float(sub.support_stability.mean()),
                     "frac_partial": float((sub.resolution_class == "partially_identifiable").mean())})
    pd.DataFrame(summ).to_csv(RESULTS / "m6_edge_network_summary.csv", index=False)
    return ed


# ================================================= E5: external transfer =====
def experiment5_external(data, ng_ref, clf, cmap, rng):
    print("E5 external sparse transfer (Talensi, Manu) ...")
    # Both external datasets measure pH, temperature and major ions, so
    # PHREEQC saturation indices are computed for them the same way as for
    # Northern Ghana (m6.load_talensi/load_manu) and used here as an
    # additional "si" evidence lifter; SI is not gated behind Sr/SiO2
    # availability, unlike the Northern Ghana Tier 3/4 ladder.
    specs = {"talensi": ("tier1_isotopes", ["iso", "si"]),
             "manu": ("tier2_fluoride", ["iso", "si"])}
    rows = []
    for name, (tier, lifters) in specs.items():
        df = data[name].reset_index(drop=True)
        edges = m6.chemistry_knn_edges(df, k=3)
        nodes = _node_lookup(df)
        cfg = m6.TIERS[tier]
        for (s, t) in edges:
            src, tgt = nodes.get(s), nodes.get(t)
            if src is None or tgt is None:
                continue
            r = m6.fit_unit(src, tgt, cfg["ions"], m6._upstream_si(tgt), clf, cmap, rng,
                            n_boot=16, lifters=lifters, evidence_row=tgt)
            if r is None:
                continue
            rows.append({"dataset": name, "native_tier": tier, "source": s, "target": t,
                         "dominant_family": r.dominant_family,
                         "dominant_process": r.dominant_process,
                         "resolution_class": r.resolution_class, "mrs": r.mrs,
                         "support_stability": r.support_stability,
                         "reactive_rmse": r.rmse, "n_support": len(r.support)})
    ext = pd.DataFrame(rows)
    ext.to_csv(RESULTS / "m6_external_transfer.csv", index=False)

    # NG reference at matched tiers (re-run NG edges at tier1 & tier2 for fair comparison)
    dry = ng_ref[ng_ref.season == "Dry"].reset_index(drop=True)
    ng_edges = m6.chemistry_knn_edges(dry, k=3)
    rng.shuffle(ng_edges)
    ng_nodes = _node_lookup(dry)
    ng_rows = []
    for tier, lifters in [("tier1_isotopes", ["iso", "si"]), ("tier2_fluoride", ["iso", "si"])]:
        cfg = m6.TIERS[tier]
        for (s, t) in ng_edges[:120]:
            src, tgt = ng_nodes.get(s), ng_nodes.get(t)
            if src is None or tgt is None:
                continue
            r = m6.fit_unit(src, tgt, cfg["ions"], m6._upstream_si(tgt), clf, cmap, rng,
                            n_boot=8, lifters=lifters, evidence_row=tgt, heldout=False)
            if r:
                ng_rows.append({"dataset": "northern_ghana_ref", "native_tier": tier,
                                "mrs": r.mrs, "support_stability": r.support_stability,
                                "resolution_class": r.resolution_class})
    ng_ref_df = pd.DataFrame(ng_rows)

    def _summ(df, label):
        return {"dataset": label, "n_edges": len(df),
                "mean_mrs": float(df.mrs.mean()),
                "mean_stability": float(df.support_stability.mean()),
                "frac_partial": float((df.resolution_class == "partially_identifiable").mean()),
                "frac_non_identifiable": float((df.resolution_class == "non_identifiable").mean())}
    summ = [_summ(ext[ext.dataset == "talensi"], "talensi (Tier1)"),
            _summ(ext[ext.dataset == "manu"], "manu (Tier2)"),
            _summ(ng_ref_df[ng_ref_df.native_tier == "tier1_isotopes"], "N.Ghana ref (Tier1)"),
            _summ(ng_ref_df[ng_ref_df.native_tier == "tier2_fluoride"], "N.Ghana ref (Tier2)")]
    pd.DataFrame(summ).to_csv(RESULTS / "m6_external_summary.csv", index=False)
    return ext


# ================================================= E6: limitation map ========
def experiment6_limitation_map(ng_pairs, ext, edge, cmap):
    print("E6 limitation map & conservative-vs-bestfit ...")
    frames = []
    frames.append(ng_pairs.assign(dataset="northern_ghana")[
        ["dataset", "dominant_family", "resolution_class"]])
    frames.append(ext[["dataset", "dominant_family", "resolution_class"]])
    allu = pd.concat(frames, ignore_index=True)
    lim = (allu.groupby(["dataset", "dominant_family", "resolution_class"]).size()
           .rename("n").reset_index())
    lim["frac"] = lim.groupby(["dataset", "dominant_family"])["n"].transform(lambda s: s / s.sum())
    lim.to_csv(RESULTS / "m6_limitation_map.csv", index=False)

    # conservative (class-level) vs best-fit (single reaction) commitment
    rows = []
    for _, r in ng_pairs.iterrows():
        support = set(str(r["support"]).split(";")) if isinstance(r["support"], str) and r["support"] else set()
        classes = {cmap[s] for s in support if s in cmap}
        # ambiguous members admitted by conservative reporting
        admitted = sum(sum(1 for v in cmap.values() if v == c) for c in classes)
        rows.append({"dataset": "northern_ghana", "well_id": r["well_id"],
                     "bestfit_n_reactions": len(support),
                     "conservative_n_admitted": admitted,
                     "resolution_class": r["resolution_class"]})
    cvb = pd.DataFrame(rows)
    cvb.to_csv(RESULTS / "m6_conservative_vs_bestfit.csv", index=False)


def main():
    t0 = time.time()
    data = m6.load_all()
    ng = data["northern_ghana"]
    clf = m6.TransferClassifier()
    cmap = m6.get_class_map()
    rng = np.random.default_rng(m6.SEED)

    experiment1_readiness(data)
    ng_pairs = experiment2_ng_transfer(ng, clf, cmap, rng)
    experiment3_tier_ablation(ng, clf, cmap, rng)
    edge = experiment4_edge_uncertainty(ng, clf, cmap, rng)
    ext = experiment5_external(data, ng, clf, cmap, rng)
    experiment6_limitation_map(ng_pairs, ext, edge, cmap)
    print(f"\nM6 field-transfer analysis complete in {time.time()-t0:.0f}s.")
    print("Results ->", RESULTS)


if __name__ == "__main__":
    main()
