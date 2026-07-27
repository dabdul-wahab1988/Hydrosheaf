"""M6 validation and robustness diagnostics for the Q1 analysis package.

A. Circularity: field tier-ablation with the evidence-lift gate ON vs OFF, plus
   gate-independent structural diagnostics per tier -> shows how much of the
   'collapse' is the gate (a conservative reporting prior) versus real structure.
B. Transport sensitivity: headline field inference under no / Cl / recharge-
   endmember transport treatments.
C. Talensi charge-balance diagnosis: what drives the ~30% imbalance.
D. Discrimination: field MRS spread and synthetic true-class spread.
E. Bootstrap 95% CIs on the tier non-identifiable fractions.
Deterministic (seed 1234).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as m6
import synthetic_reaction_truth_model as rx

RES = m6.RESULTS_DIR
RES.mkdir(parents=True, exist_ok=True)

# M6 evidence tier -> reaction-basis panel for structural diagnostics
TIER_TO_PANEL = {"tier0_majors": "majors", "tier1_isotopes": "majors",
                 "tier2_fluoride": "plus_F", "tier3_sr_sio2": "plus_Sr_SiO2",
                 "tier4_full_metadata": "plus_Sr_SiO2"}


# ---- A + E: gate on/off, structural, CIs ------------------------------------
def gate_and_structural(ng, clf, cmap, rng):
    print("A. gate on/off + structural + CIs ...")
    pairs = m6.seasonal_well_pairs(ng)
    rows = []
    for w, (wet, dry) in pairs.items():
        for tier in m6.TIER_ORDER:
            cfg = m6.TIERS[tier]
            si = m6._upstream_si(dry) if "si" in cfg["lifters"] else None
            r = m6.fit_unit(wet, dry, cfg["ions"], si, clf, cmap, rng, n_boot=10,
                            lifters=cfg["lifters"], evidence_row=dry)
            if r is None:
                continue
            rows.append({"well": w, "tier": tier,
                         "gated_class": r.resolution_class,
                         "ungated_class": r.base_resolution_class,
                         "dominant_family": r.dominant_family})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "m6_field_gate_perwell.csv", index=False)

    B = 400
    out = []
    for tier in m6.TIER_ORDER:
        sub = df[df.tier == tier]
        gated = (sub.gated_class == "non_identifiable").mean()
        ungated = (sub.ungated_class == "non_identifiable").mean()
        # bootstrap CI on gated fraction (resample wells)
        vals = (sub.gated_class == "non_identifiable").to_numpy().astype(float)
        boot = [np.mean(rng.choice(vals, len(vals), replace=True)) for _ in range(B)]
        lo, hi = np.percentile(boot, [2.5, 97.5])
        st = rx.diagnostics(rx.PANELS[TIER_TO_PANEL[tier]])
        out.append({"tier": tier, "gated_frac_non_identifiable": float(gated),
                    "gated_ci_lo": float(lo), "gated_ci_hi": float(hi),
                    "ungated_frac_non_identifiable": float(ungated),
                    "rank": st["rank"], "nullity": st["nullity"],
                    "silicate_coherence": rx.silicate_identifiability(rx.PANELS[TIER_TO_PANEL[tier]]),
                    "n_equivalence_classes": st["n_equivalence_classes"]})
    pd.DataFrame(out).to_csv(RES / "m6_field_gate_structural.csv", index=False)
    return df


# ---- B: transport sensitivity ------------------------------------------------
def transport_sensitivity(ng, clf, cmap, rng):
    print("B. transport sensitivity ...")
    pairs = m6.seasonal_well_pairs(ng)
    cfg = m6.TIERS["tier4_full_metadata"]
    # per-region recharge endmember = minimum ionic-load sample (no
    # independent aquifer-type classification exists for these boreholes;
    # see DECISIONS.md)
    ng2 = ng.copy()
    ng2["load"] = ng2.apply(m6._ionic_load, axis=1)
    recharge = {region: g.loc[g["load"].idxmin()].to_dict()
                for region, g in ng2.groupby("Region")}
    rows = []
    for treat in ["none", "cl", "recharge"]:
        for w, (wet, dry) in pairs.items():
            src = recharge[dry["Region"]] if treat == "recharge" else wet
            tcorr = treat != "none"
            r = m6.fit_unit(src, dry, cfg["ions"], m6._upstream_si(dry), clf, cmap, rng,
                            n_boot=8, lifters=cfg["lifters"], evidence_row=dry,
                            heldout=False, transport_correct=tcorr)
            if r is None:
                continue
            rows.append({"treatment": treat, "well": w,
                         "dominant_family": r.dominant_family,
                         "resolution_class": r.resolution_class, "mrs": r.mrs})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "m6_transport_sensitivity_perwell.csv", index=False)
    summ = (df.groupby(["treatment", "dominant_family"]).size().rename("n").reset_index())
    summ["frac"] = summ.groupby("treatment")["n"].transform(lambda s: s / s.sum())
    summ.to_csv(RES / "m6_transport_sensitivity.csv", index=False)
    # agreement of dominant family vs Cl treatment (reference)
    ref = df[df.treatment == "cl"].set_index("well")["dominant_family"]
    agree = []
    for treat in ["none", "cl", "recharge"]:
        cur = df[df.treatment == treat].set_index("well")["dominant_family"]
        common = ref.index.intersection(cur.index)
        agree.append({"treatment": treat,
                      "dominant_family_agreement_vs_cl": float((ref.loc[common] == cur.loc[common]).mean()),
                      "mean_mrs": float(df[df.treatment == treat]["mrs"].mean())})
    pd.DataFrame(agree).to_csv(RES / "m6_transport_agreement.csv", index=False)


# ---- C: Talensi CBE diagnosis ------------------------------------------------
def talensi_diagnosis():
    print("C. Talensi CBE diagnosis ...")
    df = m6.load_talensi()
    def cbe_of(row):
        return m6.charge_balance_error(row)
    base = df.apply(cbe_of, axis=1)
    cat = df.apply(lambda r: sum(max(r.get(i, 0) or 0, 0) * m6.CHARGE[i]
                                 for i in ["Ca", "Mg", "Na", "K"]), axis=1)
    an = df.apply(lambda r: sum(abs(r.get(i, 0) or 0) * abs(m6.CHARGE[i])
                  for i in ["HCO3", "Cl", "SO4", "NO3"] if pd.notna(r.get(i))), axis=1)
    # HCO3 is already in mmol/L (== meq/L for a monovalent anion) after harmonisation
    hco3_meq = df["HCO3"]
    # If alkalinity had been reported as mg/L CaCO3 (not HCO3), the true meq would be
    # LARGER (61.0/50.0x), which worsens the anion excess -> rules that hypothesis out.
    an_if_caco3 = an - hco3_meq + hco3_meq * (61.0168 / 50.04)
    cbe_caco3 = 100 * (cat - an_if_caco3) / (cat + an_if_caco3)
    missing_cation_meq = (an - cat)  # cation deficit that would balance the sample
    rows = [{
        "n": len(df),
        "median_cbe_pct": float(base.median()),
        "median_cation_meq": float(cat.median()),
        "median_anion_meq": float(an.median()),
        "anion_excess_median_meq": float((an - cat).median()),
        "hco3_share_of_anions_pct": float((hco3_meq / an * 100).median()),
        "cl_share_of_anions_pct": float((df["Cl"] / an * 100).median()),
        "median_cbe_if_HCO3_were_CaCO3": float(cbe_caco3.median()),
        "median_missing_cation_meq_to_balance": float(missing_cation_meq.median()),
        "interpretation": ("Persistent anion excess (~2.1 meq/L); HCO3 the largest anion. "
                           "Treating alkalinity as CaCO3 worsens the balance, so the "
                           "imbalance points to an unmeasured cation pool (no trace-metal "
                           "cations reported for this artisanal-mining catchment) or "
                           "cation under-recovery. Talensi is screening-only regardless."),
    }]
    pd.DataFrame(rows).to_csv(RES / "m6_talensi_cbe_diagnosis.csv", index=False)


# ---- D: discrimination -------------------------------------------------------
def discrimination():
    print("D. MRS discrimination ...")
    fp = pd.read_csv(RES / "m6_ng_field_pairs.csv")
    syn = pd.read_csv(RES / "m6_synthetic_class_distribution.csv")
    rows = [{
        "source": "Northern Ghana field (MRS)",
        "mean": float(fp.mrs.mean()), "std": float(fp.mrs.std()),
        "min": float(fp.mrs.min()), "max": float(fp.mrs.max()),
        "iqr": float(fp.mrs.quantile(0.75) - fp.mrs.quantile(0.25)),
    }]
    pd.DataFrame(rows).to_csv(RES / "m6_mrs_discrimination.csv", index=False)
    # synthetic class spread (proves discrimination): share of each true class
    tot = syn.groupby("true_class")["n"].sum()
    (tot / tot.sum()).rename("frac").reset_index().to_csv(
        RES / "m6_synthetic_class_spread.csv", index=False)


def main():
    data = m6.load_all()
    ng = data["northern_ghana"]
    clf = m6.TransferClassifier()
    cmap = m6.get_class_map()
    rng = np.random.default_rng(m6.SEED)
    gate_and_structural(ng, clf, cmap, rng)
    transport_sensitivity(ng, clf, cmap, rng)
    talensi_diagnosis()
    discrimination()
    print("Reviewer-response analyses complete ->", RES)


if __name__ == "__main__":
    main()
