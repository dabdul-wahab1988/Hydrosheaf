"""M6 synthetic-field validation with KNOWN ground truth.

Forward-generates synthetic groundwater wells from a dilute recharge endmember by
applying KNOWN mineral/exchange reactions with KNOWN extents (plus SiO2/Sr
signatures, evaporation and analytical noise), then runs the SAME M6 inference at
each evidence tier and measures GROUND-TRUTH reaction recovery.

Purpose (answers the Q1 review):
  * Non-circular identifiability: recovery is scored against known truth and
    against gate-independent structural diagnostics, NOT the evidence-lift rule.
  * Validation: shows the tier-collapse corresponds to a REAL loss of recovering
    the true reactions, and that removing SiO2/Sr genuinely degrades silicate/
    carbonate recovery.
  * Discrimination: the ground-truth identifiability class spans the full range
    (identifiable -> non-identifiable), demonstrating the framework discriminates.
Deterministic (seed 1234).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synthetic_reaction_truth_model as rx

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
SEED = 1234

# dilute recharge endmember (mmol/L)
X0 = {"Ca": 0.30, "Mg": 0.15, "Na": 0.35, "K": 0.05, "HCO3": 1.0, "Cl": 0.20,
      "SO4": 0.08, "NO3": 0.05, "F": 0.005, "Fe": 0.001, "PO4": 0.002,
      "SiO2": 0.15, "Sr": 0.002}

ARCHETYPES = {
    "silicate": ["albite", "k_feldspar", "biotite", "anorthite"],
    "carbonate": ["calcite", "dolomite"],
    "evaporite": ["halite", "gypsum"],
    "ion_exchange": ["NaCa_exch", "CaNa_exch"],
    "nitrate": ["NO3src"],
    "mixed": ["albite", "calcite", "NaCa_exch", "gypsum", "NO3src", "dolomite"],
}
ARCH_W = {"silicate": 0.34, "carbonate": 0.16, "evaporite": 0.12,
          "ion_exchange": 0.12, "nitrate": 0.10, "mixed": 0.16}
# archetype -> the reaction family a correct dominant-process call should return
ARCH_TO_FAMILY = {"silicate": "silicate", "carbonate": "carbonate",
                  "evaporite": "evaporite", "ion_exchange": "ion_exchange",
                  "nitrate": "anthropogenic"}

# panels ordered by information content
PANELS = [("T0_majors", rx.PANELS["majors"]),
          ("T2_plus_F", rx.PANELS["plus_F"]),
          ("T3_plus_Sr_SiO2", rx.PANELS["plus_Sr_SiO2"])]

REQ = {}  # required major species per reaction (SiO2/Sr are diagnostic, never required)
for r in rx.REACTIONS:
    REQ[r.label] = {s for s, c in r.stoich.items()
                    if s not in ("SiO2", "Sr") and abs(c) >= 0.5}


def penalty_for(panel):
    ionset = set(panel)
    return np.array([1.0 if REQ[r.label] <= ionset else 40.0 for r in rx.REACTIONS])


def make_well(rng):
    arch = rng.choice(list(ARCH_W), p=list(ARCH_W.values()))
    pool = ARCHETYPES[arch]
    k = min(len(pool), rng.integers(1, 4))
    active = list(rng.choice(pool, size=k, replace=False))
    z = {lab: float(rng.uniform(0.25, 1.5)) for lab in active}
    x = {s: X0[s] for s in rx.SPECIES}
    for lab, ext in z.items():
        for s, c in rx.REACTIONS[rx.LABELS.index(lab)].stoich.items():
            x[s] = x.get(s, 0.0) + ext * c
    f = float(rng.uniform(1.0, 2.0))  # evaporative concentration
    x = {s: max(v * f, 0.0) for s, v in x.items()}
    x = {s: v * (1 + rng.normal(0, 0.03)) for s, v in x.items()}  # analytical noise
    return arch, set(active), x


def main():
    rng = np.random.default_rng(SEED)
    N = 360
    rows = []
    for w in range(N):
        arch, truth, x = make_well(rng)
        # Cl-estimated evaporation factor (as in the field pipeline)
        f_est = x["Cl"] / X0["Cl"] if X0["Cl"] > 0 else 1.0
        f_est = float(np.clip(f_est, 0.2, 6.0))
        for pname, panel in PANELS:
            cmap = rx.equivalence_classes(panel)
            resid = np.array([x[s] - f_est * X0[s] for s in panel])
            # map panel residual onto full species residual for fit()
            full = {s: (x[s] - f_est * X0[s]) if s in panel else 0.0 for s in rx.SPECIES}
            r_full = np.array([full[s] for s in panel])
            fit = rx.fit(r_full, panel, lambda_l1=0.01, lambda_l2=1e-3,
                         penalty_scales=penalty_for(panel))
            pred = rx.support(fit["extents"])
            p, rc, f1 = rx.prf1(truth, pred)
            cf1 = rx.class_f1(truth, pred, cmap)
            # process-family level recovery (what M6 actually reports)
            true_fams = {rx.FAMILY_OF[l] for l in truth}
            pred_fams = {rx.FAMILY_OF[l] for l in pred}
            fam_f1 = rx.prf1(true_fams, pred_fams)[2]
            # dominant-family correctness (by summed |extent|)
            fe = {}
            for l in pred:
                fe[rx.FAMILY_OF[l]] = fe.get(rx.FAMILY_OF[l], 0.0) + abs(
                    fit["extents"][rx.LABELS.index(l)])
            dom_pred = max(fe, key=fe.get) if fe else "none"
            tclass = rx.true_resolution_class(truth, pred, cmap)
            d = rx.diagnostics(panel)
            # extent error on truly-active reactions
            ext_err = np.mean([abs(fit["extents"][rx.LABELS.index(l)] - 1.0) * 0
                               for l in truth]) if truth else 0.0
            rows.append({
                "well": w, "archetype": arch, "panel": pname,
                "n_true": len(truth), "phase_precision": p, "phase_recall": rc,
                "phase_f1": f1, "class_f1": cf1, "family_f1": fam_f1,
                "dominant_correct": int(dom_pred in true_fams if arch == "mixed"
                                        else dom_pred == ARCH_TO_FAMILY[arch]),
                "true_class": tclass,
                "rank": d["rank"], "nullity": d["nullity"],
                "mutual_coherence": d["mutual_coherence"],
                "silicate_coherence": rx.silicate_identifiability(panel),
                "n_equivalence_classes": d["n_equivalence_classes"],
            })
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m6_synthetic_validation.csv", index=False)

    # aggregate: recovery F1 by archetype x panel (the non-circular headline)
    agg = (df.groupby(["archetype", "panel"])
           .agg(mean_phase_f1=("phase_f1", "mean"),
                mean_class_f1=("class_f1", "mean"),
                mean_family_f1=("family_f1", "mean"),
                dominant_accuracy=("dominant_correct", "mean"),
                frac_non_identifiable=("true_class",
                    lambda s: float((s == "non_identifiable").mean())),
                n=("well", "count")).reset_index())
    agg.to_csv(RESULTS / "m6_synthetic_recovery_by_tier.csv", index=False)

    # structural identifiability by panel (gate-independent)
    struct = (df.groupby("panel")
              .agg(rank=("rank", "first"), nullity=("nullity", "first"),
                   silicate_coherence=("silicate_coherence", "first"),
                   n_equivalence_classes=("n_equivalence_classes", "first"),
                   mean_phase_f1=("phase_f1", "mean"),
                   frac_non_identifiable=("true_class",
                       lambda s: float((s == "non_identifiable").mean()))).reset_index())
    struct.to_csv(RESULTS / "m6_synthetic_structural.csv", index=False)

    # true-class distribution by panel (discrimination evidence)
    dist = (df.groupby(["panel", "true_class"]).size().rename("n").reset_index())
    dist.to_csv(RESULTS / "m6_synthetic_class_distribution.csv", index=False)

    # headline print
    print("Ground-truth PROCESS-FAMILY F1 by archetype x tier (non-circular):")
    piv = agg.pivot(index="archetype", columns="panel", values="mean_family_f1")
    print(piv.round(3).to_string())
    print("\nDominant-family accuracy by archetype x tier:")
    print(agg.pivot(index="archetype", columns="panel",
                    values="dominant_accuracy").round(3).to_string())
    print("\nExact-phase F1 by archetype x tier:")
    print(agg.pivot(index="archetype", columns="panel",
                    values="mean_phase_f1").round(3).to_string())
    print("\nStructural identifiability by tier:")
    print(struct.round(3).to_string(index=False))
    print("\nTrue-class distribution by tier:")
    print(dist.pivot(index="panel", columns="true_class", values="n").to_string())


if __name__ == "__main__":
    main()
