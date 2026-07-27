"""Export auxiliary frames for M6 figures that the main analysis does not emit:
sample coordinates + aquifer (map), the evidence-tier ladder definition, and
Piper/Gibbs hydrochemical-context data for supplementary figures.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as m6

OUT = m6.RESULTS_DIR
OUT.mkdir(parents=True, exist_ok=True)


def export_coords(data):
    rows = []
    for name, df in data.items():
        d = df.copy()
        if name == "northern_ghana":
            d = d[d.season == "Dry"]
        for _, r in d.iterrows():
            rows.append({"dataset": name, "sample_id": r["sample_id"],
                         "latitude": r.get("Latitude"), "longitude": r.get("Longitude"),
                         "aquifer": r.get("Aquifer_Type")})
    pd.DataFrame(rows).to_csv(OUT / "m6_map_coordinates.csv", index=False)


def export_tier_ladder():
    rows = []
    caps = {"Tier 0": "majors", "Tier 1": "+ isotopes", "Tier 2": "+ fluoride",
            "Tier 3": "+ Sr / SiO2", "Tier 4": "+ saturation indices"}
    attain = {"northern_ghana": 4, "manu": 2, "talensi": 1}
    for i, (tier, desc) in enumerate(caps.items()):
        for ds, cap in attain.items():
            rows.append({"tier": tier, "tier_index": i, "description": desc,
                         "dataset": ds, "attained": int(i <= cap)})
    pd.DataFrame(rows).to_csv(OUT / "m6_tier_ladder.csv", index=False)


def export_hydrochem_context(data):
    """Piper (cation/anion %) and Gibbs (TDS vs Na/(Na+Ca)) context, all datasets."""
    rows = []
    for name, df in data.items():
        for _, r in df.iterrows():
            ca, mg, na, k = [r.get(i, np.nan) for i in ["Ca", "Mg", "Na", "K"]]
            hco3, cl, so4 = [r.get(i, np.nan) for i in ["HCO3", "Cl", "SO4"]]
            catmeq = np.array([ca * 2, mg * 2, na, k])
            anmeq = np.array([hco3, cl, so4 * 2])
            cs, as_ = np.nansum(catmeq), np.nansum(anmeq)
            if cs <= 0 or as_ <= 0:
                continue
            tds = np.nansum([r.get(i, 0) * m6.m5.MOLAR_MASS_G_MOL[i]
                             for i in ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3"]])
            rows.append({
                "dataset": name, "sample_id": r["sample_id"],
                "aquifer": r.get("Aquifer_Type"),
                "ca_pct": 100 * ca * 2 / cs, "mg_pct": 100 * mg * 2 / cs,
                "nak_pct": 100 * (na + k) / cs,
                "hco3_pct": 100 * hco3 / as_, "cl_pct": 100 * cl / as_,
                "so4_pct": 100 * so4 * 2 / as_,
                "tds_mgL": tds,
                "na_ratio": na / (na + ca) if (na + ca) > 0 else np.nan,
                "cl_ratio": cl / (cl + hco3) if (cl + hco3) > 0 else np.nan,
            })
    pd.DataFrame(rows).to_csv(OUT / "m6_hydrochem_context.csv", index=False)


def main():
    data = m6.load_all()
    export_coords(data)
    export_tier_ladder()
    export_hydrochem_context(data)
    print("Auxiliary figure data exported to", OUT)


if __name__ == "__main__":
    main()
