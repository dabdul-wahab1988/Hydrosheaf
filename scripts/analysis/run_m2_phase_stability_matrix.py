"""M2 per-mineral Phase Stability Matrix (PSI) with GEOLOGICALLY-CORRECT dictionaries.

ACCURACY FIX. Figure 7 previously read M6/m6_phase_stability_index.csv, which applied a
single generic mineral list (INCLUDING gypsum) to both Ghana sites. For Talensi
(crystalline Birimian basement, artisanal gold mining) gypsum is geologically impossible
- there are no evaporite deposits - so the M6 run produced a false-positive gypsum
PSI=0.99, contradicting the pyrite-oxidation mining signature and the paper text.

This script computes the per-mineral inclusion probability (PSI) directly from the M2
field pipeline using SITE-SPECIFIC, lithology-correct dictionaries:
  - Lower Anayari (Manu): calcite, dolomite, gypsum, albite, halite, fluorite
  - Talensi:              calcite, dolomite, pyrite_oxidation_aerobic, albite, halite   (NO gypsum)
Exchange and nitrogen-cycling terms are enabled for both. PSI per mineral = mean
Monte-Carlo inclusion probability across that site's edges (matching run_edge_psi.py's
MC settings). Output schema matches the file Figure 7 consumes: site, mineral, psi, mean, std.

Writes: M2/m2_benchmark/results/m2_phase_stability_index.csv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf import Config
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.uncertainty.sensitivity import analyze_sensitivity_mc

IONS = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe']

SITES = {
    "Manu": dict(
        csv="LowerAnayari/manu.csv", id_col="Sample ID",
        minerals=["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"],
    ),
    "Talensi": dict(
        csv="Talensi_MiningArea/talensi.csv", id_col="Code",
        # NO gypsum: crystalline basement mining area has no evaporites; SO4 -> pyrite.
        minerals=["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"],
    ),
}


def parse_val(v):
    if isinstance(v, str):
        v = v.strip()
        if v in ("", "-", "nd", "ND", "<") or v.startswith("<"):
            return 0.0005
    try:
        return float(v)
    except Exception:
        return 0.0


def load_samples(csv, id_col):
    candidates = [ROOT / "data" / "FieldData" / csv, ROOT / "data" / csv]
    path = next((c for c in candidates if c.is_file()), candidates[0])
    df = pd.read_csv(path)
    samples = {}
    for _, row in df.iterrows():
        sid = str(row.get(id_col))
        samples[sid] = {ion: mgL_to_mmolL(parse_val(row.get(ion, 0.0)), ion) for ion in IONS}
    return samples


def run_site(site, cfg, n_trials=30):
    samples = load_samples(cfg["csv"], cfg["id_col"])
    config = Config(
        ion_order=IONS, weights=[1.0] * 10, conservative_weights=[0.01] * 10,
        isotope_enabled=False, active_minerals=cfg["minerals"],
        exchange_enabled=True, honest_modeling=True, geologic_bias="crystalline",
    )
    config.sensitivity_analysis_enabled = True
    field = pd.read_csv(ROOT / "M2" / "m2_benchmark" / "results" / "field_discovery_results.csv")
    edges = field[field["edge_id"].str.startswith(site)]
    incl, ext = {}, {}
    n = 0
    for _, row in edges.iterrows():
        u, v = str(row["u"]), str(row["v"])
        if u not in samples or v not in samples:
            continue
        base = {"x_u": [samples[u][i] for i in IONS], "x_v": [samples[v][i] for i in IONS],
                "config": config, "obs_v": samples[v]}
        rep = analyze_sensitivity_mc(fit_edge, base, config, n_trials=n_trials,
                                     concentration_error_pct=0.05)
        for k, p in rep.inclusion_probabilities.items():
            incl.setdefault(k, []).append(float(p))
        for k, m in rep.extent_means.items():
            ext.setdefault(k, []).append(abs(float(m)))
        n += 1
    print(f"  {site}: {n} edges processed")
    rows = []
    for mineral, probs in incl.items():
        rows.append({"site": site, "mineral": mineral,
                     "psi": float(np.mean(probs)),
                     "mean": float(np.mean(ext.get(mineral, [0.0]))),
                     "std": float(np.std(probs))})
    return rows


def main():
    # analyze_sensitivity_mc draws from the global numpy RNG, so the PSI matrix
    # (and Figure 7) drifts between runs unless the seed is fixed. Same seed as
    # scripts/analysis/run_edge_psi.py.
    np.random.seed(20260820)
    all_rows = []
    for site, cfg in SITES.items():
        print(f"Processing {site} (minerals: {cfg['minerals']})...")
        all_rows.extend(run_site(site, cfg))
    out = pd.DataFrame(all_rows)
    save = ROOT / "M2" / "m2_benchmark" / "results" / "m2_phase_stability_index.csv"
    out.to_csv(save, index=False)
    print(f"\nSaved geologically-correct PSI matrix -> {save}")
    for site in SITES:
        top = out[out.site == site].sort_values("psi", ascending=False).head(6)
        print(f"\n{site} top PSI:\n{top[['mineral','psi']].to_string(index=False)}")


if __name__ == "__main__":
    main()
