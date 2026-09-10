"""Revision analysis 3 (CAGEO-D-26-00847): PSI separation for degenerate reaction pairs.

Answers Reviewer 1 Major 15: does the process-stability index separate the
near-degenerate pairs (calcite/dolomite; CaNa_exch/NaCa_exch) or assign
similar stability to both members?

For a sample of field edges (20 per site) we record the per-reaction Monte
Carlo inclusion probability for the four reactions of interest and report:
  - co-activation rate (both members activated in the same trial),
  - PSI per member and the PSI difference,
  - which member dominates.

Output: M2/m2_benchmark/results/revision/psi_pair_separation.csv
"""
from __future__ import annotations

import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf import Config  # noqa: E402
from hydrosheaf.data.units import mgL_to_mmolL  # noqa: E402
from hydrosheaf.uncertainty.sensitivity import analyze_sensitivity_mc  # noqa: E402
from hydrosheaf.inference.edge_fit import fit_edge  # noqa: E402

OUT_DIR = PROJECT_ROOT / "M2" / "m2_benchmark" / "results" / "revision"

PAIRS = [("calcite", "dolomite"), ("CaNa_exch", "NaCa_exch")]
IONS = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"]


def parse_val(val):
    if isinstance(val, str) and "<" in val:
        return 0.0005
    try:
        return float(val)
    except (TypeError, ValueError):
        return 0.0


def load_samples(csv_file: Path) -> dict[str, dict[str, float]]:
    df = pd.read_csv(csv_file)
    samples = {}
    for _, row in df.iterrows():
        sid = str(row.get("Sample ID", row.get("Code")))
        samples[sid] = {ion: mgL_to_mmolL(parse_val(row.get(ion, 0.0)), ion) for ion in IONS}
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit PSI separation for selected degenerate reaction pairs.")
    parser.add_argument(
        "--max-edges-per-site",
        type=int,
        default=None,
        help="Optional cap per site. The default audits every retained canonical field edge.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260820,
        help="Base seed for the NumPy perturbations used by analyze_sensitivity_mc.",
    )
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    field = pd.read_csv(PROJECT_ROOT / "M2" / "m2_benchmark" / "results" / "field_discovery_results.csv")
    sites = [
        ("Manu", PROJECT_ROOT / "data" / "FieldData" / "LowerAnayari" / "manu.csv",
         ["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"]),
        ("Talensi", PROJECT_ROOT / "data" / "FieldData" / "Talensi_MiningArea" / "talensi.csv",
         ["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"]),
    ]
    rows = []
    for site_index, (site_name, csv_file, active_minerals) in enumerate(sites):
        samples = load_samples(csv_file)
        config = Config(
            ion_order=IONS,
            weights=[1.0] * 10,
            conservative_weights=[0.01] * 10,
            isotope_enabled=False,
            active_minerals=active_minerals,
            exchange_enabled=True,
            honest_modeling=True,
            geologic_bias="crystalline",
        )
        site_edges = field[field["edge_id"].str.startswith(site_name)]
        if args.max_edges_per_site is not None:
            site_edges = site_edges.head(args.max_edges_per_site)
        for edge_index, (_, row) in enumerate(site_edges.iterrows()):
            u, v = str(row["u"]), str(row["v"])
            if u not in samples or v not in samples:
                continue
            x_u = [samples[u][ion] for ion in IONS]
            x_v = [samples[v][ion] for ion in IONS]
            config.sensitivity_analysis_enabled = True
            # analyze_sensitivity_mc currently consumes NumPy's module-level
            # RNG. Seed each edge with a stable offset so this diagnostic is
            # reproducible independently of previous script/process state.
            np.random.seed(args.seed + site_index * 100_000 + edge_index)
            report = analyze_sensitivity_mc(
                fit_edge,
                {"x_u": x_u, "x_v": x_v, "config": config, "obs_v": samples[v]},
                config,
                n_trials=30,
                concentration_error_pct=0.05,
            )
            inc = report.inclusion_probabilities
            means = report.extent_means
            for a, b in PAIRS:
                pa = float(inc.get(a, 0.0))
                pb = float(inc.get(b, 0.0))
                ma = float(means.get(a, 0.0))
                mb = float(means.get(b, 0.0))
                both = (np.sign(ma) != 0) and (np.sign(mb) != 0) and (np.sign(ma) == np.sign(mb))
                rows.append(
                    {
                        "site": site_name,
                        "edge_id": row["edge_id"],
                        "pair": f"{a}~{b}",
                        "psi_a": pa,
                        "psi_b": pb,
                        "psi_diff_abs": abs(pa - pb),
                        "co_activated_both_positive": bool(both),
                        "mean_extent_a": ma,
                        "mean_extent_b": mb,
                        "n_trials": 30,
                        "seed": args.seed + site_index * 100_000 + edge_index,
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "psi_pair_separation.csv", index=False)
    print(out.groupby(["site", "pair"])[["psi_a", "psi_b", "psi_diff_abs"]].mean().round(3).to_string())
    print(f"Wrote {OUT_DIR / 'psi_pair_separation.csv'}")


if __name__ == "__main__":
    main()
