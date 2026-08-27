import pandas as pd
import numpy as np
import sys
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf import Config
from hydrosheaf.data.units import mgL_to_mmolL
from hydrosheaf.uncertainty.sensitivity import analyze_sensitivity_mc
from hydrosheaf.inference.edge_fit import fit_edge

def parse_val(val):
    if isinstance(val, str) and '<' in val:
        return 0.0005
    try:
        return float(val)
    except:
        return 0.0

def resolve_field_csv(csv_file: str) -> Path:
    """Resolve the field sample CSV under the canonical data locations.

    Revision fix (CAGEO-D-26-00847): the submitted pipeline expected
    ``data/<site>/<file>`` while the canonical field inputs live under
    ``data/FieldData/<site>/<file>`` (see run_m2_field_benchmarks.py).
    Both locations are tried so the script works regardless of layout.
    """
    candidates = [
        ROOT / "data" / "FieldData" / csv_file,
        ROOT / "data" / csv_file,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"field sample CSV not found; tried: {[str(c) for c in candidates]}"
    )


def process_site(site_name, csv_file, active_minerals, n_trials=30, max_edges=None):
    df_samples = pd.read_csv(resolve_field_csv(csv_file))
    samples = {}
    ions = ['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4', 'NO3', 'F', 'Fe']

    for _, row in df_samples.iterrows():
        sid = str(row.get('Sample ID', row.get('Code')))
        s = {}
        for ion in ions:
            raw_val = row.get(ion, 0.0)
            s[ion] = mgL_to_mmolL(parse_val(raw_val), ion)
        samples[sid] = s

    config = Config(
        ion_order=ions,
        weights=[1.0]*10,
        conservative_weights=[0.01]*10,
        isotope_enabled=False,
        active_minerals=active_minerals,
        exchange_enabled=True,
        honest_modeling=True,
        geologic_bias="crystalline"
    )

    field = pd.read_csv(ROOT / "M2" / "m2_benchmark" / "results" / "field_discovery_results.csv")
    # ARCHITECTURAL SHIFT: Process ALL edges for the full network PSI visualization
    site_edges = field[field["edge_id"].str.startswith(site_name)].copy()

    results = []
    print(f"Calculating PSI for all {len(site_edges)} edges in {site_name}...")
    if max_edges is not None:
        site_edges = site_edges.head(max_edges)
        print(f"  limited to first {len(site_edges)} edges for this run.")
    for _, row in site_edges.iterrows():
        u, v = str(row["u"]), str(row["v"])
        if u in samples and v in samples:
            x_u = [samples[u][ion] for ion in ions]
            x_v = [samples[v][ion] for ion in ions]

            base_inputs = {
                "x_u": x_u,
                "x_v": x_v,
                "config": config,
                "obs_v": samples[v]
            }

            config.sensitivity_analysis_enabled = True
            report = analyze_sensitivity_mc(
                fit_edge,
                base_inputs,
                config,
                n_trials=n_trials,
                concentration_error_pct=0.05
            )

            # Find dominant family
            best_fam = "Conservative"
            best_psi = 0.5
            max_mean = 0.0
            dom_proc = None

            valid_means = {k: abs(v) for k, v in report.extent_means.items() if "exch" not in k and "input" not in k and "denit" not in k and "NO3src" not in k}
            if valid_means:
                dom_proc = max(valid_means, key=valid_means.get)
                max_mean = valid_means[dom_proc]
                best_psi = report.inclusion_probabilities.get(dom_proc, 0.0)

            def get_family(p):
                if not p: return "Conservative"
                p = p.lower()
                if any(x in p for x in ["albite", "anorthite"]): return "Plagioclase"
                if any(x in p for x in ["k_feldspar", "microcline"]): return "K-Feldspar"
                if any(x in p for x in ["biotite", "chlorite", "enstatite", "diopside"]): return "Ferromagnesian"
                if any(x in p for x in ["calcite", "dolomite", "magnesite", "aragonite"]): return "Carbonates"
                if any(x in p for x in ["gypsum", "halite", "sylvite", "fluorite", "apatite"]): return "Evaporites"
                if any(x in p for x in ["no3", "potash", "nitrification", "road_salt"]): return "Anthropogenic"
                if any(x in p for x in ["denit", "pyrite", "sulfate_reduction", "iron_reduction"]): return "Redox"
                return "Conservative"

            if max_mean > 0.01:
                best_fam = get_family(dom_proc)

            results.append({
                "edge_id": row["edge_id"],
                "u": u,
                "v": v,
                "family": best_fam,
                "psi": best_psi
            })

    return results

def main():
    parser = argparse.ArgumentParser(description="Generate M2 field edge process-stability probabilities.")
    parser.add_argument("--trials", type=int, default=30, help="Monte Carlo trials per edge.")
    parser.add_argument("--max-edges-per-site", type=int, default=None, help="Optional smoke-test edge limit per site.")
    parser.add_argument("--seed", type=int, default=20260820,
                        help="RNG seed. analyze_sensitivity_mc draws from the global "
                             "numpy RNG, so PSI values -- and hence the reported family "
                             "distribution -- drift between runs unless this is fixed.")
    args = parser.parse_args()
    np.random.seed(args.seed)

    manu_res = process_site(
        "Manu",
        "LowerAnayari/manu.csv",
        ["calcite", "dolomite", "gypsum", "albite", "halite", "fluorite"],
        n_trials=args.trials,
        max_edges=args.max_edges_per_site,
    )
    tal_res = process_site(
        "Talensi",
        "Talensi_MiningArea/talensi.csv",
        ["calcite", "dolomite", "pyrite_oxidation_aerobic", "albite", "halite"],
        n_trials=args.trials,
        max_edges=args.max_edges_per_site,
    )

    out_df = pd.DataFrame(manu_res + tal_res)
    save_path = ROOT / "M2" / "m2_benchmark" / "results" / "top_edges_psi.csv" # Keep name for compatibility
    out_df.to_csv(save_path, index=False)
    print(f"Saved full network PSI results to {save_path}")

if __name__ == "__main__":
    main()
