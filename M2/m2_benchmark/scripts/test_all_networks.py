import sys
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from M2.m2_benchmark.scripts.run_e4b_public_chem_validation import build_public_samples, ensure_source_files, candidate_edges, read_table
from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.config import default_config

def test_networks():
    data_dir = ensure_source_files()
    sites = read_table(data_dir, "Table_1_site_list_v4.txt")
    networks = sites["Network_name"].unique()
    print(f"Total networks: {len(networks)}")
    
    for nw in networks:
        try:
            samples, _ = build_public_samples(data_dir, nw, 30)
            if len(samples) < 10:
                continue
            edges = candidate_edges(samples, 5)
            config = default_config()
            config.ion_order = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "F"]
            config.weights = [1.0] * 8
            config.active_minerals = ["calcite", "dolomite", "gypsum", "halite"]
            config.latent_endmembers_enabled = True
            config.latent_endmembers_count = 2
            results, _ = fit_network_pipeline(samples.to_dict("records"), edges, config)
            print(f"{nw}: {len(results)} edges (samples: {len(samples)})")
        except Exception as e:
            # print(f"{nw}: Failed - {e}")
            pass

if __name__ == "__main__":
    test_networks()
