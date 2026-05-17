import pandas as pd
import numpy as np
from pathlib import Path
from hydrosheaf.inference.network_fit import infer_edges, fit_network
from hydrosheaf.config import default_config

def run():
    # 1. Setup
    repo_root = Path(".").resolve()
    result_dir = repo_root / "M2" / "m2_benchmark" / "results"
    result_dir.mkdir(parents=True, exist_ok=True)
    
    # 2. Manu (Lower Anayari)
    print("Running Manu (Lower Anayari) Field Benchmark...")
    manu_df = pd.read_csv("manu.csv")
    manu_df = manu_df.rename(columns={
        "Sample ID": "site_id",
        "X coordinate": "lon",
        "Y coordinate": "lat",
        "Elevation": "elevation",
        "Temp": "temp_c"
    })
    
    config_m = default_config()
    config_m.geologic_bias = "crystalline"
    config_m.missing_policy = "impute_zero" # Critical fix
    config_m.active_minerals = ["calcite", "dolomite", "albite", "halite", "pyrite_oxidation_aerobic"]
    
    # Simple spatial inference for Manu
    samples_m = manu_df.to_dict('records')
    edges_m = infer_edges(samples_m, max_neighbors=3, method="simple")
    print(f"  Found {len(edges_m)} edges for Manu.")
    
    res_m = fit_network(samples_m, edges_m, config=config_m)
    print(f"  Successfully fitted {len(res_m)} edges for Manu.")
    
    # 3. Talensi
    print("Running Talensi Mining Area Field Benchmark...")
    talensi_df = pd.read_csv("talensi.csv")
    talensi_df = talensi_df.rename(columns={
        "Code": "site_id",
        "Longitude": "lon",
        "Latitude": "lat",
        "Elevation": "elevation",
        "Temp": "temp_c"
    })
    
    config_t = default_config()
    config_t.geologic_bias = "crystalline"
    config_t.missing_policy = "impute_zero" # Critical fix
    config_t.active_minerals = ["calcite", "dolomite", "albite", "halite", "pyrite_oxidation_aerobic"]
    
    # Simple spatial inference for Talensi
    samples_t = talensi_df.to_dict('records')
    edges_t = infer_edges(samples_t, max_neighbors=3, method="simple")
    print(f"  Found {len(edges_t)} edges for Talensi.")
    
    res_t = fit_network(samples_t, edges_t, config=config_t)
    print(f"  Successfully fitted {len(res_t)} edges for Talensi.")
    
    # 4. Combine and Save
    print("Saving combined field discovery results...")
    all_rows = []
    
    for res in res_m:
        row = {"edge_id": f"Manu_{res.edge_id}", "u": f"Manu_{res.u}", "v": f"Manu_{res.v}", 
               "transport_model": res.transport_model, "chemistry_r2": res.chemistry_r2}
        for label, extent in zip(res.z_labels, res.z_extents):
            row[f"extent_{label}"] = extent
        all_rows.append(row)
        
    for res in res_t:
        row = {"edge_id": f"Talensi_{res.edge_id}", "u": f"Talensi_{res.u}", "v": f"Talensi_{res.v}", 
               "transport_model": res.transport_model, "chemistry_r2": res.chemistry_r2}
        for label, extent in zip(res.z_labels, res.z_extents):
            row[f"extent_{label}"] = extent
        all_rows.append(row)
        
    if not all_rows:
        print("WARNING: No edges were fitted. Check your coordinate renaming and elevation data.")
        
    pd.DataFrame(all_rows).to_csv(result_dir / "field_discovery_results.csv", index=False)
    print(f"Done. Wrote {len(all_rows)} edges to {result_dir / 'field_discovery_results.csv'}")

if __name__ == "__main__":
    run()
