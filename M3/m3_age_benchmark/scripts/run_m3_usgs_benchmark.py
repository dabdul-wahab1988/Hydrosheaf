"""
USGS Groundwater Age Benchmark for M3.
Evaluates Hydrosheaf's multi-tracer performance on National and Regional datasets.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import math

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.nuclear.joint_lpm import fit_lpm_model, SUPPORTED_LPM_MODELS
from hydrosheaf.log import get_logger

logger = get_logger("m3.usgs_benchmark")

# Benchmark Paths
BENCHMARK_DIR = Path(__file__).resolve().parents[1]
# Real data paths from M2 external
M2_USGS_DATA = REPO_ROOT / "M2" / "m2_benchmark" / "external" / "usgs_age" / "input" / "DataForNationalGroundwaterAge_1_1"
RESULT_DIR = BENCHMARK_DIR / "results"

import networkx as nx
from hydrosheaf.nuclear.network_aging import infer_network_ages_bayesian
from hydrosheaf.graph.build import infer_edges_from_coordinates
from hydrosheaf.nuclear.nuclides import TRITIUM

def load_usgs_national_dataset():
    """
    Loads and joins real USGS National Public-Supply Aquifer dataset.
    """
    tracers_path = M2_USGS_DATA / "Table_3_Tracers.txt"
    ages_path = M2_USGS_DATA / "Table_2_Ages.txt"
    sites_path = M2_USGS_DATA / "Table_1_Sites.txt"
    
    if not tracers_path.exists() or not ages_path.exists() or not sites_path.exists():
        logger.warning("Real USGS National dataset not found. Generating synthetic benchmark.")
        return generate_synthetic_usgs_proxy()
    
    print(f"Loading real USGS data from {M2_USGS_DATA}...")
    
    # Load tracers (tab-separated)
    df_tracers = pd.read_csv(tracers_path, sep="\t", na_values="na")
    # Load ages
    df_ages = pd.read_csv(ages_path, sep="\t", na_values="na")
    # Load sites
    df_sites = pd.read_csv(sites_path, sep="\t", na_values="na")
    
    # Join tracers and ages on SampleID
    df = pd.merge(df_tracers, df_ages[["SampleID", "Rpt_TotAge_yrs", "LPM_Name"]], on="SampleID", how="inner")
    # Join sites for coordinates and depth
    df = pd.merge(df, df_sites[["SampleID", "LatDD83", "LongDD83", "Depth_m", "StudyUnit"]], on="SampleID", how="inner")
    
    # Map column names
    mapping = {
        "SampleID": "site_id",
        "3H_TU": "tritium_TU",
        "3He_trit_TU": "he3_trit_TU",
        "SF6_pptv": "sf6_pptv",
        "CFC-11_pptv": "cfc11_pptv",
        "CFC-12_pptv": "cfc12_pptv",
        "CFC-113_pptv": "cfc113_pptv",
        "14C_pmC": "c14_pmc",
        "4He_ccpg": "he4_ccpg",
        "Rpt_TotAge_yrs": "reference_age_years",
        "LatDD83": "lat",
        "LongDD83": "lon",
        "Depth_m": "depth_m"
    }
    df.rename(columns=mapping, inplace=True)
    
    if "StudyUnit_x" in df.columns:
        df.rename(columns={"StudyUnit_x": "StudyUnit"}, inplace=True)
    elif "StudyUnit_y" in df.columns:
        df.rename(columns={"StudyUnit_y": "StudyUnit"}, inplace=True)
    elif "StudyUnit" not in df.columns:
        # Fallback: check for case-insensitive match
        for col in df.columns:
            if col.lower() == "studyunit":
                df.rename(columns={col: "StudyUnit"}, inplace=True)
                break

    # Parse SampleDate
    df["sample_year"] = pd.to_datetime(df["SampleDate"], errors="coerce").dt.year
    df.dropna(subset=["sample_year"], inplace=True)
    
    return df

def run_graph_benchmark(df):
    """
    Performs graph-regularized age inference on clusters of sites.
    """
    logger.info("Starting Graph-Regularized M3 benchmark...")
    graph_results = []
    
    # Group by StudyUnit for regional clusters
    units = df["StudyUnit"].unique()
    for unit in units[:3]: # Process first 3 units for benchmark demo
        unit_df = df[df["StudyUnit"] == unit].copy()
        if len(unit_df) < 5: continue
        
        print(f"\nProcessing Study Unit: {unit} ({len(unit_df)} wells)")
        
        # 1. Build Graph
        # Use simple coordinate-based inference
        sample_dicts = unit_df.to_dict("records")
        # Rename lat/lon to easting/northing for build_edges (crude approx)
        for s in sample_dicts:
            s["easting"] = s["lon"]
            s["northing"] = s["lat"]
            s["screen_depth"] = s.get("depth_m", 10.0)
            
        edges = infer_edges_from_coordinates(sample_dicts, max_neighbors=2)
        graph = nx.DiGraph()
        for edge in edges:
            graph.add_edge(edge.u, edge.v, length_m=edge.attrs.get("distance_m", 1000.0))
            
        if graph.number_of_edges() == 0:
            print(f"  Skipping {unit}: No edges inferred.")
            continue
            
        # 2. Run Bayesian Network Aging
        obs_map = {row["site_id"]: row["tritium_TU"] for _, row in unit_df.iterrows() if not math.isnan(row["tritium_TU"])}
        sig_map = {s: max(0.5, v*0.1) for s, v in obs_map.items()}
        
        try:
            # We'll use 3H as the primary nuclide for this demo
            net_ages = infer_network_ages_bayesian(
                graph, 
                obs_map, 
                sig_map, 
                sample_date=2015.0, 
                nuclide=TRITIUM,
                n_samples=500, # Faster sampling for benchmark
                n_chains=2
            )
            
            for site_id, res in net_ages.items():
                if site_id.startswith("_"): continue
                ref = unit_df[unit_df["site_id"] == site_id]["reference_age_years"].iloc[0]
                ref_val = _parse_age(ref)
                
                graph_results.append({
                    "study_unit": unit,
                    "site_id": site_id,
                    "ref_age": ref_val,
                    "est_age_graph": res["mean_age_years"],
                    "error_graph": abs(res["mean_age_years"] - ref_val) if not math.isnan(ref_val) else np.nan
                })
        except Exception as e:
            print(f"  Error in Bayesian inference for {unit}: {e}")
            
    if graph_results:
        res_df = pd.DataFrame(graph_results)
        res_df.to_csv(RESULT_DIR / "m3_graph_regularized_results.csv", index=False)
        print(f"\nGraph Regularization Complete. Saved to {RESULT_DIR}")
        print(f"Avg Error (Graph-Regularized): {res_df['error_graph'].mean():.2f} yrs")

def generate_synthetic_usgs_proxy():
    """Generates a synthetic proxy for the USGS dataset to test the M3 pipeline."""
    np.random.seed(42)
    n_samples = 50
    data = {
        "site_id": [f"USGS-{i:04d}" for i in range(n_samples)],
        "sample_year": np.random.choice([2010, 2012, 2015, 2017], n_samples),
        "tritium_TU": np.random.gamma(2, 2, n_samples),
        "sf6_pptv": np.random.gamma(3, 1, n_samples),
        "c14_pmc": np.random.uniform(50, 100, n_samples),
        "kr85_pptv": np.random.gamma(2, 10, n_samples),
        "ar39_pmc": np.random.uniform(70, 95, n_samples),
        "reference_age_years": np.random.lognormal(log(50), 1, n_samples)
    }
    return pd.DataFrame(data)

def _parse_age(val):
    if pd.isna(val): return np.nan
    s = str(val).strip()
    if not s: return np.nan
    # Take first part before space or parenthesis
    part = s.split()[0].split('(')[0]
    try:
        return float(part)
    except ValueError:
        return np.nan

def run_benchmark(df):
    results = []
    logger.info(f"Starting M3 benchmark over {len(df)} samples...")
    
    # Take a subset if too large for quick profiling, or run all for final
    df_run = df.head(100) # Process first 100 for now
    
    for i, row in df_run.iterrows():
        if i % 10 == 0:
            print(f"Processing sample {i}/{len(df_run)}...")
        obs = {
            "tritium_TU": row.get("tritium_TU"),
            "he3_trit_TU": row.get("he3_trit_TU"),
            "sf6_pptv": row.get("sf6_pptv"),
            "cfc11_pptv": row.get("cfc11_pptv"),
            "cfc12_pptv": row.get("cfc12_pptv"),
            "cfc113_pptv": row.get("cfc113_pptv"),
            "c14_pmc": row.get("c14_pmc"),
            "kr85_pptv": row.get("kr85_pptv"),
            "ar39_pmc": row.get("ar39_pmc"),
            "he4_ccpg": row.get("he4_ccpg")
        }
        sample_year = float(row["sample_year"])
        ref_age = _parse_age(row.get("reference_age_years"))
        
        # Test 1: Single-Tracer Baseline (Tritium)
        fit_3h = fit_lpm_model({"tritium_TU": obs["tritium_TU"]}, sample_year=sample_year, model="EM")
        
        # Test 2: Multi-Tracer Hydrosheaf
        # We use GA (Gamma) or rank models by AIC
        fit_multi = fit_lpm_model(obs, sample_year=sample_year, model="GA")
        
        results.append({
            "site_id": row["site_id"],
            "ref_age": ref_age,
            "est_age_3h": fit_3h.parameters.get("mean_age_years", np.nan),
            "est_age_multi": fit_multi.parameters.get("mean_age_years", np.nan),
            "multi_model": fit_multi.model,
            "n_tracers": fit_multi.n_tracers,
            "error_3h": abs(fit_3h.parameters.get("mean_age_years", 0) - ref_age) if not math.isnan(ref_age) else np.nan,
            "error_multi": abs(fit_multi.parameters.get("mean_age_years", 0) - ref_age) if not math.isnan(ref_age) else np.nan
        })
    
    res_df = pd.DataFrame(results)
    res_df.to_csv(RESULT_DIR / "m3_usgs_benchmark_results.csv", index=False)
    
    # Summary Stats
    logger.info("M3 Benchmark Complete.")
    logger.info(f"Avg Error (Tritium-only): {res_df['error_3h'].mean():.2f} yrs")
    logger.info(f"Avg Error (Multi-Tracer):  {res_df['error_multi'].mean():.2f} yrs")
    
def log(x): return math.log(x)

if __name__ == "__main__":
    df = load_usgs_national_dataset()
    
    # 1. Run Point-wise Multi-Tracer Benchmark
    run_benchmark(df)
    
    # 2. Run Graph-Regularized Benchmark
    run_graph_benchmark(df)
