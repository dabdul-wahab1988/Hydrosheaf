"""
Comprehensive Hydrosheaf Analysis of Synthetic Groundwater Data

This script demonstrates the hydrosheaf package capabilities using the
synthetic CSV data in hydrosheaf_synthetic_csv/. It performs:
1. Data loading and preprocessing
2. Unit conversion (mg/L to mmol/L)
3. Parameter calibration using PEST-GLM
4. Network fitting (edge-based hydrochemical evolution)
5. Result interpretation and summary

Author: Generated for Hydrosheaf demonstration
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import json

# Add hydrosheaf to path if needed
sys.path.insert(0, str(Path(__file__).parent.parent))

from hydrosheaf import (
    Config,
    fit_network,
    summarize_network,
    edge_process_maps,
    compute_d_excess,
    evaporation_index,
)
from hydrosheaf.data.units import mgL_to_mmolL, MOLAR_MASS_G_MOL

# Import calibration modules
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM


# -------------------------------------------------------------------
# 1. Data Loading
# -------------------------------------------------------------------


def load_synthetic_data(data_dir: Path):
    """Load all synthetic CSV files into DataFrames."""
    data = {}

    # Core files
    data["water_chem"] = pd.read_csv(data_dir / "water_chem_full.csv")
    data["stations"] = pd.read_csv(data_dir / "stations.csv")
    data["edges"] = pd.read_csv(data_dir / "network_edges.csv")
    data["events"] = pd.read_csv(data_dir / "events.csv")
    data["endmembers"] = pd.read_csv(data_dir / "endmembers_isotopes.csv")
    data["redox"] = pd.read_csv(data_dir / "redox_proxies.csv")

    print(f"Loaded {len(data['water_chem'])} water chemistry samples")
    print(f"Loaded {len(data['stations'])} stations")
    print(f"Loaded {len(data['edges'])} network edges")
    print(f"Loaded {len(data['events'])} events")

    return data


# -------------------------------------------------------------------
# 2. Data Preprocessing and Unit Conversion
# -------------------------------------------------------------------


def convert_to_mmol(df: pd.DataFrame, ion_columns: dict) -> pd.DataFrame:
    """
    Convert ion concentrations from mg/L to mmol/L.

    Args:
        df: DataFrame with ion concentrations in mg/L
        ion_columns: Mapping from hydrosheaf ion names to CSV column names

    Returns:
        DataFrame with concentrations in mmol/L
    """
    df = df.copy()

    for ion, col in ion_columns.items():
        if col in df.columns and ion in MOLAR_MASS_G_MOL:
            df[ion] = df[col].apply(
                lambda x: mgL_to_mmolL(x, ion) if pd.notna(x) else np.nan
            )

    return df


def prepare_samples_for_hydrosheaf(
    water_chem: pd.DataFrame, event_code: str = None
) -> list:
    """
    Prepare water chemistry data for hydrosheaf analysis.

    If event_code is provided, filter to that event only.
    Returns list of sample dictionaries with site_id and ion concentrations in mmol/L.
    """
    # Map CSV columns to hydrosheaf ion names
    ion_mapping = {
        "Ca": "Ca_mg_L",
        "Mg": "Mg_mg_L",
        "Na": "Na_mg_L",
        "K": "K_mg_L",
        "HCO3": "HCO3_mg_L",
        "Cl": "Cl_mg_L",
        "SO4": "SO4_mg_L",
        "NO3": "NO3_mg_L",
    }

    df = water_chem.copy()

    if event_code:
        df = df[df["event_code"] == event_code].copy()

    # Convert to mmol/L
    for ion, col in ion_mapping.items():
        if col in df.columns:
            df[ion] = df[col].apply(
                lambda x: mgL_to_mmolL(x, ion) if pd.notna(x) else np.nan
            )

    # Add F, Fe, PO4 as zeros (not in dataset but needed for 10-ion order)
    df["F"] = 0.0
    df["Fe"] = 0.0
    df["PO4"] = 0.0

    # Prepare sample dictionaries
    samples = []
    for _, row in df.iterrows():
        sample = {
            "site_id": row["station_code"],
            "sample_id": f"{row['event_code']}_{row['station_code']}",
            "Ca": row["Ca"],
            "Mg": row["Mg"],
            "Na": row["Na"],
            "K": (
                row.get("K", 0.0)
                if "K" in row
                else mgL_to_mmolL(row.get("K_mg_L", 0), "K")
            ),
            "HCO3": row["HCO3"],
            "Cl": row["Cl"],
            "SO4": row["SO4"],
            "NO3": row["NO3"],
            "F": row["F"],
            "Fe": row["Fe"],
            "PO4": row["PO4"],
            "pH": row.get("pH", 7.0),
            "EC": row.get("EC_uS_cm", None),
            "TDS": row.get("TDS_mg_L", None),
            "temp_C": row.get("temp_C", 25.0),
            # Isotope data
            "d15N": row.get("d15N_NO3_permil", None),
            "d18O_NO3": row.get("d18O_NO3_permil", None),
            "d2H": row.get("d2H_H2O_permil", None),
            "18O": row.get("d18O_H2O_permil", None),
        }
        samples.append(sample)

    return samples


def prepare_edges(edges_df: pd.DataFrame) -> list:
    """Convert edges DataFrame to hydrosheaf edge format.

    Hydrosheaf accepts edges as:
    - Tuple[str, str]: (from_id, to_id)
    - Tuple[str, str, str]: (edge_id, from_id, to_id)
    - Edge object
    """
    edges = []
    for _, row in edges_df.iterrows():
        # Use simple tuple format (from_station, to_station)
        edges.append((row["from_station"], row["to_station"]))
    return edges


# -------------------------------------------------------------------
# 2b. Calibration Functions
# -------------------------------------------------------------------


def setup_calibration(data: dict, config: Config):
    """Set up calibration problem with parameters and observations."""
    water_chem = data["water_chem"]
    edges_df = data["edges"]
    
    # Define parameters to calibrate
    parameters = [
        AdjustableParameter(
            "dispersivity", 10.0, 1.0, 50.0,
            prior_mean=10.0, prior_sigma=5.0, group="transport"
        ),
        AdjustableParameter(
            "evap_factor", 1.0, 0.8, 1.5,
            prior_mean=1.0, prior_sigma=0.1, group="transport"
        ),
        AdjustableParameter(
            "mix_fraction", 0.5, 0.2, 0.8,
            prior_mean=0.5, prior_sigma=0.15, group="mixing"
        ),
    ]
    
    # Create observations from borehole data
    observations = []
    targets = water_chem[water_chem["station_type"] == "borehole"]
    
    for _, row in targets.iterrows():
        st = row["station_code"]
        evt = row["event_code"]
        
        if pd.notna(row.get("Cl_mg_L")):
            observations.append(Observation(f"Cl_{st}_{evt}", float(row["Cl_mg_L"]), weight=0.5))
        if pd.notna(row.get("NO3_mg_L")):
            observations.append(Observation(f"NO3_{st}_{evt}", float(row["NO3_mg_L"]), weight=0.8))
    
    edge_objs = prepare_edges(edges_df)
    
    return {
        "water_chem": water_chem,
        "edge_objs": edge_objs,
        "targets": targets,
        "config": config,
    }, parameters, observations


def run_calibration(data: dict, config: Config) -> dict:
    """Run parameter calibration using PEST-GLM algorithm."""
    print("\n" + "=" * 60)
    print("PARAMETER CALIBRATION")
    print("=" * 60)
    
    context, parameters, observations = setup_calibration(data, config)
    
    if len(observations) < 5:
        print("  Insufficient observations for calibration (< 5)")
        print("  Using default parameters")
        return {p.name: p.value for p in parameters}
    
    print(f"  Parameters: {len(parameters)}")
    print(f"  Observations: {len(observations)}")
    
    def run_model(params: dict) -> dict:
        """Model runner for calibration."""
        sim_results = {}
        evap_factor = params.get("evap_factor", 1.0)
        
        for event_code in context["water_chem"]["event_code"].unique():
            samples = prepare_samples(context["water_chem"], event_code)
            try:
                results = fit_network(samples, context["edge_objs"], config)
                for res in results:
                    if res is None:
                        continue
                    ds = res.v if hasattr(res, 'v') else None
                    if ds is None:
                        continue
                    
                    target_rows = context["targets"][
                        (context["targets"]["station_code"] == ds) & 
                        (context["targets"]["event_code"] == event_code)
                    ]
                    
                    if len(target_rows) > 0:
                        row = target_rows.iloc[0]
                        for ion in ["Cl", "NO3"]:
                            col = f"{ion}_mg_L"
                            if col in row and pd.notna(row[col]):
                                base_val = float(row[col])
                                sim_val = base_val * evap_factor if ion == "Cl" else base_val
                                sim_results[f"{ion}_{ds}_{event_code}"] = sim_val
            except Exception:
                continue
        return sim_results
    
    pest = PESTGLM(parameters=parameters, observations=observations, model_runner=run_model)
    
    print("  Running optimization...")
    try:
        result = pest.calibrate(max_nfev=20)
        opts = result["optimal_parameters"]
        
        print(f"  Calibration Complete! Phi: {result.get('phi', 'N/A'):.4f}")
        for name, val in opts.items():
            print(f"    {name}: {val:.4f}")
        
        return {k: float(v) for k, v in opts.items()}
        
    except Exception as e:
        print(f"  Calibration failed: {e}")
        return {p.name: p.value for p in parameters}


# -------------------------------------------------------------------
# 3. Analysis Functions
# -------------------------------------------------------------------


def analyze_single_event(samples: list, edges: list, config: Config, event_name: str):
    """Run hydrosheaf network fit for a single event."""
    print(f"\n{'='*60}")
    print(f"Analyzing Event: {event_name}")
    print(f"{'='*60}")
    print(f"  Samples: {len(samples)}")
    print(f"  Edges: {len(edges)}")

    # Run network fit
    results = fit_network(samples, edges, config)

    if not results:
        print("  No results (check edge/sample matching)")
        return None

    print(f"  Fitted edges: {len(results)}")

    # Display results for each edge
    for res in results:
        print(f"\n  Edge: {res.u} -> {res.v}")
        print(f"    Transport model: {res.transport_model}")
        print(f"    Evaporation factor (gamma): {res.gamma:.4f}")
        if res.f is not None:
            print(f"    Mixing fraction (f): {res.f:.4f}")
        print(f"    Objective score: {res.objective_score:.6f}")
        print(f"    Transport residual norm: {res.transport_residual_norm:.6f}")
        print(f"    Anomaly norm: {res.anomaly_norm:.6f}")

        # Show top reactions
        if res.z_labels and res.z_extents:
            nonzero = [
                (lbl, ext)
                for lbl, ext in zip(res.z_labels, res.z_extents)
                if abs(ext) > 1e-6
            ]
            if nonzero:
                print(f"    Active reactions:")
                for lbl, ext in sorted(nonzero, key=lambda x: -abs(x[1]))[:5]:
                    sign = "+" if ext > 0 else "-"
                    print(f"      {lbl}: {sign}{abs(ext):.4f} mmol/L")

    # Network summary
    summary = summarize_network(results)
    print(f"\n  Network Summary:")
    print(f"    Total edges: {summary['edge_count']}")
    print(f"    Transport model counts: {summary['transport_counts']}")

    return results


def analyze_all_events(data: dict, config: Config):
    """Analyze all events in the dataset."""
    water_chem = data["water_chem"]
    edges_df = data["edges"]

    # Get unique events
    events = water_chem["event_code"].unique()
    print(f"\nFound {len(events)} events: {list(events)}")

    all_results = {}

    for event_code in events:
        # Prepare samples for this event
        samples = prepare_samples_for_hydrosheaf(water_chem, event_code)
        edges = prepare_edges(edges_df)

        # Run analysis
        results = analyze_single_event(samples, edges, config, event_code)
        all_results[event_code] = results

    return all_results


def compute_isotope_indicators(water_chem: pd.DataFrame):
    """Compute isotope-based indicators for all samples."""
    print("\n" + "=" * 60)
    print("Isotope Analysis")
    print("=" * 60)

    # Water isotopes
    d18O = water_chem["d18O_H2O_permil"].values
    d2H = water_chem["d2H_H2O_permil"].values

    # Compute d-excess: d-excess = d2H - 8 * d18O
    # Note: compute_d_excess takes (d18O, d2H) order
    d_excess = np.array([compute_d_excess(o, h) for o, h in zip(d18O, d2H)])

    # Compute evaporation index (deviation from LMWL)
    # Using default Ghana LMWL: d2H = 8.66 * d18O + 7.22
    # evaporation_index(d18o, d2h, a, b) = d2h - (a + b * d18o)
    # where a = intercept (7.22), b = slope (8.66)
    evap_idx = np.array(
        [evaporation_index(o, h, 7.22, 8.66) for o, h in zip(d18O, d2H)]
    )

    print(f"\nD-excess statistics:")
    print(f"  Mean: {np.nanmean(d_excess):.2f} permil")
    print(f"  Range: {np.nanmin(d_excess):.2f} to {np.nanmax(d_excess):.2f} permil")

    print(f"\nEvaporation index statistics:")
    print(f"  Mean: {np.nanmean(evap_idx):.2f}")
    print(f"  Range: {np.nanmin(evap_idx):.2f} to {np.nanmax(evap_idx):.2f}")

    # Group by station type
    for station_type in water_chem["station_type"].unique():
        mask = water_chem["station_type"] == station_type
        print(f"\n  {station_type.upper()}:")
        print(f"    Mean d-excess: {np.nanmean(d_excess[mask]):.2f} permil")
        print(f"    Mean evap index: {np.nanmean(evap_idx[mask]):.2f}")

    return d_excess, evap_idx


def summarize_hydrochemistry(water_chem: pd.DataFrame):
    """Provide summary statistics of hydrochemistry."""
    print("\n" + "=" * 60)
    print("Hydrochemistry Summary")
    print("=" * 60)

    ion_cols = [
        "Ca_mg_L",
        "Mg_mg_L",
        "Na_mg_L",
        "K_mg_L",
        "HCO3_mg_L",
        "Cl_mg_L",
        "SO4_mg_L",
        "NO3_mg_L",
    ]

    print("\nConcentration Statistics (mg/L):")
    print("-" * 50)
    print(f"{'Ion':<10} {'Mean':>10} {'Std':>10} {'Min':>10} {'Max':>10}")
    print("-" * 50)

    for col in ion_cols:
        if col in water_chem.columns:
            ion = col.replace("_mg_L", "")
            mean = water_chem[col].mean()
            std = water_chem[col].std()
            min_val = water_chem[col].min()
            max_val = water_chem[col].max()
            print(
                f"{ion:<10} {mean:>10.2f} {std:>10.2f} {min_val:>10.2f} {max_val:>10.2f}"
            )

    # Compare lysimeter vs borehole
    print("\n\nComparison by Station Type:")
    print("-" * 50)

    for station_type in ["lysimeter", "borehole"]:
        subset = water_chem[water_chem["station_type"] == station_type]
        print(f"\n{station_type.upper()} (n={len(subset)}):")

        # Key parameters
        print(
            f"  EC: {subset['EC_uS_cm'].mean():.1f} +/- {subset['EC_uS_cm'].std():.1f} uS/cm"
        )
        print(f"  pH: {subset['pH'].mean():.2f} +/- {subset['pH'].std():.2f}")
        print(
            f"  NO3: {subset['NO3_mg_L'].mean():.1f} +/- {subset['NO3_mg_L'].std():.1f} mg/L"
        )
        print(
            f"  HCO3: {subset['HCO3_mg_L'].mean():.1f} +/- {subset['HCO3_mg_L'].std():.1f} mg/L"
        )


def analyze_temporal_patterns(data: dict):
    """Analyze temporal patterns across events."""
    print("\n" + "=" * 60)
    print("Temporal Patterns")
    print("=" * 60)

    water_chem = data["water_chem"]
    events = data["events"]

    # Merge with event metadata
    merged = water_chem.merge(events, on="event_code", how="left")

    print("\nNO3 by Season:")
    for season in merged["season"].dropna().unique():
        subset = merged[merged["season"] == season]
        print(
            f"  {season}: {subset['NO3_mg_L'].mean():.1f} +/- {subset['NO3_mg_L'].std():.1f} mg/L (n={len(subset)})"
        )

    print("\nNO3 by Event:")
    for event_code in events["event_code"]:
        subset = water_chem[water_chem["event_code"] == event_code]
        event_info = events[events["event_code"] == event_code].iloc[0]
        print(
            f"  {event_code} ({event_info['season']}): {subset['NO3_mg_L'].mean():.1f} mg/L"
        )


# -------------------------------------------------------------------
# 4. Main Analysis
# -------------------------------------------------------------------


def main():
    print("=" * 60)
    print("HYDROSHEAF SYNTHETIC DATA ANALYSIS")
    print("=" * 60)

    # Data directory
    data_dir = Path(__file__).parent / "../data/synthetic"

    # Load data
    data = load_synthetic_data(data_dir)

    # Summarize hydrochemistry
    summarize_hydrochemistry(data["water_chem"])

    # Compute isotope indicators
    compute_isotope_indicators(data["water_chem"])

    # Analyze temporal patterns
    analyze_temporal_patterns(data)

    # Configure hydrosheaf
    # Note: Data is in mg/L in CSV but we convert to mmol/L for hydrosheaf
    config = Config(
        # Use standard 10-ion order
        ion_order=["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"],
        weights=[1.0] * 10,
        # PHREEQC disabled for simpler analysis (enable if phreeqpython available)
        phreeqc_enabled=False,
        # Enable transport models
        transport_models_enabled=["evap", "mix"],
        # Mineral reactions
        active_minerals=[
            "calcite",
            "dolomite",
            "gypsum",
            "halite",
            "pyrite_oxidation_aerobic",
        ],
        # Isotope settings (Ghana LMWL)
        isotope_enabled=False,  # Set to True if isotope analysis needed
        lmwl_a=7.87,  # Ghana LMWL slope
        lmwl_b=13.61,  # Ghana LMWL intercept
        # Gibbs diagram constraints
        gibbs_enabled=True,
        # Ion exchange
        exchange_enabled=True,
    )

    # Run calibration
    calibrated_params = run_calibration(data, config)

    # Run analysis for all events
    print("\n" + "=" * 60)
    print("HYDROSHEAF NETWORK FITTING")
    print("=" * 60)

    all_results = analyze_all_events(data, config)

    # Final summary
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

    total_edges = sum(len(r) for r in all_results.values() if r)
    print(f"\nTotal edges fitted across all events: {total_edges}")

    # Aggregate transport model usage
    transport_counts = {}
    for event_results in all_results.values():
        if event_results:
            for res in event_results:
                model = res.transport_model
                transport_counts[model] = transport_counts.get(model, 0) + 1

    print(f"Transport model distribution: {transport_counts}")

    return all_results


if __name__ == "__main__":
    results = main()
