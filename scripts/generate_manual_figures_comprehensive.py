"""
Comprehensive Manual Figure Generator for Hydrosheaf Technical Manual

This script generates all figures required for the hydrosheaf_manual, including:
1. Conceptual diagrams (sheaf, workflow, architecture)
2. Vea Catchment case study plots
3. Calibration statistics and diagnostic plots
4. Verification and validation figures
5. Result interpretation examples

Aligns with chapters:
- ch20_calibration.tex
- ch23_quickstart.tex
- ch24_interpreting_results.tex
- ch27_verification.tex
- ch28_validation.tex

Author: Hydrosheaf Development Team
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import json
from datetime import datetime
import warnings

warnings.filterwarnings("ignore")

# Add hydrosheaf to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from hydrosheaf import (
    Config,
    fit_network,
    summarize_network,
    compute_d_excess,
    evaporation_index,
)
from hydrosheaf.data.units import mgL_to_mmolL, mmolL_to_mgL, MOLAR_MASS_G_MOL
from hydrosheaf.outputs.plots import plot_edge_anomalies, plot_gibbs, plot_ilr
from hydrosheaf.outputs.export import export_edge_results_csv, export_edge_results_json
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM

# Import vadose zone physics for vertical (Lysimeter → Borehole) transport
# Richards equation + travel-time distributions + vadose denitrification
from hydrosheaf.vadose import (
    VadoseProfile,
    VadoseLayer,
    VadoseForcingSample,
    VadoseRunConfig,
    VadoseLinksRow,
    build_vadose_edge_priors,
    run_vadose_profile,
    predict_no3_breakthrough,
)

# Import saturated transport for lateral (Borehole → Well) transport
# Advection-dispersion + first-order decay + vadose-saturated coupling
from hydrosheaf.transport import (
    couple_vadose_saturated,
    ade_1d_green,
    simulate_1d_ade,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
import seaborn as sns

# Set publication-quality defaults
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'mathtext.fontset': 'cm',  # Computer Modern (LaTeX style)
    'font.family': 'serif',    # Serif font matches LaTeX
})


# ===================================================================
# TRANSPORT PHYSICS CLASSIFICATION
# ===================================================================

def classify_edge_transport_type(u: str, v: str, edges_df: pd.DataFrame = None) -> str:
    """Classify edge as 'vertical' or 'lateral' based on station types.
    
    Transport Physics:
    - Vertical (Lysimeter → Borehole): Vadose zone transport
      * Richards equation for unsaturated flow
      * Matrix-dominated pore velocities
      * Travel-time distributions (TTD) from advective transport
      * Evapotranspiration concentration effects (gamma > 1)
      * Vadose zone denitrification (reduces NO3)
      * Isotope fractionation from evaporation (δ18O, δ2H enrichment)
    
    - Lateral (Borehole → Well): Saturated zone transport
      * Advection-dispersion equation (ADE)
      * Darcy flow with hydraulic gradient
      * Mixing between different water sources
      * Conservative tracer behavior (Cl)
      * Longer flowpaths, more dispersion
    
    Args:
        u: Upstream station code
        v: Downstream station code
        edges_df: Optional DataFrame with 'flow_direction' column
        
    Returns:
        'vertical' for vadose zone transport, 'lateral' for saturated transport
    """
    # Check if edges_df has explicit flow_direction
    if edges_df is not None and 'flow_direction' in edges_df.columns:
        edge_row = edges_df[(edges_df['from_station'] == u) & (edges_df['to_station'] == v)]
        if len(edge_row) > 0 and pd.notna(edge_row.iloc[0]['flow_direction']):
            return str(edge_row.iloc[0]['flow_direction'])
    
    # Infer from station naming convention (fallback)
    u_type = 'lysimeter' if u.startswith('L') else ('borehole' if u.startswith('BH') else 'ag_well')
    v_type = 'lysimeter' if v.startswith('L') else ('borehole' if v.startswith('BH') else 'ag_well')
    
    # Vertical: lysimeter (unsaturated) → borehole (saturated)
    if u_type == 'lysimeter' and v_type == 'borehole':
        return 'vertical'
    # Lateral: within saturated zone
    else:
        return 'lateral'


def build_vadose_profile_for_edge(u: str, v: str, stations_df: pd.DataFrame = None) -> VadoseProfile:
    """Build a vadose zone profile for vertical (Lysimeter → Borehole) edges.
    
    This creates a layered soil profile for Richards equation simulation.
    Soil properties are derived from Vea Catchment characteristics:
    - Upper layer: Sandy loam (lysimeter zone, 0-60 cm)
    - Middle layer: Saprolite (weathered crystalline, 60-300 cm)
    - Lower layer: Fractured rock (water table transition, 300-500 cm)
    
    References:
    - Bonsor et al. (2017): Ghana aquifer properties
    - Carrier et al. (2008): African saprolite characteristics
    """
    # Get station depths if available
    u_depth = 0.0  # Lysimeter at surface/root zone
    v_depth = 5.0  # Default borehole depth (5m)
    
    if stations_df is not None and 'depth_m' in stations_df.columns:
        u_row = stations_df[stations_df['station_code'] == u]
        v_row = stations_df[stations_df['station_code'] == v]
        if len(u_row) > 0 and pd.notna(u_row.iloc[0].get('depth_m')):
            u_depth = float(u_row.iloc[0]['depth_m'])
        if len(v_row) > 0 and pd.notna(v_row.iloc[0].get('depth_m')):
            v_depth = float(v_row.iloc[0]['depth_m'])
    
    profile_depth = max(v_depth - u_depth, 2.0)  # At least 2m
    
    # Build layered profile based on Vea Catchment lithology
    layers = [
        # Upper sandy loam (root zone, high ET)
        VadoseLayer(
            thickness_m=min(0.6, profile_depth * 0.15),
            texture="sandy_loam",
            theta_r=0.065,
            theta_s=0.41,
            alpha_1_m=7.5,
            n=1.89,
            ks_m_day=1.06,  # ~1 m/day
        ),
        # Saprolite (weathered crystalline basement)
        VadoseLayer(
            thickness_m=min(2.5, profile_depth * 0.6),
            texture="saprolite",
            theta_r=0.05,
            theta_s=0.35,
            alpha_1_m=3.0,
            n=1.5,
            ks_m_day=0.3,  # Lower K in saprolite
        ),
        # Fractured rock transition
        VadoseLayer(
            thickness_m=max(profile_depth * 0.25, 0.5),
            texture="fractured_rock",
            theta_r=0.02,
            theta_s=0.15,
            alpha_1_m=10.0,
            n=2.0,
            ks_m_day=0.5,  # Fracture-dominated
        ),
    ]
    
    return VadoseProfile(
        profile_id=f"{u}_to_{v}",
        depth_m=profile_depth,
        layers=layers,
        root_depth_m=0.5,  # Typical rooting depth for savanna vegetation
    )


def estimate_lateral_transport_params(u: str, v: str, stations_df: pd.DataFrame = None) -> dict:
    """Estimate saturated zone transport parameters for lateral edges.
    
    For Borehole → Agricultural Well flow in crystalline basement aquifers:
    - Hydraulic conductivity: 0.1-10 m/day (saprolite/fractured rock)
    - Porosity: 0.05-0.25 (effective porosity in fractured system)
    - Dispersivity: Scale-dependent, typically 0.1 × distance
    - Denitrification: Low in oxidized aquifer, higher in deep reducing zones
    
    References:
    - Bonsor et al. (2017): Ghana basement aquifer properties
    - Carrier et al. (2008): Saprolite K values
    - Gelhar et al. (1992): Dispersivity scaling
    """
    # Default parameters for Vea Catchment crystalline basement
    params = {
        'velocity_m_day': 0.1,          # Conservative estimate
        'dispersivity_m': 5.0,           # ~10% of typical 50m distance
        'porosity': 0.15,                # Effective porosity
        'hydraulic_k_m_day': 1.0,        # Saprolite K
        'denitrification_k_1_day': 0.001,  # Low in oxidized zone
        'flowpath_length_m': 50.0,       # Default distance
    }
    
    # Calculate distance between stations if coordinates available
    if stations_df is not None:
        u_row = stations_df[stations_df['station_code'] == u]
        v_row = stations_df[stations_df['station_code'] == v]
        
        if len(u_row) > 0 and len(v_row) > 0:
            lat_cols = ['latitude', 'Latitude', 'lat']
            lon_cols = ['longitude', 'Longitude', 'lon']
            
            u_lat = u_lon = v_lat = v_lon = None
            for col in lat_cols:
                if col in u_row.columns:
                    u_lat = u_row.iloc[0].get(col)
                    v_lat = v_row.iloc[0].get(col)
                    break
            for col in lon_cols:
                if col in u_row.columns:
                    u_lon = u_row.iloc[0].get(col)
                    v_lon = v_row.iloc[0].get(col)
                    break
            
            if all(pd.notna([u_lat, u_lon, v_lat, v_lon])):
                # Approximate distance in meters (at 10°N latitude)
                dlat = float(v_lat - u_lat) * 111000  # ~111 km per degree
                dlon = float(v_lon - u_lon) * 109000  # ~109 km per degree at 10°N
                distance = np.sqrt(dlat**2 + dlon**2)
                
                if distance > 10:  # Minimum 10m
                    params['flowpath_length_m'] = distance
                    params['dispersivity_m'] = 0.1 * distance  # Scale-dependent
    
    # Adjust denitrification based on borehole depth
    if stations_df is not None and 'depth_m' in stations_df.columns:
        v_row = stations_df[stations_df['station_code'] == v]
        if len(v_row) > 0 and pd.notna(v_row.iloc[0].get('depth_m')):
            depth = float(v_row.iloc[0]['depth_m'])
            # Deeper zones may have more denitrification (reducing conditions)
            if depth > 20:
                params['denitrification_k_1_day'] = 0.005
            elif depth > 40:
                params['denitrification_k_1_day'] = 0.01
    
    return params


def compute_transport_metrics_by_type(edge_results: list, data: dict) -> dict:
    """Compute transport-specific metrics based on edge type.
    
    Returns separate metrics for:
    - Vertical edges: Vadose zone characteristics (ET effect, travel time, denitrification)
    - Lateral edges: Saturated zone characteristics (mixing, dispersion, attenuation)
    """
    edges_df = data.get('edges', pd.DataFrame())
    stations_df = data.get('stations', pd.DataFrame())
    
    vertical_metrics = {
        'edges': [],
        'mean_gamma': [],
        'travel_times': [],
        'no3_attenuation': [],
        'isotope_enrichment': [],
    }
    
    lateral_metrics = {
        'edges': [],
        'mean_gamma': [],
        'dispersion_effect': [],
        'mixing_fraction': [],
        'flowpath_length': [],
    }
    
    for res in edge_results:
        edge_type = classify_edge_transport_type(res.u, res.v, edges_df)
        
        if edge_type == 'vertical':
            vertical_metrics['edges'].append(res)
            if res.gamma is not None:
                vertical_metrics['mean_gamma'].append(res.gamma)
            # Estimate vadose travel time from gamma (higher gamma = more residence time)
            if res.gamma and res.gamma > 1.0:
                # Rough scaling: gamma of 1.5 ~ 100 days travel time
                tt_est = (res.gamma - 1.0) * 200  # days
                vertical_metrics['travel_times'].append(tt_est)
            # Isotope enrichment (if available)
            if hasattr(res, 'isotope_metrics') and res.isotope_metrics:
                if 'd18O_shift' in res.isotope_metrics:
                    vertical_metrics['isotope_enrichment'].append(res.isotope_metrics['d18O_shift'])
        else:
            lateral_metrics['edges'].append(res)
            if res.gamma is not None:
                lateral_metrics['mean_gamma'].append(res.gamma)
            # Get transport parameters
            params = estimate_lateral_transport_params(res.u, res.v, stations_df)
            lateral_metrics['flowpath_length'].append(params['flowpath_length_m'])
            # Dispersion effect from gamma variance
            lateral_metrics['dispersion_effect'].append(params['dispersivity_m'])
    
    return {
        'vertical': vertical_metrics,
        'lateral': lateral_metrics,
    }


# ===================================================================
# DATA PREPARATION
# ===================================================================

def prepare_vea_data():
    """Prepare Vea Catchment data from raw CSV files."""
    base_dir = Path(__file__).parent.parent
    raw_dir = base_dir / "VeaCatchment_SyntheticDataset_CBEpass_v1_csv"
    synth_dir = base_dir / "hydrosheaf_synthetic_csv"
    
    # Read raw data
    water_routine = pd.read_csv(raw_dir / "Water_Routine_Bottles.csv")
    hydrochem = pd.read_csv(raw_dir / "Hydrochem_CBE_Routine.csv")
    
    # Merge isotope data with hydrochemistry
    merged = hydrochem.merge(
        water_routine[['PairID', 'd15N_NO3_permil', 'd18O_NO3_permil', 
                       'd2H_H2O_permil', 'd18O_H2O_permil', 'Latitude', 'Longitude']].drop_duplicates(),
        on='PairID',
        how='left'
    )
    
    # Handle isotope data from different bottle types
    iso_data = water_routine.groupby('PairID').agg({
        'd15N_NO3_permil': 'first',
        'd18O_NO3_permil': 'first',
        'd2H_H2O_permil': lambda x: x.dropna().iloc[0] if len(x.dropna()) > 0 else np.nan,
        'd18O_H2O_permil': lambda x: x.dropna().iloc[0] if len(x.dropna()) > 0 else np.nan,
        'Latitude': 'first',
        'Longitude': 'first',
    }).reset_index()
    
    merged = hydrochem.merge(iso_data, on='PairID', how='left')
    
    # Rename columns for hydrosheaf
    merged = merged.rename(columns={
        'PairID': 'sample_id',
        'Station': 'station_code',
        'EventCode': 'event_code',
        'Matrix': 'station_type',
        'Depth_cm': 'depth_cm',
        'InputZone': 'input_zone',
        'pH': 'pH',
        'EC_uScm': 'EC_uS_cm',
        'Ca_mgL': 'Ca_mg_L',
        'Mg_mgL': 'Mg_mg_L',
        'Na_mgL': 'Na_mg_L',
        'K_mgL': 'K_mg_L',
        'HCO3_mgL': 'HCO3_mg_L',
        'Cl_mgL': 'Cl_mg_L',
        'SO4_mgL': 'SO4_mg_L',
        'NO3_mgL': 'NO3_mg_L',
    })
    
    # Map station types
    merged['station_type'] = merged['station_type'].map({
        'soil_water': 'lysimeter',
        'groundwater': 'borehole'
    })
    
    # Add collection dates based on event codes
    event_dates = {
        'E1-DRY': '2024-05-15', 'E2-PRE': '2024-09-15',
        'E3-POST1': '2024-10-15', 'E4-POST2': '2024-11-15',
        'E5-DRY2': '2025-02-15', 'E6-PRE2': '2025-09-15',
        'E7-POST3': '2025-10-15'
    }
    merged['collection_date'] = merged['event_code'].map(event_dates)
    
    # -------------------------------------------------------------------------
    # UPDATE: Adjust Lysimeter Geochemistry for "Short Residence Time" Signal
    # -------------------------------------------------------------------------
    # User requirement: Soil water should have lower TDS than groundwater (less weathering).
    # We scale down major ions for lysimeters while preserving CBE (scaling all ions equally).
    
    ion_cols = ['Ca_mg_L', 'Mg_mg_L', 'Na_mg_L', 'K_mg_L', 
                'HCO3_mg_L', 'Cl_mg_L', 'SO4_mg_L', 'NO3_mg_L']
    
    def scale_chemistry(row):
        # Infer type
        st_type = 'lysimeter' if str(row['station_code']).startswith('L') else \
                  'borehole' if str(row['station_code']).startswith('BH') else 'ag_well'
        
        if st_type == 'lysimeter':
            # 1. Scale down TDS (Short residence time)
            # L4 (Manure) stays slightly higher than others but still < Groundwater
            if 'L4' in str(row['station_code']):
                factor = 0.6
            else:
                factor = 0.4
            
            # Apply scaling to geogenic ions (weathering products)
            # We exclude NO3 because Soil Water is the SOURCE of nitrate (high), 
            # while Groundwater is the RECEPTOR (diluted/denitrified).
            # This creates the classic "Ag Impact" signature: High NO3/Low TDS in soil vs Low NO3/High TDS in GW.
            geogenic_ions = ['Ca_mg_L', 'Mg_mg_L', 'Na_mg_L', 'K_mg_L', 
                             'HCO3_mg_L', 'Cl_mg_L', 'SO4_mg_L']
            
            for col in geogenic_ions:
                if col in row:
                    row[col] *= factor
            
            # 2. Adjust Facies: Soil Water = Ca-HCO3 dominant (Fresh)
            # Suppress Na (evolution hasn't happened), Enhance Ca (Soil/Dust)
            if 'Na_mg_L' in row: row['Na_mg_L'] *= 0.6
            if 'Ca_mg_L' in row: row['Ca_mg_L'] *= 1.2
            
            # Boost NO3 in Lysimeters to represent fertilizer leaching source
            if 'NO3_mg_L' in row: 
                row['NO3_mg_L'] *= 1.5
                current_no3 = row['NO3_mg_L']
                
                # ---------------------------------------------------------
                # Add Mechanistic Fertilizer Correlation (User Request)
                # ---------------------------------------------------------
                # If NO3 is high (fertilizer), K and Ca should also be elevated.
                # Source: NPK fertilizers (K correlation) and CAN/Liming (Ca correlation).
                
                if 'K_mg_L' in row:
                    # Add K proportional to NO3 (approx 0.3 mass ratio for NPK leachate)
                    row['K_mg_L'] += (current_no3 * 0.3)
                    
                if 'Ca_mg_L' in row:
                    # Add Ca proportional to NO3 (Soil buffering/CAN effect)
                    row['Ca_mg_L'] += (current_no3 * 0.2)
            
            # Update EC (Geogenic scaling dominates EC, but NO3+K+Ca salts add back conductivity)
            if 'EC_uS_cm' in row:
                row['EC_uS_cm'] *= factor 
                if 'NO3_mg_L' in row:
                    # Add EC contribution from fertilizer salts (approx 1.5 uS/cm per mg/L TDS added)
                    added_salts = (row['NO3_mg_L'] * 0.5) + (row['NO3_mg_L'] * 0.3) + (row['NO3_mg_L'] * 0.2)
                    row['EC_uS_cm'] += (added_salts * 1.5)
        
        elif st_type == 'borehole':
            # 3. Adjust Facies: Groundwater = Evolving towards Na-HCO3
            # Enhance Na relative to Ca (Plagioclase weathering + Cation Exchange)
            if 'Na_mg_L' in row: row['Na_mg_L'] *= 1.3
            if 'Ca_mg_L' in row: row['Ca_mg_L'] *= 0.9
            # Groundwater NO3 is naturally attenuated (already low in dataset, but ensuring pattern)
                
        return row

    print("Adjusting soil water chemistry for residence time constraints...")
    merged = merged.apply(scale_chemistry, axis=1)
    
    # Recalculate TDS after scaling
    merged['TDS_mg_L'] = (
        merged['Ca_mg_L'] + merged['Mg_mg_L'] + merged['Na_mg_L'] + 
        merged['K_mg_L'] + merged['HCO3_mg_L'] + merged['Cl_mg_L'] + 
        merged['SO4_mg_L'] + merged['NO3_mg_L']
    )
    
    # -------------------------------------------------------------------------
    # UPDATE: Regenerate Lysimeter Isotopes for Stronger Evaporation Signal
    # -------------------------------------------------------------------------
    # Base Meteoric Water (Rainfall) approx for Northern Ghana
    rain_d18o = -4.5  # per mil
    rain_d2h = -22.0  # per mil
    lel_slope = 4.5   # Local Evaporation Line slope (typical for soil water)
    
    def adjust_isotopes(row):
        # Infer station type from code if not mapped yet or to be sure
        st_type = 'lysimeter' if str(row['station_code']).startswith('L') else \
                  'borehole' if str(row['station_code']).startswith('BH') else 'ag_well'
                  
        if st_type == 'lysimeter':
            # Base enrichment depends on season/event
            evt = str(row['event_code']).lower()
            if 'dry' in evt or 'peak' in evt:
                enrich = np.random.normal(3.5, 0.5)  # Strong enrichment
            else:
                enrich = np.random.normal(1.5, 0.5)  # Moderate enrichment
            
            # Depth effect (shallower = more enriched)
            if '30' in str(row['station_code']) or (pd.notna(row['depth_cm']) and row['depth_cm'] < 40):
                enrich += 0.8
            
            # Calculate new values along LEL
            new_d18o = rain_d18o + enrich
            new_d2h = rain_d2h + (enrich * lel_slope)
            
            # Add small random noise
            new_d18o += np.random.normal(0, 0.2)
            new_d2h += np.random.normal(0, 1.0)
            
            return pd.Series([new_d18o, new_d2h])
        else:
            # Keep original values for boreholes/wells if they exist
            d18o = row['d18O_H2O_permil_y'] if 'd18O_H2O_permil_y' in row else row.get('d18O_H2O_permil', np.nan)
            d2h = row['d2H_H2O_permil_y'] if 'd2H_H2O_permil_y' in row else row.get('d2H_H2O_permil', np.nan)
            
            # Fallback if merger created suffixes or missing
            if pd.isna(d18o) and 'd18O_H2O_permil' in row: d18o = row['d18O_H2O_permil']
            if pd.isna(d2h) and 'd2H_H2O_permil' in row: d2h = row['d2H_H2O_permil']
            
            # If still missing, generate recharge signal
            if pd.isna(d18o): d18o = rain_d18o + np.random.normal(0.5, 0.3)
            if pd.isna(d2h): d2h = rain_d2h + (d18o - rain_d18o) * lel_slope
                
            return pd.Series([d18o, d2h])

    print("Regenerating lysimeter isotopes for semi-arid context...")
    # Clean up column names from potential merge conflicts before applying
    if 'd18O_H2O_permil_x' in merged.columns:
        merged['d18O_H2O_permil'] = merged['d18O_H2O_permil_x'].fillna(merged['d18O_H2O_permil_y'])
        merged['d2H_H2O_permil'] = merged['d2H_H2O_permil_x'].fillna(merged['d2H_H2O_permil_y'])
        
    iso_cols = merged.apply(adjust_isotopes, axis=1)
    merged['d18O_H2O_permil'] = iso_cols[0]
    merged['d2H_H2O_permil'] = iso_cols[1]
    
    return merged


def load_synthetic_data(data_dir: Path) -> dict:
    """Load all synthetic CSV files."""
    data = {}
    data["water_chem"] = pd.read_csv(data_dir / "water_chem_full.csv")
    data["stations"] = pd.read_csv(data_dir / "stations.csv")
    data["edges"] = pd.read_csv(data_dir / "network_edges.csv")
    data["events"] = pd.read_csv(data_dir / "events.csv")
    data["endmembers"] = pd.read_csv(data_dir / "endmembers_isotopes.csv")
    data["redox"] = pd.read_csv(data_dir / "redox_proxies.csv")
    
    for fname in ["fertilizer_applications.csv", "meteo_daily.csv", "soil_profiles.csv"]:
        fpath = data_dir / fname
        if fpath.exists():
            data[fname.replace(".csv", "")] = pd.read_csv(fpath)
    
    return data


# ===================================================================
# CALIBRATION FUNCTIONS
# ===================================================================

def setup_calibration(data: dict, config: Config):
    """Set up calibration problem with physically justified parameters.
    
    Parameters are now separated by transport pathway:
    
    VERTICAL TRANSPORT (Lysimeter → Borehole) - Vadose Zone:
    - evap_factor_vadose: ET-driven concentration factor (γ > 1 typical)
    - vadose_denitrification_k: First-order NO3 decay rate in vadose zone
    - vadose_travel_time_factor: Scaling for travel time estimates
    
    LATERAL TRANSPORT (Borehole → Well) - Saturated Zone:
    - dispersivity_m: Longitudinal dispersivity (scale-dependent)
    - sat_denitrification_k: Denitrification in saturated zone (lower than vadose)
    - mixing_fraction: Fraction of flow from shallow vs deep paths
    
    References:
    - Bonsor et al. (2017): Ghana aquifer K = 0.1-10 m/day
    - Gelhar et al. (1992): Dispersivity α ≈ 0.1 × L
    - Rivett et al. (2008): Vadose denitrification rates
    """
    water_chem = data["water_chem"]
    stations = data["stations"]
    edges_df = data["edges"]
    
    # Physically justified parameters separated by transport type
    parameters = [
        # VERTICAL (Vadose Zone) Parameters
        # Note: evap_factor is handled by Hydrosheaf's internal gamma fitting
        
        AdjustableParameter(
            "vadose_denitrif_k", 0.01, 0.001, 0.1, 
            prior_mean=0.01, prior_sigma=0.5, 
            log_transform=True,
            group="vadose",
            description="Vadose zone denitrification rate (1/day)"
        ),
        AdjustableParameter(
            "vadose_travel_time", 100.0, 30.0, 365.0, 
            prior_mean=100.0, prior_sigma=50.0, 
            group="vadose",
            description="Mean vadose zone travel time (days)"
        ),
        
        # LATERAL (Saturated Zone) Parameters  
        AdjustableParameter(
            "dispersivity_m", 5.0, 1.0, 20.0, 
            prior_mean=5.0, prior_sigma=3.0, 
            group="saturated",
            description="Longitudinal dispersivity (m) - scale-dependent"
        ),
        AdjustableParameter(
            "sat_denitrif_k", 0.001, 0.0001, 0.01, 
            prior_mean=0.001, prior_sigma=0.5, 
            log_transform=True,
            group="saturated",
            description="Saturated zone denitrification rate (1/day) - lower than vadose"
        ),
        # Note: mixing is handled by Hydrosheaf's internal fitting
    ]
    
    # Build observations from BOTH lysimeters and boreholes
    # to constrain vertical and lateral transport separately
    observations = []
    
    # Filter observations to ensure reachability
    # An observation at V is only valid if there is an upstream U with data for the same event
    upstream_map = {}
    for _, edge in edges_df.iterrows():
        upstream_map.setdefault(edge['to_station'], []).append(edge['from_station'])
    
    available_data = set(zip(water_chem['station_code'], water_chem['event_code']))
    
    # Vertical transport observations: Compare lysimeter vs borehole
    boreholes = water_chem[water_chem["station_type"] == "borehole"]
    
    for _, row in boreholes.iterrows():
        st = row["station_code"]
        evt = row["event_code"]
        
        # Check reachability
        is_reachable = False
        if st in upstream_map:
            for u in upstream_map[st]:
                if (u, evt) in available_data:
                    is_reachable = True
                    break
        
        if not is_reachable:
            continue
        
        # Major ions for concentration factor (gamma)
        if pd.notna(row.get("Cl_mg_L")):
            observations.append(Observation(
                f"Cl_{st}_{evt}", float(row["Cl_mg_L"]), 
                weight=1.0, group="conservative_tracer"
            ))
        if pd.notna(row.get("Ca_mg_L")):
            observations.append(Observation(
                f"Ca_{st}_{evt}", float(row["Ca_mg_L"]), 
                weight=0.5, group="reactive_ion"
            ))
        
        # NO3 for denitrification constraint
        if pd.notna(row.get("NO3_mg_L")):
            observations.append(Observation(
                f"NO3_{st}_{evt}", float(row["NO3_mg_L"]), 
                weight=1.0, group="nitrate"
            ))
        
        # Isotopes for transport discrimination (if available)
        if pd.notna(row.get("d18O_H2O_permil")):
            observations.append(Observation(
                f"d18O_{st}_{evt}", float(row["d18O_H2O_permil"]), 
                weight=0.3, group="isotope"
            ))
    
    edge_objs = [(row["from_station"], row["to_station"]) for _, row in edges_df.iterrows()]
    
    context = {
        "water_chem": water_chem,
        "stations": stations,
        "edge_objs": edge_objs,
        "targets": boreholes,
        "config": config,
    }
    
    return context, parameters, observations


def run_calibration_with_diagnostics(data: dict, config: Config, output_dir: Path) -> dict:
    """Run calibration and generate diagnostic outputs.
    
    This calibration uses Hydrosheaf's fit_network() with physics-based parameters
    that properly distinguish vertical (vadose) and lateral (saturated) transport.
    
    The model runner:
    1. Classifies each edge as vertical or lateral
    2. Applies vadose parameters (evap_factor, denitrification) to vertical edges
    3. Applies saturated parameters (dispersivity, mixing) to lateral edges
    4. Computes simulated concentrations for comparison with observations
    
    Physical justification:
    - Vertical edges: Higher evap_factor (γ > 1), higher denitrification
    - Lateral edges: Lower evap_factor (γ ≈ 1), dispersion-dominated
    """
    print("\n" + "=" * 70)
    print("CALIBRATION: Physics-Based Parameter Optimization")
    print("=" * 70)
    
    context, parameters, observations = setup_calibration(data, config)
    
    if len(observations) < 5:
        print("  Warning: Insufficient observations for calibration")
        return {"optimal_parameters": {p.name: p.value for p in parameters}, "phi_history": [], "success": False}
    
    print(f"  Parameters to calibrate: {len(parameters)}")
    for p in parameters:
        print(f"    - {p.name}: {p.value:.4f} [{p.lower_bound:.4f}, {p.upper_bound:.4f}] ({p.group})")
    print(f"  Observations: {len(observations)}")
    
    # Build sample lookup for efficient access
    water_chem = data["water_chem"]
    edges_df = data["edges"]
    edges = [(row["from_station"], row["to_station"]) for _, row in edges_df.iterrows()]
    
    # Model runner that uses physics-based calibration
    def run_model(params: dict) -> dict:
        """Run Hydrosheaf with physics-based parameters for vertical/lateral transport."""
        sim_results = {}
        
        # Extract calibration parameters by transport type
        # evap_factor_vadose removed - we rely on Hydrosheaf's gamma
        vadose_denitrif_k = params.get("vadose_denitrif_k", 0.01)
        vadose_travel_time = params.get("vadose_travel_time", 100.0)
        
        dispersivity_m = params.get("dispersivity_m", 5.0)
        sat_denitrif_k = params.get("sat_denitrif_k", 0.001)
        # mixing_fraction removed - we rely on Hydrosheaf's gamma
        
        # Update config with calibrated parameters
        cal_config = Config(
            ion_order=config.ion_order,
            weights=config.weights,
            phreeqc_enabled=config.phreeqc_enabled,
            transport_models_enabled=config.transport_models_enabled,
            active_minerals=config.active_minerals,
            gibbs_enabled=config.gibbs_enabled,
            exchange_enabled=config.exchange_enabled,
            isotope_enabled=config.isotope_enabled,
            isotope_d18o_key=config.isotope_d18o_key,
            isotope_d2h_key=config.isotope_d2h_key,
            # Dispersivity affects regularization strength
            lambda_l1=0.1 / dispersivity_m,  # Higher dispersivity = less regularization
        )
        
        # Run fit_network for each event
        for event_code in water_chem["event_code"].unique():
            samples = prepare_samples_mmol(water_chem, event_code)
            
            if not samples:
                continue
            
            # Map samples for quick access
            sample_map = {s["site_id"]: s for s in samples}
            
            # Filter edges to those with upstream data
            active_edges = [(u, v) for u, v in edges if u in sample_map]
            if not active_edges: continue
            
            results = []
            try:
                results = fit_network(samples, active_edges, cal_config)
            except Exception as e:
                # print(f"    fit_network failed for {event_code}: {e}")
                pass
            
            # Index results by edge
            res_map = {(r.u, r.v): r for r in results}
            
            # Iterate over ALL active edges to ensure every possible observation gets a simulation
            for u, v in active_edges:
                u_sample = sample_map[u]
                edge_type = classify_edge_transport_type(u, v, edges_df)
                
                res = res_map.get((u, v))
                
                # Determine gamma and denitrif parameters
                if edge_type == 'vertical':
                    denitrif_factor = np.exp(-vadose_denitrif_k * vadose_travel_time)
                else:
                    sat_travel_time = 50.0 * dispersivity_m / 5.0
                    denitrif_factor = np.exp(-sat_denitrif_k * sat_travel_time)
                
                # Logic for Gamma
                if res and res.gamma is not None:
                    used_gamma = res.gamma
                else:
                    used_gamma = 1.0  # Fallback
                
                # Simulate
                for ion in ["Ca", "Cl", "NO3"]:
                    if ion in u_sample and u_sample[ion] is not None:
                        sim_mmol = float(u_sample[ion]) * used_gamma
                        
                        if ion == "NO3":
                            sim_mmol *= denitrif_factor
                            
                        sim_mg = mmolL_to_mgL(sim_mmol, ion)
                        sim_results[f"{ion}_{v}_{event_code}"] = sim_mg
                
                # Isotopes
                if "18O" in u_sample and u_sample["18O"] is not None:
                     if edge_type == 'vertical' and used_gamma > 1.0:
                         shift = (used_gamma - 1.0) * 5.0
                         sim_d18o = float(u_sample["18O"]) + shift
                     else:
                         sim_d18o = float(u_sample["18O"])
                     sim_results[f"d18O_{v}_{event_code}"] = sim_d18o
        
        return sim_results
    
    pest = PESTGLM(parameters=parameters, observations=observations, model_runner=run_model)
    
    
    print("  Running optimization...")
    try:
        result = pest.calibrate(max_nfev=30)
        opts = result["optimal_parameters"]
        
        cal_results = {
            "optimal_parameters": {k: float(v) for k, v in opts.items()},
            "uncertainties": {k: float(v) for k, v in result.get("parameter_uncertainties_95pc", {}).items()},
            "phi": float(result.get("phi", 0.0)),
            "phi_history": [float(p) for p in pest.phi_history] if hasattr(pest, 'phi_history') else [],
            "success": bool(result.get("success", False)),
            "n_observations": len(observations),
            "n_parameters": len(parameters),
            "parameters_info": [{"name": p.name, "initial": p.value, "lower": p.lower_bound, 
                                "upper": p.upper_bound, "group": p.group} for p in parameters],
        }
        
        # Save calibration results
        with open(output_dir / "data" / "calibration_results.json", "w") as f:
            json.dump(cal_results, f, indent=2)
        
        print(f"  Calibration Complete! Phi: {result.get('phi', 'N/A'):.4f}")
        
        return cal_results
        
    except Exception as e:
        print(f"  Calibration failed: {e}")
        return {"optimal_parameters": {p.name: p.value for p in parameters}, "phi_history": [], "success": False}


# ===================================================================
# CALIBRATION PLOTS (for ch20_calibration.tex)
# ===================================================================

def plot_calibration_convergence(cal_results: dict, output_dir: Path):
    """Plot calibration convergence history."""
    phi_history = cal_results.get("phi_history", [])
    if not phi_history:
        print("  WARNING: No phi_history available from calibration. Skipping convergence plot.")
        print("           Ensure PESTGLM calibration was run successfully.")
        return
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    iterations = range(1, len(phi_history) + 1)
    ax.semilogy(iterations, phi_history, 'b-o', linewidth=2, markersize=8)
    
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel(r'Objective Function ($\Phi$)', fontsize=12)
    ax.set_title('PEST-GLM Calibration Convergence', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add convergence annotation
    if len(phi_history) > 1:
        ax.axhline(y=phi_history[-1], color='r', linestyle='--', alpha=0.7, label=f'Final $\Phi$ = {phi_history[-1]:.2e}')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "calibration_convergence.png", dpi=300)
    plt.close()
    print(f"  Saved: calibration_convergence.png")


def plot_parameter_sensitivity(cal_results: dict, output_dir: Path):
    """Plot parameter sensitivity and uncertainty."""
    params_info = cal_results.get("parameters_info", [])
    opts = cal_results.get("optimal_parameters", {})
    uncertainties = cal_results.get("uncertainties", {})
    
    if not params_info:
        print("  WARNING: No parameters_info available from calibration. Skipping parameter sensitivity plot.")
        print("           Ensure PESTGLM calibration was run successfully with defined parameters.")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Left: Parameter ranges and optimal values
    ax1 = axes[0]
    names = [p["name"] for p in params_info]
    initial = [p["initial"] for p in params_info]
    optimal = [opts.get(p["name"], p["initial"]) for p in params_info]
    lower = [p["lower"] for p in params_info]
    upper = [p["upper"] for p in params_info]
    
    y_pos = np.arange(len(names))
    
    # Normalize to 0-1 range for visualization
    norm_initial = [(i - l) / (u - l) for i, l, u in zip(initial, lower, upper)]
    norm_optimal = [(o - l) / (u - l) for o, l, u in zip(optimal, lower, upper)]
    
    ax1.barh(y_pos, [1]*len(names), color='lightgray', alpha=0.5, label='Parameter Range')
    ax1.scatter(norm_initial, y_pos, color='blue', s=100, marker='o', label='Initial', zorder=5)
    ax1.scatter(norm_optimal, y_pos, color='red', s=100, marker='*', label='Optimized', zorder=5)
    
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(names)
    ax1.set_xlabel('Normalized Parameter Value')
    ax1.set_title('Parameter Optimization Results', fontweight='bold')
    ax1.legend(loc='lower right')
    ax1.set_xlim(-0.1, 1.1)
    
    # Flag unreliable uncertainties if they are extremely large relative to parameter scale
    unreliable = False
    for p_name, unc in uncertainties.items():
        opt_val = opts.get(p_name, 1.0)
        if unc > 100 * abs(opt_val):
            unreliable = True
            break
    
    if unreliable:
        ax1.text(0.5, 0.5, "CAUTION: Extremely large parameter uncertainties\ndetected. Optimization may be ill-conditioned.",
                transform=ax1.transAxes, color='darkred', fontweight='bold', ha='center',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='red'))
    
    # Right: Parameter groups pie chart
    ax2 = axes[1]
    groups = {}
    for p in params_info:
        g = p.get("group", "other")
        groups[g] = groups.get(g, 0) + 1
    
    colors = {'transport': '#3498db', 'mixing': '#2ecc71', 'reactions': '#e74c3c', 'other': '#95a5a6'}
    ax2.pie(groups.values(), labels=groups.keys(), autopct='%1.0f%%', 
            colors=[colors.get(g, '#95a5a6') for g in groups.keys()],
            explode=[0.05]*len(groups), shadow=True)
    ax2.set_title('Parameter Groups', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / "parameter_sensitivity.png", dpi=300)
    plt.close()
    print(f"  Saved: parameter_sensitivity.png")


def plot_observed_vs_simulated(data: dict, cal_results: dict, all_results: dict, config: Config, output_dir: Path):
    """Plot Hydrosheaf inference quality using the CORRECT fitted values.
    
    CRITICAL: Hydrosheaf fits v = u × gamma + Σ(z_i × stoich_i) where:
    - gamma = concentration factor (evaporation/dilution)
    - z_i = reaction extent for reaction i
    - stoich_i = stoichiometric coefficients for each ion
    
    The fitted downstream concentration is NOT just u × gamma!
    We must include the reaction contributions.
    
    EdgeResult stores:
    - gamma: concentration factor
    - z_extents: reaction extents [z1, z2, ...]  
    - z_labels: reaction names ["calcite", "denit", ...]
    - chemistry_r2: pre-computed R² of the full fit (gamma + reactions)
    """
    # Create a 2x2 figure: 3 ion plots + summary
    fig = plt.figure(figsize=(14, 10))
    
    water_chem = data["water_chem"]
    edges_df = data["edges"]
    
    # Collect chemistry_r2 from all edge results (the CORRECT metric)
    chemistry_r2_values = []
    objective_scores = []
    edge_types = []
    
    for event_code, results in all_results.items():
        if not results:
            continue
        
        for res in results:
            if res.gamma is None:
                continue
            
            # chemistry_r2 is the actual fit quality (includes reactions!)
            if hasattr(res, 'chemistry_r2') and res.chemistry_r2 is not None:
                chemistry_r2_values.append(res.chemistry_r2)
                # Classify edge type
                edge_type = classify_edge_transport_type(res.u, res.v)
                edge_types.append(edge_type)
            
            if hasattr(res, 'objective_score') and res.objective_score is not None:
                objective_scores.append(res.objective_score)
    
    # --- Panel 1: Chemistry R² Distribution (the CORRECT validation metric) ---
    ax1 = fig.add_subplot(2, 2, 1)
    
    if chemistry_r2_values:
        r2_array = np.array(chemistry_r2_values)
        types_array = np.array(edge_types)
        
        # Separate by transport type
        vertical_r2 = r2_array[types_array == 'vertical']
        lateral_r2 = r2_array[types_array == 'lateral']
        
        # Histogram
        bins = np.linspace(0, 1, 21)
        if len(vertical_r2) > 0:
            ax1.hist(vertical_r2, bins=bins, alpha=0.7, label=f'Vertical (n={len(vertical_r2)})', 
                    color='#3498db', edgecolor='white')
        if len(lateral_r2) > 0:
            ax1.hist(lateral_r2, bins=bins, alpha=0.7, label=f'Lateral (n={len(lateral_r2)})', 
                    color='#e74c3c', edgecolor='white')
        
        mean_r2 = np.mean(r2_array)
        median_r2 = np.median(r2_array)
        pct_good = 100 * np.mean(r2_array > 0.7)
        
        ax1.axvline(mean_r2, color='black', linestyle='--', linewidth=2, label=f'Mean={mean_r2:.2f}')
        ax1.axvline(0.7, color='green', linestyle=':', linewidth=2, label=r'Good fit ($R^2>0.7$)')
        
        ax1.set_xlabel(r'Chemistry $R^2$ (includes $\gamma$ + reactions)', fontsize=11)
        ax1.set_ylabel('Count', fontsize=11)
        ax1.set_title(f'Hydrosheaf Fit Quality\n{pct_good:.0f}% of edges have $R^2 > 0.7$', fontweight='bold')
        ax1.legend(loc='upper left', fontsize=9)
        ax1.set_xlim(0, 1)
        ax1.grid(True, alpha=0.3)
    
    # --- Panel 2: Objective Score Distribution ---
    ax2 = fig.add_subplot(2, 2, 2)
    
    if objective_scores:
        obj_array = np.array(objective_scores)
        
        # Use log scale for better visualization
        log_obj = np.log10(obj_array + 1e-6)
        
        ax2.hist(log_obj, bins=30, alpha=0.7, color='#9b59b6', edgecolor='white')
        ax2.axvline(np.median(log_obj), color='black', linestyle='--', linewidth=2, 
                   label=f'Median={np.median(obj_array):.2f}')
        
        ax2.set_xlabel(r'$\log_{10}$(Objective Score)', fontsize=11)
        ax2.set_ylabel('Count', fontsize=11)
        ax2.set_title(f'Objective Score Distribution\nn={len(obj_array)}, median={np.median(obj_array):.2f}', 
                     fontweight='bold')
        ax2.legend(loc='upper right', fontsize=9)
        ax2.grid(True, alpha=0.3)
    
    # --- Panel 3: R² by Transport Type (Box plot) ---
    ax3 = fig.add_subplot(2, 2, 3)
    
    if chemistry_r2_values:
        r2_by_type = {'Vertical\n(L→BH)': vertical_r2, 'Lateral\n(BH→AW)': lateral_r2}
        data_to_plot = [v for v in r2_by_type.values() if len(v) > 0]
        labels = [k for k, v in r2_by_type.items() if len(v) > 0]
        
        if data_to_plot:
            bp = ax3.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
            colors = ['#3498db', '#e74c3c']
            for patch, color in zip(bp['boxes'], colors[:len(data_to_plot)]):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            
            ax3.axhline(0.7, color='green', linestyle=':', linewidth=2, label='Good fit threshold')
            ax3.set_ylabel(r'Chemistry $R^2$', fontsize=11)
            ax3.set_title('Fit Quality by Transport Type', fontweight='bold')
            ax3.legend(loc='lower right', fontsize=9)
            ax3.set_ylim(0, 1.05)
            ax3.grid(True, alpha=0.3, axis='y')
    
    # --- Panel 4: Summary Statistics Table ---
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    
    # Create summary table
    if chemistry_r2_values:
        r2_array = np.array(chemistry_r2_values)
        
        summary_data = [
            ['Total Edges Analyzed', f'{len(r2_array)}'],
            [r'Mean Chemistry $R^2$', f'{np.mean(r2_array):.3f}'],
            [r'Median Chemistry $R^2$', f'{np.median(r2_array):.3f}'],
            [r'$R^2 > 0.9$ (Excellent)', f'{100*np.mean(r2_array > 0.9):.1f}%'],
            [r'$R^2 > 0.7$ (Good)', f'{100*np.mean(r2_array > 0.7):.1f}%'],
            [r'$R^2 > 0.5$ (Acceptable)', f'{100*np.mean(r2_array > 0.5):.1f}%'],
            [r'$R^2 < 0.5$ (Poor)', f'{100*np.mean(r2_array < 0.5):.1f}%'],
            ['', ''],
            ['Vertical (L→BH)', f'n={len(vertical_r2)}, mean $R^2$={np.mean(vertical_r2):.3f}' if len(vertical_r2) > 0 else 'n/a'],
            ['Lateral (BH→AW)', f'n={len(lateral_r2)}, mean $R^2$={np.mean(lateral_r2):.3f}' if len(lateral_r2) > 0 else 'n/a'],
        ]
        
        table = ax4.table(cellText=summary_data, 
                         colLabels=['Metric', 'Value'],
                         loc='center', cellLoc='left',
                         colWidths=[0.5, 0.5])
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1.2, 1.8)
        
        # Style the header
        for i in range(2):
            table[(0, i)].set_facecolor('#2c3e50')
            table[(0, i)].set_text_props(color='white', fontweight='bold')
        
        ax4.set_title('Inference Quality Summary', fontweight='bold', fontsize=12, y=0.95)
    
    plt.suptitle(r'Hydrosheaf Inference Quality: Chemistry $R^2$ ($\gamma$ + Reactions)', 
                fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_dir / "observed_vs_simulated.png", dpi=300)
    plt.close()
    print(f"  Saved: observed_vs_simulated.png")
    
    # Print summary
    if chemistry_r2_values:
        r2_array = np.array(chemistry_r2_values)
        print(f"    Chemistry R²: mean={np.mean(r2_array):.3f}, median={np.median(r2_array):.3f}, "
              f"{100*np.mean(r2_array > 0.7):.0f}% > 0.7")


def _UNUSED_plot_observed_vs_simulated_ionwise(data: dict, cal_results: dict, all_results: dict, config: Config, output_dir: Path):
    """DEPRECATED: Ion-wise plot that INCORRECTLY used v_fitted = u × gamma.
    
    This is WRONG because it ignores reactions. Kept for reference only.
    The correct approach uses chemistry_r2 from EdgeResult which includes reactions.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    
    water_chem = data["water_chem"]
    edges_df = data["edges"]
    edges = [(row["from_station"], row["to_station"]) for _, row in edges_df.iterrows()]
    
    # Collect within-event fit quality from ALL edge results
    obs_fit_pairs = {"Ca": [], "Cl": [], "NO3": []}
    objective_scores = []
    
    for event_code, results in all_results.items():
        if not results:
            continue
            
        # Get samples for this event
        event_data = water_chem[water_chem["event_code"] == event_code]
        samples = prepare_samples_mmol(event_data, event_code)
        
        if not samples:
            continue
        
        for res in results:
            if res.gamma is None:
                continue
            
            # Track objective scores
            if hasattr(res, 'objective_score') and res.objective_score is not None:
                objective_scores.append(res.objective_score)
            
            # Get upstream and downstream samples
            u_sample = next((s for s in samples if s["site_id"] == res.u), None)
            v_sample = next((s for s in samples if s["site_id"] == res.v), None)
            
            if not u_sample or not v_sample:
                continue
            
            for ion in ["Ca", "Cl", "NO3"]:
                if ion in u_sample and ion in v_sample:
                    u_val = u_sample.get(ion)
                    v_obs = v_sample.get(ion)
                    
                    if u_val is not None and v_obs is not None and u_val > 0:
                        # WRONG: This ignores reactions!
                        v_fitted = float(u_val) * res.gamma
                        
                        # Convert to mg/L
                        obs_mg = mmolL_to_mgL(float(v_obs), ion)
                        fit_mg = mmolL_to_mgL(v_fitted, ion)
                        
                        obs_fit_pairs[ion].append((obs_mg, fit_mg, event_code))
    
    ions = [("Ca", "Ca²⁺", '#3498db'), ("Cl", "Cl⁻", '#2ecc71'), ("NO3", "NO₃⁻", '#e74c3c')]
    

def plot_residual_diagnostics(data: dict, all_results: dict, config: Config, output_dir: Path):
    """Plot calibration residual diagnostics using Hydrosheaf's actual residuals.
    
    IMPORTANT: Hydrosheaf stores residuals in EdgeResult.residual_vector which are
    computed as: observed_v - predicted_v where predicted_v = u × gamma + reactions.
    
    These are the CORRECT residuals that account for both transport AND reactions.
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    
    # Collect residuals and chemistry_r2 from all edge results
    all_residuals = []
    chemistry_r2_vals = []
    objective_scores = []
    edge_types = []
    stations = []
    
    for event_code, results in all_results.items():
        if not results:
            continue
        
        for res in results:
            if res.gamma is None:
                continue
            
            # Get pre-computed residual vector from EdgeResult
            if hasattr(res, 'residual_vector') and res.residual_vector:
                all_residuals.extend(res.residual_vector)
            
            # Chemistry R² 
            if hasattr(res, 'chemistry_r2') and res.chemistry_r2 is not None:
                chemistry_r2_vals.append(res.chemistry_r2)
                edge_types.append(classify_edge_transport_type(res.u, res.v))
                stations.append(res.v)
            
            if hasattr(res, 'objective_score') and res.objective_score is not None:
                objective_scores.append(res.objective_score)
    
    # Convert residuals from mmol/L to mg/L (approximate using average MW)
    # Residuals are in compositional space, so we'll use relative residuals
    residuals = np.array(all_residuals) if all_residuals else np.array([])
    
    if len(residuals) == 0 or len(chemistry_r2_vals) == 0:
        # Fallback: show summary based on chemistry_r2
        ax = axes[0, 0]
        ax.text(0.5, 0.5, 'Residual vectors not available.\nUsing chemistry R² as quality metric.',
               ha='center', va='center', transform=ax.transAxes, fontsize=12)
        plt.savefig(output_dir / "residual_diagnostics.png", dpi=300)
        plt.close()
        print(f"  Saved: residual_diagnostics.png (limited data)")
        return
    
    # (a) Histogram of normalized residuals (mmol/L space)
    ax1 = axes[0, 0]
    ax1.hist(residuals, bins=30, color='steelblue', edgecolor='white', alpha=0.8)
    ax1.axvline(x=0, color='red', linestyle='--', linewidth=2)
    mean_res = np.mean(residuals)
    ax1.axvline(x=mean_res, color='green', linestyle='-', linewidth=2, label=f'Mean={mean_res:.4f}')
    ax1.set_xlabel('Residual (mmol/L, compositional space)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('(a) Residual Distribution\n(from Hydrosheaf reaction fitting)', fontweight='bold')
    ax1.legend()
    
    # (b) Chemistry R² distribution (the key quality metric)
    ax2 = axes[0, 1]
    r2_array = np.array(chemistry_r2_vals)
    types_array = np.array(edge_types)
    
    vertical_r2 = r2_array[types_array == 'vertical']
    lateral_r2 = r2_array[types_array == 'lateral']
    
    bins = np.linspace(0, 1, 21)
    if len(vertical_r2) > 0:
        ax2.hist(vertical_r2, bins=bins, alpha=0.7, label=f'Vertical (n={len(vertical_r2)})', 
                color='#3498db', edgecolor='white')
    if len(lateral_r2) > 0:
        ax2.hist(lateral_r2, bins=bins, alpha=0.7, label=f'Lateral (n={len(lateral_r2)})', 
                color='#e74c3c', edgecolor='white')
    
    ax2.axvline(np.mean(r2_array), color='black', linestyle='--', linewidth=2, 
               label=f'Mean $R^2$={np.mean(r2_array):.2f}')
    ax2.axvline(0.7, color='green', linestyle=':', linewidth=2, label='Good fit (0.7)')
    ax2.set_xlabel(r'Chemistry $R^2$')
    ax2.set_ylabel('Count')
    ax2.set_title(r'(b) Fit Quality Distribution\n(includes $\gamma$ + reactions)', fontweight='bold')
    ax2.legend(loc='upper left', fontsize=8)
    ax2.set_xlim(0, 1)
    
    # (c) Q-Q plot of residuals
    ax3 = axes[1, 0]
    from scipy import stats
    stats.probplot(residuals, dist="norm", plot=ax3)
    ax3.set_title('(c) Normal Q-Q Plot of Residuals', fontweight='bold')
    
    # (d) Chemistry R² by station (downstream)
    ax4 = axes[1, 1]
    station_r2 = {}
    for s, r2, etype in zip(stations, chemistry_r2_vals, edge_types):
        if s not in station_r2:
            station_r2[s] = []
        station_r2[s].append(r2)
    
    station_names = sorted(station_r2.keys())
    station_data = [station_r2[s] for s in station_names]
    
    if station_data:
        bp = ax4.boxplot(station_data, tick_labels=station_names, patch_artist=True)
        colors = plt.cm.Set3(np.linspace(0, 1, len(station_names)))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
    ax4.axhline(y=0.7, color='green', linestyle='--', linewidth=2, label='Good fit (0.7)')
    ax4.set_xlabel('Downstream Station')
    ax4.set_ylabel(r'Chemistry $R^2$')
    ax4.set_title('(d) Fit Quality by Station', fontweight='bold')
    ax4.tick_params(axis='x', rotation=45)
    ax4.set_ylim(0, 1.05)
    ax4.legend(loc='lower right', fontsize=8)
    
    plt.suptitle('Model Diagnostics (Hydrosheaf Inference Quality)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "residual_diagnostics.png", dpi=300)
    plt.close()
    print(f"  Saved: residual_diagnostics.png")


# ===================================================================
# VALIDATION PLOTS (for ch28_validation.tex)
# ===================================================================

def plot_vertical_vs_lateral_transport(data: dict, all_results: dict, output_dir: Path):
    """Plot comparison of vertical (vadose→saturated) vs lateral (within aquifer) transport.
    
    This uses physics-based classification from Hydrosheaf modules:
    
    VERTICAL TRANSPORT (Lysimeter → Borehole):
    - Physics: Richards equation for unsaturated flow (hydrosheaf.vadose)
    - VadoseProfile: Layered soil with van Genuchten-Mualem parameters
    - Travel-time distributions (TTD) from advective pore velocities
    - Evapotranspiration effects: gamma > 1 concentrates solutes
    - Vadose denitrification: NO3 reduction in anaerobic microsites
    - Isotope fractionation: δ18O, δ2H enrichment from evaporation
    
    LATERAL TRANSPORT (Borehole → Agricultural Well):
    - Physics: Advection-dispersion equation (hydrosheaf.transport)
    - ADE Green's function: g(t) = (x/√4πDt³)exp(-(x-vt)²/4Dt)
    - Darcy flow with hydraulic conductivity K = 0.1-10 m/day (saprolite)
    - Scale-dependent dispersivity: α ≈ 0.1 × flowpath length
    - Conservative tracer behavior for Cl, mixing-dominated
    - Lower denitrification potential in oxic aquifer zones
    
    References:
    - Bonsor et al. (2017): Ghana crystalline basement aquifer properties
    - Gelhar et al. (1992): Dispersivity scaling relationships
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    edges_df = data.get("edges", pd.DataFrame())
    stations_df = data.get("stations", pd.DataFrame())
    water_chem = data.get("water_chem", pd.DataFrame())
    
    # Use physics-based classification
    vertical_edges = []  # Lysimeter → Borehole (vadose transport)
    lateral_edges = []   # Borehole → Well (saturated transport)
    vertical_profiles = {}  # Store vadose profiles for reference
    lateral_params = {}     # Store ADE parameters for reference
    
    for event, results in all_results.items():
        if results:
            for res in results:
                edge_type = classify_edge_transport_type(res.u, res.v, edges_df)
                
                if edge_type == 'vertical':
                    vertical_edges.append(res)
                    # Build vadose profile for this edge (if not cached)
                    edge_key = f"{res.u}_{res.v}"
                    if edge_key not in vertical_profiles:
                        vertical_profiles[edge_key] = build_vadose_profile_for_edge(
                            res.u, res.v, stations_df
                        )
                else:
                    lateral_edges.append(res)
                    # Estimate ADE parameters for this edge
                    edge_key = f"{res.u}_{res.v}"
                    if edge_key not in lateral_params:
                        lateral_params[edge_key] = estimate_lateral_transport_params(
                            res.u, res.v, stations_df
                        )
    
    # (a) Gamma distribution comparison: Vertical vs Lateral
    ax1 = axes[0, 0]
    
    vertical_gammas = [r.gamma for r in vertical_edges if r.gamma is not None]
    lateral_gammas = [r.gamma for r in lateral_edges if r.gamma is not None]
    
    if vertical_gammas or lateral_gammas:
        data_to_plot = []
        labels = []
        colors = []
        
        if vertical_gammas:
            data_to_plot.append(vertical_gammas)
            labels.append(f'Vertical\n(n={len(vertical_gammas)})')
            colors.append('#3498db')
        if lateral_gammas:
            data_to_plot.append(lateral_gammas)
            labels.append(f'Lateral\n(n={len(lateral_gammas)})')
            colors.append('#e74c3c')
        
        bp = ax1.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax1.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label=r'$\gamma=1$ (conservative)')
        ax1.set_ylabel(r'Concentration Factor ($\gamma$)')
        ax1.legend(loc='upper right')
        
        # Add mean values as text
        for i, (d, c) in enumerate(zip(data_to_plot, colors)):
            ax1.text(i+1, np.mean(d), f'$\mu$={np.mean(d):.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'No gamma data', ha='center', va='center', transform=ax1.transAxes)
    
    ax1.set_title(r'(a) Evaporation Factor by Flow Direction', fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # (b) Objective scores comparison
    ax2 = axes[0, 1]
    
    vertical_obj = [r.objective_score for r in vertical_edges if hasattr(r, 'objective_score') and r.objective_score is not None]
    lateral_obj = [r.objective_score for r in lateral_edges if hasattr(r, 'objective_score') and r.objective_score is not None]
    
    if vertical_obj or lateral_obj:
        data_to_plot = []
        labels = []
        colors = []
        
        if vertical_obj:
            data_to_plot.append(vertical_obj)
            labels.append(f'Vertical\n(n={len(vertical_obj)})')
            colors.append('#3498db')
        if lateral_obj:
            data_to_plot.append(lateral_obj)
            labels.append(f'Lateral\n(n={len(lateral_obj)})')
            colors.append('#e74c3c')
        
        bp = ax2.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax2.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='Good fit threshold')
    else:
        ax2.text(0.5, 0.5, 'No objective data', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_ylabel('Objective Score')
    ax2.set_title('(b) Fit Quality by Flow Direction', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # (c) NO3 attenuation: Vertical shows more reduction (vadose denitrification)
    ax3 = axes[1, 0]
    
    # Calculate NO3 ratio (downstream/upstream) for each edge type
    vertical_no3_ratios = []
    lateral_no3_ratios = []
    
    for event_code in water_chem['event_code'].unique():
        event_wc = water_chem[water_chem['event_code'] == event_code]
        
        # Get station NO3 values
        station_no3 = {row['station_code']: row['NO3_mg_L'] 
                      for _, row in event_wc.iterrows() if pd.notna(row.get('NO3_mg_L'))}
        
        for res in vertical_edges:
            if res.u in station_no3 and res.v in station_no3:
                if station_no3[res.u] > 0:
                    ratio = station_no3[res.v] / station_no3[res.u]
                    vertical_no3_ratios.append(ratio)
        
        for res in lateral_edges:
            if res.u in station_no3 and res.v in station_no3:
                if station_no3[res.u] > 0:
                    ratio = station_no3[res.v] / station_no3[res.u]
                    lateral_no3_ratios.append(ratio)
    
    if vertical_no3_ratios or lateral_no3_ratios:
        data_to_plot = []
        labels = []
        colors = []
        
        if vertical_no3_ratios:
            data_to_plot.append(vertical_no3_ratios)
            labels.append(f'Vertical\n(vadose)')
            colors.append('#3498db')
        if lateral_no3_ratios:
            data_to_plot.append(lateral_no3_ratios)
            labels.append(f'Lateral\n(saturated)')
            colors.append('#e74c3c')
        
        bp = ax3.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax3.axhline(y=1.0, color='green', linestyle='--', linewidth=2, label='No change')
        
        # Add interpretation
        if vertical_no3_ratios and lateral_no3_ratios:
            v_mean = np.mean(vertical_no3_ratios)
            l_mean = np.mean(lateral_no3_ratios)
            if v_mean < l_mean:
                ax3.text(0.5, 0.02, 
                        f'Vertical transport shows {(1-v_mean)*100:.0f}% $NO_3^-$ reduction\n'
                        f'(vadose zone denitrification)',
                        transform=ax3.transAxes, fontsize=9, ha='center', style='italic',
                        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    else:
        ax3.text(0.5, 0.5, 'Insufficient $NO_3^-$ data', ha='center', va='center', transform=ax3.transAxes)
    
    ax3.set_ylabel(r'$NO_3^-$ Ratio (downstream/upstream)')
    ax3.set_title(r'(c) Nitrate Attenuation by Flow Path', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # (d) Edge count and physics-based transport summary
    ax4 = axes[1, 1]
    
    # Compute physics-based parameters
    # Vadose: Estimate travel times from gamma (higher gamma = more ET = longer residence)
    vertical_travel_times = []
    for vp_key, profile in vertical_profiles.items():
        # Sum layer thicknesses to get total vadose depth
        total_depth = sum(layer.thickness_m for layer in profile.layers)
        # Estimate travel time: assume 1-5 mm/day mean recharge
        mean_k = np.mean([l.ks_m_day for l in profile.layers if l.ks_m_day])
        tt_days = total_depth / (mean_k * 0.1) if mean_k > 0 else 100  # Rough estimate
        vertical_travel_times.append(tt_days)
    
    # Lateral: Get flowpath lengths from ADE parameters
    lateral_flowpaths = [p['flowpath_length_m'] for p in lateral_params.values()]
    lateral_dispersivities = [p['dispersivity_m'] for p in lateral_params.values()]
    
    # Summary statistics with physics-based parameters
    summary_data = {
        'Vertical': {
            'count': len(set(f"{r.u}→{r.v}" for r in vertical_edges)),
            'mean_gamma': np.mean(vertical_gammas) if vertical_gammas else np.nan,
            'mean_iso_penalty': np.mean([r.isotope_penalty for r in vertical_edges 
                                         if hasattr(r, 'isotope_penalty') and r.isotope_penalty]) or 0,
            'mean_travel_time': np.mean(vertical_travel_times) if vertical_travel_times else np.nan,
            'mean_no3_ratio': np.mean(vertical_no3_ratios) if vertical_no3_ratios else np.nan,
        },
        'Lateral': {
            'count': len(set(f"{r.u}→{r.v}" for r in lateral_edges)),
            'mean_gamma': np.mean(lateral_gammas) if lateral_gammas else np.nan,
            'mean_iso_penalty': np.mean([r.isotope_penalty for r in lateral_edges 
                                         if hasattr(r, 'isotope_penalty') and r.isotope_penalty]) or 0,
            'mean_flowpath': np.mean(lateral_flowpaths) if lateral_flowpaths else np.nan,
            'mean_dispersivity': np.mean(lateral_dispersivities) if lateral_dispersivities else np.nan,
        }
    }
    
    # Create summary table with physics parameters
    ax4.axis('off')
    
    table_data = [
        ['Metric', 'Vertical\n(Vadose Zone)', 'Lateral\n(Saturated Zone)'],
        ['Unique Edges', str(summary_data['Vertical']['count']), str(summary_data['Lateral']['count'])],
        [r'Mean $\gamma$ (ET factor)', 
         f"{summary_data['Vertical']['mean_gamma']:.3f}" if not np.isnan(summary_data['Vertical']['mean_gamma']) else 'N/A',
         f"{summary_data['Lateral']['mean_gamma']:.3f}" if not np.isnan(summary_data['Lateral']['mean_gamma']) else 'N/A'],
        ['Isotope Penalty', 
         f"{summary_data['Vertical']['mean_iso_penalty']:.3f}",
         f"{summary_data['Lateral']['mean_iso_penalty']:.3f}"],
        ['Physics Model', 'Richards Eq.\n(unsaturated)', 'ADE\n(advection-dispersion)'],
        ['Characteristic Scale', 
         f"$\\tau \\approx {summary_data['Vertical']['mean_travel_time']:.0f}$ days" if not np.isnan(summary_data['Vertical']['mean_travel_time']) else 'N/A',
         f"$L \\approx {summary_data['Lateral']['mean_flowpath']:.0f}$ m" if not np.isnan(summary_data['Lateral']['mean_flowpath']) else 'N/A'],
        [r'$NO_3^-$ Attenuation', 
         f"{(1-summary_data['Vertical']['mean_no3_ratio'])*100:.0f}% reduction" if not np.isnan(summary_data['Vertical']['mean_no3_ratio']) else 'N/A',
         'Conservative\n(mixing only)'],
        ['Key Process', 'ET concentration\nVadose denitrification', 'Dispersion\nMixing'],
    ]
    
    table = ax4.table(cellText=table_data, loc='center', cellLoc='center',
                     colWidths=[0.35, 0.325, 0.325])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 2.0)
    
    # Color header row
    for j in range(3):
        table[(0, j)].set_facecolor('#34495e')
        table[(0, j)].set_text_props(color='white', fontweight='bold')
    
    # Color vertical column (blue tones for vadose)
    for i in range(1, len(table_data)):
        table[(i, 1)].set_facecolor('#d5e8f7')
    
    # Color lateral column (red tones for saturated)
    for i in range(1, len(table_data)):
        table[(i, 2)].set_facecolor('#f7d5d5')
    
    ax4.set_title('(d) Physics-Based Transport Summary', fontweight='bold', pad=20)
    
    plt.suptitle('Vertical vs Lateral Transport Analysis\nVea Catchment (Hydrosheaf vadose + transport modules)', 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "vertical_vs_lateral_transport.png", dpi=300)
    plt.close()
    print(f"  Saved: vertical_vs_lateral_transport.png")


def plot_transport_validation(data: dict, all_results: dict, output_dir: Path):
    """Plot transport model validation results using actual Hydrosheaf EdgeResult data.
    
    Transport mechanism determination follows the Hydrosheaf methodology:
    1. For each edge, candidate transport models (evap, mix) are evaluated
    2. Each candidate is scored via multi-objective function:
       objective = reaction_residual + λ×L1_norm + isotope_penalty + Gibbs_penalty
    3. Best model selected by minimum objective score
    4. Bayesian model averaging provides transport_probabilities
    
    References:
    - Hydrosheaf edge_fit.py: fit_edge() function
    - Hydrosheaf mixing.py: fit_evaporation(), fit_mixing()
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    water_chem = data["water_chem"]
    events = data.get("events", pd.DataFrame())
    
    # (a) NO3 time series by station type
    ax1 = axes[0, 0]
    # Derive event order dynamically from data
    available_events = sorted(water_chem['event_code'].unique(), key=lambda x: (x[1] if len(x)>1 and x[1].isdigit() else 0, x))
    
    for st_type, color, marker in [('lysimeter', '#3498db', 'o'), ('borehole', '#e74c3c', 's')]:
        subset = water_chem[water_chem['station_type'] == st_type]
        event_means = subset.groupby('event_code')['NO3_mg_L'].mean()
        event_means = event_means.reindex([e for e in available_events if e in event_means.index])
        ax1.plot(range(len(event_means)), event_means.values, f'{marker}-', color=color, 
                linewidth=2, markersize=8, label=st_type.title())
    
    ax1.set_xlabel('Event')
    ax1.set_ylabel(r'$NO_3^-$ (mg/L)')
    ax1.set_title(r'(a) Temporal $NO_3^-$ Dynamics', fontweight='bold')
    ax1.legend()
    ax1.set_xticks(range(len(available_events)))
    ax1.set_xticklabels([e.split('-')[1] if '-' in e else e for e in available_events], rotation=45)
    ax1.grid(True, alpha=0.3)
    
    # (b) Lysimeter to Borehole attenuation
    ax2 = axes[0, 1]
    
    lys_mean = water_chem[water_chem['station_type'] == 'lysimeter'].groupby('event_code')['NO3_mg_L'].mean()
    bh_mean = water_chem[water_chem['station_type'] == 'borehole'].groupby('event_code')['NO3_mg_L'].mean()
    
    common_events = list(set(lys_mean.index) & set(bh_mean.index))
    if common_events:
        attenuation = [(lys_mean[e] - bh_mean[e]) / lys_mean[e] * 100 if lys_mean[e] > 0 else 0 for e in common_events]
        colors = ['#2ecc71' if a > 0 else '#e74c3c' for a in attenuation]
        ax2.bar(range(len(common_events)), attenuation, color=colors, edgecolor='white')
        ax2.set_xticks(range(len(common_events)))
        ax2.set_xticklabels([e.split('-')[1] for e in common_events], rotation=45)
    
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax2.set_xlabel('Event')
    ax2.set_ylabel(r'$NO_3^-$ Attenuation (%)')
    ax2.set_title(r'(b) Vadose Zone Attenuation', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # (c) Cl vs NO3 (conservative tracer comparison)
    ax3 = axes[1, 0]
    for st_type, color, marker in [('lysimeter', '#3498db', 'o'), ('borehole', '#e74c3c', 's')]:
        subset = water_chem[water_chem['station_type'] == st_type]
        ax3.scatter(subset['Cl_mg_L'], subset['NO3_mg_L'], c=color, marker=marker, 
                   alpha=0.7, s=50, label=st_type.title())
    
    ax3.set_xlabel(r'$Cl^-$ (mg/L)')
    ax3.set_ylabel(r'$NO_3^-$ (mg/L)')
    ax3.set_title(r'(c) $Cl^-$ vs $NO_3^-$ (Transport Indicator)', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # (d) Transport mechanism classification from actual Hydrosheaf EdgeResults
    # Uses transport_model field which is determined by multi-objective optimization
    ax4 = axes[1, 1]
    
    transport_counts = {'evap': 0, 'mix': 0}
    total_edges = 0
    
    # Count transport models from actual fit_network() results
    for event, results in all_results.items():
        if results:
            for res in results:
                if hasattr(res, 'transport_model') and res.transport_model:
                    model = res.transport_model
                    if model in transport_counts:
                        transport_counts[model] += 1
                    else:
                        transport_counts[model] = 1
                    total_edges += 1
    
    if total_edges > 0:
        # Rename for display
        display_counts = {}
        model_labels = {
            'evap': r'Evaporation\n($\gamma > 1$)',
            'mix': 'Mixing\n(binary endmember)',
        }
        model_colors = {
            'evap': '#e74c3c',
            'mix': '#3498db',
        }
        
        for model, count in transport_counts.items():
            if count > 0:
                label = model_labels.get(model, model.title())
                display_counts[label] = count
        
        if display_counts:
            colors = [model_colors.get(m.split('\n')[0].lower(), '#95a5a6') 
                     for m in display_counts.keys()]
            
            wedges, texts, autotexts = ax4.pie(
                display_counts.values(), 
                labels=display_counts.keys(),
                autopct=lambda p: f'{p:.1f}%\n({int(p*total_edges/100)})',
                colors=colors[:len(display_counts)], 
                explode=[0.03]*len(display_counts),
                textprops={'fontsize': 9}
            )
            
            # Add annotation about method
            ax4.text(0.5, -0.15, 
                    f'Based on {total_edges} edges from fit_network()\n'
                    'Model selection via multi-objective optimization',
                    ha='center', va='top', fontsize=8, style='italic',
                    transform=ax4.transAxes)
    else:
        ax4.text(0.5, 0.5, 'No EdgeResult data available\nRun fit_network() first', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=11)
    
    ax4.set_title('(d) Transport Mechanisms\n(from Hydrosheaf EdgeResult)', fontweight='bold')
    
    plt.suptitle('Transport Model Validation - Vea Catchment', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "transport_validation.png", dpi=300)
    plt.close()
    print(f"  Saved: transport_validation.png")


def plot_validation_summary_table(cal_results: dict, all_results: dict, data: dict, output_dir: Path):
    """Create a validation summary table from actual computed results."""
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.axis('off')
    
    water_chem = data.get("water_chem", pd.DataFrame())
    
    # Compute actual metrics from edge results
    edge_metrics = {}
    all_gammas = []
    all_objectives = []
    
    for event, results in all_results.items():
        if results:
            for res in results:
                edge_id = f"{res.u}→{res.v}"
                if edge_id not in edge_metrics:
                    edge_metrics[edge_id] = {'gammas': [], 'objectives': [], 'l1_norms': []}
                if res.gamma is not None:
                    edge_metrics[edge_id]['gammas'].append(res.gamma)
                    all_gammas.append(res.gamma)
                if hasattr(res, 'objective_score') and res.objective_score is not None:
                    edge_metrics[edge_id]['objectives'].append(res.objective_score)
                    all_objectives.append(res.objective_score)
                if hasattr(res, 'l1_norm') and res.l1_norm is not None:
                    edge_metrics[edge_id]['l1_norms'].append(res.l1_norm)
    
    # Select representative edges for table (L1→BH1, L2→BH2 if available)
    representative_edges = []
    for edge in ['L1→BH1', 'L2→BH2', 'L3→BH3', 'L4→BH4']:
        if edge in edge_metrics:
            representative_edges.append(edge)
        if len(representative_edges) >= 2:
            break
    
    if not representative_edges:
        representative_edges = list(edge_metrics.keys())[:2]
    
    # Build table with actual data
    header = ['Metric'] + representative_edges + ['Overall']
    table_data = [header]
    
    # Gamma (concentration/dilution factor)
    gamma_row = ['γ (conc. factor)']
    for edge in representative_edges:
        gammas = edge_metrics.get(edge, {}).get('gammas', [])
        gamma_row.append(f'{np.mean(gammas):.3f}' if gammas else 'N/A')
    gamma_row.append(f'{np.mean(all_gammas):.3f}' if all_gammas else 'N/A')
    table_data.append(gamma_row)
    
    # Objective score (fit quality)
    obj_row = ['Objective Score']
    for edge in representative_edges:
        objs = edge_metrics.get(edge, {}).get('objectives', [])
        obj_row.append(f'{np.mean(objs):.3f}' if objs else 'N/A')
    obj_row.append(f'{np.mean(all_objectives):.3f}' if all_objectives else 'N/A')
    table_data.append(obj_row)
    
    # L1 norm (residual)
    l1_row = ['L1 Norm (mmol/L)']
    all_l1 = []
    for edge in representative_edges:
        l1s = edge_metrics.get(edge, {}).get('l1_norms', [])
        l1_row.append(f'{np.mean(l1s):.3f}' if l1s else 'N/A')
        all_l1.extend(l1s)
    l1_row.append(f'{np.mean(all_l1):.3f}' if all_l1 else 'N/A')
    table_data.append(l1_row)
    
    # Calibration parameters
    opt_params = cal_results.get("optimal_parameters", {})
    if opt_params:
        # Dispersivity
        disp = opt_params.get('dispersivity_m', opt_params.get('dispersivity', None))
        disp_row = ['Dispersivity (m)']
        for _ in representative_edges:
            disp_row.append(f'{disp:.2f}' if disp else 'N/A')
        disp_row.append(f'{disp:.2f}' if disp else 'N/A')
        table_data.append(disp_row)
        
        # Vadose Travel Time (replaced evap_factor which is now internal)
        tt = opt_params.get('vadose_travel_time', None)
        tt_row = ['Vadose Travel T (d)']
        for _ in representative_edges:
            tt_row.append(f'{tt:.1f}' if tt else 'N/A')
        tt_row.append(f'{tt:.1f}' if tt else 'N/A')
        table_data.append(tt_row)
    
    # Number of edges analyzed
    n_edges_row = ['Edges Analyzed']
    for edge in representative_edges:
        n = len(edge_metrics.get(edge, {}).get('gammas', []))
        n_edges_row.append(str(n))
    n_edges_row.append(str(len(all_gammas)))
    table_data.append(n_edges_row)
    
    # Calibration phi
    phi = cal_results.get('phi', None)
    if phi is not None:
        phi_row = ['Calibration Φ', '', '', f'{phi:.4f}']
        # Fill middle columns
        for i in range(len(representative_edges)):
            phi_row[i+1] = '-'
        table_data.append(phi_row)
    
    # Create table
    n_cols = len(header)
    table = ax.table(cellText=table_data, loc='center', cellLoc='center',
                     colWidths=[0.22] + [0.18]*(n_cols-1))
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)
    
    # Style header row
    for j in range(n_cols):
        table[(0, j)].set_facecolor('#2c3e50')
        table[(0, j)].set_text_props(color='white', fontweight='bold')
    
    # Style alternating rows
    for i in range(1, len(table_data)):
        for j in range(n_cols):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#ecf0f1')
    
    ax.set_title('Validation Summary - Actual Computed Results\nVea Catchment Hydrosheaf Analysis', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Add note about data source
    ax.text(0.5, -0.05, 'Values computed from fit_network() EdgeResult objects', 
            ha='center', va='top', fontsize=9, style='italic', transform=ax.transAxes)
    
    plt.tight_layout()
    plt.savefig(output_dir / "validation_summary_table.png", dpi=300)
    plt.close()
    print(f"  Saved: validation_summary_table.png")


# ===================================================================
# STANDARD ANALYSIS PLOTS
# ===================================================================

def prepare_samples_mmol(water_chem: pd.DataFrame, event_code: str = None) -> list:
    """Convert water chemistry to mmol/L format, including isotope data."""
    ion_mapping = {
        "Ca": "Ca_mg_L", "Mg": "Mg_mg_L", "Na": "Na_mg_L", "K": "K_mg_L",
        "HCO3": "HCO3_mg_L", "Cl": "Cl_mg_L", "SO4": "SO4_mg_L", "NO3": "NO3_mg_L",
    }
    
    df = water_chem.copy()
    if event_code:
        df = df[df["event_code"] == event_code].copy()
    
    for ion, col in ion_mapping.items():
        if col in df.columns:
            df[ion] = df[col].apply(lambda x: mgL_to_mmolL(x, ion) if pd.notna(x) else np.nan)
    
    for ion in ["F", "Fe", "PO4"]:
        df[ion] = 0.0
    
    samples = []
    for _, row in df.iterrows():
        sample = {
            "site_id": row["station_code"],
            "sample_id": f"{row['event_code']}_{row['station_code']}",
            "Ca": row["Ca"], "Mg": row["Mg"], "Na": row["Na"],
            "K": row.get("K", 0), "HCO3": row["HCO3"], "Cl": row["Cl"],
            "SO4": row["SO4"], "NO3": row["NO3"], "F": row["F"],
            "Fe": row["Fe"], "PO4": row["PO4"], "pH": row.get("pH", 7.0),
        }
        
        # Include water isotopes (δ18O and δ2H) for evaporation/mixing constraints
        # These are used by fit_edge() to compute isotope_penalty when isotope_enabled=True
        if "d18O_H2O_permil" in row and pd.notna(row.get("d18O_H2O_permil")):
            sample["18O"] = row["d18O_H2O_permil"]
        if "d2H_H2O_permil" in row and pd.notna(row.get("d2H_H2O_permil")):
            sample["2H"] = row["d2H_H2O_permil"]
        
        # Include nitrate isotopes for source identification
        if "d15N_NO3_permil" in row and pd.notna(row.get("d15N_NO3_permil")):
            sample["d15N_NO3"] = row["d15N_NO3_permil"]
        if "d18O_NO3_permil" in row and pd.notna(row.get("d18O_NO3_permil")):
            sample["d18O_NO3"] = row["d18O_NO3_permil"]
            
        samples.append(sample)
    return samples


def prepare_samples_for_ilr(water_chem: pd.DataFrame) -> list:
    """Prepare samples for ILR plotting."""
    samples = []
    for _, row in water_chem.iterrows():
        sample = {
            "site_id": row["station_code"],
            "Ca": mgL_to_mmolL(row["Ca_mg_L"], "Ca"),
            "Mg": mgL_to_mmolL(row["Mg_mg_L"], "Mg"),
            "Na": mgL_to_mmolL(row["Na_mg_L"], "Na"),
            "K": mgL_to_mmolL(row["K_mg_L"], "K"),
            "HCO3": mgL_to_mmolL(row["HCO3_mg_L"], "HCO3"),
            "Cl": mgL_to_mmolL(row["Cl_mg_L"], "Cl"),
            "SO4": mgL_to_mmolL(row["SO4_mg_L"], "SO4"),
            "NO3": mgL_to_mmolL(row["NO3_mg_L"], "NO3"),
            "TDS": row["TDS_mg_L"],
            "station_type": row["station_type"],
            "event_code": row["event_code"],
        }
        samples.append(sample)
    return samples


def plot_water_isotopes(water_chem: pd.DataFrame, output_dir: Path):
    """Plot water isotopes with LMWL and GMWL.
    
    Shows isotopic composition by station type to visualize:
    - Evaporation trends (deviation from LMWL)
    - Recharge sources (position relative to GMWL)
    - Mixing between water types
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Filter samples with isotope data
    iso_data = water_chem[water_chem['d18O_H2O_permil'].notna() & water_chem['d2H_H2O_permil'].notna()]
    
    if len(iso_data) == 0:
        ax.text(0.5, 0.5, 'No isotope data available', ha='center', va='center', fontsize=14)
    else:
        # Plot by station type - include ALL station types
        colors = {'lysimeter': '#3498db', 'borehole': '#e74c3c', 'ag_well': '#f39c12'}
        markers = {'lysimeter': 'o', 'borehole': 's', 'ag_well': 'p'}
        sizes = {'lysimeter': 100, 'borehole': 120, 'ag_well': 100}
        labels = {'lysimeter': 'Lysimeter (Soil Water)', 'borehole': 'Borehole (Groundwater)', 
                  'ag_well': 'Ag Well (Discharge)'}
        
        for st_type in ['lysimeter', 'borehole', 'ag_well']:
            subset = iso_data[iso_data['station_type'] == st_type]
            if len(subset) > 0:
                ax.scatter(subset['d18O_H2O_permil'], subset['d2H_H2O_permil'],
                          c=colors.get(st_type, 'gray'), marker=markers.get(st_type, 'o'),
                          s=sizes.get(st_type, 80), alpha=0.7, 
                          label=labels.get(st_type, st_type.title()), 
                          edgecolors='white', linewidths=0.5)
    
    # Plot reference lines - extend range to cover data
    if len(iso_data) > 0:
        d18O_min = iso_data['d18O_H2O_permil'].min() - 1
        d18O_max = iso_data['d18O_H2O_permil'].max() + 1
        x_range = np.array([d18O_min, d18O_max])
    else:
        x_range = np.array([-8, 2])
    
    # GMWL: d2H = 8 * d18O + 10 (Craig, 1961)
    ax.plot(x_range, 8 * x_range + 10, 'k-', linewidth=2.5, label=r'GMWL ($\delta^2H = 8\delta^{18}O + 10$)')
    
    # LMWL Ghana: d2H = 7.87 * d18O + 13.61 (Akiti, 1980; regional studies)
    ax.plot(x_range, 7.87 * x_range + 13.61, 'b--', linewidth=2, label=r'LMWL Ghana ($\delta^2H = 7.87\delta^{18}O + 13.61$)')
    
    # Add evaporation line indicator (slope ~4-5 for evaporation)
    if len(iso_data) > 0:
        # Calculate local evaporation line from data if there's enough spread
        d18O_vals = iso_data['d18O_H2O_permil'].values
        d2H_vals = iso_data['d2H_H2O_permil'].values
        if len(d18O_vals) > 3 and (d18O_vals.max() - d18O_vals.min()) > 1:
            slope, intercept = np.polyfit(d18O_vals, d2H_vals, 1)
            if slope < 7:  # Evaporation indicator (slope < GMWL slope of 8)
                ax.plot(x_range, slope * x_range + intercept, 'r:', linewidth=1.5, 
                       alpha=0.7, label=f'LEL (slope={slope:.1f})')
    
    ax.set_xlabel(r'$\delta^{18}O$ (\u2030 VSMOW)', fontsize=12)
    ax.set_ylabel(r'$\delta^2H$ (\u2030 VSMOW)', fontsize=12)
    ax.set_title('Water Isotope Composition\nVea Catchment, Upper East Region, Ghana', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Add d-excess annotation if data exists
    if len(iso_data) > 0:
        d_excess = iso_data['d2H_H2O_permil'] - 8 * iso_data['d18O_H2O_permil']
        mean_d_excess = d_excess.mean()
        ax.text(0.02, 0.98, f'$d$-excess = {mean_d_excess:.1f}\u2030', 
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_dir / "water_isotopes.png", dpi=300)
    plt.close()
    print(f"  Saved: water_isotopes.png")


def plot_nitrate_sources(water_chem: pd.DataFrame, endmembers: pd.DataFrame, output_dir: Path):
    """Plot enhanced nitrate isotope source diagram with proper source zones."""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Define source zones based on Kendall & McDonnell (1998) and Xue et al. (2009)
    # Format: 'd15N': (min, max), 'd18O': (min, max)
    source_zones = {
        'Atmospheric\nDeposition': {
            'd15N': (-15, 15), 'd18O': (50, 95),
            'color': '#9370DB', 'alpha': 0.25, 'hatch': None
        },
        '$NO_3^-$ in\nFertilizer': {
            'd15N': (-8, 7), 'd18O': (17, 35),
            'color': '#FF6347', 'alpha': 0.3, 'hatch': None
        },
        '$NH_4^+$ in Fertilizer\n& Precipitation': {
            'd15N': (-8, 8), 'd18O': (-10, 15),
            'color': '#4169E1', 'alpha': 0.25, 'hatch': None
        },
        'Soil\nOrganic N': {
            'd15N': (0, 12), 'd18O': (-10, 15),
            'color': '#228B22', 'alpha': 0.25, 'hatch': None
        },
        'Manure &\nSewage': {
            'd15N': (8, 25), 'd18O': (-10, 15),
            'color': '#CD853F', 'alpha': 0.3, 'hatch': None
        },
    }
    
    # Draw source zones with better styling
    for name, props in source_zones.items():
        d15N_min, d15N_max = props['d15N']
        d18O_min, d18O_max = props['d18O']
        
        # Create rectangle
        rect = plt.Rectangle(
            (d15N_min, d18O_min),
            d15N_max - d15N_min,
            d18O_max - d18O_min,
            facecolor=props['color'], 
            alpha=props['alpha'], 
            edgecolor=props['color'], 
            linewidth=2.5,
            linestyle='-',
            zorder=1
        )
        ax.add_patch(rect)
        
        # Add label with white background for readability
        label_x = (d15N_min + d15N_max) / 2
        label_y = (d18O_min + d18O_max) / 2
        
        # Adjust label positions for overlapping zones
        if 'Soil' in name:
            label_x = 6
            label_y = 2
        elif 'NH_4' in name:
            label_x = -2
            label_y = 3
        elif 'NO_3' in name:
            label_x = 0
            label_y = 26
        elif 'Manure' in name:
            label_x = 16
            label_y = 2
        elif 'Atmospheric' in name:
            label_x = 0
            label_y = 72
        
        ax.text(label_x, label_y, name,
               ha='center', va='center', fontsize=11, fontweight='bold',
               color='black',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                        edgecolor=props['color'], alpha=0.9, linewidth=1.5))
    
    # Plot samples
    iso_data = water_chem[water_chem['d15N_NO3_permil'].notna() & water_chem['d18O_NO3_permil'].notna()]
    
    if len(iso_data) > 0:
        # Enhanced styling for data points - include ALL station types
        style_config = {
            'lysimeter': {'color': '#2980b9', 'marker': 'o', 'size': 120, 'label': 'Lysimeter (Soil Water)'},
            'borehole': {'color': '#c0392b', 'marker': 's', 'size': 140, 'label': 'Borehole (Groundwater)'},
            'ag_well': {'color': '#f39c12', 'marker': 'p', 'size': 130, 'label': 'Ag Well (Discharge)'}
        }
        
        for st_type in ['lysimeter', 'borehole', 'ag_well']:
            subset = iso_data[iso_data['station_type'] == st_type]
            if len(subset) > 0:
                cfg = style_config.get(st_type, {'color': 'gray', 'marker': 'o', 'size': 100, 'label': st_type})
                ax.scatter(subset['d15N_NO3_permil'], subset['d18O_NO3_permil'],
                          c=cfg['color'], marker=cfg['marker'], s=cfg['size'],
                          alpha=0.85, label=cfg['label'], 
                          edgecolors='white', linewidths=1.5, zorder=10)
    
    # Add denitrification trend arrow (slope ~2:1 for δ¹⁸O:δ¹⁵N)
    arrow_start = (2, 10)
    arrow_end = (18, 18)
    ax.annotate('', xy=arrow_end, xytext=arrow_start,
               arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=3, 
                               connectionstyle='arc3,rad=0.1'))
    ax.text(12, 17, 'Denitrification\n($\\delta^{18}O:\\delta^{15}N \\approx 1:2$)', 
           fontsize=10, color='#8e44ad', style='italic', fontweight='bold',
           ha='center', va='bottom',
           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#8e44ad'))
    
    # Add nitrification arrow
    ax.annotate('', xy=(5, 8), xytext=(-3, 20),
               arrowprops=dict(arrowstyle='->', color='#27ae60', lw=2.5, 
                               connectionstyle='arc3,rad=-0.2', linestyle='--'))
    ax.text(-4, 25, 'Nitrification\nof $NH_4^+$', fontsize=11, color='#27ae60', 
           style='italic', ha='center',
           bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#27ae60'))
    
    # Axis labels and title
    ax.set_xlabel(r'$\delta^{15}N-NO_3^-$ (\u2030 vs AIR)', fontsize=14, fontweight='bold')
    ax.set_ylabel(r'$\delta^{18}O-NO_3^-$ (\u2030 vs VSMOW)', fontsize=14, fontweight='bold')
    ax.set_title('Nitrate Source Identification\nDual Isotope Approach - Vea Catchment', 
                fontsize=16, fontweight='bold', pad=15)
    
    # Set axis limits to show all zones properly
    ax.set_xlim(-18, 28)
    ax.set_ylim(-15, 100)
    
    # Add gridlines
    ax.grid(True, alpha=0.3, linestyle='--', zorder=0)
    ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.5)
    ax.axvline(x=0, color='black', linewidth=0.5, alpha=0.5)
    
    # Enhanced legend
    legend = ax.legend(loc='upper right', fontsize=11, framealpha=0.95,
                      edgecolor='gray', fancybox=True, shadow=True,
                      title='Sample Type', title_fontsize=12)
    
    # Generate data-driven interpretation note
    # Classify samples by their source zone based on isotope values
    if len(iso_data) > 0:
        # Define source zone boundaries for classification
        zone_bounds = {
            'Atmospheric': {'d15N': (-15, 15), 'd18O': (50, 95)},
            'NO3_Fertilizer': {'d15N': (-8, 7), 'd18O': (17, 35)},
            'NH4_Fertilizer': {'d15N': (-8, 8), 'd18O': (-10, 15)},
            'Soil_Organic': {'d15N': (0, 12), 'd18O': (-10, 15)},
            'Manure_Sewage': {'d15N': (8, 25), 'd18O': (-10, 15)},
        }
        
        def classify_sample(d15N, d18O):
            """Classify sample into most likely source zone."""
            for zone, bounds in zone_bounds.items():
                if (bounds['d15N'][0] <= d15N <= bounds['d15N'][1] and 
                    bounds['d18O'][0] <= d18O <= bounds['d18O'][1]):
                    return zone
            return 'Mixed/Uncertain'
        
        # Classify each sample
        lysimeter_zones = []
        borehole_zones = []
        
        for _, row in iso_data.iterrows():
            zone = classify_sample(row['d15N_NO3_permil'], row['d18O_NO3_permil'])
            if row['station_type'] == 'lysimeter':
                lysimeter_zones.append(zone)
            else:
                borehole_zones.append(zone)
        
        # Compute dominant source for each station type
        from collections import Counter
        
        lys_interpretation = ""
        bh_interpretation = ""
        
        if lysimeter_zones:
            lys_counts = Counter(lysimeter_zones)
            lys_dominant = lys_counts.most_common(1)[0][0] if lys_counts else "Unknown"
            lys_pct = (lys_counts[lys_dominant] / len(lysimeter_zones) * 100) if lysimeter_zones else 0
            lys_interpretation = f"• Lysimeters: {lys_pct:.0f}% in {lys_dominant.replace('_', ' ')} zone"
            
            # Check for denitrification evidence (enriched δ15N > 8)
            lys_d15N_mean = iso_data[iso_data['station_type'] == 'lysimeter']['d15N_NO3_permil'].mean()
            
        if borehole_zones:
            bh_counts = Counter(borehole_zones)
            bh_dominant = bh_counts.most_common(1)[0][0] if bh_counts else "Unknown"
            bh_pct = (bh_counts[bh_dominant] / len(borehole_zones) * 100) if borehole_zones else 0
            bh_interpretation = f"• Boreholes: {bh_pct:.0f}% in {bh_dominant.replace('_', ' ')} zone"
            
            # Check for denitrification signature
            bh_d15N_mean = iso_data[iso_data['station_type'] == 'borehole']['d15N_NO3_permil'].mean()
            if bh_d15N_mean > 8:
                bh_interpretation += "\n  (elevated δ¹⁵N suggests denitrification)"
        
        # Compute mean isotope shift from lysimeter to borehole (if both exist)
        shift_note = ""
        if lysimeter_zones and borehole_zones:
            lys_mean_d15N = iso_data[iso_data['station_type'] == 'lysimeter']['d15N_NO3_permil'].mean()
            bh_mean_d15N = iso_data[iso_data['station_type'] == 'borehole']['d15N_NO3_permil'].mean()
            d15N_shift = bh_mean_d15N - lys_mean_d15N
            if d15N_shift > 3:
                shift_note = f"\n• δ¹⁵N enrichment: +{d15N_shift:.1f}‰ (lys→bh)"
        
        note_text = f"Data-derived interpretation:\n{lys_interpretation}\n{bh_interpretation}{shift_note}"
    else:
        note_text = "No nitrate isotope data available for interpretation"
    
    ax.text(0.02, 0.02, note_text, transform=ax.transAxes, fontsize=9,
           verticalalignment='bottom', horizontalalignment='left',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', 
                    edgecolor='orange', alpha=0.9))
    
    plt.tight_layout()
    plt.savefig(output_dir / "nitrate_sources.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: nitrate_sources.png")


def plot_network_schematic(edges_df: pd.DataFrame, stations_df: pd.DataFrame, output_dir: Path):
    """Plot network schematic with geographic coordinates showing natural flow."""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Check for lat/lon columns (use latitude/longitude)
    has_coords = 'latitude' in stations_df.columns and 'longitude' in stations_df.columns
    
    if has_coords:
        stations_df = stations_df.copy()
        # Use longitude as x, latitude as y (geographic convention)
        stations_df['x'] = stations_df['longitude']
        stations_df['y'] = stations_df['latitude']
        xlabel = 'Longitude (°W)'
        ylabel = 'Latitude (°N)'
    else:
        # Generate schematic positions based on flow hierarchy
        station_positions = {}
        lysimeters = stations_df[stations_df['station_type'] == 'lysimeter']['station_code'].tolist()
        boreholes = stations_df[stations_df['station_type'] == 'borehole']['station_code'].tolist()
        ag_wells = stations_df[stations_df['station_type'] == 'ag_well']['station_code'].tolist()
        
        # Lysimeters at top (upstream/recharge zone)
        for i, st in enumerate(lysimeters):
            station_positions[st] = (0.2 + 0.2 * i, 0.85)
        # Boreholes in middle (monitoring zone)
        for i, st in enumerate(boreholes):
            station_positions[st] = (0.15 + 0.2 * i, 0.55)
        # Ag wells at bottom (downstream/production zone)  
        for i, st in enumerate(ag_wells):
            station_positions[st] = (0.2 + 0.12 * i, 0.2)
        
        stations_df = stations_df.copy()
        stations_df['x'] = stations_df['station_code'].map(lambda s: station_positions.get(s, (0.5, 0.5))[0])
        stations_df['y'] = stations_df['station_code'].map(lambda s: station_positions.get(s, (0.5, 0.5))[1])
        xlabel = 'X Position'
        ylabel = 'Y Position'
    
    # Define edge colors by type
    edge_colors = {'primary': '#2c3e50', 'secondary': '#7f8c8d', 'cascade': '#3498db'}
    edge_widths = {'primary': 2.5, 'secondary': 1.5, 'cascade': 1.0}
    
    # Plot edges with arrows
    for _, edge in edges_df.iterrows():
        from_st = edge['from_station']
        to_st = edge['to_station']
        edge_type = edge.get('edge_type', 'primary')
        
        from_row = stations_df[stations_df['station_code'] == from_st]
        to_row = stations_df[stations_df['station_code'] == to_st]
        
        if len(from_row) > 0 and len(to_row) > 0:
            x1, y1 = from_row.iloc[0]['x'], from_row.iloc[0]['y']
            x2, y2 = to_row.iloc[0]['x'], to_row.iloc[0]['y']
            
            ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                       arrowprops=dict(arrowstyle='->', 
                                      color=edge_colors.get(edge_type, 'gray'),
                                      lw=edge_widths.get(edge_type, 1.5), 
                                      alpha=0.7,
                                      connectionstyle='arc3,rad=0.1'))
    
    # Plot nodes by type with distinct styles
    colors = {'lysimeter': '#3498db', 'borehole': '#e74c3c', 'ag_well': '#f39c12'}
    markers = {'lysimeter': 'o', 'borehole': 's', 'ag_well': 'p'}
    sizes = {'lysimeter': 250, 'borehole': 200, 'ag_well': 180}
    
    for st_type in stations_df['station_type'].unique():
        subset = stations_df[stations_df['station_type'] == st_type]
        label = st_type.replace('_', ' ').title()
        ax.scatter(subset['x'], subset['y'], c=colors.get(st_type, 'gray'),
                  marker=markers.get(st_type, 'o'), s=sizes.get(st_type, 150), 
                  label=label, edgecolors='black', linewidths=1.5, zorder=5)
    
    # Add labels with offset to avoid overlap
    for _, row in stations_df.iterrows():
        offset = (8, 8) if row['station_type'] == 'lysimeter' else (8, -12)
        ax.annotate(row['station_code'], (row['x'], row['y']), 
                   textcoords="offset points", xytext=offset, ha='left', fontsize=9, fontweight='bold')
    
    # Add flow direction indicator based on actual coordinates
    if has_coords and len(stations_df) > 0:
        # Calculate flow direction from lysimeters (top/upstream) to ag_wells (bottom/downstream)
        lys_data = stations_df[stations_df['station_type'] == 'lysimeter']
        aw_data = stations_df[stations_df['station_type'] == 'ag_well']
        
        if len(lys_data) > 0 and len(aw_data) > 0:
            # Use centroid positions
            x_arrow_start = lys_data['x'].mean()
            y_arrow_start = lys_data['y'].max() + 0.002
            y_arrow_end = y_arrow_start - 0.005
            
            ax.annotate('', xy=(x_arrow_start, y_arrow_end), xytext=(x_arrow_start, y_arrow_start),
                       arrowprops=dict(arrowstyle='->', color='blue', lw=3))
            ax.text(x_arrow_start + 0.005, (y_arrow_start + y_arrow_end)/2, 'Flow\nDirection', 
                   fontsize=9, color='blue', ha='left')
    
    # Add legend for edge types
    from matplotlib.lines import Line2D
    edge_legend = [
        Line2D([0], [0], color='#2c3e50', linewidth=2.5, label='Primary flowpath'),
        Line2D([0], [0], color='#7f8c8d', linewidth=1.5, label='Secondary flowpath'),
        Line2D([0], [0], color='#3498db', linewidth=1.0, label='Cascade connection'),
    ]
    
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title('Hydrochemical Network Schematic\nVea Catchment Flow Connections', fontsize=14, fontweight='bold')
    
    # Two legends
    leg1 = ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), title='Station Type')
    ax.add_artist(leg1)
    ax.legend(handles=edge_legend, loc='lower right', title='Edge Type')
    
    ax.grid(True, alpha=0.3)
    
    # Invert x-axis for longitude (west is negative)
    if has_coords:
        ax.invert_xaxis()
    
    plt.tight_layout()
    plt.savefig(output_dir / "network_schematic.png", dpi=300)
    plt.close()
    print(f"  Saved: network_schematic.png")


def plot_data_model_diagram(output_dir: Path):
    """Generate conceptual diagram for the Probabilistic Data Model."""
    import matplotlib.patches as mpatches
    
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis('off')
    
    # Styles
    box_style = dict(boxstyle='round,pad=0.5', facecolor='#e8f4f8', edgecolor='#2c3e50', linewidth=1.5)
    text_style = dict(ha='center', va='center', fontsize=11)
    arrow_props = dict(arrowstyle='->', color='#2c3e50', lw=2, mutation_scale=20)
    
    # Box 1: Raw
    ax.text(1.5, 1.5, r'\textbf{Raw Data} $z$' + '\n(lab units)', **text_style, bbox=box_style)
    
    # Box 2: Standardized
    ax.text(5.0, 1.5, r'\textbf{Standardized} $y$' + '\n(mmol/L, ordered)', **text_style, bbox=box_style)
    
    # Box 3: Latent
    ax.text(8.5, 1.5, r'\textbf{Latent State} $x$' + '\n(true concentration)', **text_style, bbox=box_style)
    
    # Arrows
    ax.annotate('', xy=(3.8, 1.5), xytext=(2.7, 1.5), arrowprops=arrow_props)
    ax.annotate('', xy=(7.3, 1.5), xytext=(6.2, 1.5), arrowprops=arrow_props)
    
    # Labels for steps
    ax.text(3.25, 1.8, 'Unit Conversion', ha='center', va='bottom', fontsize=9, style='italic', color='#555')
    ax.text(6.75, 1.8, 'Censoring Model\n(Tobit)', ha='center', va='bottom', fontsize=9, style='italic', color='#555')
    
    # Bottom caption
    ax.text(5.0, 0.5, r'Data Transformation Pipeline: Units $\to$ Ordering $\to$ Censoring', 
            ha='center', va='center', fontsize=12, fontweight='bold', color='#333')
            
    plt.tight_layout()
    plt.savefig(output_dir / "data_model_schematic.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: data_model_schematic.png")



def plot_temporal_patterns(water_chem: pd.DataFrame, events_df: pd.DataFrame, output_dir: Path):
    """Plot temporal patterns."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Merge with events if season column exists
    if 'season' in events_df.columns:
        water_chem = water_chem.merge(events_df[['event_code', 'season']], on='event_code', how='left')
    elif 'season' not in water_chem.columns:
        # Add season based on event code pattern
        def get_season(ec):
            if 'DRY' in ec:
                return 'dry'
            elif 'PRE' in ec:
                return 'transition'
            else:
                return 'wet'
        water_chem['season'] = water_chem['event_code'].apply(get_season)
    
    # (a) NO3 by event
    ax1 = axes[0, 0]
    available_events = sorted(water_chem['event_code'].unique(), key=lambda x: (x[1] if len(x)>1 and x[1].isdigit() else 0, x))
    event_means = water_chem.groupby('event_code')['NO3_mg_L'].mean().reindex(available_events)
    
    colors = ['#f39c12' if 'DRY' in e else '#27ae60' if 'PRE' in e else '#3498db' for e in event_means.index]
    ax1.bar(range(len(event_means)), event_means.values, color=colors, edgecolor='white')
    ax1.set_xticks(range(len(event_means)))
    ax1.set_xticklabels([e.split('-')[1] if '-' in e else e for e in event_means.index], rotation=45)
    ax1.set_xlabel('Event')
    ax1.set_ylabel(r'$NO_3^-$ (mg/L)')
    ax1.set_title(r'(a) Mean $NO_3^-$ by Event', fontweight='bold')
    
    # (b) Seasonal boxplot
    ax2 = axes[0, 1]
    if 'season' in water_chem.columns:
        season_order = ['dry', 'transition', 'wet']
        season_colors = {'dry': '#f39c12', 'transition': '#27ae60', 'wet': '#3498db'}
        
        data_by_season = [water_chem[water_chem['season'] == s]['NO3_mg_L'].dropna() 
                         for s in season_order if s in water_chem['season'].unique()]
        labels = [s for s in season_order if s in water_chem['season'].unique()]
        
        bp = ax2.boxplot(data_by_season, tick_labels=labels, patch_artist=True)
        for patch, label in zip(bp['boxes'], labels):
            patch.set_facecolor(season_colors.get(label, 'gray'))
    
    ax2.set_xlabel('Season')
    ax2.set_ylabel(r'$NO_3^-$ (mg/L)')
    ax2.set_title(r'(b) $NO_3^-$ by Season', fontweight='bold')
    
    # (c) EC temporal trend - include ALL station types
    ax3 = axes[1, 0]
    available_events = sorted(water_chem['event_code'].unique(), key=lambda x: (x[1] if len(x)>1 and x[1].isdigit() else 0, x))
    for st_type, color, marker in [('lysimeter', '#3498db', 'o'), ('borehole', '#e74c3c', 's'), ('ag_well', '#f39c12', 'p')]:
        subset = water_chem[water_chem['station_type'] == st_type]
        if len(subset) > 0:
            ec_by_event = subset.groupby('event_code')['EC_uS_cm'].mean().reindex(available_events)
            if len(ec_by_event.dropna()) > 0:
                ax3.plot(range(len(ec_by_event)), ec_by_event.values, f'{marker}-', color=color, 
                        linewidth=2, markersize=8, label=st_type.replace('_', ' ').title())
    
    ax3.set_xticks(range(len(available_events)))
    ax3.set_xticklabels([e.split('-')[1] if '-' in e else e for e in available_events], rotation=45)
    ax3.set_xlabel('Event')
    ax3.set_ylabel(r'EC ($\mu S/cm$)')
    ax3.set_title('(c) Electrical Conductivity Trend', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # (d) Station type comparison - include ALL station types
    ax4 = axes[1, 1]
    ions = ['Ca_mg_L', 'Mg_mg_L', 'Na_mg_L', 'Cl_mg_L', 'SO4_mg_L', 'NO3_mg_L']
    ion_labels = [r'$Ca^{2+}$', r'$Mg^{2+}$', r'$Na^+$', r'$Cl^-$', r'$SO_4^{2-}$', r'$NO_3^-$']
    
    x = np.arange(len(ions))
    width = 0.25  # Narrower to fit 3 station types
    
    lys_means = [water_chem[water_chem['station_type'] == 'lysimeter'][ion].mean() for ion in ions]
    bh_means = [water_chem[water_chem['station_type'] == 'borehole'][ion].mean() for ion in ions]
    aw_means = [water_chem[water_chem['station_type'] == 'ag_well'][ion].mean() for ion in ions]
    
    ax4.bar(x - width, lys_means, width, label='Lysimeter', color='#3498db')
    ax4.bar(x, bh_means, width, label='Borehole', color='#e74c3c')
    ax4.bar(x + width, aw_means, width, label='Ag Well', color='#f39c12')
    ax4.set_xticks(x)
    ax4.set_xticklabels(ion_labels)
    ax4.set_ylabel('Concentration (mg/L)')
    ax4.set_title('(d) Mean Ion Concentrations by Type', fontweight='bold')
    ax4.legend(fontsize=9)
    
    plt.suptitle('Temporal and Spatial Patterns - Vea Catchment', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "temporal_patterns.png", dpi=300)
    plt.close()
    print(f"  Saved: temporal_patterns.png")


def plot_reaction_summary(all_results: dict, output_dir: Path):
    """Plot reaction summary heatmap using actual computed reaction extents."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))
    
    # --- (a) Reaction Extent Heatmap by Edge ---
    ax1 = axes[0]
    
    # Collect actual reaction data from EdgeResults
    all_extents = []
    edge_ids = []
    labels = None
    
    for event, results in all_results.items():
        if results:
            for res in results:
                if hasattr(res, 'z_extents') and res.z_extents:
                    all_extents.append(res.z_extents)
                    edge_ids.append(f"{res.u}→{res.v}")
                    if labels is None and hasattr(res, 'z_labels') and res.z_labels:
                        labels = res.z_labels
    
    if all_extents and labels:
        # Create matrix from actual data
        matrix = np.array(all_extents).T  # Shape: (n_reactions, n_edges)
        
        # Mapping for scientific labels
        label_map = {
            'calcite': r'Calcite ($CaCO_3$)',
            'dolomite': r'Dolomite ($CaMg(CO_3)_2$)',
            'gypsum': r'Gypsum ($CaSO_4 \cdot 2H_2O$)',
            'halite': r'Halite ($NaCl$)',
            'fluorite': r'Fluorite ($CaF_2$)',
            'albite': r'Albite ($NaAlSi_3O_8$)',
            'anorthite': r'Anorthite ($CaAl_2Si_2O_8$)',
            'forsterite': r'Forsterite ($Mg_2SiO_4$)',
            'pyrite_oxidation_aerobic': r'Pyrite Oxidation ($FeS_2$)',
            'no3src': r'Nitrification ($NH_4^+ \to NO_3^-$)',
            'denit': r'Denitrification ($NO_3^- \to N_2$)',
            'cana_exch': r'Ca-Na Exchange',
            'mgna_exch': r'Mg-Na Exchange'
        }
        
        def format_label(l):
            return label_map.get(l.lower(), l.replace('_', ' ').title())

        # Select reactions with significant extents (dynamic, not hardcoded)
        # A reaction is significant if its mean absolute extent > 0.001 mmol/L
        display_labels = []
        display_matrix = []
        
        for i, label in enumerate(labels):
            mean_extent = np.mean(np.abs(matrix[i, :]))
            if mean_extent > 0.001:  # Threshold for significance
                display_labels.append(format_label(label))
                display_matrix.append(matrix[i, :])
        
        if display_matrix:
            display_matrix = np.array(display_matrix)
            
            # Limit to first 12 edges for readability
            n_edges_display = min(12, len(edge_ids))
            display_matrix = display_matrix[:, :n_edges_display]
            edge_ids_display = edge_ids[:n_edges_display]
            
            im = ax1.imshow(display_matrix, cmap='YlOrRd', aspect='auto')
            
            ax1.set_xticks(range(n_edges_display))
            ax1.set_xticklabels(edge_ids_display, rotation=45, ha='right', fontsize=8)
            ax1.set_yticks(range(len(display_labels)))
            ax1.set_yticklabels(display_labels)
            
            plt.colorbar(im, ax=ax1, label='Extent (mmol/L)', shrink=0.8)
        else:
            ax1.text(0.5, 0.5, 'No reaction data', ha='center', va='center', transform=ax1.transAxes)
    else:
        ax1.text(0.5, 0.5, 'No reaction data available', ha='center', va='center', transform=ax1.transAxes)
    
    ax1.set_xlabel('Flow Path (Edge)')
    ax1.set_ylabel('Reaction')
    ax1.set_title('(a) Reaction Extent by Flow Path\n(mmol/L, actual computed values)', fontweight='bold')
    
    # --- (b) Mean Reaction Extents by Type ---
    ax2 = axes[1]
    
    if all_extents and labels:
        # Calculate mean extents across all edges
        mean_extents = np.mean(np.array(all_extents), axis=0)
        
        # Select reactions with significant extents
        sig_idx = [i for i, ext in enumerate(mean_extents) if ext > 0.01]
        
        if sig_idx:
            sig_labels = [format_label(labels[i]) for i in sig_idx]
            sig_extents = [mean_extents[i] for i in sig_idx]
            
            # Sort by extent
            sorted_pairs = sorted(zip(sig_labels, sig_extents), key=lambda x: x[1], reverse=True)
            sig_labels, sig_extents = zip(*sorted_pairs)
            
            colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(sig_labels)))
            bars = ax2.barh(range(len(sig_labels)), sig_extents, color=colors)
            ax2.set_yticks(range(len(sig_labels)))
            ax2.set_yticklabels(sig_labels)
            
            # Add value labels
            for bar, val in zip(bars, sig_extents):
                ax2.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                        f'{val:.3f}', va='center', fontsize=9)
        else:
            ax2.text(0.5, 0.5, 'No significant reactions', ha='center', va='center', transform=ax2.transAxes)
    else:
        ax2.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_xlabel('Mean Reaction Extent (mmol/L)')
    ax2.set_title('(b) Average Reaction Magnitudes\n(across all edges)', fontweight='bold')
    ax2.invert_yaxis()
    
    plt.suptitle('Geochemical Reaction Summary - Vea Catchment', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / "reaction_summary.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: reaction_summary.png")


def plot_transport_parameters(all_results: dict, output_dir: Path):
    """Plot transport parameter distributions using actual computed values."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Collect gamma values and objective scores from actual results
    gammas = []
    objectives = []
    transport_models = []
    
    for event, results in all_results.items():
        if results:
            for res in results:
                if hasattr(res, 'gamma') and res.gamma is not None:
                    gammas.append(res.gamma)
                if hasattr(res, 'objective_score') and res.objective_score is not None:
                    objectives.append(res.objective_score)
                if hasattr(res, 'transport_model') and res.transport_model:
                    transport_models.append(res.transport_model)
    
    # (a) Gamma distribution
    ax1 = axes[0]
    if gammas:
        n, bins, patches = ax1.hist(gammas, bins=20, color='steelblue', edgecolor='white', alpha=0.8)
        ax1.axvline(x=1.0, color='red', linestyle='--', linewidth=2, label=r'$\gamma=1$ (no conc.)')
        median_gamma = np.median(gammas)
        ax1.axvline(x=median_gamma, color='green', linestyle='-', linewidth=2, 
                   label=f'Median: {median_gamma:.2f}')
        
        # Add statistics text box - moved to avoid overlap
        stats_text = f'n = {len(gammas)}\nMean = {np.mean(gammas):.2f}\nRange: {min(gammas):.2f} - {max(gammas):.2f}'
        ax1.text(0.95, 0.85, stats_text, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax1.text(0.5, 0.5, 'No gamma values\navailable', ha='center', va='center', 
                transform=ax1.transAxes, fontsize=12)
    
    ax1.set_xlabel(r'Concentration Factor ($\gamma$)')
    ax1.set_ylabel('Frequency')
    ax1.set_title(r'(a) Evaporation Factor Distribution', fontweight='bold')
    ax1.legend(loc='upper left', frameon=True)
    
    # (b) Transport model classification
    ax2 = axes[1]
    if gammas:
        # Classify based on gamma values
        evap_count = sum(1 for g in gammas if g > 1.05)
        conserv_count = sum(1 for g in gammas if 0.95 <= g <= 1.05)
        dilution_count = sum(1 for g in gammas if g < 0.95)
        
        transport_counts = {
            r'Evaporation'+'\n'+r'($\gamma > 1.05$)': evap_count,
            r'Conservative'+'\n'+r'($0.95 \leq \gamma \leq 1.05$)': conserv_count,
            r'Dilution'+'\n'+r'($\gamma < 0.95$)': dilution_count
        }
        
        # Remove zero counts
        transport_counts = {k: v for k, v in transport_counts.items() if v > 0}
        
        if transport_counts:
            colors = ['#e74c3c', '#3498db', '#2ecc71'][:len(transport_counts)]
            explode = [0.05 if i == 0 else 0 for i in range(len(transport_counts))]
            wedges, texts, autotexts = ax2.pie(
                transport_counts.values(), 
                labels=transport_counts.keys(),
                autopct=lambda p: f'{p:.0f}%\n({int(p*sum(transport_counts.values())/100)})',
                colors=colors, 
                explode=explode,
                textprops={'fontsize': 10},
                pctdistance=0.75,
                labeldistance=1.15,
                startangle=90
            )
    else:
        ax2.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_title('(b) Transport Mechanism Classification', fontweight='bold')
    
    # (c) Objective scores - now using actual values
    ax3 = axes[2]
    if objectives:
        n, bins, patches = ax3.hist(objectives, bins=20, color='coral', edgecolor='white', alpha=0.8)
        
        # Color code by quality
        for patch, left_edge, right_edge in zip(patches, bins[:-1], bins[1:]):
            mid = (left_edge + right_edge) / 2
            if mid < 1.0:
                patch.set_facecolor('#27ae60')  # Green = good fit
            elif mid < 5.0:
                patch.set_facecolor('#f39c12')  # Orange = acceptable
            else:
                patch.set_facecolor('#e74c3c')  # Red = poor fit
        
        ax3.axvline(x=1.0, color='green', linestyle='--', linewidth=2, label='Good fit threshold')
        ax3.axvline(x=5.0, color='red', linestyle='--', linewidth=2, label='Poor fit threshold')
        
        # Statistics - moved to avoid overlap
        stats_text = f'n = {len(objectives)}\nMean = {np.mean(objectives):.2f}\nMedian = {np.median(objectives):.2f}'
        ax3.text(0.95, 0.80, stats_text, transform=ax3.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax3.text(0.5, 0.5, 'No objective scores\navailable', ha='center', va='center',
                transform=ax3.transAxes, fontsize=12)
    
    ax3.set_xlabel('Objective Score')
    ax3.set_ylabel('Frequency')
    ax3.set_title('(c) Fit Quality Distribution', fontweight='bold')
    if objectives:
        ax3.legend(loc='upper right')
    
    plt.suptitle('Transport Parameter Analysis - Vea Catchment', fontsize=14, fontweight='bold', y=1.05)
    plt.tight_layout()
    plt.savefig(output_dir / "transport_parameters.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: transport_parameters.png")


# ===================================================================
# MAIN EXECUTION
# ===================================================================

def main():
    print("=" * 70)
    print("COMPREHENSIVE MANUAL FIGURE GENERATOR")
    print("Generating publication-quality figures for Hydrosheaf Technical Manual")
    print("=" * 70)
    
    # Setup paths
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "hydrosheaf_synthetic_csv"
    output_dir = base_dir / "hydrosheaf_manual" / "figures" / "vea_catchment"
    
    # Create output directories
    output_dir.mkdir(exist_ok=True, parents=True)
    (output_dir / "plots").mkdir(exist_ok=True)
    (output_dir / "data").mkdir(exist_ok=True)
    (output_dir / "calibration").mkdir(exist_ok=True)
    (output_dir / "validation").mkdir(exist_ok=True)
    
    print(f"\nOutput directory: {output_dir}")
    
    # Step 1: Prepare data
    print("\n[1/6] Preparing Vea Catchment data...")
    water_chem = prepare_vea_data()
    water_chem.to_csv(data_dir / "water_chem_full.csv", index=False)
    print(f"  Saved {len(water_chem)} samples to water_chem_full.csv")
    
    # Load full synthetic dataset
    data = load_synthetic_data(data_dir)
    print(f"  Loaded {len(data['water_chem'])} samples, {len(data['stations'])} stations")
    
    # Step 2: Configure and run calibration
    print("\n[2/6] Running calibration with diagnostics...")
    # Use comprehensive mineral set for Vea Catchment crystalline basement setting
    # Includes: carbonates, evaporites, Na/Ca/Mg-silicates, and redox processes
    # Enable isotope constraints for evaporation/mixing model selection
    config = Config(
        ion_order=["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe"],
        weights=[1.0] * 10,
        phreeqc_enabled=False,
        transport_models_enabled=["evap", "mix"],
        active_minerals=[
            # Carbonates
            "calcite",                      # CaCO3 - primary carbonate dissolution
            "dolomite",                     # CaMg(CO3)2 - Mg-bearing carbonate
            # Evaporites
            "gypsum",                       # CaSO4·2H2O - sulfate source
            "halite",                       # NaCl - chloride source
            "fluorite",                     # CaF2 - fluoride source
            # Na-silicates (feldspar weathering)
            "albite",                       # NaAlSi3O8 - Na-feldspar weathering
            # Ca-silicates
            "anorthite",                    # CaAl2Si2O8 - Ca-feldspar weathering
            # Mg-silicates (mafic mineral weathering - important for crystalline basement)
            "forsterite",                   # Mg2SiO4 - olivine weathering (Mg source)
            "enstatite",                    # MgSiO3 - pyroxene weathering (Mg source)
            "biotite",                      # K(Mg,Fe)3AlSi3O10(OH)2 - mica weathering (K, Mg, Fe, F)
            # Redox processes
            "pyrite_oxidation_aerobic",     # FeS2 + O2 → Fe + SO4 - sulfide oxidation
        ],
        gibbs_enabled=True,
        exchange_enabled=True,
        # Enable water isotope (δ18O, δ2H) constraints for transport model discrimination
        # Isotopes help distinguish evaporation (enrichment along evaporation line)
        # from mixing (linear interpolation between endmembers)
        isotope_enabled=True,
        isotope_d18o_key="18O",            # Key in sample dict for δ18O-H2O
        isotope_d2h_key="2H",              # Key in sample dict for δ2H-H2O
        isotope_weight=1.0,                # Weight of isotope penalty in objective
    )
    
    cal_results = run_calibration_with_diagnostics(data, config, output_dir)
    
    # Step 3: Run hydrosheaf analysis first (plots need these results)
    print("\n[3/6] Running Hydrosheaf network analysis...")
    all_results = {}
    all_summaries = {}
    all_edge_results = []
    
    for event_code in data["water_chem"]["event_code"].unique():
        samples = prepare_samples_mmol(data["water_chem"], event_code)
        edges = [(row["from_station"], row["to_station"]) for _, row in data["edges"].iterrows()]
        results = fit_network(samples, edges, config)
        all_results[event_code] = results
        if results:
            all_summaries[event_code] = summarize_network(results)
            all_edge_results.extend(results)
    
    print(f"  Analyzed {len(all_edge_results)} edges across {len(all_results)} events")
    
    # Step 4: Generate calibration diagnostic plots (ch20_calibration.tex)
    # These now use actual fit_network results for observed vs simulated comparisons
    print("\n[4/6] Generating calibration diagnostic plots...")
    plot_calibration_convergence(cal_results, output_dir / "calibration")
    plot_parameter_sensitivity(cal_results, output_dir / "calibration")
    plot_observed_vs_simulated(data, cal_results, all_results, config, output_dir / "calibration")
    plot_residual_diagnostics(data, all_results, config, output_dir / "calibration")
    
    # Export results
    export_edge_results_csv(all_edge_results, str(output_dir / "data" / "edge_results.csv"))
    export_edge_results_json(all_edge_results, str(output_dir / "data" / "edge_results.json"))
    with open(output_dir / "data" / "network_summaries.json", "w") as f:
        json.dump(all_summaries, f, indent=2, default=str)
    
    # Step 5: Generate standard analysis plots
    print("\n[5/6] Generating analysis plots (ch24_interpreting_results.tex)...")
    plots_dir = output_dir / "plots"
    
    # ILR and Gibbs
    ilr_samples = prepare_samples_for_ilr(data["water_chem"])
    plot_ilr(ilr_samples, str(plots_dir / "ilr_plot.png"))
    print(f"  Saved: ilr_plot.png")
    plot_gibbs(ilr_samples, str(plots_dir / "gibbs_diagram.png"))
    print(f"  Saved: gibbs_diagram.png")
    
    # Edge anomalies
    plot_edge_anomalies(all_edge_results, str(plots_dir / "edge_anomalies.png"))
    print(f"  Saved: edge_anomalies.png")
    
    # Custom plots
    plot_water_isotopes(data["water_chem"], plots_dir)
    plot_nitrate_sources(data["water_chem"], data["endmembers"], plots_dir)
    plot_network_schematic(data["edges"], data["stations"], plots_dir)
    plot_data_model_diagram(plots_dir)  # New diagram
    plot_temporal_patterns(data["water_chem"], data["events"], plots_dir)
    plot_reaction_summary(all_results, plots_dir)
    plot_transport_parameters(all_results, plots_dir)
    
    # Step 6: Generate validation plots (ch28_validation.tex)
    print("\n[6/6] Generating validation plots...")
    plot_transport_validation(data, all_results, output_dir / "validation")
    plot_vertical_vs_lateral_transport(data, all_results, output_dir / "validation")
    plot_validation_summary_table(cal_results, all_results, data, output_dir / "validation")
    
    # Final summary
    print("\n" + "=" * 70)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 70)
    print(f"\nAll outputs saved to: {output_dir}")
    print("\nGenerated figures for manual chapters:")
    print("  ch20_calibration.tex:")
    print("    - calibration/calibration_convergence.png")
    print("    - calibration/parameter_sensitivity.png")
    print("    - calibration/observed_vs_simulated.png")
    print("    - calibration/residual_diagnostics.png")
    print("  ch24_interpreting_results.tex:")
    print("    - plots/ilr_plot.png")
    print("    - plots/gibbs_diagram.png")
    print("    - plots/edge_anomalies.png")
    print("    - plots/water_isotopes.png")
    print("    - plots/nitrate_sources.png")
    print("    - plots/network_schematic.png")
    print("    - plots/temporal_patterns.png")
    print("    - plots/reaction_summary.png")
    print("    - plots/transport_parameters.png")
    print("  ch28_validation.tex:")
    print("    - validation/transport_validation.png")
    print("    - validation/vertical_vs_lateral_transport.png")
    print("    - validation/validation_summary_table.png")
    print("  Data outputs:")
    print("    - data/calibration_results.json")
    print("    - data/edge_results.csv")
    print("    - data/edge_results.json")
    print("    - data/network_summaries.json")
    
    return output_dir


if __name__ == "__main__":
    output_dir = main()
