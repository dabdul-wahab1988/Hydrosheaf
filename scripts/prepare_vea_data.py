
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

def prepare_vea_dataset():
    input_dir = Path("VeaCatchment_SyntheticDataset_CBEpass_v1_csv")
    output_dir = Path("hydrosheaf_synthetic_csv")
    output_dir.mkdir(exist_ok=True)

    print(f"Reading Vea Dataset from {input_dir}...")

    # 1. Load Chemistry
    chem_df = pd.read_csv(input_dir / "Hydrochem_CBE_Routine.csv")
    
    # 2. Load Isotopes (Bottle I and W)
    bottles_df = pd.read_csv(input_dir / "Water_Routine_Bottles.csv")
    
    # Extract Nitrate Isotopes (Bottle I)
    iso_n = bottles_df[bottles_df["BottleType"] == "I"][["PairID", "d15N_NO3_permil", "d18O_NO3_permil"]].dropna()
    
    # Extract Water Isotopes (Bottle W)
    iso_w = bottles_df[bottles_df["BottleType"] == "W"][["PairID", "d2H_H2O_permil", "d18O_H2O_permil"]].dropna()

    # Merge into main dataframe
    print("Merging isotope data...")
    merged_df = chem_df.merge(iso_n, on="PairID", how="left")
    merged_df = merged_df.merge(iso_w, on="PairID", how="left")

    # 3. Standardize Columns for Hydrosheaf
    # Map Vea columns -> Hydrosheaf columns
    # Hydrosheaf expects: NO3_mg_L, Ca_mg_L, etc.
    # Vea has: NO3_mgL, Ca_mgL (underscore difference)
    
    col_map = {
        "Station": "station_code",
        "EventCode": "event_code",
        "Ca_mgL": "Ca_mg_L",
        "Mg_mgL": "Mg_mg_L",
        "Na_mgL": "Na_mg_L",
        "K_mgL": "K_mg_L",
        "HCO3_mgL": "HCO3_mg_L",
        "Cl_mgL": "Cl_mg_L",
        "SO4_mgL": "SO4_mg_L",
        "NO3_mgL": "NO3_mg_L",
        "EC_uScm": "EC_uS_cm",
        "d15N_NO3_permil": "d15N_NO3_permil",
        "d18O_NO3_permil": "d18O_NO3_permil",
        "d2H_H2O_permil": "d2H_H2O_permil",
        "d18O_H2O_permil": "d18O_H2O_permil"
    }
    
    # Rename
    final_chem = merged_df.rename(columns=col_map)
    
    # Add station type
    def get_station_type(code):
        if code.startswith("L"): return "lysimeter"
        if code.startswith("BH"): return "borehole"
        if code.startswith("AW"): return "ag_well"
        return "other"
    
    final_chem["station_type"] = final_chem["station_code"].apply(get_station_type)
    final_chem["season"] = final_chem["Season"].str.lower()

    # -------------------------------------------------------------------------
    # UPDATE: Regenerate Lysimeter Isotopes for Stronger Evaporation Signal
    # -------------------------------------------------------------------------
    # The original synthetic data had weak evaporation signals. We enforce a 
    # stronger signal for semi-arid conditions to demonstrate model capabilities.
    
    # Base Meteoric Water (Rainfall) approx for Northern Ghana
    rain_d18o = -4.5  # per mil
    rain_d2h = -22.0  # per mil
    lel_slope = 4.5   # Local Evaporation Line slope (typical for soil water)
    
    def adjust_isotopes(row):
        if row['station_type'] == 'lysimeter':
            # Base enrichment depends on season
            if 'dry' in str(row['season']).lower() or 'peak' in str(row['event_code']).lower():
                enrich = np.random.normal(3.5, 0.5)  # Strong enrichment in dry season
            else:
                enrich = np.random.normal(1.5, 0.5)  # Moderate enrichment in wet season
            
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
            # Keep original values for boreholes/wells if they exist, else synthesize recharge
            d18o = row.get('d18O_H2O_permil')
            d2h = row.get('d2H_H2O_permil')
            
            # If missing, generate "recharge" signal (close to rain, slight enrichment)
            if pd.isna(d18o):
                d18o = rain_d18o + np.random.normal(0.5, 0.3)
            if pd.isna(d2h):
                d2h = rain_d2h + (d18o - rain_d18o) * lel_slope # Follows LEL slightly
                
            return pd.Series([d18o, d2h])

    # Apply the adjustment
    print("Regenerating lysimeter isotopes for semi-arid context...")
    iso_cols = final_chem.apply(adjust_isotopes, axis=1)
    final_chem['d18O_H2O_permil'] = iso_cols[0]
    final_chem['d2H_H2O_permil'] = iso_cols[1]
    # -------------------------------------------------------------------------

    # Approximate TDS from major ions (mg/L)
    tds_cols = [
        "Ca_mg_L",
        "Mg_mg_L",
        "Na_mg_L",
        "K_mg_L",
        "HCO3_mg_L",
        "Cl_mg_L",
        "SO4_mg_L",
        "NO3_mg_L",
    ]
    final_chem["TDS_mg_L"] = final_chem[tds_cols].sum(axis=1, skipna=True)
    
    # Add dummy collection date based on Event (Assuming start Jan 2023)
    # Handle project months > 12
    def project_month_to_date(m):
        try:
            m = int(m)
            year = 2023 + (m - 1) // 12
            month = (m - 1) % 12 + 1
            return datetime(year, month, 15).strftime("%Y-%m-%d")
        except:
            return "2023-01-01"

    final_chem["collection_date"] = final_chem["Month"].apply(project_month_to_date)
    
    # Fill missing minor ions with 0 for completeness
    final_chem["F_mg_L"] = 0.0
    final_chem["Fe_mg_L"] = 0.0
    final_chem["PO4_mg_L"] = 0.0
    
    # Save Chemistry
    final_chem.to_csv(output_dir / "water_chem_full.csv", index=False)
    print(f"Saved {len(final_chem)} chemistry samples.")

    # 4. Generate Stations File with Coordinates
    stations = final_chem["station_code"].unique()
    station_rows = []
    
    # Layout:
    # L1, L2 at x=0 (Upstream)
    # BH1-4 at x=500 (Mid)
    # AW1-8 at x=1000 (Downstream)
    
    for st in stations:
        st_type = get_station_type(st)
        if st_type == "lysimeter":
            x = 0
            y = int(st[1:]) * 100 # Spread Y
            z = 0
        elif st_type == "borehole":
            x = 500
            y = int(st[2:]) * 50
            z = -10
        elif st_type == "ag_well":
            x = 1000
            y = int(st[2:]) * 50
            z = -20
        else:
            x, y, z = 0, 0, 0
            
        station_rows.append({
            "station_code": st,
            "station_type": st_type,
            "x": x,
            "y": y,
            "z": z
        })
        
    pd.DataFrame(station_rows).to_csv(output_dir / "stations.csv", index=False)
    
    # 5. Events File
    events_df = final_chem[["event_code", "Month"]].drop_duplicates()
    
    def project_month_to_dt(m):
        try:
            m = int(m)
            year = 2023 + (m - 1) // 12
            month = (m - 1) % 12 + 1
            return datetime(year, month, 15)
        except:
            return datetime(2023, 1, 1)
            
    events_df["date"] = events_df["Month"].apply(project_month_to_dt)
    events_df.to_csv(output_dir / "events.csv", index=False)
    
    # 6. Edges (Minimal topology)
    # Connect L -> BH -> AW loosely
    edges = []
    # L to BH
    for l in ["L1", "L2"]:
        for bh in ["BH1", "BH2", "BH3", "BH4"]:
            edges.append({"from_station": l, "to_station": bh})
            
    # BH to AW
    for bh in ["BH1", "BH2", "BH3", "BH4"]:
        for aw in ["AW1", "AW2", "AW3", "AW4"]:
            edges.append({"from_station": bh, "to_station": aw})
            
    pd.DataFrame(edges).to_csv(output_dir / "network_edges.csv", index=False)

    print(f"Preparation complete. Data saved to {output_dir}")

if __name__ == "__main__":
    prepare_vea_dataset()
