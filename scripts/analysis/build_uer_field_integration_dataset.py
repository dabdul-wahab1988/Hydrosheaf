"""Build complete, harmonized, and organized field integration dataset for UER (Northern Ghana New).

Combines:
1. Native groundwater chemistry and isotope measurements from 'compiled UER data_new.xlsx'.
2. Audited WGS84 coordinates and Ghana Geological Survey polygon joins.
3. Completed 30m DEM elevations (SRTM 30m / ASTER 30m from 'uer_elevations_dem.csv').
4. Projected UTM Zone 30N coordinates (EPSG:32630).
5. Molar (mmol/L) and equivalent (meq/L) concentrations, charge balance errors, and hydrochemical facies.
6. Isotopic diagnostic indicators (d-excess, tritium recharge regime, nitrate source classification).
7. Auxiliary rainfall and monitoring-well time series.

Outputs:
- data/FieldData/derived/uer_field_integration_dataset.csv
- data/FieldData/NorthernGhanaNew/compiled UER data_new_completed.xlsx
"""

from __future__ import annotations

import datetime
from pathlib import Path
import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter
import pandas as pd
import pyproj

ROOT = Path(__file__).resolve().parents[2]
RAW_XLSX = ROOT / "data" / "FieldData" / "NorthernGhanaNew" / "compiled UER data_new.xlsx"
JOIN_CSV = ROOT / "outputs" / "2026-09-05_northern_ghana_geology_join" / "NorthernGhanaNew_geology_join.csv"
DEM_CSV = ROOT / "data" / "FieldData" / "derived" / "uer_elevations_dem.csv"

OUT_DIR_DERIVED = ROOT / "data" / "FieldData" / "derived"
OUT_CSV = OUT_DIR_DERIVED / "uer_field_integration_dataset.csv"
OUT_XLSX = ROOT / "data" / "FieldData" / "NorthernGhanaNew" / "compiled UER data_new_completed.xlsx"

# Molar masses (g/mol) and valences
MOLAR_MASS = {
    "Ca": (40.078, 2),
    "Mg": (24.305, 2),
    "Na": (22.98977, 1),
    "K": (39.0983, 1),
    "HCO3": (61.0168, 1),
    "Cl": (35.453, 1),
    "SO4": (96.06, 2),
    "NO3": (62.0049, 1),
    "F": (18.9984, 1),
}


def classify_nitrate_source(row: pd.Series) -> str:
    n15 = row.get("d15N_NO3_permil_air")
    o18 = row.get("d18O_NO3_permil_VSMOW")
    no3_mg = row.get("no3_mg_L")
    if pd.isna(n15) or pd.isna(o18):
        if pd.notna(no3_mg) and no3_mg < 10.0:
            return "Low Nitrate (<10 mg/L) - unmeasured isotopes"
        return "Unmeasured"

    if n15 > 8.5:
        return "Manure / Septic Waste (d15N > +8.5 permil)"
    elif 3.0 <= n15 <= 8.5 and o18 < 15.0:
        return "Soil Organic N / Mixed Agricultural (d15N +3.0 to +8.5 permil)"
    elif n15 < 3.0 and o18 > 15.0:
        return "Synthetic Fertilizer / Atmospheric (high d18O, low d15N)"
    else:
        return "Mixed Agricultural / Denitrified"


def classify_tritium(tu: float | None) -> str:
    if pd.isna(tu):
        return "Unmeasured"
    if tu >= 1.0:
        return "Modern Recharge (TU >= 1.0)"
    elif tu <= 0.8:
        return "Sub-modern / Pre-bomb (TU <= 0.8)"
    else:
        return "Mixed / Intermediate (0.8 < TU < 1.0)"


def main() -> None:
    raise RuntimeError(
        "Retired: this builder reads the superseded NorthernGhanaNew folder and "
        "recreates a derived UER mirror from an old source. Use the approved "
        "UERdata/compiled UER data_new_completed.xlsx source instead."
    )
    print(f"Loading {JOIN_CSV.name}...")
    df_join = pd.read_csv(JOIN_CSV)
    print(f"Loading {DEM_CSV.name}...")
    df_dem = pd.read_csv(DEM_CSV)

    # Merge DEM data
    df = df_join.merge(
        df_dem[["sample_no", "elev_srtm30m", "elev_aster30m", "elevation_completed_m", "elevation_source"]],
        on="sample_no",
        how="left",
    )

    # Project to UTM Zone 30N (EPSG:32630)
    transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:32630", always_xy=True)
    utm_e, utm_n = transformer.transform(
        df["longitude_candidate_dd"].values,
        df["latitude_candidate_dd"].values,
    )
    df["utm_easting_m"] = np.round(utm_e, 2)
    df["utm_northing_m"] = np.round(utm_n, 2)

    # Identifiers
    df["node_id"] = [f"NGN_{int(no):03d}" for no in df["sample_no"]]
    df["site_id"] = df["node_id"]
    type_map = {"BH": "Borehole", "HDW": "Hand-Dug Well", "SW": "Surface Water"}
    df["sample_type_desc"] = df["sample_type"].map(type_map).fillna("Unknown")

    # Molar concentrations (mmol/L)
    for ion, (mass, _) in MOLAR_MASS.items():
        src_col = f"{ion.lower()}_mg_L"
        if src_col in df.columns:
            df[f"{ion.lower()}_mmol_L"] = np.round(df[src_col] / mass, 5)

    # Equivalent concentrations (meq/L)
    for ion, (mass, z) in MOLAR_MASS.items():
        src_col = f"{ion.lower()}_mg_L"
        if src_col in df.columns:
            df[f"{ion.lower()}_meq_L"] = np.round((df[src_col] / mass) * z, 5)

    # Recalculate Sum Cations, Sum Anions, CBE (%)
    sum_cat = (
        df["ca_meq_L"].fillna(0)
        + df["mg_meq_L"].fillna(0)
        + df["na_meq_L"].fillna(0)
        + df["k_meq_L"].fillna(0)
    )
    sum_an = (
        df["hco3_meq_L"].fillna(0)
        + df["cl_meq_L"].fillna(0)
        + df["so4_meq_L"].fillna(0)
        + df["no3_meq_L"].fillna(0)
        + df["f_meq_L"].fillna(0)
    )
    cbe = np.where((sum_cat + sum_an) > 0, ((sum_cat - sum_an) / (sum_cat + sum_an)) * 100.0, np.nan)
    df["sum_cations_meq_L"] = np.round(sum_cat, 4)
    df["sum_anions_meq_L"] = np.round(sum_an, 4)
    df["cbe_percent"] = np.round(cbe, 2)
    df["cbe_class"] = np.where(
        np.abs(cbe) <= 5.0,
        "ACCEPTABLE (CBE <= 5%)",
        np.where(np.abs(cbe) <= 10.0, "CAUTION (5% < CBE <= 10%)", "REJECT (CBE > 10%)"),
    )

    # Coordinate quality descriptor matching hydrosheaf contract
    df["coordinate_quality"] = np.where(
        df["coordinate_assignment_method"] == "range",
        "range_reordered_source_labels",
        "decimal_or_dms_native",
    )

    # Water Isotopes and d-excess
    d18 = df["d18O_permil"]
    d2 = df["d2H_permil"]
    df["d_excess_permil"] = np.round(d2 - 8.0 * d18, 2)

    # Tritium classification
    df["recharge_era_tritium"] = [classify_tritium(val) for val in df["tritium_TU"]]

    # Nitrate source classification
    df["nitrate_source_candidate"] = [classify_nitrate_source(row) for _, row in df.iterrows()]

    # Hydrosheaf model attributes
    df["dataset"] = "northern_ghana_new"
    df["flow_head_proxy_m"] = df["elevation_completed_m"]
    df["head_inference_tier"] = "C"

    # Select and order canonical fields
    canonical_columns = [
        "node_id",
        "site_id",
        "sample_no",
        "community",
        "sample_type",
        "sample_type_desc",
        "latitude_candidate_dd",
        "longitude_candidate_dd",
        "utm_easting_m",
        "utm_northing_m",
        "elevation_completed_m",
        "elevation_source",
        "elevation_m",
        "elev_srtm30m",
        "elev_aster30m",
        "coordinate_quality",
        "pH",
        "ec_uS_cm",
        "tds_mg_L",
        "ca_mg_L",
        "mg_mg_L",
        "na_mg_L",
        "k_mg_L",
        "hco3_mg_L",
        "cl_mg_L",
        "so4_mg_L",
        "no3_mg_L",
        "f_mg_L",
        "ca_mmol_L",
        "mg_mmol_L",
        "na_mmol_L",
        "k_mmol_L",
        "hco3_mmol_L",
        "cl_mmol_L",
        "so4_mmol_L",
        "no3_mmol_L",
        "f_mmol_L",
        "sum_cations_meq_L",
        "sum_anions_meq_L",
        "cbe_percent",
        "cbe_class",
        "cation_facies",
        "anion_facies",
        "facies_label",
        "d18O_permil",
        "d2H_permil",
        "d_excess_permil",
        "tritium_TU",
        "recharge_era_tritium",
        "d15N_NO3_permil_air",
        "d18O_NO3_permil_VSMOW",
        "nitrate_source_candidate",
        "geology_join_status",
        "geology_stratigraphic_unit",
        "geology_stratigraphic_formation",
        "geology_symbol",
        "geology_tectonic_domain",
        "geology_metamorphic_grade",
        "geology_boundary_distance_m",
        "dataset",
        "flow_head_proxy_m",
        "head_inference_tier",
    ]

    rename_map = {
        "latitude_candidate_dd": "latitude_dd",
        "longitude_candidate_dd": "longitude_dd",
        "elevation_completed_m": "elevation_m",
        "elevation_m": "elevation_surveyed_m",
    }

    # Prepare export dataframe
    df_out = df[canonical_columns].rename(columns=rename_map)

    # Save canonical CSV
    OUT_DIR_DERIVED.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(OUT_CSV, index=False)
    print(f"Wrote canonical CSV to {OUT_CSV} ({len(df_out)} rows, {len(df_out.columns)} cols)")

    # Read auxiliary tables from native workbook
    print(f"Reading native auxiliary tables from {RAW_XLSX.name}...")
    df_rain = pd.read_excel(RAW_XLSX, sheet_name="rain")
    df_rain = df_rain.dropna(how="all").copy()
    # Drop completely empty columns
    df_rain = df_rain.dropna(axis=1, how="all")

    df_mw = pd.read_excel(RAW_XLSX, sheet_name="monitoring wells")
    df_mw = df_mw.dropna(how="all").copy()

    # Create Data Dictionary
    dict_rows = [
        ("node_id", "Text", "-", "Unique graph node identifier for spatial and topological network modeling (NGN_001..237)"),
        ("site_id", "Text", "-", "Site identifier (matches node_id in NorthernGhanaNew)"),
        ("sample_no", "Integer", "-", "Sequential sample number as designated in original survey (1..237)"),
        ("community", "Text", "-", "Local community, settlement, or well installation name"),
        ("sample_type", "Code", "-", "BH = Borehole, HDW = Hand-Dug Well, SW = Surface Water"),
        ("sample_type_desc", "Text", "-", "Full descriptive name of sample source"),
        ("latitude_dd", "Float", "degrees North", "Audited WGS84 decimal latitude (re-ordered from reversed native column)"),
        ("longitude_dd", "Float", "degrees West", "Audited WGS84 decimal longitude (signed negative degrees West)"),
        ("utm_easting_m", "Float", "meters", "Projected UTM Zone 30N Easting coordinate (EPSG:32630)"),
        ("utm_northing_m", "Float", "meters", "Projected UTM Zone 30N Northing coordinate (EPSG:32630)"),
        ("elevation_m", "Float", "m a.s.l.", "Completed surface elevation: surveyed where available, SRTM 30m DEM otherwise"),
        ("elevation_source", "Text", "-", "Source of primary elevation value ('surveyed' or 'dem_srtm30m')"),
        ("elevation_surveyed_m", "Float", "m a.s.l.", "Surveyed ground elevation recorded in native table (2 values present)"),
        ("elev_srtm30m", "Float", "m a.s.l.", "NASA SRTM 1-Arc-Second (30m) global digital elevation model"),
        ("elev_aster30m", "Float", "m a.s.l.", "METI/NASA ASTER GDEM v3 (30m) digital elevation model"),
        ("coordinate_quality", "Text", "-", "Quality and normalization method applied to raw coordinate text"),
        ("pH", "Float", "pH units", "Field measured groundwater pH at sampling temperature"),
        ("ec_uS_cm", "Float", "uS/cm", "Electrical conductivity normalized to 25 deg C"),
        ("tds_mg_L", "Float", "mg/L", "Total Dissolved Solids measured in field/laboratory"),
        ("ca_mg_L", "Float", "mg/L", "Dissolved Calcium concentration (mg/L)"),
        ("mg_mg_L", "Float", "mg/L", "Dissolved Magnesium concentration (mg/L)"),
        ("na_mg_L", "Float", "mg/L", "Dissolved Sodium concentration (mg/L)"),
        ("k_mg_L", "Float", "mg/L", "Dissolved Potassium concentration (mg/L)"),
        ("hco3_mg_L", "Float", "mg/L", "Alkalinity / Bicarbonate concentration (mg/L)"),
        ("cl_mg_L", "Float", "mg/L", "Dissolved Chloride concentration (mg/L)"),
        ("so4_mg_L", "Float", "mg/L", "Dissolved Sulfate concentration (mg/L)"),
        ("no3_mg_L", "Float", "mg/L", "Dissolved Nitrate concentration (mg/L)"),
        ("f_mg_L", "Float", "mg/L", "Dissolved Fluoride concentration (mg/L; 188 observed, 49 missing)"),
        ("ca_mmol_L", "Float", "mmol/L", "Dissolved Calcium in millimoles per liter"),
        ("mg_mmol_L", "Float", "mmol/L", "Dissolved Magnesium in millimoles per liter"),
        ("na_mmol_L", "Float", "mmol/L", "Dissolved Sodium in millimoles per liter"),
        ("k_mmol_L", "Float", "mmol/L", "Dissolved Potassium in millimoles per liter"),
        ("hco3_mmol_L", "Float", "mmol/L", "Bicarbonate in millimoles per liter"),
        ("cl_mmol_L", "Float", "mmol/L", "Chloride in millimoles per liter"),
        ("so4_mmol_L", "Float", "mmol/L", "Sulfate in millimoles per liter"),
        ("no3_mmol_L", "Float", "mmol/L", "Nitrate in millimoles per liter"),
        ("f_mmol_L", "Float", "mmol/L", "Fluoride in millimoles per liter"),
        ("sum_cations_meq_L", "Float", "meq/L", "Total cation charge equivalents: 2*Ca + 2*Mg + Na + K"),
        ("sum_anions_meq_L", "Float", "meq/L", "Total anion charge equivalents: HCO3 + Cl + 2*SO4 + NO3 + F"),
        ("cbe_percent", "Float", "%", "Normalized Charge Balance Error: (Cat - An)/(Cat + An) * 100"),
        ("cbe_class", "Text", "-", "Standard geochemical reliability classification of charge balance error"),
        ("cation_facies", "Text", "-", "Dominant cation classification (Ca, Na, Mg, Mixed)"),
        ("anion_facies", "Text", "-", "Dominant anion classification (HCO3, Cl, SO4, Mixed)"),
        ("facies_label", "Text", "-", "Hydrochemical water facies classification (e.g. Ca-HCO3, Na-HCO3)"),
        ("d18O_permil", "Float", "permil VSMOW", "Stable oxygen-18 isotope delta value of water"),
        ("d2H_permil", "Float", "permil VSMOW", "Stable hydrogen-2 (deuterium) isotope delta value of water"),
        ("d_excess_permil", "Float", "permil", "Dansgaard Deuterium excess: d = d2H - 8 * d18O"),
        ("tritium_TU", "Float", "TU", "Tritium activity in Tritium Units (1 TU = 1 3H per 10^18 1H atoms)"),
        ("recharge_era_tritium", "Text", "-", "Qualitative modern vs sub-modern recharge regime based on tritium"),
        ("d15N_NO3_permil_air", "Float", "permil Air", "Nitrate nitrogen-15 isotope ratio"),
        ("d18O_NO3_permil_VSMOW", "Float", "permil VSMOW", "Nitrate oxygen-18 isotope ratio"),
        ("nitrate_source_candidate", "Text", "-", "Candidate nitrate origin inferred from dual nitrate isotopes and concentration"),
        ("geology_join_status", "Text", "-", "Spatial join status against Ghana Geological Survey polygon layer"),
        ("geology_stratigraphic_unit", "Text", "-", "Mapped stratigraphic group/unit (e.g. Birimian Supergroup, Plutonic Suites)"),
        ("geology_stratigraphic_formation", "Text", "-", "Specific geological formation name"),
        ("geology_symbol", "Text", "-", "Geological map symbol code (e.g. tmht, gskf, gvbm)"),
        ("geology_tectonic_domain", "Text", "-", "Regional tectonic framework / crustal domain"),
        ("geology_metamorphic_grade", "Text", "-", "Metamorphic grade of host crystalline basement rock"),
        ("geology_boundary_distance_m", "Float", "meters", "Distance from sample coordinate to nearest geological polygon boundary"),
        ("dataset", "Text", "-", "Package name in Hydrosheaf framework ('northern_ghana_new')"),
        ("flow_head_proxy_m", "Float", "m a.s.l.", "Topographic elevation used as Tier C hydraulic head proxy"),
        ("head_inference_tier", "Code", "-", "HydroSheaf head inference tier ('C' = topography proxy)"),
    ]
    df_dict = pd.DataFrame(dict_rows, columns=["Variable", "Data Type", "Units", "Description"])

    # Build multi-tab Excel workbook with professional styling
    print(f"Building styled Excel workbook at {OUT_XLSX.name}...")
    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        df_out.to_excel(writer, sheet_name="GW_Field_Integration", index=False)
        
        # Molar and charge view
        molar_cols = [
            "node_id", "community", "sample_type_desc", "elevation_m", "pH",
            "ca_mmol_L", "mg_mmol_L", "na_mmol_L", "k_mmol_L",
            "hco3_mmol_L", "cl_mmol_L", "so4_mmol_L", "no3_mmol_L", "f_mmol_L",
            "sum_cations_meq_L", "sum_anions_meq_L", "cbe_percent", "cbe_class", "facies_label"
        ]
        df_out[molar_cols].to_excel(writer, sheet_name="GW_Molar_Charge", index=False)
        
        df_rain.to_excel(writer, sheet_name="rain", index=False)
        df_mw.to_excel(writer, sheet_name="monitoring wells", index=False)
        df_dict.to_excel(writer, sheet_name="Data_Dictionary", index=False)

    # Style openpyxl workbook
    wb = openpyxl.load_workbook(OUT_XLSX)
    header_fill = PatternFill(start_color="1F497D", end_color="1F497D", fill_type="solid")
    header_font = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
    thin_border = Border(
        left=Side(style="thin", color="D9D9D9"),
        right=Side(style="thin", color="D9D9D9"),
        top=Side(style="thin", color="D9D9D9"),
        bottom=Side(style="thin", color="D9D9D9"),
    )

    for sheetname in wb.sheetnames:
        ws = wb[sheetname]
        ws.views.sheetView[0].showGridLines = True
        for col_idx, col in enumerate(ws.iter_cols(min_row=1, max_row=1), start=1):
            cell = col[0]
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

        # Auto-adjust column widths
        for col in ws.columns:
            max_len = max(len(str(cell.value or "")) for cell in col)
            col_letter = get_column_letter(col[0].column)
            ws.column_dimensions[col_letter].width = min(max(max_len + 3, 11), 50)

    wb.save(OUT_XLSX)
    print(f"Successfully generated {OUT_XLSX} with {len(wb.sheetnames)} styled sheets.")


if __name__ == "__main__":
    main()
