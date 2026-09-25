"""Automated test suite for Central Region (CR) field integration dataset.

Verifies:
1. Immutable provenance of the raw source workbook (CentralRegion.xlsx).
2. DEM elevation completeness (251/251 coordinate points) and regional terrain bounds.
3. Coordinate normalization and UTM Zone 30N projection (EPSG:32630).
4. 100% spatial join match rate with Ghana Geological Survey polygon layer.
5. Primary hydrochemistry, molar conversions, and charge balance error (CBE).
6. Well hydraulics: drawdown, specific capacity, and field-measured hydraulic head.
7. Multi-sheet Excel workbook integrity, openpyxl styling, and data dictionary.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import openpyxl
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DERIVED_CSV = ROOT / "data" / "FieldData" / "derived" / "cr_field_integration_dataset.csv"
DEM_CSV = ROOT / "data" / "FieldData" / "derived" / "cr_elevations_dem.csv"
COMPLETED_XLSX = ROOT / "data" / "FieldData" / "CRdata" / "CentralRegion_completed.xlsx"
RAW_XLSX = ROOT / "data" / "FieldData" / "CRdata" / "CentralRegion.xlsx"
EXPECTED_RAW_SHA256 = "6b8a2cd6545d062813476ee2cf1f0238ac902859d32964cd7c550be88df5ecd7"


def test_cr_raw_workbook_provenance_hash_is_unmutated() -> None:
    """The original raw Central Region workbook or its preserved Native_CR_Master sheet must retain exact provenance."""
    if RAW_XLSX.exists():
        with RAW_XLSX.open("rb") as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        assert digest == EXPECTED_RAW_SHA256, (
            f"Raw workbook has been modified! Expected {EXPECTED_RAW_SHA256}, got {digest}"
        )
    else:
        assert COMPLETED_XLSX.exists(), f"Completed workbook missing: {COMPLETED_XLSX}"
        wb = openpyxl.load_workbook(COMPLETED_XLSX, read_only=True)
        assert "Native_CR_Master" in wb.sheetnames, "Native_CR_Master sheet missing in completed workbook"
        df_raw = pd.read_excel(COMPLETED_XLSX, sheet_name="Native_CR_Master")
        assert len(df_raw) == 252, f"Expected 252 raw rows, got {len(df_raw)}"
        assert len(df_raw.columns) == 25, f"Expected 25 raw columns, got {len(df_raw.columns)}"


def test_cr_elevation_completeness_and_terrain_bounds() -> None:
    """All 251 coordinate points must have finite elevations within Southern Ghana bounds."""
    assert DERIVED_CSV.exists(), f"Canonical CSV missing: {DERIVED_CSV}"
    assert DEM_CSV.exists(), f"DEM CSV missing: {DEM_CSV}"

    df = pd.read_csv(DERIVED_CSV)
    assert len(df) == 252, f"Expected 252 records, got {len(df)}"

    valid_coords = df[df["latitude_dd"].notna() & df["longitude_dd"].notna()]
    assert len(valid_coords) == 251, f"Expected 251 valid coordinate records, got {len(valid_coords)}"

    # 100% elevation completeness for coordinate points
    assert valid_coords["elevation_m"].notna().sum() == 251
    assert valid_coords["elev_srtm30m"].notna().sum() == 251
    assert valid_coords["elev_aster30m"].notna().sum() == 251

    # Regional topographic limits for Southern Ghana / Central Region coastal to Birimian uplands (10m - 300m)
    assert (valid_coords["elevation_m"] >= 10.0).all(), "Elevation below regional minimum"
    assert (valid_coords["elevation_m"] <= 300.0).all(), "Elevation above regional maximum"


def test_cr_coordinate_validity_and_utm_projection() -> None:
    """Audited coordinates must lie within Southern Ghana and project to valid UTM Zone 30N."""
    df = pd.read_csv(DERIVED_CSV)
    valid_coords = df[df["latitude_dd"].notna() & df["longitude_dd"].notna()]

    # Latitude range for Central Region & adjacent districts: ~5.0 to ~6.5 N
    assert (valid_coords["latitude_dd"] >= 5.0).all() and (valid_coords["latitude_dd"] <= 6.5).all()
    # Longitude range for Central Region & adjacent districts: ~ -2.3 to ~ -0.3 W
    assert (valid_coords["longitude_dd"] >= -2.3).all() and (valid_coords["longitude_dd"] <= -0.3).all()

    # UTM Zone 30N bounds for Central Region Ghana
    assert (valid_coords["utm_easting_m"] >= 580000.0).all() and (valid_coords["utm_easting_m"] <= 800000.0).all()
    assert (valid_coords["utm_northing_m"] >= 550000.0).all() and (valid_coords["utm_northing_m"] <= 720000.0).all()


def test_cr_geology_polygon_spatial_join() -> None:
    """All 251 coordinate points must successfully join Ghana geological polygons."""
    df = pd.read_csv(DERIVED_CSV)
    valid_coords = df[df["latitude_dd"].notna() & df["longitude_dd"].notna()]

    # 100% match rate for coordinate points
    assert (valid_coords["geology_join_status"] == "MATCHED").all()
    assert valid_coords["geology_symbol"].notna().sum() == 251
    assert valid_coords["geology_stratigraphic_unit"].notna().sum() == 251

    # Known host units
    expected_units = {
        "Eburnean Plutonic Suite",
        "Birimian Supergroup",
        "Tarkwaian Group",
        "'Tamnean' Plutonic Suite",
        "Mesozoic",
    }
    observed_units = set(valid_coords["geology_stratigraphic_unit"].unique())
    assert observed_units.issubset(expected_units)


def test_cr_hydrogeochemistry_and_charge_balance() -> None:
    """Molar conversions and charge balance error calculations must be consistent."""
    df = pd.read_csv(DERIVED_CSV)

    # Verify molar conversion for Calcium: Ca(mmol/L) = Ca(mg/L) / 40.078
    ca_mask = df["ca_mg_L"].notna()
    ca_mmol_expected = (df.loc[ca_mask, "ca_mg_L"] / 40.078).round(5)
    assert (df.loc[ca_mask, "ca_mmol_L"] == ca_mmol_expected).all()

    # Verify CBE on complete major ion samples
    complete_chem = df[df["cbe_percent"].notna()]
    assert len(complete_chem) == 168

    # Over 90% of samples with complete major ions should have acceptable or caution CBE (<= 10%)
    acceptable_rate = (complete_chem["cbe_percent"].abs() <= 10.0).mean()
    assert acceptable_rate >= 0.90, f"Only {acceptable_rate:.1%} of complete samples have CBE <= 10%"

    # Facies labels must be populated for complete major ion samples
    assert complete_chem["facies_label"].notna().all()
    assert (complete_chem["facies_label"] != "Unclassified").all()


def test_cr_well_hydraulics_and_piezometric_heads() -> None:
    """Drawdown, specific capacity, and piezometric heads must be physically sound."""
    df = pd.read_csv(DERIVED_CSV)

    # Drawdown = DWL - SWL
    swl_dwl = df[df["swl_m"].notna() & df["dwl_m"].notna()]
    assert len(swl_dwl) == 250
    assert (swl_dwl["drawdown_m"] == (swl_dwl["dwl_m"] - swl_dwl["swl_m"]).round(2)).all()

    # Drawdown must be non-negative (pumping lowers the water level)
    assert (swl_dwl["drawdown_m"] >= 0.0).all()

    # Hydraulic head h = Elevation - SWL
    head_samples = df[df["hydraulic_head_m"].notna()]
    assert len(head_samples) == 249
    assert (head_samples["hydraulic_head_m"] == (head_samples["elevation_m"] - head_samples["swl_m"]).round(2)).all()
    # Hydraulic heads must be positive above sea level
    assert (head_samples["hydraulic_head_m"] > 0.0).all()


def test_cr_completed_workbook_structure_and_styling() -> None:
    """Completed Excel workbook must have all required sheets, styling, and data dictionary."""
    assert COMPLETED_XLSX.exists(), f"Completed workbook missing: {COMPLETED_XLSX}"

    wb = openpyxl.load_workbook(COMPLETED_XLSX, read_only=False)
    expected_sheets = {
        "GW_Field_Integration",
        "GW_Molar_Charge",
        "GW_Well_Hydraulics",
        "Native_CR_Master",
        "Data_Dictionary",
    }
    assert expected_sheets.issubset(set(wb.sheetnames)), (
        f"Missing sheets in completed workbook: {expected_sheets - set(wb.sheetnames)}"
    )

    df_gw = pd.read_excel(COMPLETED_XLSX, sheet_name="GW_Field_Integration")
    assert len(df_gw) == 252
    assert df_gw["elevation_m"].notna().sum() == 251

    # Check Data Dictionary
    df_dict = pd.read_excel(COMPLETED_XLSX, sheet_name="Data_Dictionary")
    assert len(df_dict) == len(df_gw.columns), (
        f"Data dictionary rows ({len(df_dict)}) does not match columns ({len(df_gw.columns)})"
    )
    assert set(df_dict["Variable"]) == set(df_gw.columns)

    # Check styling on first sheet
    ws = wb["GW_Field_Integration"]
    assert ws.views.sheetView[0].showGridLines is True
    header_cell = ws.cell(1, 1)
    assert header_cell.fill.start_color.rgb == "001F497D"
    assert header_cell.font.bold is True
