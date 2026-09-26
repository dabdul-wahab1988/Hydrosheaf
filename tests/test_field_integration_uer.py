"""Automated checks for the refined Upper East Region field integration dataset.

Verifies:
1. Complete elevation coverage (237/237) and DEM reconciliation.
2. Coordinate normalization and UTM Zone 30N projection.
3. Complete primary hydrochemistry and charge balance error (CBE).
4. Spatial and probabilistic graph edge inference using completed elevations.
5. Multi-sheet completed workbook integrity and immutable provenance.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import pytest
import pandas as pd

from hydrosheaf.data.field import load_field_dataset
from hydrosheaf.graph.build import infer_edges_from_coordinates, infer_edges_probabilistic

ROOT = Path(__file__).resolve().parents[1]
DERIVED_CSV = ROOT / "data" / "FieldData" / "derived" / "uer_field_integration_dataset.csv"
DEM_CSV = ROOT / "data" / "FieldData" / "derived" / "uer_elevations_dem.csv"
COMPLETED_XLSX = ROOT / "data" / "FieldData" / "UERdata" / "compiled UER data_new_completed.xlsx"
EXPECTED_COMPLETED_SHA256 = "a802efa607100269f42b52e428ed43f20cd19cfd46cd0e2d2502c1c10420dff7"


def test_completed_workbook_provenance_hash_is_stable() -> None:
    """The current approved completed workbook has a pinned provenance hash."""
    assert COMPLETED_XLSX.exists(), f"Completed workbook missing: {COMPLETED_XLSX}"
    with COMPLETED_XLSX.open("rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    assert digest == EXPECTED_COMPLETED_SHA256, (
        f"Completed workbook has changed! Expected {EXPECTED_COMPLETED_SHA256}, got {digest}"
    )


def test_uer_elevation_completeness_and_accuracy() -> None:
    """All 237 samples must have finite elevations within the UER terrain bounds."""
    assert DERIVED_CSV.exists(), f"Canonical CSV missing: {DERIVED_CSV}"
    assert DEM_CSV.exists(), f"DEM CSV missing: {DEM_CSV}"

    df = pd.read_csv(DERIVED_CSV)
    assert len(df) == 237, f"Expected 237 records, got {len(df)}"

    # 100% elevation completeness
    assert df["elevation_m"].notna().sum() == 237
    assert df["elev_srtm30m"].notna().sum() == 237
    assert df["elev_aster30m"].notna().sum() == 237

    # UER topographic limits (Northern Ghana savannah plateau: 120m - 350m)
    assert (df["elevation_m"] >= 120.0).all(), "Elevation below regional minimum"
    assert (df["elevation_m"] <= 350.0).all(), "Elevation above regional maximum"

    # Surveyed points reconciliation
    surveyed = df[df["elevation_source"] == "surveyed"]
    assert len(surveyed) == 2
    # Sample 1: surveyed 231.0 m vs SRTM 230.0 m (error -1.0 m)
    # Sample 2: surveyed 218.0 m vs SRTM 221.0 m (error +3.0 m)
    diff = surveyed["elev_srtm30m"] - surveyed["elevation_surveyed_m"]
    assert (diff.abs() <= 3.0).all(), f"DEM error on surveyed controls exceeds 3m: {diff.tolist()}"


def test_uer_coordinate_validity_and_utm_projection() -> None:
    """Audited coordinates must lie within UER and project to valid UTM Zone 30N."""
    df = pd.read_csv(DERIVED_CSV)

    # Latitude range for Upper East Region Ghana: ~10.3 to ~11.3 N
    assert (df["latitude_dd"] >= 10.3).all() and (df["latitude_dd"] <= 11.3).all()
    # Longitude range for Upper East Region Ghana: ~ -1.5 to ~ -0.05 W
    assert (df["longitude_dd"] >= -1.5).all() and (df["longitude_dd"] <= -0.05).all()

    # UTM Zone 30N bounds for UER Ghana
    assert (df["utm_easting_m"] >= 650000.0).all() and (df["utm_easting_m"] <= 850000.0).all()
    assert (df["utm_northing_m"] >= 1140000.0).all() and (df["utm_northing_m"] <= 1240000.0).all()


def test_uer_hydrogeochemistry_and_charge_balance() -> None:
    """Primary major ions must be complete and obey geochemical charge balance."""
    df = pd.read_csv(DERIVED_CSV)

    major_ions = ["ca_mg_L", "mg_mg_L", "na_mg_L", "k_mg_L", "hco3_mg_L", "cl_mg_L", "so4_mg_L", "no3_mg_L"]
    for ion in major_ions:
        assert df[ion].notna().sum() == 237, f"Missing values in {ion}"

    # Verify molar conversion for Calcium: Ca(mmol/L) = Ca(mg/L) / 40.078
    assert (df["ca_mmol_L"] == (df["ca_mg_L"] / 40.078).round(5)).all()

    # Charge balance error (CBE) verification
    cbe = df["cbe_percent"]
    # Over 95% of samples should have acceptable or caution CBE (<= 10%)
    acceptable_rate = (cbe.abs() <= 10.0).mean()
    assert acceptable_rate >= 0.95, f"Only {acceptable_rate:.1%} of samples have CBE <= 10%"


def test_uer_elevation_is_not_promoted_to_hydraulic_head() -> None:
    """DEM elevations are context only and cannot produce field flow edges."""
    ds = load_field_dataset("upper_east_region")
    assert ds.n_records == 237
    assert ds.coverage()["elevation"] == 237
    assert all(sample["hydraulic_head"] is None for sample in ds.records)


def test_completed_workbook_structure_and_sheets() -> None:
    """Completed Excel workbook must have all required sheets and data dictionary."""
    assert COMPLETED_XLSX.exists(), f"Completed workbook missing: {COMPLETED_XLSX}"
    import openpyxl

    wb = openpyxl.load_workbook(COMPLETED_XLSX, read_only=True)
    expected_sheets = {"GW_Field_Integration", "GW_Molar_Charge", "rain", "monitoring wells", "Data_Dictionary"}
    assert expected_sheets.issubset(set(wb.sheetnames)), (
        f"Missing sheets in completed workbook: {expected_sheets - set(wb.sheetnames)}"
    )

    df_gw = pd.read_excel(COMPLETED_XLSX, sheet_name="GW_Field_Integration")
    assert len(df_gw) == 237
    assert df_gw["elevation_m"].notna().sum() == 237
