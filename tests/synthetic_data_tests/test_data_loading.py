import unittest
import pandas as pd
from pathlib import Path
from hydrosheaf.config import Config, DEFAULT_ION_ORDER
from hydrosheaf.data.units import mgL_to_mmolL

SYNTHETIC_DIR = Path(__file__).parents[2] / "hydrosheaf_synthetic_csv"

class DataLoadingTests(unittest.TestCase):
    def test_paths_exist(self):
        """Verify that the synthetic data directory and key files exist."""
        self.assertTrue(SYNTHETIC_DIR.exists(), f"Synthetic data dir not found: {SYNTHETIC_DIR}")
        files = [
            "water_chem_full.csv",
            "stations.csv",
            "network_edges.csv",
            "borehole_heads.csv"
        ]
        for f in files:
            p = SYNTHETIC_DIR / f
            self.assertTrue(p.exists(), f"File not found: {f}")

    def test_load_water_chem(self):
        """Validate water chemistry data structure and content."""
        df = pd.read_csv(SYNTHETIC_DIR / "water_chem_full.csv")
        
        # Check basic dimensions
        self.assertGreater(len(df), 0, "Water chemistry file is empty")
        
        # Check required columns (base)
        base_cols = ["event_code", "station_code", "collection_date"]
        for col in base_cols:
            self.assertIn(col, df.columns, f"Missing base column: {col}")
            
        # Check ions (allow for _mg_L suffix)
        # F, Fe, PO4 are missing in synthetic data, so we only enforce major ions
        major_ions = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3"]
        for ion in major_ions:
            # Check for ion OR ion_mg_L
            col = ion if ion in df.columns else f"{ion}_mg_L"
            self.assertIn(col, df.columns, f"Missing column for {ion}")
            
        # Check value ranges for pH
        self.assertTrue(df["pH"].between(3, 11).all(), "pH values out of realistic range")
        
        # Check for numeric concentrations
        for ion in DEFAULT_ION_ORDER:
            col = ion if ion in df.columns else f"{ion}_mg_L"
            if col in df.columns:
                self.assertTrue(pd.to_numeric(df[col], errors='coerce').notnull().all(), f"Non-numeric values in {col}")

    def test_load_stations(self):
        """Validate stations metadata."""
        df = pd.read_csv(SYNTHETIC_DIR / "stations.csv")
        self.assertEqual(len(df), 6, "Expected 6 stations")
        self.assertIn("station_code", df.columns)
        self.assertIn("station_type", df.columns)
        
        types = set(df["station_type"].unique())
        self.assertTrue({"lysimeter", "borehole"}.issubset(types))

    def test_unit_conversion(self):
        """Test conversion utility used for loading data."""
        # Test unit conversion function directly
        val_mg = 40.078 # Ca molar mass
        val_mmol = mgL_to_mmolL(val_mg, "Ca")
        self.assertAlmostEqual(val_mmol, 1.0, places=3)

if __name__ == "__main__":
    unittest.main()
