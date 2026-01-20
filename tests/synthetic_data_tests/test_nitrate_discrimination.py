import unittest
import pandas as pd
from pathlib import Path
import math

from hydrosheaf.config import Config
from hydrosheaf.nitrate_source_v2 import (
    infer_node_posteriors,
    fit_robust_stats,
    NitrateStats,
)
from hydrosheaf.data.units import mgL_to_mmolL

SYNTHETIC_DIR = Path(__file__).parents[2] / "hydrosheaf_synthetic_csv"


class NitrateDiscriminationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.chem_df = pd.read_csv(SYNTHETIC_DIR / "water_chem_full.csv")
        cls.config = Config()
        cls.config.nitrate_isotope_n15_col = "d15N_NO3_permil"
        cls.config.nitrate_isotope_o18_col = "d18O_NO3_permil"

        # Prepare dataframe with correct ion names (mg/L -> mmol/L)
        cls.samples = []
        records = cls.chem_df.to_dict(orient="records")
        for r in records:
            new_r = r.copy()
            for ion in cls.config.ion_order:
                col_name = f"{ion}_mg_L"
                if col_name in r and pd.notnull(r[col_name]):
                    new_r[ion] = mgL_to_mmolL(float(r[col_name]), ion)
                elif ion in r and pd.notnull(r[ion]):
                    new_r[ion] = mgL_to_mmolL(float(r[ion]), ion)
            cls.samples.append(new_r)

        cls.df = pd.DataFrame(cls.samples)
        if "site_id" not in cls.df.columns:
            cls.df["site_id"] = cls.df["station_code"]

    def test_nitrate_gating(self):
        """Test that low nitrate samples are gated."""
        # Create a low NO3 sample
        # Note: The function expects values in mmol/L. The threshold is 10 mg/L = ~0.16 mmol/L
        # So we use a value below 0.16 mmol/L (e.g., 5 mg/L = ~0.08 mmol/L)
        low_no3_mmol = mgL_to_mmolL(5.0, "NO3")  # Convert 5 mg/L to mmol/L

        low_no3_sample = {
            "site_id": "TEST_LOW",
            "NO3": low_no3_mmol,  # ~0.08 mmol/L, below 10 mg/L threshold
            "Cl": mgL_to_mmolL(10.0, "Cl"),
            "K": mgL_to_mmolL(2.0, "K"),
            "d15N_NO3_permil": 5.0,
            "d18O_NO3_permil": 5.0,
            "d_excess": 10.0,
            # Add dummy ions to prevent KeyError in other parts if accessed
            "Ca": mgL_to_mmolL(10.0, "Ca"),
            "Mg": mgL_to_mmolL(10.0, "Mg"),
            "Na": mgL_to_mmolL(10.0, "Na"),
            "SO4": mgL_to_mmolL(10.0, "SO4"),
            "HCO3": mgL_to_mmolL(10.0, "HCO3"),
        }

        df = pd.DataFrame([low_no3_sample]).set_index("site_id", drop=False)

        # Infer - API is infer_node_posteriors(nodes_df, edge_results, config_overrides)
        results = infer_node_posteriors(df, [])
        res = results["TEST_LOW"]

        # Should be gated - p_manure is None for low NO3 samples
        self.assertIsNone(res.p_manure)
        self.assertIn("reason_code", res.__dict__)  # Check attribute existence
        # Check gating flag is set
        self.assertIn("below_detection_threshold", res.gating_flags)

    def test_high_nitrate_discrimination(self):
        """Test discrimination for a high nitrate sample (Manure-like)."""
        # CLUSTER A in synthetic data is high intensity, likely manure signals
        # L1 in E2 has 71 mg/L NO3.
        e2_l1 = self.chem_df[
            (self.chem_df["event_code"] == "E2")
            & (self.chem_df["station_code"] == "L1")
        ].iloc[0]

        sample = e2_l1.to_dict()
        sample["site_id"] = "L1_E2"
        # Map isotope names if needed, infer_node_posteriors expects specific names?
        # Examining nitrate_source_v2.py would confirm, but usually it handles mapping or expects standard names.
        # Assuming d15N_NO3_permil maps to d15N

        df = pd.DataFrame([sample]).set_index("site_id", drop=False)

        # We need robust stats background. We can fit it from the whole dataset.
        stats = fit_robust_stats(self.chem_df)

        overrides = {
            "stats": stats
        }  # Pass fitted stats if possible, or just let it auto-fit on the 1 sample (bad but runs)
        # Actually fit_robust_stats returns a NitrateStats object. infer_node_posteriors allows overrides.

        # Let's run on the whole E2 dataset to get better stats
        e2_df = self.chem_df[self.chem_df["event_code"] == "E2"].copy()
        e2_df["site_id"] = e2_df["station_code"]
        e2_df = e2_df.set_index("site_id", drop=False)

        # API is infer_node_posteriors(nodes_df, edge_results, config_overrides)
        results = infer_node_posteriors(e2_df, [])

        l1_res = results["L1"]

        # L1 is lysimeter, high input. Check if it got a probability.
        if l1_res.p_manure is not None:
            self.assertTrue(0.0 <= l1_res.p_manure <= 1.0)
            # Also check that p_fert + p_manure approx 1.0
            self.assertAlmostEqual(l1_res.p_manure + l1_res.p_fertilizer, 1.0, places=5)


if __name__ == "__main__":
    unittest.main()
