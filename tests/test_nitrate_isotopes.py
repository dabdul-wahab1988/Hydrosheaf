
import math
import unittest
from pathlib import Path

import pandas as pd
from hydrosheaf.nitrate_source_v2 import infer_node_posteriors, NitrateSourceResult
from hydrosheaf.models import nitrate_isotopes

class TestNitrateIsotopes(unittest.TestCase):
    def test_endmembers_loading(self):
        """Verify we can load the default JSON."""
        sources = nitrate_isotopes.load_isotope_endmembers()
        self.assertTrue(len(sources) >= 2)
        manure = next((s for s in sources if s.name == "Manure"), None)
        self.assertIsNotNone(manure)
        self.assertAlmostEqual(manure.d15N_mean, 15.0)

    def test_mixing_prob_manure(self):
        """Test a clear manure sample (High N15, Low O18)."""
        sources = nitrate_isotopes.load_isotope_endmembers()
        # Manure: d15N=15, d18O=5
        sample = nitrate_isotopes.IsotopeSample(d15N=15.0, d18O=5.0)
        probs = nitrate_isotopes.compute_isotope_prob(sample, sources)
        
        self.assertTrue(probs["Manure"] > 0.8, f"Manure should be dominant, got {probs}")

    def test_mixing_prob_fertilizer(self):
        """Test a clear fertilizer sample (Low N15, High O18)."""
        sources = nitrate_isotopes.load_isotope_endmembers()
        # Fertilizer (NO3): d15N=0, d18O=20
        sample = nitrate_isotopes.IsotopeSample(d15N=0.0, d18O=20.0)
        probs = nitrate_isotopes.compute_isotope_prob(sample, sources)
        
        self.assertTrue(probs["Fertilizer"] > 0.8, f"Fertilizer should be dominant, got {probs}")

    def test_integration_with_inference(self):
        """Test that infer_node_posteriors uses the isotope logic."""
        # 1. Create a mock dataframe with isotope columns
        df = pd.DataFrame([
            {
                "site_id": "Well_A", 
                "NO3": 50.0, 
                "Cl": 1.0, 
                "d15N": 15.0, 
                "d18O_NO3": 5.0,
                "d2H": -10.0, 
                "d18O": -3.0
            },
            {
                "site_id": "Well_B",
                "NO3": 50.0,
                "Cl": 1.0,
                # Missing isotopes -> Fallback
                "d2H": -10.0,
                "d18O": -3.0
            },
            # Add background samples to shift median
            {"site_id": "Bg_1", "NO3": 1.0, "Cl": 1.0, "d2H": -10, "d18O": -3},
            {"site_id": "Bg_2", "NO3": 2.0, "Cl": 1.0, "d2H": -10, "d18O": -3},
            {"site_id": "Bg_3", "NO3": 1.5, "Cl": 1.0, "d2H": -10, "d18O": -3},
        ])
        df.set_index("site_id", inplace=True)
        
        # 2. Run inference
        results = infer_node_posteriors(df, edge_results=[])
        
        # 3. Check Well_A (Isotope)
        res_a = results["Well_A"]
        self.assertEqual(res_a.reason_code, "Dual Isotope Mixing")
        self.assertGreater(res_a.p_manure, 0.8)
        self.assertIn("dual_isotope_priority", res_a.gating_flags)
        
        # 4. Check Well_B (Fallback)
        res_b = results["Well_B"]
        self.assertEqual(res_b.reason_code, "Hydrochemical Ratios (No Isotopes)")
        # With High NO3/Cl (50) vs Median (~1.5), Z-score > 0
        # High Ratio => Fertilizer => Low p_manure
        self.assertLess(res_b.p_manure, 0.5)

if __name__ == "__main__":
    unittest.main()
