"""Tests for CoDA SBP-ilr coordinate system."""

import math
import unittest
from hydrosheaf.coda_sbp import ilr_from_sbp, geometric_mean

class CodaSbpTests(unittest.TestCase):
    def test_geometric_mean(self):
        val = geometric_mean([1.0, 4.0])
        self.assertAlmostEqual(val, 2.0)
        
    def test_ilr_basic_validity(self):
        # 7 ions: Ca, Mg, Na, K, HCO3, Cl, SO4
        sample = {
            "Ca": 10.0, "Mg": 5.0, "Na": 100.0, "K": 20.0,
            "HCO3": 200.0, "Cl": 50.0, "SO4": 30.0
        }
        ilr, valid = ilr_from_sbp(sample)
        self.assertTrue(valid)
        self.assertEqual(len(ilr), 6)
        
    def test_ilr_missing_ion(self):
        sample = {"Ca": 10.0} # Missing others
        ilr, valid = ilr_from_sbp(sample)
        self.assertFalse(valid)
        self.assertIsNone(ilr)
        
    def test_ilr_zero_handling(self):
        sample = {
            "Ca": 0.0, "Mg": 5.0, "Na": 100.0, "K": 20.0,
            "HCO3": 200.0, "Cl": 50.0, "SO4": 30.0
        }
        # Should replace 0 with epsilon and return valid
        ilr, valid = ilr_from_sbp(sample)
        self.assertTrue(valid)
        # ilr1 = sqrt(4*3/(4+3)) * ln(g(Cat)/g(An))
        # 0.0 replaced by 1e-12 -> geometric mean will be very small
        # log will be very negative
        
    def test_ilr_scale_invariance(self):
        # Multiplying all concentrations by constant k should yield SAME ilr
        sample_a = {
            "Ca": 10.0, "Mg": 5.0, "Na": 100.0, "K": 20.0,
            "HCO3": 200.0, "Cl": 50.0, "SO4": 30.0
        }
        sample_b = {k: v * 5.5 for k, v in sample_a.items()}
        
        ilr_a, _ = ilr_from_sbp(sample_a)
        ilr_b, _ = ilr_from_sbp(sample_b)
        
        for a_val, b_val in zip(ilr_a, ilr_b):
            self.assertAlmostEqual(a_val, b_val, places=6)

if __name__ == "__main__":
    unittest.main()
