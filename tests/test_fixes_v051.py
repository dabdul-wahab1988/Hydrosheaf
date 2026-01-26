
import unittest
import math
from hydrosheaf.data.units import MOLAR_MASS_G_MOL
from hydrosheaf.data.schema import parse_numeric
from hydrosheaf.models.reactions import fit_reactions
from hydrosheaf.models.mixing import fit_evaporation, fit_mixing

class TestFixesV051(unittest.TestCase):
    def test_molar_weights_precision(self):
        """Verify high-precision molar weights."""
        self.assertAlmostEqual(MOLAR_MASS_G_MOL["Ca"], 40.078, places=3)
        self.assertAlmostEqual(MOLAR_MASS_G_MOL["SO4"], 96.064, places=3)
        self.assertAlmostEqual(MOLAR_MASS_G_MOL["Cl"], 35.453, places=3)

    def test_schema_negative_rejection(self):
        """Verify negative values are rejected."""
        self.assertIsNone(parse_numeric("-5.0", "half"))
        self.assertIsNone(parse_numeric("-1e-5", "value"))
        self.assertEqual(parse_numeric("5.0", "half"), 5.0)

    def test_reaction_fitting_stability(self):
        """Verify reaction fitting with singular matrix."""
        # Two identical reactions -> singular Gram matrix
        reaction_matrix = [
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0]
        ]
        residual = [2.0, 0.0, 0.0]
        weights = [1.0, 1.0, 1.0]
        
        # Should not raise error and should converge
        result = fit_reactions(residual, reaction_matrix, weights, lambda_l1=0.0)
        self.assertTrue(result.converged)
        # Check coefficients are finite
        for c in result.extents:
            self.assertTrue(math.isfinite(c))

    def test_evaporation_bounds(self):
        """Verify evaporation gamma upper bound."""
        x_u = [1.0, 1.0]
        x_v = [10000.0, 10000.0] # 10,000x concentration
        weights = [1.0, 1.0]
        
        gamma, _, _ = fit_evaporation(x_u, x_v, weights)
        self.assertLessEqual(gamma, 1000.0)
        self.assertGreaterEqual(gamma, 1.0)

    def test_mixing_stability(self):
        """Verify mixing with zero denominator."""
        x_u = [1.0]
        x_v = [1.0]
        x_end = [1.0] # Same as start, so delta is 0
        weights = [1.0]
        
        f, _, _ = fit_mixing(x_u, x_v, x_end, weights)
        self.assertEqual(f, 0.0)

if __name__ == '__main__':
    unittest.main()
