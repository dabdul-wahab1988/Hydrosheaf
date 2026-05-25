"""Unit tests for null-model package (Phase 0-1 assumption calibration)."""

import unittest

from hydrosheaf.config import Config
from hydrosheaf.null_models.chemistry import chemistry_null_score
from hydrosheaf.null_models.lithology import lithology_null_score
from hydrosheaf.null_models.endmembers import endmember_null_score
from hydrosheaf.null_models import compute_null_penalty


class ChemistryNullModelTests(unittest.TestCase):
    def test_chemistry_similar_scores_high(self):
        """Two wells with near-identical major-ion chemistry should get high null score."""
        config = Config()
        sample_a = {"Ca": 2.0, "Mg": 1.0, "Na": 5.0, "Cl": 3.0, "SO4": 0.5}
        sample_b = {"Ca": 2.1, "Mg": 0.9, "Na": 5.1, "Cl": 3.1, "SO4": 0.55}
        score, flags = chemistry_null_score(sample_a, sample_b, config)
        # Very similar -> high null score
        self.assertGreater(score, 0.5, f"Expected high null score, got {score}")
        self.assertIn("null_chemistry_similar", flags)

    def test_different_chemistry_scores_low(self):
        """Two wells with very different chemistry should get low null score."""
        config = Config()
        sample_a = {"Ca": 20.0, "Mg": 10.0, "Na": 50.0, "Cl": 3.0, "SO4": 0.5}
        sample_b = {"Ca": 0.5, "Mg": 0.1, "Na": 2.0, "Cl": 30.0, "SO4": 15.0}
        score, flags = chemistry_null_score(sample_a, sample_b, config)
        self.assertLess(score, 0.5, f"Expected low null score, got {score}")

    def test_missing_ions_handled(self):
        """Partially missing ion data should not crash."""
        config = Config()
        sample_a = {"Ca": 2.0}
        sample_b = {"Ca": 2.1, "Mg": 0.9, "Na": 5.1}
        score, flags = chemistry_null_score(sample_a, sample_b, config)
        # Should not crash; Ca only is similar -> some null score
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    def test_empty_samples_zero_score(self):
        """Completely empty samples should return zero null score."""
        config = Config()
        score, flags = chemistry_null_score({}, {}, config)
        self.assertEqual(score, 0.0)


class LithologyNullModelTests(unittest.TestCase):
    def test_common_lithology_raises_null(self):
        """Same aquifer layer should increase null score."""
        config = Config()
        sample_a = {"aquifer_layer": "Upper Aquifer"}
        sample_b = {"aquifer_layer": "Upper Aquifer"}
        score, flags = lithology_null_score(sample_a, sample_b, config)
        self.assertGreater(score, 0.0)
        self.assertIn("null_common_lithology", flags)

    def test_different_lithology_no_null(self):
        """Different aquifer layers should not raise null flag."""
        config = Config()
        sample_a = {"aquifer_layer": "Upper Aquifer"}
        sample_b = {"aquifer_layer": "Bedrock"}
        score, flags = lithology_null_score(sample_a, sample_b, config)
        self.assertEqual(score, 0.0)
        self.assertNotIn("null_common_lithology", flags)

    def test_missing_lithology_zero(self):
        """No lithology info should return zero."""
        config = Config()
        score, flags = lithology_null_score({}, {}, config)
        self.assertEqual(score, 0.0)


class EndmemberNullModelTests(unittest.TestCase):
    def test_similar_isotopes_yields_null_score(self):
        """Wells with nearly identical isotopes should get null flag."""
        config = Config()
        config.isotope_d18o_key = "18O"
        config.isotope_d2h_key = "2H"
        sample_a = {"18O": -5.0, "2H": -30.0}
        sample_b = {"18O": -5.0, "2H": -30.0}
        score, flags = endmember_null_score(sample_a, sample_b, config)
        # Identical isotopes -> very plausible shared recharge
        self.assertGreater(score, 0.0)
        self.assertIn("null_shared_recharge", flags)

    def test_different_isotopes_low_null(self):
        """Very different isotopes should not suggest shared recharge."""
        config = Config()
        config.isotope_d18o_key = "18O"
        config.isotope_d2h_key = "2H"
        sample_a = {"18O": -2.0, "2H": -10.0}
        sample_b = {"18O": -12.0, "2H": -80.0}
        score, flags = endmember_null_score(sample_a, sample_b, config)
        self.assertLess(score, 0.5)

    def test_missing_isotopes_zero(self):
        """No isotope data should return zero endmember null."""
        config = Config()
        score, flags = endmember_null_score({}, {}, config)
        self.assertEqual(score, 0.0)

    def test_anthropogenic_flag(self):
        """Both wells with high nitrate should flag common anthropogenic source."""
        config = Config()
        sample_a = {"NO3": 50.0, "Cl": 20.0}
        sample_b = {"NO3": 45.0, "Cl": 22.0}
        score, flags = endmember_null_score(sample_a, sample_b, config)
        self.assertIn("null_common_anthropogenic", flags)


class CombinedNullPenaltyTests(unittest.TestCase):
    def test_combined_returns_bounded_score(self):
        """Combined null penalty should be in [0, 1]."""
        config = Config()
        config.null_model_enabled = True
        sample_a = {
            "Ca": 2.0, "Mg": 1.0, "Na": 5.0, "Cl": 3.0,
            "18O": -5.0, "2H": -30.0,
            "aquifer_layer": "Upper Aquifer",
        }
        sample_b = {
            "Ca": 2.1, "Mg": 0.9, "Na": 5.1, "Cl": 3.1,
            "18O": -5.0, "2H": -30.0,
            "aquifer_layer": "Upper Aquifer",
        }
        score, flags = compute_null_penalty(sample_a, sample_b, config)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    def test_combined_no_crash_on_empty(self):
        """Combined penalty should not crash on empty samples."""
        config = Config()
        score, flags = compute_null_penalty({}, {}, config)
        self.assertEqual(score, 0.0)


if __name__ == "__main__":
    unittest.main()
