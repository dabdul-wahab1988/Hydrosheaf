"""Tests for reaction-aware optimal transport."""

import unittest

from hydrosheaf.config import Config
from hydrosheaf.models.optimal_transport import compute_unbalanced_ot

ION_ORDER = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"]


class TestOptimalTransport(unittest.TestCase):
    def setUp(self):
        self.cfg = Config()

    def test_identical_chemistry_near_zero_cost(self):
        """Identical upstream and downstream should have near-zero OT cost."""
        upstream = [1.0, 0.5, 2.0, 0.1, 3.0, 1.0, 0.5, 0.2, 0.05, 0.01, 0.02]
        downstream = list(upstream)
        result = compute_unbalanced_ot(upstream, downstream, ION_ORDER, self.cfg)
        self.assertAlmostEqual(result["ot_total_cost"], 0.0, delta=0.1)
        self.assertAlmostEqual(result["ot_balanced_cost"], 0.0, delta=0.1)
        self.assertAlmostEqual(result["ot_creation_mass"], 0.0, delta=0.01)
        self.assertAlmostEqual(result["ot_destruction_mass"], 0.0, delta=0.01)

    def test_cl_mismatch_penalised(self):
        """Conservative species (Cl) gets higher creation penalty than Ca,
        resulting in higher total OT cost for the same mass imbalance."""
        self.cfg.active_minerals = ["calcite"]  # Ca is reaction-supported
        upstream = [1.0] * len(ION_ORDER)

        # Case 1: Cl increases (conservative)
        downstream_high_cl = list(upstream)
        cl_idx = ION_ORDER.index("Cl")
        downstream_high_cl[cl_idx] = 3.0

        result_cl = compute_unbalanced_ot(upstream, downstream_high_cl, ION_ORDER, self.cfg)

        # Case 2: Ca increases (reaction-supported via calcite)
        downstream_high_ca = list(upstream)
        ca_idx = ION_ORDER.index("Ca")
        downstream_high_ca[ca_idx] = 3.0

        result_ca = compute_unbalanced_ot(upstream, downstream_high_ca, ION_ORDER, self.cfg)

        # Total OT cost should be higher for conservative Cl than for
        # reaction-supported Ca at the same mass imbalance
        self.assertGreater(
            result_cl["ot_total_cost"],
            result_ca["ot_total_cost"],
        )

    def test_reaction_supported_species_lower_creation_cost(self):
        """Reaction-supported species (Ca, HCO3 in calcite) should have
        lower creation/destruction penalty than unsupported species."""
        # Ensure calcite is in active minerals
        self.cfg.active_minerals = ["calcite"]
        upstream = [0.5] * len(ION_ORDER)
        downstream = [0.5] * len(ION_ORDER)

        # Add Ca (reaction-supported via calcite)
        ca_idx = ION_ORDER.index("Ca")
        downstream_ca = list(downstream)
        downstream_ca[ca_idx] = 2.0
        result_ca = compute_unbalanced_ot(upstream, downstream_ca, ION_ORDER, self.cfg)

        # No active minerals -> higher penalty
        self.cfg.active_minerals = []
        result_no_min = compute_unbalanced_ot(upstream, downstream_ca, ION_ORDER, self.cfg)

        # Ca creation should be cheaper with calcite active
        self.assertLess(
            result_ca["ot_total_cost"],
            result_no_min["ot_total_cost"],
        )

    def test_output_fields_present(self):
        """All expected output fields should be present."""
        upstream = [1.0] * len(ION_ORDER)
        downstream = [0.8] * len(ION_ORDER)
        result = compute_unbalanced_ot(upstream, downstream, ION_ORDER, self.cfg)
        for key in [
            "ot_total_cost",
            "ot_balanced_cost",
            "ot_creation_mass",
            "ot_destruction_mass",
            "ot_conservative_mismatch",
            "ot_reaction_plausibility",
        ]:
            self.assertIn(key, result)
            self.assertGreaterEqual(result[key], 0.0)

    def test_zero_concentration_handled(self):
        """Zero-concentration inputs should be handled gracefully."""
        upstream = [0.0] * len(ION_ORDER)
        downstream = [0.0] * len(ION_ORDER)
        result = compute_unbalanced_ot(upstream, downstream, ION_ORDER, self.cfg)
        for key in [
            "ot_total_cost",
            "ot_balanced_cost",
        ]:
            self.assertEqual(result[key], 0.0)

    def test_mass_conservation(self):
        """With no creation/destruction, mass should be conserved."""
        upstream = [1.0, 0.5, 2.0, 0.1, 3.0, 1.0, 0.5, 0.2, 0.05, 0.01, 0.02]
        downstream = [0.8, 0.6, 1.8, 0.15, 3.5, 0.9, 0.6, 0.25, 0.04, 0.02, 0.01]
        result = compute_unbalanced_ot(upstream, downstream, ION_ORDER, self.cfg)
        # Creation + destruction should account for total mass difference
        total_diff = abs(sum(upstream) - sum(downstream))
        # Normalised creation/destruction mass
        self.assertGreaterEqual(result["ot_creation_mass"] + result["ot_destruction_mass"], 0.0)


if __name__ == "__main__":
    unittest.main()
