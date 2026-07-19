"""Tests for dynamic mineral library."""

import unittest
from hydrosheaf.data.minerals import MINERAL_LIBRARY, get_mineral_stoich
from hydrosheaf.config import Config
from hydrosheaf.models.reactions import build_reaction_dictionary


class MineralLibraryTests(unittest.TestCase):
    def test_library_contents(self) -> None:
        """Verify standard minerals are present with correct keys."""
        self.assertIn("calcite", MINERAL_LIBRARY)
        self.assertIn("pyrite_oxidation_aerobic", MINERAL_LIBRARY)
        self.assertIn("enstatite", MINERAL_LIBRARY)
        self.assertIn("diopside", MINERAL_LIBRARY)
        self.assertIn("sulfate_reduction", MINERAL_LIBRARY)
        self.assertIn("iron_reduction", MINERAL_LIBRARY)
        # Check specific stoichiometry
        calcite = MINERAL_LIBRARY["calcite"]
        self.assertEqual(calcite["Ca"], 1)
        self.assertEqual(calcite["HCO3"], 2)

    def test_get_mineral_stoich(self) -> None:
        """Verify normalization lookup."""
        stoich = get_mineral_stoich("Calcite")
        self.assertEqual(stoich["Ca"], 1)

        with self.assertRaises(ValueError):
            get_mineral_stoich("NonExistentMineral")

    def test_dynamic_reaction_building(self) -> None:
        """Verify dictionary builder respects active_minerals config."""
        config = Config()
        config.active_minerals = ["calcite", "halite"]
        config.exchange_enabled = False  # disable for clarity

        matrix, labels, mineral_mask, *_ = build_reaction_dictionary(config)

        self.assertIn("calcite", labels)
        self.assertIn("halite", labels)
        self.assertNotIn("gypsum", labels)

        # Check standard reactions always present
        self.assertIn("NO3src", labels)
        self.assertIn("denit", labels)

    def test_pyrite_redox_stoichiometry(self) -> None:
        """Verify Pyrite redox modes."""
        aerobic = get_mineral_stoich("pyrite_oxidation_aerobic")
        self.assertEqual(aerobic["SO4"], 2)
        self.assertEqual(aerobic["Fe"], 1)

        denit = get_mineral_stoich("pyrite_oxidation_denit")
        self.assertEqual(denit["SO4"], 2)
        self.assertEqual(denit["NO3"], -2.8)
        self.assertEqual(denit["Fe"], 1)

    def test_reducing_processes_require_indicator_ions_and_redox_gate(self):
        reducing = Config(
            active_minerals=[],
            measured_ions=["SO4", "HCO3", "Fe", "NO3"],
            exchange_enabled=False,
        )
        _, labels, _, _ = build_reaction_dictionary(
            reducing,
            sample={"SO4": 0.4, "HCO3": 2.0, "Fe": 0.02, "NO3": 0.01},
        )
        self.assertIn("sulfate_reduction", labels)
        self.assertIn("iron_reduction", labels)

        oxic = Config(
            active_minerals=[],
            measured_ions=["SO4", "HCO3", "Fe", "NO3"],
            exchange_enabled=False,
        )
        _, oxic_labels, _, _ = build_reaction_dictionary(
            oxic,
            sample={"SO4": 0.4, "HCO3": 2.0, "Fe": 0.001, "NO3": 1.0},
        )
        self.assertNotIn("sulfate_reduction", oxic_labels)
        self.assertNotIn("iron_reduction", oxic_labels)

        missing_sulfate = Config(
            active_minerals=[],
            measured_ions=["HCO3", "Fe", "NO3"],
            exchange_enabled=False,
        )
        _, missing_labels, _, _ = build_reaction_dictionary(
            missing_sulfate,
            sample={"HCO3": 2.0, "Fe": 0.02, "NO3": 0.01},
        )
        self.assertNotIn("sulfate_reduction", missing_labels)


if __name__ == "__main__":
    unittest.main()
