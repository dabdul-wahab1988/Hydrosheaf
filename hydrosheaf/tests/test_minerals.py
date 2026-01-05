"""Tests for dynamic mineral library."""

import unittest
from hydrosheaf.data.minerals import MINERAL_LIBRARY, get_mineral_stoich, Stoich
from hydrosheaf.config import Config
from hydrosheaf.models.reactions import build_reaction_dictionary

class MineralLibraryTests(unittest.TestCase):
    def test_library_contents(self) -> None:
        """Verify standard minerals are present with correct keys."""
        self.assertIn("calcite", MINERAL_LIBRARY)
        self.assertIn("pyrite_oxidation_aerobic", MINERAL_LIBRARY)
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
        config.exchange_enabled = False # disable for clarity
        
        matrix, labels, mineral_mask = build_reaction_dictionary(config)
        
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

if __name__ == "__main__":
    unittest.main()
