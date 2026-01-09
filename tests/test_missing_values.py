import unittest

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import fit_network


class MissingValuesTests(unittest.TestCase):
    def test_missing_required_skips_edge(self):
        samples = [
            {
                "sample_id": "s1",
                "site_id": "A",
                "Ca": 1.0,
                "Mg": 1.0,
                "Na": 1.0,
                "HCO3": 1.0,
                "Cl": 1.0,
                "SO4": 1.0,
                "NO3": 1.0,
                "F": 1.0,
                "EC": 10.0,
                "TDS": 20.0,
                "pH": 7.0,
            },
            {
                "sample_id": "s2",
                "site_id": "B",
                "Ca": 1.0,
                "Mg": 1.0,
                "Na": 1.0,
                "HCO3": 1.0,
                "Cl": None,
                "SO4": 1.0,
                "NO3": 1.0,
                "F": 1.0,
                "EC": 11.0,
                "TDS": 22.0,
                "pH": 7.1,
            },
        ]
        results = fit_network(samples, [("A", "B")], Config())
        self.assertEqual(results, [])

    def test_missing_required_impute_zero(self):
        samples = [
            {
                "sample_id": "s1",
                "site_id": "A",
                "Ca": 1.0,
                "Mg": 1.0,
                "Na": 1.0,
                "HCO3": 1.0,
                "Cl": 1.0,
                "SO4": 1.0,
                "NO3": 1.0,
                "F": 1.0,
                "EC": 10.0,
                "TDS": 20.0,
                "pH": 7.0,
            },
            {
                "sample_id": "s2",
                "site_id": "B",
                "Ca": 1.0,
                "Mg": 1.0,
                "Na": 1.0,
                "HCO3": 1.0,
                "Cl": None,
                "SO4": 1.0,
                "NO3": 1.0,
                "F": 1.0,
                "EC": 11.0,
                "TDS": 22.0,
                "pH": 7.1,
            },
        ]
        results = fit_network(samples, [("A", "B")], Config(missing_policy="impute_zero"))
        self.assertEqual(len(results), 1)


if __name__ == "__main__":
    unittest.main()
