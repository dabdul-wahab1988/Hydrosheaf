import unittest

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import edge_process_maps, fit_network, summarize_network


class NetworkSummaryTests(unittest.TestCase):
    def test_summarize_and_maps(self):
        samples = [
            {
                "sample_id": "s1",
                "site_id": "A",
                "Ca": 1.0,
                "Mg": 1.0,
                "Na": 0.0,
                "HCO3": 1.0,
                "Cl": 0.0,
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
                "Ca": 1.2,
                "Mg": 1.2,
                "Na": 0.0,
                "HCO3": 1.2,
                "Cl": 0.0,
                "SO4": 1.2,
                "NO3": 1.2,
                "F": 1.2,
                "EC": 12.0,
                "TDS": 24.0,
                "pH": 7.1,
            },
        ]
        results = fit_network(samples, [("A", "B")], Config(lambda_sparse=0.0, missing_policy="impute_zero"))
        summary = summarize_network(results)
        self.assertEqual(summary["edge_count"], 1)
        maps = edge_process_maps(results)
        self.assertEqual(len(maps["transport_likelihoods"]), 1)
        self.assertEqual(len(maps["reaction_intensity"]), 1)


if __name__ == "__main__":
    unittest.main()
