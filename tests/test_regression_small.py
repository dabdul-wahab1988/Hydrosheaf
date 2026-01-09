import json
import pathlib
import unittest

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import fit_network


class RegressionSmallNetworkTests(unittest.TestCase):
    def test_small_network_fixture(self):
        fixtures_dir = pathlib.Path(__file__).parent / "fixtures"
        samples = json.loads((fixtures_dir / "small_network_samples.json").read_text())
        expected = json.loads((fixtures_dir / "small_network_expected.json").read_text())

        config = Config(lambda_sparse=0.0, missing_policy="impute_zero")
        results = fit_network(samples, [("A", "B")], config)

        self.assertEqual(len(results), 1)
        result = results[0]

        self.assertEqual(result.edge_id, expected["edge_id"])
        self.assertEqual(result.transport_model, expected["transport_model"])
        self.assertAlmostEqual(result.gamma, expected["gamma"], places=6)
        self.assertAlmostEqual(
            result.transport_residual_norm,
            expected["transport_residual_norm"],
            places=6,
        )
        self.assertAlmostEqual(result.anomaly_norm, expected["anomaly_norm"], places=6)
        self.assertAlmostEqual(result.objective_score, expected["objective_score"], places=6)

        z_map = dict(zip(result.z_labels, result.z_extents))
        for label, value in expected["z"].items():
            self.assertAlmostEqual(z_map[label], value, places=6)
        for label, value in z_map.items():
            if label not in expected["z"]:
                self.assertAlmostEqual(value, 0.0, places=6)


if __name__ == "__main__":
    unittest.main()
