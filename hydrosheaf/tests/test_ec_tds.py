import unittest

from hydrosheaf.config import Config
from hydrosheaf.inference.network_fit import predict_node_ec_tds
from hydrosheaf.models.ec_tds import calibrate_ec_tds


class EcTdsTests(unittest.TestCase):
    def test_calibrate_ec_tds(self):
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
                "EC": 17.0,
                "TDS": 26.0,
                "pH": 7.0,
            },
            {
                "sample_id": "s2",
                "site_id": "B",
                "Ca": 2.0,
                "Mg": 2.0,
                "Na": 2.0,
                "HCO3": 2.0,
                "Cl": 2.0,
                "SO4": 2.0,
                "NO3": 2.0,
                "F": 2.0,
                "EC": 33.0,
                "TDS": 50.0,
                "pH": 7.2,
            },
        ]
        config = Config(missing_policy="impute_zero")
        calibrate_ec_tds(samples, config)
        self.assertAlmostEqual(config.ec_model[0], 2.0, places=6)
        self.assertAlmostEqual(config.ec_model[1], 1.0, places=6)
        self.assertAlmostEqual(config.tds_model[0], 3.0, places=6)
        self.assertAlmostEqual(config.tds_model[1], 2.0, places=6)

    def test_predict_node_ec_tds(self):
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
                "EC": 17.0,
                "TDS": 26.0,
                "pH": 7.0,
            }
        ]
        config = Config(ec_model=(2.0, 1.0), tds_model=(3.0, 2.0), missing_policy="impute_zero")
        rows = predict_node_ec_tds(samples, config)
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]["ec_pred"], 17.0, places=6)
        self.assertAlmostEqual(rows[0]["tds_pred"], 26.0, places=6)


if __name__ == "__main__":
    unittest.main()
