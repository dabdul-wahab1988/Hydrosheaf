import unittest

from hydrosheaf import Config, infer_edges


class SheafTopologyTests(unittest.TestCase):
    def test_sheaf_prefers_isotope_cluster(self):
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0,
                "2H": -30.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.1,
                "18O": -5.1,
                "2H": -30.5,
            },
            {
                "site_id": "C",
                "lat": 0.01,
                "lon": 0.0,
                "head_meas": 85.0,
                "Cl": 1.0,
                "18O": -10.0,
                "2H": -70.0,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.1
        config.sheaf_weight_head_prior = 0.2
        config.sheaf_weight_isotope = 1.0
        config.sheaf_weight_cl = 0.5

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->B", edge_ids)
        self.assertNotIn("A->C", edge_ids)

    def test_evaporation_gate_allows_enrichment(self):
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0,
                "2H": -30.0,
                "screen_depth": 5.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 2.0,
                "18O": -2.0,
                "2H": -15.0,
                "screen_depth": 3.0,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.1
        config.sheaf_weight_head_prior = 0.2
        config.sheaf_weight_isotope = 1.0
        config.sheaf_weight_cl = 0.5

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        self.assertEqual(len(edges), 1)
        pi_evap = edges[0].attrs.get("sheaf_pi_evap")
        self.assertIsNotNone(pi_evap)
        self.assertGreater(float(pi_evap), 0.2)

    def test_missing_isotopes_fallback(self):
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.2,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.1

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        self.assertEqual(len(edges), 1)
        flags = edges[0].attrs.get("sheaf_flags") or ""
        self.assertIn("iso_missing", flags)

    def test_depth_gating_limits_evaporation(self):
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0,
                "2H": -30.0,
                "screen_depth": 5.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 2.0,
                "18O": -2.0,
                "2H": -15.0,
                "screen_depth": 150.0,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.1

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        self.assertEqual(len(edges), 1)
        pi_evap = edges[0].attrs.get("sheaf_pi_evap")
        self.assertIsNotNone(pi_evap)
        self.assertLess(float(pi_evap), 0.4)


if __name__ == "__main__":
    unittest.main()
