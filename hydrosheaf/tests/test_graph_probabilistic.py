import unittest

from hydrosheaf.graph.build import infer_edges_probabilistic


class GraphProbabilisticTests(unittest.TestCase):
    def test_probabilistic_edges(self):
        samples = [
            {"site_id": "A", "lat": 0.0, "lon": 0.0, "head_meas": 100.0},
            {"site_id": "B", "lat": 0.0, "lon": 0.01, "head_meas": 90.0},
            {"site_id": "C", "lat": 1.0, "lon": 1.0, "head_meas": 80.0},
        ]
        edges = infer_edges_probabilistic(
            samples,
            radius_km=5.0,
            max_neighbors=2,
            p_min=0.7,
            sigma_meas=0.5,
            sigma_dtw=1.0,
            sigma_elev=1.0,
            sigma_topo=10.0,
            gradient_min=0.0,
            depth_mismatch=20.0,
        )
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->B", edge_ids)
        self.assertNotIn("B->A", edge_ids)


if __name__ == "__main__":
    unittest.main()
