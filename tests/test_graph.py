import unittest

from hydrosheaf.graph.build import infer_edges_from_coordinates


class GraphInferTests(unittest.TestCase):
    def test_infer_edges_downhill(self):
        samples = [
            {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0},
            {"site_id": "B", "lat": 0.0, "lon": 0.5, "elevation": 5.0},
            {"site_id": "C", "lat": 1.0, "lon": 0.0, "elevation": 1.0},
        ]
        edges = infer_edges_from_coordinates(samples, max_neighbors=1, allow_uphill=False)
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->B", edge_ids)
        self.assertIn("B->C", edge_ids)
        self.assertNotIn("C->A", edge_ids)

    def test_infer_edges_flow_to(self):
        samples = [
            {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0, "flow_to": "C"},
            {"site_id": "B", "lat": 0.0, "lon": 0.5, "elevation": 5.0},
            {"site_id": "C", "lat": 1.0, "lon": 0.0, "elevation": 1.0},
        ]
        edges = infer_edges_from_coordinates(samples, max_neighbors=1, allow_uphill=False)
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->C", edge_ids)

    def test_infer_edges_head_key(self):
        samples = [
            {"site_id": "A", "lat": 0.0, "lon": 0.0, "hydraulic_head": 10.0},
            {"site_id": "B", "lat": 0.0, "lon": 0.5, "hydraulic_head": 5.0},
        ]
        edges = infer_edges_from_coordinates(samples, max_neighbors=1, allow_uphill=False)
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->B", edge_ids)


if __name__ == "__main__":
    unittest.main()
