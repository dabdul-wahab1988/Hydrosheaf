import unittest
import numpy as np
from hydrosheaf import Config, infer_edges
from hydrosheaf.sheaf.topology_refine import (
    NodeIsotopeInfo,
    _edge_age_cost,
    refine_edges_with_sheaf,
)
from hydrosheaf.graph.types import Edge

class SheafStabilityTests(unittest.TestCase):
    @staticmethod
    def _age_node(node_id, age, sigma=2.0, identifiable=True):
        return NodeIsotopeInfo(
            node_id=node_id,
            d18o=None,
            d2h=None,
            d_excess=None,
            evap_index=None,
            cl=None,
            depth_m=None,
            p_evap=0.0,
            age_years=age,
            age_sigma_years=sigma,
            age_identifiable=identifiable,
        )

    def test_age_likelihood_ranks_matching_travel_increment(self):
        upstream = self._age_node("U", 10.0)
        downstream = self._age_node("V", 20.0)
        matching, _ = _edge_age_cost(
            upstream,
            downstream,
            expected_travel_years=10.0,
            travel_sigma_years=2.0,
        )
        mismatching, flags = _edge_age_cost(
            upstream,
            downstream,
            expected_travel_years=40.0,
            travel_sigma_years=2.0,
        )
        self.assertLess(matching, mismatching)
        self.assertIn("age_travel_mismatch", flags)

    def test_unidentified_age_is_neutral(self):
        upstream = self._age_node("U", 10.0, identifiable=False)
        downstream = self._age_node("V", 20.0)
        cost, flags = _edge_age_cost(upstream, downstream)
        self.assertEqual(cost, 0.0)
        self.assertIn("age_unidentified", flags)

    def test_age_direction_likelihood_penalizes_reversal(self):
        upstream = self._age_node("U", 10.0)
        older = self._age_node("V", 20.0)
        younger = self._age_node("W", 0.0)
        forward_cost, _ = _edge_age_cost(
            upstream,
            older,
            expected_travel_years=10.0,
            travel_sigma_years=5.0,
        )
        reverse_cost, flags = _edge_age_cost(
            upstream,
            younger,
            expected_travel_years=10.0,
            travel_sigma_years=5.0,
        )
        self.assertLess(forward_cost, reverse_cost)
        self.assertIn("age_reversal", flags)

    def test_bistable_network_convergence(self):
        """
        Test a scenario where two parallel paths are nearly identical.
        The algorithm should converge to a stable solution (either splitting weight
        or stably choosing one), rather than oscillating or diverging.
        """
        # A -> B, with two potential paths (Direct A->B via different candidate definitions
        # or effectively identical candidates).
        # To make it "bistable", we provide two edges that are both plausible but
        # slightly conflicting in terms of which one explains the data "better" depending on 
        # small fluctuations.
        
        # Actually, let's create a diamond graph: S -> (A, B) -> T.
        # S is source, T is target.
        # Path S->A->T and S->B->T.
        # If we only look at local fit, both might look good.
        
        samples = [
            # Source
            {"site_id": "S", "head_meas": 100.0, "Cl": 10.0, "18O": -5.0, "2H": -30.0, "lat": 0.0, "lon": 0.0},
            # Intermediate A
            {"site_id": "A", "head_meas": 90.0, "Cl": 12.0, "18O": -4.8, "2H": -29.0, "lat": 0.0, "lon": 0.01},
            # Intermediate B (Identical to A initially)
            {"site_id": "B", "head_meas": 90.0, "Cl": 12.0, "18O": -4.8, "2H": -29.0, "lat": 0.01, "lon": 0.0},
            # Target
            {"site_id": "T", "head_meas": 80.0, "Cl": 14.4, "18O": -4.6, "2H": -28.0, "lat": 0.01, "lon": 0.01},
        ]

        config = Config()
        config.edge_radius_km = 10.0
        config.edge_max_neighbors = 2  # Allow multiple parents if weights allow
        config.sheaf_max_iter = 10     # Enough iterations to see oscillation if present
        config.sheaf_weight_global = 2.0 # Strong feedback from global solve
        
        # We manually construct candidates to ensure the topology exists
        candidates = [
            Edge(u="S", v="A", edge_id="S->A"),
            Edge(u="S", v="B", edge_id="S->B"),
            Edge(u="A", v="T", edge_id="A->T"),
            Edge(u="B", v="T", edge_id="B->T"),
        ]

        # In a hard selection logic, if A->T and B->T are competing for "best parent of T"
        # and max_neighbors=1, it might flip flop if the global solve shifts the residuals slightly.
        # Here we test if it runs without error and produces a result. 
        # The key for "Soft Selection" is that we should see stable weights or a stable selection.
        
        selected_edges = refine_edges_with_sheaf(samples, candidates, config)
        
        # Check basic validity
        self.assertTrue(len(selected_edges) > 0)
        
        # With soft selection (to be implemented), we expect 'sheaf_weight' or similar
        # to be present and potentially < 1.0 for competing edges.
        edge_map = {e.edge_id: e for e in selected_edges}
        
        # Just ensure we didn't lose everything
        self.assertIn("S->A", edge_map)
        self.assertIn("S->B", edge_map)
        
    def test_soft_weighting_properties(self):
        """
        Verify that edges get weights between 0 and 1, and that
        bad edges get low weights but aren't necessarily hard-deleted immediately
        if we use a soft approach.
        """
        samples = [
            {"site_id": "U", "head_meas": 100.0, "Cl": 10.0, "lat": 0.0, "lon": 0.0},
            {"site_id": "V_good", "head_meas": 90.0, "Cl": 10.1, "lat": 0.0, "lon": 0.01}, # Compatible
            {"site_id": "V_bad", "head_meas": 90.0, "Cl": 500.0, "lat": 0.01, "lon": 0.0}, # Incompatible
        ]
        
        config = Config()
        # Enable soft selection mode (we will add this flag/logic)
        # config.sheaf_soft_selection = True # No longer needed, implicit
        
        candidates = [
            Edge(u="U", v="V_good", edge_id="Good"),
            Edge(u="U", v="V_bad", edge_id="Bad"),
        ]
        
        selected = refine_edges_with_sheaf(samples, candidates, config)
        edge_map = {e.edge_id: e for e in selected}
        
        if "Bad" in edge_map:
            # If the bad edge is kept, it should have a very low weight
            attrs = edge_map["Bad"].attrs
            # We assume we'll store the soft weight in 'sheaf_weight'
            weight = attrs.get("sheaf_weight", 1.0)
            self.assertLess(weight, 0.5, "Bad edge should have low weight in soft selection")
            
        if "Good" in edge_map:
             attrs = edge_map["Good"].attrs
             weight = attrs.get("sheaf_weight", 0.0)
             self.assertGreater(weight, 0.5, "Good edge should have high weight")

if __name__ == "__main__":
    unittest.main()
