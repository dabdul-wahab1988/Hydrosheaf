import unittest
import numpy as np
from hydrosheaf import Config
from hydrosheaf.sheaf.directed_section import solve_directed_section, DirectedEdgeMap
from hydrosheaf.graph.types import Edge

class NegativeConcentrationTests(unittest.TestCase):
    def test_negative_concentration_unconstrained(self):
        """
        Verify that the standard linear solver produces physically impossible 
        negative concentrations when reaction offsets (precipitation) are large.
        """
        # Graph: U -> V
        # U is fixed (obs), V is unknown.
        node_ids = ["U", "V"]
        
        # U has concentration 10.0
        # V has NO observation (None)
        node_obs = {
            "U": [10.0],
            "V": None 
        }
        
        # Edge U->V with alpha=1.0 (conservative transport)
        # But Offset = -20.0 (Strong precipitation/decay)
        # Weight = 1.0
        edge = Edge(u="U", v="V", edge_id="e1")
        edge_map = DirectedEdgeMap(
            edge=edge,
            alpha=1.0,
            offset=[-20.0], # Removes 20 units
            weight=1.0,
            objective=0.0,
            transport_model="evap",
            endmember_id=None,
            residual_norm=0.0
        )
        
        # Current solver (exact linear solve)
        # Equation: x_v - (1.0 * x_u + (-20)) = 0  => x_v = x_u - 20
        # If x_u is fixed at 10 (via strong obs weight), x_v should be -10.
        
        # Note: solve_directed_section adds obs_weight to diagonal.
        # For U: mat[0,0] += obs_weight. vec[0] += obs_weight * 10.
        # For V: No obs, so just the edge constraint.
        
        results = solve_directed_section(
            node_ids,
            [edge_map],
            node_obs,
            obs_weight=1000.0, # Strong observation at U
            non_negative=False # Explicitly disable constraint if added later
        )
        
        val_u = results["U"][0]
        val_v = results["V"][0]
        
        self.assertAlmostEqual(val_u, 10.0, places=2)
        # Expect negative result
        self.assertLess(val_v, 0.0)
        self.assertAlmostEqual(val_v, -10.0, places=2)

    def test_negative_concentration_constrained(self):
        """
        Verify that the new constrained solver prevents negative concentrations.
        """
        node_ids = ["U", "V"]
        node_obs = {"U": [10.0], "V": None}
        
        edge = Edge(u="U", v="V", edge_id="e1")
        edge_map = DirectedEdgeMap(
            edge=edge,
            alpha=1.0,
            offset=[-20.0],
            weight=1.0,
            objective=0.0,
            transport_model="evap",
            endmember_id=None,
            residual_norm=0.0
        )
        
        # This argument 'non_negative=True' is what we need to implement
        try:
            results = solve_directed_section(
                node_ids,
                [edge_map],
                node_obs,
                obs_weight=1000.0,
                non_negative=True
            )
        except TypeError:
            self.skipTest("non_negative parameter not yet implemented")
            
        val_v = results["V"][0]
        
        # Should be clamped to 0.0
        self.assertGreaterEqual(val_v, 0.0)
        self.assertLess(val_v, 0.1) # Close to zero

if __name__ == "__main__":
    unittest.main()
