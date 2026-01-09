
import unittest
from hydrosheaf.graph3d.constraints import check_hydraulic_feasibility
from hydrosheaf.graph3d.types_3d import Node3D

class Test3DConstraints(unittest.TestCase):
    
    def test_hydraulic_head_check(self):
        """
        Test hydraulic feasibility logic.
        Flow should go from High Head -> Low Head.
        """
        # Uphill flow (Impossible)
        u_low = Node3D("U", x=0, y=0, z=0, elevation_m=100, hydraulic_head=90.0)
        v_high = Node3D("V", x=100, y=0, z=0, elevation_m=100, hydraulic_head=95.0)
        
        is_feasible, prob, grad = check_hydraulic_feasibility(u_low, v_high)
        self.assertFalse(is_feasible, "Flow from 90m to 95m head should be impossible")
        self.assertEqual(prob, 0.0)
        self.assertTrue(grad < 0)
        
        # Downhill flow (Possible)
        u_high = Node3D("U", x=0, y=0, z=0, elevation_m=100, hydraulic_head=95.0)
        v_low = Node3D("V", x=100, y=0, z=0, elevation_m=100, hydraulic_head=90.0)
        
        is_feasible, prob, grad = check_hydraulic_feasibility(u_high, v_low)
        self.assertTrue(is_feasible)
        self.assertEqual(prob, 1.0) # Strong gradient > default min
        self.assertTrue(grad > 0)
        
        # Missing data (Permissive) - With Bayesian fallback
        # U has high elevation (1000m), V low (100m). Prob should be ~1.0
        u_mtn = Node3D("U", x=0, y=0, z=0, elevation_m=1000, hydraulic_head=None)
        v_val = Node3D("V", x=100, y=0, z=0, elevation_m=100, hydraulic_head=None)
        
        is_feasible, prob, _ = check_hydraulic_feasibility(u_mtn, v_val)
        self.assertTrue(is_feasible)
        self.assertGreater(prob, 0.99)
        
        # Uncertain case: Equal elevation (100m vs 100m)
        # Prob should be 0.5 (could go either way)
        u_flat = Node3D("U", x=0, y=0, z=0, elevation_m=100, hydraulic_head=None)
        v_flat = Node3D("V", x=100, y=0, z=0, elevation_m=100, hydraulic_head=None)
        
        is_feasible, prob, _ = check_hydraulic_feasibility(u_flat, v_flat)
        self.assertTrue(is_feasible)
        self.assertAlmostEqual(prob, 0.5, delta=0.01)
        
        # Impossible case (Bayesian): U (100m) -> V (200m)
        # Flowing 100m uphill is very unlikely
        u_low_elv = Node3D("U", x=0, y=0, z=0, elevation_m=100, hydraulic_head=None)
        v_high_elv = Node3D("V", x=100, y=0, z=0, elevation_m=200, hydraulic_head=None)
        
        is_feasible, prob, _ = check_hydraulic_feasibility(u_low_elv, v_high_elv)
        # Should be rejected (P < 0.1)
        self.assertFalse(is_feasible)
        self.assertLess(prob, 0.01)
