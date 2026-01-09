
import math
import unittest
from hydrosheaf.graph3d.distance import (
    haversine_distance, 
    compute_3d_distance, 
    compute_screen_overlap,
    classify_edge_type
)
from hydrosheaf.graph3d.layers import (
    get_aquitard_probability,
    compute_layer_probability
)
from hydrosheaf.graph3d.types_3d import Node3D, LayeredAquiferSystem

class Test3DFlowNetworksAccuracy(unittest.TestCase):
    """
    Test accuracy of 3D Flow Networks extension from a mathematical hydrogeology perspective.
    """

    def test_haversine_distance(self):
        """
        Verify Haversine distance calculation against known values.
        """
        # Point 1: Greenwich (51.4779 N, 0.0000 E)
        # Point 2: Eiffel Tower (48.8584 N, 2.2945 E)
        # Approximate distance: 343 km
        
        # Test 1 degree longitude at Equator
        # Should be exactly (pi * R) / 180
        # 6371 * pi / 180 = 111.19 km
        dist_eq = haversine_distance(0, 0, 0, 1)
        self.assertAlmostEqual(dist_eq, 111.195, delta=0.1)
        
        # Test 1 degree latitude
        dist_lat = haversine_distance(0, 0, 1, 0)
        self.assertAlmostEqual(dist_lat, 111.195, delta=0.1)
        
        # Zero distance
        self.assertAlmostEqual(haversine_distance(10, 10, 10, 10), 0.0)

    def test_anisotropic_distance(self):
        """
        Verify anisotropic 3D distance calculation.
        d_3d = sqrt(d_xy^2 + (d_z / alpha)^2)
        """
        # Setup nodes using Euclidean coordinates for simplicity
        node_a = Node3D("A", x=0, y=0, z=0, elevation_m=100)
        node_b = Node3D("B", x=300, y=400, z=100, elevation_m=100)
        
        # d_xy = sqrt(300^2 + 400^2) = 500 m
        # d_z = 100 m
        
        # Case 1: Isotropic (alpha = 1.0)
        d_3d_iso = math.sqrt(500**2 + 100**2) # sqrt(250000 + 10000) = sqrt(260000) ~= 509.9
        
        _, _, calc_iso = compute_3d_distance(node_a, node_b, anisotropy_factor=1.0, use_haversine=False)
        self.assertAlmostEqual(calc_iso, d_3d_iso)
        
        # Case 2: Anisotropic (alpha = 0.1) -> Vertical distance magnified 10x
        # d_z_scaled = 100 / 0.1 = 1000
        d_3d_aniso = math.sqrt(500**2 + 1000**2) # sqrt(250000 + 1000000) = sqrt(1250000) ~= 1118.0
        
        _, _, calc_aniso = compute_3d_distance(node_a, node_b, anisotropy_factor=0.1, use_haversine=False)
        self.assertAlmostEqual(calc_aniso, d_3d_aniso)

    def test_screen_overlap(self):
        """
        Verify logic for screened interval overlap.
        """
        # Case 1: Partial Overlap
        # A: 30-50, B: 40-60. Overlap 40-50 (10m).
        node_a = Node3D("A", x=0, y=0, z=0, elevation_m=0, screen_top=30, screen_bottom=50)
        node_b = Node3D("B", x=0, y=0, z=0, elevation_m=0, screen_top=40, screen_bottom=60)
        
        overlap_m, overlap_frac = compute_screen_overlap(node_a, node_b)
        self.assertEqual(overlap_m, 10.0)
        # Shorter screen is 20m vs 20m. Frac = 10/20 = 0.5
        self.assertEqual(overlap_frac, 0.5)
        
        # Case 2: No Overlap
        # A: 30-40, B: 50-60.
        node_c = Node3D("C", x=0, y=0, z=0, elevation_m=0, screen_top=50, screen_bottom=60)
        overlap_m, overlap_frac = compute_screen_overlap(node_a, node_c)
        self.assertEqual(overlap_m, 0.0)
        self.assertEqual(overlap_frac, 0.0)
        
        # Case 3: Complete Inclusion
        # A: 30-60 (30m), D: 40-50 (10m). Overlap 10m.
        node_d = Node3D("D", x=0, y=0, z=0, elevation_m=0, screen_top=40, screen_bottom=50)
        node_e = Node3D("E", x=0, y=0, z=0, elevation_m=0, screen_top=30, screen_bottom=60)
        
        overlap_m, overlap_frac = compute_screen_overlap(node_d, node_e)
        self.assertEqual(overlap_m, 10.0)
        # Shorter screen is D (10m). Frac = 10/10 = 1.0
        self.assertEqual(overlap_frac, 1.0)

    def test_edge_classification(self):
        """
        Verify geometric classification of 3D edges.
        """
        # Horizontal flow: High d_xy / d_z, same layer
        # ratio 100/1 = 100 > 10
        self.assertEqual(
            classify_edge_type(100.0, 1.0, same_layer=True), 
            "horizontal"
        )
        
        # Vertical leakage: Low d_xy / d_z, diff layer
        # ratio 1/100 = 0.01 < 0.1
        self.assertEqual(
            classify_edge_type(1.0, 100.0, same_layer=False),
            "vertical_leakage"
        )
        
        # Oblique: Intermediate
        # ratio 50/50 = 1.0
        self.assertEqual(
            classify_edge_type(50.0, 50.0, same_layer=False),
            "oblique"
        )

    def test_layer_probability(self):
        """
        Verify compound probability calculation for crossing multiple layers.
        """
        # Setup multi-layer system
        # Layer 1 -> 2 (p=0.5)
        # Layer 2 -> 3 (p=0.1)
        # Layer 1 -> 3 should be 0.5 * 0.1 = 0.05
        
        system = LayeredAquiferSystem(
            n_layers=3,
            layer_names=["L1", "L2", "L3"],
            layer_tops=[0, 10, 20],
            layer_bottoms=[10, 20, 30],
            aquitard_p=[0.5, 0.1],
            anisotropy_factors=[1, 1, 1]
        )
        
        # Direct neighbors
        self.assertAlmostEqual(get_aquitard_probability(system, 1, 2), 0.5)
        self.assertAlmostEqual(get_aquitard_probability(system, 2, 3), 0.1)
        
        # Skip one layer
        self.assertAlmostEqual(get_aquitard_probability(system, 1, 3), 0.05)
        
        # Same layer
        self.assertEqual(get_aquitard_probability(system, 1, 1), 1.0)

