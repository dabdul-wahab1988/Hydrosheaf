import pytest
from hydrosheaf.graph3d.build_3d import infer_edges_3d_probabilistic
from hydrosheaf.graph3d.types_3d import Node3D
from hydrosheaf.config import Config

class TestStratigraphicPriority:
    def test_prefer_same_layer_distant_over_cross_layer_close(self):
        """
        Verify that the new algorithm prefers a distant neighbor in the SAME layer
        over a very close neighbor in a DIFFERENT layer, if the cross-layer quota is filled
        or if probability favors the horizontal connection due to anisotropy.
        """
        config = Config()
        config.edge_max_neighbors = 2
        config.edge_max_neighbors_primary = 2 # Prioritize 2 same-layer
        config.edge_max_neighbors_secondary = 0 # Disallow cross-layer for this test
        config.vertical_anisotropy = 0.1 # Strong penalty for vertical distance
        config.edge_p_min = 0.000001 # Extremely low threshold
        config.edge_radius_km = 500.0 # Make sure 200km is within radius
        
        # Source Node: Layer 1
        src = Node3D(node_id="src", x=0, y=0, z=50, elevation_m=100, hydraulic_head=100, aquifer_layer=1)
        
        # Candidate 1: Close but Layer 2 (Vertical leakage candidate)
        # Delta Z = 10m. Effective dist = 10 / 0.1 = 100m.
        c1 = Node3D(node_id="c1_cross", x=0, y=0, z=60, elevation_m=100, hydraulic_head=99, aquifer_layer=2)
        
        # Candidate 2: Distant but Layer 1 (Horizontal flow candidate)
        # Delta X = 200m. Effective dist = 200m.
        # Even though 200m > 10m physical distance, the "Same Layer" priority should win
        # because secondary quota is 0.
        c2 = Node3D(node_id="c2_same", x=200, y=0, z=50, elevation_m=100, hydraulic_head=98, aquifer_layer=1)
        
        nodes = [src, c1, c2]
        edges = infer_edges_3d_probabilistic(nodes, config, layer_system=None, use_haversine=False)
        
        # Debugging
        print(f"\nTotal Edges Found: {len(edges)}")
        for e in edges:
            print(f"Edge: {e.u} -> {e.v} | Prob: {e.prob_combined} | Type: {e.edge_type}")
        
        # Should only select the source->c2 edge
        selected_targets = [e.v for e in edges if e.u == "src"]
        
        assert "c2_same" in selected_targets
        assert "c1_cross" not in selected_targets

    def test_secondary_quota_allows_leakage(self):
        """
        Verify that we can still get leakage connections if we allow secondary quota.
        """
        config = Config()
        config.edge_max_neighbors_primary = 1
        config.edge_max_neighbors_secondary = 1 # Allow 1 leakage
        config.edge_p_min = 0.000001
        config.edge_radius_km = 500.0
        
        src = Node3D(node_id="src", x=0, y=0, z=50, elevation_m=100, hydraulic_head=100, aquifer_layer=1)
        c1 = Node3D(node_id="c1_cross", x=0, y=0, z=60, elevation_m=100, hydraulic_head=99, aquifer_layer=2)
        c2 = Node3D(node_id="c2_same", x=200, y=0, z=50, elevation_m=100, hydraulic_head=98, aquifer_layer=1)
        
        nodes = [src, c1, c2]
        edges = infer_edges_3d_probabilistic(nodes, config, use_haversine=False)
        
        selected_targets = [e.v for e in edges if e.u == "src"]
        
        # Should have both: 1 primary (c2) and 1 secondary (c1)
        assert "c2_same" in selected_targets
        assert "c1_cross" in selected_targets
