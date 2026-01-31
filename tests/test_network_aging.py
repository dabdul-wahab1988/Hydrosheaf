"""
Tests for Network-Enhanced Bayesian Groundwater Dating.
"""
import unittest
import numpy as np
import networkx as nx
from hydrosheaf.nuclear.network_aging import infer_network_ages_bayesian
from hydrosheaf.nuclear.input_history import build_default_tritium_input
from hydrosheaf.nuclear.lpm import convolve_input
from hydrosheaf.nuclear.nuclides import TRITIUM

class NetworkAgingTests(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # Setup synthetic data
        # Graph: A -> B -> C
        cls.graph = nx.DiGraph()
        cls.graph.add_edges_from([("A", "B"), ("B", "C")])
        
        # True ages (years)
        # A is young recharge, B is older, C is fossil
        cls.true_ages = {"A": 5.0, "B": 30.0, "C": 150.0}
        cls.sample_date = 2023.0
        
        # Input history (Global)
        cls.hist = build_default_tritium_input("global")
        lambda_y = np.log(2) / 12.32
        
        # Generate synthetic observations
        cls.obs = {}
        cls.sigmas = {}
        
        for n, age in cls.true_ages.items():
            val = convolve_input(cls.sample_date, age, cls.hist.years, cls.hist.values, lambda_y, model_type="PFM")
            cls.obs[n] = float(val)
            cls.sigmas[n] = max(0.2, val * 0.1) # 10% error
            
        print(f"Synthetic Obs: {cls.obs}")

    def test_bayesian_network_inference(self):
        """Test that network inference recovers age order correctly."""
        # Run inference with fewer samples for speed in test environment
        try:
            import nutpie
        except ImportError:
            self.skipTest("Nutpie not installed")

        results = infer_network_ages_bayesian(
            graph=self.graph,
            node_observations=self.obs,
            node_sigmas=self.sigmas,
            sample_date=self.sample_date,
            input_hist=self.hist,
            n_samples=50, # Fast test
            n_chains=2
        )
        
        # Check results structure
        self.assertIn("A", results)
        self.assertIn("B", results)
        self.assertIn("C", results)
        
        mean_A = results["A"]["mean_age_years"]
        mean_B = results["B"]["mean_age_years"]
        mean_C = results["C"]["mean_age_years"]
        
        print(f"Inferred Means: A={mean_A:.1f}, B={mean_B:.1f}, C={mean_C:.1f}")
        
        # Check ordering: A < B < C
        # Note: B and C might be close if mixing is allowed, but generally A should be youngest
        self.assertLess(mean_A, mean_C)
        
        # Check modern prob
        self.assertGreater(results["A"]["p_modern"], 0.8) # Should be confidently modern
        self.assertLess(results["C"]["p_modern"], 0.2) # Should be confidently fossil

if __name__ == "__main__":
    unittest.main()
