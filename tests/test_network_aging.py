"""
Tests for Network-Enhanced Bayesian Groundwater Dating.
"""

import unittest
import numpy as np
import networkx as nx
from hydrosheaf.nuclear.network_aging import (
    TracerObservationSet,
    infer_network_ages_bayesian,
)
from hydrosheaf.nuclear.input_history import (
    InputHistory,
    build_default_tritium_input,
)
from hydrosheaf.nuclear.lpm import convolve_input
from hydrosheaf.nuclear.nuclides import ARGON39


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
            val = convolve_input(
                cls.sample_date,
                age,
                cls.hist.years,
                cls.hist.values,
                lambda_y,
                model_type="PFM",
            )
            cls.obs[n] = float(val)
            cls.sigmas[n] = max(0.2, val * 0.1)  # 10% error

        print(f"Synthetic Obs: {cls.obs}")

    def test_bayesian_network_inference(self):
        """Test that network inference recovers age order correctly."""
        # SMC does not depend on Numba and is robust to separated bomb-peak modes.
        try:
            import pymc

            _ = pymc
        except Exception as exc:
            self.skipTest(f"PyMC unavailable in this environment: {exc}")

        results = infer_network_ages_bayesian(
            graph=self.graph,
            node_observations=self.obs,
            node_sigmas=self.sigmas,
            sample_date=self.sample_date,
            input_hist=self.hist,
            n_samples=300,
            n_chains=4,
            sampler="smc",
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

        # This 300-draw SMC run is a structural smoke test, not a
        # convergence-qualified posterior.  Check that diagnostics are finite
        # and that the API's explicit convergence flag applies the published
        # thresholds; do not mislabel a short run as converged by weakening the
        # production criterion.
        self.assertGreater(results["A"]["p_modern"], 0.5)
        diagnostics = results["_diagnostics"]
        self.assertEqual(diagnostics["sampler"], "smc")
        self.assertEqual(diagnostics["divergences"], 0)
        self.assertTrue(np.isfinite(diagnostics["age_r_hat_max"]))
        self.assertGreater(diagnostics["age_ess_bulk_min"], 0.0)
        self.assertGreater(diagnostics["age_ess_tail_min"], 0.0)
        expected_converged = bool(
            diagnostics["age_r_hat_max"] <= 1.01
            and diagnostics["age_ess_bulk_min"] >= 400.0
            and diagnostics["age_ess_tail_min"] >= 400.0
            and diagnostics["divergences"] == 0
        )
        self.assertEqual(diagnostics["converged"], expected_converged)
        self.assertLess(
            diagnostics["posterior_predictive"]["3H"]["standardized_rmse"],
            2.0,
        )

    def test_cycles_are_rejected_before_sampling(self):
        graph = nx.DiGraph([("A", "B"), ("B", "A")])
        with self.assertRaisesRegex(ValueError, "directed acyclic"):
            infer_network_ages_bayesian(
                graph=graph,
                node_observations={"A": 1.0, "B": 0.5},
                node_sigmas={"A": 0.1, "B": 0.1},
                sample_date=2025.0,
                n_samples=10,
                n_chains=1,
            )

    def test_exact_grid_multi_tracer_panel_converges(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(["young", "old"])
        ages = {"young": 5.0, "old": 80.0}
        tritium_history = InputHistory(
            np.asarray([1850.0, 2030.0]),
            np.asarray([6.0, 6.0]),
        )
        lambda_3h = np.log(2.0) / 12.32
        tritium = {node: 6.0 * np.exp(-lambda_3h * age) for node, age in ages.items()}
        argon = {
            node: 100.0 * np.exp(-np.log(2.0) * age / 269.0)
            for node, age in ages.items()
        }
        result = infer_network_ages_bayesian(
            graph,
            tritium,
            {node: 0.12 for node in ages},
            sample_date=2025.0,
            input_hist=tritium_history,
            additional_tracers=[
                TracerObservationSet(
                    name="39Ar",
                    nuclide=ARGON39,
                    observations=argon,
                    sigmas={node: 1.0 for node in ages},
                    sample_date=2025.0,
                )
            ],
            n_samples=500,
            n_chains=4,
            sampler="grid",
            max_age_years=150.0,
        )
        diagnostics = result["_diagnostics"]
        self.assertTrue(diagnostics["converged"])
        self.assertEqual(diagnostics["divergences"], 0)
        self.assertLess(abs(result["young"]["mean_age_years"] - 5.0), 5.0)
        self.assertLess(abs(result["old"]["mean_age_years"] - 80.0), 10.0)

    def test_exact_grid_rejects_network_edges(self):
        graph = nx.DiGraph([("A", "B")])
        with self.assertRaisesRegex(ValueError, "edge-free"):
            infer_network_ages_bayesian(
                graph,
                {"A": 5.0, "B": 2.0},
                {"A": 0.2, "B": 0.2},
                sample_date=2025.0,
                n_samples=20,
                n_chains=2,
                sampler="grid",
            )


if __name__ == "__main__":
    unittest.main()
