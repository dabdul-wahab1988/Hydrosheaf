import unittest
import sys
from unittest import mock

from hydrosheaf.graph.build import infer_edges_probabilistic
from hydrosheaf.graph.head_inference import infer_heads_bayesian_mcmc


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

    def test_probabilistic_edges_bayesian_topography(self):
        samples = [
            {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 100.0},
            {"site_id": "B", "lat": 0.0, "lon": 0.01, "elevation": 90.0},
        ]
        edges = infer_edges_probabilistic(
            samples,
            radius_km=5.0,
            max_neighbors=1,
            p_min=0.6,
            head_key="head_meas",
            dtw_key="dtw",
            elevation_key="elevation",
            sigma_meas=0.5,
            sigma_dtw=1.0,
            sigma_elev=1.0,
            sigma_topo=5.0,
            gradient_min=0.0,
            depth_mismatch=20.0,
            head_inference="bayesian",
            dtw_prior_mu=5.0,
            dtw_prior_sigma=5.0,
        )
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->B", edge_ids)

    @unittest.skipIf(sys.platform.startswith("win"), "PyMC MCMC causes heap corruption on Windows CI")
    def test_probabilistic_edges_bayesian_mcmc_fallback(self):
        samples = [
            {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 100.0},
            {"site_id": "B", "lat": 0.0, "lon": 0.01, "elevation": 90.0},
        ]
        # PyMC is not guaranteed available in all environments. This call should
        # still work by falling back to the closed-form solver.
        edges = infer_edges_probabilistic(
            samples,
            radius_km=5.0,
            max_neighbors=1,
            p_min=0.6,
            head_key="head_meas",
            dtw_key="dtw",
            elevation_key="elevation",
            sigma_meas=0.5,
            sigma_dtw=1.0,
            sigma_elev=1.0,
            sigma_topo=5.0,
            gradient_min=0.0,
            depth_mismatch=20.0,
            head_inference="bayesian_mcmc",
            dtw_prior_mu=5.0,
            dtw_prior_sigma=5.0,
            mcmc_draws=50,
            mcmc_chains=1,
        )
        edge_ids = {edge.edge_id for edge in edges}
        self.assertIn("A->B", edge_ids)

    def test_head_mcmc_falls_back_to_linear_when_pymc_unavailable(self):
        samples = [
            {"site_id": "A", "elevation": 100.0},
            {"site_id": "B", "elevation": 90.0},
        ]

        with mock.patch(
            "hydrosheaf.graph.head_inference._load_pymc",
            side_effect=ImportError("pymc unavailable"),
        ):
            posterior = infer_heads_bayesian_mcmc(
                samples,
                node_id_key="site_id",
                head_key="head_meas",
                dtw_key="dtw",
                elevation_key="elevation",
                sigma_meas=0.5,
                sigma_dtw=1.0,
                sigma_elev=1.0,
                sigma_topo=5.0,
                dtw_prior_mu=5.0,
                dtw_prior_sigma=5.0,
                head_prior_mu=0.0,
                head_prior_sigma=1000.0,
                mcmc_draws=10,
                mcmc_chains=1,
            )

        self.assertEqual(set(posterior.node_ids), {"A", "B"})
        self.assertEqual(len(posterior.head_mean), 2)



if __name__ == "__main__":
    unittest.main()
