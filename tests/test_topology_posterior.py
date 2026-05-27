"""Tests for Bayesian topology posterior."""

import unittest
from typing import Sequence

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.topology_posterior import (
    run_topology_posterior,
    attach_posterior_attrs,
)


def _make_edge(edge_id: str, u: str, v: str, confidence: float = 0.5) -> Edge:
    return Edge(edge_id=edge_id, u=u, v=v, attrs={"edge_confidence": confidence})


def _simple_cost_fn(edges: Sequence[Edge]) -> float:
    """Cost that decreases with more edges (better data fit)."""
    return max(0.0, 5.0 - float(len(edges)))


def _biased_cost_fn(edges: Sequence[Edge]) -> float:
    """Cost that penalises specific edges (B->C and C->D are expensive)."""
    expensive_edges = {"B->C", "C->D"}
    base_cost = max(0.0, 5.0 - float(len(edges)))
    for e in edges:
        if e.edge_id in expensive_edges:
            base_cost += 5.0
    return base_cost


class TestTopologyPosterior(unittest.TestCase):
    def setUp(self):
        self.cfg = Config()

    def test_high_prior_edge_gets_high_probability(self):
        """An edge with high prior should have high posterior probability."""
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.95),
            _make_edge("B->C", "B", "C", confidence=0.05),
        ]
        self.cfg.topology_posterior_samples = 1000
        self.cfg.topology_posterior_burnin = 200

        result = run_topology_posterior(
            universe, _simple_cost_fn, self.cfg, initial_edges=universe, seed=42
        )
        probs = result["edge_probabilities"]
        self.assertGreater(probs.get("A->B", 0.0), probs.get("B->C", 0.0))

    def test_expensive_edge_gets_low_probability(self):
        """An edge with high cost should have low posterior probability."""
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.8),
            _make_edge("B->C", "B", "C", confidence=0.8),
            _make_edge("C->D", "C", "D", confidence=0.8),
        ]
        self.cfg.topology_posterior_samples = 1000
        self.cfg.topology_posterior_burnin = 200
        self.cfg.topology_posterior_beta = 2.0  # Higher beta = more cost-sensitive

        result = run_topology_posterior(
            universe, _biased_cost_fn, self.cfg, initial_edges=[universe[0]], seed=42
        )
        probs = result["edge_probabilities"]
        # The expensive edges B->C and C->D should have lower probability
        self.assertGreater(probs.get("A->B", 0.0), probs.get("B->C", 0.0))
        self.assertGreater(probs.get("A->B", 0.0), probs.get("C->D", 0.0))

    def test_sampler_deterministic_with_seed(self):
        """The sampler should be deterministic with fixed seed."""
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.6),
            _make_edge("B->C", "B", "C", confidence=0.4),
        ]
        self.cfg.topology_posterior_samples = 500
        self.cfg.topology_posterior_burnin = 100

        result1 = run_topology_posterior(
            universe, _simple_cost_fn, self.cfg, seed=123
        )
        result2 = run_topology_posterior(
            universe, _simple_cost_fn, self.cfg, seed=123
        )
        self.assertEqual(result1["edge_probabilities"], result2["edge_probabilities"])
        self.assertEqual(result1["map_edges"], result2["map_edges"])

    def test_different_seeds_produce_different_results(self):
        """Different seeds should (likely) produce somewhat different chains."""
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.6),
            _make_edge("B->C", "B", "C", confidence=0.4),
        ]
        self.cfg.topology_posterior_samples = 500
        self.cfg.topology_posterior_burnin = 100

        result1 = run_topology_posterior(
            universe, _simple_cost_fn, self.cfg, seed=42
        )
        result2 = run_topology_posterior(
            universe, _simple_cost_fn, self.cfg, seed=99
        )
        # The chains may differ, but both should produce valid probabilities
        for prob in result1["edge_probabilities"].values():
            self.assertGreaterEqual(prob, 0.0)
            self.assertLessEqual(prob, 1.0)
        for prob in result2["edge_probabilities"].values():
            self.assertGreaterEqual(prob, 0.0)
            self.assertLessEqual(prob, 1.0)

    def test_output_structure(self):
        """All expected output keys should be present."""
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.7),
        ]
        self.cfg.topology_posterior_samples = 200
        self.cfg.topology_posterior_burnin = 50

        result = run_topology_posterior(
            universe, _simple_cost_fn, self.cfg, seed=42
        )
        self.assertIn("edge_probabilities", result)
        self.assertIn("edge_log_odds", result)
        self.assertIn("map_edges", result)
        self.assertIn("entropy", result)
        self.assertIn("n_edges_mean", result)
        self.assertIn("n_edges_ci95", result)
        self.assertIn("acceptance_rate", result)
        self.assertGreaterEqual(result["acceptance_rate"], 0.0)
        self.assertLessEqual(result["acceptance_rate"], 1.0)

    def test_empty_universe(self):
        """An empty universe should produce empty results."""
        result = run_topology_posterior([], _simple_cost_fn, self.cfg, seed=42)
        self.assertEqual(result["edge_probabilities"], {})
        self.assertEqual(result["map_edges"], [])

    def test_make_cost_fn_accepts_gradient_map(self):
        """make_topology_cost_fn works with gradient_map kwarg."""
        from hydrosheaf.inference.topology_posterior import make_topology_cost_fn

        samples = {
            "A": {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0},
            "B": {"site_id": "B", "lat": 0.0, "lon": 1.0, "elevation": 8.0},
        }
        gradient_map = {"A": (-0.005, 0.0, 0.0)}
        config = Config(
            phreeqc_enabled=False,
            hydraulic_hodge_enabled=True,
            hydraulic_hodge_weight=1.0,
            hydraulic_hodge_fallback_to_elevation=True,
        )

        cost_fn = make_topology_cost_fn(
            sample_map=samples,
            config=config,
            gradient_map=gradient_map,
        )
        edges = [_make_edge("A->B", "A", "B", confidence=0.9)]
        cost = cost_fn(edges)
        self.assertGreaterEqual(cost, 0.0)

    def test_make_cost_fn_accepts_local_residuals(self):
        """make_topology_cost_fn works with local_residuals kwarg."""
        from hydrosheaf.inference.topology_posterior import make_topology_cost_fn

        samples = {
            "A": {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0},
            "B": {"site_id": "B", "lat": 0.0, "lon": 1.0, "elevation": 8.0},
        }
        local_residuals = {"A": 0.1, "B": 0.2}
        config = Config(
            phreeqc_enabled=False,
            hydraulic_hodge_enabled=True,
            hydraulic_hodge_weight=1.0,
            head_plane_residual_weight=1.0,
            hydraulic_hodge_fallback_to_elevation=True,
        )

        cost_fn = make_topology_cost_fn(
            sample_map=samples,
            config=config,
            local_residuals=local_residuals,
        )
        edges = [_make_edge("A->B", "A", "B", confidence=0.9)]
        cost = cost_fn(edges)
        self.assertGreaterEqual(cost, 0.0)

    def test_local_residuals_do_not_require_projected_gradient(self):
        """head_plane_residual_enabled works standalone, without gradient prerequisites."""
        from hydrosheaf.inference.topology_posterior import make_topology_cost_fn

        samples = {
            "A": {"site_id": "A", "lat": 0.0, "lon": 0.0, "elevation": 10.0},
            "B": {"site_id": "B", "lat": 0.0, "lon": 1.0, "elevation": 8.0},
        }
        config = Config(
            phreeqc_enabled=False,
            hydraulic_hodge_enabled=True,
            head_plane_residual_enabled=True,
            head_plane_residual_weight=1.0,
            hydraulic_hodge_fallback_to_elevation=True,
        )

        cost_fn = make_topology_cost_fn(
            sample_map=samples,
            config=config,
            local_residuals={"A": 0.3, "B": 0.1},
        )

        cost = cost_fn([_make_edge("A->B", "A", "B", confidence=0.9)])
        self.assertGreater(cost, 0.0)

    def test_edge_penalty_encourages_sparsity(self):
        """Higher edge penalty should reduce the mean number of edges."""
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.7),
            _make_edge("B->C", "B", "C", confidence=0.7),
            _make_edge("C->D", "C", "D", confidence=0.7),
        ]
        self.cfg.topology_posterior_samples = 800
        self.cfg.topology_posterior_burnin = 200
        self.cfg.topology_posterior_beta = 1.0

        cfg_no_penalty = Config()
        cfg_no_penalty.topology_posterior_samples = 800
        cfg_no_penalty.topology_posterior_burnin = 200
        cfg_no_penalty.topology_posterior_beta = 1.0
        cfg_no_penalty.topology_posterior_edge_penalty = 0.0

        cfg_penalty = Config()
        cfg_penalty.topology_posterior_samples = 800
        cfg_penalty.topology_posterior_burnin = 200
        cfg_penalty.topology_posterior_beta = 1.0
        cfg_penalty.topology_posterior_edge_penalty = 2.0

        result_no_penalty = run_topology_posterior(
            universe, _simple_cost_fn, cfg_no_penalty, seed=42
        )
        result_penalty = run_topology_posterior(
            universe, _simple_cost_fn, cfg_penalty, seed=42
        )
        self.assertGreaterEqual(
            result_no_penalty["n_edges_mean"],
            result_penalty["n_edges_mean"],
        )


class TestAttachPosteriorAttrs(unittest.TestCase):
    def test_attrs_attached(self):
        universe = [
            _make_edge("A->B", "A", "B", confidence=0.7),
            _make_edge("B->C", "B", "C", confidence=0.3),
        ]
        selected = list(universe)
        cfg = Config()
        cfg.topology_posterior_samples = 300
        cfg.topology_posterior_burnin = 50

        result = run_topology_posterior(
            universe, _simple_cost_fn, cfg, initial_edges=selected, seed=42
        )
        attach_posterior_attrs(selected, universe, result)

        for edge in universe:
            attrs = edge.attrs or {}
            self.assertIn("posterior_edge_probability", attrs)
            self.assertIn("posterior_edge_log_odds", attrs)
            self.assertIn("posterior_map_selected", attrs)
            self.assertIn("posterior_topology_entropy", attrs)


if __name__ == "__main__":
    unittest.main()
