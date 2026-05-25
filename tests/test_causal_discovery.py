"""Tests for causal discovery layer."""

import unittest

import numpy as np

from hydrosheaf.config import Config
from hydrosheaf.causal.discovery import (
    CausalResult,
    assess_edge_causality,
    compute_causal_support,
    _common_driver_screen,
)


class TestCausalDiscovery(unittest.TestCase):
    def setUp(self):
        self.cfg = Config()

    def test_static_data_returns_insufficient_data(self):
        """Single-sample data with no temporal history returns insufficient_data."""
        upstream = {"Cl": 1.0, "Ca": 2.0}
        downstream = {"Cl": 1.2, "Ca": 1.8}
        result = assess_edge_causality(upstream, downstream, config=self.cfg)
        self.assertEqual(result.status, "insufficient_data")
        self.assertEqual(result.method, "insufficient_data")
        self.assertEqual(result.support_score, 0.0)

    def test_static_data_runs_common_driver_screen(self):
        """Common driver screen runs even with static data; results go to
        diagnostics but confounded_score is zero when temporal data is
        insufficient to distinguish causal flow from shared drivers."""
        upstream = {
            "Cl": 1.0, "Ca": 2.0,
            "aquifer_unit": "GravelAquifer",
            "screen_depth": 50.0,
            "elevation": 100.0,
        }
        downstream = {
            "Cl": 1.2, "Ca": 1.8,
            "aquifer_unit": "GravelAquifer",
            "screen_depth": 55.0,
            "elevation": 105.0,
        }
        result = assess_edge_causality(upstream, downstream, config=self.cfg)
        # Driver screen still runs — results in diagnostics
        self.assertIn("same_aquifer", result.diagnostics)
        self.assertTrue(result.diagnostics["same_aquifer"])
        self.assertIn("common_driver_flags", result.diagnostics)
        # But confounded_score is 0 because static data is insufficient
        self.assertEqual(result.confounded_score, 0.0)
        self.assertEqual(result.status, "insufficient_data")

    def test_synthetic_temporal_a_to_b_higher_support(self):
        """Synthetic temporal data: A -> B should have higher support for A->B
        than for B->A, using the lagged tracer correlation."""
        rng = np.random.RandomState(42)
        n = 20
        # Upstream A drives downstream B with a lag of 1
        a_data = rng.normal(10, 2, n)
        b_data = np.zeros(n)
        for i in range(1, n):
            b_data[i] = 0.7 * a_data[i - 1] + rng.normal(0.0, 0.5)

        # Build full ion vectors so partial correlation across species works
        ion_order = list(self.cfg.ion_order)
        cl_idx = ion_order.index("Cl")

        def _make_samples(data, site):
            samples = []
            for i in range(n):
                vec = {ion: 1.0 for ion in ion_order}  # baseline values
                vec["Cl"] = float(data[i])
                vec["site_id"] = site
                vec["date"] = f"2020-01-{i+1:02d}"
                samples.append(vec)
            return samples

        upstream_history = _make_samples(a_data, "A")
        downstream_history = _make_samples(b_data, "B")

        result_ab = assess_edge_causality(
            upstream_history[0], downstream_history[0],
            upstream_history=upstream_history,
            downstream_history=downstream_history,
            config=self.cfg,
        )

        # Forward direction A->B should find support via lagged correlation
        lag_support = result_ab.diagnostics.get("lag_support", 0.0)
        self.assertGreater(lag_support, 0.0)

        # Reverse direction B->A should have lower lagged support
        result_ba = assess_edge_causality(
            downstream_history[0], upstream_history[0],
            upstream_history=downstream_history,
            downstream_history=upstream_history,
            config=self.cfg,
        )
        ba_lag_support = result_ba.diagnostics.get("lag_support", 0.0)
        self.assertGreater(lag_support, ba_lag_support)

    def test_common_source_produces_confounding(self):
        """Two wells in the same aquifer with similar depth/elevation
        should have driver flags in diagnostics.  confounded_score is zero
        for insufficient_data because static snapshots cannot distinguish
        causal flow from shared-driver explanations."""
        upstream = {
            "Cl": 1.0, "aquifer_unit": "Sandstone",
            "screen_depth": 30.0, "elevation": 200.0,
        }
        downstream = {
            "Cl": 0.9, "aquifer_unit": "Sandstone",
            "screen_depth": 32.0, "elevation": 205.0,
        }
        result = assess_edge_causality(upstream, downstream, config=self.cfg)
        # Driver info is recorded in diagnostics for transparency
        self.assertIn("shared_aquifer", result.diagnostics.get("common_driver_flags", []))
        # But confounded_score is zero — insufficient temporal data
        self.assertEqual(result.confounded_score, 0.0)

    def test_different_aquifer_no_confounding(self):
        """Wells in different aquifers should not flag confounding."""
        upstream = {
            "Cl": 1.0, "aquifer_unit": "Sandstone",
            "screen_depth": 30.0, "elevation": 200.0,
        }
        downstream = {
            "Cl": 2.0, "aquifer_unit": "Limestone",
            "screen_depth": 80.0, "elevation": 150.0,
        }
        result = assess_edge_causality(upstream, downstream, config=self.cfg)
        # Different aquifer, depth, elevation -> no confounding
        self.assertEqual(result.confounded_score, 0.0)

    def test_compute_causal_support_returns_all_fields(self):
        """Convenience wrapper returns all expected dict fields."""
        upstream = {"Cl": 1.0, "Ca": 2.0}
        downstream = {"Cl": 1.2, "Ca": 1.8}
        result = compute_causal_support(upstream, downstream, config=self.cfg)
        for key in [
            "causal_support_score",
            "causal_confounded_score",
            "causal_p_value",
            "causal_method",
            "causal_n_observations",
            "causal_status",
        ]:
            self.assertIn(key, result)

    def test_min_observations_respected(self):
        """Insufficient history (< causal_min_observations) returns insufficient_data."""
        self.cfg.causal_min_observations = 10
        upstream_history = [
            {"site_id": "A", "Cl": 1.0, "date": f"2020-01-{i+1:02d}"}
            for i in range(5)  # Only 5 observations, < 10 required
        ]
        downstream_history = list(upstream_history)

        result = assess_edge_causality(
            upstream_history[0], downstream_history[0],
            upstream_history=upstream_history,
            downstream_history=downstream_history,
            config=self.cfg,
        )
        self.assertEqual(result.status, "insufficient_data")


class TestCommonDriverScreen(unittest.TestCase):
    def test_no_shared_attributes(self):
        cfg = Config()
        u = {"aquifer_unit": "A", "screen_depth": 10.0, "elevation": 100.0}
        v = {"aquifer_unit": "B", "screen_depth": 50.0, "elevation": 200.0}
        info = _common_driver_screen(u, v, cfg)
        self.assertEqual(info["confound_score"], 0.0)

    def test_all_shared_attributes(self):
        cfg = Config()
        u = {"aquifer_unit": "Same", "screen_depth": 25.0, "elevation": 150.0}
        v = {"aquifer_unit": "Same", "screen_depth": 28.0, "elevation": 155.0}
        info = _common_driver_screen(u, v, cfg)
        self.assertGreater(info["confound_score"], 0.5)
        self.assertEqual(len(info["common_driver_flags"]), 3)


if __name__ == "__main__":
    unittest.main()
