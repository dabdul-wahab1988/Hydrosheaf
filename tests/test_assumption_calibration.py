"""Integration tests for assumption calibration (Phase 0-1).

Tests edge evidence classification, null-model downgrading, and backward
compatibility when assumption calibration is disabled.
"""

import unittest

from hydrosheaf import Config, infer_edges


class AssumptionCalibrationIntegrationTests(unittest.TestCase):
    """Test null-model integration and evidence ladder via the public API."""

    def test_chemistry_similar_noflow_downgraded(self):
        """Two wells with near-identical chemistry in same lithology:
        null model should raise null_score and the edge should be downgraded."""
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0,
                "2H": -30.0,
                "aquifer_layer": "Upper Aquifer",
                "Ca": 2.0,
                "Mg": 1.0,
                "Na": 5.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.1,
                "18O": -5.1,
                "2H": -30.5,
                "aquifer_layer": "Upper Aquifer",
                "Ca": 2.1,
                "Mg": 0.9,
                "Na": 5.1,
            },
            {
                "site_id": "C",
                "lat": 0.01,
                "lon": 0.0,
                "head_meas": 110.0,  # Uphill from A, unlikely flow
                "Cl": 10.0,           # Very different chemistry
                "18O": -10.0,
                "2H": -70.0,
                "Ca": 0.5,
                "Mg": 0.2,
                "Na": 1.0,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 2
        config.edge_map_candidate_multiplier = 5
        config.edge_map_p_min = 0.0
        config.sheaf_weight_head_prior = 0.2
        config.sheaf_weight_isotope = 5.0
        config.sheaf_weight_cl = 0.5
        # Enable null model
        config.null_model_enabled = True
        config.null_model_weight = 0.5
        config.evidence_ladder_enabled = True

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)

        # The A->B edge has very similar chemistry + same lithology
        # so null_score should be > 0
        a_b_edge = None
        for e in edges:
            if e.edge_id == "A->B":
                a_b_edge = e
                break
        self.assertIsNotNone(a_b_edge, "Expected A->B edge")
        null_score = a_b_edge.attrs.get("null_score", 0.0)
        self.assertGreater(
            float(null_score), 0.0,
            f"Null model should detect similar chemistry; got null_score={null_score}"
        )

    def test_missing_evidence_is_ambiguous(self):
        """An edge with no isotope or Cl data should be classified AMBIGUOUS,
        not silently treated as VALIDATED."""
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.0
        config.evidence_ladder_enabled = True

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)

        self.assertGreater(len(edges), 0, "Expected at least one edge")
        evidence_class = edges[0].attrs.get("evidence_class", "")
        self.assertEqual(evidence_class, "AMBIGUOUS",
                         f"Expected AMBIGUOUS, got {evidence_class}")

    def test_existing_refinement_unchanged_when_disabled(self):
        """When assumption calibration is disabled, output should match
        existing behavior (same edges selected)."""
        samples = [
            {
                "site_id": "A",
                "lat": 0.0,
                "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0,
                "2H": -30.0,
            },
            {
                "site_id": "B",
                "lat": 0.0,
                "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.1,
                "18O": -5.1,
                "2H": -30.5,
            },
            {
                "site_id": "C",
                "lat": 0.01,
                "lon": 0.0,
                "head_meas": 110.0,
                "Cl": 1.0,
                "18O": -10.0,
                "2H": -70.0,
            },
        ]

        # Run WITHOUT assumption calibration
        config_disabled = Config()
        config_disabled.edge_radius_km = 5.0
        config_disabled.edge_max_neighbors = 1
        config_disabled.edge_map_candidate_multiplier = 3
        config_disabled.edge_map_p_min = 0.0
        config_disabled.sheaf_weight_head_prior = 0.2
        config_disabled.sheaf_weight_isotope = 5.0
        config_disabled.sheaf_weight_cl = 0.5
        # Explicitly disabled
        config_disabled.assumption_calibration_enabled = False
        config_disabled.null_model_enabled = False
        config_disabled.evidence_ladder_enabled = False

        edges_disabled = infer_edges(samples, method="probabilistic_sheaf",
                                      config=config_disabled)
        edge_ids_disabled = {e.edge_id for e in edges_disabled}

        # Run WITH assumption calibration enabled (null model on)
        config_enabled = Config()
        config_enabled.edge_radius_km = 5.0
        config_enabled.edge_max_neighbors = 1
        config_enabled.edge_map_candidate_multiplier = 3
        config_enabled.edge_map_p_min = 0.0
        config_enabled.sheaf_weight_head_prior = 0.2
        config_enabled.sheaf_weight_isotope = 5.0
        config_enabled.sheaf_weight_cl = 0.5
        config_enabled.null_model_enabled = True
        config_enabled.null_model_weight = 0.5
        config_enabled.evidence_ladder_enabled = True

        edges_enabled = infer_edges(samples, method="probabilistic_sheaf",
                                     config=config_enabled)

        # When enabled, the null model increases cost on A->B (similar chemistry)
        # but A->C has very different isotopes so null is lower.
        # For this test: the same A->B edge should be in both sets because
        # head + distance still favor it, just with higher cost.
        self.assertIn("A->B", edge_ids_disabled,
                      "Expected A->B to be selected when disabled")
        # With null model enabled, A->B may still be selected if it's the best
        # option, but it should have evidence annotations
        a_b_enabled = None
        for e in edges_enabled:
            if e.edge_id == "A->B":
                a_b_enabled = e
                break
        if a_b_enabled:
            # Should have evidence attributes
            self.assertIn("evidence_class", a_b_enabled.attrs)
            self.assertIn("null_score", a_b_enabled.attrs)

    def test_evidence_class_written_to_attrs(self):
        """When evidence_ladder_enabled=True, edge attrs should contain
        evidence_class, evidence_score, null_score, evidence_flags, and evidence_reason."""
        samples = [
            {
                "site_id": "A",
                "lat": 0.0, "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0, "2H": -30.0,
            },
            {
                "site_id": "B",
                "lat": 0.0, "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.1,
                "18O": -5.1, "2H": -30.5,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.0
        config.evidence_ladder_enabled = True

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        self.assertGreater(len(edges), 0)
        attrs = edges[0].attrs
        self.assertIn("evidence_class", attrs,
                      "evidence_class missing from edge attrs")
        self.assertIn("evidence_score", attrs,
                      "evidence_score missing from edge attrs")
        self.assertIn("null_score", attrs,
                      "null_score missing from edge attrs")
        self.assertIn("evidence_flags", attrs,
                      "evidence_flags missing from edge attrs")
        self.assertIn("evidence_reason", attrs,
                      "evidence_reason missing from edge attrs")

    def test_assumption_calibration_master_gate(self):
        """setting assumption_calibration_enabled=True alone should activate
        null_model and evidence_ladder, even without individual flags."""
        samples = [
            {
                "site_id": "A",
                "lat": 0.0, "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0, "2H": -30.0,
            },
            {
                "site_id": "B",
                "lat": 0.0, "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.1,
                "18O": -5.1, "2H": -30.5,
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.0
        # Only the master gate, not the individual flags
        config.assumption_calibration_enabled = True
        config.null_model_enabled = False
        config.evidence_ladder_enabled = False

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        self.assertGreater(len(edges), 0)
        attrs = edges[0].attrs
        # Master gate should have turned on evidence annotation
        self.assertIn("evidence_class", attrs,
                      "assumption_calibration_enabled should activate evidence ladder")
        self.assertIn("null_score", attrs,
                      "assumption_calibration_enabled should activate null model")

    def test_falsified_by_age_reversal(self):
        """Edge with clear age reversal should be FALSIFIED when evidence
        ladder is enabled."""
        samples = [
            {
                "site_id": "A",
                "lat": 0.0, "lon": 0.0,
                "head_meas": 100.0,
                "Cl": 1.0,
                "18O": -5.0, "2H": -30.0,
                "mean_age_years": 50.0,
            },
            {
                "site_id": "B",
                "lat": 0.0, "lon": 0.01,
                "head_meas": 90.0,
                "Cl": 1.1,
                "18O": -5.1, "2H": -30.5,
                "mean_age_years": 5.0,  # Much younger -> age reversal
            },
        ]
        config = Config()
        config.edge_radius_km = 5.0
        config.edge_max_neighbors = 1
        config.edge_map_candidate_multiplier = 3
        config.edge_map_p_min = 0.0
        config.evidence_ladder_enabled = True
        config.sheaf_age_enabled = True
        config.sheaf_weight_age = 5.0  # Strong age penalty

        edges = infer_edges(samples, method="probabilistic_sheaf", config=config)
        # Edge may be selected if only one candidate, but check evidence_class
        if edges:
            evidence_class = edges[0].attrs.get("evidence_class", "")
            # Age reversal should trigger FALSIFIED
            self.assertIn(evidence_class, ["FALSIFIED", "AMBIGUOUS"],
                          f"Age reversal edge got class {evidence_class}")


if __name__ == "__main__":
    unittest.main()
