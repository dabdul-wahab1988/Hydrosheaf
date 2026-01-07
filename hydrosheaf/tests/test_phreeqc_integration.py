"""
Comprehensive Integration Tests for Hydrosheaf PHREEQC Module

This test suite validates the full PHREEQC integration workflow including:
1. Database loading and SI calculations
2. Reaction constraint building from SI values  
3. Full network fitting with PHREEQC constraints
4. Edge-level reaction fitting with bounds
"""

import unittest
from pathlib import Path

from hydrosheaf.config import Config, default_config
from hydrosheaf.phreeqc.runner import run_phreeqc
from hydrosheaf.phreeqc.constraints import build_edge_bounds
from hydrosheaf.models.reactions import build_reaction_dictionary, fit_reactions
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.inference.network_fit import fit_network, summarize_network


class PhreeqcIntegrationTests(unittest.TestCase):
    """Test PHREEQC integration with actual database loading and SI calculations."""

    @classmethod
    def setUpClass(cls):
        """Create test samples once for all tests."""
        cls.samples = [
            {
                "sample_id": "GW001",
                "pH": 7.5,
                "temp_c": 25.0,
                "Ca": 80.0,
                "Mg": 25.0,
                "Na": 45.0,
                "Cl": 35.0,
                "SO4": 120.0,
                "HCO3": 250.0,
                "NO3": 10.0,
                "F": 0.5,
            },
            {
                "sample_id": "GW002",
                "pH": 6.8,
                "temp_c": 22.0,
                "Ca": 40.0,
                "Mg": 15.0,
                "Na": 20.0,
                "Cl": 25.0,
                "SO4": 60.0,
                "HCO3": 150.0,
                "NO3": 5.0,
                "F": 0.3,
            },
            {
                "sample_id": "GW003",
                "pH": 8.2,
                "temp_c": 28.0,
                "Ca": 120.0,
                "Mg": 40.0,
                "Na": 80.0,
                "Cl": 90.0,
                "SO4": 200.0,
                "HCO3": 320.0,
                "NO3": 2.0,
                "F": 0.8,
            },
        ]
        cls.config = default_config()
        cls.config.missing_policy = "impute_zero"

    def test_phreeqc_calculates_si_values(self):
        """Test that PHREEQC calculates saturation indices for all minerals."""
        results = run_phreeqc(self.samples, self.config)
        
        for sample_id in ["GW001", "GW002", "GW003"]:
            result = results[sample_id]
            self.assertTrue(result["phreeqc_ok"], f"PHREEQC failed for {sample_id}")
            self.assertIsNone(result["skipped_reason"])
            
            # Check all SI values are calculated (not None)
            self.assertIsNotNone(result["ionic_strength"])
            self.assertIsNotNone(result["si_calcite"])
            self.assertIsNotNone(result["si_dolomite"])
            self.assertIsNotNone(result["si_gypsum"])
            self.assertIsNotNone(result["si_halite"])
            self.assertIsNotNone(result["si_fluorite"])
            
            # SI values should be finite numbers
            self.assertTrue(-50 < result["si_calcite"] < 50)
            self.assertTrue(-50 < result["si_dolomite"] < 50)

    def test_ionic_strength_reasonable(self):
        """Test that ionic strength values are reasonable for freshwater."""
        results = run_phreeqc(self.samples, self.config)
        
        for sample_id in ["GW001", "GW002", "GW003"]:
            I = results[sample_id]["ionic_strength"]
            # Freshwater typically has I between 0.001 and 0.1
            self.assertGreater(I, 0.001)
            self.assertGreater(I, 0.001)
            self.assertLess(I, 1.0)

    def test_si_calcite_varies_with_ph(self):
        """Test that calcite SI increases with pH (more supersaturated)."""
        results = run_phreeqc(self.samples, self.config)
        
        # GW002 (pH 6.8) should have lower SI calcite than GW003 (pH 8.2)
        si_low_ph = results["GW002"]["si_calcite"]
        si_high_ph = results["GW003"]["si_calcite"]
        self.assertLess(si_low_ph, si_high_ph)


class ConstraintsBuildingTests(unittest.TestCase):
    """Test PHREEQC constraints building from SI values."""

    def test_dissolution_only_constraint(self):
        """Test that undersaturated minerals get dissolution-only constraint."""
        config = default_config()
        config.si_threshold_tau = 0.2
        
        # Mock PHREEQC results with undersaturated calcite
        phreeqc_results = {
            "A": {"phreeqc_ok": True, "si_calcite": -1.5, "si_dolomite": -2.0, 
                  "si_gypsum": -2.5, "si_halite": -8.0, "si_fluorite": -2.0},
            "B": {"phreeqc_ok": True, "si_calcite": -1.0, "si_dolomite": -1.5,
                  "si_gypsum": -2.0, "si_halite": -7.5, "si_fluorite": -1.5},
        }
        
        edges = [("A", "B")]
        matrix, labels, mineral_mask = build_reaction_dictionary(config)
        bounds = build_edge_bounds(phreeqc_results, edges, labels, mineral_mask, config)
        
        edge_bounds = bounds["A->B"]
        self.assertEqual(edge_bounds["constraints_active"]["calcite"], "dissolution_only")
        self.assertEqual(edge_bounds["constraints_active"]["halite"], "dissolution_only")
        
        # lb should be 0.0 for dissolution-only (no precipitation)
        calcite_idx = labels.index("calcite")
        self.assertEqual(edge_bounds["lb"][calcite_idx], 0.0)

    def test_precipitation_only_constraint(self):
        """Test that supersaturated downstream minerals get precipitation-only constraint."""
        config = default_config()
        config.si_threshold_tau = 0.2
        
        # Mock PHREEQC results with supersaturated calcite at v
        phreeqc_results = {
            "A": {"phreeqc_ok": True, "si_calcite": 0.0, "si_dolomite": -0.5,
                  "si_gypsum": -2.0, "si_halite": -8.0, "si_fluorite": -2.0},
            "B": {"phreeqc_ok": True, "si_calcite": 1.5, "si_dolomite": 1.0,
                  "si_gypsum": -1.5, "si_halite": -7.0, "si_fluorite": -1.0},
        }
        
        edges = [("A", "B")]
        matrix, labels, mineral_mask = build_reaction_dictionary(config)
        bounds = build_edge_bounds(phreeqc_results, edges, labels, mineral_mask, config)
        
        edge_bounds = bounds["A->B"]
        self.assertEqual(edge_bounds["constraints_active"]["calcite"], "precipitation_only")
        
        # ub should be 0.0 for precipitation-only (no dissolution)
        calcite_idx = labels.index("calcite")
        self.assertEqual(edge_bounds["ub"][calcite_idx], 0.0)

    def test_free_constraint_near_equilibrium(self):
        """Test that near-equilibrium minerals get free constraint."""
        config = default_config()
        config.si_threshold_tau = 0.2
        
        # Mock PHREEQC results with calcite near equilibrium
        phreeqc_results = {
            "A": {"phreeqc_ok": True, "si_calcite": 0.1, "si_dolomite": 0.0,
                  "si_gypsum": -0.1, "si_halite": -8.0, "si_fluorite": -2.0},
            "B": {"phreeqc_ok": True, "si_calcite": -0.1, "si_dolomite": 0.1,
                  "si_gypsum": 0.0, "si_halite": -7.5, "si_fluorite": -1.5},
        }
        
        edges = [("A", "B")]
        matrix, labels, mineral_mask = build_reaction_dictionary(config)
        bounds = build_edge_bounds(phreeqc_results, edges, labels, mineral_mask, config)
        
        edge_bounds = bounds["A->B"]
        self.assertEqual(edge_bounds["constraints_active"]["calcite"], "free")


class ReactionFittingTests(unittest.TestCase):
    """Test reaction fitting with various configurations."""

    def test_fit_multiple_reactions(self):
        """Test fitting multiple concurrent reactions."""
        config = default_config()
        matrix, labels, _ = build_reaction_dictionary(config)
        
        calcite_idx = labels.index("calcite")
        halite_idx = labels.index("halite")
        
        # Create residual from calcite + halite dissolution
        calcite = matrix[calcite_idx]
        halite = matrix[halite_idx]
        z_calcite = 0.5
        z_halite = 0.3
        residual = [z_calcite * c + z_halite * h for c, h in zip(calcite, halite)]
        
        fit = fit_reactions(
            residual, matrix, 
            weights=[1.0] * 8,
            lambda_l1=0.0
        )
        
        self.assertAlmostEqual(fit.extents[calcite_idx], z_calcite, places=2)
        self.assertAlmostEqual(fit.extents[halite_idx], z_halite, places=2)
        self.assertTrue(fit.converged)

    def test_fit_with_bounds(self):
        """Test reaction fitting respects bounds constraints."""
        config = default_config()
        matrix, labels, _ = build_reaction_dictionary(config)
        
        calcite_idx = labels.index("calcite")
        
        # Residual that would require calcite dissolution
        residual = [1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]  # Ca and HCO3
        
        # But set precipitation-only constraint (ub=0)
        lb = [-float("inf")] * len(labels)
        ub = [float("inf")] * len(labels)
        ub[calcite_idx] = 0.0  # No dissolution allowed
        
        fit = fit_reactions(
            residual, matrix,
            weights=[1.0] * 8,
            lambda_l1=0.0,
            lb=lb, ub=ub
        )
        
        # Calcite extent should be <= 0 (precipitation only)
        self.assertLessEqual(fit.extents[calcite_idx], 0.0 + 1e-6)

    def test_l1_regularization_sparsity(self):
        """Test that L1 regularization promotes sparsity."""
        config = default_config()
        matrix, labels, _ = build_reaction_dictionary(config)
        
        # Residual that could be explained by multiple reactions
        residual = [0.5, 0.3, 0.2, 1.0, 0.1, 0.4, 0.1, 0.05]
        
        # Fit without regularization
        fit_no_reg = fit_reactions(
            residual, matrix,
            weights=[1.0] * 8,
            lambda_l1=0.0
        )
        
        # Fit with strong regularization
        fit_reg = fit_reactions(
            residual, matrix,
            weights=[1.0] * 8,
            lambda_l1=1.0
        )
        
        # Regularized fit should have smaller L1 norm (more sparse)
        self.assertLess(fit_reg.l1_norm, fit_no_reg.l1_norm)


class NetworkIntegrationTests(unittest.TestCase):
    """Test full network fitting with PHREEQC integration."""

    def test_network_fit_with_phreeqc(self):
        """Test complete network fitting workflow with PHREEQC."""
        samples = [
            {
                "site_id": "A",
                "pH": 7.0,
                "temp_c": 25.0,
                "Ca": 50.0,
                "Mg": 15.0,
                "Na": 25.0,
                "Cl": 20.0,
                "SO4": 60.0,
                "HCO3": 180.0,
                "NO3": 5.0,
                "F": 0.3,
            },
            {
                "site_id": "B",
                "pH": 7.5,
                "temp_c": 25.0,
                "Ca": 65.0,
                "Mg": 20.0,
                "Na": 30.0,
                "Cl": 25.0,
                "SO4": 80.0,
                "HCO3": 220.0,
                "NO3": 3.0,
                "F": 0.4,
            },
        ]
        
        config = default_config()
        config.phreeqc_enabled = True
        config.missing_policy = "impute_zero"
        
        results = fit_network(samples, [("A", "B")], config)
        
        self.assertEqual(len(results), 1)
        result = results[0]
        
        self.assertEqual(result.edge_id, "A->B")
        self.assertIn(result.transport_model, ["evap", "mix"])
        self.assertIsNotNone(result.objective_score)
        self.assertEqual(len(result.z_extents), len(result.z_labels))

    def test_network_summary(self):
        """Test network summary statistics."""
        samples = [
            {"site_id": "A", "pH": 7.0, "Ca": 50, "Mg": 15, "Na": 25, 
             "Cl": 20, "SO4": 60, "HCO3": 180, "NO3": 5, "F": 0.3},
            {"site_id": "B", "pH": 7.5, "Ca": 65, "Mg": 20, "Na": 30,
             "Cl": 25, "SO4": 80, "HCO3": 220, "NO3": 3, "F": 0.4},
            {"site_id": "C", "pH": 7.8, "Ca": 80, "Mg": 25, "Na": 40,
             "Cl": 35, "SO4": 100, "HCO3": 260, "NO3": 2, "F": 0.5},
        ]
        
        config = default_config()
        config.phreeqc_enabled = False  # Faster for this test
        config.missing_policy = "impute_zero"
        
        results = fit_network(samples, [("A", "B"), ("B", "C")], config)
        summary = summarize_network(results)
        
        self.assertEqual(summary["edge_count"], 2)
        self.assertIn("evap", summary["transport_counts"] or {})
        self.assertIn("reaction_means", summary)
        self.assertIn("reaction_nonzero", summary)


class EdgeFitDetailedTests(unittest.TestCase):
    """Detailed edge fitting tests with various configurations."""

    def test_edge_fit_returns_all_fields(self):
        """Test that edge fit returns all expected fields."""
        config = default_config()
        
        x_u = [1.0, 0.5, 0.3, 0.5, 0.1, 0.6, 0.2, 0.1]
        x_v = [1.3, 0.7, 0.4, 0.7, 0.15, 0.8, 0.15, 0.12]
        
        result = fit_edge(x_u, x_v, config, edge_id="test", u="A", v="B")
        
        # Check all expected fields exist
        self.assertEqual(result.edge_id, "test")
        self.assertEqual(result.u, "A")
        self.assertEqual(result.v, "B")
        self.assertIn(result.transport_model, ["evap", "mix"])
        self.assertIsNotNone(result.z_extents)
        self.assertIsNotNone(result.z_labels)
        self.assertIsNotNone(result.transport_residual_norm)
        self.assertIsNotNone(result.anomaly_norm)
        self.assertIsNotNone(result.objective_score)
        self.assertIsNotNone(result.transport_probabilities)
        self.assertAlmostEqual(sum(result.transport_probabilities.values()), 1.0, places=6)

    def test_edge_fit_with_constraints(self):
        """Test edge fitting with PHREEQC-derived constraints."""
        config = default_config()
        matrix, labels, _ = build_reaction_dictionary(config)
        
        x_u = [1.0, 0.5, 0.3, 0.5, 0.1, 0.6, 0.2, 0.1]
        x_v = [1.5, 0.8, 0.4, 0.8, 0.15, 0.9, 0.18, 0.12]
        
        # Simulate bounds from PHREEQC (dissolution only for calcite)
        calcite_idx = labels.index("calcite")
        lb = [0.0] * len(labels)
        ub = [float("inf")] * len(labels)
        
        bounds = {
            "lb": lb,
            "ub": ub,
            "constraints_active": {"calcite": "dissolution_only"},
            "phreeqc_ok": True,
        }
        
        result = fit_edge(x_u, x_v, config, edge_id="test", u="A", v="B", bounds=bounds)
        
        # Calcite should only dissolve (positive or zero extent)
        self.assertGreaterEqual(result.z_extents[calcite_idx], -1e-6)


class DatabaseCompatibilityTests(unittest.TestCase):
    """Test different PHREEQC databases."""

    def test_all_databases_work(self):
        """Test that all available databases load and calculate SI."""
        db_dir = Path(__file__).parent.parent / "databases"
        databases = ["phreeqc.dat", "wateq4f.dat", "minteq.dat"]
        
        sample = [{
            "sample_id": "test",
            "pH": 7.0,
            "temp_c": 25.0,
            "Ca": 50.0,
            "Mg": 15.0,
            "Na": 25.0,
            "Cl": 20.0,
            "SO4": 60.0,
            "HCO3": 180.0,
        }]
        
        for db_name in databases:
            db_path = db_dir / db_name
            if not db_path.exists():
                self.skipTest(f"Database {db_name} not found")
            
            config = Config(phreeqc_database=str(db_path))
            results = run_phreeqc(sample, config)
            
            msg = f"Database {db_name} failed"
            self.assertTrue(results["test"]["phreeqc_ok"], msg)
            self.assertIsNotNone(results["test"]["si_calcite"], msg)


if __name__ == "__main__":
    unittest.main(verbosity=2)
