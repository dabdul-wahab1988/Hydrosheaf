import unittest
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
from hydrosheaf.config import Config
from hydrosheaf.nitrate_source_v2 import infer_node_posteriors
from hydrosheaf.models import nitrate_isotopes


class TestNitrateIsotopes(unittest.TestCase):
    def test_endmembers_loading(self):
        """Verify we can load the default JSON."""
        sources = nitrate_isotopes.load_isotope_endmembers()
        self.assertTrue(len(sources) >= 2)
        manure = next((s for s in sources if s.name == "Manure"), None)
        self.assertIsNotNone(manure)
        self.assertAlmostEqual(manure.d15N_mean, 15.0)

    def test_mixing_prob_manure(self):
        """Test a clear manure sample (High N15, Low O18)."""
        sources = nitrate_isotopes.load_isotope_endmembers()
        # Manure: d15N=15, d18O=5
        sample = nitrate_isotopes.IsotopeSample(d15N=15.0, d18O=5.0)
        probs = nitrate_isotopes.compute_isotope_prob(sample, sources)

        self.assertTrue(
            probs["Manure"] > 0.8, f"Manure should be dominant, got {probs}"
        )

    def test_mixing_prob_fertilizer(self):
        """Test a clear fertilizer sample (Low N15, High O18)."""
        sources = nitrate_isotopes.load_isotope_endmembers()
        # Fertilizer (NO3): d15N=0, d18O=20
        sample = nitrate_isotopes.IsotopeSample(d15N=0.0, d18O=20.0)
        probs = nitrate_isotopes.compute_isotope_prob(sample, sources)

        self.assertTrue(
            probs["Fertilizer"] > 0.8, f"Fertilizer should be dominant, got {probs}"
        )

    def test_covariance_likelihood_prefers_matching_correlation(self):
        """Correlation structure should influence source posterior probability."""
        sources = [
            nitrate_isotopes.SourceIsotopes(
                name="PosCorr",
                d15N_mean=0.0,
                d15N_std=2.0,
                d18O_mean=0.0,
                d18O_std=2.0,
                d15N_d18O_corr=0.9,
            ),
            nitrate_isotopes.SourceIsotopes(
                name="NegCorr",
                d15N_mean=0.0,
                d15N_std=2.0,
                d18O_mean=0.0,
                d18O_std=2.0,
                d15N_d18O_corr=-0.9,
            ),
        ]
        same_sign = nitrate_isotopes.compute_isotope_prob(
            nitrate_isotopes.IsotopeSample(d15N=1.5, d18O=1.5), sources
        )
        opposite_sign = nitrate_isotopes.compute_isotope_prob(
            nitrate_isotopes.IsotopeSample(d15N=1.5, d18O=-1.5), sources
        )
        self.assertGreater(same_sign["PosCorr"], same_sign["NegCorr"])
        self.assertGreater(opposite_sign["NegCorr"], opposite_sign["PosCorr"])

    def test_endmember_loader_reads_covariance_fields(self):
        """JSON parser should map optional covariance and correlation fields."""
        payload = {
            "sources": {
                "WithCorr": {
                    "d15N": {"mean": 5.0, "std": 2.0},
                    "d18O": {"mean": 8.0, "std": 4.0},
                    "correlation": {"d15N_d18O": 0.5},
                },
                "WithCov": {
                    "d15N": {"mean": 3.0, "std": 1.0},
                    "d18O": {"mean": 7.0, "std": 2.0},
                    "covariance": {"d15N_d18O": 0.8},
                },
            }
        }
        with tempfile.TemporaryDirectory() as tmp_dir:
            db_path = Path(tmp_dir) / "endmembers.json"
            with open(db_path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle)

            sources = nitrate_isotopes.load_isotope_endmembers(db_path)
            by_name = {s.name: s for s in sources}
            self.assertAlmostEqual(by_name["WithCorr"].covariance_d15N_d18O(), 4.0)
            self.assertAlmostEqual(by_name["WithCov"].covariance_d15N_d18O(), 0.8)

    def test_process_prior_denitrification_trajectory(self):
        """Denitrification trajectory should favor sources aligned with enrichment ray."""
        sample = nitrate_isotopes.IsotopeSample(d15N=10.0, d18O=5.0)
        sources = [
            nitrate_isotopes.SourceIsotopes(
                name="Aligned",
                d15N_mean=5.0,
                d15N_std=2.0,
                d18O_mean=2.0,
                d18O_std=2.0,
            ),
            nitrate_isotopes.SourceIsotopes(
                name="OffRay",
                d15N_mean=5.0,
                d15N_std=2.0,
                d18O_mean=12.0,
                d18O_std=2.0,
            ),
        ]
        prior_probs, flags, diagnostics = nitrate_isotopes.compute_process_prior_probs(
            sample=sample,
            sources=sources,
            denitrification_extent=1.0,
        )
        self.assertGreater(prior_probs["Aligned"], prior_probs["OffRay"])
        self.assertIn("denitrification_trajectory_match", flags)
        self.assertIn("denitrification_best_score", diagnostics)

    def test_process_prior_nitrification_pathway(self):
        """Nitrification pathway should prefer sources near expected d18O_NO3."""
        sample = nitrate_isotopes.IsotopeSample(d15N=4.0, d18O=4.0)
        sources = [
            nitrate_isotopes.SourceIsotopes(
                name="SoilLike",
                d15N_mean=5.0,
                d15N_std=2.0,
                d18O_mean=4.0,
                d18O_std=2.0,
            ),
            nitrate_isotopes.SourceIsotopes(
                name="AtmosphericLike",
                d15N_mean=0.0,
                d15N_std=2.0,
                d18O_mean=40.0,
                d18O_std=6.0,
            ),
        ]
        prior_probs, flags, diagnostics = nitrate_isotopes.compute_process_prior_probs(
            sample=sample,
            sources=sources,
            water_d18O=-6.0,
        )
        self.assertGreater(prior_probs["SoilLike"], prior_probs["AtmosphericLike"])
        self.assertIn("nitrification_pathway_match", flags)
        self.assertIn("nitrification_expected_d18O", diagnostics)

    def test_qc_low_identifiability_flag(self):
        """QC should flag weak source separation in near-uniform posteriors."""
        sources = [
            nitrate_isotopes.SourceIsotopes("A", 5.0, 2.0, 5.0, 2.0),
            nitrate_isotopes.SourceIsotopes("B", 5.0, 2.0, 5.0, 2.0),
        ]
        sample = nitrate_isotopes.IsotopeSample(d15N=5.0, d18O=5.0)
        diagnostics, flags = nitrate_isotopes.compute_isotope_qc_diagnostics(
            sample=sample,
            sources=sources,
            source_probs={"A": 0.5, "B": 0.5},
            qc_config={"min_top_probability": 0.60, "min_top_margin": 0.2},
        )
        self.assertIn("qc_low_identifiability", flags)
        self.assertIn("normalized_entropy", diagnostics)

    def test_qc_posterior_predictive_mismatch_flag(self):
        """QC should flag observations that are very unlikely under inferred source cloud."""
        sources = [
            nitrate_isotopes.SourceIsotopes("A", 0.0, 1.0, 0.0, 1.0),
            nitrate_isotopes.SourceIsotopes("B", 1.0, 1.0, 1.0, 1.0),
        ]
        sample = nitrate_isotopes.IsotopeSample(d15N=18.0, d18O=18.0)
        diagnostics, flags = nitrate_isotopes.compute_isotope_qc_diagnostics(
            sample=sample,
            sources=sources,
            source_probs={"A": 1.0, "B": 0.0},
            qc_config={"min_tail_probability": 0.05},
        )
        self.assertIn("qc_posterior_predictive_mismatch", flags)
        self.assertLess(diagnostics["posterior_predictive_tail_probability"], 0.05)

    def test_integration_with_inference(self):
        """Test that infer_node_posteriors uses the isotope logic."""
        # 1. Create a mock dataframe with isotope columns
        df = pd.DataFrame(
            [
                {
                    "site_id": "Well_A",
                    "NO3": 50.0,
                    "Cl": 1.0,
                    "d15N": 15.0,
                    "d18O_NO3": 5.0,
                    "d2H": -10.0,
                    "d18O": -3.0,
                },
                {
                    "site_id": "Well_B",
                    "NO3": 50.0,
                    "Cl": 1.0,
                    # Missing isotopes -> Fallback
                    "d2H": -10.0,
                    "d18O": -3.0,
                },
                # Add background samples to shift median
                {"site_id": "Bg_1", "NO3": 1.0, "Cl": 1.0, "d2H": -10, "d18O": -3},
                {"site_id": "Bg_2", "NO3": 2.0, "Cl": 1.0, "d2H": -10, "d18O": -3},
                {"site_id": "Bg_3", "NO3": 1.5, "Cl": 1.0, "d2H": -10, "d18O": -3},
            ]
        )
        df.set_index("site_id", inplace=True)

        # 2. Run inference
        results = infer_node_posteriors(df, edge_results=[])

        # 3. Check Well_A (Isotope)
        res_a = results["Well_A"]
        self.assertEqual(res_a.reason_code, "Dual Isotope Mixing")
        self.assertGreater(res_a.p_manure, 0.8)
        self.assertIn("dual_isotope_priority", res_a.gating_flags)

        # 4. Check Well_B (Fallback)
        res_b = results["Well_B"]
        self.assertEqual(res_b.reason_code, "Hydrochemical Ratios (No Isotopes)")
        # With High NO3/Cl (50) vs Median (~1.5), Z-score > 0
        # High Ratio => Fertilizer => Low p_manure
        self.assertLess(res_b.p_manure, 0.5)

    def test_integration_process_constraints_flag(self):
        """Inference should expose process-constraint gating flags when applied."""
        df = pd.DataFrame(
            [
                {
                    "site_id": "Well_A",
                    "NO3": 50.0,
                    "Cl": 1.0,
                    "d15N": 12.0,
                    "d18O_NO3": 7.0,
                    "d18O": -2.0,
                    "d2H": -10.0,
                },
                {"site_id": "Bg_1", "NO3": 1.0, "Cl": 1.0, "d2H": -10, "d18O": -3},
            ]
        ).set_index("site_id", drop=False)

        edge_results = [
            SimpleNamespace(v="Well_A", z_labels=["denit"], z_extents=[1.0]),
        ]
        results = infer_node_posteriors(df, edge_results=edge_results)
        res = results["Well_A"]

        self.assertEqual(res.reason_code, "Dual Isotope Mixing")
        self.assertIn("process_constraints_applied", res.gating_flags)
        self.assertIsNotNone(res.diagnostics)
        self.assertIn("top_probability", res.diagnostics)

    def test_integration_hierarchical_mcmc_flag(self):
        """Hierarchical MCMC mode should annotate isotope results with a hierarchical flag."""
        df = pd.DataFrame(
            [
                {
                    "site_id": "Well_A",
                    "NO3": 50.0,
                    "Cl": 1.0,
                    "d15N": 12.0,
                    "d18O_NO3": 7.0,
                    "d18O": -2.0,
                },
                {
                    "site_id": "Well_B",
                    "NO3": 45.0,
                    "Cl": 1.2,
                    "d15N": 5.0,
                    "d18O_NO3": 15.0,
                    "d18O": -3.0,
                },
                {"site_id": "Bg_1", "NO3": 1.0, "Cl": 1.0, "d18O": -3.0},
            ]
        ).set_index("site_id", drop=False)

        cfg = Config()
        cfg.isotope_mcmc_enabled = True
        cfg.isotope_mcmc_hierarchical_enabled = True
        cfg.isotope_mcmc_n_samples = 50
        cfg.isotope_mcmc_n_chains = 2
        cfg.isotope_mcmc_warmup = 20

        results = infer_node_posteriors(df, edge_results=[], config=cfg)
        res_a = results["Well_A"]
        res_b = results["Well_B"]

        self.assertEqual(res_a.reason_code, "MCMC Bayesian Isotope Mixing")
        self.assertEqual(res_b.reason_code, "MCMC Bayesian Isotope Mixing")
        self.assertIn("hierarchical_mcmc", res_a.gating_flags)
        self.assertIn("hierarchical_mcmc", res_b.gating_flags)


if __name__ == "__main__":
    unittest.main()
