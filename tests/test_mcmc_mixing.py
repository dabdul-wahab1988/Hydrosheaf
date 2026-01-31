"""
Tests for MCMC Bayesian Isotope Mixing Model.
"""

import unittest
import sys

from hydrosheaf.models.nitrate_isotopes import IsotopeSample, SourceIsotopes
from hydrosheaf.models.nitrate_isotopes_mcmc import (
    MCMCMixingResult,
    check_pymc_available,
    run_mcmc_mixing,
    run_mcmc_mixing_batch,
    summarize_mcmc_results,
)


class MCMCMixingTests(unittest.TestCase):
    """Tests for MCMC isotope mixing."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        # Define source endmembers with realistic isotope values
        cls.sources = [
            SourceIsotopes(
                name="Manure",
                d15N_mean=15.0,
                d15N_std=3.0,
                d18O_mean=5.0,
                d18O_std=2.0,
            ),
            SourceIsotopes(
                name="Fertilizer",
                d15N_mean=0.0,
                d15N_std=2.0,
                d18O_mean=20.0,
                d18O_std=3.0,
            ),
            SourceIsotopes(
                name="Soil_N",
                d15N_mean=5.0,
                d15N_std=2.0,
                d18O_mean=2.0,
                d18O_std=2.0,
            ),
        ]

        # Sample close to manure signature
        cls.manure_sample = IsotopeSample(d15N=14.0, d18O=6.0)

        # Sample close to fertilizer signature
        cls.fertilizer_sample = IsotopeSample(d15N=1.0, d18O=18.0)

        # Mixed sample
        cls.mixed_sample = IsotopeSample(d15N=8.0, d18O=10.0)

    def test_pymc_availability_check(self):
        """Test that PyMC availability check works."""
        available = check_pymc_available()
        self.assertIsInstance(available, bool)

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_mixing_manure_dominant(self):
        """Test MCMC mixing with manure-like sample."""
        result = run_mcmc_mixing(
            sample=self.manure_sample,
            sources=self.sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        self.assertIsInstance(result, MCMCMixingResult)
        self.assertIn("Manure", result.source_fractions)
        self.assertIn("Fertilizer", result.source_fractions)

        # Manure should be dominant
        self.assertGreater(result.source_fractions["Manure"], 0.4)

        # Fractions should sum to 1
        total = sum(result.source_fractions.values())
        self.assertAlmostEqual(total, 1.0, places=2)

        # CI should be valid
        self.assertIn("Manure", result.ci_lower)
        self.assertIn("Manure", result.ci_upper)
        self.assertLess(result.ci_lower["Manure"], result.source_fractions["Manure"])
        self.assertGreater(result.ci_upper["Manure"], result.source_fractions["Manure"])

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_mixing_fertilizer_dominant(self):
        """Test MCMC mixing with fertilizer-like sample."""
        result = run_mcmc_mixing(
            sample=self.fertilizer_sample,
            sources=self.sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        # Fertilizer should be dominant
        self.assertGreater(result.source_fractions["Fertilizer"], 0.4)

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_mixing_diagnostics(self):
        """Test that MCMC diagnostics are computed."""
        result = run_mcmc_mixing(
            sample=self.mixed_sample,
            sources=self.sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        # R-hat should be computed
        self.assertIn("Manure", result.r_hat)

        # ESS should be computed
        self.assertIn("Manure", result.ess_bulk)

        # Posterior samples should be available
        self.assertIsNotNone(result.posterior_samples)
        self.assertEqual(result.posterior_samples.shape[1], len(self.sources))

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_convergence_check(self):
        """Test convergence checking logic."""
        result = run_mcmc_mixing(
            sample=self.mixed_sample,
            sources=self.sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        # With enough samples, should converge
        # R-hat should be close to 1
        for name, rh in result.r_hat.items():
            self.assertLess(rh, 1.1, f"R-hat for {name} too high: {rh}")

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_batch_processing(self):
        """Test batch processing of multiple samples."""
        samples = [self.manure_sample, self.fertilizer_sample, self.mixed_sample]

        results = run_mcmc_mixing_batch(
            samples=samples,
            sources=self.sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        self.assertEqual(len(results), 3)
        for result in results:
            self.assertIsInstance(result, MCMCMixingResult)

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_summary(self):
        """Test summarization of multiple MCMC results."""
        samples = [self.manure_sample, self.fertilizer_sample]

        results = run_mcmc_mixing_batch(
            samples=samples,
            sources=self.sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        summary = summarize_mcmc_results(results)

        self.assertIn("n_samples", summary)
        self.assertEqual(summary["n_samples"], 2)
        self.assertIn("mean_fractions", summary)
        self.assertIn("Manure", summary["mean_fractions"])

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_mcmc_with_two_sources(self):
        """Test MCMC works with minimum 2 sources."""
        two_sources = self.sources[:2]  # Just Manure and Fertilizer

        if not check_pymc_available():
            self.skipTest("PyMC not available")

        result = run_mcmc_mixing(
            sample=self.mixed_sample,
            sources=two_sources,
            n_samples=50,
            n_chains=2,
            warmup=20,
            random_seed=42,
        )


        self.assertEqual(len(result.source_fractions), 2)
        total = sum(result.source_fractions.values())
        self.assertAlmostEqual(total, 1.0, places=2)

    @unittest.skipIf(not check_pymc_available(), "PyMC not available")
    def test_insufficient_sources_raises(self):
        """Test that single source raises error."""
        if not check_pymc_available():
            self.skipTest("PyMC not available")

        with self.assertRaises(ValueError):
            run_mcmc_mixing(
                sample=self.mixed_sample,
                sources=[self.sources[0]],  # Only one source
                n_samples=50,
            )



class MCMCFallbackTests(unittest.TestCase):
    """Tests for MCMC fallback behavior when PyMC not available."""

    def test_import_without_pymc(self):
        """Test module imports without errors even if PyMC not available."""
        # This should not raise even if PyMC is not installed
        from hydrosheaf.models.nitrate_isotopes_mcmc import check_pymc_available

        available = check_pymc_available()
        self.assertIsInstance(available, bool)


if __name__ == "__main__":
    unittest.main()
