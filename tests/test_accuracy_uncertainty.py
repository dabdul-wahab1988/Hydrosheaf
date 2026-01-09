
import unittest
import numpy as np
from hydrosheaf.uncertainty.bootstrap import compute_bca_ci
from hydrosheaf.uncertainty.bayesian import compute_r_hat, compute_ess

class TestUncertaintyAccuracy(unittest.TestCase):
    """
    Test accuracy of Uncertainty Quantification extension from a statistical perspective.
    """

    def test_bca_ci_symmetric(self):
        """
        Test BCa confidence interval for a symmetric distribution (Normal).
        For Normal(0, 1), 95% CI should be approx [-1.96, 1.96].
        """
        np.random.seed(42)
        # Generate 1000 samples from Normal(0, 1)
        samples = np.random.normal(0, 1, 1000)
        point_estimate = 0.0
        
        ci_low, ci_high = compute_bca_ci(samples, point_estimate, alpha=0.05)
        
        # Check proximity to theoretical values
        self.assertAlmostEqual(ci_low, -1.96, delta=0.2)
        self.assertAlmostEqual(ci_high, 1.96, delta=0.2)

    def test_bca_ci_skewed(self):
        """
        Test BCa CI for a skewed distribution (LogNormal).
        BCa should handle skewness better than simple percentile.
        """
        np.random.seed(42)
        # LogNormal(0, 1)
        samples = np.random.lognormal(0, 1, 2000)
        point_estimate = np.exp(0.5) # Theoretical mean ? No, mode is 1. Median is 1. Mean is exp(0.5).
        # Let's just use the sample median as point estimate for stability
        point_estimate = np.median(samples)
        
        ci_low, ci_high = compute_bca_ci(samples, point_estimate, alpha=0.05)
        
        # Just verify it returns valid bounds within the range
        self.assertTrue(ci_low < point_estimate < ci_high)
        self.assertTrue(ci_low > 0) # LogNormal is positive

    def test_r_hat_convergence(self):
        """
        Test Gelman-Rubin R-hat statistic.
        """
        np.random.seed(42)
        
        # Case 1: Converged chains (all from same distribution)
        n_chains = 4
        n_samples = 1000
        chains_converged = np.random.normal(0, 1, (n_chains, n_samples))
        
        r_hat = compute_r_hat(chains_converged)
        
        # Should be very close to 1.0 (typically < 1.01)
        self.assertTrue(0.99 <= r_hat < 1.01, f"Converged R-hat {r_hat} should be ~1.0")
        
        # Case 2: Diverged chains (different means)
        chains_diverged = np.zeros((n_chains, n_samples))
        chains_diverged[0, :] = np.random.normal(0, 0.1, n_samples)
        chains_diverged[1, :] = np.random.normal(10, 0.1, n_samples) # Way off
        chains_diverged[2, :] = np.random.normal(0, 0.1, n_samples)
        chains_diverged[3, :] = np.random.normal(0, 0.1, n_samples)
        
        r_hat_div = compute_r_hat(chains_diverged)
        
        # Should be large
        self.assertTrue(r_hat_div > 1.5, f"Diverged R-hat {r_hat_div} should be > 1.5")

    def test_ess_calculation(self):
        """
        Test Effective Sample Size (ESS) calculation.
        """
        np.random.seed(42)
        
        # Case 1: Independent samples (White noise)
        # ESS should be close to N
        N = 1000
        samples_indep = np.random.normal(0, 1, N)
        
        ess = compute_ess(samples_indep)
        self.assertTrue(N * 0.8 < ess < N * 1.2, f"Independent ESS {ess} should be close to {N}")
        
        # Case 2: Highly correlated samples (Random walk / AR1)
        # x_t = 0.9 * x_{t-1} + noise
        samples_corr = [0.0]
        for _ in range(N-1):
            samples_corr.append(0.9 * samples_corr[-1] + np.random.normal(0, 1))
        samples_corr = np.array(samples_corr)
        
        ess_corr = compute_ess(samples_corr)
        
        # ESS should be much lower than N
        self.assertTrue(ess_corr < N * 0.2, f"Correlated ESS {ess_corr} should be low")

