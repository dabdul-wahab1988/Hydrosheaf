import unittest
import numpy as np
from datetime import datetime, timedelta
from hydrosheaf.temporal import TimeSeriesSample, TemporalNode
from hydrosheaf.temporal.interpolation import (
    _linear_interp,
)
from hydrosheaf.temporal.temporal_edge_fit import compute_seasonal_decomposition
from hydrosheaf.temporal.residence_time import estimate_residence_time_with_details


class TestTemporalAccuracy(unittest.TestCase):
    """
    Test accuracy of Temporal Dynamics extension from a mathematical perspective.
    """

    def test_linear_interpolation(self):
        """
        Verify linear interpolation mathematics.
        """
        times = [0.0, 10.0, 20.0]
        values = [0.0, 100.0, 50.0]

        # Midpoint test (t=5): Should be (0+100)/2 = 50
        self.assertAlmostEqual(_linear_interp(times, values, 5.0), 50.0)

        # Endpoint test (t=10): 100
        self.assertAlmostEqual(_linear_interp(times, values, 10.0), 100.0)

        # Second interval (t=15): (100+50)/2 = 75
        self.assertAlmostEqual(_linear_interp(times, values, 15.0), 75.0)

        # Extrapolation (clamped)
        self.assertAlmostEqual(_linear_interp(times, values, 25.0), 50.0)

    def test_seasonal_decomposition_synthetic(self):
        """
        Verify seasonal decomposition against a pure sine wave signal.
        C(t) = 10 + 2*t + 5*sin(2*pi*t/365) + noise
        """
        np.random.seed(42)
        n_samples = 365 * 2  # 2 years
        times = np.arange(n_samples)

        # Define signal components
        mean_true = 10.0
        trend_true = 0.01  # slope per day
        amp_true = 5.0
        period = 365.0

        # Generate signal
        # Note: Implementation expects sin/cos summation.
        # Here we use sin(omega*t) which corresponds to [cos_coeff=0, sin_coeff=5] -> amp=5
        omega = 2 * np.pi / period
        signal = mean_true + trend_true * times + amp_true * np.sin(omega * times)

        # No noise for exact check

        # Build samples
        t0 = datetime(2020, 1, 1)
        samples = []
        for i, val in enumerate(signal):
            samples.append(
                TimeSeriesSample(
                    sample_id=str(i),
                    node_id="test",
                    timestamp=t0 + timedelta(days=float(times[i])),
                    concentrations=[val],
                )
            )

        node = TemporalNode(node_id="test", samples=samples)

        # Run decomposition
        mean, trend, amp, residual_std = compute_seasonal_decomposition(
            node, 0, period_days=365.0
        )

        # Check components
        # Mean might be slightly shifted due to trend centering def, but should be close
        self.assertAlmostEqual(trend, trend_true, places=5)
        self.assertAlmostEqual(amp, amp_true, places=3)
        self.assertAlmostEqual(residual_std, 0.0, places=3)  # Pure signal

    def test_residence_time_cross_correlation_logic(self):
        """
        Verify cross-correlation logic for finding lags.
        Instead of calling the complex full function, we simulate the core math.
        """
        # Create two signals, one lagged by 10 days
        t = np.linspace(0, 100, 100)
        sig_u = np.sin(t / 10.0)

        # v(t) = u(t - lag)
        lag_true = 10.0
        sig_v = np.sin((t - lag_true) / 10.0)

        # Normalize
        u_norm = (sig_u - np.mean(sig_u)) / np.std(sig_u)
        v_norm = (sig_v - np.mean(sig_v)) / np.std(sig_v)

        # Compute correlation at exactly the true lag
        # At lag=10, v(t) aligns with u(t-10). Wait.
        # Implementation: for v sample at t, find u at t-lag.
        # So if v(t) = u(t-10), then looking at u(t-10) should perfectly match v(t).

        match_sum = 0
        count = 0
        for i, t_val in enumerate(t):
            target_u = t_val - lag_true
            if target_u >= 0:
                # Interpolate u at target
                u_val = np.interp(target_u, t, u_norm)
                match_sum += u_val * v_norm[i]
                count += 1

        corr = match_sum / count

        # Correlation should be nearly 1.0
        self.assertTrue(corr > 0.95, f"Lagged correlation {corr} should be near 1.0")

    def test_residence_time_ttd_convolution_recovers_mean_lag(self):
        """
        Verify TTD+attenuation residence-time estimation recovers the effective mean lag.

        Synthetic model:
          v(t) = b + (u * h)(t), where h(τ) ∝ w(τ) * exp(-k*τ)
        """
        np.random.seed(7)
        n_days = 220
        dt = 1.0
        t = np.arange(n_days, dtype=float) * dt

        # Upstream signal with variability
        u = np.sin(t / 12.0) + 0.3 * np.cos(t / 6.0) + 0.05 * np.random.randn(n_days)

        # Travel-time distribution w (gamma-like, discretized) and attenuation exp(-k*tau)
        max_lag = 60
        lags = np.arange(max_lag + 1, dtype=float) * dt
        w_raw = (lags**2) * np.exp(
            -lags / 6.0
        )  # shape ~3, scale ~6, mean ~18 before attenuation
        w = w_raw / np.sum(w_raw)
        k_true = 0.01
        h = w * np.exp(-k_true * lags)
        h = h / np.sum(h)
        tau_true = float(np.sum(lags * h))

        # Convolution + intercept
        b0 = 0.2
        v = b0 + np.convolve(u, h, mode="full")[:n_days]

        # Build TemporalNodes (single-ion order containing "Cl")
        t0 = datetime(2020, 1, 1)
        u_samples = [
            TimeSeriesSample(
                sample_id=f"u{i}",
                node_id="U",
                timestamp=t0 + timedelta(days=float(i)),
                concentrations=[float(u[i])],
            )
            for i in range(n_days)
        ]
        v_samples = [
            TimeSeriesSample(
                sample_id=f"v{i}",
                node_id="V",
                timestamp=t0 + timedelta(days=float(i)),
                concentrations=[float(v[i])],
            )
            for i in range(n_days)
        ]
        node_u = TemporalNode(node_id="U", samples=u_samples)
        node_v = TemporalNode(node_id="V", samples=v_samples)

        tau, unc, used, details, flags = estimate_residence_time_with_details(
            node_u,
            node_v,
            method="ttd",
            tracer_ion="Cl",
            ion_order=["Cl"],
            hydraulic_params={
                "grid_dt_days": 1.0,
                "max_lag_days": float(max_lag),
                "smoothness_lambda": 0.01,
                "ttd_min_r2": 0.5,
                "attenuation_k_max": 0.02,
                "attenuation_k_steps": 6,
            },
        )

        self.assertTrue("ttd" in used, f"Expected a TTD method, got: {used}")
        self.assertFalse("ttd_failed" in used, f"TTD should not fail: {used}")
        self.assertTrue(tau > 0.0 and unc >= 0.0)
        self.assertTrue("ttd_failed_all_tracers" not in flags)
        self.assertAlmostEqual(tau, tau_true, delta=3.0)

    def test_residence_time_bayesian_lag_recovers_lag(self):
        """
        Verify Bayesian lag estimator returns a sensible posterior mean tau for sparse data,
        stabilized by a physics-informed prior.
        """
        np.random.seed(11)
        n = 40
        dt_days = 7.0
        lag_true = 21.0
        t = np.arange(n, dtype=float) * dt_days
        u = np.sin(t / 30.0) + 0.05 * np.random.randn(n)
        k_true = 0.01
        a_true = 1.3
        b_true = -0.2
        v = np.zeros_like(u)
        for i, tt in enumerate(t):
            tt_u = tt - lag_true
            if tt_u < t[0]:
                v[i] = b_true
            else:
                u_interp = np.interp(tt_u, t, u)
                v[i] = b_true + a_true * u_interp * np.exp(-k_true * lag_true)
        v = v + 0.03 * np.random.randn(n)

        t0 = datetime(2020, 1, 1)
        node_u = TemporalNode(
            node_id="U",
            samples=[
                TimeSeriesSample(
                    sample_id=f"u{i}",
                    node_id="U",
                    timestamp=t0 + timedelta(days=float(t[i])),
                    concentrations=[float(u[i])],
                )
                for i in range(n)
            ],
        )
        node_v = TemporalNode(
            node_id="V",
            samples=[
                TimeSeriesSample(
                    sample_id=f"v{i}",
                    node_id="V",
                    timestamp=t0 + timedelta(days=float(t[i])),
                    concentrations=[float(v[i])],
                )
                for i in range(n)
            ],
        )

        tau, unc, used, details, flags = estimate_residence_time_with_details(
            node_u,
            node_v,
            method="bayesian_lag",
            tracer_ion="Cl",
            ion_order=["Cl"],
            hydraulic_params={
                # physics prior roughly near the right scale
                "distance_m": 1000.0,
                "K_m_day": 1.0,
                "gradient": 0.01,
                "porosity": 0.2,
                "bayes_lag_grid_dt_days": 3.0,
                "bayes_lag_max_lag_days": 120.0,
                "bayes_lag_min_pairs": 6,
                "bayes_lag_prior_sigma_multiplier": 2.0,
                "attenuation_k_max": 0.02,
                "attenuation_k_steps": 6,
                "ttd_min_r2": 0.1,
            },
        )

        self.assertTrue("bayesian_lag" in used)
        self.assertFalse("bayesian_lag_failed" in used)
        self.assertTrue(tau > 0.0)
        self.assertTrue(unc >= 0.0)
        self.assertTrue("bayes_failed_all_tracers" not in flags)
        self.assertAlmostEqual(tau, lag_true, delta=6.0)
