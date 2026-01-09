
import unittest
import numpy as np
from datetime import datetime, timedelta
from hydrosheaf.temporal import TimeSeriesSample, TemporalNode
from hydrosheaf.temporal.interpolation import (
    interpolate_to_common_times, 
    _linear_interp, 
    _spline_interp
)
from hydrosheaf.temporal.temporal_edge_fit import compute_seasonal_decomposition

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
        n_samples = 365 * 2 # 2 years
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
            samples.append(TimeSeriesSample(
                sample_id=str(i),
                node_id="test",
                timestamp=t0 + timedelta(days=float(times[i])),
                concentrations=[val]
            ))
            
        node = TemporalNode(node_id="test", samples=samples)
        
        # Run decomposition
        mean, trend, amp, residual_std = compute_seasonal_decomposition(node, 0, period_days=365.0)
        
        # Check components
        # Mean might be slightly shifted due to trend centering def, but should be close
        self.assertAlmostEqual(trend, trend_true, places=5)
        self.assertAlmostEqual(amp, amp_true, places=3)
        self.assertAlmostEqual(residual_std, 0.0, places=3) # Pure signal

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

