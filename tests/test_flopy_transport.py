"""
Tests for FloPy 1D Transport Model.
"""

import unittest
from unittest import mock
import numpy as np

from hydrosheaf.transport.flopy_1d import (
    TransportResult,
    check_flopy_available,
    analytical_1d_transport,
    run_analytical_1d_transport,
)
from hydrosheaf.transport.coupling import (
    VadoseCouplingResult,
    couple_vadose_saturated,
    estimate_saturated_travel_time,
    estimate_denitrification_attenuation,
)
from hydrosheaf.config import Config


class AnalyticalTransportTests(unittest.TestCase):
    """Tests for analytical 1D transport solutions."""

    def test_analytical_solution_zero_time(self):
        """Test analytical solution at t=0."""
        c = analytical_1d_transport(
            x=10.0,
            t=0.0,
            c0=1.0,
            v=0.1,
            D=0.1,
        )
        self.assertEqual(c, 0.0)

    def test_analytical_solution_far_field(self):
        """Test that concentration is low far from source at early times."""
        c = analytical_1d_transport(
            x=100.0,  # Far from source
            t=10.0,  # Early time
            c0=1.0,
            v=1.0,
            D=1.0,
        )
        # Front position is v*t = 10m, x=100m is 90m ahead
        # With 6*sqrt(D*t) = 6*sqrt(10) ≈ 19m, x is beyond threshold
        self.assertLessEqual(c, 0.01)

    def test_analytical_solution_concentration_bounded(self):
        """Test that concentration doesn't exceed source."""
        c = analytical_1d_transport(
            x=10.0,
            t=1000.0,  # Long time
            c0=1.0,
            v=0.1,
            D=0.1,
        )
        self.assertLessEqual(c, 1.0)
        self.assertGreaterEqual(c, 0.0)

    def test_analytical_solution_with_decay(self):
        """Test that decay reduces concentration."""
        c_no_decay = analytical_1d_transport(
            x=50.0,
            t=500.0,
            c0=1.0,
            v=0.1,
            D=0.1,
            k=0.0,
        )

        c_with_decay = analytical_1d_transport(
            x=50.0,
            t=500.0,
            c0=1.0,
            v=0.1,
            D=0.1,
            k=0.001,
        )

        self.assertLess(c_with_decay, c_no_decay)

    def test_analytical_breakthrough_curve(self):
        """Test analytical breakthrough curve generation."""
        times = np.linspace(0, 1000, 100)

        result = run_analytical_1d_transport(
            times=times,
            distance_m=50.0,
            velocity_m_day=0.1,
            dispersivity_m=1.0,
            decay_rate_1_day=0.0,
            source_concentration=1.0,
        )

        self.assertIsInstance(result, TransportResult)
        self.assertEqual(len(result.times), 100)
        self.assertEqual(len(result.concentrations), 100)
        self.assertTrue(result.success)

        # Concentration should increase and then plateau
        self.assertEqual(result.concentrations[0], 0.0)
        self.assertGreater(result.concentrations[-1], 0.5)


class VadoseCouplingTests(unittest.TestCase):
    """Tests for vadose-saturated coupling."""

    def test_coupling_basic(self):
        """Test basic vadose-saturated coupling."""
        # Create simple vadose output
        vadose_times = np.linspace(0, 365, 50)
        vadose_concentration = np.sin(vadose_times / 365 * np.pi) * 0.5 + 0.5

        result = couple_vadose_saturated(
            vadose_times=vadose_times,
            vadose_concentration=vadose_concentration,
            use_analytical=True,
            aquifer_length_m=100.0,
            hydraulic_k_m_day=1.0,
            porosity=0.25,
            dispersivity_m=1.0,
            denitrification_k_1_day=0.0,
            head_gradient=0.01,
        )

        self.assertIsInstance(result, VadoseCouplingResult)
        self.assertTrue(result.success)
        self.assertIsNotNone(result.combined_times)
        self.assertIsNotNone(result.combined_concentration)

    def test_coupling_with_decay(self):
        """Test that denitrification reduces output concentration via analytical function."""
        # Test directly with analytical function at steady state
        from hydrosheaf.transport.flopy_1d import analytical_1d_transport

        # Long time, concentration should reach steady state
        t = 1000.0
        x = 50.0
        v = 0.1
        D = 1.0
        c0 = 1.0

        # Without decay
        c_no_decay = analytical_1d_transport(x=x, t=t, c0=c0, v=v, D=D, k=0.0)

        # With decay
        c_with_decay = analytical_1d_transport(x=x, t=t, c0=c0, v=v, D=D, k=0.001)

        # Decay should reduce concentration
        self.assertLess(c_with_decay, c_no_decay)

    def test_coupling_with_config(self):
        """Test coupling using Config object."""
        config = Config()
        config.aquifer_thickness_m = 15.0
        config.aquifer_porosity = 0.3
        config.dispersivity_m = 2.0
        config.denitrification_k_1_day = 0.005

        vadose_times = np.linspace(0, 365, 50)
        vadose_concentration = np.ones(50) * 1.0

        result = couple_vadose_saturated(
            vadose_times=vadose_times,
            vadose_concentration=vadose_concentration,
            config=config,
            use_analytical=True,
        )

        self.assertTrue(result.success)

    def test_travel_time_estimation(self):
        """Test saturated zone travel time estimation."""
        travel_time = estimate_saturated_travel_time(
            distance_m=100.0,
            hydraulic_k_m_day=1.0,
            head_gradient=0.01,
            porosity=0.25,
        )

        # Expected: 100 / (1.0 * 0.01 / 0.25) = 100 / 0.04 = 2500 days
        self.assertAlmostEqual(travel_time, 2500.0, places=0)

    def test_denitrification_attenuation(self):
        """Test denitrification attenuation calculation."""
        # No decay
        atten = estimate_denitrification_attenuation(
            travel_time_days=100.0,
            decay_rate_1_day=0.0,
        )
        self.assertEqual(atten, 1.0)

        # With decay
        atten = estimate_denitrification_attenuation(
            travel_time_days=100.0,
            decay_rate_1_day=0.01,
        )
        # exp(-0.01 * 100) = exp(-1) ≈ 0.368
        self.assertAlmostEqual(atten, np.exp(-1), places=3)


class FloPyAvailabilityTests(unittest.TestCase):
    """Tests for FloPy availability checking."""

    def test_flopy_availability_check(self):
        """Test that FloPy availability check works."""
        available = check_flopy_available()
        self.assertIsInstance(available, bool)

    def test_build_model_raises_when_flopy_unavailable(self):
        """Guard rail: transport model build should fail clearly when FloPy is unavailable."""
        from hydrosheaf.transport import flopy_1d

        with mock.patch.object(flopy_1d, "FLOPY_AVAILABLE", False):
            with self.assertRaises(ImportError):
                flopy_1d.build_1d_transport_model()

    def test_run_model_raises_when_flopy_unavailable(self):
        """Guard rail: transport run should fail clearly when FloPy is unavailable."""
        from hydrosheaf.transport import flopy_1d

        with mock.patch.object(flopy_1d, "FLOPY_AVAILABLE", False):
            with self.assertRaises(ImportError):
                flopy_1d.run_1d_transport(mf=object(), mt=object(), params={})


class TransportResultTests(unittest.TestCase):
    """Tests for TransportResult dataclass."""

    def test_transport_result_creation(self):
        """Test TransportResult dataclass."""
        result = TransportResult(
            times=np.array([0, 1, 2]),
            concentrations=np.array([0, 0.5, 0.8]),
            n_cells=50,
            velocity_m_day=0.1,
        )

        self.assertEqual(len(result.times), 3)
        self.assertEqual(result.n_cells, 50)
        self.assertTrue(result.success)


class ConfigValidationTests(unittest.TestCase):
    """Tests for transport-related config validation."""

    def test_valid_transport_config(self):
        """Test valid transport configuration."""
        config = Config()
        config.aquifer_thickness_m = 10.0
        config.aquifer_porosity = 0.25
        config.dispersivity_m = 1.0
        config.denitrification_k_1_day = 0.001

        # Should not raise
        config.validate()

    def test_invalid_porosity(self):
        """Test invalid porosity raises error."""
        config = Config()
        config.aquifer_porosity = 1.5  # Invalid: > 1

        with self.assertRaises(ValueError):
            config.validate()

    def test_invalid_thickness(self):
        """Test invalid aquifer thickness raises error."""
        config = Config()
        config.aquifer_thickness_m = -10.0  # Invalid: negative

        with self.assertRaises(ValueError):
            config.validate()

    def test_invalid_decay_rate(self):
        """Test invalid decay rate raises error."""
        config = Config()
        config.denitrification_k_1_day = -0.01  # Invalid: negative

        with self.assertRaises(ValueError):
            config.validate()


if __name__ == "__main__":
    unittest.main()
