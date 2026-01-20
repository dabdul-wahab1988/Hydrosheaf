import unittest
from unittest import mock

from hydrosheaf.api import (
    auto_disable_missing_modules,
    fit_network_pipeline,
    fit_network_with_priors,
)
from hydrosheaf.config import Config


class AutoDisableMissingModulesTests(unittest.TestCase):
    def test_disables_all_when_data_missing(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A"}]

        updated = auto_disable_missing_modules(samples, config)

        self.assertFalse(updated.phreeqc_enabled)
        self.assertFalse(updated.isotope_enabled)
        self.assertFalse(updated.nitrate_source_enabled)

    def test_keeps_isotopes_when_present(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A", "18O": -5.0, "2H": -35.0}]

        updated = auto_disable_missing_modules(samples, config)

        self.assertFalse(updated.phreeqc_enabled)
        self.assertTrue(updated.isotope_enabled)
        self.assertFalse(updated.nitrate_source_enabled)

    def test_keeps_phreeqc_and_nitrate_when_present(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A", "pH": 7.2, "NO3": 1.5}]

        updated = auto_disable_missing_modules(samples, config)

        self.assertTrue(updated.phreeqc_enabled)
        self.assertFalse(updated.isotope_enabled)
        self.assertTrue(updated.nitrate_source_enabled)

    def test_fit_network_with_priors_auto_disable(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A"}, {"site_id": "B"}]
        edges = [("A", "B")]
        captured = {}

        def _fake_fit_network(samples_arg, edges_arg, config_arg, **kwargs):
            captured["config"] = config_arg
            return []

        with mock.patch("hydrosheaf.api.fit_network", side_effect=_fake_fit_network):
            fit_network_with_priors(samples, edges, config)

        self.assertIn("config", captured)
        self.assertFalse(captured["config"].phreeqc_enabled)
        self.assertFalse(captured["config"].isotope_enabled)
        self.assertFalse(captured["config"].nitrate_source_enabled)

    def test_fit_network_pipeline_auto_disable(self):
        config = Config(
            phreeqc_enabled=True, isotope_enabled=True, nitrate_source_enabled=True
        )
        samples = [{"site_id": "A"}, {"site_id": "B"}]
        edges = [("A", "B")]
        captured = {}

        def _fake_fit_network(samples_arg, edges_arg, config_arg, **kwargs):
            captured["config"] = config_arg
            return []

        with mock.patch("hydrosheaf.api.fit_network", side_effect=_fake_fit_network):
            results, extras = fit_network_pipeline(samples, edges, config)

        self.assertEqual(results, [])
        self.assertIn("config", captured)
        self.assertFalse(captured["config"].phreeqc_enabled)
        self.assertFalse(captured["config"].isotope_enabled)
        self.assertFalse(captured["config"].nitrate_source_enabled)
        self.assertIn("edges", extras)

    def test_fit_network_with_priors_strict_requires_phreeqc_major_ions(self):
        config = Config(phreeqc_enabled=True, strict_input_validation=True)
        samples = [{"site_id": "A", "pH": 7.1}, {"site_id": "B", "pH": 7.3}]
        edges = [("A", "B")]

        with mock.patch("hydrosheaf.api.fit_network") as mocked_fit:
            with self.assertRaises(ValueError):
                fit_network_with_priors(samples, edges, config)

        mocked_fit.assert_not_called()


if __name__ == "__main__":
    unittest.main()
