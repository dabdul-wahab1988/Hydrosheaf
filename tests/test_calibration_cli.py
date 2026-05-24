import unittest
from unittest import mock

from hydrosheaf.calibration.cli import _resolve_internal_parameters
from hydrosheaf.calibration.definitions import AdjustableParameter
from hydrosheaf.calibration.config import load_calibration_config

import tempfile
import yaml


class _DummyProblem:
    def __init__(self, parameters):
        self._parameters = parameters

    def get_parameters(self):
        return self._parameters


class CalibrationCliParameterResolutionTests(unittest.TestCase):
    def test_prefers_adapter_parameters_and_overrides_matching_config(self):
        problem_params = [
            AdjustableParameter(
                name="dispersivity",
                value=1.0,
                lower_bound=0.01,
                upper_bound=100.0,
                log_transform=True,
                prior_mean=1.0,
                prior_sigma=1.0,
            ),
            AdjustableParameter(name="decay", value=0.05, lower_bound=0.0, upper_bound=1.0),
        ]
        config_params = [
            AdjustableParameter(
                name="dispersivity",
                value=2.5,
                lower_bound=0.1,
                upper_bound=50.0,
                log_transform=False,
                prior_mean=2.0,
                prior_sigma=0.25,
            ),
            AdjustableParameter(name="unknown_param", value=9.9),
        ]
        logger = mock.MagicMock()

        resolved = _resolve_internal_parameters(
            _DummyProblem(problem_params), config_params, logger
        )

        self.assertEqual([param.name for param in resolved], ["dispersivity", "decay"])
        self.assertEqual(resolved[0].value, 2.5)
        self.assertEqual(resolved[0].lower_bound, 0.1)
        self.assertEqual(resolved[0].upper_bound, 50.0)
        self.assertFalse(resolved[0].log_transform)
        self.assertEqual(resolved[0].prior_mean, 2.0)
        self.assertEqual(resolved[0].prior_sigma, 0.25)
        self.assertEqual(resolved[1].value, 0.05)

        # Ensure original adapter parameters are not mutated in place.
        self.assertEqual(problem_params[0].value, 1.0)

        self.assertTrue(logger.warning.called)

    def test_falls_back_to_config_parameters_when_adapter_has_none(self):
        config_params = [
            AdjustableParameter(name="k", value=1.0, lower_bound=0.0, upper_bound=2.0)
        ]
        logger = mock.MagicMock()

        resolved = _resolve_internal_parameters(
            _DummyProblem([]), config_params, logger
        )

        self.assertEqual(len(resolved), 1)
        self.assertEqual(resolved[0].name, "k")
        self.assertEqual(resolved[0].value, 1.0)
        self.assertIsNot(resolved[0], config_params[0])
        self.assertTrue(logger.warning.called)

    def test_fixed_and_tied_fields_preserved_in_resolution(self):
        """Verify that fixed and tied_to fields propagate during parameter resolution."""
        problem_params = [
            AdjustableParameter(name="p1", value=1.0, lower_bound=0.0, upper_bound=10.0),
            AdjustableParameter(name="p2", value=2.0, lower_bound=0.0, upper_bound=10.0),
        ]
        config_params = [
            AdjustableParameter(name="p1", value=1.5, lower_bound=0.0, upper_bound=10.0, fixed=True),
            AdjustableParameter(name="p2", value=2.5, lower_bound=0.0, upper_bound=10.0, tied_to="p1"),
        ]
        logger = mock.MagicMock()

        resolved = _resolve_internal_parameters(
            _DummyProblem(problem_params), config_params, logger
        )

        self.assertEqual(len(resolved), 2)
        self.assertTrue(resolved[0].fixed)
        self.assertIsNone(resolved[0].tied_to)
        self.assertFalse(resolved[1].fixed)
        self.assertEqual(resolved[1].tied_to, "p1")

    def test_yaml_fixed_tied_parsing(self):
        """Verify that fixed and tied_to are parsed from calibration YAML."""
        config_data = {
            "calibration": {
                "type": "transport",
                "settings": {"n_workers": 1, "max_iterations": 10},
                "parameters": [
                    {"name": "dispersivity", "initial": 1.0, "bounds": [0.1, 10.0], "log": True},
                    {"name": "decay", "initial": 0.0, "bounds": [0.0, 1.0], "fixed": True},
                    {"name": "velocity", "initial": 0.1, "bounds": [0.01, 1.0], "tied_to": "dispersivity"},
                ]
            }
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            yaml.dump(config_data, f)
            f.flush()
            config = load_calibration_config(f.name)

        names = {p.name: p for p in config.parameters}
        self.assertIn("dispersivity", names)
        self.assertFalse(names["dispersivity"].fixed)
        self.assertIsNone(names["dispersivity"].tied_to)
        self.assertTrue(names["decay"].fixed)
        self.assertIsNone(names["decay"].tied_to)
        self.assertFalse(names["velocity"].fixed)
        self.assertEqual(names["velocity"].tied_to, "dispersivity")


if __name__ == "__main__":
    unittest.main()


class LoadCalibrationJsonTests(unittest.TestCase):
    def test_loads_direct_and_aliased_parameters(self):
        import tempfile
        import json
        from hydrosheaf.config import Config

        data = {
            "optimal_parameters": {
                "dispersivity": 2.5,
                "porosity": 0.3,
                "unknown_adapter_param": 999.0,
            }
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            f.flush()
            config = Config()
            config.dispersivity_m = 1.0  # current value before loading
            config.aquifer_porosity = 0.2
            config.load_from_calibration_json(f.name)

        self.assertAlmostEqual(config.dispersivity_m, 2.5)
        self.assertAlmostEqual(config.aquifer_porosity, 0.3)
        # unknown_adapter_param should be silently ignored
        self.assertFalse(hasattr(config, "unknown_adapter_param"))

    def test_loads_ensemble_result_via_posterior_means(self):
        import tempfile
        import json
        from hydrosheaf.config import Config

        data = {
            "success": True,
            "optimal_parameters": {},  # empty dict
            "posterior_parameters": {
                "porosity": [0.18, 0.22, 0.20],
                "decay": [0.001, 0.002, 0.0015],
            }
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            f.flush()
            config = Config()
            config.aquifer_porosity = 0.15
            config.denitrification_k_1_day = 0.0005
            config.load_from_calibration_json(f.name)

        # 0.18 + 0.22 + 0.20 = 0.60 / 3 = 0.20
        self.assertAlmostEqual(config.aquifer_porosity, 0.20, places=5)
        # 0.001 + 0.002 + 0.0015 = 0.0045 / 3 = 0.0015
        self.assertAlmostEqual(config.denitrification_k_1_day, 0.0015, places=5)

    def test_loads_kinetic_logk_mineral_rate(self):
        import tempfile
        import json
        from hydrosheaf.config import Config

        data = {
            "optimal_parameters": {
                "log_k_calcite": -6.5,
            }
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(data, f)
            f.flush()
            config = Config()
            config.mineral_rate_constants = {"calcite": -6.0}
            config.load_from_calibration_json(f.name)

        self.assertAlmostEqual(config.mineral_rate_constants["calcite"], -6.5)
