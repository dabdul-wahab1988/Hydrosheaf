import unittest
from unittest import mock

from hydrosheaf.calibration.cli import _resolve_internal_parameters
from hydrosheaf.calibration.definitions import AdjustableParameter


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


if __name__ == "__main__":
    unittest.main()
