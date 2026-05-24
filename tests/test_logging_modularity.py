from unittest import mock

import pytest

from hydrosheaf.log import get_logger
from hydrosheaf.calibration.config import CalibrationConfig
from hydrosheaf.calibration.factory import build_calibration_problem
from hydrosheaf.calibration.cli import setup_vadose_adapter


def test_get_logger_does_not_double_prefix_package_names():
    assert get_logger("api").name == "hydrosheaf.api"
    assert get_logger("hydrosheaf.api").name == "hydrosheaf.api"
    assert get_logger("hydrosheaf.nuclear.invert").name == "hydrosheaf.nuclear.invert"


def test_calibration_cli_reexports_factory_setup_helpers():
    assert callable(setup_vadose_adapter)


def test_build_calibration_problem_rejects_unknown_type():
    logger = mock.MagicMock()
    config = CalibrationConfig(problem_type="unknown")

    with pytest.raises(ValueError, match="Unknown problem type"):
        build_calibration_problem(config, logger)

    logger.error.assert_called()
