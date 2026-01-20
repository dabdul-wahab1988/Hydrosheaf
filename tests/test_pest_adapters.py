"""
Tests for new PESTGLM Adapters (Vadose & Transport).
"""

import numpy as np
from unittest.mock import MagicMock, patch

from hydrosheaf.calibration.glm import PESTGLM
from hydrosheaf.calibration.adapters import (
    TransportCalibrationAdapter,
    TransportExperiment,
)


def test_transport_adapter_structure():
    """Test TransportAdapter param generation."""
    exp = TransportExperiment(
        id="col1",
        times=[10, 20, 30],
        observed_concentrations=[0.1, 0.5, 0.9],
        distance_m=10.0,
    )

    adapter = TransportCalibrationAdapter(
        experiments=[exp], params_to_fit=["dispersivity", "velocity"]
    )

    params = adapter.get_parameters()
    names = {p.name for p in params}

    assert "dispersivity" in names
    assert "velocity" in names
    assert "decay" not in names

    obs = adapter.get_observations()
    assert len(obs) == 3
    assert obs[0].name == "col1_0"
    assert obs[0].value == 0.1


@patch("hydrosheaf.transport.flopy_1d.run_analytical_1d_transport")
def test_transport_adapter_execution(mock_run):
    """Test TransportAdapter model execution."""
    # Setup mock return
    # Needs to return a TransportResult object (simulated)
    mock_result = MagicMock()
    mock_result.concentrations = np.array([0.15, 0.55, 0.95])
    mock_run.return_value = mock_result

    exp = TransportExperiment(
        id="col1",
        times=[10, 20, 30],
        observed_concentrations=[0.1, 0.5, 0.9],
        distance_m=10.0,
    )

    adapter = TransportCalibrationAdapter(
        experiments=[exp], params_to_fit=["dispersivity"]
    )

    pest = PESTGLM.from_problem(adapter)

    # Eval
    residuals = pest._objective_function(pest.x0)

    # Check
    assert mock_run.called
    # Residual = Sim - Obs.
    # 0.15 - 0.1 = 0.05
    assert np.isclose(residuals[0], 0.05)
