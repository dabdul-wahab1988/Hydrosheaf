import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import shutil
import pytest
from unittest.mock import MagicMock, patch

# Ensure we can import hydrosheaf from local source
sys.path.insert(0, os.getcwd())

from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.adapters import (
    KineticCalibrationAdapter,
    KineticExperiment
)
from hydrosheaf.reactive_transport import KineticParameters
from hydrosheaf.data.units import mgL_to_mmolL

# --- Test Data ---
def test_unit_conversion():
    """Test that mg/L observations are correctly converted to mmol/L."""
    
    # Setup Experiment in mg/L
    # Ca: 40.078 mg/mmol. 40.078 mg/L = 1.0 mmol/L
    exp = KineticExperiment(
        id="test_units",
        initial_solution={"Ca": 40.078},
        residence_time_days=1.0,
        reaction_extents=[0.0],
        reaction_labels=["calcite"],
        observations={"Ca": 80.156}, # Target: 2.0 mmol/L
        units="mg/L"
    )
    
    # Setup Adapter
    config = MagicMock()
    config.ion_order = ["Ca", "C"]
    
    adapter = KineticCalibrationAdapter(
        base_params={"calcite": KineticParameters("calcite", 1e-6, 1.0)},
        experiments=[exp],
        config=config
    )
    
    # Check conversion
    processed_exp = adapter.experiments[0]
    
    # Initial solution
    assert processed_exp.units == "mmol/L"
    assert np.isclose(processed_exp.initial_solution["Ca"], 1.0)
    
    # Observations
    pest_obs = adapter.get_observations()
    assert len(pest_obs) == 1
    assert pest_obs[0].name == "test_units_Ca"
    assert np.isclose(pest_obs[0].value, 2.0)

def test_parameter_scoping():
    """Test global vs edge-specific parameter resolution."""
    
    exp1 = KineticExperiment(
        id="exp_global",
        edge_id="edge_A",
        initial_solution={},
        residence_time_days=1.0,
        reaction_extents=[],
        reaction_labels=[],
        observations={},
        units="mmol/L"
    )
    
    exp2 = KineticExperiment(
        id="exp_local",
        edge_id="edge_B",
        initial_solution={},
        residence_time_days=1.0,
        reaction_extents=[],
        reaction_labels=[],
        observations={},
        units="mmol/L"
    )
    
    base_params = {"calcite": KineticParameters("calcite", 1e-6, 1.0)}
    
    adapter = KineticCalibrationAdapter(
        base_params=base_params,
        experiments=[exp1, exp2],
        config=MagicMock()
    )
    
    # Case 1: Global Only
    # Note: param_values dict contains REAL values (already transformed back from log space if needed)
    # The parameter name is "log_k_calcite" but the value is k (e.g. 1e-5).
    param_values = {"log_k_calcite": 1e-5} 
    
    # Test internal resolution
    p1 = adapter._resolve_parameters(exp1, param_values)
    assert np.isclose(p1["calcite"].rate_constant, 1e-5), f"p1 mismatch: {p1['calcite'].rate_constant}"
    
    p2 = adapter._resolve_parameters(exp2, param_values)
    assert np.isclose(p2["calcite"].rate_constant, 1e-5), f"p2 mismatch: {p2['calcite'].rate_constant}"
    
    # Case 2: Edge Specific Override
    param_values["log_k_calcite__edge_B"] = 1e-4
    
    p1 = adapter._resolve_parameters(exp1, param_values)
    assert np.isclose(p1["calcite"].rate_constant, 1e-5), f"p1 override mismatch: {p1['calcite'].rate_constant}" # Unchanged
    
    p2 = adapter._resolve_parameters(exp2, param_values)
    assert np.isclose(p2["calcite"].rate_constant, 1e-4), f"p2 override mismatch: {p2['calcite'].rate_constant}" # Overridden

if __name__ == "__main__":
    # Run tests manually
    try:
        test_unit_conversion()
        print("Unit conversion test passed.")
        test_parameter_scoping()
        print("Parameter scoping test passed.")
    except AssertionError as e:
        print(f"Test failed: {e}")
    except Exception as e:
        print(f"Error: {e}")
