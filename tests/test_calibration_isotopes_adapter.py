"""
Tests for Water Isotope Mixing Adapter (PEST calibration).
"""

import pytest
import numpy as np
from hydrosheaf.calibration.adapters_iso import WaterIsotopeMixingAdapter, WaterEndmember
from hydrosheaf.calibration.definitions import AdjustableParameter

def test_softmax_fractions():
    """Test that mixing fractions sum to 1 and respect logits."""
    
    endmembers = [
        WaterEndmember("Rain", -5.0, -30.0),
        WaterEndmember("River", -7.0, -45.0),
        WaterEndmember("Sea", 0.0, 0.0)
    ]
    
    observations = {
        "S1": {"18O": -2.5, "2H": -15.0} # roughly 50% Rain, 50% Sea
    }
    
    group_map = {"S1": "GroupA"}
    
    adapter = WaterIsotopeMixingAdapter(
        endmembers=endmembers,
        observations=observations,
        group_map=group_map
    )
    
    # 3 Endmembers -> 2 Thetas (last is 0)
    # If theta1=0, theta2=0 -> all 3 logits are 0 -> fractions 1/3, 1/3, 1/3
    params = {
        "theta_GroupA_Rain": 0.0,
        "theta_GroupA_River": 0.0
    }
    
    res = adapter.run_model(params)
    
    # Expected d18O = (-5 + -7 + 0)/3 = -4.0
    assert np.isclose(res["S1_18O"], -4.0)
    
    # Test large theta for Sea (implicit third param)
    # If theta_Rain = -10, theta_River = -10, theta_Sea (implicit) = 0
    # exp(-10) is tiny. exp(0) = 1.
    # Fraction of Sea -> 1.0
    
    params = {
        "theta_GroupA_Rain": -10.0,
        "theta_GroupA_River": -10.0
    }
    res = adapter.run_model(params)
    assert np.isclose(res["S1_18O"], 0.0, atol=1e-3)


def test_parameter_generation():
    """Test that correct parameters are generated for multiple groups."""
    endmembers = [
        WaterEndmember("A", 0, 0),
        WaterEndmember("B", 0, 0)
    ]
    
    # 2 groups
    group_map = {"S1": "G1", "S2": "G2"}
    obs = {"S1": {}, "S2": {}}
    
    adapter = WaterIsotopeMixingAdapter(endmembers, obs, group_map)
    params = adapter.get_parameters()
    
    names = {p.name for p in params}
    # For G1: theta_G1_A (B is implicit)
    # For G2: theta_G2_A (B is implicit)
    assert "theta_G1_A" in names
    assert "theta_G2_A" in names
    assert len(names) == 2


def test_prediction_correctness():
    """Test exact mixing calculation."""
    # Endmembers: A(0,0), B(10,100)
    endmembers = [
        WaterEndmember("A", 0.0, 0.0),
        WaterEndmember("B", 10.0, 100.0)
    ]
    
    obs = {"S1": {"18O": 5.0, "2H": 50.0}} # Target: 50/50 mix
    group_map = {"S1": "G1"}
    
    adapter = WaterIsotopeMixingAdapter(endmembers, obs, group_map)
    
    # To get 50/50, logits must be equal.
    # theta_A = x, theta_B(implicit) = 0.
    # So we set theta_A = 0.
    
    params = {"theta_G1_A": 0.0}
    res = adapter.run_model(params)
    
    assert np.isclose(res["S1_18O"], 5.0)
    assert np.isclose(res["S1_2H"], 50.0)
    
    # To get 73% A, 27% B (logistic sigma(1) approx 0.73)
    # theta_A = 1.0, theta_B = 0.0
    # exp(1)/(exp(1)+1) = 2.718/3.718 = 0.731
    # Pred = 0.731*0 + 0.269*10 = 2.69
    
    params = {"theta_G1_A": 1.0}
    res = adapter.run_model(params)
    
    f_A = np.exp(1.0) / (np.exp(1.0) + 1.0)
    f_B = 1.0 - f_A
    expected_18O = f_A*0 + f_B*10
    
    assert np.isclose(res["S1_18O"], expected_18O)
