"""
Tests for the Hydrosheaf Nuclear module.
"""
import pytest
import numpy as np
from hydrosheaf.nuclear.nuclides import TRITIUM, get_nuclide
from hydrosheaf.nuclear.input_history import build_default_tritium_input
from hydrosheaf.nuclear.lpm import convolve_input, piston_flow_model, exponential_model
from hydrosheaf.nuclear.invert import infer_age_from_tracer

def test_nuclide_properties():
    t3 = get_nuclide("3H")
    assert t3 is not None
    assert t3.symbol == "3H"
    assert abs(t3.half_life_years - 12.32) < 0.01

def test_lpm_convolution_pfm():
    # Test Piston Flow Model convolution
    # If age is 0, output should match input at sample time
    hist = build_default_tritium_input()
    sample_year = 2020.0
    val_in = hist.interpolate([sample_year])[0][0]
    
    # Predict with tau=0
    # lambda = ln(2)/12.32
    lambda_y = np.log(2) / 12.32
    
    val_out = convolve_input(sample_year, 0.0, hist.years, hist.values, lambda_y, model_type="PFM")
    assert abs(val_out - val_in) < 1.0 # Should be very close

    # Predict with tau=12.32 years (1 half life)
    # Output should be Input(2020 - 12.32) * 0.5
    t_recharge = sample_year - 12.32
    val_recharge = hist.interpolate([t_recharge])[0][0]
    expected = val_recharge * 0.5
    
    val_out_aged = convolve_input(sample_year, 12.32, hist.years, hist.values, lambda_y, model_type="PFM")
    assert abs(val_out_aged - expected) < 1.0

def test_inverse_inference_synthetic():
    # Create synthetic observation
    true_age = 25.0 # years
    sample_year = 2023.0
    hist = build_default_tritium_input()
    lambda_y = np.log(2) / TRITIUM.half_life_years
    
    # Forward model PFM
    obs_val = convolve_input(sample_year, true_age, hist.years, hist.values, lambda_y, model_type="PFM")
    obs_sigma = 1.0 # 1 TU uncertainty
    
    # Invert
    result = infer_age_from_tracer(
        obs_val, 
        obs_sigma, 
        sample_year, 
        TRITIUM, 
        input_hist=hist, 
        model="PFM",
        tau_range=(0, 100)
    )
    
    print(f"True Age: {true_age}, Inferred MAP: {result['tau_map_years']}")
    
    # Check if true age is within CI or close to MAP
    # Note: PFM with Tritium can be multi-modal.
    # But for recent years (post-bomb peak washout), monotonic decay is often observed if input is steady.
    # The default input history has a bomb peak. 2023 - 25 = 1998. Post-peak.
    # However, 2023 - 60 = 1963 (peak). So very old water might look like young water?
    # Actually, 1963 water would have decayed by 2^-(60/12.32) approx 2^-5 = 1/32.
    # 3000 TU / 32 approx 100 TU.
    # Modern rain is ~5-10 TU.
    # So 1963 water is still "hotter" than modern water.
    # So there shouldn't be too much ambiguity between 25 year old water and 60 year old water.
    
    assert abs(result["tau_map_years"] - true_age) < 5.0 # Allow some grid error
    assert result["tau_ci_low_years"] <= true_age + 2.0
    assert result["tau_ci_high_years"] >= true_age - 2.0

if __name__ == "__main__":
    test_nuclide_properties()
    test_lpm_convolution_pfm()
    test_inverse_inference_synthetic()
    print("All nuclear tests passed!")
