"""Test for M3 advanced nuclear tracers and Gamma TTD."""

import pytest
import numpy as np
from hydrosheaf.nuclear.joint_lpm import fit_lpm_model, predict_lpm_tracers

def test_gamma_ttd_prediction():
    # Mean age 100 years, shape 1.0 (approaches Exponential)
    params = {"mean_age_years": 100.0, "shape": 1.0}
    sample_year = 2020.0
    
    # Predict Ar-39 (Half-life 269 years)
    # For a purely exponential model with mean 100, 
    # expected = C0 / (1 + lambda * tau) = 100 / (1 + (ln2/269) * 100)
    # lambda = 0.002576
    # expected = 100 / (1 + 0.2576) approx 79.5
    
    res = predict_lpm_tracers("GA", params, sample_year, ["39Ar"])
    assert "39Ar" in res
    assert 70.0 < res["39Ar"] < 85.0

def test_kr85_joint_fit():
    # Sample with modern Kr-85
    obs = {
        "kr85_pptv": 50.0, # Moderately modern (atmospheric equivalent)
        "tritium_TU": 5.0
    }
    sample_year = 2015.0
    
    fit = fit_lpm_model(obs, sample_year=sample_year, model="EM")
    assert fit.converged
    assert fit.n_tracers == 2
    assert "85Kr" in [r.tracer for r in fit.residuals]

if __name__ == "__main__":
    test_gamma_ttd_prediction()
    test_kr85_joint_fit()
    print("M3 Advanced Nuclear tests passed!")
