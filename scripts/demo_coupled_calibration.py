import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import shutil

# Ensure we can import hydrosheaf
sys.path.append(os.getcwd())

from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.adapters import (
    GenericFunctionAdapter, 
    CompositeCalibrationAdapter,
    KineticCalibrationAdapter,
    KineticExperiment
)
from hydrosheaf.calibration.pestpp.runner import run_pestpp
from hydrosheaf.reactive_transport import KineticParameters

# --- Mocking PHREEQC for Demo ---
from unittest.mock import MagicMock
import hydrosheaf.calibration.adapters as adapters

def mock_phreeqc_runner(initial_solution, kinetics_block, residence_time_days, config, **kwargs):
    import re
    # Look for: -parms <k> <A>
    match = re.search(r"-parms\s+([\d\.eE\+\-]+)\s+([\d\.eE\+\-]+)", kinetics_block)
    
    k = 1e-6 # Default
    A = 1.0
    if match:
        k = float(match.group(1))
        A = float(match.group(2))
        
    ca_0 = initial_solution.get("Ca", 0.0)
    rate = k * A * 1e5 # Arbitrary scaling to get concentration changes
    ca_final = ca_0 + rate * residence_time_days
    
    return {
        "success": True,
        "final_composition": [ca_final, 0, 0, 0, 0, 0, 0, 0, 0, 0], # Ca is first
        "time_series": {},
        "si_series": {}
    }

# Apply the mock
adapters.run_phreeqc_kinetic = mock_phreeqc_runner

# --- Transport Model (Physics) ---
from demo_models import coupled_model

def run_coupled_demo():
    print("=== Hydrosheaf Coupled Physics-Chemistry Calibration ===")
    print("Goal: Constrain groundwater velocity using both hydraulic and chemical data.")
    
    # Parameters
    # Velocity: True 0.5
    # Log K: True -5.0 (k=1e-5)
    params = [
        AdjustableParameter("velocity", 0.2, 0.01, 2.0),
        AdjustableParameter("log_k_calcite", -6.0, -8.0, -4.0)
    ]
    
    # Observations
    # Physics: 200 days
    # Chem: 1e-5 * 200 * 1e3 = 2.0 mg/L (example)
    obs = [
        Observation("obs_travel_time", 200.0, weight=1.0),
        Observation("obs_calcium", 2.0, weight=5.0) # High weight on chem
    ]
    
    prob_coupled = GenericFunctionAdapter(coupled_model, params, obs)
    
    # 3. Run PEST++ IES
    work_dir = Path("demo_coupled_work")
    if work_dir.exists():
        shutil.rmtree(work_dir)
        
    print(f"Running in {work_dir.absolute()}...")
    
    result = run_pestpp(
        problem=prob_coupled,
        engine="pestpp-ies",
        work_dir=str(work_dir),
        case_name="coupled",
        max_nfev=3,
        pestpp_options={
            "ies_num_reals": 15,
            "ies_no_noise": "true"
        }
    )
    
    if result["success"]:
        print("\n=== CALIBRATION SUCCESS ===")
        print(f"Calibrated Velocity: {result['optimal_parameters']['velocity']:.4f} m/d (True: 0.5)")
        print(f"Calibrated Log K:    {result['optimal_parameters']['log_k_calcite']:.4f} (True: -5.0)")
        
        unc = result["parameter_uncertainty"]
        print(f"\nUncertainty (Sigma):")
        print(f"Velocity Sigma: {unc.get('velocity', 0):.4f}")
        print(f"Log K Sigma:    {unc.get('log_k_calcite', 0):.4f}")
        
    else:
        print("\nFAILED")
        print(result.get("message"))

if __name__ == "__main__":
    run_coupled_demo()
