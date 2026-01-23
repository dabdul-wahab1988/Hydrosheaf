import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import shutil

# Ensure we can import hydrosheaf
sys.path.append(os.getcwd())

from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.adapters import GenericFunctionAdapter, CompositeCalibrationAdapter
from hydrosheaf.calibration.pestpp.runner import run_pestpp, get_executable_path

def vadose_model(params):
    """
    Simulates recharge = precipitation * factor.
    """
    # params: 'recharge_factor'
    rf = params.get("recharge_factor", 0.5)
    precip = 100.0 # mm/yr
    recharge = precip * rf
    return {"obs_recharge": recharge}

def transport_model(params):
    """
    Simulates concentration = source / (1 + decay * time).
    Depends on 'recharge_factor' implicitly (conceptual link)? 
    Let's say C_source depends on recharge (dilution).
    C_source = Mass / Recharge.
    """
    rf = params.get("recharge_factor", 0.5)
    decay = params.get("decay", 0.1)
    
    # Physics link
    recharge_vol = rf * 100.0
    mass_loading = 1000.0
    c_source = mass_loading / (recharge_vol + 1e-6)
    
    # Observations at T=10
    t = 10.0
    c_out = c_source * np.exp(-decay * t)
    
    return {"obs_conc_t10": c_out}

def run_demo():
    print("=== Hydrosheaf Advanced Calibration Demo: PEST++ IES ===")
    
    # 1. Setup Parameters
    # Shared parameter: recharge_factor
    p_rf = AdjustableParameter("recharge_factor", 0.5, 0.1, 1.0, prior_mean=0.5, prior_sigma=0.2)
    # Transport only: decay
    p_decay = AdjustableParameter("decay", 0.1, 0.0, 1.0, prior_mean=0.1, prior_sigma=0.05)
    
    # 2. Setup Observations (Synthetic Truth)
    # True RF = 0.6 -> Recharge = 60. C_source = 1000/60 = 16.66
    # True Decay = 0.15 -> C_out = 16.66 * exp(-1.5) = 16.66 * 0.223 = 3.71
    obs_r = Observation("obs_recharge", 60.0, weight=1.0)
    obs_c = Observation("obs_conc_t10", 3.71, weight=10.0) # Higher weight
    
    # 3. Create Adapters
    # Note: GenericAdapter usually takes a fixed list of params.
    # We pass the relevant ones.
    
    # Vadose problem sees 'recharge_factor'
    prob1 = GenericFunctionAdapter(vadose_model, [p_rf], [obs_r])
    
    # Transport problem sees 'recharge_factor' and 'decay'
    prob2 = GenericFunctionAdapter(transport_model, [p_rf, p_decay], [obs_c])
    
    # 4. Composite
    master_prob = CompositeCalibrationAdapter([prob1, prob2])
    
    # 5. Run PEST++ IES
    work_dir = Path("demo_ies_work")
    if work_dir.exists():
        shutil.rmtree(work_dir)
    
    print(f"Running in {work_dir.absolute()}...")
    
    # PEST++ IES options would normally go into the PST control data
    # Standard 'pestpp-ies' runs with defaults (n_ensemble etc)
    # We can rely on defaults for this demo.
    
    result = run_pestpp(
        problem=master_prob,
        engine="pestpp-ies",
        work_dir=str(work_dir),
        case_name="demo",
        max_nfev=3, # Low iterations for demo speed
        pestpp_options={
            "ies_num_reals": 10, # Small ensemble for speed
            "ies_no_noise": "true" # Simplified for demo
        }
    )
    
    if result["success"]:
        print("\nSUCCESS!")
        print("Optimal Parameters (Mean of Posterior):")
        print(result["optimal_parameters"])
        print("\nUncertainty (Std Dev):")
        print(result["parameter_uncertainty"])
        
        # Check ensemble results
        ens_info = result.get("ensemble_results")
        if ens_info:
            print(f"\nEnsemble file: {ens_info['parameters_file']}")
            df = pd.read_csv(ens_info['parameters_file'])
            print("\nFirst 5 Realizations:")
            print(df.head())
            
    else:
        print("\nFAILED.")
        print(result.get("message"))
        print("STDOUT:", result.get("stdout"))
        print("STDERR:", result.get("stderr"))

if __name__ == "__main__":
    run_demo()
