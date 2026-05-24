"""
Advanced verification tests for PEST++ runner, control file writing, and output parsing.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.pestpp.runner import (
    PestControlFile,
    PestParameter,
    PestObservation,
    PestRunner,
    parse_pest_residuals,
    parse_pest_sensitivities,
    parse_pest_identifiability,
    parse_ies_csv,
    generate_sweep_csv
)

class GlobalMockProblem:
    def __init__(self, params=None, obs=None):
        self.params = params or []
        self.obs = obs or []
    def get_parameters(self):
        return self.params
    def get_observations(self):
        return self.obs

def test_pest_control_file_prior_info(tmp_path):
    """Verify prior-information equations generation from parameter prior definitions."""
    pst = PestControlFile()
    pst.add_parameter(PestParameter("p1", "none", "factor", 10.0, 1.0, 100.0, prior_mean=12.0, prior_sigma=2.0))
    pst.add_parameter(PestParameter("p2", "log", "factor", 5.0, 0.1, 50.0, prior_mean=6.0, prior_sigma=0.3))
    
    pst_file = tmp_path / "test_prior.pst"
    pst.write(pst_file)
    
    content = pst_file.read_text()
    assert "* prior information" in content
    assert "pi_p1 1.0 * p1 = 12.0 0.5 prior" in content
    assert "pi_p2 1.0 * log(p2) =" in content
    # log10(6) approx 0.77815
    assert "0.77815" in content
    # weight = 1.0 / 0.3 approx 3.3333
    assert "3.333" in content
    assert "prior" in pst.observation_groups
    assert "prior" in content.split("* observation groups")[1].split("*")[0]

def test_pest_control_file_fixed_and_tied(tmp_path):
    """Verify parameter fixed and tied transformations are formatted correctly."""
    p_fixed = AdjustableParameter("p_fixed", 1.5, fixed=True)
    p_tied = AdjustableParameter("p_tied", 3.0, tied_to="p_parent")
    
    p1 = PestParameter.from_adjustable_parameter(p_fixed)
    p2 = PestParameter.from_adjustable_parameter(p_tied)
    
    assert p1.trans == "fixed"
    assert p2.trans == "tied"
    assert p2.tied_to == "p_parent"
    
    pst = PestControlFile()
    pst.add_parameter(p1)
    pst.add_parameter(p2)
    
    pst_file = tmp_path / "test_fixed.pst"
    pst.write(pst_file)
    content = pst_file.read_text()
    
    assert "p_fixed fixed factor 1.5" in content
    assert "p_tied tied factor 3.0" in content
    assert "p_tied tied factor 3.0 -1e+30 1e+30 pargp 1.0 0.0 p_parent" in content

def test_pest_control_file_multiple_files(tmp_path):
    """Verify writing multiple template and instruction files."""
    pst = PestControlFile()
    pst.tpl_files = ["t1.tpl", "t2.tpl"]
    pst.ins_files = ["i1.ins", "i2.ins"]
    pst.model_input_files = ["in1.txt", "in2.txt"]
    pst.model_output_files = ["out1.txt", "out2.txt"]
    
    pst_file = tmp_path / "test_multi.pst"
    pst.write(pst_file)
    
    content = pst_file.read_text()
    assert "2 2 single" in content  # ntplfle=2, ninsfle=2
    assert "0 0 1 0 1" in content
    assert "t1.tpl in1.txt" in content
    assert "t2.tpl in2.txt" in content
    assert "i1.ins out1.txt" in content
    assert "i2.ins out2.txt" in content

def test_glm_parser(tmp_path):
    """Verify parsing GLM residuals, sensitivities, covariance path, and identifiability."""
    # Write mock residual (.rei) file
    rei_content = """
                    REI File
 Name  Group  Measured  Simulated  Residual  Weight
 o1    obsgp  10.0      9.5        -0.5      2.0
 o2    obsgp  20.0      20.2       0.2       1.0
 o3    prior  1.0       1.1        0.1       10.0
"""
    rei_path = tmp_path / "case.rei"
    rei_path.write_text(rei_content)
    
    # Write mock sensitivity (.sen) file
    sen_content = """
PARAMETER NAME      COMPOSITE SENSITIVITY
p1                  15.5
p2                  0.05
"""
    sen_path = tmp_path / "case.sen"
    sen_path.write_text(sen_content)
    
    # Write mock identifiability (.id) file
    id_content = """
p1                  0.95
p2                  0.12
"""
    id_path = tmp_path / "case.id"
    id_path.write_text(id_content)
    
    # Check Residual Parser
    residuals, simulated, phi_by_group = parse_pest_residuals(rei_path)
    assert len(residuals) == 3
    assert residuals["o1"] == -0.5
    assert simulated["o2"] == 20.2
    assert np.isclose(phi_by_group["obsgp"], (-0.5 * 2.0)**2 + (0.2 * 1.0)**2)
    assert np.isclose(phi_by_group["prior"], (0.1 * 10.0)**2)
    
    # Check Sensitivity Parser
    sens = parse_pest_sensitivities(sen_path)
    assert sens["p1"] == 15.5
    assert sens["p2"] == 0.05
    
    # Check Identifiability Parser
    ident = parse_pest_identifiability(id_path)
    assert ident["p1"] == 0.95
    assert ident["p2"] == 0.12

def test_ies_parser(tmp_path):
    """Verify parsing IES ensemble CSV files and phi history."""
    # Create mock prior parameter ensemble CSV
    prior_df = pd.DataFrame({
        "real_name": ["r0", "r1"],
        "p1": [1.0, 1.2],
        "p2": [20.0, 18.0]
    })
    prior_path = tmp_path / "case.0.par.csv"
    prior_df.to_csv(prior_path, index=False)
    
    # Create mock posterior parameter ensemble CSV
    post_df = pd.DataFrame({
        "realname": ["r0", "r1"],
        "p1": [1.5, 1.6],
        "p2": [10.0, 11.0]
    })
    post_path = tmp_path / "case.2.par.csv"
    post_df.to_csv(post_path, index=False)
    
    # Create mock phi.actual.csv
    phi_df = pd.DataFrame({
        "iteration": [0, 1, 2],
        "r0": [100.0, 50.0, 5.0],
        "r1": [120.0, 60.0, 7.0]
    })
    phi_path = tmp_path / "case.phi.actual.csv"
    phi_df.to_csv(phi_path, index=False)
    
    # Verify Parser
    prior_par = parse_ies_csv(prior_path)
    assert prior_par is not None
    assert prior_par["p1"] == [1.0, 1.2]
    
    # Test runner output mapping mock
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0), AdjustableParameter("p2", 20.0)],
        obs=[Observation("o1", 5.0)]
    )

    runner = PestRunner(
        problem=prob,
        engine="pestpp-ies",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=2
    )
    
    # Call parser logic by mocking run dependencies
    # We create mock files in the work_dir first
    (tmp_path / "case.pst").write_text("dummy pst")
    
    # Mock subprocess run to avoid executing binary
    import subprocess
    orig_run = subprocess.run
    try:
        subprocess.run = lambda *args, **kwargs: None
        result = runner.run()
    finally:
        subprocess.run = orig_run
        
    assert result["success"]
    assert np.isclose(result["phi"], 6.0)  # (5.0 + 7.0) / 2
    assert result["optimal_parameters"]["p1"] == 1.55  # (1.5 + 1.6) / 2
    assert result["optimal_parameters"]["p2"] == 10.5
    assert result["phi_history"] == [110.0, 55.0, 6.0]

def test_sen_parser(tmp_path):
    """Verify Sobol and Morris parser works correctly."""
    # Mock Sobol file
    sobol_content = """parameter,first_order,total_order
p1,0.3,0.85
p2,0.1,0.05
"""
    (tmp_path / "case.sobol.csv").write_text(sobol_content)
    
    # Mock Morris file
    morris_content = """name,mean,std
p1,4.5,1.2
p2,0.2,0.01
"""
    (tmp_path / "case.morris.csv").write_text(morris_content)
    
    prob = GlobalMockProblem(params=[], obs=[])
        
    runner = PestRunner(
        problem=prob,
        engine="pestpp-sen",
        work_dir=str(tmp_path),
        case_name="case"
    )
    
    import subprocess
    orig_run = subprocess.run
    try:
        subprocess.run = lambda *args, **kwargs: None
        result = runner.run()
    finally:
        subprocess.run = orig_run
        
    assert result["success"]
    assert result["first_order_sobol"]["p1"] == 0.3
    assert result["total_sobol"]["p1"] == 0.85
    assert result["morris_elementary_effects"]["p1"] == [4.5, 1.2]
    assert result["ranked_importance"] == ["p1", "p2"]
    assert result["recommended_calibratable_subset"] == ["p1"]

def test_sweep_generator(tmp_path):
    """Verify that generate_sweep_csv produces the expected shape and headers."""
    params = [
        AdjustableParameter("p1", 5.0, lower_bound=1.0, upper_bound=10.0),
        AdjustableParameter("p2", 0.1, lower_bound=0.01, upper_bound=1.0, log_transform=True)
    ]
    
    filepath = tmp_path / "sweep.csv"
    generate_sweep_csv(filepath, params, n_runs=10, method="latin_hypercube")
    
    assert filepath.exists()
    df = pd.read_csv(filepath)
    assert len(df) == 10
    assert list(df.columns) == ["run_id", "p1", "p2"]
    assert (df["p1"] >= 1.0).all() and (df["p1"] <= 10.0).all()
    assert (df["p2"] >= 0.01).all() and (df["p2"] <= 1.0).all()
