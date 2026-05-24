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
    parse_sen_msn,
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
    """Verify Sobol, Morris, and real PESTPP-SEN .msn parsing works correctly."""
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

    msn_content = """parameter_name,n_samples,sen_mean,sen_mean_abs,sen_std_dev
p3,4,-2.0,2.0,0.4
"""
    (tmp_path / "case.msn").write_text(msn_content)

    msn = parse_sen_msn(tmp_path / "case.msn")
    assert msn["p3"]["sen_mean_abs"] == 2.0
    
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
    assert result["morris_elementary_effects"]["p3"] == [2.0, 0.4]
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


def test_parallel_spawn_flags(tmp_path):
    """Verify that n_workers > 1 triggers manager/agent /h flags."""
    from unittest.mock import patch, MagicMock

    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0), AdjustableParameter("p2", 20.0)],
        obs=[Observation("o1", 5.0)]
    )

    runner = PestRunner(
        problem=prob,
        engine="pestpp-ies",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=1,
        n_workers=3,
    )

    # Create dummy PST so prepare() succeeds
    (tmp_path / "case.pst").write_text("dummy")

    # Capture both Popen calls (manager + agents)
    popen_calls = []

    def mock_popen(cmd, **kwargs):
        mock_proc = MagicMock()
        mock_proc.pid = 12345
        mock_proc.returncode = 0
        mock_proc.poll.return_value = 0
        popen_calls.append(cmd)
        return mock_proc

    with patch("subprocess.Popen", side_effect=mock_popen):
        result = runner.run()

    # Manager should be first call with /h :<port>
    assert popen_calls, "No Popen calls were made"
    manager_cmd = popen_calls[0]
    assert manager_cmd[0].endswith("pestpp-ies") or "pestpp-ies" in str(manager_cmd[0])
    assert "/h" in manager_cmd
    # Port should be somewhere after /h
    port_arg = manager_cmd[manager_cmd.index("/h") + 1]
    assert port_arg.startswith(":")

    # Should have 1 manager + 3 agents
    assert len(popen_calls) == 4, f"Expected 4 Popen calls, got {len(popen_calls)}"
    for agent_cmd in popen_calls[1:]:
        assert "/h" in agent_cmd
        host_port = agent_cmd[agent_cmd.index("/h") + 1]
        assert host_port.startswith("localhost:")


def test_ies_covariance_and_localizer_generated(tmp_path):
    """Verify auto-generated parcov, obscov, and localizer files for IES."""
    from hydrosheaf.calibration.pestpp.runner import (
        write_parcov_matrix, write_obscov_matrix, write_localizer_matrix
    )

    params = [
        AdjustableParameter("p1", 5.0, lower_bound=1.0, upper_bound=10.0),
        AdjustableParameter("p2", 0.1, lower_bound=0.01, upper_bound=1.0, log_transform=True)
    ]
    obs = [
        Observation("o1", 10.0, weight=2.0),
        Observation("o2", 20.0, weight=1.0)
    ]

    parcov_file = tmp_path / "parcov.csv"
    write_parcov_matrix(parcov_file, params)
    content = parcov_file.read_text()
    assert "p1" in content
    assert "p2" in content
    # Check diagonal non-zero
    lines = content.strip().splitlines()
    header = lines[0].split(",")
    assert header[1] == "p1"
    assert header[2] == "p2"

    obscov_file = tmp_path / "obscov.csv"
    write_obscov_matrix(obscov_file, obs)
    content = obscov_file.read_text()
    assert "o1" in content
    assert "o2" in content

    loc_file = tmp_path / "localizer.csv"
    write_localizer_matrix(loc_file, params, obs)
    content = loc_file.read_text()
    lines = content.strip().splitlines()
    assert lines[0].split(",")[1] == "p1"
    assert lines[1] == "o1,1.0,1.0"
    assert lines[2] == "o2,1.0,1.0"


def test_ies_prepare_covariance_pass_passthrough(tmp_path):
    """Verify that parcov/obscov/localizer paths are passed into pestpp_options when specified in options."""
    import subprocess

    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0, lower_bound=0.1, upper_bound=10.0),
                AdjustableParameter("p2", 2.0, lower_bound=0.1, upper_bound=10.0)],
        obs=[Observation("o1", 5.0, weight=1.0)]
    )

    # Provide explicit covariance paths
    runner = PestRunner(
        problem=prob,
        engine="pestpp-ies",
        work_dir=str(tmp_path),
        case_name="case",
        pestpp_options={
            "parcov": "user_parcov.csv",
            "obscov": "user_obscov.csv",
            "ies_localizer": "user_loc.csv",
            "ies_restart_parameter_ensemble": "restart_par.csv",
            "ies_restart_observation_ensemble": "restart_obs.csv",
        }
    )
    pst = runner.prepare()

    # Verify paths appear in written .pst file
    pst_path = tmp_path / "case.pst"
    assert pst_path.exists()
    content = pst_path.read_text()
    assert "++parcov(user_parcov.csv)" in content
    assert "++obscov(user_obscov.csv)" in content
    assert "++ies_localizer(user_loc.csv)" in content
    assert "++ies_restart_parameter_ensemble(restart_par.csv)" in content
    assert "++ies_restart_observation_ensemble(restart_obs.csv)" in content

    # Run mocked
    orig_run = subprocess.run
    try:
        subprocess.run = lambda *args, **kwargs: None
        result = runner.run()
    finally:
        subprocess.run = orig_run

    assert result["success"]
    assert result.get("parcov_path") == "user_parcov.csv"
    assert result.get("obscov_path") == "user_obscov.csv"
    assert result.get("localizer_path") == "user_loc.csv"
    assert result.get("restart_from") == "restart_par.csv"


def test_ies_forecast_summary_parsing(tmp_path):
    """Verify that forecast summaries are computed when forecasts option is set."""
    import subprocess

    prior_df = pd.DataFrame({
        "real_name": ["r0", "r1"],
        "p1": [1.0, 1.2],
    })
    prior_df.to_csv(tmp_path / "case.0.par.csv", index=False)

    # Posterior with forecast observations
    post_df = pd.DataFrame({
        "realname": ["r0", "r1"],
        "p1": [1.5, 1.6],
        "fc1": [10.0, 12.0],
        "fc2": [5.0, 6.0],
    })
    post_df.to_csv(tmp_path / "case.2.obs.csv", index=False)

    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0)],
        obs=[Observation("fc1", 10.0), Observation("fc2", 5.0)]
    )

    runner = PestRunner(
        problem=prob,
        engine="pestpp-ies",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=2,
        pestpp_options={"forecasts": "fc1,fc2"}
    )

    orig_run = subprocess.run
    try:
        subprocess.run = lambda *args, **kwargs: None
        result = runner.run()
    finally:
        subprocess.run = orig_run

    assert result["success"]
    assert "posterior_forecast_summaries" in result
    assert "fc1" in result["posterior_forecast_summaries"]
    assert "fc2" in result["posterior_forecast_summaries"]
    fc1 = result["posterior_forecast_summaries"]["fc1"]
    assert np.isclose(fc1["mean"], 11.0)
    assert np.isclose(fc1["std"], 1.0)


def test_sen_workflow_defaults(tmp_path):
    """Verify SEN prepare() does not emit rejected legacy option defaults."""
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0), AdjustableParameter("p2", 2.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-sen",
        work_dir=str(tmp_path),
        case_name="case",
    )
    pst = runner.prepare()
    pst_path = tmp_path / "case.pst"
    assert pst_path.exists()
    content = pst_path.read_text()
    assert "++sen_method" not in content
    assert "++sen_num_samples" not in content
    assert "++sen_morris_delta" not in content


def test_sen_legacy_options_are_normalized(tmp_path):
    """Verify old Hydrosheaf SEN option names are translated to PEST++ GSA names."""
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0), AdjustableParameter("p2", 2.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-sen",
        work_dir=str(tmp_path),
        case_name="case",
        pestpp_options={
            "sen_method": "morris",
            "sen_num_samples": 12,
            "sen_morris_delta": 0.2,
        },
    )
    runner.prepare()
    content = (tmp_path / "case.pst").read_text()
    assert "++gsa_method(morris)" in content
    assert "++gsa_morris_r(12)" in content
    assert "++gsa_morris_delta(0.2)" in content
    assert "++sen_method" not in content


def test_swp_workflow_defaults(tmp_path):
    """Verify SWP prepare() generates sweep CSV and sets output CSV."""
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 5.0, lower_bound=1.0, upper_bound=10.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-swp",
        work_dir=str(tmp_path),
        case_name="case",
        pestpp_options={"sweep_n_runs": 5}
    )
    pst = runner.prepare()
    sweep_file = tmp_path / "case.swp.in.csv"
    assert sweep_file.exists()
    df = pd.read_csv(sweep_file)
    assert len(df) == 5

    pst_path = tmp_path / "case.pst"
    content = pst_path.read_text()
    assert "++sweep_parameter_csv_file" in content
    assert "++sweep_output_csv_file" in content
    assert "++sweep_n_runs" not in content
    assert "++sweep_sampler" not in content


def test_mou_workflow_defaults(tmp_path):
    """Verify MOU prepare() sets supported population and generator defaults."""
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-mou",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=20,
    )
    pst = runner.prepare()
    pst_path = tmp_path / "case.pst"
    content = pst_path.read_text()
    assert "++mou_population_size" in content
    assert "++mou_generator" in content
    assert "++mou_max_generations" not in content
    assert "++opt_risk" in content


def test_mou_legacy_options_are_normalized(tmp_path):
    """Verify old MOU/OPT aliases are normalized to supported PEST++ option names."""
    prob = GlobalMockProblem(
        params=[
            AdjustableParameter("p1", 1.0, group="dv"),
            AdjustableParameter("p2", 2.0, group="uncertain"),
        ],
        obs=[Observation("obj_cost", 5.0, group="less_than_obj")]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-mou",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=20,
        pestpp_options={
            "mou_max_generations": 3,
            "dec_var_groups": "dec_var p1",
            "risk": 0.8,
        },
    )
    pst = runner.prepare()
    content = (tmp_path / "case.pst").read_text()
    assert pst.control_data.noptmax == 3
    assert "++opt_dec_var_groups(dv)" in content
    assert "++opt_risk(0.8)" in content
    assert "++mou_max_generations" not in content
    assert "++dec_var_groups" not in content


def test_da_workflow_defaults(tmp_path):
    """Verify DA prepare() avoids unsupported cycle-count aliases."""
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-da",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=3,
    )
    pst = runner.prepare()
    pst_path = tmp_path / "case.pst"
    content = pst_path.read_text()
    assert "++da_num_cycles" not in content


def test_da_internal_state_options_are_not_written(tmp_path):
    """Verify Hydrosheaf-only DA state helpers generate files without invalid ++ options."""
    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0, lower_bound=0.1, upper_bound=2.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-da",
        work_dir=str(tmp_path),
        case_name="case",
        pestpp_options={
            "ies_num_reals": 3,
            "da_state_names": ["head_1", "head_2"],
            "da_state_initial_values": [10.0, 11.0],
        },
    )
    runner.prepare()
    content = (tmp_path / "case.pst").read_text()
    assert (tmp_path / "case.state.0.par.csv").exists()
    assert "++da_state_parameter_ensemble(case.state.0.par.csv)" in content
    assert "++da_state_names" not in content
    assert "++da_state_initial_values" not in content


def test_da_real_output_pattern_parser(tmp_path):
    """Verify DA parser recognizes real PESTPP-DA cycle/global output names."""
    import subprocess

    pd.DataFrame({
        "real_name": ["r0", "r1"],
        "p1": [1.0, 2.0],
    }).to_csv(tmp_path / "case.0.0.par.csv", index=False)
    pd.DataFrame({
        "real_name": ["r0", "r1"],
        "o1": [4.0, 5.0],
    }).to_csv(tmp_path / "case.0.0.obs.csv", index=False)
    pd.DataFrame({
        "iteration": [0, 1],
        "r0": [10.0, 4.0],
        "r1": [12.0, 6.0],
    }).to_csv(tmp_path / "case.global.phi.actual.csv", index=False)

    prob = GlobalMockProblem(
        params=[AdjustableParameter("p1", 1.0)],
        obs=[Observation("o1", 5.0)]
    )
    runner = PestRunner(
        problem=prob,
        engine="pestpp-da",
        work_dir=str(tmp_path),
        case_name="case",
        max_nfev=1,
    )

    orig_run = subprocess.run
    try:
        subprocess.run = lambda *args, **kwargs: None
        result = runner.run()
    finally:
        subprocess.run = orig_run

    assert result["success"]
    assert result["phi_history"] == [11.0, 5.0]
    assert result["cycles"][0]["cycle"] == 0
    assert result["cycles"][0]["iteration"] == 0
    assert result["posterior_parameters"]["p1"] == [1.0, 2.0]
