import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, patch

from hydrosheaf.graph.types import Edge
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.adapters import (
    TopologyCalibrationObservation,
    TopologyOuterLoopCalibrator,
    VadoseCalibrationAdapter,
    NitrateSourceCalibrationAdapter,
    NitrateSourceCalibrationObservation,
    AgeTemporalCalibrationAdapter,
    AgeTemporalExperiment,
)
from hydrosheaf.calibration.glm import PESTGLM


def test_outer_loop_topology_search():
    # 1. Setup candidate edges
    edge_a_b = Edge("A_B", "A", "B", attrs={"edge_confidence": 0.9})
    edge_b_c = Edge("B_C", "B", "C", attrs={"edge_confidence": 0.8})
    edge_c_a = Edge("C_A", "C", "A", attrs={"edge_confidence": 0.7})  # Creates a cycle
    
    # Observations
    obs = [
        TopologyCalibrationObservation("A_B", 1.0),
        TopologyCalibrationObservation("B_C", 1.0),
        TopologyCalibrationObservation("C_A", 0.0),
    ]
    
    # Samples with heads and layers
    samples = [
        {"site_id": "A", "hydraulic_head": 10.0, "aquifer_layer": 0, "Ca": 1.0, "Mg": 0.5, "Na": 0.2, "K": 0.1, "HCO3": 1.5, "Cl": 0.5, "SO4": 0.3, "NO3": 0.1, "F": 0.01, "Fe": 0.01, "PO4": 0.01},
        {"site_id": "B", "hydraulic_head": 8.0, "aquifer_layer": 0, "Ca": 1.0, "Mg": 0.5, "Na": 0.2, "K": 0.1, "HCO3": 1.5, "Cl": 0.5, "SO4": 0.3, "NO3": 0.1, "F": 0.01, "Fe": 0.01, "PO4": 0.01},
        {"site_id": "C", "hydraulic_head": 6.0, "aquifer_layer": 1, "Ca": 1.0, "Mg": 0.5, "Na": 0.2, "K": 0.1, "HCO3": 1.5, "Cl": 0.5, "SO4": 0.3, "NO3": 0.1, "F": 0.01, "Fe": 0.01, "PO4": 0.01},
    ]

    calibrator = TopologyOuterLoopCalibrator(
        samples_df=samples,
        candidate_edges=[edge_a_b, edge_b_c, edge_c_a],
        observations=obs,
        max_iterations=3,
        max_neighbors=1,
    )
    
    # Cycle check
    assert calibrator.check_constraints([edge_a_b, edge_b_c]) is True
    # If we add C_A, it creates a cycle A->B->C->A, so check_constraints must return False
    assert calibrator.check_constraints([edge_a_b, edge_b_c, edge_c_a]) is False

    # Run discrete search
    res = calibrator.search()
    assert res["success"] is True
    # Cycle should be broken, so C_A should not be selected
    assert edge_c_a not in res["selected_edges"]


def test_huber_robust_loss():
    # Simple linear problem with an outlier
    # y = 2x
    # Data points: x = [1, 2, 3], y = [2, 4, 20] (outlier)
    def model_runner(params):
        x_val = params["coeff"]
        return {
            "obs1": 1.0 * x_val,
            "obs2": 2.0 * x_val,
            "obs3": 3.0 * x_val,
        }

    params = [AdjustableParameter("coeff", 1.0, 0.0, 10.0)]
    obs = [
        Observation("obs1", 2.0, weight=1.0),
        Observation("obs2", 4.0, weight=1.0),
        Observation("obs3", 20.0, weight=1.0),  # Outlier
    ]

    # Run with linear loss
    pest_linear = PESTGLM(params, obs, model_runner, loss="linear")
    res_linear = pest_linear.calibrate(max_nfev=20)
    coeff_linear = res_linear["optimal_parameters"]["coeff"]

    # Run with Huber loss
    pest_huber = PESTGLM(params, obs, model_runner, loss="huber")
    res_huber = pest_huber.calibrate(max_nfev=20)
    coeff_huber = res_huber["optimal_parameters"]["coeff"]

    # Huber should be less affected by the outlier at obs3, so coeff should be closer to 2.0 than linear
    assert abs(coeff_huber - 2.0) < abs(coeff_linear - 2.0)


def test_vadose_adapter_parameters():
    # Setup dummy VadoseProfile with layers
    from dataclasses import dataclass
    
    @dataclass
    class DummyLayer:
        ks_m_day: float = 1.0
        alpha_1_m: float = 0.5
        n: float = 1.5

    @dataclass
    class DummyProfile:
        layers: list

    @dataclass
    class DummyConfig:
        kc: float = 1.0
        preferential_flow_fraction: float = 0.1
        ttd_default_cv: float = 0.7

    profile = DummyProfile(layers=[DummyLayer()])
    config = DummyConfig()
    
    adapter = VadoseCalibrationAdapter(
        profile=profile,
        forcing=[],
        observations=[],
        config=config,
        layers_to_fit=[0],
        params_to_fit=["ks_L0", "alpha_L0", "n_L0", "kc", "preferential_flow_fraction", "ttd_cv"],
    )

    params = adapter.get_parameters()
    param_names = [p.name for p in params]
    assert "ks_L0" in param_names
    assert "alpha_L0" in param_names
    assert "n_L0" in param_names
    assert "kc" in param_names
    assert "preferential_flow_fraction" in param_names
    assert "ttd_cv" in param_names


def test_nitrate_adapter_process_constraints():
    nodes_df = pd.DataFrame(
        [{"site_id": "N1", "NO3": 1.0, "Cl": 0.1, "K": 0.1, "PO4": 0.01, "Fe": 0.001}]
    )
    
    adapter = NitrateSourceCalibrationAdapter(
        nodes_df=nodes_df,
        observations=[NitrateSourceCalibrationObservation("N1", 0.7)],
        params_to_fit=[
            "denitrification_min_extent",
            "min_top_probability",
            "nitrate_source_min_mg_L",
        ],
    )
    
    params = adapter.get_parameters()
    param_names = [p.name for p in params]
    assert "denitrification_min_extent" in param_names
    assert "min_top_probability" in param_names
    assert "nitrate_source_min_mg_L" in param_names


def test_age_temporal_node_pairs():
    # Test age temporal adapter with node pair files
    # Mock data times and values
    u_times = [0.0, 10.0, 20.0, 30.0]
    u_vals = [10.0, 15.0, 20.0, 10.0]
    v_times = [10.0, 20.0, 30.0, 40.0]
    v_vals = [5.0, 7.5, 10.0, 5.0]

    exp = AgeTemporalExperiment(
        id="test_exp",
        node_u_times=u_times,
        node_u_values=u_vals,
        node_v_times=v_times,
        node_v_values=v_vals,
        default_tau_days=10.0,
    )

    adapter = AgeTemporalCalibrationAdapter(
        experiments=[exp],
        params_to_fit=["tau", "ttd_cv", "decay"],
        base_decay_rate_1_day=0.0693147,  # half life 10 days
        base_ttd_cv=0.1,
    )

    params = adapter.get_parameters()
    param_names = [p.name for p in params]
    assert "tau_days" in param_names
    assert "ttd_cv" in param_names
    assert "decay_rate_1_day" in param_names

    # Run model simulation
    res = adapter.run_model({"tau_days": 10.0, "ttd_cv": 0.01, "decay_rate_1_day": 0.0693147})
    # For CV close to 0, it behaves like pure advection with decay: C_v(t) = C_u(t - tau) * exp(-decay * tau)
    # At t=20.0, C_u(10.0) = 15.0. Decay = exp(-0.0693147 * 10) = 0.5. Expected C_v(20.0) = 15.0 * 0.5 = 7.5
    assert np.isclose(res["test_exp_v_1"], 7.5, rtol=1e-2)


def test_cli_write_template(capsys):
    from hydrosheaf.calibration.cli import main
    import sys

    with patch.object(sys, "argv", ["hydrosheaf-cal", "--write-template", "composite"]):
        exit_code = main()
        assert exit_code == 0
        captured = capsys.readouterr()
        assert "calibration:" in captured.out
        assert "type: composite" in captured.out


def test_cli_dry_run(tmp_path, capsys):
    from hydrosheaf.calibration.cli import main
    import sys
    import yaml

    config_data = {
        "calibration": {
            "type": "transport",
            "settings": {
                "engine": "internal",
                "max_iterations": 10,
            },
            "parameters": [
                {"name": "coeff", "initial": 1.0, "bounds": [0.0, 5.0]}
            ]
        }
    }
    
    config_file = tmp_path / "cal_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    with patch.object(sys, "argv", ["hydrosheaf-cal", str(config_file), "--dry-run"]):
        exit_code = main()
        assert exit_code == 0
        captured = capsys.readouterr()
        assert "DRY RUN: RESOLVED CALIBRATION TARGETS" in captured.out
        assert "coeff" in captured.out
