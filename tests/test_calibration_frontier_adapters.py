from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from hydrosheaf.calibration.adapters import (
    AgeTemporalCalibrationAdapter,
    AgeTemporalExperiment,
    GenericFunctionAdapter,
    KineticCalibrationAdapter,
    KineticExperiment,
    NitrateSourceCalibrationAdapter,
    NitrateSourceCalibrationObservation,
    TopologyCalibrationAdapter,
    TopologyCalibrationObservation,
)
from hydrosheaf.calibration.cli import setup_vadose_adapter
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM
from hydrosheaf.graph.types import Edge
from hydrosheaf.reactive_transport import KineticParameters


def test_kinetic_adapter_supports_hierarchical_scopes_and_extent_parameters():
    exp_a = KineticExperiment(
        id="exp_a",
        edge_id="edge_A",
        geological_layer="layer_1",
        initial_solution={},
        residence_time_days=1.0,
        reaction_extents=[0.0],
        reaction_labels=["calcite"],
        observations={},
    )
    exp_b = KineticExperiment(
        id="exp_b",
        edge_id="edge_B",
        geological_layer="layer_2",
        initial_solution={},
        residence_time_days=1.0,
        reaction_extents=[0.0],
        reaction_labels=["calcite"],
        observations={},
    )
    adapter = KineticCalibrationAdapter(
        base_params={"calcite": KineticParameters("calcite", 1e-6, 1.0)},
        experiments=[exp_a, exp_b],
        config=SimpleNamespace(ion_order=["Ca"]),
        params_to_fit=["calcite:k:per_layer", "calcite:extent:per_edge"],
    )

    names = {param.name for param in adapter.get_parameters()}
    assert "log_k_calcite__layer_1" in names
    assert "log_k_calcite__layer_2" in names
    assert "extent_calcite__edge_A" in names
    assert "extent_calcite__edge_B" in names

    values = {
        "log_k_calcite__layer_1": 1e-5,
        "log_k_calcite__layer_2": 1e-4,
        "extent_calcite__edge_A": 2.5,
        "extent_calcite__edge_B": 5.0,
    }
    assert np.isclose(adapter._resolve_parameters(exp_a, values)["calcite"].rate_constant, 1e-5)
    assert np.isclose(adapter._resolve_parameters(exp_b, values)["calcite"].rate_constant, 1e-4)
    assert np.isclose(adapter._resolve_extents(exp_a, values)[0], 2.5)
    assert np.isclose(adapter._resolve_extents(exp_b, values)[0], 5.0)


def test_age_temporal_adapter_can_fit_global_tau():
    adapter = AgeTemporalCalibrationAdapter(
        experiments=[AgeTemporalExperiment(id="age1", observed_age_days=20.0, default_tau_days=5.0)],
        params_to_fit=["tau"],
    )
    pest = PESTGLM.from_problem(adapter)

    result = pest.calibrate(max_nfev=30)

    assert result["success"]
    assert np.isclose(result["optimal_parameters"]["tau_days"], 20.0, rtol=1e-3)


def test_age_temporal_adapter_uses_decay_for_tracer_age():
    adapter = AgeTemporalCalibrationAdapter(
        experiments=[
            AgeTemporalExperiment(
                id="tritium",
                observed_age_days=10.0,
                tracer_initial=100.0,
                tracer_observed=50.0,
            )
        ],
        params_to_fit=["decay"],
        base_decay_rate_1_day=0.069314718,
    )

    result = adapter.run_model({"decay_rate_1_day": 0.069314718})

    assert np.isclose(result["tritium_age_days"], 10.0, rtol=1e-5)


def test_topology_adapter_calibrates_soft_edge_presence():
    edge = Edge(edge_id="A_B", u="A", v="B", attrs={"edge_confidence": 0.25})
    adapter = TopologyCalibrationAdapter(
        candidate_edges=[edge],
        observations=[TopologyCalibrationObservation("A_B", 1.0)],
    )

    params = adapter.get_parameters()
    result = adapter.run_model({params[0].name: 6.0})

    assert result["topology_A_B"] > 0.99
    assert adapter.selected_edges({params[0].name: 6.0}) == [edge]


def test_nitrate_adapter_passes_calibrated_overrides_to_inference():
    nodes_df = pd.DataFrame(
        [{"site_id": "N1", "NO3": 1.0, "Cl": 0.1, "K": 0.1, "PO4": 0.01, "Fe": 0.001}]
    )
    adapter = NitrateSourceCalibrationAdapter(
        nodes_df=nodes_df,
        observations=[NitrateSourceCalibrationObservation("N1", 0.7)],
        params_to_fit=[
            "prior_logit",
            "w1_no3_cl",
            "denitrification_min_extent",
            "min_top_probability",
        ],
    )

    with patch("hydrosheaf.nitrate_source_v2.infer_node_posteriors") as mock_infer:
        mock_infer.return_value = {
            "N1": SimpleNamespace(p_manure=0.7, p_fertilizer=0.3)
        }
        result = adapter.run_model(
            {
                "prior_logit": 0.0,
                "w1_no3_cl": 2.0,
                "denitrification_min_extent": 0.2,
                "min_top_probability": 0.8,
            }
        )

    overrides = mock_infer.call_args.kwargs["config_overrides"]
    assert np.isclose(overrides["prior_prob"], 0.5)
    assert overrides["weights"]["w1_no3_cl"] == 2.0
    assert overrides["isotope_process_constraints"]["denitrification_min_extent"] == 0.2
    assert overrides["isotope_qc"]["min_top_probability"] == 0.8
    assert result["nitrate_N1_p_manure"] == 0.7


def test_setup_vadose_adapter_loads_files(tmp_path):
    profile = tmp_path / "profile.json"
    profile.write_text(
        """
{
  "profile_id": "test",
  "depth_m": 1.0,
  "layers": [
    {"thickness_m": 1.0, "theta_r": 0.05, "theta_s": 0.45, "alpha_1_m": 1.0, "n": 1.5, "ks_m_day": 0.2}
  ]
}
""".strip(),
        encoding="utf-8",
    )
    forcing = tmp_path / "forcing.csv"
    forcing.write_text(
        "timestamp,P_mm_day,ET0_mm_day\n2020-01-01,1.0,0.2\n2020-01-02,0.0,0.2\n",
        encoding="utf-8",
    )
    observations = tmp_path / "obs.csv"
    observations.write_text(
        "timestamp,depth_m,theta\n2020-01-01,0.5,0.25\n",
        encoding="utf-8",
    )
    config = SimpleNamespace(
        adapter_settings={
            "profile_file": str(profile),
            "forcing_file": str(forcing),
            "layers_to_fit": [0],
        },
        observations_file=str(observations),
        model_config_file=None,
    )

    adapter = setup_vadose_adapter(config)

    assert [p.name for p in adapter.get_parameters()] == ["ks_L0", "alpha_L0"]
    assert adapter.get_observations()[0].name == "obs_0_theta"


def test_glm_uses_svd_covariance_for_rank_deficient_jacobian():
    def model(params):
        return {"obs": params["x"] + params["y"]}

    adapter = GenericFunctionAdapter(
        model,
        [
            AdjustableParameter("x", 0.0, -10.0, 10.0),
            AdjustableParameter("y", 0.0, -10.0, 10.0),
        ],
        [Observation("obs", 1.0)],
    )
    result = PESTGLM.from_problem(adapter).calibrate(max_nfev=20)

    assert result["success"]
    assert result["covariance_method"] == "svd_pseudoinverse"
