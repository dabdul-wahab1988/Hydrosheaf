import math

import pytest

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.graph3d.build_3d import build_network_3d
from hydrosheaf.inference.network_fit import estimate_edge_residence_time_days
from hydrosheaf.inference.topology_posterior import run_topology_posterior
from hydrosheaf.isotopes import craig_gordon_enrichment
from hydrosheaf.models import nitrate_isotopes
from hydrosheaf.nuclear.input_history import build_default_tritium_input
from hydrosheaf.nuclear.lpm import convolve_input
from hydrosheaf.phreeqc.constraints import build_edge_bounds
from hydrosheaf.temporal.convolution import convolve_reactive_ttd
from hydrosheaf.vadose.ttd import mixture_ttd_from_series


def test_evaporation_does_not_fit_dilution():
    cfg = Config(transport_models_enabled=["evap"], lambda_sparse=0.0)
    x_u = [1.0 + 0.1 * i for i in range(len(cfg.ion_order))]
    x_v = [0.5 * value for value in x_u]

    result = fit_edge(x_u, x_v, cfg, edge_id="A->B", u="A", v="B")

    assert result.transport_model == "evap"
    assert result.gamma == pytest.approx(1.0)


def test_transport_objective_includes_transport_residual():
    cfg = Config(transport_models_enabled=["evap"], lambda_sparse=0.0)
    x_u = [1.0 + 0.1 * i for i in range(len(cfg.ion_order))]
    x_v = [1.2 * value for value in x_u]

    result = fit_edge(x_u, x_v, cfg, edge_id="A->B", u="A", v="B")
    entry = result.candidate_scores[0]

    assert entry["transport_residual_norm"] >= 0.0
    assert entry["combined_residual_norm"] >= entry["transport_residual_norm"]
    assert entry["objective_score"] >= entry["combined_residual_norm"]
    assert entry["information_score"] is not None


def test_topology_posterior_recovers_independent_edge_prior():
    cfg = Config()
    cfg.topology_posterior_samples = 12000
    cfg.topology_posterior_burnin = 1000
    cfg.topology_posterior_beta = 0.0
    cfg.topology_posterior_edge_penalty = 0.0
    universe = [
        Edge(edge_id=f"E{i}", u="A", v=f"B{i}", attrs={"edge_confidence": 0.8})
        for i in range(3)
    ]

    result = run_topology_posterior(
        universe, lambda edges: 0.0, cfg, initial_edges=[], seed=7
    )

    for prob in result["edge_probabilities"].values():
        assert prob == pytest.approx(0.8, abs=0.03)


def test_phreeqc_allows_dissolution_to_saturation():
    cfg = Config(phreeqc_enabled=True)
    edge = Edge(edge_id="A->B", u="A", v="B", attrs={})
    phreeqc_results = {
        "A": {"phreeqc_ok": True, "si_calcite": -1.0},
        "B": {"phreeqc_ok": True, "si_calcite": 0.0},
    }

    bounds = build_edge_bounds(
        phreeqc_results,
        [edge],
        labels=["calcite"],
        mineral_mask=[True],
        config=cfg,
    )["A->B"]

    assert bounds["constraints_active"]["calcite"] == "dissolution_only"
    assert bounds["lb"][0] == 0.0
    assert math.isinf(bounds["ub"][0])


def test_requested_bootstrap_uncertainty_is_attached_to_edge_result():
    cfg = Config(
        transport_models_enabled=["evap"],
        uncertainty_method="bootstrap",
        bootstrap_n_resamples=20,
        bootstrap_ci_method="bca",
        lambda_sparse=0.0,
    )
    x_u = [1.0 + 0.2 * i for i in range(len(cfg.ion_order))]
    x_v = [
        1.15 * value + (0.02 if idx % 2 else -0.015)
        for idx, value in enumerate(x_u)
    ]

    result = fit_edge(x_u, x_v, cfg, edge_id="A->B", u="A", v="B")

    assert result.uncertainty is not None
    assert result.uncertainty_method == "bootstrap"
    assert result.gamma_std is not None
    assert result.gamma_ci_low is not None
    assert result.gamma_ci_high is not None


def test_temporal_convolution_uses_lagged_input_series():
    c_out = convolve_reactive_ttd(
        999.0,
        [0.5, 0.5],
        [0.0, 5.0],
        0.0,
        input_times=[0.0, 10.0],
        input_values=[0.0, 10.0],
        t_sample=10.0,
    )

    assert c_out == pytest.approx(7.5)


def test_isotope_mixture_uncertainty_uses_fraction_squared_weights():
    sources = [
        nitrate_isotopes.SourceIsotopes("A", 0.0, 1.0, 0.0, 2.0),
        nitrate_isotopes.SourceIsotopes("B", 10.0, 1.0, 0.0, 2.0),
    ]

    moments = nitrate_isotopes._summary_moments({"A": 0.5, "B": 0.5}, sources)

    assert moments["mean_d15N"] == pytest.approx(5.0)
    assert moments["var_d15N"] == pytest.approx(0.5)
    assert moments["var_d18O"] == pytest.approx(2.0)


def test_3d_edges_propagate_zero_head_delta_h():
    cfg = Config(edge_radius_km=1.0, edge_p_min=0.0, edge_max_neighbors=2)
    samples = [
        {"site_id": "A", "x": 0.0, "y": 0.0, "z": 10.0, "elevation": 100.0, "head": 0.0},
        {"site_id": "B", "x": 100.0, "y": 0.0, "z": 10.0, "elevation": 100.0, "head": -1.0},
    ]

    network = build_network_3d(samples, cfg, use_haversine=False)
    edge = next(e for e in network.edges if e.edge_id == "A->B")

    assert edge.delta_h == pytest.approx(1.0)
    assert edge.horizontal_gradient == pytest.approx(0.01)


def test_residence_time_prefers_delta_h_over_vertical_gradient_fallback():
    cfg = Config(residence_time_hydraulic_k=1.0, residence_time_porosity=0.1)
    tau_mean, _ = estimate_edge_residence_time_days(
        {"distance_3d_m": 100.0, "delta_h": 1.0, "vertical_gradient": 1000.0},
        cfg,
    )

    assert tau_mean == pytest.approx(1000.0)


def test_dual_domain_ttd_summary_matches_effective_kernel():
    _, _, single = mixture_ttd_from_series(
        [100.0] * 10,
        [1.0] * 10,
        grid_dt_days=1.0,
        max_lag_days=200.0,
        cv=0.1,
        preferential_flow_fraction=0.0,
    )
    _, _, dual = mixture_ttd_from_series(
        [100.0] * 10,
        [1.0] * 10,
        grid_dt_days=1.0,
        max_lag_days=200.0,
        cv=0.1,
        preferential_flow_fraction=0.2,
        preferential_velocity_factor=10.0,
    )

    assert dual.mean_days < single.mean_days
    assert dual.mean_days == pytest.approx(82.0, abs=2.0)


def test_lpm_prebomb_baseline_uses_input_history_value():
    hist = build_default_tritium_input()

    out = convolve_input(
        1900.0,
        0.0,
        hist.years,
        hist.values,
        0.0,
        model_type="PFM",
    )

    assert out == pytest.approx(float(hist.values[0]))


def test_craig_gordon_enrichment_is_not_noop():
    out = craig_gordon_enrichment(
        delta_L=-5.0,
        delta_A=-20.0,
        h=0.5,
        epsilon_eq=9.8,
        epsilon_k=14.2,
    )

    assert out != pytest.approx(-5.0)
    assert out == pytest.approx((-5.0 - 0.5 * -20.0 - 24.0) / 0.5142)
