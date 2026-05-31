import math

import pytest

from hydrosheaf.config import Config
from hydrosheaf.graph.types import Edge
from hydrosheaf.inference.edge_fit import fit_edge
from hydrosheaf.inference.topology_posterior import run_topology_posterior
from hydrosheaf.phreeqc.constraints import build_edge_bounds


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
