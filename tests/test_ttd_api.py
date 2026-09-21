"""Public API integration tests for the explicit network TTD workflow."""

import networkx as nx
import numpy as np

from hydrosheaf import (
    NodeMixingSpecification,
    NodeTracerPanel,
    TracerObservation,
    build_advection_dispersion_operator,
    build_forward_system,
    build_local_recharge_distribution,
    build_uniform_ttd_grid,
    fit_ttd_network,
)


def test_fit_ttd_network_builds_forward_systems_and_solves():
    grid = build_uniform_ttd_grid(max_age_years=20.0, dt_years=1.0)
    graph = nx.DiGraph([("A", "B")])
    operator = build_advection_dispersion_operator(
        ("A", "B"), grid, delta_tau_years=1.0, dispersion=0.0
    )

    g_a = np.exp(-grid.taus / 4.0)
    g_a /= np.sum(g_a)
    g_b = operator.transport(g_a)
    template_a = NodeTracerPanel(
        "A",
        2024.0,
        (
            TracerObservation("3H", 1.0, 0.2),
            TracerObservation("14C", 50.0, 1.0),
        ),
        head_m=100.0,
    )
    template_b = NodeTracerPanel(
        "B",
        2024.0,
        (
            TracerObservation("3H", 1.0, 0.2),
            TracerObservation("14C", 50.0, 1.0),
        ),
        head_m=90.0,
    )
    system_a = build_forward_system(template_a, grid)
    system_b = build_forward_system(template_b, grid)
    c_a = system_a.predict(g_a)
    c_b = system_b.predict(g_b)

    panels = [
        NodeTracerPanel(
            "A",
            2024.0,
            (
                TracerObservation("3H", c_a[0], 0.2),
                TracerObservation("14C", c_a[1], 1.0),
            ),
            head_m=100.0,
        ),
        NodeTracerPanel(
            "B",
            2024.0,
            (
                TracerObservation("3H", c_b[0], 0.2),
                TracerObservation("14C", c_b[1], 1.0),
            ),
            head_m=90.0,
        ),
    ]

    result = fit_ttd_network(
        graph,
        panels,
        grid,
        transport_operators={("A", "B"): operator},
        mixing_specs={
            "B": NodeMixingSpecification(
                node_id="B",
                local_fraction=0.0,
                upstream_weights={"A": 1.0},
                recharge_distribution=build_local_recharge_distribution(grid),
            )
        },
        lambda_smoothness=0.01,
        lambda_graph=2.0,
    )

    assert result.status == "ESTIMATED"
    assert result.node_g("A") is not None
    assert result.node_g("B") is not None
    assert np.sum(result.node_g("B")) == 1.0
