"""Unit tests for non-parametric single-node and network TTD inverse solver."""

import networkx as nx
import numpy as np
import pytest

from hydrosheaf.nuclear.ttd_diagnostics import (
    CODE_INSUFFICIENT_TRACERS,
    CODE_GRAPH_CYCLES,
    CODE_HEAD_CONFLICT,
    CODE_MISSING_FORWARD_SYSTEM,
    CODE_MISSING_MIXING_SPEC,
    CODE_TRACER_CONFLICT,
)
from hydrosheaf.nuclear.ttd_grid import build_uniform_ttd_grid
from hydrosheaf.nuclear.ttd_kernel_builder import (
    NodeTracerPanel,
    TracerObservation,
    build_forward_system,
)
from hydrosheaf.nuclear.ttd_network_solver import (
    solve_network_ttd,
    solve_single_node_ttd,
)
from hydrosheaf.nuclear.ttd_transport import (
    NodeMixingSpecification,
    build_advection_dispersion_operator,
    build_local_recharge_distribution,
)


def test_single_node_exponential_recovery():
    grid = build_uniform_ttd_grid(max_age_years=60.0, dt_years=1.0)
    # Synthetic ground truth: Exponential model with mean age = 15 years
    tau_true = 15.0
    g_true = (1.0 / tau_true) * np.exp(-grid.taus / tau_true)
    g_true = g_true / np.sum(g_true)

    # Build forward panel with 3H, SF6, 14C
    dummy_panel = NodeTracerPanel(
        node_id="W1",
        sample_year=2024.0,
        observations=(
            TracerObservation("3H", 1.0, 0.2),
            TracerObservation("SF6", 1.0, 0.2),
            TracerObservation("14C", 1.0, 1.0),
        ),
    )
    sys = build_forward_system(dummy_panel, grid)

    # Synthetic observations: c_synth = A @ g_true
    c_synth = sys.predict(g_true)
    panel_synth = NodeTracerPanel(
        node_id="W1",
        sample_year=2024.0,
        observations=(
            TracerObservation("3H", c_synth[0], 0.2),
            TracerObservation("SF6", c_synth[1], 0.2),
            TracerObservation("14C", c_synth[2], 1.5),
        ),
    )
    sys_synth = build_forward_system(panel_synth, grid)

    res = solve_single_node_ttd(sys_synth, grid, lambda_smoothness=0.01)

    assert res.status == "ESTIMATED"
    assert np.sum(res.g) == pytest.approx(1.0, abs=1e-5)
    # Recovered young-water fraction should be close to ground truth
    yf_true = grid.young_water_fraction(g_true)
    assert res.young_water_fraction == pytest.approx(yf_true, abs=0.08)
    assert res.chi_squared_per_observation < 1.0


def test_network_branch_merge_inversion():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=1.0)

    # Graph: U1 -> V, U2 -> V (branch & merge)
    graph = nx.DiGraph()
    graph.add_edge("U1", "V")
    graph.add_edge("U2", "V")

    # Heads: U1=100m, U2=95m, V=80m (valid hydraulic flow)
    heads = {"U1": 100.0, "U2": 95.0, "V": 80.0}

    # True distributions:
    # U1: Young water (mean 5 yr)
    g_u1_true = np.exp(-grid.taus / 5.0)
    g_u1_true /= np.sum(g_u1_true)

    # U2: Intermediate water (mean 20 yr)
    g_u2_true = np.exp(-grid.taus / 20.0)
    g_u2_true /= np.sum(g_u2_true)

    # Operators: U1 -> V (+5 yr), U2 -> V (+5 yr)
    op1 = build_advection_dispersion_operator(("U1", "V"), grid, delta_tau_years=5.0, dispersion=0.02)
    op2 = build_advection_dispersion_operator(("U2", "V"), grid, delta_tau_years=5.0, dispersion=0.02)
    transport_ops = {("U1", "V"): op1, ("U2", "V"): op2}

    # Downstream V true arrival: 50% U1 routed, 50% U2 routed
    g_v_true = 0.5 * op1.transport(g_u1_true) + 0.5 * op2.transport(g_u2_true)

    # Forward systems
    p_u1 = NodeTracerPanel("U1", 2024.0, (TracerObservation("3H", 1.0, 0.2), TracerObservation("14C", 1.0, 1.0)), head_m=100.0)
    p_u2 = NodeTracerPanel("U2", 2024.0, (TracerObservation("3H", 1.0, 0.2), TracerObservation("14C", 1.0, 1.0)), head_m=95.0)
    p_v = NodeTracerPanel("V", 2024.0, (TracerObservation("3H", 1.0, 0.2), TracerObservation("14C", 1.0, 1.0)), head_m=80.0)

    s_u1 = build_forward_system(p_u1, grid)
    s_u2 = build_forward_system(p_u2, grid)
    s_v = build_forward_system(p_v, grid)

    # Replace with synthetic observations
    c_u1 = s_u1.predict(g_u1_true)
    c_u2 = s_u2.predict(g_u2_true)
    c_v = s_v.predict(g_v_true)

    fwd_systems = {
        "U1": build_forward_system(NodeTracerPanel("U1", 2024.0, (TracerObservation("3H", c_u1[0], 0.2), TracerObservation("14C", c_u1[1], 1.0)), head_m=100.0), grid),
        "U2": build_forward_system(NodeTracerPanel("U2", 2024.0, (TracerObservation("3H", c_u2[0], 0.2), TracerObservation("14C", c_u2[1], 1.0)), head_m=95.0), grid),
        "V": build_forward_system(NodeTracerPanel("V", 2024.0, (TracerObservation("3H", c_v[0], 0.2), TracerObservation("14C", c_v[1], 1.0)), head_m=80.0), grid),
    }

    net_res = solve_network_ttd(
        graph,
        fwd_systems,
        transport_ops,
        grid,
        mixing_specs={
            "V": NodeMixingSpecification(
                node_id="V",
                local_fraction=0.0,
                upstream_weights={"U1": 0.5, "U2": 0.5},
                recharge_distribution=build_local_recharge_distribution(grid),
            )
        },
        lambda_smoothness=0.01,
        lambda_graph=2.0,
        node_heads=heads,
    )

    assert net_res.status == "ESTIMATED"
    for nid in ["U1", "U2", "V"]:
        g_rec = net_res.node_g(nid)
        assert g_rec is not None
        assert np.sum(g_rec) == pytest.approx(1.0, abs=1e-4)

    # Downstream mean age at V should exceed U1 mean age
    mtts = net_res.mean_transit_times()
    assert mtts["V"] > mtts["U1"]


def test_abstention_on_graph_cycle():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    # Graph with cycle: A -> B -> C -> A
    graph = nx.DiGraph()
    graph.add_edge("A", "B")
    graph.add_edge("B", "C")
    graph.add_edge("C", "A")

    fwd_systems = {
        "A": build_forward_system(NodeTracerPanel("A", 2024.0, (TracerObservation("3H", 2.0, 0.2),)), grid),
        "B": build_forward_system(NodeTracerPanel("B", 2024.0, (TracerObservation("3H", 2.0, 0.2),)), grid),
        "C": build_forward_system(NodeTracerPanel("C", 2024.0, (TracerObservation("3H", 2.0, 0.2),)), grid),
    }

    res = solve_network_ttd(graph, fwd_systems, {}, grid)
    assert res.status == "ABSTAIN"
    assert CODE_GRAPH_CYCLES in res.abstention_reasons


def test_abstention_on_adverse_heads():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    graph = nx.DiGraph([("A", "B")])
    # Adverse heads: A=50m, B=100m (water cannot flow uphill by 50m)
    heads = {"A": 50.0, "B": 100.0}

    fwd_systems = {
        "A": build_forward_system(NodeTracerPanel("A", 2024.0, (TracerObservation("3H", 2.0, 0.2),)), grid),
        "B": build_forward_system(NodeTracerPanel("B", 2024.0, (TracerObservation("3H", 2.0, 0.2),)), grid),
    }

    res = solve_network_ttd(graph, fwd_systems, {}, grid, node_heads=heads)
    assert res.status == "ABSTAIN"
    assert CODE_HEAD_CONFLICT in res.abstention_reasons


def test_abstention_on_tracer_contradiction():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    # Contradictory panel: Modern SF6 (4.5 pptv) with 0 tritium (0.0 TU) and dead 14C (2 pmc)
    panel = NodeTracerPanel(
        "ErrWell",
        2024.0,
        observations=(
            TracerObservation("SF6", 4.5, 0.2),
            TracerObservation("3H", 0.05, 0.05),
            TracerObservation("14C", 2.0, 0.5),
        ),
    )
    sys = build_forward_system(panel, grid)
    res = solve_single_node_ttd(sys, grid)
    assert res.status == "ABSTAIN"
    assert CODE_TRACER_CONFLICT in res.abstention_reasons


def test_network_abstains_on_missing_forward_system():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    graph = nx.DiGraph([("A", "B")])
    panel = NodeTracerPanel(
        "A",
        2024.0,
        (TracerObservation("3H", 2.0, 0.2), TracerObservation("14C", 50.0, 1.0)),
    )
    system = build_forward_system(panel, grid)
    op = build_advection_dispersion_operator(
        ("A", "B"), grid, delta_tau_years=1.0, dispersion=0.0
    )

    result = solve_network_ttd(
        graph,
        {"A": system},
        {("A", "B"): op},
        grid,
    )

    assert result.status == "ABSTAIN"
    assert CODE_MISSING_FORWARD_SYSTEM in result.abstention_reasons


def test_network_abstains_on_node_tracer_contradiction():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    graph = nx.DiGraph()
    graph.add_node("ErrWell")
    system = build_forward_system(
        NodeTracerPanel(
            "ErrWell",
            2024.0,
            (
                TracerObservation("SF6", 4.5, 0.2),
                TracerObservation("3H", 0.05, 0.05),
                TracerObservation("14C", 2.0, 0.5),
            ),
        ),
        grid,
    )

    result = solve_network_ttd(graph, {"ErrWell": system}, {}, grid)

    assert result.status == "ABSTAIN"
    assert CODE_TRACER_CONFLICT in result.abstention_reasons


def test_network_abstains_on_missing_mixing_specification():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    graph = nx.DiGraph([("A", "B")])
    systems = {
        node_id: build_forward_system(
            NodeTracerPanel(
                node_id,
                2024.0,
                (
                    TracerObservation("3H", 2.0, 0.2),
                    TracerObservation("14C", 50.0, 1.0),
                ),
            ),
            grid,
        )
        for node_id in ("A", "B")
    }
    op = build_advection_dispersion_operator(
        ("A", "B"), grid, delta_tau_years=1.0, dispersion=0.0
    )

    result = solve_network_ttd(graph, systems, {("A", "B"): op}, grid)

    assert result.status == "ABSTAIN"
    assert CODE_MISSING_MIXING_SPEC in result.abstention_reasons


def test_single_node_abstains_with_one_tracer():
    grid = build_uniform_ttd_grid(max_age_years=30.0, dt_years=1.0)
    system = build_forward_system(
        NodeTracerPanel(
            "A",
            2024.0,
            (TracerObservation("3H", 2.0, 0.2),),
        ),
        grid,
    )

    result = solve_single_node_ttd(system, grid)

    assert result.status == "ABSTAIN"
    assert CODE_INSUFFICIENT_TRACERS in result.abstention_reasons
