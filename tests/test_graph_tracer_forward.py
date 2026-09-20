"""Tests for causal dynamic graph transport and convolution forward operator."""

from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.nuclear.dynamic_edge_kernel import DynamicEdgeKernel
from hydrosheaf.nuclear.graph_tracer_forward import (
    convolve_causal_dynamic_kernel,
    simulate_dynamic_graph_transport,
    topological_sort,
    validate_candidate_graph,
)


def test_topological_sort():
    nodes = ("C", "A", "B", "R")
    edges = (("R", "A"), ("A", "B"), ("B", "C"), ("A", "C"))
    order = topological_sort(nodes, edges)
    assert order.index("R") < order.index("A")
    assert order.index("A") < order.index("B")
    assert order.index("B") < order.index("C")
    assert order.index("A") < order.index("C")

    # Cyclic graph must fail
    cyclic_edges = (("A", "B"), ("B", "A"))
    with pytest.raises(ValueError, match="cycles"):
        topological_sort(("A", "B"), cyclic_edges)


def test_convolve_causal_dynamic_kernel_pure_delay():
    T = 20
    source = np.sin(np.linspace(0, 4 * np.pi, T))
    L = 5
    lags = np.arange(L, dtype=float)

    # Pure delay of 2 steps: kernel has 1.0 at lag index 2
    k_mat = np.zeros((T, L), dtype=float)
    k_mat[:, 2] = 1.0

    out = convolve_causal_dynamic_kernel(source, k_mat, lags)
    # Output at t should equal source[t - 2] for t >= 2
    for t in range(2, T):
        assert pytest.approx(out[t]) == source[t - 2]


def test_convolution_uses_declared_physical_lag_units():
    source = np.arange(12, dtype=float)
    kernel = np.zeros((12, 3), dtype=float)
    kernel[:, 1] = 1.0
    out = convolve_causal_dynamic_kernel(source, kernel, np.array([0.0, 7.0, 14.0]), step_days=7.0)
    np.testing.assert_allclose(out[1:], source[:-1])
    with pytest.raises(ValueError, match="align"):
        convolve_causal_dynamic_kernel(source, kernel, np.array([0.0, 5.0, 10.0]), step_days=7.0)


def test_candidate_graph_boundary_gate_is_data_structural():
    assert validate_candidate_graph(("R", "A", "B"), (("R", "A"), ("A", "B")))[0]
    ok, reason = validate_candidate_graph(("R", "A", "B"), (("A", "R"), ("B", "A")))
    assert ok is False
    assert "incoming" in reason


def test_simulate_dynamic_graph_transport_mass_conservation():
    T = 30
    times = np.arange(T, dtype=float)
    lags = np.array([0.0, 1.0, 2.0, 3.0])
    nodes = ("R", "A", "B")
    edges = (("R", "A"), ("A", "B"))

    # Uniform kernel on both edges
    k_vals = np.full((2, T, 4), 0.25)
    dek = DynamicEdgeKernel(
        edge_ids=("R->A", "A->B"),
        time_grid=times,
        lag_grid=lags,
        values=k_vals,
    )

    # Constant input at R = 10.0
    local_inputs = {"R": np.full(T, 10.0), "A": np.full(T, 0.0), "B": np.full(T, 0.0)}
    # Pure transport: rho_A = 0.0, pi_RA = 1.0; rho_B = 0.0, pi_AB = 1.0
    rho = {"R": 1.0, "A": 0.0, "B": 0.0}
    pi = {"R->A": 1.0, "A->B": 1.0}

    signals = simulate_dynamic_graph_transport(
        nodes=nodes,
        edges=edges,
        time_grid=times,
        dynamic_kernel=dek,
        local_inputs=local_inputs,
        local_recharge_fractions=rho,
        edge_mixing_fractions=pi,
    )

    # After warm-up (t >= 6), constant 10.0 should be preserved everywhere
    np.testing.assert_allclose(signals["R"], 10.0)
    np.testing.assert_allclose(signals["A"][4:], 10.0)
    np.testing.assert_allclose(signals["B"][7:], 10.0)


def test_simplex_violation_rejection():
    T = 10
    times = np.arange(T, dtype=float)
    lags = np.array([0.0, 1.0])
    nodes = ("R", "A")
    edges = (("R", "A"),)
    k_vals = np.full((1, T, 2), 0.5)
    dek = DynamicEdgeKernel(edge_ids=("R->A",), time_grid=times, lag_grid=lags, values=k_vals)

    # Simplex sum = rho_A (0.3) + pi_RA (0.8) = 1.1 != 1.0
    with pytest.raises(ValueError, match="Simplex conservation violation"):
        simulate_dynamic_graph_transport(
            nodes=nodes,
            edges=edges,
            time_grid=times,
            dynamic_kernel=dek,
            local_inputs={"R": np.ones(T)},
            local_recharge_fractions={"R": 1.0, "A": 0.3},
            edge_mixing_fractions={"R->A": 0.8},
            check_simplex=True,
        )
