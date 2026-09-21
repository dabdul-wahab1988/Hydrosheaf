"""Unit tests for edge transport operators and node mixing algebra."""

import numpy as np
import pytest

from hydrosheaf.nuclear.ttd_grid import build_uniform_ttd_grid
from hydrosheaf.nuclear.ttd_transport import (
    EdgeTransportOperator,
    build_advection_dispersion_operator,
    build_local_recharge_distribution,
    NodeMixingSpecification,
)


def test_transport_operator_properties():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=1.0)
    op = build_advection_dispersion_operator(
        edge_id=("Well_A", "Well_B"),
        grid=grid,
        delta_tau_years=5.0,
        dispersion=0.02,
        length_m=200.0,
    )

    assert op.edge_id == ("Well_A", "Well_B")
    assert op.n_bins == grid.n_bins
    assert op.delta_tau_years == 5.0

    # 1. Exact column stochasticity: sum of each column must be 1.0
    col_sums = np.sum(op.matrix, axis=0)
    assert np.allclose(col_sums, 1.0)

    # 2. Causality: no mass can travel backwards in time (tau_arrival < tau_source)
    for k in range(grid.n_bins):
        tau_src = grid.taus[k]
        illegal_indices = np.where(grid.taus < tau_src - 1e-6)[0]
        if len(illegal_indices) > 0:
            assert np.all(op.matrix[illegal_indices, k] == 0.0)


def test_transport_propagation_shift():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=1.0)
    delta_tau = 10.0
    op = build_advection_dispersion_operator(
        edge_id=("U", "V"),
        grid=grid,
        delta_tau_years=delta_tau,
        dispersion=1e-5,  # Approximating pure piston flow
    )

    # Upstream mass at tau = 5.0 yr (bin index 5)
    g_up = np.zeros(grid.n_bins)
    g_up[5] = 1.0

    g_down = op.transport(g_up)

    # Downstream mass should peak at 5 + 10 = 15 yr (bin index 15)
    assert np.sum(g_down) == pytest.approx(1.0)
    assert np.argmax(g_down) == 15
    assert g_down[15] == pytest.approx(1.0)


def test_transport_operator_rejects_noncausal_matrix():
    grid = build_uniform_ttd_grid(max_age_years=5.0, dt_years=1.0)
    matrix = np.eye(grid.n_bins)
    matrix[0, -1] = 0.5
    matrix[-1, -1] = 0.5

    with pytest.raises(ValueError, match="causality"):
        EdgeTransportOperator(
            edge_id=("U", "V"),
            matrix=matrix,
            delta_tau_years=1.0,
        )


def test_local_recharge_distribution():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=0.5)
    r_v = build_local_recharge_distribution(grid, mean_recharge_age_years=2.0)

    assert np.sum(r_v) == pytest.approx(1.0)
    # Recharge should be young: majority of mass within 10 years
    young_mass = np.sum(r_v[grid.taus <= 10.0])
    assert young_mass > 0.95


def test_node_mixing_composite():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=1.0)

    # Two upstream nodes U1 and U2
    g_u1 = np.zeros(grid.n_bins)
    g_u1[10] = 1.0  # age 10
    g_u2 = np.zeros(grid.n_bins)
    g_u2[20] = 1.0  # age 20

    r_v = build_local_recharge_distribution(grid, mean_recharge_age_years=1.0)

    # Operators with travel times: U1 -> V (+5 yr), U2 -> V (+10 yr)
    op1 = build_advection_dispersion_operator(("U1", "V"), grid, delta_tau_years=5.0, dispersion=1e-5)
    op2 = build_advection_dispersion_operator(("U2", "V"), grid, delta_tau_years=10.0, dispersion=1e-5)

    mixing_spec = NodeMixingSpecification(
        node_id="V",
        local_fraction=0.2,
        upstream_weights={"U1": 0.5, "U2": 0.3},
        recharge_distribution=r_v,
    )

    g_v = mixing_spec.composite_distribution(
        upstream_distributions={"U1": g_u1, "U2": g_u2},
        transport_operators={("U1", "V"): op1, ("U2", "V"): op2},
    )

    # Check total probability mass conservation
    assert np.sum(g_v) == pytest.approx(1.0)

    # U1 arrival at 10 + 5 = 15 yr should carry mass 0.5
    assert g_v[15] == pytest.approx(0.5)

    # U2 arrival at 20 + 10 = 30 yr should carry mass 0.3
    assert g_v[30] == pytest.approx(0.3)
