"""Tests for dynamic edge-kernel data structures and regularizers."""

from __future__ import annotations

import numpy as np
import pytest

from hydrosheaf.nuclear.dynamic_edge_kernel import (
    DynamicEdgeKernel,
    build_harmonic_basis_matrix,
    build_lag_curvature_matrix,
    build_phase_basis_matrix,
    build_temporal_smoothness_matrix,
)


def test_dynamic_edge_kernel_invariants_and_properties():
    edges = ("R->A", "A->B")
    time_grid = np.linspace(0.0, 100.0, 11)
    lag_grid = np.array([0.0, 7.0, 14.0, 28.0, 60.0])
    n_edges = len(edges)
    n_times = len(time_grid)
    n_lags = len(lag_grid)

    # Valid uniform kernel tensor
    values = np.full((n_edges, n_times, n_lags), 1.0 / n_lags)
    dek = DynamicEdgeKernel(
        edge_ids=edges,
        time_grid=time_grid,
        lag_grid=lag_grid,
        values=values,
        kernel_mode="time",
    )

    assert dek.n_edges == 2
    assert dek.n_output_times == 11
    assert dek.n_lags == 5
    assert dek.edge_index("A->B") == 1

    # Extract kernel at integer index and floating time
    k0 = dek.kernel_at("R->A", 0)
    assert len(k0) == 5
    np.testing.assert_allclose(np.sum(k0), 1.0, atol=1e-5)

    k_interp = dek.kernel_at("R->A", 25.0)
    assert len(k_interp) == 5
    np.testing.assert_allclose(np.sum(k_interp), 1.0, atol=1e-5)

    # Young water fraction and mean transit time
    fy = dek.young_water_fraction(cutoff_days=14.0, edge_id="R->A")
    assert len(fy) == 11
    # Lags <= 14 are 0, 7, 14 (3 out of 5 lags, each 0.2 -> 0.6)
    np.testing.assert_allclose(fy, 0.6, atol=1e-5)

    mu = dek.mean_transit_time("R->A")
    assert len(mu) == 11
    expected_mu = float(np.mean(lag_grid))
    np.testing.assert_allclose(mu, expected_mu, atol=1e-5)


def test_dynamic_edge_kernel_validation_rejections():
    edges = ("A->B",)
    times = np.array([0.0, 10.0])
    lags = np.array([0.0, 5.0])

    # 1. Unnormalized values must fail
    bad_vals = np.array([[[0.5, 0.2], [0.5, 0.2]]])
    with pytest.raises(ValueError, match="normalization violation"):
        DynamicEdgeKernel(edge_ids=edges, time_grid=times, lag_grid=lags, values=bad_vals)

    # 2. Negative values must fail
    neg_vals = np.array([[[-0.1, 1.1], [0.5, 0.5]]])
    with pytest.raises(ValueError, match="must be non-negative"):
        DynamicEdgeKernel(edge_ids=edges, time_grid=times, lag_grid=lags, values=neg_vals)

    # 3. Non-monotonic lag grid must fail
    bad_lags = np.array([10.0, 5.0])
    good_vals = np.array([[[0.5, 0.5], [0.5, 0.5]]])
    with pytest.raises(ValueError, match="strictly monotonically increasing"):
        DynamicEdgeKernel(edge_ids=edges, time_grid=times, lag_grid=bad_lags, values=good_vals)


def test_from_stationary_embedding():
    edges = ("R->A", "A->B")
    lags = (0.0, 7.0, 14.0)
    times = (0.0, 10.0, 20.0, 30.0)
    stationary = {
        "R->A": [0.6, 0.3, 0.1],
        "A->B": [0.1, 0.8, 0.1],
    }
    dek = DynamicEdgeKernel.from_stationary(edges, lags, stationary, times)
    assert dek.kernel_mode == "stationary"
    assert dek.n_output_times == 4
    for t_idx in range(4):
        np.testing.assert_allclose(dek.kernel_at("R->A", t_idx), [0.6, 0.3, 0.1])
        np.testing.assert_allclose(dek.kernel_at("A->B", t_idx), [0.1, 0.8, 0.1])


def test_regularization_matrices():
    lags = [0.0, 7.0, 14.0, 28.0]
    d2 = build_lag_curvature_matrix(lags)
    assert d2.shape == (2, 4)
    # A linear lag ramp f(tau) = a * tau + b has zero 2nd derivative
    linear_ramp = 2.0 * np.array(lags) + 3.0
    np.testing.assert_allclose(d2 @ linear_ramp, [0.0, 0.0], atol=1e-12)

    times = np.array([0.0, 10.0, 20.0, 30.0])
    dt_mat = build_temporal_smoothness_matrix(times)
    assert dt_mat.shape == (3, 4)
    constant_signal = np.full(4, 5.0)
    np.testing.assert_allclose(dt_mat @ constant_signal, 0.0, atol=1e-12)

    # Phase basis
    b_phase = build_phase_basis_matrix(times, season_period=20.0, n_phases=2)
    assert b_phase.shape == (4, 2)
    # Rows sum to 1
    np.testing.assert_allclose(np.sum(b_phase, axis=1), 1.0)

    # Harmonic basis
    b_harm = build_harmonic_basis_matrix(times, season_period=20.0, n_harmonics=1)
    assert b_harm.shape == (4, 3)  # [1, cos, sin]
