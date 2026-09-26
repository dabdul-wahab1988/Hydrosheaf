"""Unit tests for multi-scale age discretization grid and operators."""

import numpy as np
import pytest

from hydrosheaf.nuclear.ttd_grid import (
    build_multiscale_ttd_grid,
    build_uniform_ttd_grid,
    build_d1_difference_matrix,
    build_d2_curvature_matrix,
    shannon_entropy,
    wasserstein_1d,
)


def test_uniform_grid_construction():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=1.0)
    assert grid.min_age == 0.0
    assert grid.max_age == 50.0
    assert grid.n_bins == 51
    # Check trapezoidal weights: endpoints are 0.5, interior are 1.0
    assert grid.weights[0] == pytest.approx(0.5)
    assert grid.weights[-1] == pytest.approx(0.5)
    assert np.allclose(grid.weights[1:-1], 1.0)
    # Total integral of dx on [0, 50] is 50.0
    assert np.sum(grid.weights) == pytest.approx(50.0)


def test_multiscale_grid_eras():
    grid = build_multiscale_ttd_grid(
        max_age_years=20000.0,
        dt_young=1.0,
        dt_holocene=100.0,
        dt_pleistocene=1000.0,
        young_cutoff_years=70.0,
        holocene_cutoff_years=11700.0,
    )
    assert grid.min_age == 0.0
    assert grid.max_age == 20000.0
    # Young window has 1-year resolution
    idx_young = np.where(grid.taus <= 70.0)[0]
    assert len(idx_young) > 60
    assert np.all(np.diff(grid.taus) > 0.0)
    # Total sum of weights equals interval length
    assert np.sum(grid.weights) == pytest.approx(20000.0)


def test_grid_statistical_functionals():
    grid = build_uniform_ttd_grid(max_age_years=300.0, dt_years=0.5)
    # Exponential distribution f(tau) = 1/20 * exp(-tau / 20)
    tau_m = 20.0
    pdf = (1.0 / tau_m) * np.exp(-grid.taus / tau_m)
    # Normalize density with trapezoidal weights
    pdf = pdf / np.sum(pdf * grid.weights)

    # Mean transit time on [0, 300] for exp(20) has negligible truncation (exp(-15) ~ 3e-7)
    mean_tau = grid.mean_transit_time(pdf, is_pdf=True)
    assert mean_tau == pytest.approx(tau_m, rel=0.005)

    # Anthropocene fraction (tau <= 70) for exp(20): 1 - exp(-70/20) = 1 - exp(-3.5) = 0.9698
    fractions = grid.age_fractions(pdf, is_pdf=True)
    assert fractions["anthropocene"] == pytest.approx(1.0 - np.exp(-70.0 / 20.0), rel=0.01)
    assert fractions["pleistocene"] == 0.0  # Max age is 300 yr


def test_difference_and_curvature_operators():
    grid = build_uniform_ttd_grid(max_age_years=10.0, dt_years=1.0)
    d1 = build_d1_difference_matrix(grid)
    d2 = build_d2_curvature_matrix(grid)

    assert d1.shape == (grid.n_bins - 1, grid.n_bins)
    assert d2.shape == (grid.n_bins - 2, grid.n_bins)

    # Linear function f(tau) = 2 * tau + 3: D1 should be 2.0, D2 should be 0.0
    f_linear = 2.0 * grid.taus + 3.0
    assert np.allclose(d1 @ f_linear, 2.0)
    assert np.allclose(d2 @ f_linear, 0.0)

    # Quadratic function f(tau) = tau^2: D2 should be 2.0 everywhere
    f_quad = grid.taus ** 2
    assert np.allclose(d2 @ f_quad, 2.0)


def test_wasserstein_and_entropy():
    grid = build_uniform_ttd_grid(max_age_years=10.0, dt_years=1.0)
    # Delta at tau=2 vs Delta at tau=5
    p1 = np.zeros(grid.n_bins)
    p2 = np.zeros(grid.n_bins)
    p1[2] = 1.0  # at tau = 2
    p2[5] = 1.0  # at tau = 5

    # W1 distance should equal |5 - 2| = 3.0
    w1 = wasserstein_1d(p1, p2, grid)
    assert w1 == pytest.approx(3.0)

    # Entropy of deterministic mass is 0
    assert shannon_entropy(p1) == pytest.approx(0.0, abs=1e-6)

    # Uniform mass over 4 bins
    p_unif = np.zeros(grid.n_bins)
    p_unif[:4] = 0.25
    assert shannon_entropy(p_unif) == pytest.approx(np.log(4.0))

    # Endpoint atoms must integrate the full interval rather than half of it.
    p0 = np.zeros(grid.n_bins)
    p10 = np.zeros(grid.n_bins)
    p0[0] = 1.0
    p10[-1] = 1.0
    assert wasserstein_1d(p0, p10, grid) == pytest.approx(10.0)
