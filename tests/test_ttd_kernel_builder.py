"""Unit tests for multi-tracer forward kernel builder."""

import numpy as np
import pytest

from hydrosheaf.nuclear.ttd_grid import build_uniform_ttd_grid
from hydrosheaf.nuclear.ttd_kernel_builder import (
    TracerObservation,
    NodeTracerPanel,
    build_forward_system,
    build_network_forward_systems,
)


def test_forward_system_construction():
    grid = build_uniform_ttd_grid(max_age_years=60.0, dt_years=1.0)
    panel = NodeTracerPanel(
        node_id="Well_01",
        sample_year=2024.0,
        observations=(
            TracerObservation("3H", 5.2, 0.5, units="TU"),
            TracerObservation("SF6", 4.1, 0.4, units="pptv"),
            TracerObservation("14C", 85.0, 3.0, units="pmc"),
        ),
        head_m=125.4,
    )

    sys = build_forward_system(panel, grid)
    assert sys.node_id == "Well_01"
    assert sys.n_tracers == 3
    assert sys.matrix.shape == (3, grid.n_bins)
    assert sys.head_m == 125.4
    assert np.all(sys.observations == [5.2, 4.1, 85.0])
    assert np.all(sys.sigmas == [0.5, 0.4, 3.0])

    # Check 14C row: at tau=0, 14C response should be 100 pmc (q=1.0)
    assert sys.matrix[2, 0] == pytest.approx(100.0)

    # Condition number should be finite
    cond = sys.condition_number()
    assert np.isfinite(cond)
    assert cond > 1.0

    # Effective rank should be at least 2
    rank = sys.effective_rank()
    assert rank >= 2.0


def test_forward_prediction_and_misfit():
    grid = build_uniform_ttd_grid(max_age_years=60.0, dt_years=1.0)
    panel = NodeTracerPanel(
        node_id="Well_02",
        sample_year=2024.0,
        observations=(
            TracerObservation("14C", 50.0, 2.0, units="pmc"),
        ),
    )
    sys = build_forward_system(panel, grid)

    # A mass distribution concentrated at tau = 0 should give 100 pmc
    g_modern = np.zeros(grid.n_bins)
    g_modern[0] = 1.0
    pred = sys.predict(g_modern)
    assert pred[0] == pytest.approx(100.0)

    # Chi-squared for modern water when observed is 50 pmc: ((100 - 50)/2)^2 = 25^2 = 625
    chi2 = sys.chi_squared(g_modern)
    assert chi2 == pytest.approx(625.0)


def test_carbon14_q_correction():
    grid = build_uniform_ttd_grid(max_age_years=60.0, dt_years=1.0)
    panel = NodeTracerPanel(
        node_id="Well_03",
        sample_year=2024.0,
        observations=(TracerObservation("14C", 70.0, 2.0, units="pmc"),),
    )

    # When q = 0.85, initial 14C response at tau=0 should be 85 pmc
    sys = build_forward_system(panel, grid, q_carbon_correction=0.85)
    assert sys.matrix[0, 0] == pytest.approx(85.0)


def test_build_network_forward_systems():
    grid = build_uniform_ttd_grid(max_age_years=50.0, dt_years=2.0)
    panels = [
        NodeTracerPanel("W1", 2024.0, (TracerObservation("3H", 4.0, 0.4),)),
        NodeTracerPanel("W2", 2024.0, (TracerObservation("3H", 2.0, 0.3), TracerObservation("14C", 60.0, 2.0))),
    ]
    net_systems = build_network_forward_systems(panels, grid)
    assert set(net_systems.keys()) == {"W1", "W2"}
    assert net_systems["W1"].n_tracers == 1
    assert net_systems["W2"].n_tracers == 2
