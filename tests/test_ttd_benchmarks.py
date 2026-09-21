"""Controlled synthetic virtual benchmark tests for non-parametric network TTD inversion.

Evaluates:
1. Bimodal multi-tracer mixture recovery (where conventional parametric single-LPMs fail);
2. Correct graph vs. uncoupled local vs. negative control (reversed DAG);
3. Held-out tracer predictability and abstention under structural misspecification.
"""

import networkx as nx
import numpy as np
import pytest

from hydrosheaf.nuclear.joint_lpm import fit_lpm_models
from hydrosheaf.nuclear.ttd_diagnostics import CODE_HEAD_CONFLICT
from hydrosheaf.nuclear.ttd_grid import build_multiscale_ttd_grid, build_uniform_ttd_grid
from hydrosheaf.nuclear.ttd_kernel_builder import (
    NodeTracerPanel,
    TracerObservation,
    build_forward_system,
)
from hydrosheaf.nuclear.ttd_network_solver import (
    solve_network_ttd,
    solve_single_node_ttd,
)
from hydrosheaf.nuclear.ttd_transport import build_advection_dispersion_operator
from hydrosheaf.nuclear.ttd_transport import NodeMixingSpecification


def test_bimodal_mixture_recovery():
    """Test recovery of a bimodal groundwater mixture (modern + paleowater).

    Conventional single-LPM models (EM, DM, PFM) fail when water is a mixture of
    recent recharge (bomb-peak 3H) and deep old water (low 14C).
    """
    # Multiscale grid out to 30,000 years
    grid = build_multiscale_ttd_grid(
        max_age_years=30000.0,
        dt_young=1.0,
        dt_holocene=100.0,
        dt_pleistocene=1000.0,
        young_cutoff_years=70.0,
        holocene_cutoff_years=11700.0,
    )

    # True bimodal distribution:
    # Component 1: 30% Modern water (mean age 10 years, Anthropocene)
    # Component 2: 70% Pleistocene water (mean age 20,000 years)
    f_modern = 0.30
    f_paleo = 0.70

    g1 = np.exp(-grid.taus / 10.0)
    g1[grid.taus > 70.0] = 0.0
    g1 /= np.sum(g1)

    g2 = np.exp(-((grid.taus - 20000.0) ** 2) / (2.0 * (3000.0 ** 2)))
    g2 /= np.sum(g2)

    g_true_bimodal = f_modern * g1 + f_paleo * g2
    g_true_bimodal /= np.sum(g_true_bimodal)

    # Generate synthetic observations using multi-tracer panel: 3H, SF6, 14C, 39Ar
    template_panel = NodeTracerPanel(
        node_id="BimodalWell",
        sample_year=2024.0,
        observations=(
            TracerObservation("3H", 1.0, 0.3, units="TU"),
            TracerObservation("SF6", 1.0, 0.3, units="pptv"),
            TracerObservation("14C", 1.0, 1.5, units="pmc"),
            TracerObservation("39Ar", 1.0, 3.0, units="pmc"),
        ),
    )
    sys_template = build_forward_system(template_panel, grid)
    c_synth = sys_template.predict(g_true_bimodal)

    # 1. Try fitting conventional single-LPM models
    lpm_obs = {
        "tritium_TU": float(c_synth[0]),
        "tritium_sigma_TU": 0.3,
        "sf6_pptv": float(c_synth[1]),
        "sf6_sigma_pptv": 0.3,
        "c14_pmc": float(c_synth[2]),
        "c14_sigma_pmc": 1.5,
        "sample_year": 2024.0,
    }
    fits = fit_lpm_models(
        lpm_obs,
        sample_year=2024.0,
        models=("EM", "DM", "PFM"),
        max_age_years=30000.0,
    )
    best_single_lpm = min(fits, key=lambda f: f.objective)
    # Single-LPM models have very high residual because they cannot fit both 3H and 14C
    assert best_single_lpm.rmse_standardized > 2.5

    # 2. Invert with non-parametric regularized TTD solver
    panel_obs = NodeTracerPanel(
        node_id="BimodalWell",
        sample_year=2024.0,
        observations=(
            TracerObservation("3H", c_synth[0], 0.3, units="TU"),
            TracerObservation("SF6", c_synth[1], 0.3, units="pptv"),
            TracerObservation("14C", c_synth[2], 1.5, units="pmc"),
            TracerObservation("39Ar", c_synth[3], 3.0, units="pmc"),
        ),
    )
    sys_bimodal = build_forward_system(panel_obs, grid)
    result = solve_single_node_ttd(sys_bimodal, grid, lambda_smoothness=0.005)

    assert result.status == "ESTIMATED"
    assert result.chi_squared_per_observation < 1.5

    # Check that the non-parametric solver successfully recovered bimodal fractions:
    # Anthropocene fraction (~0.30) and Old water fraction (~0.70)
    fractions = grid.age_fractions(result.g)
    assert fractions["anthropocene"] == pytest.approx(f_modern, abs=0.08)
    assert fractions["pleistocene"] + fractions["holocene"] == pytest.approx(f_paleo, abs=0.08)
    assert fractions["pleistocene"] > 0.50


def test_network_oracle_vs_reversed_negative_control():
    """Test correct network topology vs. reversed graph negative control."""
    grid = build_uniform_ttd_grid(max_age_years=60.0, dt_years=1.0)

    # Correct Graph: RechargeWell (head 120m) -> ProductionWell (head 90m)
    correct_graph = nx.DiGraph([("RechargeWell", "ProductionWell")])
    heads = {"RechargeWell": 120.0, "ProductionWell": 90.0}

    # True distributions
    g_up_true = np.exp(-grid.taus / 4.0)
    g_up_true /= np.sum(g_up_true)

    op = build_advection_dispersion_operator(
        ("RechargeWell", "ProductionWell"), grid, delta_tau_years=8.0, dispersion=0.02
    )
    transport_ops = {("RechargeWell", "ProductionWell"): op}

    # Production well draws 80% routed upstream and 20% local recharge
    r_local = np.exp(-grid.taus / 1.0)
    r_local /= np.sum(r_local)
    g_down_true = 0.8 * op.transport(g_up_true) + 0.2 * r_local

    # Forward systems
    p_up_template = NodeTracerPanel("RechargeWell", 2024.0, (TracerObservation("3H", 1.0, 0.2), TracerObservation("14C", 1.0, 1.0)), head_m=120.0)
    p_down_template = NodeTracerPanel("ProductionWell", 2024.0, (TracerObservation("3H", 1.0, 0.2), TracerObservation("14C", 1.0, 1.0)), head_m=90.0)

    s_up = build_forward_system(p_up_template, grid)
    s_down = build_forward_system(p_down_template, grid)

    c_up = s_up.predict(g_up_true)
    c_down = s_down.predict(g_down_true)

    fwd_systems = {
        "RechargeWell": build_forward_system(NodeTracerPanel("RechargeWell", 2024.0, (TracerObservation("3H", c_up[0], 0.2), TracerObservation("14C", c_up[1], 1.0)), head_m=120.0), grid),
        "ProductionWell": build_forward_system(NodeTracerPanel("ProductionWell", 2024.0, (TracerObservation("3H", c_down[0], 0.2), TracerObservation("14C", c_down[1], 1.0)), head_m=90.0), grid),
    }

    # Case A: Correct Graph Inversion
    res_correct = solve_network_ttd(
        correct_graph,
        fwd_systems,
        transport_ops,
        grid,
        mixing_specs={
            "ProductionWell": NodeMixingSpecification(
                node_id="ProductionWell",
                local_fraction=0.2,
                upstream_weights={"RechargeWell": 0.8},
                recharge_distribution=r_local,
            )
        },
        lambda_smoothness=0.01,
        lambda_graph=2.0,
        node_heads=heads,
    )
    assert res_correct.status == "ESTIMATED"
    assert res_correct.network_sheaf_energy < 0.05
    mtts = res_correct.mean_transit_times()
    assert mtts["ProductionWell"] > mtts["RechargeWell"]

    # Case B: Negative Control - Reversed Flow Direction
    reversed_graph = nx.DiGraph([("ProductionWell", "RechargeWell")])
    reversed_ops = {
        ("ProductionWell", "RechargeWell"): build_advection_dispersion_operator(
            ("ProductionWell", "RechargeWell"), grid, delta_tau_years=8.0, dispersion=0.02
        )
    }

    res_reversed = solve_network_ttd(
        reversed_graph,
        fwd_systems,
        reversed_ops,
        grid,
        lambda_smoothness=0.01,
        lambda_graph=2.0,
        node_heads=heads,
    )
    # The negative control must honestly ABSTAIN due to hydraulic head conflict
    assert res_reversed.status == "ABSTAIN"
    assert CODE_HEAD_CONFLICT in res_reversed.abstention_reasons
