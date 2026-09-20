"""Tests for HydroSheaf's regularized, truth-blind graph TTD inversion engine."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pytest

from hydrosheaf.benchmarks.independent_particle_ttd import (
    ParticleTTDConfig,
    generate_independent_particle_ttd_case,
    score_particle_ttd_recovery,
)
from hydrosheaf.benchmarks.ttd_graph_dynamic import (
    DynamicTTDGraphConfig,
    evaluate_dynamic_ttd_recovery,
    generate_dynamic_ttd_case,
)
from hydrosheaf.benchmarks.ttd_graph_static import (
    generate_static_ttd_graph_case,
    score_static_ttd_graph_submission,
)
from hydrosheaf.nuclear.graph_ttd_inversion import (
    GraphTTDInversionConfig,
    solve_dynamic_virtual_benchmark,
    solve_node_ttd_inversion,
    solve_particle_virtual_benchmark,
    solve_static_virtual_benchmark,
    verify_truth_blindness,
)


def test_truth_blindness_violation_rejection():
    """Verify that truth-blind gatekeeping strictly rejects sealed truth tokens."""
    # Clean dictionary must pass
    clean_payload = {"observations": [1.0, 2.0], "metadata": {"node": "A", "public_truth_policy": "sealed"}}
    verify_truth_blindness(clean_payload)

    # Dict containing forbidden key must raise ValueError
    with pytest.raises(ValueError, match="Truth-blindness violation"):
        verify_truth_blindness({"true_edges": [("R", "A")]})

    with pytest.raises(ValueError, match="Truth-blindness violation"):
        verify_truth_blindness({"edge_kernels": {"R->A": [0.5, 0.5]}})

    with pytest.raises(ValueError, match="Truth-blindness violation"):
        verify_truth_blindness({"ground_truth": 42.0})

    # Nested dataclass with truth attribute must be rejected
    @dataclass
    class LeakyPayload:
        node_id: str
        node_ttd_masses: list[float]

    leaky = LeakyPayload(node_id="M", node_ttd_masses=[0.1, 0.9])
    with pytest.raises(ValueError, match="Truth-blindness violation"):
        verify_truth_blindness(leaky)


def test_node_inversion_mass_conservation_and_recovery():
    """Verify exact mass conservation (sum to 1) and recovery of mixing weights."""
    np.random.seed(42)
    t_max = 200
    times = np.arange(t_max)
    lags = (0, 7, 14, 28, 60, 90)

    # Synthetic continuous local input and parent signals
    local_input = np.sin(2 * np.pi * times / 60.0) + 0.1 * np.random.randn(t_max)
    parent_a = np.cos(2 * np.pi * times / 90.0) + 0.1 * np.random.randn(t_max)

    # True system: 30% local recharge (lag 0), 70% from parent A with lag 14
    true_target = np.zeros(t_max, dtype=float)
    for t in range(t_max):
        loc_val = local_input[t]
        p_val = parent_a[t - 14] if t >= 14 else parent_a[0]
        true_target[t] = 0.30 * loc_val + 0.70 * p_val

    # Sample calibration observations at subset of time indices
    obs_idx = np.sort(np.random.choice(times[30:], size=40, replace=False))
    obs_vals = true_target[obs_idx] + 0.02 * np.random.randn(len(obs_idx))

    cfg = GraphTTDInversionConfig(
        lag_grid_days=lags,
        young_water_cutoff_days=30,
        regularization_smoothness=0.01,
        regularization_ridge=1e-4,
        regularization_edge_sparsity=0.01,
    )

    res = solve_node_ttd_inversion(
        node_id="test_node",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=local_input,
        parent_series_map={"parent_a": parent_a},
        config=cfg,
    )

    assert res.status == "ESTIMATED"
    assert res.r2 is not None and res.r2 > 0.85

    # Check strict mass conservation: rho + pi_A = 1.0
    total_mass = (res.local_recharge_fraction or 0.0) + sum(res.upstream_weights.values())
    assert np.isclose(total_mass, 1.0, atol=1e-5)

    # Weights must be non-negative
    assert (res.local_recharge_fraction or 0.0) >= 0.0
    assert all(w >= 0.0 for w in res.upstream_weights.values())

    # Check recovery close to true mixing fractions (30% local, 70% parent A)
    assert np.isclose(res.local_recharge_fraction or 0.0, 0.30, atol=0.15)
    assert np.isclose(res.upstream_weights["parent_a"], 0.70, atol=0.15)


def test_second_order_curvature_smoothness_effect():
    """Verify that curvature penalty D_2 reduces high-frequency lag variance."""
    np.random.seed(99)
    t_max = 150
    times = np.arange(t_max)
    lags = (0, 7, 14, 21, 28, 35, 42)

    parent_signal = np.sin(2 * np.pi * times / 45.0) + 0.05 * np.random.randn(t_max)
    # Target is convolution with smooth bell curve around lag 14
    weights_true = np.array([0.05, 0.20, 0.50, 0.20, 0.05, 0.0, 0.0])
    target = np.zeros(t_max, dtype=float)
    for t in range(t_max):
        for lag, w in zip(lags, weights_true):
            idx = t - lag
            val = parent_signal[idx] if idx >= 0 else parent_signal[0]
            target[t] += w * val

    obs_idx = np.sort(np.random.choice(times[50:], size=40, replace=False))
    obs_vals = target[obs_idx] + 0.10 * np.random.randn(len(obs_idx))

    # Unregularized fit (smoothness = 0.0)
    cfg_unsmooth = GraphTTDInversionConfig(
        lag_grid_days=lags,
        regularization_smoothness=0.0,
        regularization_ridge=1e-4,
    )
    res_unsmooth = solve_node_ttd_inversion(
        node_id="N",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=None,
        parent_series_map={"P": parent_signal},
        config=cfg_unsmooth,
    )

    # Curvature regularized fit (smoothness = 0.20)
    cfg_smooth = GraphTTDInversionConfig(
        lag_grid_days=lags,
        regularization_smoothness=0.20,
        regularization_ridge=1e-4,
    )
    res_smooth = solve_node_ttd_inversion(
        node_id="N",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=None,
        parent_series_map={"P": parent_signal},
        config=cfg_smooth,
    )

    assert res_unsmooth.status == "ESTIMATED"
    assert res_smooth.status == "ESTIMATED"

    # Compute second difference norm of the recovered kernel
    k_unsmooth = np.asarray(res_unsmooth.edge_kernels["P"])
    k_smooth = np.asarray(res_smooth.edge_kernels["P"])

    d2_unsmooth = np.diff(k_unsmooth, n=2)
    d2_smooth = np.diff(k_smooth, n=2)

    # Curvature regularization must suppress second difference norm
    assert np.linalg.norm(d2_smooth) <= np.linalg.norm(d2_unsmooth) + 1e-4


def test_false_edge_pruning():
    """Verify that spurious uncorrelated upstream parents are pruned."""
    np.random.seed(123)
    t_max = 200
    times = np.arange(t_max)
    lags = (0, 14, 28, 60)

    # True parent A
    parent_a = np.sin(2 * np.pi * times / 50.0)
    # Spurious parent B: uncorrelated high-frequency noise
    parent_b = np.random.randn(t_max)

    # Target is driven 100% by Parent A
    target = np.zeros(t_max, dtype=float)
    for t in range(t_max):
        idx = t - 14
        target[t] = parent_a[idx] if idx >= 0 else parent_a[0]

    obs_idx = np.sort(np.random.choice(times[40:], size=45, replace=False))
    obs_vals = target[obs_idx] + 0.05 * np.random.randn(len(obs_idx))

    cfg = GraphTTDInversionConfig(
        lag_grid_days=lags,
        regularization_edge_sparsity=0.03,
        edge_selection_threshold=0.05,
    )

    res = solve_node_ttd_inversion(
        node_id="downstream",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=None,
        parent_series_map={"parent_a": parent_a, "parent_b": parent_b},
        config=cfg,
    )

    assert res.status == "ESTIMATED"
    # True parent weight must dominate
    assert res.upstream_weights["parent_a"] > 0.85
    # False parent weight must be pruned below threshold
    assert res.upstream_weights["parent_b"] < cfg.edge_selection_threshold


def test_sharp_linear_programming_bounds():
    """Verify Highs LP uncertainty bounds strictly enclose point estimate and scale with residual cone."""
    np.random.seed(7)
    t_max = 120
    times = np.arange(t_max)
    lags = (0, 7, 14, 28, 60)

    p_series = np.sin(2 * np.pi * times / 40.0) + 0.05 * np.random.randn(t_max)
    target = np.roll(p_series, 7)

    obs_idx = np.sort(np.random.choice(times[30:], size=35, replace=False))
    obs_vals = target[obs_idx] + 0.04 * np.random.randn(len(obs_idx))

    cfg_tight = GraphTTDInversionConfig(
        lag_grid_days=lags,
        residual_cone_tolerance=0.03,
    )
    res_tight = solve_node_ttd_inversion(
        node_id="test",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=None,
        parent_series_map={"p": p_series},
        config=cfg_tight,
    )

    cfg_wide = GraphTTDInversionConfig(
        lag_grid_days=lags,
        residual_cone_tolerance=0.15,
    )
    res_wide = solve_node_ttd_inversion(
        node_id="test",
        target_times=obs_idx,
        target_values=obs_vals,
        local_input_series=None,
        parent_series_map={"p": p_series},
        config=cfg_wide,
    )

    assert res_tight.status == "ESTIMATED"
    assert res_tight.young_water_interval is not None
    assert res_tight.young_water_fraction is not None

    y_low, y_high = res_tight.young_water_interval
    # Bounds must be valid probabilities
    assert 0.0 <= y_low <= 1.0
    assert 0.0 <= y_high <= 1.0
    # Must strictly enclose the point estimate
    assert y_low <= res_tight.young_water_fraction <= y_high

    # Wider tolerance cone must yield equal or wider interval
    w_low, w_high = res_wide.young_water_interval  # type: ignore
    assert (w_high - w_low) >= (y_high - y_low) - 1e-6


def test_honest_calibrated_abstention():
    """Verify calibrated abstention under low sample count, constant signal, or poor fit."""
    # 1. Insufficient sample count (N = 5 < 8)
    res_low_n = solve_node_ttd_inversion(
        node_id="N_low",
        target_times=[1, 2, 3, 4, 5],
        target_values=[1.0, 1.2, 1.1, 1.3, 1.2],
        local_input_series=np.ones(10),
        parent_series_map={},
    )
    assert res_low_n.status == "ABSTAIN"
    assert res_low_n.reason == "insufficient_calibration_samples"

    # 2. No signal variation (target is flat line)
    res_flat = solve_node_ttd_inversion(
        node_id="N_flat",
        target_times=list(range(20)),
        target_values=[5.0] * 20,
        local_input_series=np.ones(30),
        parent_series_map={},
    )
    assert res_flat.status == "ABSTAIN"
    assert res_flat.reason == "no_target_variation"

    # 3. Uncorrelated noise (R^2 << min_r2)
    np.random.seed(888)
    t_len = 50
    noise_target = np.random.randn(t_len)
    signal_input = np.sin(np.linspace(0, 10, t_len))
    res_noise = solve_node_ttd_inversion(
        node_id="N_noise",
        target_times=list(range(t_len)),
        target_values=list(noise_target),
        local_input_series=signal_input,
        parent_series_map={},
        config=GraphTTDInversionConfig(lag_grid_days=(0, 1), min_r2=0.30),
    )
    assert res_noise.status == "ABSTAIN"
    assert "poor_fit_r2" in (res_noise.reason or "")


def test_solve_static_virtual_benchmark_integration():
    """Verify that solve_static_virtual_benchmark executes truth-blindly and scores cleanly."""
    case = generate_static_ttd_graph_case(seed=17)
    submission = solve_static_virtual_benchmark(case.observations)

    assert submission.method_id == "hydrosheaf_regularized_graph_ttd_inversion_v1"
    assert len(submission.held_out_predictions) > 0
    assert any(
        interval.status == "ESTIMATED"
        for interval in submission.node_intervals.values()
    )

    # Score with the sealed benchmark scorer
    score = score_static_ttd_graph_submission(case.truth, case.observations, submission)
    metrics = score.to_dict()["metrics"]

    # Provenance check: submission did NOT access truth
    assert metrics["provenance"]["truth_used_by_submission"] is False
    # Check that topology and held-out prediction metrics are populated
    assert "f1" in metrics["topology"]["candidate_graph"]
    assert "mae" in metrics["held_out_prediction"]
    assert "conditional_coverage" in metrics["interval_recovery"]


def test_solve_particle_virtual_benchmark_integration():
    """Verify that solve_particle_virtual_benchmark executes on particle tracking cases."""
    # 1. Nominal reference case: should estimate and verify commitment
    ref_case = generate_independent_particle_ttd_case(ParticleTTDConfig(scenario="nominal", seed=17))
    ref_recovery = solve_particle_virtual_benchmark(ref_case.observations)

    assert ref_recovery.status == "ESTIMATED"
    assert ref_recovery.sealed_commitment == ref_case.observations.sealed_commitment
    assert len(ref_recovery.predictions) > 0

    ref_score = score_particle_ttd_recovery(ref_case, ref_recovery)
    assert ref_score.commitment_verified is True
    assert ref_score.held_out_rmse is not None
    assert ref_score.young_water_absolute_error is not None

    # 2. Sparse calibration case: should appropriately abstain
    sparse_case = generate_independent_particle_ttd_case(ParticleTTDConfig(scenario="sparse", seed=17))
    sparse_recovery = solve_particle_virtual_benchmark(sparse_case.observations)

    assert sparse_recovery.status == "ABSTAIN"
    sparse_score = score_particle_ttd_recovery(sparse_case, sparse_recovery)
    assert sparse_score.appropriate_abstention is True


def test_solve_dynamic_virtual_benchmark_integration():
    """Verify solve_dynamic_virtual_benchmark across nominal and stress scenarios."""
    # 1. Nominal scenario
    truth_nom, obs_nom, _ = generate_dynamic_ttd_case(
        DynamicTTDGraphConfig(scenario="nominal", seed=20260917)
    )
    verify_truth_blindness(obs_nom)

    recs_nom = solve_dynamic_virtual_benchmark(obs_nom)
    assert len(recs_nom) == len(obs_nom.candidate_edges)

    eval_nom = evaluate_dynamic_ttd_recovery(truth_nom, obs_nom, recs_nom)
    assert eval_nom.truth_commitment_verified is True
    s_nom = eval_nom.summary
    assert s_nom["justified_recoveries"] == 2  # R->A and A->B
    assert s_nom["correct_abstentions"] == 2  # B->C and A->C multi-path mixing
    assert s_nom["unsupported_point_estimates"] == 0
    assert s_nom["identifiable_recovery_misses"] == 0
    assert s_nom["false_abstentions"] == 0
    assert s_nom["mean_heldout_forecast_r2"] is not None
    assert s_nom["mean_heldout_forecast_r2"] > 0.98

    # 2. Sparse sampling scenario
    truth_sp, obs_sp, _ = generate_dynamic_ttd_case(
        DynamicTTDGraphConfig(scenario="sparse_sampling", seed=20260917)
    )
    recs_sp = solve_dynamic_virtual_benchmark(obs_sp)
    eval_sp = evaluate_dynamic_ttd_recovery(truth_sp, obs_sp, recs_sp)
    assert eval_sp.summary["correct_abstentions"] == 4
    assert eval_sp.summary["unsupported_point_estimates"] == 0

    # 3. Wrong forcing scenario
    truth_wf, obs_wf, _ = generate_dynamic_ttd_case(
        DynamicTTDGraphConfig(scenario="wrong_forcing", seed=20260917)
    )
    recs_wf = solve_dynamic_virtual_benchmark(obs_wf)
    eval_wf = evaluate_dynamic_ttd_recovery(truth_wf, obs_wf, recs_wf)
    assert eval_wf.summary["justified_recoveries"] == 1  # A->B
    assert eval_wf.summary["correct_abstentions"] == 3  # R->A, B->C, A->C
    assert eval_wf.summary["unsupported_point_estimates"] == 0

    # 4. Reversed graph topology scenario
    truth_rev, obs_rev, _ = generate_dynamic_ttd_case(
        DynamicTTDGraphConfig(scenario="reversed_graph", seed=20260917)
    )
    recs_rev = solve_dynamic_virtual_benchmark(obs_rev)
    eval_rev = evaluate_dynamic_ttd_recovery(truth_rev, obs_rev, recs_rev)
    assert eval_rev.summary["correct_abstentions"] == 4
    assert eval_rev.summary["unsupported_point_estimates"] == 0
