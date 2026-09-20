"""Focused tests for the truth-sealed dynamic graph-TTD virtual benchmark."""

from dataclasses import replace
import json

import numpy as np

from hydrosheaf.benchmarks.ttd_graph_dynamic import (
    SCENARIOS,
    DynamicTTDGraphConfig,
    evaluate_dynamic_ttd_recovery,
    generate_dynamic_ttd_case,
    recover_dynamic_ttd_baseline,
    run_dynamic_ttd_stress_suite,
)


def _config(**overrides):
    defaults = {
        "seed": 19,
        "n_steps": 156,
        "min_pairs_per_phase": 8,
        "min_identified_phases": 2,
    }
    defaults.update(overrides)
    return DynamicTTDGraphConfig(**defaults)


def test_dynamic_case_is_deterministic_and_keeps_truth_out_of_observations():
    config = _config(seed=41, n_steps=116)
    truth_a, observations_a, manifest_a = generate_dynamic_ttd_case(config)
    truth_b, observations_b, manifest_b = generate_dynamic_ttd_case(config)

    assert truth_a.truth_commitment == observations_a.sealed_case_commitment
    assert truth_a.truth_commitment == manifest_a["sealed_case_commitment"]
    assert truth_a.truth_commitment == truth_b.truth_commitment
    assert observations_a.observation_digest == observations_b.observation_digest
    assert manifest_a == manifest_b
    assert not hasattr(observations_a, "edge_kernels")
    assert not hasattr(observations_a, "node_signals")
    assert manifest_a["sealed_payload_in_manifest"] is False
    assert manifest_a["validation_scope"] == "controlled_synthetic_only"
    json.dumps(manifest_a, sort_keys=True)

    for node in ("A", "B", "C"):
        np.testing.assert_array_equal(
            observations_a.node_observations[node], observations_b.node_observations[node]
        )
        assert np.isnan(observations_a.node_observations[node]).any()
        assert (~observations_a.node_observation_masks[node]).any()


def test_dynamic_edge_kernels_are_normalized_causal_and_time_varying():
    truth, _, _ = generate_dynamic_ttd_case(_config(seed=7))
    lags = np.arange(next(iter(truth.edge_kernels.values())).shape[1], dtype=float)

    for kernel in truth.edge_kernels.values():
        np.testing.assert_allclose(kernel.sum(axis=1), 1.0, atol=1.0e-12)
        assert np.all(kernel >= 0.0)
        # The forward operator has support only at non-negative lag indices.
        assert kernel.shape[1] == len(lags)

    means = truth.edge_kernels["A->B"] @ lags
    assert np.ptp(means) > 1.0
    assert np.ptp(truth.local_recharge_fractions["B"]) > 0.05
    assert np.ptp(truth.direct_path_fraction) > 0.15


def test_truth_blind_baseline_recovers_identifiable_targets_and_forecasts_holdout():
    config = _config(seed=19)
    truth, observations, _ = generate_dynamic_ttd_case(config)
    recoveries = recover_dynamic_ttd_baseline(observations)
    by_edge = {result.edge_id: result for result in recoveries}

    assert by_edge["R->A"].status == "RECOVERED"
    assert by_edge["A->B"].status == "RECOVERED"
    assert by_edge["R->A"].forecast_n > 0
    assert by_edge["R->A"].forecast_r2 is not None

    score = evaluate_dynamic_ttd_recovery(
        truth,
        observations,
        recoveries,
        season_period_steps=config.season_period_steps,
        n_phase_bins=config.n_phase_bins,
    )
    assert score.truth_commitment_verified is True
    assert score.per_edge["R->A"]["classification"] == (
        "IDENTIFIABLE_RECOVERY_WITHIN_TOLERANCE"
    )
    assert score.summary["justified_recoveries"] >= 1
    assert score.summary["heldout_forecast_n"] > 0


def test_holdout_target_values_do_not_change_fitted_lags():
    config = _config(seed=23)
    _, observations, _ = generate_dynamic_ttd_case(config)
    baseline = {
        result.edge_id: result for result in recover_dynamic_ttd_baseline(observations)
    }["R->A"]

    altered_nodes = {
        name: np.array(values, copy=True)
        for name, values in observations.node_observations.items()
    }
    altered_nodes["A"][observations.training_end_step :] += 1000.0
    altered = replace(observations, node_observations=altered_nodes)
    altered_result = {
        result.edge_id: result for result in recover_dynamic_ttd_baseline(altered)
    }["R->A"]

    assert altered_result.phase_lag_steps == baseline.phase_lag_steps
    assert altered_result.mean_lag_steps == baseline.mean_lag_steps
    assert altered_result.forecast_rmse != baseline.forecast_rmse


def test_sparse_and_declared_stress_scenarios_exercise_abstention_and_protocols():
    sparse_config = _config(scenario="sparse_sampling", seed=51, n_steps=124)
    truth, observations, _ = generate_dynamic_ttd_case(sparse_config)
    sparse_recoveries = recover_dynamic_ttd_baseline(observations)
    sparse_score = evaluate_dynamic_ttd_recovery(
        truth,
        observations,
        sparse_recoveries,
        season_period_steps=sparse_config.season_period_steps,
        n_phase_bins=sparse_config.n_phase_bins,
    )
    assert all(result.status == "ABSTAIN" for result in sparse_recoveries)
    assert sparse_score.summary["correct_abstentions"] == len(sparse_recoveries)

    wrong_forcing_truth, wrong_forcing_obs, _ = generate_dynamic_ttd_case(
        _config(scenario="wrong_forcing", seed=52, n_steps=124)
    )
    assert not np.allclose(wrong_forcing_truth.true_forcing, wrong_forcing_obs.forcing)

    suite = run_dynamic_ttd_stress_suite(_config(seed=70, n_steps=112))
    assert set(suite) == set(SCENARIOS)
    assert ("A", "C") not in suite["wrong_topology"].observations.candidate_edges
    assert suite["unmodelled_local_recharge"].observations.local_recharge_available is False
    assert suite["unmodelled_local_recharge"].observations.local_recharge_inputs == {}
    assert "mis-specified_recharge_forcing" in suite["wrong_forcing"].manifest[
        "declared_stressors"
    ]
    assert suite["reversed_graph"].observations.candidate_edges == tuple(
        (target, source) for source, target in (("R", "A"), ("A", "B"), ("B", "C"), ("A", "C"))
    )
    assert all(
        edge not in (("R", "A"), ("A", "B"), ("B", "C"), ("A", "C"))
        for edge in suite["random_graph"].observations.candidate_edges
    )
    assert ("B", "C") not in suite["edge_removed_graph"].observations.candidate_edges
