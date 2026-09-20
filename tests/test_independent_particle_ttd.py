"""Regression tests for the independent particle graph-TTD generator family."""

from __future__ import annotations

import ast
from dataclasses import fields, replace
from pathlib import Path

import numpy as np
import pytest

import hydrosheaf.benchmarks.independent_particle_ttd as particle
from hydrosheaf.benchmarks.independent_particle_ttd import (
    ParticleTTDConfig,
    assert_independent_particle_generator,
    assert_particle_observation_view_is_truth_blind,
    generate_independent_particle_ttd_case,
    independent_particle_generator_provenance,
    particle_ttd_public_manifest,
    recover_particle_ttd_baseline,
    run_independent_particle_ttd_stress_suite,
    score_particle_ttd_recovery,
)


def _config(scenario: str = "nominal", *, seed: int = 417) -> ParticleTTDConfig:
    return ParticleTTDConfig(
        seed=seed,
        scenario=scenario,
        n_particles=700,
        n_steps=160,
        max_lag_steps=28,
        minimum_calibration_points=12,
    )


def test_generator_has_no_hydrosheaf_import_and_declares_independence() -> None:
    """The second family cannot silently start sharing inverse implementation code."""

    assert_independent_particle_generator()
    provenance = independent_particle_generator_provenance()
    assert provenance["imports_hydrosheaf"] is False
    assert provenance["independent_from_hydrosheaf_inference"] is True
    assert provenance["forward_representation"] == "particle_paths_then_causal_histogram"

    tree = ast.parse(Path(particle.__file__).read_text(encoding="utf-8"))
    imported_modules = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.append(node.module)
    assert not [name for name in imported_modules if name.startswith("hydrosheaf")]


def test_particle_case_is_deterministic_and_separates_latent_truth_from_public_view() -> None:
    first = generate_independent_particle_ttd_case(_config(seed=88))
    second = generate_independent_particle_ttd_case(_config(seed=88))

    assert first.truth.commitment == second.truth.commitment
    assert first.observations.sealed_commitment == second.observations.sealed_commitment
    assert first.observations.forcing_samples == second.observations.forcing_samples
    assert first.observations.calibration_samples == second.observations.calibration_samples
    assert first.observations.future_targets == second.observations.future_targets
    assert first.truth.true_edges == second.truth.true_edges
    for edge_id, kernel in first.truth.edge_kernels.items():
        other = second.truth.edge_kernels[edge_id]
        assert np.array_equal(kernel.particle_travel_times_days, other.particle_travel_times_days)
        assert np.array_equal(kernel.mass, other.mass)

    # The inference-facing object has no direct latent/truth field and its
    # manifest contains only counts, candidates, provenance, and commitment.
    public_field_names = {item.name.lower() for item in fields(first.observations)}
    assert not {name for name in public_field_names if name.startswith(("true_", "truth_"))}
    manifest = particle_ttd_public_manifest(first.observations)
    assert "future_targets" in manifest
    assert "node_responses" not in manifest
    assert "node_age_particles" not in manifest


def test_particle_kernels_conserve_mass_and_are_causal() -> None:
    case = generate_independent_particle_ttd_case(_config(seed=103))
    for kernel in case.truth.edge_kernels.values():
        assert np.all(kernel.particle_travel_times_days >= 0.0)
        assert np.all(kernel.lag_days >= 0.0)
        assert np.all(np.diff(kernel.lag_days) >= 0.0)
        assert np.all(kernel.mass >= 0.0)
        assert np.isclose(kernel.mass.sum(), 1.0, atol=1e-12)
        # An impulse at virtual time zero must reproduce a causal kernel;
        # no response appears at a negative lag because there is no such bin.
        impulse = np.zeros(kernel.mass.size, dtype=float)
        impulse[0] = 1.0
        assert np.allclose(np.convolve(impulse, kernel.mass, mode="full")[: impulse.size], kernel.mass)
    assert np.all(case.truth.node_age_particles["well_terminal"] >= 0.0)


def test_public_observation_view_is_blind_and_future_holdout_has_no_values() -> None:
    case = generate_independent_particle_ttd_case(_config(seed=157))
    observations = case.observations
    assert_particle_observation_view_is_truth_blind(observations)
    with pytest.raises(ValueError, match="forbidden latent field"):
        assert_particle_observation_view_is_truth_blind({"truth_kernel": [0.2, 0.8]})

    tampered = replace(
        observations,
        provenance={**dict(observations.provenance), "true_mean_lag": 999.0},
    )
    with pytest.raises(ValueError, match="forbidden latent field"):
        recover_particle_ttd_baseline(tampered)

    latest_calibration = max(row.time_index for row in observations.calibration_samples)
    assert observations.future_targets
    assert all(target.time_index > latest_calibration for target in observations.future_targets)
    assert all(not hasattr(target, "value") for target in observations.future_targets)


def test_truth_blind_baseline_forecasts_future_holdout_and_scores_post_hoc() -> None:
    case = generate_independent_particle_ttd_case(_config(seed=41))
    recovery = recover_particle_ttd_baseline(case.observations)
    assert recovery.status == "ESTIMATED"
    assert 0.0 <= float(recovery.young_water_fraction) <= 1.0
    assert float(recovery.mean_lag_days) >= 0.0
    assert np.isclose(sum(recovery.fitted_mass), 1.0, atol=1e-12)
    expected_targets = {
        target.time_index
        for target in case.observations.future_targets
        if target.tracer == case.observations.primary_tracer
    }
    assert {prediction.time_index for prediction in recovery.predictions} == expected_targets

    score = score_particle_ttd_recovery(case, recovery)
    assert score.commitment_verified is True
    assert score.n_held_out == len(expected_targets)
    assert score.held_out_rmse is not None
    # It is a deliberately simple coarse baseline, not an oracle.  These loose
    # bounds prove recovery of an interpretable signal without coding to truth.
    assert score.young_water_absolute_error is not None
    assert score.young_water_absolute_error < 0.40
    assert score.mean_lag_absolute_error_days is not None
    assert score.mean_lag_absolute_error_days < 95.0


def test_sparse_stress_case_explicitly_abstains() -> None:
    case = generate_independent_particle_ttd_case(_config("sparse", seed=251))
    recovery = recover_particle_ttd_baseline(case.observations)
    assert recovery.status == "ABSTAIN"
    assert recovery.reason == "insufficient_calibration_observations"
    assert recovery.predictions == ()
    score = score_particle_ttd_recovery(case, recovery)
    assert score.appropriate_abstention is True
    assert score.young_water_absolute_error is None


def test_predeclared_stress_suite_includes_all_required_negative_controls() -> None:
    runs = run_independent_particle_ttd_stress_suite(
        seed=703,
        n_particles=450,
        n_steps=150,
        max_lag_steps=24,
    )
    by_scenario = {run.case.observations.scenario: run for run in runs}
    assert set(by_scenario) == {
        "nominal",
        "sparse",
        "wrong_forcing",
        "local_recharge_misspecification",
        "wrong_topology",
    }
    assert by_scenario["wrong_forcing"].case.observations.provenance["visible_forcing_mode"] == "phase_shifted_attenuated"
    assert by_scenario["local_recharge_misspecification"].case.observations.provenance["local_forcing_exposed"] is False
    wrong_topology = by_scenario["wrong_topology"].case
    assert ("merge", "well_terminal") in wrong_topology.truth.true_edges
    assert ("merge", "well_terminal") not in wrong_topology.observations.candidate_edges
    assert ("well_terminal", "merge") in wrong_topology.observations.candidate_edges
    assert by_scenario["sparse"].recovery.status == "ABSTAIN"
