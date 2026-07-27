from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_m8_confirmatory.py"
SPEC = importlib.util.spec_from_file_location("m8_confirmatory_for_test", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

bootstrap_median_ci = MODULE.bootstrap_median_ci
fixed_designs = MODULE.fixed_designs
information_diagnostics = MODULE.information_diagnostics
observation_sigma = MODULE.observation_sigma
whitened_log_jacobian = MODULE.whitened_log_jacobian


def test_fixed_design_contract_has_sixteen_equal_size_designs():
    designs = fixed_designs()
    assert len(designs) == 16
    assert all(len(times) == 4 for times in designs.values())
    assert "very_early" in designs
    assert "very_late" in designs


def test_log_jacobian_detects_exact_product_ridge():
    truth = {"k": 2.0, "A": 5.0}

    def forward(params):
        product = params["k"] * params["A"]
        return np.array([product, 2.0 * product, product**2])

    jacobian = whitened_log_jacobian(
        forward,
        truth,
        sigma=np.ones(3),
        names=("k", "A"),
        step=1e-5,
    )
    diagnostics = information_diagnostics(jacobian)
    assert np.allclose(jacobian[:, 0], jacobian[:, 1], rtol=1e-8, atol=1e-8)
    assert diagnostics["rank"] == 1
    assert diagnostics["sensitivity_cosine"] > 0.999999


def test_independent_surface_area_measurement_restores_rank():
    jacobian = np.array([[1.0, 1.0], [2.0, 2.0]])
    augmented = np.vstack([jacobian, [0.0, 10.0]])
    assert information_diagnostics(jacobian)["rank"] == 1
    assert information_diagnostics(augmented)["rank"] == 2


def test_observation_sigma_has_absolute_floor():
    sigma = observation_sigma(np.array([0.0, 0.1, 1.0]))
    assert np.allclose(sigma, [0.01, 0.01, 0.03])


def test_bootstrap_interval_is_deterministic():
    values = np.arange(1.0, 10.0)
    first = bootstrap_median_ci(values, seed=123, samples=200)
    second = bootstrap_median_ci(values, seed=123, samples=200)
    assert first == second
    assert first[1] <= first[0] <= first[2]
