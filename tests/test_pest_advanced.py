"""
Advanced Verification Suite for PESTGLM.
Tests parallel consistency, uncertainty quantification accuracy, and boundary handling.
"""

import pytest
import numpy as np

from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM
from hydrosheaf.calibration.adapters import GenericFunctionAdapter

# --- Helpers ---


def get_linear_problem(n_params=5, noise_level=0.0):
    """
    Creates a linear problem y = sum(p_i * x)
    """
    # True params: all 1.0
    {f"p_{i}": 1.0 for i in range(n_params)}

    # Obs: y = sum(p) for different coefficients?
    # Simple setup: Obs i depends only on p_i (diagonal Jacobian)
    # y_i = p_i * 1.0

    params = [AdjustableParameter(f"p_{i}", 0.5, -10, 10) for i in range(n_params)]
    obs = []

    for i in range(n_params):
        val = 1.0
        if noise_level > 0:
            val += np.random.normal(0, noise_level)
        obs.append(Observation(f"obs_{i}", val))

    def model(p):
        res = {}
        for i in range(n_params):
            res[f"obs_{i}"] = p[f"p_{i}"]
        return res

    return GenericFunctionAdapter(model, params, obs)


def get_rosenbrock_problem():
    """
    Standard optimization benchmark.
    f(x,y) = (a-x)^2 + b(y-x^2)^2
    Global min at (a, a^2). Usually a=1, b=100.
    Min at (1, 1).

    We formulate this as least squares:
    r1 = a - x
    r2 = sqrt(b) * (y - x^2)
    Targets = 0.
    """
    params = [
        AdjustableParameter("x", 0.0, -5, 5),
        AdjustableParameter("y", 0.0, -5, 5),
    ]

    obs = [Observation("r1", 0.0), Observation("r2", 0.0)]  # Matches r1  # Matches r2

    def model(p):
        x = p["x"]
        y = p["y"]
        # r1 = 1 - x
        # r2 = 10 * (y - x^2)
        return {"r1": 1.0 - x, "r2": 10.0 * (y - x**2)}

    return GenericFunctionAdapter(model, params, obs)


# --- Tests ---


def test_parallel_consistency():
    """
    Verify that Serial and Parallel modes produce IDENTICAL results
    for a non-trivial problem (Rosenbrock).
    """
    problem = get_rosenbrock_problem()

    # Serial Run
    # Seed randomness just in case (though least_squares is deterministic)
    pest_serial = PESTGLM.from_problem(problem, n_workers=1)
    res_serial = pest_serial.calibrate(max_nfev=200)

    # Parallel Run
    pest_parallel = PESTGLM.from_problem(problem, n_workers=2, worker_type="thread")
    res_parallel = pest_parallel.calibrate(max_nfev=200)

    # 1. Check Success
    assert res_serial["success"], "Serial optimization failed"
    assert res_parallel["success"], "Parallel optimization failed"

    # 2. Check Optimal Values (Approaching 1.0, 1.0)
    opt_s = res_serial["optimal_parameters"]
    opt_p = res_parallel["optimal_parameters"]

    print(f"Serial: {opt_s}")
    print(f"Parallel: {opt_p}")

    assert np.isclose(opt_s["x"], 1.0, atol=1e-4)
    assert np.isclose(opt_s["y"], 1.0, atol=1e-4)

    # 3. Check Consistency between modes
    # Allow small tolerance because finite difference step sizes might vary slightly
    # between scipy internal Jacobian and our custom parallel Jacobian.
    assert np.allclose(opt_s["x"], opt_p["x"], atol=1e-3), "Serial/Parallel x mismatch"
    assert np.allclose(opt_s["y"], opt_p["y"], atol=1e-3), "Serial/Parallel y mismatch"


def test_uncertainty_quantification():
    """
    Verify that reported uncertainty matches statistical theory.
    Problem: y = p * x.
    We have 100 observations of p=2.0 with noise sigma=0.1.
    Standard Error of Mean (SEM) for p should be sigma / sqrt(n).
    SEM = 0.1 / 10 = 0.01.
    95% Conf Interval approx +/- 1.96 * SEM = +/- 0.0196.
    """
    n_obs = 100
    true_p = 2.0
    noise_sigma = 0.1

    np.random.seed(42)

    params = [AdjustableParameter("p", 1.5, 0, 10)]
    obs = []

    for i in range(n_obs):
        # x is 1.0 for simplicity, so y = p + noise
        val = true_p + np.random.normal(0, noise_sigma)
        obs.append(Observation(f"o_{i}", val))

    def model(p):
        return {f"o_{i}": p["p"] for i in range(n_obs)}

    pest = PESTGLM.from_problem(GenericFunctionAdapter(model, params, obs), n_workers=1)
    res = pest.calibrate()

    est_p = res["optimal_parameters"]["p"]
    unc_p = res["parameter_uncertainties_95pc"]["p"]  # This is 1.96 * sigma_p
    sigma_p = unc_p / 1.96

    print(f"Estimated p: {est_p:.4f}")
    print(f"Estimated sigma_p: {sigma_p:.4f}")
    print(f"Theoretical sigma: {noise_sigma/np.sqrt(n_obs):.4f}")

    # Check Parameter Accuracy
    assert np.isclose(est_p, true_p, atol=2 * 0.01)  # Within 2 std errors

    # Check Uncertainty Accuracy
    # Note: PESTGLM calculates sigma^2 based on *actual* residuals (Phi).
    # Since we added noise 0.1, the residuals should reflect that.
    theoretical_sigma = noise_sigma / np.sqrt(n_obs)  # 0.01

    # Allow 20% error in variance estimation (finite sample)
    assert np.isclose(sigma_p, theoretical_sigma, rtol=0.20)


def test_log_bounds_handling():
    """
    Test that log-transformed parameters strictly respect bounds.
    """
    # Problem: Minimize (x - 50)^2
    # Bounds: [100, 200].
    # Optimal solution should be clamped at 100 (closest to 50).
    # Log Bounds: [2, 2.3]. Target 50 is log 1.7 (out of bounds).

    params = [AdjustableParameter("x", 150.0, 100.0, 200.0, log_transform=True)]
    obs = [Observation("obs", 50.0)]  # Target 50

    def model(p):
        return {"obs": p["x"]}

    pest = PESTGLM.from_problem(GenericFunctionAdapter(model, params, obs))
    res = pest.calibrate()

    opt_x = res["optimal_parameters"]["x"]

    print(f"Target: 50. Bounds: [100, 200]. Result: {opt_x}")

    # Should be effectively 100.0
    assert np.isclose(opt_x, 100.0, atol=1e-3)
    assert opt_x >= 100.0 - 1e-6  # strictly respects lower bound


def test_parallel_failure_recovery():
    """
    Test how parallel mode handles a model run that returns garbage/missing data.
    """
    params = [AdjustableParameter("x", 1.0), AdjustableParameter("y", 1.0)]
    obs = [Observation("o1", 1.0), Observation("o2", 1.0)]

    def buggy_model(p):
        # Simulate a crash for certain parameters
        if p["x"] > 1.000001:
            # When computing Jacobian, x will be perturbed to 1.000001
            # Let's return empty dict to simulate crash
            return {}
        return {"o1": p["x"], "o2": p["y"]}

    pest = PESTGLM.from_problem(
        GenericFunctionAdapter(buggy_model, params, obs),
        n_workers=2,
        worker_type="thread",
    )

    # The GenericFunctionAdapter doesn't handle failures gracefully by default (returns KeyError if key missing).
    # But PESTGLM._get_residuals uses .get(name, -9999.0).
    # So if model returns {}, residuals will be huge (-9999 - value).
    # This should drive the optimizer away from that region or fail, but not crash code.

    try:
        # We perform one function eval to calculate Jacobian
        # Base x=[1,1]. Perturbed x=[1+eps, 1].
        # Buggy model fails on perturbation.
        # Residuals for perturbation become huge.
        # Jacobian derivative becomes huge.
        pest.calibrate(max_nfev=2)
    except Exception as e:
        pytest.fail(f"Optimizer crashed on model failure: {e}")

    # We expect it to run, even if results are bad
    assert True
