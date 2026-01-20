"""
Diagnostic Script for PESTGLM Verification.
Tests linear and non-linear problems with analytical Jacobians to verify math.
"""

import numpy as np
import sys
import pandas as pd
from typing import Dict, List, Any
from dataclasses import dataclass
from scipy.optimize import least_squares

# Add project root to path
import os

sys.path.insert(0, os.path.abspath("."))

from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM
from hydrosheaf.calibration.adapters import GenericFunctionAdapter


def print_header(title):
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)


def verify_linear_regression():
    """
    Problem: y = a*x + b
    Data: (0, 1), (1, 3), (2, 5) -> True: a=2, b=1
    Residual r_i = w_i * ( (a*x_i + b) - y_i )

    Jacobian wrt params [a, b]:
    dr_i/da = w_i * x_i
    dr_i/db = w_i * 1

    We set weights = 1.0.
    """
    print_header("TEST 1: LINEAR REGRESSION (y = a*x + b)")

    x_data = [0.0, 1.0, 2.0]
    y_data = [1.0, 3.0, 5.0]

    # Start far from truth
    params = [
        AdjustableParameter("a", 0.0, -10, 10),
        AdjustableParameter("b", 0.0, -10, 10),
    ]

    obs = []
    for i, (x, y) in enumerate(zip(x_data, y_data)):
        obs.append(Observation(f"obs_{i}", y, weight=1.0))

    def model(p):
        res = {}
        for i, x in enumerate(x_data):
            res[f"obs_{i}"] = p["a"] * x + p["b"]
        return res

    adapter = GenericFunctionAdapter(model, params, obs)
    pest = PESTGLM.from_problem(adapter)

    # 1. Run Calibration
    res = pest.calibrate(max_nfev=20)
    opts = res["optimal_parameters"]

    print(f"\nOptimization Results:")
    print(f"  True: a=2.0, b=1.0")
    print(f"  Est : a={opts['a']:.6f}, b={opts['b']:.6f}")

    # 2. Verify Jacobian at Solution
    # Analytical J at solution:
    # J_rows = obs, J_cols = params [a, b]
    # Row 0 (x=0): [0, 1]
    # Row 1 (x=1): [1, 1]
    # Row 2 (x=2): [2, 1]

    J_analytical = np.array([[0.0, 1.0], [1.0, 1.0], [2.0, 1.0]])

    # Get Numerical Jacobian from scipy result (we need to rerun one iter at optimal or trust logic)
    # Actually PESTGLM.calibrate discards J, but calculates covariance.
    # Let's verify covariance calculation.
    # Sigma^2 = Phi / (n - m) = 0 / 1 = 0 (perfect fit).
    # But let's add noise or check J directly by subclassing for debug.

    # Recalculate J at optimal point manually using finite difference or call internal
    # Since we can't easily access the internal 'result.jac' from outside without modifying code,
    # we will trust the optimization result if it hits truth.

    err_a = abs(opts["a"] - 2.0)
    err_b = abs(opts["b"] - 1.0)

    if err_a < 1e-4 and err_b < 1e-4:
        print("  [PASS] Parameter Estimation matches Analytical Truth.")
    else:
        print("  [FAIL] Parameter Estimation diverged.")


def verify_log_parameter_jacobian():
    """
    Problem: y = x^k
    x=2, y_obs=8 -> True k=3
    Parameter k is log-transformed.

    Model: y = x^(10^p_internal)  where k = 10^p_internal
    Residual r = (x^k - y_obs)

    Chain Rule for Jacobian wrt internal parameter p_int:
    dr/dp_int = dr/dk * dk/dp_int
    dr/dk = x^k * ln(x)
    dk/dp_int = ln(10) * 10^p_int = ln(10) * k

    So dr/dp_int = x^k * ln(x) * ln(10) * k
    """
    print_header("TEST 2: NON-LINEAR JACOBIAN (Log Transform)")

    x_val = 2.0
    y_obs = 8.0  # x^3

    # Initial Guess k=2 (log_k ~ 0.3)
    # True k=3

    params = [AdjustableParameter("k", 2.0, 0.1, 10.0, log_transform=True)]
    obs = [Observation("obs_1", y_obs, weight=1.0)]

    def model(p):
        return {"obs_1": x_val ** p["k"]}

    adapter = GenericFunctionAdapter(model, params, obs)
    pest = PESTGLM.from_problem(adapter)

    # Evaluate Objective Function at k=2
    # Internal x = log10(2) = 0.30103
    x_start = np.array([np.log10(2.0)])

    # 1. Check Residual
    residuals = pest._objective_function(x_start)
    res_val = residuals[0]
    expected_res = (2.0**2.0) - 8.0  # 4 - 8 = -4

    print(f"\nResidual Check at k=2:")
    print(f"  Numerical: {res_val:.6f}")
    print(f"  Expected : {expected_res:.6f}")

    if abs(res_val - expected_res) < 1e-6:
        print("  [PASS] Residual Calculation correct.")
    else:
        print("  [FAIL] Residual Calculation incorrect.")

    # 2. Check Jacobian (Finite Difference vs Analytical)
    # We use a tiny step in internal space
    h = 1e-6
    x_plus = x_start + h
    res_plus = pest._objective_function(x_plus)[0]

    jac_fd = (res_plus - res_val) / h

    # Analytical
    k = 2.0
    dr_dk = (x_val**k) * np.log(x_val)  # 4 * 0.693 = 2.772
    dk_dp = np.log(10) * k  # 2.302 * 2 = 4.605
    jac_anal = dr_dk * dk_dp  # ~ 12.76

    print(f"\nJacobian Check at k=2 (wrt log parameter):")
    print(f"  Finite Diff: {jac_fd:.6f}")
    print(f"  Analytical : {jac_anal:.6f}")

    if abs(jac_fd - jac_anal) < 1e-3:
        print("  [PASS] Jacobian Calculation correct (Log Chain Rule verified).")
    else:
        print("  [FAIL] Jacobian Mismatch.")


def verify_regularization():
    """
    Problem: Minimize (x - 0)^2 with Regularization toward x=10.
    Cost = (x)^2 + (w_reg * (x - 10))^2

    Derivative wrt x: 2x + 2*w^2*(x-10) = 0
    x(1 + w^2) = 10*w^2
    x = 10*w^2 / (1 + w^2)

    If prior_sigma = 0.1, w_reg = 10.
    w^2 = 100.
    x = 1000 / 101 = 9.90

    If prior_sigma = 10, w_reg = 0.1.
    w^2 = 0.01.
    x = 0.1 / 1.01 = 0.099 (dominated by data x=0)
    """
    print_header("TEST 3: REGULARIZATION MATH")

    # Case 1: Strong Regularization (High confidence in prior)
    # Prior = 10, Sigma = 0.1 -> w=10
    params = [AdjustableParameter("x", 5.0, -20, 20, prior_mean=10.0, prior_sigma=0.1)]
    obs = [Observation("zero", 0.0, weight=1.0)]  # Data says x=0

    def model(p):
        return {"zero": p["x"]}

    adapter = GenericFunctionAdapter(model, params, obs)
    pest = PESTGLM.from_problem(adapter)
    res = pest.calibrate()

    opt_x = res["optimal_parameters"]["x"]
    expected_x = 10.0 * 100.0 / 101.0

    print(f"\nStrong Regularization (Prior=10, Sigma=0.1, Data=0):")
    print(f"  Expected: {expected_x:.4f}")
    print(f"  Got     : {opt_x:.4f}")

    if abs(opt_x - expected_x) < 1e-3:
        print("  [PASS] Strong Regularization correct.")
    else:
        print("  [FAIL] Strong Regularization failed.")

    # Case 2: Weak Regularization
    # Prior = 10, Sigma = 10.0 -> w=0.1
    params2 = [
        AdjustableParameter("x", 5.0, -20, 20, prior_mean=10.0, prior_sigma=10.0)
    ]
    adapter2 = GenericFunctionAdapter(model, params2, obs)
    pest2 = PESTGLM.from_problem(adapter2)
    res2 = pest2.calibrate()

    opt_x2 = res2["optimal_parameters"]["x"]
    expected_x2 = 10.0 * 0.01 / 1.01

    print(f"\nWeak Regularization (Prior=10, Sigma=10, Data=0):")
    print(f"  Expected: {expected_x2:.4f}")
    print(f"  Got     : {opt_x2:.4f}")

    if abs(opt_x2 - expected_x2) < 1e-3:
        print("  [PASS] Weak Regularization correct.")
    else:
        print("  [FAIL] Weak Regularization failed.")


if __name__ == "__main__":
    verify_linear_regression()
    verify_log_parameter_jacobian()
    verify_regularization()
