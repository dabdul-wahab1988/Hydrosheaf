"""
Performance Test for Parallel PESTGLM.
"""

import time
import numpy as np
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM
from hydrosheaf.calibration.adapters import GenericFunctionAdapter


def slow_model(params):
    """A model that simulates a slow process (0.1s delay)."""
    import time

    time.sleep(0.1)  # Simulate workload
    # We used "x_0" etc in the params definition
    return {"obs": params["x_0"] * 2.0}  # Just use one param for dummy output


def benchmark_parallel():
    print("\n" + "=" * 80)
    print("PARALLEL BENCHMARK (Simulated Workload)")
    print("=" * 80)

    # Setup Problem with 10 parameters
    # Total Jacobian cost = 10 runs + 1 base run = 11 runs * 0.1s = 1.1s minimum in serial.
    params = [AdjustableParameter(f"x_{i}", 1.0) for i in range(10)]
    obs = [Observation("obs", 2.0)]

    # We need a model that uses all params to be realistic, but for timing
    # the structure doesn't matter, only the sleep.

    # Serial Run
    print("\nRunning SERIAL (1 worker)...")
    adapter = GenericFunctionAdapter(slow_model, params, obs)
    pest_serial = PESTGLM.from_problem(adapter, n_workers=1)

    t0 = time.time()
    # Force one Jacobian calculation
    x0 = pest_serial.x0
    J_serial = pest_serial._calculate_jacobian(x0)
    t_serial = time.time() - t0

    print(f"Serial Time: {t_serial:.4f} s")

    # Parallel Run (4 workers)
    print("\nRunning PARALLEL (4 workers)...")
    pest_parallel = PESTGLM.from_problem(adapter, n_workers=4, worker_type="thread")

    t0 = time.time()
    J_parallel = pest_parallel._calculate_jacobian(x0)
    t_parallel = time.time() - t0

    print(f"Parallel Time: {t_parallel:.4f} s")

    speedup = t_serial / t_parallel
    print(f"\nSpeedup: {speedup:.2f}x")

    if speedup > 1.5:
        print("[PASS] Parallelism is working effectively.")
    else:
        print("[WARN] Speedup low (overhead might dominate for this small test).")


if __name__ == "__main__":
    benchmark_parallel()
