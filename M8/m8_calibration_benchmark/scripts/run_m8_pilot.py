"""M8 pilot: does calibration success imply parameter identifiability?

The M3-M7 programme established, for *mechanisms*, that a low residual does not
imply a correct answer.  M8 asks the same question of *parameters*, which is the
one place in this thesis where a non-circular positive result is available:
synthetic truth is known exactly, so recovery can be measured directly rather
than against another model's output.

Three designs are run through the same PEST-style Gauss-Levenberg-Marquardt
engine (`hydrosheaf.calibration.glm.PESTGLM`), each at four noise levels and
many replicates:

``identifiable``
    Two parameters, each with an observation that constrains it independently.
    Full-rank Jacobian expected.

``correlated``
    Two parameters that enter the model only through their product.  Every
    observation constrains the product, none constrains either factor.  This is
    structural (not statistical) non-identifiability: the residual can be driven
    to zero along a ridge of equivalent solutions.

``overparameterised``
    More parameters than observations.  Rank-deficient by construction.

Reported per run: whether the optimiser reported success, the final objective
(phi), the true parameter error, the engine's own 95% uncertainty, whether that
interval covers the truth, and which covariance path the engine had to take.
``covariance_method`` is the key diagnostic: a fall back from ``inverse`` to
``svd_pseudoinverse`` means J^T J was singular, i.e. the engine itself detected
rank deficiency.  Nothing in standard calibration practice reports this.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hydrosheaf.calibration.adapters import GenericFunctionAdapter
from hydrosheaf.calibration.definitions import AdjustableParameter, Observation
from hydrosheaf.calibration.glm import PESTGLM

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

SEED = 20260727
NOISE_LEVELS = (0.0, 0.01, 0.03, 0.10)
N_REPLICATES = 40


# --- three designs -----------------------------------------------------------

def design_identifiable(truth):
    """y1 = a, y2 = b  -- each observation constrains one parameter."""
    a, b = truth["a"], truth["b"]
    def model(p):
        return {"y1": p["a"], "y2": p["b"], "y3": p["a"] + p["b"]}
    return model, {"y1": a, "y2": b, "y3": a + b}


def design_correlated(truth):
    """y = a*b only -- the product is identified, the factors are not."""
    a, b = truth["a"], truth["b"]
    def model(p):
        ab = p["a"] * p["b"]
        return {"y1": ab, "y2": 2.0 * ab, "y3": 0.5 * ab}
    return model, {"y1": a * b, "y2": 2 * a * b, "y3": 0.5 * a * b}


def design_overparameterised(truth):
    """Four parameters, two observations."""
    a, b, c, d = truth["a"], truth["b"], truth["c"], truth["d"]
    def model(p):
        return {"y1": p["a"] + p["b"] + p["c"] + p["d"],
                "y2": p["a"] - p["b"] + p["c"] - p["d"]}
    return model, {"y1": a + b + c + d, "y2": a - b + c - d}


DESIGNS = {
    "identifiable": (design_identifiable, {"a": 2.0, "b": 5.0}),
    "correlated": (design_correlated, {"a": 2.0, "b": 5.0}),
    "overparameterised": (design_overparameterised, {"a": 2.0, "b": 5.0, "c": 1.0, "d": 3.0}),
}


def run_case(design_name, noise, rep, rng):
    builder, truth = DESIGNS[design_name]
    model, clean_obs = builder(truth)

    noisy = {}
    for k, v in clean_obs.items():
        sigma = abs(v) * noise
        noisy[k] = v + (rng.normal(0.0, sigma) if sigma > 0 else 0.0)

    # Start away from truth so recovery is a real test, not a no-op.
    params = [
        AdjustableParameter(name=n, value=float(rng.uniform(0.2, 8.0)),
                            lower_bound=0.01, upper_bound=50.0)
        for n in truth
    ]
    obs = [Observation(name=k, value=v, weight=1.0) for k, v in noisy.items()]
    problem = GenericFunctionAdapter(model, params, obs)

    pest = PESTGLM.from_problem(problem)
    res = pest.calibrate(max_nfev=200)

    opt = res["optimal_parameters"]
    unc = res["parameter_uncertainties_95pc"]
    rows = []
    for name, true_val in truth.items():
        est = float(opt.get(name, np.nan))
        halfwidth = float(unc.get(name, np.nan))
        err = est - true_val
        rows.append({
            "design": design_name, "noise": noise, "replicate": rep,
            "parameter": name, "true_value": true_val, "estimate": est,
            "abs_error": abs(err), "rel_error": abs(err) / abs(true_val),
            "uncertainty_95pc_halfwidth": halfwidth,
            "covers_truth": bool(np.isfinite(halfwidth) and abs(err) <= halfwidth),
            "phi": float(res["phi"]),
            "optimiser_success": bool(res["success"]),
            "covariance_method": res["covariance_method"],
            "n_parameters": len(truth), "n_observations": len(obs),
        })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--replicates", type=int, default=N_REPLICATES)
    args = ap.parse_args(argv)

    rng = np.random.default_rng(SEED)
    rows = []
    for design in DESIGNS:
        for noise in NOISE_LEVELS:
            for rep in range(args.replicates):
                rows.extend(run_case(design, noise, rep, rng))
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m8_pilot_parameter_recovery.csv", index=False)

    summary = (df.groupby(["design", "noise"])
                 .agg(n=("abs_error", "size"),
                      median_rel_error=("rel_error", "median"),
                      mean_phi=("phi", "mean"),
                      optimiser_success_rate=("optimiser_success", "mean"),
                      coverage_95pc=("covers_truth", "mean"),
                      median_uncertainty=("uncertainty_95pc_halfwidth", "median"),
                      pseudoinverse_fraction=("covariance_method",
                                              lambda s: (s == "svd_pseudoinverse").mean()))
                 .reset_index())
    summary.to_csv(RESULTS / "m8_pilot_summary.csv", index=False)

    pd.set_option("display.width", 200)
    print(summary.to_string(index=False))

    (RESULTS / "m8_pilot_manifest.json").write_text(json.dumps({
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED, "replicates": args.replicates,
        "noise_levels": list(NOISE_LEVELS), "designs": list(DESIGNS),
        "engine": "hydrosheaf.calibration.glm.PESTGLM",
    }, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
