"""M8 extension: the identifiability threshold, and optimal experimental design.

Two additions to the pilot.

**Part 1 - continuous conditioning.**  The pilot contrasted three discrete designs.
A contrast tells a practitioner that identifiability matters; a *threshold* tells
them when to worry.  Here a single parameter epsilon controls how ill-conditioned
the problem is:

    y1 = a + b
    y2 = a + (1 + eps) * b
    y3 = a + (1 + 2*eps) * b

As eps -> 0 the two parameters become indistinguishable and the Fisher
information matrix approaches singularity, continuously.  Sweeping eps locates
the condition number at which parameter recovery and interval coverage collapse.

**Part 2 - optimal experimental design.**  The Fisher information matrix
``FIM = J^T W J`` is already implicit in the Jacobian PESTGLM computes, but
Hydrosheaf currently exposes no design machinery.  Its eigen-decomposition names
the unidentifiable direction, so a next-measurement can be chosen to attack it:

    D-optimality : maximise det(FIM)          - shrinks the confidence ellipsoid volume
    A-optimality : minimise trace(FIM^-1)     - shrinks average parameter variance
    E-optimality : maximise lambda_min(FIM)   - shrinks the WORST-determined direction

E-optimality is the natural criterion for this thesis because the worst-determined
direction *is* the non-identifiability.  Each criterion selects a next observation
from a candidate pool; the selection is then *verified* by actually adding the
observation, re-calibrating, and measuring whether recovery improved, against a
random-selection control and an exhaustive oracle.
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
TRUTH = {"a": 2.0, "b": 5.0}
NOISE = 0.03


# --- shared helpers ----------------------------------------------------------

def numeric_jacobian(model, params, names, eps=1e-6):
    base = model(params)
    keys = list(base)
    J = np.zeros((len(keys), len(names)))
    for j, n in enumerate(names):
        q = dict(params)
        q[n] = params[n] + eps
        pert = model(q)
        for i, k in enumerate(keys):
            J[i, j] = (pert[k] - base[k]) / eps
    return J


def fim_diagnostics(J):
    """Condition number, eigenvalues and the least-determined direction."""
    FIM = J.T @ J
    w, v = np.linalg.eigh(FIM)
    w = np.clip(w, 0.0, None)
    lam_min, lam_max = float(w.min()), float(w.max())
    cond = lam_max / lam_min if lam_min > 0 else np.inf
    with np.errstate(divide="ignore", invalid="ignore"):
        det = float(np.linalg.det(FIM))
        trace_inv = float(np.trace(np.linalg.pinv(FIM, rcond=1e-15)))
    return {
        "cond_fim": cond, "lambda_min": lam_min, "lambda_max": lam_max,
        "det_fim": det, "trace_fim_inv": trace_inv,
        "null_direction": v[:, 0].tolist(),
    }


def calibrate(model, truth, obs_defs, rng, noise=NOISE, max_nfev=200):
    clean = model(truth)
    obs_vals = {k: clean[k] + rng.normal(0.0, abs(clean[k]) * noise) for k in obs_defs}
    params = [AdjustableParameter(name=n, value=float(rng.uniform(0.2, 8.0)),
                                  lower_bound=0.01, upper_bound=50.0) for n in truth]
    obs = [Observation(name=k, value=v, weight=1.0) for k, v in obs_vals.items()]
    res = PESTGLM.from_problem(GenericFunctionAdapter(model, params, obs)).calibrate(max_nfev=max_nfev)
    opt, unc = res["optimal_parameters"], res["parameter_uncertainties_95pc"]
    errs, covers = [], []
    for n, tv in truth.items():
        e = abs(float(opt.get(n, np.nan)) - tv) / abs(tv)
        hw = float(unc.get(n, np.nan))
        errs.append(e)
        covers.append(bool(np.isfinite(hw) and abs(float(opt.get(n, np.nan)) - tv) <= hw))
    return {
        "median_rel_error": float(np.median(errs)),
        "coverage": float(np.mean(covers)),
        "phi": float(res["phi"]),
        "covariance_method": res["covariance_method"],
    }


# --- Part 1: conditioning sweep ---------------------------------------------

def make_conditioned_model(eps):
    def model(p):
        return {"y1": p["a"] + p["b"],
                "y2": p["a"] + (1.0 + eps) * p["b"],
                "y3": p["a"] + (1.0 + 2.0 * eps) * p["b"]}
    return model


def part1_conditioning(replicates, rng):
    rows = []
    for eps in np.logspace(0, -6, 13):
        model = make_conditioned_model(float(eps))
        diag = fim_diagnostics(numeric_jacobian(model, TRUTH, list(TRUTH)))
        for rep in range(replicates):
            out = calibrate(model, TRUTH, ["y1", "y2", "y3"], rng)
            rows.append({"epsilon": float(eps), "replicate": rep,
                         "cond_fim": diag["cond_fim"], "lambda_min": diag["lambda_min"],
                         **{k: v for k, v in out.items()}})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m8_conditioning_sweep.csv", index=False)
    summ = (df.groupby(["epsilon", "cond_fim"])
              .agg(n=("median_rel_error", "size"),
                   median_rel_error=("median_rel_error", "median"),
                   coverage_95pc=("coverage", "mean"),
                   mean_phi=("phi", "mean"),
                   pseudoinverse_fraction=("covariance_method",
                                           lambda s: (s == "svd_pseudoinverse").mean()))
              .reset_index().sort_values("cond_fim"))
    summ.to_csv(RESULTS / "m8_conditioning_summary.csv", index=False)
    return summ


# --- Part 2: optimal experimental design ------------------------------------

CANDIDATES = {
    "measure_a_directly":   lambda p: p["a"],
    "measure_b_directly":   lambda p: p["b"],
    "contrast_a_minus_b":   lambda p: p["a"] - p["b"],
    "ratio_like_2a_minus_b": lambda p: 2.0 * p["a"] - p["b"],
    "redundant_a_plus_b":   lambda p: p["a"] + p["b"],
}


def part2_oed(replicates, rng, eps=1e-3):
    base_model = make_conditioned_model(eps)
    base_names = ["y1", "y2", "y3"]

    # Score each candidate by the three classical design criteria.
    scores = []
    for cname, cfun in CANDIDATES.items():
        def augmented(p, cfun=cfun):
            out = base_model(p)
            out["y_new"] = cfun(p)
            return out
        d = fim_diagnostics(numeric_jacobian(augmented, TRUTH, list(TRUTH)))
        scores.append({"candidate": cname, "D_det_fim": d["det_fim"],
                       "A_trace_fim_inv": d["trace_fim_inv"],
                       "E_lambda_min": d["lambda_min"], "cond_fim": d["cond_fim"]})
    sdf = pd.DataFrame(scores)
    picks = {
        "E_optimal": sdf.loc[sdf.E_lambda_min.idxmax(), "candidate"],
        "D_optimal": sdf.loc[sdf.D_det_fim.idxmax(), "candidate"],
        "A_optimal": sdf.loc[sdf.A_trace_fim_inv.idxmin(), "candidate"],
        "worst_choice": sdf.loc[sdf.E_lambda_min.idxmin(), "candidate"],
    }

    # Verify by actually adding the observation and re-calibrating.
    strategies = dict(picks)
    strategies["random"] = None  # drawn per replicate
    rows = []
    for strat, chosen in strategies.items():
        for rep in range(replicates):
            cname = chosen if chosen is not None else rng.choice(list(CANDIDATES))
            cfun = CANDIDATES[cname]
            def augmented(p, cfun=cfun):
                out = base_model(p)
                out["y_new"] = cfun(p)
                return out
            out = calibrate(augmented, TRUTH, base_names + ["y_new"], rng)
            rows.append({"strategy": strat, "chosen_observation": cname,
                         "replicate": rep, **out})
    # No-addition baseline
    for rep in range(replicates):
        out = calibrate(base_model, TRUTH, base_names, rng)
        rows.append({"strategy": "no_new_measurement", "chosen_observation": "-",
                     "replicate": rep, **out})

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m8_oed_verification.csv", index=False)
    sdf.to_csv(RESULTS / "m8_oed_candidate_scores.csv", index=False)
    summ = (df.groupby("strategy")
              .agg(n=("median_rel_error", "size"),
                   median_rel_error=("median_rel_error", "median"),
                   coverage_95pc=("coverage", "mean"))
              .reset_index().sort_values("median_rel_error"))
    summ.to_csv(RESULTS / "m8_oed_summary.csv", index=False)
    return sdf, summ, picks


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--replicates", type=int, default=30)
    args = ap.parse_args(argv)
    rng = np.random.default_rng(SEED)
    pd.set_option("display.width", 200)

    print("=" * 92)
    print("PART 1 - identifiability threshold (continuous conditioning sweep)")
    print("=" * 92)
    s1 = part1_conditioning(args.replicates, rng)
    print(s1.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    print()
    print("=" * 92)
    print("PART 2 - optimal experimental design (eps = 1e-3, near-unidentifiable)")
    print("=" * 92)
    sdf, s2, picks = part2_oed(args.replicates, rng)
    print("candidate scores:")
    print(sdf.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()
    print("selections:", picks)
    print()
    print("verification by re-calibration:")
    print(s2.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    (RESULTS / "m8_conditioning_oed_manifest.json").write_text(json.dumps({
        "run_utc": datetime.now(timezone.utc).isoformat(), "seed": SEED,
        "replicates": args.replicates, "noise": NOISE, "truth": TRUTH,
        "oed_selections": picks,
    }, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
