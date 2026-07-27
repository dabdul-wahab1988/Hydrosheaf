"""Locked M8 transport OED and PHREEQC kinetics confirmation.

The script deliberately writes new immutable provenance runs and never reads the
exploratory ``results/`` directory. See ``M8_CONFIRMATORY_PROTOCOL.md``.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import io
import json
import logging
import math
import platform
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hydrosheaf.calibration.adapters import (  # noqa: E402
    KineticCalibrationAdapter,
    KineticExperiment,
    TransportCalibrationAdapter,
    TransportExperiment,
)
from hydrosheaf.calibration.glm import PESTGLM  # noqa: E402
from hydrosheaf.config import Config  # noqa: E402
from hydrosheaf.reactive_transport import KineticParameters  # noqa: E402

RUN_ID = "RUN-M8-CONFIRM-20260728-01"
SEED = 2026072801
BOOTSTRAP_SEED = 2026072802
TRUTH = {"dispersivity": 2.0, "decay": 0.005}
DISTANCE_M = 10.0
VELOCITY_M_DAY = 0.1
SOURCE_CONCENTRATION = 1.0
NOISE_REL = 0.03
NOISE_ABS_FLOOR = 0.01
PRIMARY_FD_STEP = 1e-4
FD_STEPS = (1e-3, 1e-4, 1e-5)
PARAMETERS = ("dispersivity", "decay")
BASE_TIMES = (150.0, 180.0, 210.0)
CANDIDATE_TIMES = (50.0, 70.0, 90.0, 105.0, 120.0, 150.0, 180.0, 240.0)


def fixed_designs() -> dict[str, tuple[float, ...]]:
    designs: dict[str, tuple[float, ...]] = {}
    for centre in (60.0, 90.0, 120.0):
        for spread in (5.0, 20.0, 45.0, 70.0):
            times = (
                max(2.0, centre - spread),
                max(3.0, centre - spread / 3.0),
                centre + spread / 3.0,
                centre + spread,
            )
            designs[f"c{int(centre)}_s{int(spread)}"] = tuple(round(v, 1) for v in times)
    designs.update(
        {
            "very_early": (5.0, 15.0, 30.0, 50.0),
            "very_late": (150.0, 180.0, 210.0, 240.0),
            "log_spread": (8.0, 25.0, 75.0, 220.0),
            "dense_front": (30.0, 45.0, 60.0, 80.0),
        }
    )
    return designs


DESIGNS = fixed_designs()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json_hash(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def git_revision() -> str:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        )
        return proc.stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def verify_lock(expected_run_id: str) -> dict:
    lock_path = PROJECT / "m8_confirmatory_protocol.lock.json"
    if not lock_path.is_file():
        raise RuntimeError(f"Protocol lock missing: {lock_path}")
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("run_id") != expected_run_id:
        raise RuntimeError(
            f"Lock run_id {lock.get('run_id')!r} does not match {expected_run_id!r}"
        )
    roots = (("sha256", PROJECT), ("repo_sha256", REPO))
    for key, root in roots:
        resolved_root = root.resolve()
        for relative, expected_hash in lock.get(key, {}).items():
            path = (root / relative).resolve()
            if not path.is_relative_to(resolved_root):
                raise RuntimeError(f"Locked path escapes its root: {relative}")
            actual = sha256_file(path)
            if actual != expected_hash:
                raise RuntimeError(
                    f"Locked file changed: {relative}: {actual} != {expected_hash}"
                )
    return lock


def observation_sigma(clean: np.ndarray) -> np.ndarray:
    clean = np.asarray(clean, dtype=float)
    return np.maximum(NOISE_ABS_FLOOR, NOISE_REL * np.abs(clean))


def transport_forward(times: Sequence[float], params: Mapping[str, float]) -> np.ndarray:
    exp = TransportExperiment(
        id="transport",
        times=list(times),
        observed_concentrations=[0.0] * len(times),
        distance_m=DISTANCE_M,
        velocity_m_day=VELOCITY_M_DAY,
        source_concentration=SOURCE_CONCENTRATION,
    )
    adapter = TransportCalibrationAdapter(
        experiments=[exp],
        params_to_fit=list(PARAMETERS),
        base_dispersivity=TRUTH["dispersivity"],
        base_decay=TRUTH["decay"],
        base_velocity=VELOCITY_M_DAY,
    )
    simulated = adapter.run_model(dict(params))
    return np.array([simulated[f"transport_{i}"] for i in range(len(times))], dtype=float)


def whitened_log_jacobian(
    forward: Callable[[Mapping[str, float]], np.ndarray],
    params: Mapping[str, float],
    sigma: np.ndarray,
    names: Sequence[str],
    step: float = PRIMARY_FD_STEP,
) -> np.ndarray:
    columns: list[np.ndarray] = []
    for name in names:
        plus = dict(params)
        minus = dict(params)
        plus[name] = params[name] * (10.0**step)
        minus[name] = params[name] / (10.0**step)
        derivative = (forward(plus) - forward(minus)) / (2.0 * step)
        columns.append(derivative / sigma)
    return np.column_stack(columns)


def information_diagnostics(jacobian: np.ndarray) -> dict[str, float | int]:
    fim = np.asarray(jacobian, dtype=float).T @ np.asarray(jacobian, dtype=float)
    eigenvalues = np.clip(np.linalg.eigvalsh(fim), 0.0, None)
    lambda_min = float(eigenvalues[0])
    lambda_max = float(eigenvalues[-1])
    tolerance = max(lambda_max * 1e-10, 1e-14)
    rank = int(np.sum(eigenvalues > tolerance))
    condition = lambda_max / lambda_min if lambda_min > 0 else math.inf
    sign, logdet = np.linalg.slogdet(fim)
    if rank == fim.shape[0]:
        covariance = np.linalg.inv(fim)
        sd = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        corr = covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])
        trace_inverse = float(np.trace(covariance))
    else:
        sd = np.full(fim.shape[0], np.inf)
        corr = math.nan
        trace_inverse = math.inf
    cosine = math.nan
    if jacobian.shape[1] == 2:
        norms = np.linalg.norm(jacobian, axis=0)
        if np.all(norms > 0):
            cosine = float(np.dot(jacobian[:, 0], jacobian[:, 1]) / np.prod(norms))
    return {
        "rank": rank,
        "lambda_min": lambda_min,
        "lambda_max": lambda_max,
        "eigen_ratio": lambda_min / lambda_max if lambda_max > 0 else 0.0,
        "condition": condition,
        "logdet_fim": float(logdet) if sign > 0 else -math.inf,
        "trace_inverse": trace_inverse,
        "sd_log10_p1": float(sd[0]),
        "sd_log10_p2": float(sd[1]),
        "parameter_correlation": float(corr),
        "sensitivity_cosine": cosine,
    }


def transport_information(times: Sequence[float], step: float = PRIMARY_FD_STEP) -> dict:
    clean = transport_forward(times, TRUTH)
    sigma = observation_sigma(clean)
    jacobian = whitened_log_jacobian(
        lambda p: transport_forward(times, p), TRUTH, sigma, PARAMETERS, step=step
    )
    result = information_diagnostics(jacobian)
    result.update(
        {
            "sd_log10_dispersivity": result.pop("sd_log10_p1"),
            "sd_log10_decay": result.pop("sd_log10_p2"),
            "fd_step": step,
        }
    )
    return result


def calibration_start(rep: int, radius: float) -> dict[str, float]:
    rng = np.random.default_rng(SEED + 100_000 + rep)
    offsets = rng.uniform(-radius, radius, size=2)
    return {
        name: TRUTH[name] * (10.0 ** float(offsets[i]))
        for i, name in enumerate(PARAMETERS)
    }


def replicate_noise(rep: int, count: int = 4) -> np.ndarray:
    return np.random.default_rng(SEED + rep).normal(size=count)


def calibrate_transport(
    times: Sequence[float],
    standard_noise: np.ndarray,
    start: Mapping[str, float],
    max_nfev: int = 200,
) -> dict:
    clean = transport_forward(times, TRUTH)
    sigma = observation_sigma(clean)
    z = np.asarray(standard_noise, dtype=float)[: len(times)]
    observed = clean + sigma * z
    exp = TransportExperiment(
        id="transport",
        times=list(times),
        observed_concentrations=observed.tolist(),
        weights=(1.0 / sigma).tolist(),
        distance_m=DISTANCE_M,
        velocity_m_day=VELOCITY_M_DAY,
        source_concentration=SOURCE_CONCENTRATION,
    )
    adapter = TransportCalibrationAdapter(
        experiments=[exp],
        params_to_fit=list(PARAMETERS),
        base_dispersivity=float(start["dispersivity"]),
        base_decay=float(start["decay"]),
        base_velocity=VELOCITY_M_DAY,
    )
    parameters = adapter.get_parameters()
    for parameter in parameters:
        parameter.prior_mean = None
        parameter.prior_sigma = None
    pest = PESTGLM(parameters, adapter.get_observations(), adapter.run_model)
    with contextlib.redirect_stdout(io.StringIO()):
        result = pest.calibrate(max_nfev=max_nfev)
    output: dict[str, float | int | bool | str] = {
        "success": bool(result["success"]),
        "phi": float(result["phi"]),
        "n_iterations": int(result["n_iterations"]),
        "covariance_method": str(result["covariance_method"]),
    }
    for name in PARAMETERS:
        truth = TRUTH[name]
        estimate = float(result["optimal_parameters"].get(name, math.nan))
        halfwidth = float(result["parameter_uncertainties_95pc"].get(name, math.nan))
        output[f"start_{name}"] = float(start[name])
        output[f"estimate_{name}"] = estimate
        output[f"abs_log10_error_{name}"] = abs(math.log10(estimate / truth)) if estimate > 0 else math.inf
        output[f"relative_error_{name}"] = abs(estimate - truth) / truth
        output[f"covered_95pc_{name}"] = bool(
            math.isfinite(halfwidth) and abs(estimate - truth) <= halfwidth
        )
        output[f"halfwidth_95pc_{name}"] = halfwidth
    return output


def bootstrap_median_ci(values: Sequence[float], seed: int, samples: int) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return math.nan, math.nan, math.nan
    rng = np.random.default_rng(seed)
    boot = np.empty(samples, dtype=float)
    for index in range(samples):
        boot[index] = np.median(rng.choice(array, size=array.size, replace=True))
    return float(np.median(array)), float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def run_design_sweep(replicates: int, bootstrap_samples: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    information = {name: transport_information(times) for name, times in DESIGNS.items()}
    rows: list[dict] = []
    for rep in range(replicates):
        noise = replicate_noise(rep)
        start = calibration_start(rep, radius=0.25)
        for design, times in DESIGNS.items():
            row = {
                "replicate": rep,
                "design": design,
                "times_days": json.dumps(times),
                **information[design],
                **calibrate_transport(times, noise, start),
            }
            rows.append(row)
    raw = pd.DataFrame(rows)

    summary_rows: list[dict] = []
    for design, group in raw.groupby("design", sort=True):
        for p_index, parameter in enumerate(PARAMETERS):
            values = group[f"abs_log10_error_{parameter}"].to_numpy(float)
            median, lower, upper = bootstrap_median_ci(
                values,
                BOOTSTRAP_SEED + p_index * 10_000 + list(sorted(DESIGNS)).index(design),
                bootstrap_samples,
            )
            summary_rows.append(
                {
                    "design": design,
                    "times_days": group["times_days"].iloc[0],
                    "parameter": parameter,
                    "n": int(len(group)),
                    "median_abs_log10_error": median,
                    "ci95_low": lower,
                    "ci95_high": upper,
                    "mean_coverage_95pc": float(group[f"covered_95pc_{parameter}"].mean()),
                    "success_rate": float(group["success"].mean()),
                    "condition": float(group["condition"].iloc[0]),
                    "lambda_min": float(group["lambda_min"].iloc[0]),
                    "sd_log10_parameter": float(
                        group[f"sd_log10_{parameter}"].iloc[0]
                    ),
                }
            )
    summary = pd.DataFrame(summary_rows)

    correlation_rows: list[dict] = []
    metric_frame = pd.DataFrame(
        [
            {
                "design": design,
                "log10_condition": math.log10(float(values["condition"])),
                "sd_log10_dispersivity": values["sd_log10_dispersivity"],
                "sd_log10_decay": values["sd_log10_decay"],
            }
            for design, values in information.items()
        ]
    )
    for rep, rep_frame in raw.groupby("replicate"):
        merged = rep_frame.merge(metric_frame, on="design", suffixes=("", "_metric"))
        for parameter in PARAMETERS:
            error = merged[f"abs_log10_error_{parameter}"].to_numpy(float)
            for metric in ("log10_condition", f"sd_log10_{parameter}"):
                rho = spearmanr(merged[metric].to_numpy(float), error).statistic
                correlation_rows.append(
                    {
                        "replicate": int(rep),
                        "parameter": parameter,
                        "metric": metric,
                        "spearman_rho": float(rho),
                    }
                )
    correlations = pd.DataFrame(correlation_rows)
    return raw, summary, correlations


def score_candidates() -> tuple[pd.DataFrame, dict[str, float], pd.DataFrame]:
    score_rows: list[dict] = []
    stability_rows: list[dict] = []
    for candidate in CANDIDATE_TIMES:
        times = tuple(sorted((*BASE_TIMES, candidate)))
        for step in FD_STEPS:
            diagnostics = transport_information(times, step=step)
            stability_rows.append({"candidate_time_days": candidate, **diagnostics})
            if step == PRIMARY_FD_STEP:
                score_rows.append(
                    {"candidate_time_days": candidate, "times_days": json.dumps(times), **diagnostics}
                )
    scores = pd.DataFrame(score_rows).sort_values("candidate_time_days")
    picks = {
        "dispersivity_target": float(
            scores.loc[scores["sd_log10_dispersivity"].idxmin(), "candidate_time_days"]
        ),
        "decay_target": float(
            scores.loc[scores["sd_log10_decay"].idxmin(), "candidate_time_days"]
        ),
        "D_optimal": float(scores.loc[scores["logdet_fim"].idxmax(), "candidate_time_days"]),
        "A_optimal": float(scores.loc[scores["trace_inverse"].idxmin(), "candidate_time_days"]),
        "E_optimal": float(scores.loc[scores["lambda_min"].idxmax(), "candidate_time_days"]),
        "balanced": float(
            scores.loc[
                scores[["sd_log10_dispersivity", "sd_log10_decay"]].max(axis=1).idxmin(),
                "candidate_time_days",
            ]
        ),
        "worst_joint": float(scores.loc[scores["trace_inverse"].idxmax(), "candidate_time_days"]),
    }
    return scores, picks, pd.DataFrame(stability_rows)


def run_oed(
    replicates: int, bootstrap_samples: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, float]]:
    scores, picks, stability = score_candidates()
    rows: list[dict] = []
    for regime, radius in (("local", 0.25), ("distant", 1.0)):
        for rep in range(replicates):
            noise = replicate_noise(10_000 + rep)
            start = calibration_start(10_000 + rep, radius=radius)
            chosen = dict(picks)
            random_rng = np.random.default_rng(SEED + 300_000 + rep)
            chosen["random"] = float(random_rng.choice(CANDIDATE_TIMES))
            for strategy, candidate in chosen.items():
                times = tuple(sorted((*BASE_TIMES, candidate)))
                rows.append(
                    {
                        "start_regime": regime,
                        "replicate": rep,
                        "strategy": strategy,
                        "candidate_time_days": candidate,
                        "times_days": json.dumps(times),
                        **calibrate_transport(times, noise, start),
                    }
                )
            rows.append(
                {
                    "start_regime": regime,
                    "replicate": rep,
                    "strategy": "no_new_measurement",
                    "candidate_time_days": math.nan,
                    "times_days": json.dumps(BASE_TIMES),
                    **calibrate_transport(BASE_TIMES, noise[:3], start),
                }
            )
    raw = pd.DataFrame(rows)

    summary_rows: list[dict] = []
    for (regime, strategy), group in raw.groupby(["start_regime", "strategy"], sort=True):
        for p_index, parameter in enumerate(PARAMETERS):
            values = group[f"abs_log10_error_{parameter}"].to_numpy(float)
            seed = BOOTSTRAP_SEED + 100_000 + p_index * 10_000 + sum(map(ord, regime + strategy))
            median, lower, upper = bootstrap_median_ci(values, seed, bootstrap_samples)
            summary_rows.append(
                {
                    "start_regime": regime,
                    "strategy": strategy,
                    "candidate_time_days": group["candidate_time_days"].dropna().iloc[0]
                    if group["candidate_time_days"].notna().any()
                    else math.nan,
                    "parameter": parameter,
                    "n": int(len(group)),
                    "median_abs_log10_error": median,
                    "ci95_low": lower,
                    "ci95_high": upper,
                    "mean_coverage_95pc": float(group[f"covered_95pc_{parameter}"].mean()),
                    "success_rate": float(group["success"].mean()),
                }
            )
    summary = pd.DataFrame(summary_rows)

    comparisons: list[dict] = []
    comparison_pairs = [
        (strategy, "no_new_measurement") for strategy in picks
    ] + [
        ("dispersivity_target", "decay_target"),
        ("decay_target", "dispersivity_target"),
    ]
    for regime in ("local", "distant"):
        frame = raw.loc[raw["start_regime"] == regime]
        for parameter_index, parameter in enumerate(PARAMETERS):
            pivot = frame.pivot(
                index="replicate", columns="strategy", values=f"abs_log10_error_{parameter}"
            )
            for left, right in comparison_pairs:
                difference = (pivot[left] - pivot[right]).to_numpy(float)
                median, lower, upper = bootstrap_median_ci(
                    difference,
                    BOOTSTRAP_SEED
                    + 200_000
                    + parameter_index * 10_000
                    + sum(map(ord, regime + left + right)),
                    bootstrap_samples,
                )
                comparisons.append(
                    {
                        "start_regime": regime,
                        "parameter": parameter,
                        "left_strategy": left,
                        "right_strategy": right,
                        "median_paired_error_difference": median,
                        "ci95_low": lower,
                        "ci95_high": upper,
                        "favours_left_if_negative": True,
                    }
                )
    return raw, summary, pd.DataFrame(comparisons), stability, picks


KINETIC_TRUTH = {"log_k_calcite": 1e-10, "log_A_calcite": 0.1}


def kinetic_adapter(experiments: Sequence[KineticExperiment]) -> KineticCalibrationAdapter:
    return KineticCalibrationAdapter(
        base_params={
            "calcite": KineticParameters(
                reaction_name="calcite",
                rate_constant=KINETIC_TRUTH["log_k_calcite"],
                surface_area=KINETIC_TRUTH["log_A_calcite"],
            )
        },
        experiments=list(experiments),
        config=Config(),
        params_to_fit=["calcite:k", "calcite:A"],
    )


def kinetic_experiments() -> list[KineticExperiment]:
    initial = {
        "pH": 6.5,
        "Ca": 0.05,
        "Mg": 0.02,
        "Na": 0.10,
        "K": 0.01,
        "Cl": 0.10,
        "SO4": 0.01,
        "NO3": 0.01,
        "HCO3": 0.10,
        "F": 0.0,
        "Fe": 0.0,
        "PO4": 0.0,
    }
    return [
        KineticExperiment(
            id=f"calcite_t{int(days)}",
            initial_solution=dict(initial),
            residence_time_days=days,
            reaction_extents=[1.0],
            reaction_labels=["calcite"],
            observations={"Ca": 0.0, "HCO3": 0.0},
            units="mmol/L",
        )
        for days in (10.0, 30.0, 100.0, 300.0)
    ]


def kinetic_forward(experiments: Sequence[KineticExperiment], params: Mapping[str, float]) -> np.ndarray:
    adapter = kinetic_adapter(experiments)
    result = adapter.run_model(dict(params))
    names = [
        f"{experiment.id}_{ion}"
        for experiment in experiments
        for ion in ("Ca", "HCO3")
    ]
    values = np.array([result[name] for name in names], dtype=float)
    if not np.all(np.isfinite(values)) or np.any(np.abs(values) >= 1e5):
        raise RuntimeError(f"Invalid PHREEQC kinetic output: {values}")
    return values


def run_kinetics() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    experiments = kinetic_experiments()
    designs = {
        "single_time_data_only": experiments[:1],
        "multi_time_data_only": experiments,
    }
    rows: list[dict] = []
    truth_outputs: dict[str, np.ndarray] = {}
    for name, selected in designs.items():
        def forward(
            params: Mapping[str, float],
            selected_experiments: Sequence[KineticExperiment] = selected,
        ) -> np.ndarray:
            return kinetic_forward(selected_experiments, params)

        clean = forward(KINETIC_TRUTH)
        truth_outputs[name] = clean
        sigma = np.full(clean.shape, 0.005, dtype=float)
        jacobian = whitened_log_jacobian(
            forward,
            KINETIC_TRUTH,
            sigma,
            ("log_k_calcite", "log_A_calcite"),
            step=PRIMARY_FD_STEP,
        )
        diagnostics = information_diagnostics(jacobian)
        rows.append(
            {
                "design": name,
                "n_concentration_observations": int(clean.size),
                "independent_surface_area_observation": False,
                **diagnostics,
            }
        )
        if name == "multi_time_data_only":
            augmented = np.vstack([jacobian, np.array([[0.0, 1.0 / 0.10]])])
            augmented_diagnostics = information_diagnostics(augmented)
            rows.append(
                {
                    "design": "multi_time_plus_surface_area",
                    "n_concentration_observations": int(clean.size),
                    "independent_surface_area_observation": True,
                    **augmented_diagnostics,
                }
            )

    base = truth_outputs["multi_time_data_only"]
    sigma = np.full(base.shape, 0.005, dtype=float)
    ridge_rows: list[dict] = []
    for factor in (0.1, 0.3, 1.0, 3.0, 10.0):
        params = {
            "log_k_calcite": KINETIC_TRUTH["log_k_calcite"] * factor,
            "log_A_calcite": KINETIC_TRUTH["log_A_calcite"] / factor,
        }
        predicted = kinetic_forward(experiments, params)
        ridge_rows.append(
            {
                "case": "constant_product_ridge",
                "factor_on_k": factor,
                "k": params["log_k_calcite"],
                "surface_area": params["log_A_calcite"],
                "product": params["log_k_calcite"] * params["log_A_calcite"],
                "max_abs_concentration_difference": float(np.max(np.abs(predicted - base))),
                "max_standardized_difference": float(np.max(np.abs((predicted - base) / sigma))),
            }
        )
    off_ridge = dict(KINETIC_TRUTH)
    off_ridge["log_k_calcite"] *= 2.0
    off_prediction = kinetic_forward(experiments, off_ridge)
    ridge_rows.append(
        {
            "case": "doubled_product_off_ridge",
            "factor_on_k": 2.0,
            "k": off_ridge["log_k_calcite"],
            "surface_area": off_ridge["log_A_calcite"],
            "product": off_ridge["log_k_calcite"] * off_ridge["log_A_calcite"],
            "max_abs_concentration_difference": float(np.max(np.abs(off_prediction - base))),
            "max_standardized_difference": float(np.max(np.abs((off_prediction - base) / sigma))),
        }
    )
    summary = pd.DataFrame(rows)
    ridge = pd.DataFrame(ridge_rows)
    data_only = summary.loc[summary["design"] == "multi_time_data_only"].iloc[0]
    constrained = summary.loc[summary["design"] == "multi_time_plus_surface_area"].iloc[0]
    ridge_max = ridge.loc[
        ridge["case"] == "constant_product_ridge", "max_standardized_difference"
    ].max()
    off_response = ridge.loc[
        ridge["case"] == "doubled_product_off_ridge", "max_standardized_difference"
    ].iloc[0]
    decision = {
        "structural_nonidentifiability_supported": bool(
            int(data_only["rank"]) == 1
            and abs(float(data_only["sensitivity_cosine"])) >= 0.999
            and float(data_only["eigen_ratio"]) <= 1e-8
            and float(ridge_max) <= 1e-5
            and float(off_response) >= 1e-3
        ),
        "surface_area_constraint_restores_rank": bool(
            int(constrained["rank"]) == 2 and float(constrained["eigen_ratio"]) > 1e-10
        ),
        "ridge_max_standardized_difference": float(ridge_max),
        "off_ridge_standardized_response": float(off_response),
    }
    return summary, ridge, decision


def transport_claim_decision(comparisons: pd.DataFrame, picks: Mapping[str, float]) -> dict:
    local = comparisons.loc[comparisons["start_regime"] == "local"]

    def comparison(parameter: str, left: str, right: str) -> pd.Series:
        match = local.loc[
            (local["parameter"] == parameter)
            & (local["left_strategy"] == left)
            & (local["right_strategy"] == right)
        ]
        if len(match) != 1:
            raise RuntimeError(f"Missing comparison {parameter}: {left} vs {right}")
        return match.iloc[0]

    disp = comparison("dispersivity", "dispersivity_target", "decay_target")
    decay = comparison("decay", "decay_target", "dispersivity_target")
    distinct = picks["dispersivity_target"] != picks["decay_target"]
    disp_supported = float(disp["ci95_high"]) < 0.0
    decay_supported = float(decay["ci95_high"]) < 0.0
    if distinct and disp_supported and decay_supported:
        status = "SUPPORTED"
    elif distinct and (disp_supported or decay_supported):
        status = "PARTIAL"
    else:
        status = "NOT_SUPPORTED"
    return {
        "target_specific_candidate_times_differ": bool(distinct),
        "dispersivity_direction_supported": bool(disp_supported),
        "decay_direction_supported": bool(decay_supported),
        "tradeoff_status": status,
        "picks": dict(picks),
    }


def make_figure(
    design_summary: pd.DataFrame,
    correlations: pd.DataFrame,
    candidate_scores: pd.DataFrame,
    oed_summary: pd.DataFrame,
    kinetics_summary: pd.DataFrame,
    output_path: Path,
) -> None:
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    colours = {"dispersivity": "#1f77b4", "decay": "#d62728"}
    for parameter in PARAMETERS:
        frame = design_summary.loc[design_summary["parameter"] == parameter]
        ax.scatter(
            np.log10(frame["condition"].to_numpy(float)),
            frame["median_abs_log10_error"].to_numpy(float),
            label=parameter,
            color=colours[parameter],
            alpha=0.85,
        )
    ax.set_xlabel("log10 condition number (whitened log-parameter FIM)")
    ax.set_ylabel("Median absolute log10 error")
    ax.set_title("A  Fixed-design recovery is parameter-specific")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)

    ax = axes[0, 1]
    ax.plot(
        candidate_scores["candidate_time_days"],
        candidate_scores["sd_log10_dispersivity"],
        "o-",
        color=colours["dispersivity"],
        label="dispersivity",
    )
    ax.plot(
        candidate_scores["candidate_time_days"],
        candidate_scores["sd_log10_decay"],
        "s-",
        color=colours["decay"],
        label="decay",
    )
    ax.set_xlabel("Added sampling time (d)")
    ax.set_ylabel("Predicted marginal SD (log10 units)")
    ax.set_title("B  Candidate value depends on target parameter")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)

    ax = axes[1, 0]
    local = oed_summary.loc[
        (oed_summary["start_regime"] == "local")
        & oed_summary["strategy"].isin(
            ["dispersivity_target", "decay_target", "balanced", "random", "no_new_measurement"]
        )
    ].copy()
    strategies = [
        "dispersivity_target",
        "decay_target",
        "balanced",
        "random",
        "no_new_measurement",
    ]
    x = np.arange(len(strategies), dtype=float)
    width = 0.36
    for offset, parameter in ((-width / 2, "dispersivity"), (width / 2, "decay")):
        frame = local.loc[local["parameter"] == parameter].set_index("strategy").reindex(strategies)
        ax.bar(
            x + offset,
            frame["median_abs_log10_error"],
            width,
            color=colours[parameter],
            label=parameter,
        )
    ax.set_xticks(x, ["disp. target", "decay target", "balanced", "random", "no new"], rotation=18)
    ax.set_ylabel("Median absolute log10 error")
    ax.set_title("C  Realised recovery under paired local starts")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)

    ax = axes[1, 1]
    labels = kinetics_summary["design"].str.replace("_", " ").tolist()
    ratios = np.maximum(kinetics_summary["eigen_ratio"].to_numpy(float), 1e-18)
    bars = ax.bar(np.arange(len(labels)), ratios, color=["#999999", "#666666", "#2ca02c"])
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=15)
    ax.set_ylabel("Smallest/largest FIM eigenvalue")
    ax.set_title("D  Kinetic k*A confounding needs external A evidence")
    ax.axhline(1e-10, color="black", linestyle="--", linewidth=0.8, label="rank tolerance")
    ax.legend(frameon=False)
    for bar, rank in zip(bars, kinetics_summary["rank"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"rank {int(rank)}", ha="center", va="bottom")
    ax.grid(axis="y", alpha=0.2)

    fig.suptitle(
        "M8 confirmatory controlled-synthetic benchmark",
        fontsize=13,
        fontweight="bold",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    fig.savefig(output_path.with_suffix(".pdf"))
    plt.close(fig)


def environment_record() -> dict:
    return {
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "git_commit": git_revision(),
    }


def append_run_ledger(record: Mapping[str, object]) -> None:
    provenance = PROJECT / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    ledger = provenance / "run_ledger.csv"
    fields = [
        "run_id",
        "status",
        "started_at_utc",
        "ended_at_utc",
        "git_commit",
        "config_hash",
        "input_manifest_hash",
        "environment_hash",
        "database_snapshot_id",
    ]
    exists = ledger.exists()
    with ledger.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not exists:
            writer.writeheader()
        writer.writerow({field: record.get(field, "") for field in fields})


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, lineterminator="\n")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=250)
    parser.add_argument("--bootstrap-samples", type=int, default=2000)
    parser.add_argument("--run-id", default=RUN_ID)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.replicates < 2:
        raise ValueError("replicates must be at least 2")
    if args.bootstrap_samples < 100:
        raise ValueError("bootstrap-samples must be at least 100")
    lock = verify_lock(args.run_id)
    if args.dry_run:
        print(json.dumps({"status": "DRY-RUN-PASS", "run_id": args.run_id, "lock": lock}, indent=2))
        return 0

    run_root = PROJECT / "provenance" / "runs" / args.run_id
    if run_root.exists():
        raise FileExistsError(f"Immutable run already exists: {run_root}")
    run_root.mkdir(parents=True)
    started = utc_now()
    environment = environment_record()
    config = {
        "run_id": args.run_id,
        "replicates": args.replicates,
        "bootstrap_samples": args.bootstrap_samples,
        "seed": SEED,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "truth": TRUTH,
        "noise_rel": NOISE_REL,
        "noise_abs_floor": NOISE_ABS_FLOOR,
        "fd_steps": FD_STEPS,
    }
    manifest = {
        "schema_version": "1.0",
        "run_id": args.run_id,
        "status": "RUNNING",
        "started_at_utc": started,
        "ended_at_utc": None,
        "config": config,
        "lock_sha256": sha256_file(PROJECT / "m8_confirmatory_protocol.lock.json"),
        "environment": environment,
        "artifacts": [],
    }
    manifest_path = run_root / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    start_time = time.perf_counter()
    logging.getLogger("hydrosheaf.calibration.glm").setLevel(logging.ERROR)

    try:
        raw_design, design_summary, correlations = run_design_sweep(
            args.replicates, args.bootstrap_samples
        )
        candidate_scores, _, fd_stability = score_candidates()
        raw_oed, oed_summary, comparisons, _, picks = run_oed(
            args.replicates, args.bootstrap_samples
        )
        kinetics_summary, kinetics_ridge, kinetics_decision = run_kinetics()
        transport_decision = transport_claim_decision(comparisons, picks)
        claim_decision = {
            "transport": transport_decision,
            "kinetics": kinetics_decision,
            "overall": (
                "TARGET_DEPENDENCE_CONFIRMED_WITH_STRUCTURAL_KINETIC_LIMIT"
                if transport_decision["tradeoff_status"] == "SUPPORTED"
                and kinetics_decision["structural_nonidentifiability_supported"]
                and kinetics_decision["surface_area_constraint_restores_rank"]
                else "SCOPED_OR_NEGATIVE_RESULT_REQUIRED"
            ),
        }

        raw_files = {
            "transport_design_replicates.csv": raw_design,
            "transport_design_correlations.csv": correlations,
            "transport_candidate_scores.csv": candidate_scores,
            "transport_fd_stability.csv": fd_stability,
            "transport_oed_replicates.csv": raw_oed,
            "transport_oed_paired_comparisons.csv": comparisons,
            "kinetics_constant_product_ridge.csv": kinetics_ridge,
        }
        for name, frame in raw_files.items():
            write_csv(frame, run_root / name)
        (run_root / "claim_decision.json").write_text(
            json.dumps(claim_decision, indent=2), encoding="utf-8"
        )
        (run_root / "environment.json").write_text(
            json.dumps(environment, indent=2), encoding="utf-8"
        )

        artifact_root = PROJECT / "manuscript" / "artifacts"
        artifact_paths = {
            "TAB-M8-TRANSPORT": artifact_root / "m8_transport_parameter_summary.csv",
            "TAB-M8-OED": artifact_root / "m8_transport_oed_summary.csv",
            "TAB-M8-KINETICS": artifact_root / "m8_kinetics_structural_summary.csv",
            "FIG-M8-CONFIRMATORY": artifact_root / "m8_confirmatory_figure.png",
        }
        write_csv(design_summary, artifact_paths["TAB-M8-TRANSPORT"])
        write_csv(oed_summary, artifact_paths["TAB-M8-OED"])
        write_csv(kinetics_summary, artifact_paths["TAB-M8-KINETICS"])
        make_figure(
            design_summary,
            correlations,
            candidate_scores,
            oed_summary,
            kinetics_summary,
            artifact_paths["FIG-M8-CONFIRMATORY"],
        )
        (artifact_root / "m8_transport_candidate_scores.csv").write_text(
            candidate_scores.to_csv(index=False, lineterminator="\n"), encoding="utf-8"
        )
        (artifact_root / "m8_transport_paired_comparisons.csv").write_text(
            comparisons.to_csv(index=False, lineterminator="\n"), encoding="utf-8"
        )
        (artifact_root / "m8_claim_decision.json").write_text(
            json.dumps(claim_decision, indent=2), encoding="utf-8"
        )

        manifest["status"] = "PASS"
        manifest["ended_at_utc"] = utc_now()
        manifest["runtime_seconds"] = time.perf_counter() - start_time
        for artifact_id, path in artifact_paths.items():
            manifest["artifacts"].append(
                {
                    "id": artifact_id,
                    "path": path.relative_to(PROJECT).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        append_run_ledger(
            {
                "run_id": args.run_id,
                "status": "PASS",
                "started_at_utc": started,
                "ended_at_utc": manifest["ended_at_utc"],
                "git_commit": environment["git_commit"],
                "config_hash": stable_json_hash(config),
                "input_manifest_hash": manifest["lock_sha256"],
                "environment_hash": stable_json_hash(environment),
                "database_snapshot_id": "NOT-APPLICABLE-TRACK-B-FILESYSTEM",
            }
        )
        print(json.dumps(claim_decision, indent=2))
        print(f"PASS {args.run_id} in {manifest['runtime_seconds']:.1f} s")
        return 0
    except Exception as exc:
        manifest["status"] = "FAIL"
        manifest["ended_at_utc"] = utc_now()
        manifest["runtime_seconds"] = time.perf_counter() - start_time
        manifest["error"] = str(exc)
        manifest["traceback"] = traceback.format_exc()
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        append_run_ledger(
            {
                "run_id": args.run_id,
                "status": "FAIL",
                "started_at_utc": started,
                "ended_at_utc": manifest["ended_at_utc"],
                "git_commit": environment["git_commit"],
                "config_hash": stable_json_hash(config),
                "input_manifest_hash": manifest["lock_sha256"],
                "environment_hash": stable_json_hash(environment),
                "database_snapshot_id": "NOT-APPLICABLE-TRACK-B-FILESYSTEM",
            }
        )
        raise


if __name__ == "__main__":
    raise SystemExit(main())
