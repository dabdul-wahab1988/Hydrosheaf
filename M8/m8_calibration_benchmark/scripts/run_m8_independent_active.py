"""Independent numerical truth and active-learning portability benchmark."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import math
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import factorized

PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hydrosheaf.calibration.active_learning import rank_next_measurements  # noqa: E402
from hydrosheaf.calibration.adapters import (  # noqa: E402
    TransportCalibrationAdapter,
    TransportExperiment,
)
from hydrosheaf.calibration.glm import PESTGLM  # noqa: E402
from hydrosheaf.graph.types import Edge  # noqa: E402

from run_m8_confirmatory import (  # noqa: E402
    BASE_TIMES,
    CANDIDATE_TIMES,
    DISTANCE_M,
    NOISE_ABS_FLOOR,
    NOISE_REL,
    PARAMETERS,
    SOURCE_CONCENTRATION,
    TRUTH,
    VELOCITY_M_DAY,
    score_candidates,
)

RUN_ID = "RUN-M8-INDEPENDENT-20260728-01"
DEVELOPMENT_REPLICATES = 80
TEST_REPLICATES = 250
BOOTSTRAP_SAMPLES = 5000
DT_DAYS = 0.05
TRUTH_CELLS = 240
GRID_CELLS = (120, 240, 480)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "UNAVAILABLE"


def _verify_lock(
    development_replicates: int,
    test_replicates: int,
    bootstrap_samples: int,
) -> dict:
    path = PROJECT / "m8_independent_active_protocol.lock.json"
    lock = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "development_replicates": int(development_replicates),
        "test_replicates": int(test_replicates),
        "bootstrap_samples": int(bootstrap_samples),
    }
    if lock.get("run_id") != RUN_ID:
        raise RuntimeError("Independent-model protocol lock has the wrong run ID.")
    for key, value in expected.items():
        if int(lock.get(key, -1)) != value:
            raise RuntimeError(f"{key} does not match the locked protocol.")
    for relative, expected_hash in lock.get("sha256", {}).items():
        file_path = (REPO / relative).resolve()
        if not file_path.is_relative_to(REPO.resolve()):
            raise RuntimeError(f"Locked path escapes repository: {relative}")
        actual = _sha256(file_path)
        if actual != expected_hash:
            raise RuntimeError(
                f"Locked file changed: {relative}: {actual} != {expected_hash}"
            )
    return lock


def numerical_transport_truth(
    times: Sequence[float],
    *,
    dispersivity: float = TRUTH["dispersivity"],
    decay: float = TRUTH["decay"],
    cells: int = TRUTH_CELLS,
    dt_days: float = DT_DAYS,
) -> np.ndarray:
    """Implicit finite-volume/upwind truth, independent of HydroSheaf transport."""
    targets = np.asarray(times, dtype=float)
    if np.any(targets < 0.0) or cells < 10 or dt_days <= 0.0:
        raise ValueError("Invalid numerical transport design.")
    dx = DISTANCE_M / int(cells)
    dispersion = float(dispersivity) * VELOCITY_M_DAY
    diffusion_number = dispersion * dt_days / (dx * dx)
    courant = VELOCITY_M_DAY * dt_days / dx
    matrix = lil_matrix((cells + 1, cells + 1), dtype=float)
    matrix[0, 0] = 1.0
    for index in range(1, cells):
        matrix[index, index - 1] = -(diffusion_number + courant)
        matrix[index, index] = 1.0 + 2.0 * diffusion_number + courant + decay * dt_days
        matrix[index, index + 1] = -diffusion_number
    matrix[cells, cells - 1] = -1.0
    matrix[cells, cells] = 1.0
    solve = factorized(matrix.tocsc())

    state = np.zeros(cells + 1, dtype=float)
    state[0] = SOURCE_CONCENTRATION
    max_time = float(np.max(targets)) if len(targets) else 0.0
    n_steps = int(math.ceil(max_time / dt_days))
    time_grid = np.arange(n_steps + 1, dtype=float) * dt_days
    outlet = np.empty(n_steps + 1, dtype=float)
    outlet[0] = state[-1]
    for step in range(1, n_steps + 1):
        rhs = state.copy()
        rhs[0] = SOURCE_CONCENTRATION
        rhs[-1] = 0.0
        state = solve(rhs)
        outlet[step] = state[-1]
    return np.interp(targets, time_grid, outlet)


def observation_sigma(clean: np.ndarray) -> np.ndarray:
    clean = np.asarray(clean, dtype=float)
    return np.maximum(NOISE_ABS_FLOOR, NOISE_REL * np.abs(clean))


def _start(rep: int, family_seed: int) -> dict[str, float]:
    offsets = np.random.default_rng(family_seed + rep).uniform(-0.25, 0.25, size=2)
    return {
        name: float(TRUTH[name] * 10.0 ** offsets[index])
        for index, name in enumerate(PARAMETERS)
    }


def calibrate_observed(
    times: Sequence[float],
    observed: np.ndarray,
    sigma: np.ndarray,
    start: Mapping[str, float],
) -> dict:
    experiment = TransportExperiment(
        id="transport",
        times=list(map(float, times)),
        observed_concentrations=np.asarray(observed, dtype=float).tolist(),
        weights=(1.0 / np.asarray(sigma, dtype=float)).tolist(),
        distance_m=DISTANCE_M,
        velocity_m_day=VELOCITY_M_DAY,
        source_concentration=SOURCE_CONCENTRATION,
    )
    adapter = TransportCalibrationAdapter(
        experiments=[experiment],
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
        result = pest.calibrate(max_nfev=200)
    output: dict[str, object] = {
        "success": bool(result["success"]),
        "phi": float(result["phi"]),
        "covariance_method": str(result["covariance_method"]),
    }
    joint = []
    for name in PARAMETERS:
        estimate = float(result["optimal_parameters"].get(name, math.nan))
        halfwidth = float(result["parameter_uncertainties_95pc"].get(name, math.nan))
        error = (
            abs(math.log10(estimate / TRUTH[name]))
            if estimate > 0.0
            else math.inf
        )
        output[f"estimate_{name}"] = estimate
        output[f"abs_log10_error_{name}"] = error
        output[f"covered_95pc_{name}"] = bool(
            math.isfinite(halfwidth) and abs(estimate - TRUTH[name]) <= halfwidth
        )
        joint.append(error)
    output["joint_abs_log10_error"] = float(np.mean(joint))
    return output


def _noise_for_times(rep: int, times: Sequence[float], family_seed: int) -> np.ndarray:
    all_times = sorted(set((*BASE_TIMES, *CANDIDATE_TIMES)))
    values = np.random.default_rng(family_seed + rep).normal(size=len(all_times))
    by_time = dict(zip(all_times, values))
    return np.asarray([by_time[float(time)] for time in times], dtype=float)


def _calibrate_design(
    rep: int,
    times: Sequence[float],
    truth_by_time: Mapping[float, float],
    *,
    family_seed: int,
) -> dict:
    clean = np.asarray([truth_by_time[float(time)] for time in times], dtype=float)
    sigma = observation_sigma(clean)
    observed = clean + sigma * _noise_for_times(rep, times, family_seed)
    return calibrate_observed(times, observed, sigma, _start(rep, family_seed + 500_000))


def _bootstrap_difference(
    raw: pd.DataFrame,
    strategy: str,
    metric: str,
    *,
    samples: int,
    seed: int,
) -> dict:
    pivot = raw.pivot(index="replicate", columns="strategy", values=metric)
    difference = (pivot[strategy] - pivot["no_new_measurement"]).to_numpy(float)
    rng = np.random.default_rng(seed)
    boot = np.empty(int(samples), dtype=float)
    for index in range(int(samples)):
        boot[index] = float(
            np.median(rng.choice(difference, size=len(difference), replace=True))
        )
    return {
        "strategy": strategy,
        "reference": "no_new_measurement",
        "metric": metric,
        "median_paired_difference": float(np.median(difference)),
        "ci95_low": float(np.quantile(boot, 0.025)),
        "ci95_high": float(np.quantile(boot, 0.975)),
    }


def heuristic_portability_test() -> tuple[dict, dict]:
    candidate_ids = [f"transport_time_{int(time)}d" for time in CANDIDATE_TIMES]
    report = {
        "variants": {
            name: {"selected_edge_ids": list(candidate_ids)}
            for name in ("baseline", "null_model_defaults", "assumption_calibrated")
        },
        "independent_validation": True,
        "manuscript_claim_allowed": False,
    }
    edges = [
        Edge(edge_id=edge_id, u="base_design", v=edge_id, attrs={})
        for edge_id in candidate_ids
    ]
    result = rank_next_measurements(
        benchmark_report=report,
        candidate_edges=edges,
        top_k=len(edges),
    )
    recommendations = sorted(
        result.get("recommendations", []),
        key=lambda item: str(item.get("edge_id", "")),
    )
    for rank, recommendation in enumerate(recommendations, start=1):
        recommendation["rank"] = rank
    result["recommendations"] = recommendations
    priorities = {float(item["priority_score"]) for item in recommendations}
    texts = [
        str(measurement).lower()
        for item in recommendations
        for measurement in item.get("recommended_measurements", [])
    ]
    explicit_transport_action = any(
        "concentration" in text and ("day" in text or "time" in text)
        for text in texts
    )
    actionable = bool(recommendations and len(priorities) > 1 and explicit_transport_action)
    audit = {
        "status": (
            "ACTIONABLE" if actionable else "NOT_ACTIONABLE_FOR_TRANSPORT_TIME_SELECTION"
        ),
        "n_recommendations": len(recommendations),
        "n_unique_priority_scores": len(priorities),
        "explicit_transport_time_action": explicit_transport_action,
        "reason": (
            "The existing heuristic consumes topology evidence classes and does not "
            "consume candidate-time Jacobians or emit a transport concentration time."
        ),
    }
    return result, audit


def _safe_output(output: Path, overwrite: bool) -> None:
    resolved = output.resolve()
    if not resolved.is_relative_to(PROJECT.resolve()):
        raise ValueError("Output must stay within the M8 benchmark directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(f"Output exists: {resolved}")
        if resolved == PROJECT.resolve():
            raise ValueError("Refusing to remove the benchmark root.")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)


def run_benchmark(
    *,
    output: Path,
    development_replicates: int = DEVELOPMENT_REPLICATES,
    test_replicates: int = TEST_REPLICATES,
    bootstrap_samples: int = BOOTSTRAP_SAMPLES,
    overwrite: bool = False,
) -> dict:
    locked = (
        development_replicates == DEVELOPMENT_REPLICATES
        and test_replicates == TEST_REPLICATES
        and bootstrap_samples == BOOTSTRAP_SAMPLES
    )
    if locked:
        _verify_lock(development_replicates, test_replicates, bootstrap_samples)
    _safe_output(output, overwrite)

    all_times = tuple(sorted(set((*BASE_TIMES, *CANDIDATE_TIMES))))
    grid_predictions = {
        cells: numerical_transport_truth(all_times, cells=cells)
        for cells in GRID_CELLS
    }
    reference = grid_predictions[480]
    reference_sigma = observation_sigma(reference)
    convergence_rows = []
    for cells in GRID_CELLS[:-1]:
        standardized = np.abs(grid_predictions[cells] - reference) / reference_sigma
        convergence_rows.append(
            {
                "cells": cells,
                "reference_cells": 480,
                "max_abs_concentration_difference": float(
                    np.max(np.abs(grid_predictions[cells] - reference))
                ),
                "max_standardized_difference": float(np.max(standardized)),
            }
        )
    convergence = pd.DataFrame(convergence_rows)
    convergence_passed = bool(
        convergence.loc[convergence["cells"] == TRUTH_CELLS, "max_standardized_difference"].iloc[0]
        <= 0.25
    )
    truth_by_time = dict(zip(all_times, grid_predictions[TRUTH_CELLS]))

    candidate_scores, analytical_picks, fd_stability = score_candidates()
    e_optimal_time = float(analytical_picks["E_optimal"])

    development_rows: list[dict] = []
    for rep in range(int(development_replicates)):
        for candidate in CANDIDATE_TIMES:
            times = tuple(sorted((*BASE_TIMES, float(candidate))))
            development_rows.append(
                {
                    "replicate": rep,
                    "candidate_time_days": float(candidate),
                    **_calibrate_design(
                        rep,
                        times,
                        truth_by_time,
                        family_seed=2026072810,
                    ),
                }
            )
    development = pd.DataFrame(development_rows)
    oracle_time = float(
        development.groupby("candidate_time_days")["joint_abs_log10_error"]
        .median()
        .idxmin()
    )

    test_rows: list[dict] = []
    for rep in range(int(test_replicates)):
        random_time = float(
            np.random.default_rng(2026072820 + 300_000 + rep).choice(CANDIDATE_TIMES)
        )
        strategies = {
            "E_optimal": e_optimal_time,
            "development_oracle": oracle_time,
            "random": random_time,
            "no_new_measurement": None,
        }
        for strategy, candidate in strategies.items():
            times = (
                BASE_TIMES
                if candidate is None
                else tuple(sorted((*BASE_TIMES, float(candidate))))
            )
            test_rows.append(
                {
                    "replicate": rep,
                    "strategy": strategy,
                    "candidate_time_days": candidate,
                    **_calibrate_design(
                        rep,
                        times,
                        truth_by_time,
                        family_seed=2026072820,
                    ),
                }
            )
    test = pd.DataFrame(test_rows)
    summary_rows = []
    for strategy, group in test.groupby("strategy", sort=True):
        row: dict[str, object] = {
            "strategy": strategy,
            "n": len(group),
            "candidate_time_days": (
                float(group["candidate_time_days"].dropna().iloc[0])
                if strategy not in {"random", "no_new_measurement"}
                else math.nan
            ),
            "median_joint_abs_log10_error": float(
                group["joint_abs_log10_error"].median()
            ),
            "success_rate": float(group["success"].mean()),
        }
        for parameter in PARAMETERS:
            row[f"median_abs_log10_error_{parameter}"] = float(
                group[f"abs_log10_error_{parameter}"].median()
            )
            row[f"coverage_95pc_{parameter}"] = float(
                group[f"covered_95pc_{parameter}"].mean()
            )
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    heuristic_raw, heuristic_audit = heuristic_portability_test()
    summary = pd.concat(
        [
            summary,
            pd.DataFrame(
                [
                    {
                        "strategy": "existing_active_learning_heuristic",
                        "n": 0,
                        "candidate_time_days": math.nan,
                        "median_joint_abs_log10_error": math.nan,
                        "success_rate": math.nan,
                        **{
                            f"median_abs_log10_error_{parameter}": math.nan
                            for parameter in PARAMETERS
                        },
                        **{
                            f"coverage_95pc_{parameter}": math.nan
                            for parameter in PARAMETERS
                        },
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    contrasts = pd.DataFrame(
        [
            _bootstrap_difference(
                test,
                strategy,
                f"abs_log10_error_{parameter}",
                samples=bootstrap_samples,
                seed=2026072830 + 100 * strategy_index + parameter_index,
            )
            for strategy_index, strategy in enumerate(
                ("E_optimal", "development_oracle", "random")
            )
            for parameter_index, parameter in enumerate(PARAMETERS)
        ]
    )
    e_contrasts = contrasts[contrasts["strategy"] == "E_optimal"]
    robust = bool(
        convergence_passed
        and len(e_contrasts) == len(PARAMETERS)
        and (e_contrasts["ci95_high"] < 0.0).all()
    )
    heuristic_supported = False
    claim_decision = {
        "numerical_convergence_gate_passed": convergence_passed,
        "independent_model_robustness_supported": robust,
        "analytical_E_optimal_time_days": e_optimal_time,
        "development_oracle_time_days": oracle_time,
        "active_learning_transport_status": heuristic_audit["status"],
        "active_learning_transport_claim_supported": heuristic_supported,
        "allowed_claim": (
            "The analytically selected E-optimal front observation improved recovery "
            "of both transport parameters when truth was generated by an independent "
            "converged numerical solver."
            if robust
            else "Independent-model robustness of the E-optimal transport result was "
            "not confirmed under the locked dual-parameter rule."
        ),
        "active_learning_consequence": (
            "The existing active-learning routine is topology/categorical decision "
            "support, not a transport experimental-design optimiser; the M8 transport "
            "performance claim is removed unless a Jacobian-aware interface is added."
        ),
        "guardrail": "Controlled-synthetic model-class evidence only; not field validation.",
    }

    development.to_csv(output / "oracle_development.csv", index=False, lineterminator="\n")
    test.to_csv(output / "locked_test_calibrations.csv", index=False, lineterminator="\n")
    summary.to_csv(output / "strategy_summary.csv", index=False, lineterminator="\n")
    contrasts.to_csv(output / "paired_bootstrap_contrasts.csv", index=False, lineterminator="\n")
    convergence.to_csv(output / "numerical_grid_convergence.csv", index=False, lineterminator="\n")
    candidate_scores.to_csv(output / "analytical_candidate_scores.csv", index=False, lineterminator="\n")
    fd_stability.to_csv(output / "analytical_fd_stability.csv", index=False, lineterminator="\n")
    (output / "active_learning_raw.json").write_text(
        json.dumps(heuristic_raw, indent=2), encoding="utf-8"
    )
    (output / "active_learning_portability.json").write_text(
        json.dumps(heuristic_audit, indent=2), encoding="utf-8"
    )
    (output / "claim_decision.json").write_text(
        json.dumps(claim_decision, indent=2), encoding="utf-8"
    )
    manifest = {
        "run_id": RUN_ID,
        "status": "PASS" if convergence_passed else "FAIL",
        "locked_protocol": locked,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "development_replicates": int(development_replicates),
        "test_replicates": int(test_replicates),
        "bootstrap_samples": int(bootstrap_samples),
        "truth_generator": "independent_implicit_finite_volume_upwind_v1",
        "calibration_model": "HydroSheaf analytical TransportCalibrationAdapter",
        "truth_cells": TRUTH_CELLS,
        "dt_days": DT_DAYS,
        "excluded_capabilities": ["temporal_series", "graph_3d", "vadose"],
        "artifacts": {},
    }
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            manifest["artifacts"][path.name] = _sha256(path)
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return {"manifest": manifest, "claim_decision": claim_decision}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path, default=PROJECT / "results" / RUN_ID
    )
    parser.add_argument("--development-replicates", type=int, default=DEVELOPMENT_REPLICATES)
    parser.add_argument("--test-replicates", type=int, default=TEST_REPLICATES)
    parser.add_argument("--bootstrap-samples", type=int, default=BOOTSTRAP_SAMPLES)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    development = 8 if args.quick else args.development_replicates
    test = 12 if args.quick else args.test_replicates
    bootstrap = 200 if args.quick else args.bootstrap_samples
    if min(development, test) < 2 or bootstrap < 100:
        raise ValueError("Insufficient replicate or bootstrap count.")
    result = run_benchmark(
        output=args.output,
        development_replicates=development,
        test_replicates=test,
        bootstrap_samples=bootstrap,
        overwrite=args.overwrite,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
