"""Two-stage prospectively locked M7.5 robust/hybrid sheaf diagnostic."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    f1_score,
    log_loss,
    roc_auc_score,
)

PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hydrosheaf.sheaf.directed_section import solve_directed_section  # noqa: E402

from independent_sheaf_graph_generator import (  # noqa: E402
    SCENARIOS,
    IndependentSectionCase,
    generate_independent_section_case,
    scenario_for_seed,
)
from run_m7_sheaf_vs_graph import (  # noqa: E402
    _edge_maps,
    _ece,
    _git_revision,
    _iteratively_reweighted_section,
    _select_threshold,
)

RUN_ID = "RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01"
DEVELOPMENT_SEEDS = tuple(range(8401, 8465))
LOCKED_TEST_SEEDS = tuple(range(8501, 8629))
QUICK_DEVELOPMENT_SEEDS = tuple(range(8701, 8717))
QUICK_TEST_SEEDS = tuple(range(8751, 8767))
BOOTSTRAP_SAMPLES = 10_000
LOO_ITERATIONS = 3
GROUP_FOLDS = 8
LAMBDA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
C_GRID = (0.1, 1.0, 10.0)

NATIVE_ARMS = (
    "edge_local",
    "identity_graph",
    "original_affine_global",
    "original_hybrid",
    "robust_affine_global",
    "robust_hybrid",
)
MODEL_ORDER = (*NATIVE_ARMS, "robust_hybrid_permuted")
CALIBRATION_REGIMES = ("separate_crossfit", "shared_calibrator")
HYBRID_ARMS = {
    "original_hybrid",
    "robust_hybrid",
    "robust_hybrid_permuted",
}

GLOBAL_COLUMN = {
    "identity_graph": "identity_section_residual",
    "original_affine_global": "original_affine_residual",
    "original_hybrid": "original_affine_residual",
    "robust_affine_global": "robust_loo_residual",
    "robust_hybrid": "robust_loo_residual",
    "robust_hybrid_permuted": "permuted_robust_loo_residual",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _seed_spec(spec: Mapping[str, object]) -> tuple[int, ...]:
    first = int(spec["first"])
    last = int(spec["last"])
    values = tuple(range(first, last + 1))
    if len(values) != int(spec["count"]):
        raise RuntimeError("Seed range and count disagree in protocol lock.")
    return values


def _verify_protocol_lock(
    development_seeds: Sequence[int],
    test_seeds: Sequence[int] | None = None,
) -> dict:
    path = PROJECT / "m7_robust_hybrid_protocol.lock.json"
    lock = json.loads(path.read_text(encoding="utf-8"))
    if lock.get("run_id") != RUN_ID:
        raise RuntimeError("M7.5 protocol lock has the wrong run ID.")
    if tuple(map(int, development_seeds)) != _seed_spec(lock["development_seeds"]):
        raise RuntimeError("Development seeds do not match the M7.5 protocol lock.")
    if test_seeds is not None and tuple(map(int, test_seeds)) != _seed_spec(
        lock["locked_test_seeds"]
    ):
        raise RuntimeError("Test seeds do not match the M7.5 protocol lock.")
    checks = {
        "bootstrap_samples": BOOTSTRAP_SAMPLES,
        "loo_iterations": LOO_ITERATIONS,
        "group_folds": GROUP_FOLDS,
        "hybrid_local_weight_grid": list(LAMBDA_GRID),
        "logistic_c_grid": list(C_GRID),
    }
    for key, actual in checks.items():
        if lock.get(key) != actual:
            raise RuntimeError(f"Frozen setting changed: {key}")
    for relative, expected in lock.get("sha256", {}).items():
        target = (REPO / relative).resolve()
        if not target.is_relative_to(REPO.resolve()):
            raise RuntimeError(f"Locked path escapes repository: {relative}")
        actual = _sha256(target)
        if actual != expected:
            raise RuntimeError(f"Protocol-bound file changed: {relative}")
    return lock


def _verify_confirmatory_lock(development_dir: Path) -> dict:
    lock_path = PROJECT / "m7_robust_hybrid_confirmatory.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("run_id") != RUN_ID or not lock.get("locked_before_test_generation"):
        raise RuntimeError("M7.5 confirmatory lock is absent or invalid.")
    targets = {
        "runner_sha256": Path(__file__).resolve(),
        "protocol_lock_sha256": PROJECT / "m7_robust_hybrid_protocol.lock.json",
        "frozen_settings_sha256": development_dir / "frozen_settings.json",
        "development_manifest_sha256": development_dir / "development_manifest.json",
    }
    for key, target in targets.items():
        if _sha256(target) != lock.get(key):
            raise RuntimeError(f"Confirmatory-bound artifact changed: {key}")
    return lock


def _safe_new_directory(path: Path, overwrite: bool) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(PROJECT.resolve()):
        raise ValueError("Output must remain inside the M7 project directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(f"Output exists: {resolved}")
        if resolved in {PROJECT.resolve(), (PROJECT / "results").resolve()}:
            raise ValueError("Refusing to remove a broad project directory.")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    return resolved


def _node_observations(case: IndependentSectionCase) -> dict[str, list[float] | None]:
    return {
        node_id: (None if value is None else [float(value)])
        for node_id, value in case.observations.items()
    }


def _held_out_residual(case: IndependentSectionCase, edge_map, maps) -> float:
    remaining = [candidate for candidate in maps if candidate.edge.edge_id != edge_map.edge.edge_id]
    section = solve_directed_section(
        case.node_ids,
        remaining,
        _node_observations(case),
        obs_weight=4.0,
        diag_eps=1.0e-8,
        non_negative=False,
    )
    x_u = float(section[edge_map.edge.u][0])
    x_v = float(section[edge_map.edge.v][0])
    return abs(float(edge_map.alpha) * x_u + float(edge_map.offset[0]) - x_v)


def _robust_loo_section(
    case: IndependentSectionCase,
    mode: str,
    *,
    permutation_seed: int | None = None,
) -> tuple[dict[str, float], dict[str, float], float]:
    maps = _edge_maps(case, mode, permutation_seed=permutation_seed)
    started = time.perf_counter()
    residuals: dict[str, float] = {}
    for _ in range(LOO_ITERATIONS):
        residuals = {
            edge_map.edge.edge_id: _held_out_residual(case, edge_map, maps)
            for edge_map in maps
        }
        scale = max(float(np.median(list(residuals.values()))), 1.0e-6)
        for edge_map in maps:
            base = float((edge_map.edge.attrs or {})["prior_edge_probability"])
            standardised = residuals[edge_map.edge.edge_id] / scale
            edge_map.weight = float(
                np.clip(base * np.exp(-0.5 * standardised * standardised), 1.0e-4, 1.0)
            )
    weights = {edge_map.edge.edge_id: float(edge_map.weight) for edge_map in maps}
    return residuals, weights, 1000.0 * (time.perf_counter() - started)


def _feature_frame(case: IndependentSectionCase) -> tuple[pd.DataFrame, dict]:
    identity, identity_ms = _iteratively_reweighted_section(case, "identity")
    original, original_ms = _iteratively_reweighted_section(case, "affine")
    robust, robust_weights, robust_ms = _robust_loo_section(case, "affine")
    permuted, permuted_weights, permuted_ms = _robust_loo_section(
        case, "permuted", permutation_seed=case.seed + 91_337
    )
    rows: list[dict] = []
    for candidate in case.candidates:
        obs_u = case.observations[candidate.u]
        obs_v = case.observations[candidate.v]
        missing = int(obs_u is None or obs_v is None)
        local = np.nan
        if not missing:
            local = abs(candidate.alpha * float(obs_u) + candidate.offset - float(obs_v))
        rows.append(
            {
                "seed": case.seed,
                "scenario": case.scenario,
                "edge_id": candidate.edge_id,
                "u": candidate.u,
                "v": candidate.v,
                "prior_probability": candidate.prior_probability,
                "prior_logit": float(
                    np.log(candidate.prior_probability / (1.0 - candidate.prior_probability))
                ),
                "local_missing": missing,
                "local_affine_residual": local,
                "identity_section_residual": identity[candidate.edge_id],
                "original_affine_residual": original[candidate.edge_id],
                "robust_loo_residual": robust[candidate.edge_id],
                "permuted_robust_loo_residual": permuted[candidate.edge_id],
                "robust_final_weight": robust_weights[candidate.edge_id],
                "permuted_robust_final_weight": permuted_weights[candidate.edge_id],
                "is_true_edge": candidate.is_true_edge,
                "is_corrupted_edge": candidate.is_corrupted_edge,
            }
        )
    identity_difference = np.nan
    if case.scenario == "identity_limit":
        identity_difference = max(abs(identity[key] - original[key]) for key in identity)
    diagnostics = {
        "seed": case.seed,
        "scenario": case.scenario,
        "identity_original_residual_max_abs_difference": float(identity_difference),
        "identity_runtime_ms": float(identity_ms),
        "original_affine_runtime_ms": float(original_ms),
        "robust_loo_runtime_ms": float(robust_ms),
        "permuted_robust_loo_runtime_ms": float(permuted_ms),
    }
    return pd.DataFrame(rows), diagnostics


def _representation(
    frame: pd.DataFrame,
    arm: str,
    local_weight: float | None,
    *,
    imputation: float | None = None,
) -> tuple[np.ndarray, float]:
    local = np.log1p(np.maximum(frame["local_affine_residual"].to_numpy(float), 0.0))
    if arm == "edge_local":
        if imputation is None:
            imputation = float(np.nanmedian(local))
        return np.where(np.isfinite(local), local, float(imputation)), float(imputation)
    global_value = np.log1p(
        np.maximum(frame[GLOBAL_COLUMN[arm]].to_numpy(float), 0.0)
    )
    if arm not in HYBRID_ARMS:
        return global_value, float(np.nanmedian(global_value))
    if local_weight is None:
        raise ValueError(f"Hybrid arm requires a local weight: {arm}")
    observed = np.isfinite(local)
    blended = global_value.copy()
    blended[observed] = (
        float(local_weight) * local[observed]
        + (1.0 - float(local_weight)) * global_value[observed]
    )
    return blended, float(np.nanmedian(blended))


def _matrix(
    frame: pd.DataFrame,
    arm: str,
    local_weight: float | None,
    state: Mapping[str, object] | None = None,
) -> tuple[np.ndarray, dict]:
    imputation = None if state is None else float(state["residual_imputation"])
    residual, residual_imputation = _representation(
        frame, arm, local_weight, imputation=imputation
    )
    raw = np.column_stack(
        [
            frame["prior_logit"].to_numpy(float),
            residual,
            frame["local_missing"].to_numpy(float),
        ]
    )
    if state is None:
        mean = raw.mean(axis=0)
        scale = raw.std(axis=0)
        scale[scale < 1.0e-12] = 1.0
    else:
        mean = np.asarray(state["mean"], dtype=float)
        scale = np.asarray(state["scale"], dtype=float)
    return (raw - mean) / scale, {
        "residual_imputation": float(residual_imputation),
        "mean": mean.astype(float).tolist(),
        "scale": scale.astype(float).tolist(),
    }


def _fit_logistic(matrix: np.ndarray, truth: np.ndarray, c_value: float) -> dict:
    model = LogisticRegression(
        C=float(c_value),
        class_weight=None,
        max_iter=2_000,
        random_state=20260729,
        solver="lbfgs",
    ).fit(matrix, truth)
    return {
        "intercept": float(model.intercept_[0]),
        "coefficients": model.coef_[0].astype(float).tolist(),
    }


def _predict_logistic(matrix: np.ndarray, model: Mapping[str, object]) -> np.ndarray:
    linear = float(model["intercept"]) + matrix @ np.asarray(
        model["coefficients"], dtype=float
    )
    probability = 1.0 / (1.0 + np.exp(-np.clip(linear, -40.0, 40.0)))
    return np.clip(probability, 1.0e-8, 1.0 - 1.0e-8)


def _fold_assignments(frame: pd.DataFrame) -> dict[int, int]:
    seeds = sorted(frame["seed"].astype(int).unique())
    return {seed: index % GROUP_FOLDS for index, seed in enumerate(seeds)}


def _crossfit(
    frame: pd.DataFrame,
    arm: str,
    local_weight: float | None,
    c_value: float,
) -> tuple[np.ndarray, float]:
    truth = frame["is_true_edge"].to_numpy(int)
    probability = np.full(len(frame), np.nan, dtype=float)
    assignments = _fold_assignments(frame)
    folds = frame["seed"].astype(int).map(assignments).to_numpy(int)
    for fold in range(GROUP_FOLDS):
        train_mask = folds != fold
        valid_mask = folds == fold
        train_matrix, state = _matrix(frame.loc[train_mask], arm, local_weight)
        model = _fit_logistic(train_matrix, truth[train_mask], c_value)
        valid_matrix, _ = _matrix(frame.loc[valid_mask], arm, local_weight, state)
        probability[valid_mask] = _predict_logistic(valid_matrix, model)
    if not np.isfinite(probability).all():
        raise RuntimeError(f"Cross-fitting failed for {arm}.")
    per_case = []
    for seed in sorted(frame["seed"].astype(int).unique()):
        mask = frame["seed"].to_numpy(int) == seed
        per_case.append(log_loss(truth[mask], probability[mask], labels=[0, 1]))
    return probability, float(np.mean(per_case))


def _fit_separate_models(
    development: pd.DataFrame,
) -> tuple[dict[str, dict], pd.DataFrame, pd.DataFrame]:
    truth = development["is_true_edge"].to_numpy(int)
    specs: dict[str, dict] = {}
    grid_rows: list[dict] = []
    oof_blocks: list[pd.DataFrame] = []
    for arm in NATIVE_ARMS:
        lambdas: Sequence[float | None] = LAMBDA_GRID if arm in HYBRID_ARMS else (None,)
        candidates = []
        for local_weight in lambdas:
            for c_value in C_GRID:
                probability, score = _crossfit(
                    development, arm, local_weight, float(c_value)
                )
                tie_lambda = 0.5 if local_weight is None else float(local_weight)
                row = {
                    "arm": arm,
                    "local_weight": local_weight,
                    "C": float(c_value),
                    "mean_case_log_loss": score,
                }
                grid_rows.append(row)
                candidates.append(
                    (
                        score,
                        abs(tie_lambda - 0.5),
                        tie_lambda,
                        abs(np.log10(float(c_value))),
                        float(c_value),
                        probability,
                        local_weight,
                    )
                )
        selected = min(candidates, key=lambda item: item[:5])
        score, _, _, _, c_value, probability, local_weight = selected
        matrix, state = _matrix(development, arm, local_weight)
        model = _fit_logistic(matrix, truth, c_value)
        threshold = _select_threshold(truth, probability)
        specs[arm] = {
            "arm": arm,
            "local_weight": local_weight,
            "C": float(c_value),
            "selection_mean_case_log_loss": float(score),
            "threshold": float(threshold),
            "state": state,
            "model": model,
        }
        block = development[["seed", "scenario", "edge_id", "is_true_edge"]].copy()
        block["calibration_regime"] = "separate_crossfit"
        block["model"] = arm
        block["probability"] = probability
        oof_blocks.append(block)
    return specs, pd.DataFrame(grid_rows), pd.concat(oof_blocks, ignore_index=True)


def _fit_shared_model(
    development: pd.DataFrame,
    separate_specs: Mapping[str, Mapping[str, object]],
) -> tuple[dict, pd.DataFrame]:
    truth = development["is_true_edge"].to_numpy(int)
    assignments = _fold_assignments(development)
    folds = development["seed"].astype(int).map(assignments).to_numpy(int)
    oof_by_arm = {arm: np.full(len(development), np.nan) for arm in NATIVE_ARMS}
    for fold in range(GROUP_FOLDS):
        train_mask = folds != fold
        valid_mask = folds == fold
        matrices = []
        states = {}
        for arm in NATIVE_ARMS:
            local_weight = separate_specs[arm]["local_weight"]
            matrix, state = _matrix(development.loc[train_mask], arm, local_weight)
            matrices.append(matrix)
            states[arm] = state
        pooled_matrix = np.vstack(matrices)
        pooled_truth = np.tile(truth[train_mask], len(NATIVE_ARMS))
        model = _fit_logistic(pooled_matrix, pooled_truth, 1.0)
        for arm in NATIVE_ARMS:
            valid_matrix, _ = _matrix(
                development.loc[valid_mask],
                arm,
                separate_specs[arm]["local_weight"],
                states[arm],
            )
            oof_by_arm[arm][valid_mask] = _predict_logistic(valid_matrix, model)
    pooled_oof = np.concatenate([oof_by_arm[arm] for arm in NATIVE_ARMS])
    pooled_truth = np.tile(truth, len(NATIVE_ARMS))
    threshold = _select_threshold(pooled_truth, pooled_oof)
    states = {}
    matrices = []
    for arm in NATIVE_ARMS:
        matrix, state = _matrix(
            development, arm, separate_specs[arm]["local_weight"]
        )
        matrices.append(matrix)
        states[arm] = state
    model = _fit_logistic(np.vstack(matrices), pooled_truth, 1.0)
    blocks = []
    for arm in NATIVE_ARMS:
        block = development[["seed", "scenario", "edge_id", "is_true_edge"]].copy()
        block["calibration_regime"] = "shared_calibrator"
        block["model"] = arm
        block["probability"] = oof_by_arm[arm]
        blocks.append(block)
    spec = {
        "C": 1.0,
        "threshold": float(threshold),
        "states": states,
        "local_weights": {
            arm: separate_specs[arm]["local_weight"] for arm in NATIVE_ARMS
        },
        "model": model,
    }
    return spec, pd.concat(blocks, ignore_index=True)


def _generate_features(seeds: Sequence[int], split: str) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    first_seed = int(min(seeds))
    frames: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    provenance: list[dict] = []
    for index, seed in enumerate(seeds, start=1):
        scenario = scenario_for_seed(int(seed), first_seed)
        case = generate_independent_section_case(int(seed), scenario)
        frame, case_diagnostics = _feature_frame(case)
        frame["split"] = split
        frames.append(frame)
        diagnostics.append(case_diagnostics | {"split": split})
        provenance.append(dict(case.provenance) | {"split": split})
        if index % 8 == 0:
            print(f"{split}: generated {index}/{len(seeds)} cases", flush=True)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(diagnostics), provenance


def _stable_csv(frame: pd.DataFrame, path: Path) -> None:
    stable = frame.copy()
    float_columns = stable.select_dtypes(include=["float32", "float64"]).columns
    stable[float_columns] = stable[float_columns].round(12)
    stable.to_csv(path, index=False, lineterminator="\n")


def run_development(output: Path, *, overwrite: bool = False, quick: bool = False) -> dict:
    seeds = QUICK_DEVELOPMENT_SEEDS if quick else DEVELOPMENT_SEEDS
    if not quick:
        _verify_protocol_lock(seeds)
    output = _safe_new_directory(output, overwrite)
    features, diagnostics, provenance = _generate_features(seeds, "development")
    separate, grid, separate_oof = _fit_separate_models(features)
    shared, shared_oof = _fit_shared_model(features, separate)
    settings = {
        "schema_version": "1.0",
        "run_id": RUN_ID,
        "development_seeds": list(map(int, seeds)),
        "separate_crossfit": separate,
        "shared_calibrator": shared,
        "feature_order": ["prior_logit", "residual_scalar", "local_missing"],
        "truth_bearing_inference_features": [],
    }
    _stable_csv(features, output / "development_edge_features.csv")
    _stable_csv(diagnostics, output / "development_case_diagnostics.csv")
    _stable_csv(grid, output / "selection_grid.csv")
    _stable_csv(
        pd.concat([separate_oof, shared_oof], ignore_index=True),
        output / "development_oof_predictions.csv",
    )
    (output / "frozen_settings.json").write_text(
        json.dumps(settings, indent=2), encoding="utf-8"
    )
    (output / "generator_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    manifest = {
        "schema_version": "1.0",
        "run_id": RUN_ID,
        "stage": "development",
        "status": "PASS",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "seeds": list(map(int, seeds)),
        "locked_test_generated": False,
        "artifacts": {},
    }
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "development_manifest.json":
            manifest["artifacts"][path.name] = _sha256(path)
    (output / "development_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest


def _predict_models(features: pd.DataFrame, settings: Mapping[str, object]) -> pd.DataFrame:
    blocks: list[pd.DataFrame] = []
    separate = settings["separate_crossfit"]
    for arm in MODEL_ORDER:
        source_arm = "robust_hybrid" if arm == "robust_hybrid_permuted" else arm
        spec = separate[source_arm]
        matrix, _ = _matrix(
            features,
            arm,
            spec["local_weight"],
            spec["state"],
        )
        probability = _predict_logistic(matrix, spec["model"])
        block = features[
            ["seed", "scenario", "edge_id", "u", "v", "is_true_edge", "is_corrupted_edge"]
        ].copy()
        block["calibration_regime"] = "separate_crossfit"
        block["model"] = arm
        block["probability"] = probability
        block["selected"] = (probability >= float(spec["threshold"])).astype(int)
        blocks.append(block)
    shared = settings["shared_calibrator"]
    for arm in MODEL_ORDER:
        source_arm = "robust_hybrid" if arm == "robust_hybrid_permuted" else arm
        matrix, _ = _matrix(
            features,
            arm,
            shared["local_weights"][source_arm],
            shared["states"][source_arm],
        )
        probability = _predict_logistic(matrix, shared["model"])
        block = features[
            ["seed", "scenario", "edge_id", "u", "v", "is_true_edge", "is_corrupted_edge"]
        ].copy()
        block["calibration_regime"] = "shared_calibrator"
        block["model"] = arm
        block["probability"] = probability
        block["selected"] = (probability >= float(shared["threshold"])).astype(int)
        blocks.append(block)
    return pd.concat(blocks, ignore_index=True)


def _case_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (seed, scenario, regime, model), group in predictions.groupby(
        ["seed", "scenario", "calibration_regime", "model"], sort=True
    ):
        truth = group["is_true_edge"].to_numpy(int)
        probability = group["probability"].to_numpy(float)
        selected = group["selected"].to_numpy(int)
        false_mask = truth == 0
        confident = (probability <= 0.20) | (probability >= 0.80)
        corrupted = group["is_corrupted_edge"].to_numpy(int)
        conflict_ap = np.nan
        if corrupted.sum() > 0 and corrupted.sum() < len(corrupted):
            conflict_ap = average_precision_score(corrupted, 1.0 - probability)
        rows.append(
            {
                "seed": int(seed),
                "scenario": scenario,
                "calibration_regime": regime,
                "model": model,
                "pr_auc": average_precision_score(truth, probability),
                "roc_auc": roc_auc_score(truth, probability),
                "brier": brier_score_loss(truth, probability),
                "log_loss": log_loss(truth, probability, labels=[0, 1]),
                "ece": _ece(truth, probability),
                "selected_f1": f1_score(truth, selected, zero_division=0),
                "false_confidence_rate": float(
                    (probability[false_mask] >= 0.80).sum() / max(false_mask.sum(), 1)
                ),
                "abstention_coverage": float(confident.mean()),
                "abstention_accuracy": float(
                    ((probability[confident] >= 0.5) == truth[confident]).mean()
                    if confident.any()
                    else np.nan
                ),
                "conflict_localisation_pr_auc": float(conflict_ap),
            }
        )
    return pd.DataFrame(rows)


def _paired_bootstrap(
    metrics: pd.DataFrame,
    regime: str,
    left: str,
    right: str,
    metric: str,
    scenario: str,
    samples: int,
    seed: int,
) -> dict:
    subset = metrics[metrics["calibration_regime"] == regime]
    if scenario != "all":
        subset = subset[subset["scenario"] == scenario]
    pivot = subset.pivot(index="seed", columns="model", values=metric).dropna()
    differences = (pivot[left] - pivot[right]).to_numpy(float)
    if len(differences) == 0:
        return {
            "scenario": scenario,
            "calibration_regime": regime,
            "left": left,
            "right": right,
            "metric": metric,
            "n_cases": 0,
            "mean_difference": np.nan,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
        }
    rng = np.random.default_rng(int(seed))
    indices = rng.integers(0, len(differences), size=(samples, len(differences)))
    boot = differences[indices].mean(axis=1)
    return {
        "scenario": scenario,
        "calibration_regime": regime,
        "left": left,
        "right": right,
        "metric": metric,
        "n_cases": len(differences),
        "mean_difference": float(differences.mean()),
        "ci95_low": float(np.quantile(boot, 0.025)),
        "ci95_high": float(np.quantile(boot, 0.975)),
    }


def _contrast_table(metrics: pd.DataFrame, samples: int) -> pd.DataFrame:
    comparisons = (
        ("robust_hybrid", "edge_local"),
        ("robust_hybrid", "robust_hybrid_permuted"),
        ("original_hybrid", "original_affine_global"),
        ("robust_affine_global", "original_affine_global"),
        ("robust_hybrid", "original_hybrid"),
        ("robust_hybrid", "robust_affine_global"),
        ("original_affine_global", "edge_local"),
    )
    metric_names = (
        "pr_auc",
        "brier",
        "log_loss",
        "ece",
        "selected_f1",
        "false_confidence_rate",
        "abstention_accuracy",
        "conflict_localisation_pr_auc",
    )
    rows = []
    counter = 0
    for regime in CALIBRATION_REGIMES:
        for scenario in ("all", *SCENARIOS):
            for left, right in comparisons:
                for metric in metric_names:
                    counter += 1
                    rows.append(
                        _paired_bootstrap(
                            metrics,
                            regime,
                            left,
                            right,
                            metric,
                            scenario,
                            samples,
                            2026072900 + counter,
                        )
                    )
    return pd.DataFrame(rows)


def _contrast(
    contrasts: pd.DataFrame,
    left: str,
    right: str,
    metric: str,
    scenario: str = "all",
    regime: str = "separate_crossfit",
) -> dict:
    row = contrasts[
        (contrasts["calibration_regime"] == regime)
        & (contrasts["scenario"] == scenario)
        & (contrasts["left"] == left)
        & (contrasts["right"] == right)
        & (contrasts["metric"] == metric)
    ]
    if len(row) != 1:
        raise KeyError((regime, scenario, left, right, metric))
    return row.iloc[0].to_dict()


def _claim_decision(
    metrics: pd.DataFrame,
    contrasts: pd.DataFrame,
    diagnostics: pd.DataFrame,
    provenance: Sequence[Mapping[str, object]],
) -> dict:
    finite = bool(np.isfinite(metrics[["pr_auc", "brier", "log_loss"]]).all().all())
    independent = all(item.get("imports_hydrosheaf") is False for item in provenance)
    identity_raw = float(
        diagnostics["identity_original_residual_max_abs_difference"].dropna().max()
    )
    execution = bool(finite and independent and identity_raw <= 1.0e-10)
    primary = {
        metric: _contrast(contrasts, "robust_hybrid", "edge_local", metric)
        for metric in ("pr_auc", "brier", "log_loss")
    }
    identity = {
        metric: _contrast(
            contrasts,
            "robust_hybrid",
            "edge_local",
            metric,
            scenario="identity_limit",
        )
        for metric in ("pr_auc", "brier", "log_loss")
    }
    permuted = {
        metric: _contrast(
            contrasts, "robust_hybrid", "robust_hybrid_permuted", metric
        )
        for metric in ("pr_auc", "brier", "log_loss")
    }
    incremental = bool(
        execution
        and primary["pr_auc"]["ci95_low"] > 0.0
        and primary["brier"]["ci95_high"] < 0.0
        and primary["log_loss"]["ci95_high"] < 0.0
        and identity["pr_auc"]["ci95_low"] >= -0.02
        and identity["brier"]["ci95_high"] <= 0.01
        and identity["log_loss"]["ci95_high"] <= 0.02
        and permuted["pr_auc"]["ci95_low"] > 0.0
        and permuted["brier"]["ci95_high"] < 0.0
    )
    mechanism = {
        "local_evidence_restoration": {
            metric: _contrast(
                contrasts, "original_hybrid", "original_affine_global", metric
            )
            for metric in ("pr_auc", "brier", "log_loss")
        },
        "loo_robustness": {
            metric: _contrast(
                contrasts, "robust_affine_global", "original_affine_global", metric
            )
            for metric in ("pr_auc", "brier", "log_loss")
        },
        "robust_hybrid_increment": {
            metric: _contrast(contrasts, "robust_hybrid", "original_hybrid", metric)
            for metric in ("pr_auc", "brier", "log_loss")
        },
    }
    favourable_scenarios = []
    for scenario in SCENARIOS:
        pr = _contrast(
            contrasts, "robust_hybrid", "edge_local", "pr_auc", scenario=scenario
        )
        ll = _contrast(
            contrasts, "robust_hybrid", "edge_local", "log_loss", scenario=scenario
        )
        if pr["ci95_low"] > 0.0 or ll["ci95_high"] < 0.0:
            favourable_scenarios.append(
                {"scenario": scenario, "pr_auc": pr, "log_loss": ll}
            )
    if incremental:
        allowed_claim = (
            "In this prospectively locked controlled-synthetic benchmark, the "
            "LOO robust hybrid sheaf supplied incremental discrimination and "
            "calibration value beyond the edge-local weighted graph and retained "
            "native-map information beyond a permuted adverse control."
        )
    elif favourable_scenarios:
        names = ", ".join(item["scenario"] for item in favourable_scenarios)
        allowed_claim = (
            "The full incremental robust-hybrid gate failed. Conditional improvements "
            f"were confined to the following prespecified strata/metrics: {names}."
        )
    else:
        allowed_claim = (
            "The full incremental robust-hybrid gate failed and no favourable "
            "scenario-specific interval excluded zero; the tested global information "
            "did not add reliable value beyond the edge-local graph in this generator."
        )
    return {
        "execution_provenance_gate_passed": execution,
        "scientific_workflow_status": "PASS" if execution else "FAIL",
        "incremental_robust_hybrid_sheaf_claim_allowed": incremental,
        "identity_original_residual_max_abs_difference": identity_raw,
        "primary_robust_hybrid_minus_edge_local": primary,
        "identity_limit_non_degradation": identity,
        "native_minus_permuted_control": permuted,
        "mechanism_diagnostics": mechanism,
        "favourable_scenario_diagnostics": favourable_scenarios,
        "allowed_claim": allowed_claim,
        "guardrail": (
            "Controlled-synthetic scalar, static, two-dimensional evidence only; not "
            "field validation or evidence for temporal, three-dimensional, vadose-zone, "
            "or general HydroSheaf superiority."
        ),
    }


def _write_publication_assets(
    metrics: pd.DataFrame, contrasts: pd.DataFrame
) -> list[Path]:
    figure_dir = PROJECT / "figures" / "publication"
    table_dir = PROJECT / "tables" / "publication"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    primary = metrics[metrics["calibration_regime"] == "separate_crossfit"]
    summary = (
        primary.groupby("model", sort=False)
        .agg(
            n_cases=("seed", "count"),
            pr_auc_mean=("pr_auc", "mean"),
            pr_auc_sd=("pr_auc", "std"),
            brier_mean=("brier", "mean"),
            log_loss_mean=("log_loss", "mean"),
            ece_mean=("ece", "mean"),
            selected_f1_mean=("selected_f1", "mean"),
        )
        .reindex(MODEL_ORDER)
        .reset_index()
    )
    table_csv = table_dir / "table9_robust_hybrid_sheaf.csv"
    table_md = table_dir / "table9_robust_hybrid_sheaf.md"
    summary.round(6).to_csv(table_csv, index=False, lineterminator="\n")
    table_md.write_text(summary.round(4).to_markdown(index=False) + "\n", encoding="utf-8")

    labels = {
        "edge_local": "Edge-local",
        "identity_graph": "Identity graph",
        "original_affine_global": "Original global",
        "original_hybrid": "Original hybrid",
        "robust_affine_global": "Robust global",
        "robust_hybrid": "Robust hybrid",
        "robust_hybrid_permuted": "Permuted control",
    }
    colours = ["#777777", "#4c78a8", "#72a0c1", "#59a14f", "#2a9d8f", "#006d77", "#e76f51"]
    overall = primary.groupby("model", sort=False).mean(numeric_only=True).reindex(MODEL_ORDER)
    x = np.arange(len(MODEL_ORDER))
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.0), constrained_layout=True)
    axes[0, 0].bar(x, overall["pr_auc"], color=colours, edgecolor="black", linewidth=0.4)
    axes[0, 0].set_xticks(x, [labels[m] for m in MODEL_ORDER], rotation=28, ha="right")
    axes[0, 0].set_ylabel("Mean case PR-AUC")
    axes[0, 0].set_title("a  Locked-test discrimination")

    axes[0, 1].bar(x, overall["log_loss"], color=colours, edgecolor="black", linewidth=0.4)
    axes[0, 1].set_xticks(x, [labels[m] for m in MODEL_ORDER], rotation=28, ha="right")
    axes[0, 1].set_ylabel("Mean log loss (lower is better)")
    axes[0, 1].set_title("b  Locked-test calibration")

    scenario_rows = contrasts[
        (contrasts["calibration_regime"] == "separate_crossfit")
        & (contrasts["left"] == "robust_hybrid")
        & (contrasts["right"] == "edge_local")
        & (contrasts["metric"] == "pr_auc")
        & (contrasts["scenario"] != "all")
    ].set_index("scenario").reindex(SCENARIOS)
    means = scenario_rows["mean_difference"].to_numpy(float)
    errors = np.vstack(
        [means - scenario_rows["ci95_low"].to_numpy(float), scenario_rows["ci95_high"].to_numpy(float) - means]
    )
    axes[1, 0].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 0].errorbar(means, np.arange(len(SCENARIOS)), xerr=errors, fmt="o", color="#006d77", capsize=3)
    axes[1, 0].set_yticks(np.arange(len(SCENARIOS)), [s.replace("_", " ") for s in SCENARIOS])
    axes[1, 0].set_xlabel("PR-AUC difference (robust hybrid - edge local)")
    axes[1, 0].set_title("c  Prespecified scenario diagnostics")

    mechanisms = [
        ("Original hybrid - global", "original_hybrid", "original_affine_global"),
        ("Robust global - original", "robust_affine_global", "original_affine_global"),
        ("Robust hybrid - original", "robust_hybrid", "original_hybrid"),
        ("Native - permuted", "robust_hybrid", "robust_hybrid_permuted"),
    ]
    rows = [
        _contrast(contrasts, left, right, "pr_auc") | {"label": label}
        for label, left, right in mechanisms
    ]
    mechanism_frame = pd.DataFrame(rows)
    means = mechanism_frame["mean_difference"].to_numpy(float)
    errors = np.vstack(
        [means - mechanism_frame["ci95_low"].to_numpy(float), mechanism_frame["ci95_high"].to_numpy(float) - means]
    )
    axes[1, 1].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 1].errorbar(means, np.arange(len(rows)), xerr=errors, fmt="o", color="#2a9d8f", capsize=3)
    axes[1, 1].set_yticks(np.arange(len(rows)), mechanism_frame["label"])
    axes[1, 1].set_xlabel("Paired PR-AUC difference")
    axes[1, 1].set_title("d  Mechanism attribution")
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.16)
    fig.suptitle("M7.5: robust and hybrid sheaf estimator diagnostic", fontsize=14, fontweight="bold")
    figure_png = figure_dir / "figure7_robust_hybrid_sheaf.png"
    figure_pdf = figure_dir / "figure7_robust_hybrid_sheaf.pdf"
    figure_tif = figure_dir / "figure7_robust_hybrid_sheaf.tif"
    fig.savefig(figure_png, dpi=600, bbox_inches="tight")
    fig.savefig(figure_pdf, bbox_inches="tight")
    fig.savefig(figure_tif, dpi=300, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")
    plt.close(fig)
    return [table_csv, table_md, figure_png, figure_pdf, figure_tif]


def run_confirmatory(
    output: Path,
    development_dir: Path,
    *,
    overwrite: bool = False,
    quick: bool = False,
) -> dict:
    seeds = QUICK_TEST_SEEDS if quick else LOCKED_TEST_SEEDS
    if not quick:
        _verify_protocol_lock(DEVELOPMENT_SEEDS, seeds)
        _verify_confirmatory_lock(development_dir)
    output = _safe_new_directory(output, overwrite)
    settings = json.loads((development_dir / "frozen_settings.json").read_text(encoding="utf-8"))
    features, diagnostics, provenance = _generate_features(seeds, "locked_test")
    predictions = _predict_models(features, settings)
    metrics = _case_metrics(predictions)
    contrasts = _contrast_table(metrics, BOOTSTRAP_SAMPLES if not quick else 500)
    claim = _claim_decision(metrics, contrasts, diagnostics, provenance)
    summary = (
        metrics.groupby(["calibration_regime", "model"], sort=False)
        .mean(numeric_only=True)
        .reset_index()
    )
    scenario_summary = (
        metrics.groupby(["scenario", "calibration_regime", "model"], sort=False)
        .mean(numeric_only=True)
        .reset_index()
    )
    for filename, frame in {
        "locked_test_edge_features.csv": features,
        "locked_test_case_diagnostics.csv": diagnostics,
        "edge_predictions.csv": predictions,
        "case_metrics.csv": metrics,
        "model_summary.csv": summary,
        "scenario_summary.csv": scenario_summary,
        "paired_bootstrap_contrasts.csv": contrasts,
    }.items():
        _stable_csv(frame, output / filename)
    (output / "generator_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    (output / "claim_decision.json").write_text(
        json.dumps(claim, indent=2), encoding="utf-8"
    )
    assets = [] if quick else _write_publication_assets(metrics, contrasts)
    manifest = {
        "schema_version": "1.0",
        "run_id": RUN_ID,
        "stage": "locked_test",
        "status": "PASS" if claim["execution_provenance_gate_passed"] else "FAIL",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "development_seeds": list(map(int, DEVELOPMENT_SEEDS if not quick else QUICK_DEVELOPMENT_SEEDS)),
        "locked_test_seeds": list(map(int, seeds)),
        "locked_test_execution_count": 1,
        "bootstrap_samples": BOOTSTRAP_SAMPLES if not quick else 500,
        "models": list(MODEL_ORDER),
        "calibration_regimes": list(CALIBRATION_REGIMES),
        "locked_test_truth_joined_only_for_scoring": True,
        "excluded_capabilities": ["temporal_series", "graph_3d", "vadose"],
        "publication_assets": [str(path.relative_to(PROJECT)) for path in assets],
        "artifacts": {},
    }
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            manifest["artifacts"][path.name] = _sha256(path)
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return {"manifest": manifest, "claim_decision": claim}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("development", "confirmatory"), required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--development-dir",
        type=Path,
        default=PROJECT / "results" / RUN_ID / "development",
    )
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = PROJECT / "results" / RUN_ID
    if args.stage == "development":
        output = args.output or root / ("quick_development" if args.quick else "development")
        result = run_development(output, overwrite=args.overwrite, quick=args.quick)
    else:
        output = args.output or root / ("quick_test" if args.quick else "locked_test")
        result = run_confirmatory(
            output,
            args.development_dir,
            overwrite=args.overwrite,
            quick=args.quick,
        )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
