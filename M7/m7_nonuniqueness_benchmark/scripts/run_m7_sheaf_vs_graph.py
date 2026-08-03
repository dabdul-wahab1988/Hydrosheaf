"""Prospectively locked M7.4 sheaf-versus-weighted-graph benchmark."""

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

from hydrosheaf.graph.types import Edge  # noqa: E402
from hydrosheaf.sheaf.directed_section import (  # noqa: E402
    DirectedEdgeMap,
    compute_edge_section_residuals,
    solve_directed_section,
)

from independent_sheaf_graph_generator import (  # noqa: E402
    SCENARIOS,
    IndependentSectionCase,
    generate_independent_section_case,
    scenario_for_seed,
)

RUN_ID = "RUN-M7-SHEAF-VS-GRAPH-20260729-01"
DEVELOPMENT_SEEDS = tuple(range(7401, 7433))
LOCKED_TEST_SEEDS = tuple(range(7501, 7565))
QUICK_DEVELOPMENT_SEEDS = tuple(range(7901, 7909))
QUICK_TEST_SEEDS = tuple(range(7951, 7959))
BOOTSTRAP_SAMPLES = 10_000
MODEL_ORDER = (
    "weighted_graph",
    "graph_laplacian",
    "affine_sheaf",
    "affine_sheaf_permuted",
)
PRIMARY_MODELS = MODEL_ORDER[:3]
RESIDUAL_COLUMNS = {
    "weighted_graph": "local_affine_residual",
    "graph_laplacian": "identity_section_residual",
    "affine_sheaf": "affine_section_residual",
    "affine_sheaf_permuted": "permuted_section_residual",
}


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
    development_seeds: Sequence[int],
    test_seeds: Sequence[int],
    bootstrap_samples: int,
) -> dict:
    lock_path = PROJECT / "m7_sheaf_vs_graph_protocol.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock.get("run_id") != RUN_ID:
        raise RuntimeError("M7.4 lock has the wrong run ID.")
    if list(map(int, development_seeds)) != lock.get("development_seeds"):
        raise RuntimeError("Requested development seeds do not match the M7.4 lock.")
    if list(map(int, test_seeds)) != lock.get("locked_test_seeds"):
        raise RuntimeError("Requested test seeds do not match the M7.4 lock.")
    if int(bootstrap_samples) != int(lock.get("bootstrap_samples", -1)):
        raise RuntimeError("Bootstrap count does not match the M7.4 lock.")
    for relative, expected in lock.get("sha256", {}).items():
        path = (REPO / relative).resolve()
        if not path.is_relative_to(REPO.resolve()):
            raise RuntimeError(f"Locked path escapes repository: {relative}")
        actual = _sha256(path)
        if actual != expected:
            raise RuntimeError(f"Locked file changed: {relative}: {actual} != {expected}")
    return lock


def _safe_output(path: Path, overwrite: bool) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(PROJECT.resolve()):
        raise ValueError("Output must stay within the M7 benchmark directory.")
    if resolved.exists():
        if not overwrite:
            raise FileExistsError(f"Output exists: {resolved}")
        if resolved == PROJECT.resolve():
            raise ValueError("Refusing to remove the M7 benchmark root.")
        shutil.rmtree(resolved)
    resolved.mkdir(parents=True)
    return resolved


def _edge_maps(
    case: IndependentSectionCase,
    mode: str,
    *,
    permutation_seed: int | None = None,
) -> list[DirectedEdgeMap]:
    candidates = list(case.candidates)
    params = [(candidate.alpha, candidate.offset) for candidate in candidates]
    if mode == "identity":
        params = [(1.0, 0.0)] * len(candidates)
    elif mode == "permuted":
        rng = np.random.default_rng(int(permutation_seed))
        params = [params[index] for index in rng.permutation(len(params))]
    elif mode != "affine":
        raise ValueError(f"Unknown map mode: {mode}")

    maps: list[DirectedEdgeMap] = []
    for candidate, (alpha, offset) in zip(candidates, params):
        edge = Edge(
            edge_id=candidate.edge_id,
            u=candidate.u,
            v=candidate.v,
            attrs={"prior_edge_probability": candidate.prior_probability},
        )
        maps.append(
            DirectedEdgeMap(
                edge=edge,
                alpha=float(alpha),
                offset=[float(offset)],
                weight=float(candidate.prior_probability),
                objective=0.0,
                transport_model=mode,
                endmember_id=None,
                residual_norm=0.0,
            )
        )
    return maps


def _iteratively_reweighted_section(
    case: IndependentSectionCase,
    mode: str,
    *,
    permutation_seed: int | None = None,
    iterations: int = 3,
) -> tuple[dict[str, float], float]:
    maps = _edge_maps(case, mode, permutation_seed=permutation_seed)
    node_obs = {
        node_id: (None if value is None else [float(value)])
        for node_id, value in case.observations.items()
    }
    started = time.perf_counter()
    residuals: dict[str, float] = {}
    for _ in range(int(iterations)):
        section = solve_directed_section(
            case.node_ids,
            maps,
            node_obs,
            obs_weight=4.0,
            diag_eps=1.0e-8,
            non_negative=False,
        )
        residuals = compute_edge_section_residuals(
            {edge_map.edge.edge_id: edge_map for edge_map in maps},
            section,
            [1.0],
        )
        scale = max(float(np.median(list(residuals.values()))), 1.0e-6)
        for edge_map in maps:
            base = float((edge_map.edge.attrs or {})["prior_edge_probability"])
            edge_map.weight = float(
                np.clip(base * np.exp(-residuals[edge_map.edge.edge_id] / scale), 1e-4, 1.0)
            )
    elapsed_ms = 1000.0 * (time.perf_counter() - started)
    return residuals, float(elapsed_ms)


def _feature_frame(case: IndependentSectionCase) -> tuple[pd.DataFrame, dict]:
    identity, graph_ms = _iteratively_reweighted_section(case, "identity")
    affine, sheaf_ms = _iteratively_reweighted_section(case, "affine")
    permuted, permuted_ms = _iteratively_reweighted_section(
        case, "permuted", permutation_seed=case.seed + 91_337
    )
    rows: list[dict] = []
    for candidate in case.candidates:
        obs_u = case.observations[candidate.u]
        obs_v = case.observations[candidate.v]
        local_missing = int(obs_u is None or obs_v is None)
        local_residual = np.nan
        if not local_missing:
            local_residual = abs(
                candidate.alpha * float(obs_u) + candidate.offset - float(obs_v)
            )
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
                "local_missing": local_missing,
                "local_affine_residual": local_residual,
                "identity_section_residual": identity[candidate.edge_id],
                "affine_section_residual": affine[candidate.edge_id],
                "permuted_section_residual": permuted[candidate.edge_id],
                "is_true_edge": candidate.is_true_edge,
                "is_corrupted_edge": candidate.is_corrupted_edge,
            }
        )
    diagnostics = {
        "seed": case.seed,
        "scenario": case.scenario,
        "identity_affine_residual_max_abs_difference": float(
            max(abs(identity[key] - affine[key]) for key in identity)
            if case.scenario == "identity_limit"
            else np.nan
        ),
        "graph_laplacian_runtime_ms": graph_ms,
        "affine_sheaf_runtime_ms": sheaf_ms,
        "permuted_sheaf_runtime_ms": permuted_ms,
    }
    return pd.DataFrame(rows), diagnostics


def _model_matrix(
    frame: pd.DataFrame,
    residual_column: str,
    *,
    median: float | None = None,
    mean: Sequence[float] | None = None,
    scale: Sequence[float] | None = None,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    residual = frame[residual_column].to_numpy(float)
    if median is None:
        median = float(np.nanmedian(residual))
    residual = np.where(np.isfinite(residual), residual, float(median))
    raw = np.column_stack(
        [
            frame["prior_logit"].to_numpy(float),
            np.log1p(np.maximum(residual, 0.0)),
            frame["local_missing"].to_numpy(float),
        ]
    )
    if mean is None:
        mean_array = raw.mean(axis=0)
    else:
        mean_array = np.asarray(mean, dtype=float)
    if scale is None:
        scale_array = raw.std(axis=0)
        scale_array[scale_array < 1.0e-12] = 1.0
    else:
        scale_array = np.asarray(scale, dtype=float)
    return (raw - mean_array) / scale_array, float(median), mean_array, scale_array


def _select_threshold(truth: np.ndarray, probability: np.ndarray) -> float:
    candidates = np.linspace(0.10, 0.90, 81)
    scores = np.asarray(
        [f1_score(truth, probability >= threshold, zero_division=0) for threshold in candidates]
    )
    best = np.flatnonzero(scores == scores.max())
    return float(candidates[best[np.argmin(abs(candidates[best] - 0.5))]])


def _fit_models(development: pd.DataFrame) -> dict[str, dict]:
    truth = development["is_true_edge"].to_numpy(int)
    fitted: dict[str, dict] = {}
    # The local comparator receives its own equal-budget calibration.  The two
    # global representations share one calibrator fitted to the stacked
    # identity- and affine-residual development rows.  This makes predictions
    # identical when restrictions are identical and isolates representation
    # rather than calibration as the only graph-versus-sheaf difference.
    matrix, median, mean, scale = _model_matrix(
        development, RESIDUAL_COLUMNS["weighted_graph"]
    )
    local_model = LogisticRegression(
        C=1.0,
        class_weight=None,
        max_iter=2_000,
        random_state=20260729,
        solver="lbfgs",
    ).fit(matrix, truth)
    local_probability = local_model.predict_proba(matrix)[:, 1]
    fitted["weighted_graph"] = {
        "model": local_model,
        "residual_column": RESIDUAL_COLUMNS["weighted_graph"],
        "residual_median": median,
        "mean": mean,
        "scale": scale,
        "threshold": _select_threshold(truth, local_probability),
    }

    identity = development.copy()
    identity["global_section_residual"] = identity["identity_section_residual"]
    affine = development.copy()
    affine["global_section_residual"] = affine["affine_section_residual"]
    pooled = pd.concat([identity, affine], ignore_index=True)
    pooled_truth = pooled["is_true_edge"].to_numpy(int)
    matrix, median, mean, scale = _model_matrix(pooled, "global_section_residual")
    global_model = LogisticRegression(
        C=1.0,
        class_weight=None,
        max_iter=2_000,
        random_state=20260729,
        solver="lbfgs",
    ).fit(matrix, pooled_truth)
    pooled_probability = global_model.predict_proba(matrix)[:, 1]
    global_threshold = _select_threshold(pooled_truth, pooled_probability)
    for model_name in ("graph_laplacian", "affine_sheaf"):
        fitted[model_name] = {
            "model": global_model,
            "residual_column": RESIDUAL_COLUMNS[model_name],
            "residual_median": median,
            "mean": mean,
            "scale": scale,
            "threshold": global_threshold,
        }
    return fitted


def _predict_models(frame: pd.DataFrame, fitted: Mapping[str, Mapping]) -> pd.DataFrame:
    outputs: list[pd.DataFrame] = []
    for model_name in MODEL_ORDER:
        calibration_name = (
            "affine_sheaf" if model_name == "affine_sheaf_permuted" else model_name
        )
        spec = fitted[calibration_name]
        matrix, _, _, _ = _model_matrix(
            frame,
            RESIDUAL_COLUMNS[model_name],
            median=float(spec["residual_median"]),
            mean=spec["mean"],
            scale=spec["scale"],
        )
        probability = spec["model"].predict_proba(matrix)[:, 1]
        block = frame[
            ["seed", "scenario", "edge_id", "u", "v", "is_true_edge", "is_corrupted_edge"]
        ].copy()
        block["model"] = model_name
        block["probability"] = np.clip(probability, 1.0e-8, 1.0 - 1.0e-8)
        block["selected"] = (
            block["probability"].to_numpy(float) >= float(spec["threshold"])
        ).astype(int)
        block["residual"] = frame[RESIDUAL_COLUMNS[model_name]].to_numpy(float)
        outputs.append(block)
    return pd.concat(outputs, ignore_index=True)


def _ece(truth: np.ndarray, probability: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    total = 0.0
    for index in range(bins):
        if index == bins - 1:
            mask = (probability >= edges[index]) & (probability <= edges[index + 1])
        else:
            mask = (probability >= edges[index]) & (probability < edges[index + 1])
        if mask.any():
            total += float(mask.mean()) * abs(
                float(probability[mask].mean()) - float(truth[mask].mean())
            )
    return float(total)


def _case_metric_rows(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (seed, scenario, model), group in predictions.groupby(
        ["seed", "scenario", "model"], sort=True
    ):
        truth = group["is_true_edge"].to_numpy(int)
        probability = group["probability"].to_numpy(float)
        selected = group["selected"].to_numpy(int)
        confident = (probability <= 0.20) | (probability >= 0.80)
        prediction = probability >= 0.50
        false_mask = truth == 0
        corrupted = group["is_corrupted_edge"].to_numpy(int)
        conflict_ap = np.nan
        if corrupted.sum() > 0 and corrupted.sum() < len(corrupted):
            conflict_ap = average_precision_score(corrupted, 1.0 - probability)
        rows.append(
            {
                "seed": int(seed),
                "scenario": scenario,
                "model": model,
                "n_edges": len(group),
                "prevalence": float(truth.mean()),
                "pr_auc": float(average_precision_score(truth, probability)),
                "roc_auc": float(roc_auc_score(truth, probability)),
                "brier": float(brier_score_loss(truth, probability)),
                "log_loss": float(log_loss(truth, probability, labels=[0, 1])),
                "ece": _ece(truth, probability),
                "selected_f1": float(f1_score(truth, selected, zero_division=0)),
                "false_confidence_rate": float(
                    ((probability[false_mask] >= 0.80).sum()) / max(int(false_mask.sum()), 1)
                ),
                "abstention_coverage": float(confident.mean()),
                "abstention_accuracy": float(
                    (prediction[confident] == truth[confident]).mean()
                    if confident.any()
                    else np.nan
                ),
                "conflict_localisation_pr_auc": float(conflict_ap),
            }
        )
    return pd.DataFrame(rows)


def _paired_bootstrap(
    metrics: pd.DataFrame,
    left: str,
    right: str,
    metric: str,
    *,
    samples: int,
    seed: int,
    scenario: str = "all",
) -> dict:
    subset = metrics if scenario == "all" else metrics[metrics["scenario"] == scenario]
    pivot = subset.pivot(index="seed", columns="model", values=metric).dropna()
    differences = (pivot[left] - pivot[right]).to_numpy(float)
    if len(differences) == 0:
        return {
            "scenario": scenario,
            "left": left,
            "right": right,
            "metric": metric,
            "n_cases": 0,
            "mean_difference": np.nan,
            "ci95_low": np.nan,
            "ci95_high": np.nan,
        }
    rng = np.random.default_rng(int(seed))
    indices = rng.integers(0, len(differences), size=(int(samples), len(differences)))
    boot = differences[indices].mean(axis=1)
    return {
        "scenario": scenario,
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
        ("affine_sheaf", "weighted_graph"),
        ("affine_sheaf", "graph_laplacian"),
        ("affine_sheaf", "affine_sheaf_permuted"),
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
    for scenario in ("all", *SCENARIOS):
        for left, right in comparisons:
            for metric in metric_names:
                counter += 1
                rows.append(
                    _paired_bootstrap(
                        metrics,
                        left,
                        right,
                        metric,
                        samples=samples,
                        seed=2026072900 + counter,
                        scenario=scenario,
                    )
                )
    return pd.DataFrame(rows)


def _contrast(
    contrasts: pd.DataFrame,
    left: str,
    right: str,
    metric: str,
    scenario: str = "all",
) -> Mapping[str, object]:
    row = contrasts[
        (contrasts["left"] == left)
        & (contrasts["right"] == right)
        & (contrasts["metric"] == metric)
        & (contrasts["scenario"] == scenario)
    ]
    if len(row) != 1:
        raise KeyError((scenario, left, right, metric))
    return row.iloc[0].to_dict()


def _claim_decision(
    metrics: pd.DataFrame,
    contrasts: pd.DataFrame,
    diagnostics: pd.DataFrame,
    generator_provenance: Sequence[Mapping[str, object]],
) -> dict:
    identity_max = float(
        diagnostics["identity_affine_residual_max_abs_difference"].dropna().max()
    )
    identity_equivalence = identity_max <= 1.0e-10
    generator_independent = all(
        item.get("imports_hydrosheaf") is False for item in generator_provenance
    )
    finite = bool(
        np.isfinite(
            metrics[["pr_auc", "brier", "log_loss", "ece", "selected_f1"]]
        ).all().all()
    )

    sw_pr = _contrast(contrasts, "affine_sheaf", "weighted_graph", "pr_auc")
    sw_brier = _contrast(contrasts, "affine_sheaf", "weighted_graph", "brier")
    sg_pr = _contrast(contrasts, "affine_sheaf", "graph_laplacian", "pr_auc")
    sg_brier = _contrast(contrasts, "affine_sheaf", "graph_laplacian", "brier")
    sp_pr = _contrast(
        contrasts, "affine_sheaf", "affine_sheaf_permuted", "pr_auc"
    )
    identity_pr = _contrast(
        contrasts,
        "affine_sheaf",
        "graph_laplacian",
        "pr_auc",
        "identity_limit",
    )
    identity_brier = _contrast(
        contrasts,
        "affine_sheaf",
        "graph_laplacian",
        "brier",
        "identity_limit",
    )
    conflict_ap = _contrast(
        contrasts,
        "affine_sheaf",
        "graph_laplacian",
        "conflict_localisation_pr_auc",
        "incompatible_cycles",
    )

    execution_gate = bool(identity_equivalence and generator_independent and finite)
    primary_increment = bool(
        execution_gate
        and float(sw_pr["ci95_low"]) > 0.0
        and float(sw_brier["ci95_high"]) < 0.0
        and float(sg_pr["ci95_low"]) > 0.0
        and float(sg_brier["ci95_high"]) < 0.0
        and float(sp_pr["ci95_low"]) > 0.0
        and float(identity_pr["ci95_low"]) >= -0.02
        and float(identity_brier["ci95_high"]) <= 0.01
    )
    conditional_conflict = bool(
        execution_gate
        and np.isfinite(float(conflict_ap["ci95_low"]))
        and float(conflict_ap["ci95_low"]) > 0.0
    )
    if primary_increment:
        allowed = (
            "In this prospectively locked controlled-synthetic benchmark, the affine "
            "sheaf supplied incremental topology-discrimination and calibration value "
            "beyond competence-matched edge-local and identity-Laplacian weighted graphs, "
            "while collapsing to the graph formulation in the identity limit."
        )
    elif conditional_conflict:
        allowed = (
            "The affine sheaf did not satisfy the full incremental-superiority rule, but "
            "it conditionally improved localisation of planted global inconsistencies in "
            "the controlled-synthetic incompatible-cycle stratum."
        )
    else:
        allowed = (
            "No incremental sheaf claim is allowed; in these controlled-synthetic cases, "
            "the tested sheaf layer should be treated as an optional global-consistency "
            "diagnostic rather than a superior topology model."
        )
    return {
        "execution_and_equivalence_gate_passed": execution_gate,
        "controlled_synthetic_incremental_sheaf_claim_allowed": primary_increment,
        "conditional_conflict_localisation_claim_allowed": conditional_conflict,
        "identity_residual_max_abs_difference": identity_max,
        "primary_contrasts": {
            "sheaf_minus_weighted_graph_pr_auc": sw_pr,
            "sheaf_minus_weighted_graph_brier": sw_brier,
            "sheaf_minus_graph_laplacian_pr_auc": sg_pr,
            "sheaf_minus_graph_laplacian_brier": sg_brier,
            "native_sheaf_minus_permuted_map_pr_auc": sp_pr,
            "identity_limit_sheaf_minus_graph_pr_auc": identity_pr,
            "identity_limit_sheaf_minus_graph_brier": identity_brier,
            "conflict_localisation_sheaf_minus_graph": conflict_ap,
        },
        "allowed_claim": allowed,
        "guardrail": (
            "The result is a model-conditioned scalar, static, two-dimensional graph "
            "capability test. It is not field validation, does not establish general "
            "HydroSheaf superiority, and excludes temporal, spatial three-dimensional, "
            "and vadose-zone capability."
        ),
    }


def _write_publication_assets(
    metrics: pd.DataFrame,
    contrasts: pd.DataFrame,
    claim: Mapping[str, object],
) -> list[Path]:
    figure_dir = PROJECT / "figures" / "publication"
    table_dir = PROJECT / "tables" / "publication"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    summary = (
        metrics.groupby("model", sort=False)
        .agg(
            n_cases=("seed", "count"),
            pr_auc_mean=("pr_auc", "mean"),
            pr_auc_sd=("pr_auc", "std"),
            brier_mean=("brier", "mean"),
            ece_mean=("ece", "mean"),
            selected_f1_mean=("selected_f1", "mean"),
            false_confidence_rate_mean=("false_confidence_rate", "mean"),
            abstention_coverage_mean=("abstention_coverage", "mean"),
            abstention_accuracy_mean=("abstention_accuracy", "mean"),
        )
        .reset_index()
    )
    table_csv = table_dir / "table8_sheaf_vs_weighted_graph.csv"
    table_md = table_dir / "table8_sheaf_vs_weighted_graph.md"
    summary.round(6).to_csv(table_csv, index=False, lineterminator="\n")
    table_md.write_text(summary.round(4).to_markdown(index=False) + "\n", encoding="utf-8")

    labels = {
        "weighted_graph": "Edge-local\nweighted graph",
        "graph_laplacian": "Identity\ngraph Laplacian",
        "affine_sheaf": "Affine\nsheaf",
        "affine_sheaf_permuted": "Permuted-map\nsheaf",
    }
    colours = ["#8c8c8c", "#4c78a8", "#2a9d8f", "#e76f51"]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.8), constrained_layout=True)
    overall = metrics.groupby("model", sort=False).mean(numeric_only=True).reindex(MODEL_ORDER)
    axes[0, 0].bar(
        np.arange(len(MODEL_ORDER)), overall["pr_auc"], color=colours, edgecolor="black", linewidth=0.5
    )
    axes[0, 0].set_xticks(np.arange(len(MODEL_ORDER)), [labels[m] for m in MODEL_ORDER])
    axes[0, 0].set_ylabel("Mean case PR-AUC")
    axes[0, 0].set_title("a  Locked-test discrimination")

    axes[0, 1].bar(
        np.arange(len(MODEL_ORDER)), overall["brier"], color=colours, edgecolor="black", linewidth=0.5
    )
    axes[0, 1].set_xticks(np.arange(len(MODEL_ORDER)), [labels[m] for m in MODEL_ORDER])
    axes[0, 1].set_ylabel("Mean Brier score (lower is better)")
    axes[0, 1].set_title("b  Locked-test calibration")

    scenario = (
        metrics.groupby(["scenario", "model"], sort=False)["pr_auc"].mean().unstack("model")
    ).reindex(SCENARIOS)
    x = np.arange(len(SCENARIOS))
    width = 0.19
    for index, model in enumerate(MODEL_ORDER):
        axes[1, 0].bar(
            x + (index - 1.5) * width,
            scenario[model],
            width,
            color=colours[index],
            label=labels[model].replace("\n", " "),
        )
    axes[1, 0].set_xticks(x, [s.replace("_", "\n") for s in SCENARIOS])
    axes[1, 0].set_ylabel("Mean case PR-AUC")
    axes[1, 0].set_title("c  Scenario-specific performance")
    axes[1, 0].legend(frameon=False, fontsize=8, ncol=2)

    selected = contrasts[
        (contrasts["scenario"] == "all")
        & (contrasts["left"] == "affine_sheaf")
        & (contrasts["right"].isin(["weighted_graph", "graph_laplacian", "affine_sheaf_permuted"]))
        & (contrasts["metric"] == "pr_auc")
    ].copy()
    selected["label"] = selected["right"].map(
        {
            "weighted_graph": "vs edge-local graph",
            "graph_laplacian": "vs graph Laplacian",
            "affine_sheaf_permuted": "vs permuted maps",
        }
    )
    y = np.arange(len(selected))
    means = selected["mean_difference"].to_numpy(float)
    errors = np.vstack(
        [means - selected["ci95_low"].to_numpy(float), selected["ci95_high"].to_numpy(float) - means]
    )
    axes[1, 1].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 1].errorbar(means, y, xerr=errors, fmt="o", color="#2a9d8f", capsize=3)
    axes[1, 1].set_yticks(y, selected["label"])
    axes[1, 1].set_xlabel("Paired PR-AUC difference (sheaf - comparator)")
    axes[1, 1].set_title("d  Case-block bootstrap contrasts")

    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.18)
    fig.suptitle(
        "M7.4: incremental contribution of affine sheaf structure",
        fontsize=14,
        fontweight="bold",
    )
    figure_png = figure_dir / "figure6_sheaf_vs_weighted_graph.png"
    figure_pdf = figure_dir / "figure6_sheaf_vs_weighted_graph.pdf"
    figure_tif = figure_dir / "figure6_sheaf_vs_weighted_graph.tif"
    fig.savefig(figure_png, dpi=600, bbox_inches="tight")
    fig.savefig(figure_pdf, bbox_inches="tight")
    fig.savefig(figure_tif, dpi=300, pil_kwargs={"compression": "tiff_lzw"}, bbox_inches="tight")
    plt.close(fig)
    return [table_csv, table_md, figure_png, figure_pdf, figure_tif]


def _serialise_models(fitted: Mapping[str, Mapping]) -> dict:
    output: dict[str, dict] = {}
    for name, spec in fitted.items():
        model = spec["model"]
        output[name] = {
            "feature_order": ["prior_logit", "log1p_residual", "local_missing"],
            "residual_column": spec["residual_column"],
            "residual_median": float(spec["residual_median"]),
            "standardisation_mean": np.asarray(spec["mean"], dtype=float).tolist(),
            "standardisation_scale": np.asarray(spec["scale"], dtype=float).tolist(),
            "intercept": model.intercept_.astype(float).tolist(),
            "coefficients": model.coef_.astype(float).tolist(),
            "threshold": float(spec["threshold"]),
            "C": float(model.C),
            "class_weight": model.class_weight,
            "solver": model.solver,
            "max_iter": int(model.max_iter),
        }
    return output


def run_benchmark(
    *,
    output: Path,
    development_seeds: Sequence[int] = DEVELOPMENT_SEEDS,
    test_seeds: Sequence[int] = LOCKED_TEST_SEEDS,
    bootstrap_samples: int = BOOTSTRAP_SAMPLES,
    overwrite: bool = False,
    write_publication_assets: bool = True,
) -> dict:
    if tuple(map(int, test_seeds)) == LOCKED_TEST_SEEDS:
        _verify_lock(development_seeds, test_seeds, bootstrap_samples)
    output = _safe_output(output, overwrite)
    if set(map(int, development_seeds)) & set(map(int, test_seeds)):
        raise ValueError("Development and locked-test seeds must be disjoint.")

    all_features: list[pd.DataFrame] = []
    diagnostics: list[dict] = []
    provenance: list[Mapping[str, object]] = []
    development_first = int(min(development_seeds))
    test_first = int(min(test_seeds))
    for split, seeds, first_seed in (
        ("development", development_seeds, development_first),
        ("locked_test", test_seeds, test_first),
    ):
        for seed in seeds:
            scenario = scenario_for_seed(int(seed), first_seed)
            case = generate_independent_section_case(int(seed), scenario)
            frame, case_diagnostics = _feature_frame(case)
            frame["split"] = split
            all_features.append(frame)
            diagnostics.append(case_diagnostics | {"split": split})
            provenance.append(dict(case.provenance) | {"split": split})

    features = pd.concat(all_features, ignore_index=True)
    development = features[features["split"] == "development"].copy()
    locked_test = features[features["split"] == "locked_test"].copy()
    fitted = _fit_models(development)
    predictions = _predict_models(locked_test, fitted)
    metrics = _case_metric_rows(predictions)
    contrast_frame = _contrast_table(metrics, int(bootstrap_samples))
    diagnostics_frame = pd.DataFrame(diagnostics)
    claim = _claim_decision(metrics, contrast_frame, diagnostics_frame, provenance)

    model_summary = (
        metrics.groupby("model", sort=False)
        .mean(numeric_only=True)
        .reset_index()
    )
    scenario_summary = (
        metrics.groupby(["scenario", "model"], sort=False)
        .mean(numeric_only=True)
        .reset_index()
    )
    stable_frames = {
        "development_edge_features.csv": development,
        "locked_test_edge_features.csv": locked_test,
        "edge_predictions.csv": predictions,
        "case_metrics.csv": metrics,
        "model_summary.csv": model_summary,
        "scenario_summary.csv": scenario_summary,
        "paired_bootstrap_contrasts.csv": contrast_frame,
        "case_diagnostics.csv": diagnostics_frame,
    }
    for filename, frame in stable_frames.items():
        stable = frame.copy()
        float_columns = stable.select_dtypes(include=["float32", "float64"]).columns
        stable[float_columns] = stable[float_columns].round(12)
        stable.to_csv(output / filename, index=False, lineterminator="\n")
    (output / "frozen_models.json").write_text(
        json.dumps(_serialise_models(fitted), indent=2), encoding="utf-8"
    )
    (output / "generator_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    (output / "claim_decision.json").write_text(
        json.dumps(claim, indent=2), encoding="utf-8"
    )

    publication_assets: list[Path] = []
    if write_publication_assets:
        publication_assets = _write_publication_assets(metrics, contrast_frame, claim)
    manifest = {
        "schema_version": "1.0",
        "run_id": RUN_ID,
        "status": "PASS" if claim["execution_and_equivalence_gate_passed"] else "FAIL",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": _git_revision(),
        "development_seeds": list(map(int, development_seeds)),
        "locked_test_seeds": list(map(int, test_seeds)),
        "scenario_order": list(SCENARIOS),
        "bootstrap_samples": int(bootstrap_samples),
        "models": list(MODEL_ORDER),
        "development_truth_used_only_for_calibration": True,
        "locked_test_truth_joined_only_for_scoring": True,
        "excluded_capabilities": ["temporal_series", "graph_3d", "vadose"],
        "artifacts": {},
        "publication_assets": [str(path.relative_to(PROJECT)) for path in publication_assets],
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
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT / "results" / RUN_ID,
    )
    parser.add_argument("--bootstrap-samples", type=int, default=BOOTSTRAP_SAMPLES)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-publication-assets", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.bootstrap_samples < 200:
        raise ValueError("bootstrap-samples must be at least 200")
    result = run_benchmark(
        output=args.output,
        development_seeds=(QUICK_DEVELOPMENT_SEEDS if args.quick else DEVELOPMENT_SEEDS),
        test_seeds=(QUICK_TEST_SEEDS if args.quick else LOCKED_TEST_SEEDS),
        bootstrap_samples=int(args.bootstrap_samples),
        overwrite=bool(args.overwrite),
        write_publication_assets=not args.no_publication_assets,
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
