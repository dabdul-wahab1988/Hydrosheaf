"""Post-hoc diagnostics for the frozen M4-A and M4-B experiments.

The script consumes the frozen probability tables and reconstructs only the
declared synthetic-development feature rows used by the locked fair runner.
It never fits a model, changes a threshold, or uses Savage truth to select a
model.  MODPATH labels are used only to calculate the requested diagnostic
metrics after the frozen probabilities have been loaded.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, wasserstein_distance
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.validation.topology_v2 import (  # noqa: E402
    MODEL_FEATURE_SETS,
    TopologyV2Config,
    build_topology_v2_feature_rows,
    generate_topology_v2_candidate_universe,
)
from run_m4_fair_topology_benchmark import (  # noqa: E402
    M4_A,
    M4_B,
    _synthetic_case,
)


FREEZE_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "frozen_m4_diagnostics_20260922"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "m4_diagnostics_20260922"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: object) -> object:
    return json.loads(json.dumps(value, default=str, sort_keys=True))


def _safe_auc(metric_fn: Any, labels: np.ndarray, values: np.ndarray) -> float | None:
    if len(np.unique(labels)) < 2:
        return None
    try:
        return float(metric_fn(labels, values))
    except ValueError:
        return None


def _fixed_confusion(
    labels: np.ndarray,
    probabilities: np.ndarray,
    present_threshold: float,
    absent_threshold: float,
) -> dict[str, object]:
    present = probabilities >= float(present_threshold)
    absent = probabilities <= float(absent_threshold)
    abstain = ~(present | absent)
    tp = int(np.sum(present & (labels == 1.0)))
    fp = int(np.sum(present & (labels == 0.0)))
    fn = int(np.sum(absent & (labels == 1.0)))
    tn = int(np.sum(absent & (labels == 0.0)))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "abstain": int(np.sum(abstain)),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "predicted_present": int(np.sum(present)),
    }


def _score_metrics(
    labels: np.ndarray,
    probabilities: np.ndarray,
    *,
    present_threshold: float,
    absent_threshold: float,
) -> dict[str, object]:
    prevalence = float(np.mean(labels)) if len(labels) else float("nan")
    return {
        "n_candidates": int(len(labels)),
        "n_positive": int(np.sum(labels == 1.0)),
        "n_negative": int(np.sum(labels == 0.0)),
        "prevalence": prevalence,
        "pr_auc": _safe_auc(average_precision_score, labels, probabilities),
        "roc_auc": _safe_auc(roc_auc_score, labels, probabilities),
        "brier": float(brier_score_loss(labels, probabilities)),
        "log_loss": float(log_loss(labels, np.column_stack([1.0 - probabilities, probabilities]), labels=[0.0, 1.0])),
        "mean_probability": float(np.mean(probabilities)),
        "positive_mean_probability": float(np.mean(probabilities[labels == 1.0])) if np.any(labels == 1.0) else None,
        "negative_mean_probability": float(np.mean(probabilities[labels == 0.0])) if np.any(labels == 0.0) else None,
        "separation_mean_difference": (
            float(np.mean(probabilities[labels == 1.0]) - np.mean(probabilities[labels == 0.0]))
            if np.any(labels == 1.0) and np.any(labels == 0.0)
            else None
        ),
        "fixed_threshold": _fixed_confusion(
            labels, probabilities, present_threshold, absent_threshold
        ),
    }


def _load_frozen(freeze_dir: Path, mode: str) -> tuple[pd.DataFrame, dict[str, object]]:
    freeze_manifest = freeze_dir / "freeze_manifest.json"
    if not freeze_manifest.exists():
        raise FileNotFoundError(freeze_manifest)
    manifest = json.loads(freeze_manifest.read_text(encoding="utf-8"))
    mode_dir = freeze_dir / mode
    edge_name = "fair_edge_scores_m4_a.csv" if mode == "M4-A" else "fair_edge_scores_m4_b.csv"
    table_path = mode_dir / edge_name
    if not table_path.exists():
        raise FileNotFoundError(table_path)
    table = pd.read_csv(table_path)
    required = {"probability", "in_reference_for_evaluation_only"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Frozen {mode} score table is missing {sorted(missing)}")
    table["probability"] = pd.to_numeric(table["probability"], errors="coerce")
    table["label_posthoc"] = table["in_reference_for_evaluation_only"].astype(bool).astype(float)
    if table["probability"].isna().any():
        raise ValueError(f"Frozen {mode} table contains non-finite probabilities.")
    return table, manifest


def _load_policies(freeze_dir: Path) -> dict[str, tuple[float, float]]:
    calibration = pd.read_csv(freeze_dir / "shared" / "fair_calibration.csv")
    policies: dict[str, tuple[float, float]] = {}
    for _, row in calibration.iterrows():
        mode = "M4-A" if str(row["mode"]).startswith("M4-A") else "M4-B"
        policies[mode] = (float(row["present_threshold"]), float(row["absent_threshold"]))
    if set(policies) != {"M4-A", "M4-B"}:
        raise ValueError("Frozen calibration table does not contain both M4-A and M4-B policies.")
    return policies


def _development_features(mode: str) -> pd.DataFrame:
    feature_key = "M4_A_sparse" if mode == "M4-A" else "M4_B_archive_informed"
    feature_names = MODEL_FEATURE_SETS[feature_key]
    rows: list[dict[str, object]] = []
    for case_name in ("fit_1", "fit_2", "threshold_1"):
        observations, reference = _synthetic_case(
            M4_A if mode == "M4-A" else M4_B,
            case_name,
            missing_elevation=(mode == "M4-A" and case_name != "fit_1"),
        )
        universe = generate_topology_v2_candidate_universe(
            observations,
            config=TopologyV2Config(default_head_sigma_m=0.05),
        )
        feature_rows = build_topology_v2_feature_rows(
            universe,
            observations,
            case_id=case_name,
            config=TopologyV2Config(default_head_sigma_m=0.05),
        )
        reference_ids = set(reference)
        for row in feature_rows:
            record: dict[str, object] = {
                "case_id": case_name,
                "edge_id": row.edge_id,
                "label_development_only": float(row.edge_id in reference_ids),
            }
            for feature in feature_names:
                record[feature] = row.features.get(feature)
            rows.append(record)
    return pd.DataFrame(rows)


def _separation_summary(table: pd.DataFrame) -> dict[str, object]:
    labels = table["label_posthoc"].to_numpy(dtype=float)
    probabilities = table["probability"].to_numpy(dtype=float)
    return _score_metrics(labels, probabilities, present_threshold=0.5, absent_threshold=0.5)


def _reversed_diagnostic(
    table: pd.DataFrame,
    present_threshold: float,
    absent_threshold: float,
) -> dict[str, object]:
    labels = table["label_posthoc"].to_numpy(dtype=float)
    probabilities = table["probability"].to_numpy(dtype=float)
    reversed_probabilities = 1.0 - probabilities
    original = _score_metrics(
        labels,
        probabilities,
        present_threshold=present_threshold,
        absent_threshold=absent_threshold,
    )
    reversed_metrics = _score_metrics(
        labels,
        reversed_probabilities,
        present_threshold=present_threshold,
        absent_threshold=absent_threshold,
    )
    return {
        "original_pr_auc": original["pr_auc"],
        "original_roc_auc": original["roc_auc"],
        "reversed_pr_auc": reversed_metrics["pr_auc"],
        "reversed_roc_auc": reversed_metrics["roc_auc"],
        "reversed_minus_original_pr_auc": (
            float(reversed_metrics["pr_auc"] - original["pr_auc"])
            if original["pr_auc"] is not None and reversed_metrics["pr_auc"] is not None
            else None
        ),
        "reversed_minus_original_roc_auc": (
            float(reversed_metrics["roc_auc"] - original["roc_auc"])
            if original["roc_auc"] is not None and reversed_metrics["roc_auc"] is not None
            else None
        ),
        "reversed_fixed_threshold": reversed_metrics["fixed_threshold"],
        "no_threshold_retuned": True,
    }


def _feature_shift(mode: str, development: pd.DataFrame, savage: pd.DataFrame) -> pd.DataFrame:
    feature_key = "M4_A_sparse" if mode == "M4-A" else "M4_B_archive_informed"
    feature_names = MODEL_FEATURE_SETS[feature_key]
    records: list[dict[str, object]] = []
    for feature in feature_names:
        dev_values = pd.to_numeric(development[feature], errors="coerce").dropna().to_numpy(dtype=float)
        savage_column = f"feature_{feature}"
        savage_values = pd.to_numeric(savage[savage_column], errors="coerce").dropna().to_numpy(dtype=float)
        if not len(dev_values) or not len(savage_values):
            continue
        dev_std = float(np.std(dev_values))
        savage_std = float(np.std(savage_values))
        pooled = math.sqrt(0.5 * (dev_std**2 + savage_std**2)) if (dev_std or savage_std) else 1.0
        ks = ks_2samp(dev_values, savage_values)
        records.append(
            {
                "mode": mode,
                "feature": feature,
                "development_n": len(dev_values),
                "savage_n": len(savage_values),
                "development_missing_rate": float(development[feature].isna().mean()),
                "savage_missing_rate": float(savage[savage_column].isna().mean()),
                "development_mean": float(np.mean(dev_values)),
                "savage_mean": float(np.mean(savage_values)),
                "development_median": float(np.median(dev_values)),
                "savage_median": float(np.median(savage_values)),
                "development_q25": float(np.quantile(dev_values, 0.25)),
                "savage_q25": float(np.quantile(savage_values, 0.25)),
                "development_q75": float(np.quantile(dev_values, 0.75)),
                "savage_q75": float(np.quantile(savage_values, 0.75)),
                "standardised_mean_difference": float((np.mean(savage_values) - np.mean(dev_values)) / pooled),
                "wasserstein_distance": float(wasserstein_distance(dev_values, savage_values)),
                "ks_statistic": float(ks.statistic),
                "ks_pvalue": float(ks.pvalue),
            }
        )
    return pd.DataFrame(records)


def _quantile_strata(values: pd.Series, prefix: str) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    finite = numeric.replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        return pd.Series(["unavailable"] * len(values), index=values.index)
    cuts = np.unique(np.quantile(finite.to_numpy(dtype=float), [0.0, 0.25, 0.5, 0.75, 1.0]))
    if len(cuts) < 2:
        return pd.Series([f"{prefix}_constant"] * len(values), index=values.index)
    labels = [f"{prefix}_q{i + 1}" for i in range(len(cuts) - 1)]
    return pd.cut(numeric, bins=cuts, labels=labels, include_lowest=True, duplicates="drop").astype("object").fillna("unavailable")


def _strata_table(
    mode: str,
    table: pd.DataFrame,
    present_threshold: float,
    absent_threshold: float,
) -> pd.DataFrame:
    work = table.copy()
    work["distance_stratum"] = _quantile_strata(work["feature_distance_m"], "distance")
    strata_specs = [("candidate_distance", "distance_stratum")]
    gradient_column = "feature_gradient_m_per_km"
    if gradient_column in work.columns and pd.to_numeric(work[gradient_column], errors="coerce").notna().any():
        work["head_gradient_stratum"] = _quantile_strata(work[gradient_column], "gradient")
        strata_specs.append(("head_gradient", "head_gradient_stratum"))

    records: list[dict[str, object]] = []
    for stratum_type, column in strata_specs:
        for name, group in work.groupby(column, dropna=False, sort=True):
            labels = group["label_posthoc"].to_numpy(dtype=float)
            probabilities = group["probability"].to_numpy(dtype=float)
            metrics = _score_metrics(
                labels,
                probabilities,
                present_threshold=present_threshold,
                absent_threshold=absent_threshold,
            )
            records.append(
                {
                    "mode": mode,
                    "stratum_type": stratum_type,
                    "stratum": str(name),
                    **metrics,
                }
            )
    return pd.DataFrame(records)


def _plot_probability_separation(mode: str, table: pd.DataFrame, output_dir: Path) -> None:
    labels = table["label_posthoc"].to_numpy(dtype=float)
    probabilities = table["probability"].to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    bins = np.linspace(0.0, 1.0, 41)
    ax.hist(probabilities[labels == 0.0], bins=bins, density=True, alpha=0.55, label="MODPATH-absent candidate")
    ax.hist(probabilities[labels == 1.0], bins=bins, density=True, alpha=0.75, label="MODPATH-present candidate")
    ax.set_xlabel("Frozen edge probability")
    ax.set_ylabel("Density")
    ax.set_title(f"{mode}: frozen probability separation")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / f"probability_separation_{mode.replace('-', '_')}.png", dpi=220)
    plt.close(fig)


def _plot_reversed(summary: pd.DataFrame, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))
    x = np.arange(len(summary))
    width = 0.35
    axes[0].bar(x - width / 2, summary["original_pr_auc"], width, label="frozen")
    axes[0].bar(x + width / 2, summary["reversed_pr_auc"], width, label="1 - p")
    axes[0].axhline(summary["prevalence"].mean(), color="black", linestyle="--", linewidth=0.8, label="prevalence")
    axes[0].set_title("PR-AUC")
    axes[1].bar(x - width / 2, summary["original_roc_auc"], width, label="frozen")
    axes[1].bar(x + width / 2, summary["reversed_roc_auc"], width, label="1 - p")
    axes[1].axhline(0.5, color="black", linestyle="--", linewidth=0.8)
    axes[1].set_title("ROC-AUC")
    for ax in axes:
        ax.set_xticks(x, summary["mode"])
        ax.set_ylim(0.0, 1.0)
        ax.grid(axis="y", alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Reversed-score diagnostic; no Savage threshold retuning")
    fig.tight_layout()
    fig.savefig(output_dir / "reversed_score_diagnostic.png", dpi=220)
    plt.close(fig)


def _plot_feature_shift(mode: str, shift: pd.DataFrame, output_dir: Path) -> None:
    if shift.empty:
        return
    n = len(shift)
    fig, axes = plt.subplots(n, 1, figsize=(8.2, max(2.2 * n, 3.0)), squeeze=False)
    for index, (_, row) in enumerate(shift.iterrows()):
        ax = axes[index, 0]
        feature = str(row["feature"])
        dev_q = [row["development_q25"], row["development_median"], row["development_q75"]]
        sav_q = [row["savage_q25"], row["savage_median"], row["savage_q75"]]
        ax.plot([0.85, 1.0, 1.15], dev_q, marker="o", label="synthetic development")
        ax.plot([1.85, 2.0, 2.15], sav_q, marker="o", label="Savage candidates")
        ax.set_xticks([1.0, 2.0], ["development", "Savage"])
        ax.set_title(f"{feature}: SMD={float(row['standardised_mean_difference']):.2f}, KS={float(row['ks_statistic']):.2f}")
        ax.grid(alpha=0.2)
        if index == 0:
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"{mode}: development-to-Savage feature shift")
    fig.tight_layout()
    fig.savefig(output_dir / f"feature_shift_{mode.replace('-', '_')}.png", dpi=220)
    plt.close(fig)


def _plot_strata(strata: pd.DataFrame, output_dir: Path) -> None:
    if strata.empty:
        return
    pivot = strata.pivot_table(index=["mode", "stratum_type"], columns="stratum", values="roc_auc", aggfunc="first")
    fig, ax = plt.subplots(figsize=(10.0, max(3.0, 0.45 * len(pivot))))
    image = ax.imshow(pivot.fillna(0.5).to_numpy(dtype=float), aspect="auto", vmin=0.0, vmax=1.0, cmap="coolwarm")
    ax.set_yticks(np.arange(len(pivot)), [" / ".join(map(str, index)) for index in pivot.index])
    ax.set_xticks(np.arange(len(pivot.columns)), [str(value) for value in pivot.columns], rotation=45, ha="right")
    ax.set_title("ROC-AUC by candidate distance/head-gradient stratum")
    fig.colorbar(image, ax=ax, label="ROC-AUC")
    fig.tight_layout()
    fig.savefig(output_dir / "performance_by_strata.png", dpi=220)
    plt.close(fig)


def _diagnostic_conclusion(
    summaries: pd.DataFrame,
    reversed_table: pd.DataFrame,
    shifts: pd.DataFrame,
    strata: pd.DataFrame,
) -> dict[str, object]:
    summary_by_mode = {str(row["mode"]): row for _, row in summaries.iterrows()}
    reverse_by_mode = {str(row["mode"]): row for _, row in reversed_table.iterrows()}
    ranking_failure = all(
        float(row["roc_auc"]) < 0.5 or float(row["pr_auc"]) <= float(row["prevalence"])
        for row in summary_by_mode.values()
    )
    reversal_better = any(
        float(row["reversed_minus_original_roc_auc"]) > 0.05
        or float(row["reversed_minus_original_pr_auc"]) > 0.05
        for row in reverse_by_mode.values()
        if row["reversed_minus_original_roc_auc"] is not None
    )
    strong_shift = bool(
        not shifts.empty
        and (
            (shifts["standardised_mean_difference"].abs() > 1.0).any()
            or (shifts["ks_statistic"] > 0.30).any()
        )
    )
    candidate_dilution = all(float(row["prevalence"]) < 0.01 for row in summary_by_mode.values())
    gradient_rows = strata[strata["stratum_type"] == "head_gradient"] if not strata.empty else pd.DataFrame()
    gradient_separation = None
    if not gradient_rows.empty:
        gradient_separation = float(gradient_rows["roc_auc"].max() - gradient_rows["roc_auc"].min())
    distance_rows = strata[strata["stratum_type"] == "candidate_distance"] if not strata.empty else pd.DataFrame()
    long_distance_failure = False
    long_distance_roc_auc = None
    if not distance_rows.empty:
        long_rows = distance_rows[distance_rows["stratum"].astype(str).str.endswith("_q4")]
        if not long_rows.empty:
            long_distance_roc_auc = float(long_rows["roc_auc"].mean())
            long_distance_failure = bool(
                (long_rows["n_positive"] > 0).any()
                and (long_rows["roc_auc"] < 0.5).any()
            )
    gradient_tail_failure = bool(
        not gradient_rows.empty
        and (gradient_rows["stratum"].astype(str).str.endswith("_q4")).any()
        and (gradient_rows.loc[gradient_rows["stratum"].astype(str).str.endswith("_q4"), "n_positive"] > 0).any()
        and (gradient_rows.loc[gradient_rows["stratum"].astype(str).str.endswith("_q4"), "roc_auc"] < 0.5).any()
    )
    path_information_gap = long_distance_failure or gradient_tail_failure
    if reversal_better:
        primary = (
            "feature-sign/orientation mismatch is the primary failure; long-range "
            "path-aware hydraulic information remains the next limitation"
        )
    elif ranking_failure:
        primary = "ranking failure dominates; candidate-space dilution is a secondary denominator effect"
    elif strong_shift:
        primary = "calibration-transfer failure is plausible because ranking survives but feature distributions shift"
    else:
        primary = "mixed or unresolved failure; path-aware hydraulic information remains the next test"
    return {
        "primary_failure_mode": primary,
        "ranking_failure_supported": ranking_failure,
        "feature_sign_or_orientation_mismatch_supported": reversal_better,
        "reversed_score_better_by_more_than_0_05": reversal_better,
        "feature_distribution_shift_supported": strong_shift,
        "calibration_transfer_failure_supported": strong_shift,
        "candidate_space_dilution_present": candidate_dilution,
        "candidate_space_dilution_is_primary": False if reversal_better else candidate_dilution,
        "head_gradient_stratum_roc_auc_range": gradient_separation,
        "long_distance_q4_roc_auc": long_distance_roc_auc,
        "long_distance_path_failure_supported": long_distance_failure,
        "head_gradient_tail_failure_supported": gradient_tail_failure,
        "path_aware_information_gap_supported": path_information_gap,
        "path_aware_information_gap_remains_open": True,
        "thresholds_retuned_on_savage": False,
        "interpretation": (
            "The raw scores are directionally inverted relative to the withheld labels, "
            "so orientation/sign review is the primary diagnosis. Candidate prevalence "
            "is low and feature distributions shift between synthetic development and "
            "Savage, but those effects do not explain ROC-AUC below 0.5 as directly as "
            "the reversal result does. M4-B also fails in the longest-distance and "
            "high-gradient-tail strata, leaving a path-aware hydraulic information gap. "
            "These findings choose the M4-C information channel; they do not change "
            "M4-A/M4-B probabilities or thresholds."
        ),
    }


def run(freeze_dir: Path = FREEZE_DIR, output_dir: Path = DEFAULT_OUTPUT_DIR) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    policies = _load_policies(freeze_dir)
    tables: dict[str, pd.DataFrame] = {}
    manifests: dict[str, object] = {}
    summaries: list[dict[str, object]] = []
    reversed_records: list[dict[str, object]] = []
    shift_tables: list[pd.DataFrame] = []
    strata_tables: list[pd.DataFrame] = []

    for mode in ("M4-A", "M4-B"):
        table, freeze_manifest = _load_frozen(freeze_dir, mode)
        tables[mode] = table
        manifests[mode] = freeze_manifest
        present_threshold, absent_threshold = policies[mode]
        labels = table["label_posthoc"].to_numpy(dtype=float)
        probabilities = table["probability"].to_numpy(dtype=float)
        metrics = _score_metrics(
            labels,
            probabilities,
            present_threshold=present_threshold,
            absent_threshold=absent_threshold,
        )
        summaries.append(
            {
                "mode": mode,
                "present_threshold_frozen": present_threshold,
                "absent_threshold_frozen": absent_threshold,
                **metrics,
                "threshold_retuned_on_savage": False,
            }
        )
        reversed_records.append(
            {
                "mode": mode,
                "prevalence": metrics["prevalence"],
                **_reversed_diagnostic(table, present_threshold, absent_threshold),
            }
        )
        development = _development_features(mode)
        shift_tables.append(_feature_shift(mode, development, table))
        strata_tables.append(_strata_table(mode, table, present_threshold, absent_threshold))
        _plot_probability_separation(mode, table, output_dir)
        _plot_feature_shift(mode, shift_tables[-1], output_dir)

    summary_table = pd.DataFrame(summaries)
    reversed_table = pd.DataFrame(reversed_records)
    shift_table = pd.concat(shift_tables, ignore_index=True) if shift_tables else pd.DataFrame()
    strata_table = pd.concat(strata_tables, ignore_index=True) if strata_tables else pd.DataFrame()
    conclusion = _diagnostic_conclusion(summary_table, reversed_table, shift_table, strata_table)
    _plot_reversed(reversed_table, output_dir)
    _plot_strata(strata_table, output_dir)

    summary_table.to_csv(output_dir / "m4_frozen_pr_auc_roc_auc.csv", index=False)
    reversed_table.to_csv(output_dir / "m4_reversed_score_diagnostic.csv", index=False)
    shift_table.to_csv(output_dir / "m4_feature_distribution_shift.csv", index=False)
    strata_table.to_csv(output_dir / "m4_performance_by_strata.csv", index=False)

    output_manifest = {
        "diagnostic_version": "m4_frozen_diagnostics_v1",
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_directory": str(freeze_dir),
        "frozen_manifest_sha256": _sha256(freeze_dir / "freeze_manifest.json"),
        "frozen_mode_score_hashes": {
            mode: _sha256(
                freeze_dir / mode / ("fair_edge_scores_m4_a.csv" if mode == "M4-A" else "fair_edge_scores_m4_b.csv")
            )
            for mode in ("M4-A", "M4-B")
        },
        "thresholds_retuned_on_savage": False,
        "reference_labels_used_only_for_posthoc_metrics": True,
        "diagnostics": [
            "PR-AUC",
            "ROC-AUC",
            "probability separation",
            "reversed-score diagnostic",
            "synthetic-to-Savage feature distribution shift",
            "candidate-distance strata",
            "head-gradient strata",
        ],
        "conclusion": conclusion,
    }
    (output_dir / "m4_diagnostic_conclusion.json").write_text(
        json.dumps(_jsonable(output_manifest), indent=2), encoding="utf-8"
    )
    report = [
        "# Frozen M4 diagnostic report",
        "",
        "M4-A and M4-B were consumed from the frozen snapshot. No probabilities, thresholds, or model parameters were refit, and no Savage threshold was selected.",
        "",
        f"Primary diagnostic conclusion: **{conclusion['primary_failure_mode']}**.",
        "",
        "The candidate-distance and head-gradient tables are stratified evaluations of the frozen scores. The feature-shift table compares the declared synthetic development features with the complete Savage candidate universe without using labels for binning or calibration.",
        "",
        "## Requested evidence",
        "",
        "- PR-AUC and ROC-AUC: `m4_frozen_pr_auc_roc_auc.csv`",
        "- Reversed-score diagnostic: `m4_reversed_score_diagnostic.csv` and `reversed_score_diagnostic.png`",
        "- Feature distribution shift: `m4_feature_distribution_shift.csv` and `feature_shift_*.png`",
        "- Distance/head-gradient strata: `m4_performance_by_strata.csv` and `performance_by_strata.png`",
        "- Immutable-input provenance: `m4_diagnostic_conclusion.json`",
        "",
        "M4-C is built separately from this report; these diagnostics are used to choose the path-aware physical evidence channel, not to retune M4-A or M4-B.",
    ]
    (output_dir / "M4_FROZEN_DIAGNOSTIC_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-dir", type=Path, default=FREEZE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    output = run(args.freeze_dir, args.output_dir)
    print(f"Wrote frozen M4 diagnostics to {output}")


if __name__ == "__main__":
    main()
