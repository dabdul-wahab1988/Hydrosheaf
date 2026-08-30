"""Regenerate estimator-diagnostic assets without changing the locked runner."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RUN = (
    PROJECT
    / "results"
    / "RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01"
    / "locked_test"
)
ADJUSTED = (
    PROJECT
    / "results"
    / "post_review_audit_20260730"
    / "m7_5_full_family_simultaneous.csv"
)
TABLES = PROJECT / "tables" / "publication"
FIGURES = PROJECT / "figures" / "publication"
PRINTED_WIDTH_IN = 7.08
MODEL_ORDER = (
    "edge_local",
    "identity_graph",
    "original_affine_global",
    "original_hybrid",
    "robust_affine_global",
    "robust_hybrid",
    "robust_hybrid_permuted",
)
SCENARIOS = (
    "identity_limit",
    "heterogeneous_affine",
    "incompatible_cycles",
    "noisy_missing",
)
LABELS = {
    "edge_local": "Edge-local",
    "identity_graph": "Identity",
    "original_affine_global": "Original global",
    "original_hybrid": "Original local-first",
    "robust_affine_global": "Robust global",
    "robust_hybrid": "Robust local-first",
    "robust_hybrid_permuted": "Permuted",
}
COLOURS = ["#777777", "#4C78A8", "#72A0C1", "#59A14F", "#2A9D8F", "#006D77", "#E76F51"]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def contrast(
    adjusted: pd.DataFrame,
    left: str,
    right: str,
    metric: str,
    *,
    scenario: str = "all",
) -> pd.Series:
    row = adjusted[
        (adjusted["calibration_regime"] == "separate_crossfit")
        & (adjusted["scenario"] == scenario)
        & (adjusted["left"] == left)
        & (adjusted["right"] == right)
        & (adjusted["metric"] == metric)
    ]
    if len(row) != 1:
        raise KeyError((scenario, left, right, metric))
    return row.iloc[0]


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(
        FIGURES / f"{stem}.tif",
        dpi=300,
        pil_kwargs={"compression": "tiff_lzw"},
        bbox_inches="tight",
    )


def main() -> int:
    configure_style()
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(RUN / "case_metrics.csv")
    contrasts = pd.read_csv(RUN / "paired_bootstrap_contrasts.csv")
    adjusted = pd.read_csv(ADJUSTED)
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
    summary.round(6).to_csv(
        TABLES / "table9_robust_hybrid_sheaf_diagnostics.csv",
        index=False,
        lineterminator="\n",
    )
    display = summary[
        [
            "model",
            "n_cases",
            "pr_auc_mean",
            "brier_mean",
            "log_loss_mean",
            "selected_f1_mean",
        ]
    ].copy()
    display["model"] = display["model"].map(
        {
            key: value.replace("\n", " ")
            for key, value in LABELS.items()
        }
    )
    display = display.rename(
        columns={
            "model": "Model",
            "n_cases": "Cases",
            "pr_auc_mean": "PR-AUC",
            "brier_mean": "Brier",
            "log_loss_mean": "Log loss",
            "selected_f1_mean": "Selected F1",
        }
    )
    display.round(6).to_csv(
        TABLES / "table9_robust_hybrid_sheaf.csv",
        index=False,
        lineterminator="\n",
    )
    (TABLES / "table9_robust_hybrid_sheaf.md").write_text(
        display.round(4).to_markdown(index=False) + "\n", encoding="utf-8"
    )
    contrasts.round(6).to_csv(
        TABLES / "tableS9_robust_hybrid_contrasts.csv",
        index=False,
        lineterminator="\n",
    )
    (TABLES / "tableS9_robust_hybrid_contrasts.md").write_text(
        contrasts.round(4).to_markdown(index=False) + "\n", encoding="utf-8"
    )

    overall = primary.groupby("model", sort=False).mean(numeric_only=True).reindex(MODEL_ORDER)
    y_models = np.arange(len(MODEL_ORDER))
    fig, axes = plt.subplots(2, 2, figsize=(PRINTED_WIDTH_IN, 5.25))
    axes[0, 0].barh(
        y_models,
        overall["pr_auc"],
        color=COLOURS,
        edgecolor="black",
        linewidth=0.4,
    )
    axes[0, 0].set_yticks(y_models, [LABELS[m] for m in MODEL_ORDER])
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_xlabel("Mean case PR-AUC")
    axes[0, 0].set_title("a  Locked-test discrimination", loc="left", fontweight="bold")

    axes[0, 1].barh(
        y_models,
        overall["log_loss"],
        color=COLOURS,
        edgecolor="black",
        linewidth=0.4,
    )
    axes[0, 1].set_yticks(y_models, [LABELS[m] for m in MODEL_ORDER])
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlabel("Mean log loss (lower is better)")
    axes[0, 1].set_title("b  Locked-test calibration", loc="left", fontweight="bold")

    scenario_rows = pd.DataFrame(
        [
            contrast(
                adjusted,
                "robust_hybrid",
                "edge_local",
                "pr_auc",
                scenario=scenario,
            )
            for scenario in SCENARIOS
        ]
    )
    means = scenario_rows["mean_difference"].to_numpy(float)
    errors = np.vstack(
        [
            means - scenario_rows["simultaneous_ci95_low"].to_numpy(float),
            scenario_rows["simultaneous_ci95_high"].to_numpy(float) - means,
        ]
    )
    axes[1, 0].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 0].errorbar(
        means,
        np.arange(len(SCENARIOS)),
        xerr=errors,
        fmt="o",
        color="#006D77",
        capsize=3,
    )
    axes[1, 0].set_yticks(
        np.arange(len(SCENARIOS)),
        [scenario.replace("_", " ") for scenario in SCENARIOS],
    )
    axes[1, 0].set_xlabel("PR-AUC difference vs edge-local")
    axes[1, 0].set_title(
        "c  Scenario diagnostics (family-wise 95%)",
        loc="left",
        fontweight="bold",
    )

    mechanisms = [
        ("Original local-first vs global", "original_hybrid", "original_affine_global"),
        ("Robust global vs original", "robust_affine_global", "original_affine_global"),
        ("Robust local-first vs original", "robust_hybrid", "original_hybrid"),
        ("Native vs permuted", "robust_hybrid", "robust_hybrid_permuted"),
    ]
    mechanism_rows = pd.DataFrame(
        [
            dict(contrast(adjusted, left, right, "pr_auc")) | {"label": label}
            for label, left, right in mechanisms
        ]
    )
    means = mechanism_rows["mean_difference"].to_numpy(float)
    errors = np.vstack(
        [
            means - mechanism_rows["simultaneous_ci95_low"].to_numpy(float),
            mechanism_rows["simultaneous_ci95_high"].to_numpy(float) - means,
        ]
    )
    axes[1, 1].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 1].errorbar(
        means,
        np.arange(len(mechanisms)),
        xerr=errors,
        fmt="o",
        color="#2A9D8F",
        capsize=3,
    )
    axes[1, 1].set_yticks(np.arange(len(mechanisms)), mechanism_rows["label"])
    axes[1, 1].set_xlabel("Paired PR-AUC difference")
    axes[1, 1].set_title(
        "d  Mechanism attribution (family-wise 95%)",
        loc="left",
        fontweight="bold",
    )
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.16)
    fig.suptitle(
        "Local-first/global-fallback estimator diagnostic",
        fontsize=11,
        fontweight="bold",
    )
    fig.subplots_adjust(
        hspace=0.48, wspace=0.55, left=0.16, right=0.98, top=0.92, bottom=0.11
    )
    save_figure(fig, "figure7_robust_hybrid_sheaf")
    plt.close(fig)

    manifest_path = FIGURES / "figure_source_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    if "Printed width" not in manifest:
        manifest["Printed width"] = ""
    if "Minimum label size" not in manifest:
        manifest["Minimum label size"] = ""
    manifest = manifest[manifest["Manuscript item"] != "Figure 7"]
    row = pd.DataFrame(
        [
            {
                "Manuscript item": "Figure 7",
                "Artifact stem": "figure7_robust_hybrid_sheaf",
                "Locked source": (
                    "RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01/locked_test/"
                    "case_metrics.csv; paired_bootstrap_contrasts.csv; "
                    "post_review_audit_20260730/m7_5_full_family_simultaneous.csv"
                ),
                "Purpose": (
                    "Fresh-seed local-first/global-fallback estimator and "
                    "family-wise mechanism diagnostic"
                ),
                "Printed width": f"{PRINTED_WIDTH_IN:.2f} in (180 mm)",
                "Minimum label size": "8.0 pt",
            }
        ]
    )
    pd.concat([manifest, row], ignore_index=True).to_csv(
        manifest_path, index=False, lineterminator="\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
