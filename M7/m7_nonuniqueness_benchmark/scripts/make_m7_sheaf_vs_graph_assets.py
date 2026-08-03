"""Render publication assets from the immutable representation-benchmark tables."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
RUN = PROJECT / "results" / "RUN-M7-SHEAF-VS-GRAPH-20260729-01"
SYSTEM_RUN = PROJECT / "provenance" / "runs" / "RUN-M7-SYSTEM-20260728-01"
ADJUSTED = (
    PROJECT
    / "results"
    / "post_review_audit_20260730"
    / "m7_4_full_family_simultaneous.csv"
)
FIGURES = PROJECT / "figures" / "publication"
TABLES = PROJECT / "tables" / "publication"
PRINTED_WIDTH_IN = 7.08
MODEL_ORDER = (
    "weighted_graph",
    "graph_laplacian",
    "affine_sheaf",
    "affine_sheaf_permuted",
)
SCENARIOS = (
    "identity_limit",
    "heterogeneous_affine",
    "incompatible_cycles",
    "noisy_missing",
)


def main() -> int:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    metrics = pd.read_csv(RUN / "case_metrics.csv")
    contrasts = pd.read_csv(RUN / "paired_bootstrap_contrasts.csv")
    adjusted = pd.read_csv(ADJUSTED)
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)

    summary = (
        metrics.groupby("model", sort=False)
        .agg(
            n_cases=("seed", "count"),
            pr_auc_mean=("pr_auc", "mean"),
            brier_mean=("brier", "mean"),
            log_loss_mean=("log_loss", "mean"),
            selected_f1_mean=("selected_f1", "mean"),
        )
        .reset_index()
    )
    summary["model"] = summary["model"].map(
        {
            "weighted_graph": "Edge-local graph",
            "graph_laplacian": "Identity Laplacian",
            "affine_sheaf": "Affine sheaf",
            "affine_sheaf_permuted": "Permuted-map sheaf",
        }
    )
    summary = summary.rename(
        columns={
            "model": "Model",
            "n_cases": "Cases",
            "pr_auc_mean": "PR-AUC",
            "brier_mean": "Brier",
            "log_loss_mean": "Log loss",
            "selected_f1_mean": "Selected F1",
        }
    )
    summary.round(6).to_csv(
        TABLES / "table8_sheaf_vs_weighted_graph.csv", index=False, lineterminator="\n"
    )
    (TABLES / "table8_sheaf_vs_weighted_graph.md").write_text(
        summary.round(4).to_markdown(index=False) + "\n", encoding="utf-8"
    )

    contrasts.round(6).to_csv(
        TABLES / "tableS7_sheaf_vs_graph_contrasts.csv",
        index=False,
        lineterminator="\n",
    )
    (TABLES / "tableS7_sheaf_vs_graph_contrasts.md").write_text(
        contrasts.round(4).to_markdown(index=False) + "\n", encoding="utf-8"
    )

    system_summary = pd.read_csv(SYSTEM_RUN / "condition_summary.csv")
    system_contrasts = pd.read_csv(SYSTEM_RUN / "paired_bootstrap_contrasts.csv")
    system_export = pd.concat(
        [
            system_summary.assign(record_type="condition_summary"),
            system_contrasts.assign(record_type="paired_contrast"),
        ],
        ignore_index=True,
        sort=False,
    )
    leading = ["record_type"] + [
        column for column in system_export.columns if column != "record_type"
    ]
    system_export[leading].round(6).to_csv(
        TABLES / "tableS8_public_pipeline_acceptance.csv",
        index=False,
        lineterminator="\n",
    )
    system_markdown = (
        "### Condition summary\n\n"
        + system_summary.round(4).to_markdown(index=False)
        + "\n\n### Paired case-block contrasts\n\n"
        + system_contrasts.round(4).to_markdown(index=False)
        + "\n"
    )
    (TABLES / "tableS8_public_pipeline_acceptance.md").write_text(
        system_markdown, encoding="utf-8"
    )

    tick_labels = {
        "weighted_graph": "Edge-local\ngraph",
        "graph_laplacian": "Identity\nLaplacian",
        "affine_sheaf": "Affine\nsheaf",
        "affine_sheaf_permuted": "Permuted-map\nsheaf",
    }
    legend_labels = {
        "weighted_graph": "Edge-local graph",
        "graph_laplacian": "Identity Laplacian",
        "affine_sheaf": "Affine sheaf",
        "affine_sheaf_permuted": "Permuted-map sheaf",
    }
    colours = ["#8c8c8c", "#4c78a8", "#2a9d8f", "#e76f51"]
    fig, axes = plt.subplots(2, 2, figsize=(PRINTED_WIDTH_IN, 5.25))
    overall = metrics.groupby("model", sort=False).mean(numeric_only=True).reindex(MODEL_ORDER)
    x_models = np.arange(len(MODEL_ORDER))

    axes[0, 0].bar(x_models, overall["pr_auc"], color=colours, edgecolor="black", linewidth=0.5)
    axes[0, 0].set_xticks(x_models, [tick_labels[m] for m in MODEL_ORDER])
    axes[0, 0].set_ylabel("Mean case PR-AUC")
    axes[0, 0].set_title("a  Locked-test discrimination", loc="left", fontweight="bold")

    axes[0, 1].bar(x_models, overall["brier"], color=colours, edgecolor="black", linewidth=0.5)
    axes[0, 1].set_xticks(x_models, [tick_labels[m] for m in MODEL_ORDER])
    axes[0, 1].set_ylabel("Mean Brier score (lower is better)")
    axes[0, 1].set_title("b  Locked-test calibration", loc="left", fontweight="bold")

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
            label=legend_labels[model],
        )
    axes[1, 0].set_xticks(
        x,
        ["Identity", "Heterog.", "Conflict", "Noisy/miss."],
    )
    axes[1, 0].set_ylabel("Mean case PR-AUC")
    axes[1, 0].set_title("c  Scenario-specific performance", loc="left", fontweight="bold")
    axes[1, 0].legend(frameon=False, ncol=2, loc="upper center")

    selected = adjusted[
        (adjusted["scenario"] == "all")
        & (adjusted["left"] == "affine_sheaf")
        & adjusted["right"].isin(
            ["weighted_graph", "graph_laplacian", "affine_sheaf_permuted"]
        )
        & (adjusted["metric"] == "pr_auc")
    ].copy()
    selected["label"] = selected["right"].map(
        {
            "weighted_graph": "vs edge-local graph",
            "graph_laplacian": "vs identity Laplacian",
            "affine_sheaf_permuted": "vs permuted maps",
        }
    )
    y = np.arange(len(selected))
    means = selected["mean_difference"].to_numpy(float)
    errors = np.vstack(
        [
            means - selected["simultaneous_ci95_low"].to_numpy(float),
            selected["simultaneous_ci95_high"].to_numpy(float) - means,
        ]
    )
    axes[1, 1].axvline(0.0, color="black", linewidth=0.8)
    axes[1, 1].errorbar(means, y, xerr=errors, fmt="o", color="#2a9d8f", capsize=3)
    axes[1, 1].set_yticks(y, selected["label"])
    axes[1, 1].set_xlabel("Paired PR-AUC difference (sheaf - comparator)")
    axes[1, 1].set_title(
        "d  Family-wise simultaneous 95% intervals",
        loc="left",
        fontweight="bold",
    )

    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", alpha=0.18)
    fig.suptitle(
        "Incremental contribution of affine sheaf structure",
        fontsize=11,
        fontweight="bold",
    )
    fig.subplots_adjust(
        hspace=0.58, wspace=0.38, left=0.10, right=0.98, top=0.92, bottom=0.10
    )
    fig.savefig(FIGURES / "figure6_sheaf_vs_weighted_graph.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "figure6_sheaf_vs_weighted_graph.pdf", bbox_inches="tight")
    fig.savefig(
        FIGURES / "figure6_sheaf_vs_weighted_graph.tif",
        dpi=300,
        pil_kwargs={"compression": "tiff_lzw"},
        bbox_inches="tight",
    )
    plt.close(fig)

    manifest_path = FIGURES / "figure_source_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    if "Printed width" not in manifest:
        manifest["Printed width"] = ""
    if "Minimum label size" not in manifest:
        manifest["Minimum label size"] = ""
    manifest = manifest[manifest["Manuscript item"] != "Figure 6"]
    figure6_row = pd.DataFrame(
        [
            {
                "Manuscript item": "Figure 6",
                "Artifact stem": "figure6_sheaf_vs_weighted_graph",
                "Locked source": (
                    "RUN-M7-SHEAF-VS-GRAPH-20260729-01/case_metrics.csv; "
                    "paired_bootstrap_contrasts.csv; "
                    "post_review_audit_20260730/m7_4_full_family_simultaneous.csv"
                ),
                "Purpose": (
                    "Competence-matched affine-sheaf versus weighted-graph "
                    "ablation with family-wise intervals"
                ),
                "Printed width": f"{PRINTED_WIDTH_IN:.2f} in (180 mm)",
                "Minimum label size": "8.0 pt",
            }
        ]
    )
    pd.concat([manifest, figure6_row], ignore_index=True).to_csv(
        manifest_path, index=False, lineterminator="\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
