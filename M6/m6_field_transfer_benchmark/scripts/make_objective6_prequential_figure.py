"""Create Q1 Figure 5 from the locked supporting field hold-forward audit."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE = (
    REPO_ROOT
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "results"
    / "supporting_validation"
)
OUT = REPO_ROOT / "M6" / "m6_field_transfer_benchmark" / "figures" / "r_publication"

BLUE = "#0072B2"
SKY = "#56B4E9"
GRAY = "#6B6B6B"
GREEN = "#009E73"
RED = "#D55E00"
BLACK = "#222222"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 6.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "grid.color": "#E6E6E6",
            "grid.linewidth": 0.5,
            "pdf.fonttype": 42,
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )


def main() -> None:
    configure_style()
    summary = pd.read_csv(SOURCE / "field_prequential_summary.csv")
    audit = json.loads(
        (SOURCE / "field_prequential_audit.json").read_text(encoding="utf-8")
    )
    overall = summary[summary["ion"] == "ALL"].set_index("method")
    methods = ["persistence", "expanding_mean_delta", "hydrosheaf_graph_ridge"]
    labels = ["Persistence", "Expanding mean Δ", "Graph ridge"]
    colors = [GRAY, SKY, BLUE]

    fig, axes = plt.subplots(2, 2, figsize=(7.08, 4.55))
    ax = axes[0, 0]
    bars = ax.bar(labels, overall.loc[methods, "mae_log1p"], color=colors, width=0.68)
    for bar, value in zip(bars, overall.loc[methods, "mae_log1p"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.009,
            f"{value:.3f}",
            ha="center",
            fontsize=6.5,
        )
    ax.set_ylim(0, 0.39)
    ax.set_ylabel("MAE (log1p concentration)")
    ax.set_title("One-step seasonal hold-forward")
    ax.tick_params(axis="x", rotation=16)
    ax.grid(axis="y")
    panel_label(ax, "a")

    ax = axes[0, 1]
    bootstrap = audit["paired_block_bootstrap"]
    entries = [
        (
            "Graph − persistence",
            bootstrap["graph_ridge_minus_persistence"],
            GREEN,
        ),
        (
            "Graph − mean Δ",
            bootstrap["graph_ridge_minus_expanding_mean_delta"],
            BLUE,
        ),
    ]
    for y, (label, row, color) in enumerate(entries[::-1]):
        value = row["mean_paired_mae_difference_log1p"]
        low = row["ci95_low"]
        high = row["ci95_high"]
        ax.errorbar(
            value,
            y,
            xerr=[[value - low], [high - value]],
            fmt="o",
            color=color,
            capsize=2,
            lw=1.1,
        )
    ax.axvline(0, color=GRAY, ls="--", lw=0.8)
    ax.set_yticks([0, 1], [entries[1][0], entries[0][0]])
    ax.set_xlabel("Paired MAE difference\n(negative favours graph ridge)")
    ax.set_title("Well-block uncertainty")
    ax.grid(axis="x")
    panel_label(ax, "b")

    ax = axes[1, 0]
    ion = summary[summary["ion"] != "ALL"].copy()
    graph = ion[ion["method"] == "hydrosheaf_graph_ridge"].set_index("ion")
    baseline = ion[ion["method"] == "expanding_mean_delta"].set_index("ion")
    common = sorted(set(graph.index) & set(baseline.index))
    delta = graph.loc[common, "mae_log1p"] - baseline.loc[common, "mae_log1p"]
    display = [value.replace("_mg_L", "").replace("SiO2", "SiO₂") for value in common]
    order = np.argsort(delta.to_numpy())
    values = delta.to_numpy()[order]
    names = np.asarray(display)[order]
    ax.barh(
        np.arange(len(values)),
        values,
        color=[GREEN if value < 0 else RED for value in values],
    )
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_yticks(np.arange(len(values)), names)
    ax.set_xlabel("Graph − expanding mean Δ MAE")
    ax.set_title("Ion-level descriptive differences")
    ax.grid(axis="x")
    panel_label(ax, "c")

    ax = axes[1, 1]
    coverage = overall.loc[methods, "coverage90"]
    bars = ax.bar(labels, coverage, color=colors, width=0.68)
    ax.axhline(0.90, color=BLACK, ls="--", lw=0.8, label="Nominal 90%")
    for bar, value in zip(bars, coverage):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            f"{value:.2f}",
            ha="center",
            fontsize=6.5,
        )
    ax.set_ylim(0.55, 1.0)
    ax.set_ylabel("Empirical 90% coverage")
    ax.set_title("Predictive interval calibration")
    ax.tick_params(axis="x", rotation=16)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="y")
    panel_label(ax, "d")

    fig.text(
        0.01,
        0.01,
        "Within-campaign chemistry prediction only; coordinates are masked and no future wet-season labels were used.",
        fontsize=6.2,
        color=GRAY,
    )
    fig.subplots_adjust(
        hspace=0.55, wspace=0.55, left=0.15, right=0.99, top=0.93, bottom=0.14
    )
    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / "figure5_field_prequential"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=600)
    fig.savefig(
        stem.with_suffix(".tif"),
        bbox_inches="tight",
        dpi=300,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    print(f"Objective 6 field figure -> {stem}")


if __name__ == "__main__":
    main()
