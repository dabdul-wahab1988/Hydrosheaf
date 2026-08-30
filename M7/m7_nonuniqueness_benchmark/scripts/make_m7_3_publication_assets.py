"""Build the manuscript-ready process-based integration figures and tables.

The graphics use a compact, colour-blind-safe style suitable for a two-column
Nature Portfolio manuscript. Every plotted value is read from the locked
process-based integration bundle or from the supporting field diagnostics.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK = REPO_ROOT / "M7" / "m7_nonuniqueness_benchmark"
RESULTS = BENCHMARK / "results" / "m7_3_locked"
FIGURES = BENCHMARK / "figures" / "publication"
TABLES = BENCHMARK / "tables" / "publication"
FIELD_RESULTS = (
    REPO_ROOT
    / "M7"
    / "m7_nonuniqueness_benchmark"
    / "results"
    / "supporting_validation"
)

WIDTH_TWO_COLUMN = 7.08
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
GRAY = "#6B6B6B"
LIGHT_GRAY = "#D9D9D9"
VERY_LIGHT = "#F3F3F3"
BLACK = "#222222"


def configure_style() -> None:
    """Apply one reproducible publication theme."""

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "normal",
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 6.8,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "grid.color": "#E6E6E6",
            "grid.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(ax: plt.Axes, label: str, *, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", bbox_inches="tight", dpi=600)
    fig.savefig(
        FIGURES / f"{stem}.tif",
        bbox_inches="tight",
        dpi=300,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / name)


def forest(
    ax: plt.Axes,
    frame: pd.DataFrame,
    *,
    labels: Sequence[str],
    colors: Sequence[str],
    xlabel: str,
    zero: bool = True,
) -> None:
    y = np.arange(len(frame))[::-1]
    estimate = frame["mean_difference"].to_numpy(float)
    low = frame["ci95_low"].to_numpy(float)
    high = frame["ci95_high"].to_numpy(float)
    if zero:
        ax.axvline(0.0, color=GRAY, lw=0.8, ls="--", zorder=0)
    for index, (yi, value, left, right, color) in enumerate(
        zip(y, estimate, low, high, colors)
    ):
        ax.errorbar(
            value,
            yi,
            xerr=[[value - left], [right - value]],
            fmt="o",
            ms=4.2,
            color=color,
            ecolor=color,
            elinewidth=1.1,
            capsize=2,
            zorder=3,
        )
        offset = 0.02 * max(1.0e-9, float(high.max() - low.min()))
        ax.text(
            right + offset,
            yi,
            f"{value:+.3g}",
            va="center",
            ha="left",
            fontsize=6.5,
            color=BLACK,
        )
    x_min = min(float(low.min()), 0.0)
    x_max = max(float(high.max()), 0.0)
    span = max(x_max - x_min, 1.0e-6)
    ax.set_xlim(x_min - 0.06 * span, x_max + 0.24 * span)
    ax.set_yticks(y, labels)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x")


def figure1_design() -> None:
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(WIDTH_TWO_COLUMN, 3.85),
        gridspec_kw={"height_ratios": [1.05, 1.0]},
    )
    for ax in axes:
        ax.set_axis_off()

    ax = axes[0]
    panel_label(ax, "a")
    ax.set_title(
        "Locked synthetic-truth benchmark",
        loc="left",
        pad=4,
    )
    boxes = [
        (0.02, "Official MODFLOW 6\n+ MODPATH 7", BLUE),
        (0.27, "Independent nonlinear\nchemistry + ³H/³⁹Ar", GREEN),
        (0.52, "Blind HydroSheaf\nH / A / C evidence", ORANGE),
        (0.77, "Locked test +\nnegative controls", PURPLE),
    ]
    for x, text, color in boxes:
        patch = FancyBboxPatch(
            (x, 0.32),
            0.19,
            0.38,
            boxstyle="round,pad=0.012,rounding_size=0.015",
            transform=ax.transAxes,
            fc=mpl.colors.to_rgba(color, 0.12),
            ec=color,
            lw=1.0,
        )
        ax.add_patch(patch)
        ax.text(
            x + 0.095,
            0.51,
            text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=7.2,
        )
    for start in (0.215, 0.465, 0.715):
        ax.add_patch(
            FancyArrowPatch(
                (start, 0.51),
                (start + 0.045, 0.51),
                transform=ax.transAxes,
                arrowstyle="-|>",
                mutation_scale=8,
                lw=0.8,
                color=GRAY,
            )
        )
    ax.text(
        0.02,
        0.12,
        "Development seeds freeze fusion coefficients; 5301–5312 remain untouched until scoring.",
        transform=ax.transAxes,
        color=GRAY,
        fontsize=6.7,
    )

    ax = axes[1]
    panel_label(ax, "b")
    ax.set_title(
        "Field interpretation is restricted by the Ghana evidence ceiling",
        loc="left",
        pad=4,
    )
    left = FancyBboxPatch(
        (0.03, 0.25),
        0.38,
        0.48,
        boxstyle="round,pad=0.015,rounding_size=0.015",
        transform=ax.transAxes,
        fc=mpl.colors.to_rgba(GREEN, 0.10),
        ec=GREEN,
        lw=1.0,
    )
    right = FancyBboxPatch(
        (0.59, 0.25),
        0.38,
        0.48,
        boxstyle="round,pad=0.015,rounding_size=0.015",
        transform=ax.transAxes,
        fc=mpl.colors.to_rgba(RED, 0.08),
        ec=RED,
        lw=1.0,
    )
    ax.add_patch(left)
    ax.add_patch(right)
    ax.text(
        0.22,
        0.64,
        "Supportable component diagnostics",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        0.22,
        0.43,
        "QC • seasonal chemistry • reaction equivalence\n"
        "measurement value • alternative-edge sensitivity",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=7,
    )
    ax.text(
        0.78,
        0.64,
        "Non-identifiable field claims",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        0.78,
        0.43,
        "Residence time • exact directed topology\n"
        "screen-resolved flow • unique reactions • field twin",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=7,
    )
    ax.add_patch(
        FancyArrowPatch(
            (0.425, 0.49),
            (0.575, 0.49),
            transform=ax.transAxes,
            arrowstyle="<|-|>",
            mutation_scale=8,
            lw=0.8,
            color=GRAY,
        )
    )
    ax.text(
        0.50,
        0.56,
        "claim boundary",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.5,
        color=GRAY,
    )
    fig.subplots_adjust(hspace=0.18, left=0.06, right=0.99, top=0.95, bottom=0.05)
    save_figure(fig, "figure1_benchmark_and_claim_design")


def figure2_evidence_integration() -> None:
    summary = read_csv("evidence_panel_summary.csv")
    contrasts = read_csv("evidence_case_bootstrap_contrasts.csv")
    native = summary[summary["condition"] == "native"].set_index("panel")
    panel_order = ["H", "A", "C", "HA", "HC", "AC", "HAC"]

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_TWO_COLUMN, 4.85))
    ax = axes[0, 0]
    colors = [GRAY, GRAY, GREEN, SKY, BLUE, BLUE, ORANGE]
    values = native.loc[panel_order, "pr_auc"]
    bars = ax.bar(panel_order, values, color=colors, width=0.72)
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=6.2,
        )
    ax.set_ylim(0, 0.56)
    ax.set_ylabel("PR-AUC")
    ax.set_title("Native evidence panels")
    ax.grid(axis="y")
    panel_label(ax, "a")

    selected = contrasts[
        (contrasts["condition"] == "native") & (contrasts["metric"] == "pr_auc")
    ].copy()
    order = [
        "native_incremental_chemistry",
        "native_incremental_hydraulics",
        "native_incremental_age",
    ]
    selected = selected.set_index("contrast").loc[order].reset_index()
    ax = axes[0, 1]
    forest(
        ax,
        selected,
        labels=["Chemistry added", "Hydraulics added", "Age added"],
        colors=[GREEN, BLUE, RED],
        xlabel="Change in PR-AUC",
    )
    ax.set_title("Incremental evidence value")
    panel_label(ax, "b")

    adverse_order = [
        "native_incremental_age",
        "permuted_age_increment",
        "permuted_hydraulic_increment",
        "joint_misspecification",
    ]
    labels = [
        "Native age",
        "Permuted age",
        "Permuted hydraulics",
        "Joint misspecification",
    ]
    colors = [ORANGE, PURPLE, BLUE, RED]
    def stress_rows(metric: str) -> pd.DataFrame:
        rows = []
        for contrast in adverse_order:
            rows.append(
                contrasts[
                    (contrasts["contrast"] == contrast)
                    & (contrasts["metric"] == metric)
                ].iloc[0]
            )
        return pd.DataFrame(rows)

    ax = axes[1, 0]
    forest(
        ax,
        stress_rows("mean_edge_entropy"),
        labels=labels,
        colors=colors,
        xlabel="Change in edge entropy\n(negative = more certain)",
    )
    ax.axvspan(ax.get_xlim()[0], 0, color=mpl.colors.to_rgba(ORANGE, 0.07), zorder=0)
    ax.set_title("Certainty shift")
    panel_label(ax, "c")

    ax = axes[1, 1]
    forest(
        ax,
        stress_rows("pr_auc"),
        labels=labels,
        colors=colors,
        xlabel="Change in PR-AUC\n(negative = degraded ranking)",
    )
    ax.axvspan(ax.get_xlim()[0], 0, color=mpl.colors.to_rgba(RED, 0.07), zorder=0)
    ax.set_yticklabels([])
    ax.set_title("Discrimination shift")
    panel_label(ax, "d")
    fig.subplots_adjust(
        hspace=0.58, wspace=0.42, left=0.17, right=0.99, top=0.94, bottom=0.11
    )
    save_figure(fig, "figure2_evidence_integration")


def figure3_topology_age() -> None:
    contrasts = read_csv("topology_age_bootstrap_contrasts.csv")
    sensitivity = read_csv("topology_age_sensitivity.csv")
    comparison_order = [
        ("informative", "correct_minus_no_topology"),
        ("tritium_only", "correct_minus_no_topology"),
        ("informative", "reversed_minus_correct"),
    ]
    labels = [
        "Correct − none\n³H + ³⁹Ar",
        "Correct − none\n³H only",
        "Reversed − correct\n³H + ³⁹Ar",
    ]

    def select_metric(metric: str) -> pd.DataFrame:
        rows = []
        for regime, contrast in comparison_order:
            rows.append(
                contrasts[
                    (contrasts["tracer_regime"] == regime)
                    & (contrasts["contrast"] == contrast)
                    & (contrasts["metric"] == metric)
                ].iloc[0]
            )
        return pd.DataFrame(rows)

    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_TWO_COLUMN, 4.95))
    ax = axes[0, 0]
    forest(
        ax,
        select_metric("age_mae_years"),
        labels=labels,
        colors=[GREEN, GREEN, RED],
        xlabel="Change in age MAE (years)",
    )
    ax.set_title("Accuracy")
    panel_label(ax, "a", x=0.0)

    ax = axes[0, 1]
    forest(
        ax,
        select_metric("mean_interval_width_years"),
        labels=labels,
        colors=[GREEN, GREEN, RED],
        xlabel="Change in 95% interval width (years)",
    )
    ax.set_title("Posterior precision")
    ax.set_yticklabels([])
    panel_label(ax, "b")

    ax = axes[1, 0]
    graph_order = ["none", "partial_true", "complete_true", "reversed"]
    graph_labels = ["None", "Partial true", "Complete true", "Reversed"]
    x = np.arange(len(graph_order))
    for offset, regime, color, marker, label in [
        (-0.12, "informative", BLUE, "o", "³H + ³⁹Ar"),
        (0.12, "tritium_only", ORANGE, "s", "³H only"),
    ]:
        subset = sensitivity[sensitivity["tracer_regime"] == regime]
        for index, condition in enumerate(graph_order):
            values = subset.loc[
                subset["graph_condition"] == condition, "importance_ess_fraction"
            ].to_numpy(float)
            jitter = np.linspace(-0.045, 0.045, len(values))
            ax.scatter(
                np.full(len(values), x[index] + offset) + jitter,
                values,
                s=12,
                color=color,
                marker=marker,
                alpha=0.72,
                edgecolors="none",
            )
            ax.plot(
                [x[index] + offset - 0.07, x[index] + offset + 0.07],
                [np.median(values), np.median(values)],
                color=BLACK,
                lw=1.1,
            )
    ax.axhline(400 / 50_000, color=RED, ls="--", lw=0.8)
    ax.text(3.42, 400 / 50_000 * 1.25, "ESS rule", color=RED, fontsize=6.2)
    ax.set_yscale("log")
    ax.set_ylim(0.0004, 1.5)
    ax.set_xticks(x, graph_labels, rotation=18, ha="right")
    ax.set_ylabel("Importance ESS fraction")
    ax.set_title("Numerical compatibility")
    ax.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color=BLUE, label="³H + ³⁹Ar"),
            Line2D([], [], marker="s", ls="", color=ORANGE, label="³H only"),
        ],
        frameon=False,
        loc="lower left",
    )
    ax.grid(axis="y", which="both")
    panel_label(ax, "c")

    ax = axes[1, 1]
    coverage = select_metric("age_95_coverage")
    forest(
        ax,
        coverage,
        labels=labels,
        colors=[GREEN, GREEN, RED],
        xlabel="Change in 95% coverage",
    )
    ax.text(
        0.98,
        0.08,
        "Reversed ³H-only estimate omitted:\nESS failed in 8/12 cases.",
        transform=ax.transAxes,
        fontsize=6.2,
        color=RED,
        va="bottom",
        ha="right",
    )
    ax.set_title("Calibration")
    ax.set_yticklabels([])
    panel_label(ax, "d")
    fig.subplots_adjust(
        hspace=0.52, wspace=0.32, left=0.15, right=0.99, top=0.93, bottom=0.12
    )
    save_figure(fig, "figure3_topology_conditions_age")


def figure4_reactions() -> None:
    summary = read_csv("reaction_nonuniqueness_summary.csv")
    summary = summary[summary["true_process"] != "ALL"].copy()
    order = [
        "carbonate_precipitation",
        "carbonate_weathering",
        "iron_reduction",
        "silicate_weathering",
        "denitrification",
        "sulfate_reduction",
    ]
    labels = [
        "Carbonate precipitation",
        "Carbonate weathering",
        "Iron reduction",
        "Silicate weathering",
        "Denitrification",
        "Sulfate reduction",
    ]
    core = summary[summary["tier"] == "core"].set_index("true_process").loc[order]
    enhanced = (
        summary[summary["tier"] == "enhanced"].set_index("true_process").loc[order]
    )
    y = np.arange(len(order))
    height = 0.34

    fig, axes = plt.subplots(
        2, 2, figsize=(WIDTH_TWO_COLUMN, 4.85), sharey=True
    )
    metrics = [
        (
            "modal_family_accuracy",
            "Modal-family recovery",
            "Accuracy",
            (0, 1.05),
        ),
        (
            "mean_true_family_probability",
            "Probability assigned to truth",
            "Probability",
            (0, 1.05),
        ),
        (
            "mean_family_support_entropy",
            "Support entropy",
            "Entropy (nats)",
            (0, 0.58),
        ),
        (
            "mean_effective_supported_families",
            "Effective support size",
            "Effective supported families",
            (0.9, 1.85),
        ),
    ]
    for index, (ax, (metric, title, xlabel, limits)) in enumerate(
        zip(axes.flat, metrics)
    ):
        ax.barh(y + height / 2, core[metric], height, color=GRAY, label="Core")
        ax.barh(
            y - height / 2,
            enhanced[metric],
            height,
            color=BLUE,
            label="Enhanced",
        )
        ax.set_xlim(*limits)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.grid(axis="x")
        panel_label(ax, chr(ord("a") + index))
        if index in (0, 2):
            ax.set_yticks(y, labels)
            ax.tick_params(axis="y", labelleft=True)
        else:
            ax.tick_params(axis="y", length=0, labelleft=False)
    axes[0, 0].invert_yaxis()
    axes[0, 0].legend(frameon=False, loc="upper right", ncol=2)
    axes[0, 0].text(
        0.03,
        0.87,
        "Carbonate:\n0/36 recovered",
        transform=axes[0, 0].transAxes,
        color=RED,
        fontsize=6.4,
        fontweight="bold",
        va="top",
    )
    axes[1, 0].text(
        0.03,
        0.93,
        "Low entropy can be\nconfidently wrong",
        transform=axes[1, 0].transAxes,
        color=RED,
        fontsize=6.5,
        va="top",
    )
    fig.subplots_adjust(
        hspace=0.48, wspace=0.30, left=0.25, right=0.99, top=0.94, bottom=0.11
    )
    save_figure(fig, "figure4_reaction_nonuniqueness")


def figure5_ghana_boundary() -> None:
    audit = json.loads(
        (RESULTS / "ghana_data_scope_audit.json").read_text(encoding="utf-8")
    )
    field_summary = pd.read_csv(FIELD_RESULTS / "field_prequential_summary.csv")
    field_audit = json.loads(
        (FIELD_RESULTS / "field_prequential_audit.json").read_text(encoding="utf-8")
    )

    availability = [
        ("Major chemistry", True),
        ("Stable water isotopes", True),
        ("Sr + SiO₂", True),
        ("Single static-water field", audit["single_static_water_level_available"]),
        (
            "Environmental age tracers",
            audit["environmental_age_tracer_panel_available"],
        ),
        ("Screen intervals", audit["screen_intervals_available"]),
        ("Repeated heads", audit["time_varying_head_series_available"]),
        (
            "Independent flow truth",
            audit["independent_field_connectivity_truth_available"],
        ),
    ]
    claims = [
        ("Data readiness / QC", True),
        ("Seasonal chemistry prediction", True),
        ("Reaction equivalence classes", True),
        ("Alternative-edge sensitivity", True),
        ("Residence-time validation", False),
        ("Exact directed topology", False),
        ("Unique field reactions", False),
        ("Validated field digital twin", False),
    ]

    fig, axes_grid = plt.subplots(2, 2, figsize=(WIDTH_TWO_COLUMN, 5.0))
    axes = list(axes_grid.flat)
    for panel, (ax, rows, title) in enumerate(
        [
            (axes[0], availability, "Observed evidence"),
            (axes[1], claims, "Defensible field claims"),
        ]
    ):
        rows = rows[::-1]
        names = [name for name, _ in rows]
        values = [value for _, value in rows]
        y = np.arange(len(rows))
        colors = [GREEN if value else RED for value in values]
        ax.barh(
            y, np.ones(len(rows)), color=[mpl.colors.to_rgba(c, 0.12) for c in colors]
        )
        for yi, name, value, color in zip(y, names, values, colors):
            ax.text(
                0.06,
                yi,
                "✓" if value else "—",
                color=color,
                va="center",
                ha="center",
                fontsize=9,
                fontweight="bold",
            )
            ax.text(
                0.14,
                yi,
                name,
                color=BLACK,
                va="center",
                ha="left",
                fontsize=6.4,
            )
        ax.set_yticks([])
        ax.set_xlim(0, 1)
        ax.set_xticks([])
        ax.set_title(title)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        panel_label(ax, chr(ord("a") + panel))

    ax = axes[2]
    coverage_labels = ["Wells", "Complete\npairs", "Seasonal\nobservations"]
    coverage_values = [
        audit["n_wells"],
        field_audit["n_complete_quantitative_pairs"],
        audit["n_hydrochemistry_rows"],
    ]
    bars = ax.bar(
        coverage_labels,
        coverage_values,
        color=[BLUE, SKY, GREEN],
        width=0.68,
    )
    for bar, value in zip(bars, coverage_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 7,
            str(value),
            ha="center",
            va="bottom",
            fontsize=6.4,
        )
    ax.set_ylim(0, 360)
    ax.set_ylabel("Count")
    ax.grid(axis="y")
    ax.set_title("Field-audit coverage")
    panel_label(ax, "c")

    ax = axes[3]
    overall = field_summary[field_summary["ion"] == "ALL"].set_index("method")
    methods = ["persistence", "expanding_mean_delta", "hydrosheaf_graph_ridge"]
    method_labels = ["Persist.", "Mean Δ", "Graph"]
    bars = ax.bar(
        method_labels,
        overall.loc[methods, "mae_log1p"],
        color=[GRAY, SKY, BLUE],
        width=0.7,
    )
    for bar, value in zip(bars, overall.loc[methods, "mae_log1p"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.012,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=5.8,
        )
    ax.set_ylabel("MAE")
    ax.set_ylim(0, 0.39)
    ax.tick_params(axis="x", rotation=10)
    ax.grid(axis="y")
    ax.set_title("Seasonal hold-forward")
    panel_label(ax, "d")
    contrast = field_audit["paired_block_bootstrap"][
        "graph_ridge_minus_expanding_mean_delta"
    ]
    ax.text(
        1.6,
        0.31,
        "Graph vs mean Δ:\n"
        f"{contrast['mean_paired_mae_difference_log1p']:+.3f} "
        f"[{contrast['ci95_low']:+.3f}, {contrast['ci95_high']:+.3f}]",
        fontsize=6.0,
        va="center",
        ha="center",
    )
    fig.subplots_adjust(
        hspace=0.48, wspace=0.34, left=0.08, right=0.99, top=0.94, bottom=0.10
    )
    save_figure(fig, "figure5_ghana_supportability_boundary_m7_only")


def markdown_table(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in frame.iterrows():
        values = ["" if pd.isna(value) else str(value) for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


# Stems that are authored by hand and must never be produced by a generator.
#
# `table1_benchmark_design` is the seven-audit design-and-claim map added in
# response to reviewer comment M7-20260729-R04 and cited as TAB-1 from both
# 01-introduction and 03-results.  This script previously wrote an unrelated
# three-column scale table to the same stem, so any re-run silently replaced the
# reviewer-requested table and the manuscript still assembled without error.
# `table2_primary_m7_3_decision` is likewise hand-authored and has no CSV.
HAND_MAINTAINED_STEMS = frozenset(
    {
        "table1_benchmark_design",
        "table2_primary_m7_3_decision",
    }
)


def write_table(frame: pd.DataFrame, stem: str, title: str) -> None:
    if stem in HAND_MAINTAINED_STEMS:
        raise RuntimeError(
            f"Refusing to overwrite hand-maintained table {stem!r}. "
            "It is authored by hand and referenced from the manuscript; "
            "a generator must write to a different stem."
        )
    TABLES.mkdir(parents=True, exist_ok=True)
    frame.to_csv(TABLES / f"{stem}.csv", index=False)
    (TABLES / f"{stem}.md").write_text(
        f"# {title}\n\n{markdown_table(frame)}\n",
        encoding="utf-8",
    )


def interval_text(row: pd.Series, digits: int = 3) -> str:
    return (
        f"{row['mean_difference']:.{digits}f} "
        f"[{row['ci95_low']:.{digits}f}, {row['ci95_high']:.{digits}f}]"
    )


def build_tables() -> None:
    manifest = json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8"))
    design = pd.DataFrame(
        [
            ("Development cases", len(manifest["development_seeds"]), "5201–5206"),
            ("Locked test cases", len(manifest["locked_test_seeds"]), "5301–5312"),
            ("Age particles", manifest["age_importance_particles"], "per case/regime"),
            (
                "Reaction bootstraps",
                manifest["reaction_bootstrap_per_case"],
                "per test case",
            ),
            (
                "Case-block bootstraps",
                manifest["paired_case_bootstrap"],
                "independent case resampling",
            ),
            ("Candidate recall", f"{manifest['candidate_recall']:.3f}", "locked test"),
        ],
        columns=["Design item", "Value", "Scope"],
    )
    # Not Table 1.  TAB-1 is the hand-authored seven-audit design-and-claim map
    # (reviewer M7-20260729-R04).  This is the M7.3 scale record derived from the
    # locked manifest; it is a provenance artifact and is not cited in the
    # manuscript body.
    write_table(
        design,
        "table1b_locked_benchmark_scale",
        "Table 1b | Locked M7.3 benchmark scale",
    )

    native = read_csv("evidence_panel_summary.csv")
    native = native[native["condition"] == "native"][
        [
            "panel",
            "pr_auc",
            "roc_auc",
            "brier",
            "log_loss",
            "mean_edge_entropy",
            "expected_calibration_error",
        ]
    ].copy()
    native.columns = [
        "Panel",
        "PR-AUC",
        "ROC-AUC",
        "Brier",
        "Log loss",
        "Edge entropy",
        "ECE",
    ]
    numeric = native.columns[1:]
    native[numeric] = native[numeric].map(lambda value: f"{value:.3f}")
    write_table(
        native,
        "table2_native_evidence_panels",
        "Table 2 | Native locked-test evidence-panel performance",
    )

    contrasts = read_csv("evidence_case_bootstrap_contrasts.csv")
    contrast_names = {
        "native_incremental_age": "Native age added",
        "native_incremental_chemistry": "Native chemistry added",
        "native_incremental_hydraulics": "Native hydraulics added",
        "permuted_age_increment": "Permuted age added",
        "permuted_hydraulic_increment": "Permuted hydraulics added",
        "joint_misspecification": "Age + hydraulics misspecified",
    }
    rows = []
    for contrast, label in contrast_names.items():
        selected = contrasts[contrasts["contrast"] == contrast].set_index("metric")
        rows.append(
            {
                "Contrast": label,
                "PR-AUC Δ [95% CI]": interval_text(selected.loc["pr_auc"]),
                "Brier Δ [95% CI]": interval_text(selected.loc["brier"], 4),
                "Log-loss Δ [95% CI]": interval_text(selected.loc["log_loss"], 4),
                "Entropy Δ [95% CI]": interval_text(
                    selected.loc["mean_edge_entropy"], 4
                ),
            }
        )
    write_table(
        pd.DataFrame(rows),
        "table3_evidence_contrasts",
        "Table 3 | Case-block evidence contrasts",
    )

    topology = read_csv("topology_age_bootstrap_contrasts.csv")
    topology = topology[
        topology["metric"].isin(
            ["age_mae_years", "age_95_coverage", "mean_interval_width_years"]
        )
    ].copy()
    topology["Tracer regime"] = topology["tracer_regime"].map(
        {"informative": "³H + ³⁹Ar", "tritium_only": "³H only"}
    )
    topology["Comparison"] = topology["contrast"].map(
        {
            "correct_minus_no_topology": "Complete true − none",
            "partial_minus_no_topology": "Partial true − none",
            "reversed_minus_correct": "Reversed − complete true",
        }
    )
    topology["Metric"] = topology["metric"].map(
        {
            "age_mae_years": "Age MAE (years)",
            "age_95_coverage": "95% coverage",
            "mean_interval_width_years": "95% interval width (years)",
        }
    )
    topology["Difference [95% CI]"] = topology.apply(interval_text, axis=1)
    write_table(
        topology[
            ["Tracer regime", "Comparison", "Metric", "Difference [95% CI]", "n_cases"]
        ].rename(columns={"n_cases": "Cases"}),
        "table4_topology_age_contrasts",
        "Table 4 | Topology-conditioned age contrasts",
    )

    reactions = read_csv("reaction_nonuniqueness_summary.csv")
    reactions = reactions[reactions["true_process"] != "ALL"].copy()
    reactions["Panel"] = reactions["tier"].map({"core": "Core", "enhanced": "Enhanced"})
    reactions["Process"] = reactions["true_process"].map(
        {
            "carbonate_precipitation": "Carbonate precip.",
            "carbonate_weathering": "Carbonate weather.",
            "denitrification": "Denitrif.",
            "iron_reduction": "Iron reduction",
            "silicate_weathering": "Silicate weather.",
            "sulfate_reduction": "Sulfate reduction",
        }
    )
    for column in [
        "modal_family_accuracy",
        "mean_true_family_probability",
        "mean_family_support_entropy",
        "mean_effective_supported_families",
    ]:
        reactions[column] = reactions[column].map(lambda value: f"{value:.3f}")
    reactions = reactions.rename(
        columns={
            "n_edges": "Edges",
            "modal_family_accuracy": "Modal accuracy",
            "mean_true_family_probability": "True-family probability",
            "mean_family_support_entropy": "Support entropy",
            "mean_effective_supported_families": "Effective families",
        }
    )
    write_table(
        reactions[
            [
                "Panel",
                "Process",
                "Edges",
                "Modal accuracy",
                "True-family probability",
                "Support entropy",
                "Effective families",
            ]
        ],
        "table5_reaction_nonuniqueness",
        "Table 5 | Reaction-family recovery and non-uniqueness",
    )

    audit = json.loads(
        (RESULTS / "ghana_data_scope_audit.json").read_text(encoding="utf-8")
    )
    scope = pd.DataFrame(
        [
            ("Major hydrochemistry", "Available", "Component inference and QC"),
            ("Stable water isotopes", "Available", "Recharge/source evidence"),
            ("Environmental age tracers", "Absent", "Residence time non-identifiable"),
            ("Screen intervals", "Absent", "Vertical connectivity non-identifiable"),
            ("Repeated head series", "Absent", "Dynamic head validation unavailable"),
            ("Coordinates", "Masked", "No site-scale connectivity truth"),
            ("Processed graph edges", "Absent", "Edge sets are self-generated, not supplied"),
            ("Intra-season sampling dates", "Absent", "Batch order is disclosed and arbitrary"),
            ("Independent aquifer-type classification", "Absent", "Stratified reporting uses region instead"),
            ("Independent reaction truth", "Absent", "Unique mechanisms unvalidated"),
        ],
        columns=["Evidence", "Status", "Defensible use"],
    )
    assert audit["n_wells"] == 160
    write_table(
        scope,
        "table6_ghana_claim_boundary",
        "Table 6 | Ghana data scope and claim boundary",
    )

    full_tables = {
        "tableS1_all_evidence_conditions": read_csv("evidence_panel_summary.csv"),
        "tableS2_case_block_contrasts": read_csv(
            "evidence_case_bootstrap_contrasts.csv"
        ),
        "tableS3_topology_age_sensitivity": read_csv("topology_age_sensitivity.csv"),
        "tableS4_reaction_edge_nonuniqueness": read_csv(
            "reaction_edge_nonuniqueness.csv"
        ),
        "tableS5_conflict_diagnostics": read_csv("evidence_conflict_summary.csv"),
    }
    for stem, frame in full_tables.items():
        title = stem.replace("_", " ").replace("tableS", "Table S").title()
        write_table(frame, stem, title)


def write_manifest() -> None:
    rows = [
        (
            "Figure 1",
            "figure1_benchmark_and_claim_design",
            "Protocol; Ghana data-scope audit",
            "Benchmark architecture and field claim boundary",
        ),
        (
            "Figure 2",
            "figure2_evidence_integration",
            "evidence_panel_summary.csv; evidence_case_bootstrap_contrasts.csv",
            "Complementarity, redundancy and false confidence",
        ),
        (
            "Figure 3",
            "figure3_topology_conditions_age",
            "topology_age_sensitivity.csv; topology_age_bootstrap_contrasts.csv",
            "Topology-conditioned age accuracy, precision and ESS",
        ),
        (
            "Figure 4",
            "figure4_reaction_nonuniqueness",
            "reaction_nonuniqueness_summary.csv",
            "Process-family recovery and carbonate failure",
        ),
        (
            "Figure 5",
            "figure5_ghana_supportability_boundary_m7_only",
            "Ghana audit; supporting prequential audit",
            "Objective 6 supportability and non-identifiability",
        ),
    ]
    frame = pd.DataFrame(
        rows,
        columns=["Manuscript item", "Artifact stem", "Locked source", "Purpose"],
    )
    FIGURES.mkdir(parents=True, exist_ok=True)
    frame.to_csv(FIGURES / "figure_source_manifest.csv", index=False)


def verify_outputs() -> None:
    for stem in [
        "figure1_benchmark_and_claim_design",
        "figure2_evidence_integration",
        "figure3_topology_conditions_age",
        "figure4_reaction_nonuniqueness",
        "figure5_ghana_supportability_boundary_m7_only",
    ]:
        for suffix in (".pdf", ".png", ".tif"):
            path = FIGURES / f"{stem}{suffix}"
            if not path.exists() or path.stat().st_size == 0:
                raise RuntimeError(f"Missing publication figure: {path}")
    # Verify the M7.3 contract by name rather than by total directory count:
    # later locked extensions (M7.4 and system-acceptance tables) legitimately
    # share this publication directory.
    expected_tables = [f"table{index}" for index in range(1, 8)] + [
        f"tableS{index}" for index in range(1, 7)
    ]
    for stem in expected_tables:
        path = TABLES / f"{stem}.csv"
        matching = list(TABLES.glob(f"{stem}_*.csv"))
        if not path.exists() and not matching:
            raise RuntimeError(f"Missing required M7.3 publication table: {stem}")


def main() -> None:
    configure_style()
    figure1_design()
    figure2_evidence_integration()
    figure3_topology_age()
    figure4_reactions()
    figure5_ghana_boundary()
    build_tables()
    write_manifest()
    verify_outputs()
    print(f"M7.3 publication figures -> {FIGURES}")
    print(f"M7.3 publication tables -> {TABLES}")


if __name__ == "__main__":
    main()
