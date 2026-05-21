"""Create Nature-style M4 figures from benchmark result tables."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Dict, Iterable, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = BENCHMARK_ROOT / "results"
FIGURE_DIR = BENCHMARK_ROOT / "figures" / "Manuscript_Ready"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
GRAY = "#666666"
LIGHT_GRAY = "#D9D9D9"

SCENARIO_LABELS = {
    "well_resolved_chain": "Well-resolved\nchain",
    "false_positive_shortcut": "False-positive\nshortcut",
    "false_negative_missing_link": "False-negative\nlink",
    "scale_mismatch_shortcuts": "Scale-mismatch\nshortcuts",
    "2.6.1_spatial_only": "Spatial only",
    "2.6.2_head_constrained": "Head-constrained",
    "2.6.3_hydrostratigraphic": "Hydrostratigraphic",
}


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _save(fig: plt.Figure, name: str) -> None:
    for suffix in (".png", ".svg"):
        out = FIGURE_DIR / f"{name}{suffix}"
        tmp = out.with_name(f"{out.stem}.tmp{out.suffix}")
        fig.savefig(tmp, dpi=300, bbox_inches="tight")
        try:
            if out.exists():
                out.unlink()
            tmp.replace(out)
        except Exception:
            time.sleep(0.5)
            alt = out.with_name(out.stem + ".new" + out.suffix)
            tmp.replace(alt)
        finally:
            if tmp.exists():
                try:
                    tmp.unlink()
                except Exception:
                    pass
    plt.close(fig)


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
        ha="left",
    )


def _display_scenario(value: str) -> str:
    return SCENARIO_LABELS.get(value, value.replace("_", " "))


def _draw_box(ax: plt.Axes, xy: Tuple[float, float], text: str, width: float, height: float, color: str) -> None:
    rect = Rectangle(xy, width, height, facecolor="white", edgecolor=color, linewidth=1.0)
    ax.add_patch(rect)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=6,
        color="#222222",
        wrap=True,
    )


def _arrow(ax: plt.Axes, start: Tuple[float, float], end: Tuple[float, float], color: str = GRAY, dashed: bool = False) -> None:
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=0.8,
        linestyle=(0, (3, 2)) if dashed else "solid",
        color=color,
    )
    ax.add_patch(patch)


def make_figure1_workflow() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")

    boxes = {
        "archive": ((0.3, 2.4), "MODPATH\narchive outputs", BLUE),
        "reference": ((2.1, 2.4), "Reference\ndirected graph", GREEN),
        "inference": ((4.0, 2.4), "Independent\nHydrosheaf inference", SKY),
        "metrics": ((6.1, 2.4), "Edge metrics\nand scale checks", PURPLE),
        "figures": ((8.1, 2.4), "Figure/table\nevidence package", GRAY),
        "priors": ((4.0, 0.7), "MODPATH-informed\nHydrosheaf priors", ORANGE),
        "guardrail": ((6.1, 0.7), "Prior mode:\nnot independent\nvalidation", VERMILLION),
    }
    for xy, text, color in boxes.values():
        _draw_box(ax, xy, text, 1.55, 0.85, color)

    _arrow(ax, (1.75, 2.83), (2.1, 2.83), BLUE)
    _arrow(ax, (3.55, 2.83), (4.0, 2.83), GREEN)
    _arrow(ax, (5.45, 2.83), (6.1, 2.83), SKY)
    _arrow(ax, (7.55, 2.83), (8.1, 2.83), PURPLE)
    _arrow(ax, (2.8, 2.4), (4.0, 1.2), ORANGE, dashed=True)
    _arrow(ax, (5.45, 1.12), (6.1, 1.12), VERMILLION, dashed=True)
    _arrow(ax, (7.55, 1.12), (8.8, 2.4), VERMILLION, dashed=True)

    ax.text(
        0.3,
        0.15,
        "Dashed branch documents MODPATH prior use separately from independent validation.",
        fontsize=6,
        color="#444444",
    )
    _save(fig, "Manuscript_Fig1_M4_Benchmark_Workflow")


def make_figure2_performance() -> None:
    df = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    if df.empty:
        return
    metric_cols = ["precision", "recall", "f1", "false_positive_rate", "false_negative_rate"]
    count_cols = ["tp", "fp", "fn", "tn"]
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(7.2, 2.8),
        gridspec_kw={"width_ratios": [1.35, 1.0, 1.0], "wspace": 0.65},
    )

    ax = axes[0]
    matrix = df[metric_cols].to_numpy(float)
    im = ax.imshow(matrix, vmin=0, vmax=1, cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(metric_cols)))
    ax.set_xticklabels(["Precision", "Recall", "F1", "FPR", "FNR"], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels([_display_scenario(value) for value in df["scenario"]])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color="white")
    ax.set_title("Agreement metrics")
    _panel_label(ax, "a")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=5.5)

    ax = axes[1]
    left = np.zeros(len(df))
    colors = {"tp": GREEN, "fp": VERMILLION, "fn": ORANGE, "tn": LIGHT_GRAY}
    for column in count_cols:
        values = pd.to_numeric(df[column], errors="coerce").fillna(0).to_numpy(float)
        ax.barh(np.arange(len(df)), values, left=left, color=colors[column], label=column.upper(), height=0.55)
        left += values
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel("Directed edges")
    ax.legend(frameon=False, ncol=2, loc="lower right", handlelength=1.0)
    ax.set_title("Edge classes")
    _panel_label(ax, "b")

    ax = axes[2]
    x = np.arange(len(df))
    width = 0.35
    ref = pd.to_numeric(df["median_reference_length"], errors="coerce")
    inf = pd.to_numeric(df["median_inferred_length"], errors="coerce")
    ax.bar(x - width / 2, ref, width, color=GRAY, label="Reference")
    ax.bar(x + width / 2, inf, width, color=BLUE, label="Inferred")
    scale = df["scale_mismatch"].astype(str).str.lower().eq("true")
    for i, is_mismatch in enumerate(scale):
        if is_mismatch:
            ymax = max(float(ref.iloc[i]), float(inf.iloc[i]))
            ax.plot(i, ymax * 1.08, marker="v", color=VERMILLION, markersize=4, clip_on=False)
    ax.set_xticks(x)
    ax.set_xticklabels([_display_scenario(value) for value in df["scenario"]], rotation=35, ha="right")
    ax.set_ylabel("Median edge length")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Scale check")
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
    _panel_label(ax, "c")

    _save(fig, "Manuscript_Fig2_M4_Independent_Topology_Performance")


def _edge_nodes(edge: str) -> Tuple[str, str]:
    left, right = edge.split("->", 1)
    return left, right


def _draw_edge(ax: plt.Axes, u: str, v: str, color: str, linestyle: str, label: str = "") -> None:
    positions = {"R": (0.0, 0.18), "A": (1.0, 0.42), "B": (2.0, 0.18), "C": (3.0, 0.42), "D": (4.0, 0.18)}
    x1, y1 = positions[u]
    x2, y2 = positions[v]
    distance = abs(x2 - x1)
    rad = 0.0 if distance <= 1.1 else 0.18
    patch = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=1.0,
        linestyle=linestyle,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(patch)


def make_figure3_edge_networks() -> None:
    edges = _read_csv(RESULT_DIR / "edge_classification.csv")
    if edges.empty:
        return
    scenarios = list(edges["scenario"].drop_duplicates())
    fig, axes = plt.subplots(1, len(scenarios), figsize=(7.2, 1.55), sharex=True, sharey=True)
    if len(scenarios) == 1:
        axes = [axes]
    positions = {"R": (0.0, 0.18), "A": (1.0, 0.42), "B": (2.0, 0.18), "C": (3.0, 0.42), "D": (4.0, 0.18)}
    style = {
        "TP": (GREEN, "solid"),
        "FP": (VERMILLION, "solid"),
        "FN": (ORANGE, (0, (3, 2))),
    }
    for ax, scenario in zip(axes, scenarios):
        subset = edges[edges["scenario"] == scenario]
        ax.set_xlim(-0.25, 4.25)
        ax.set_ylim(-0.05, 0.68)
        ax.axis("off")
        for _, row in subset.iterrows():
            label = str(row["classification"])
            if label == "TN":
                continue
            u, v = _edge_nodes(str(row["edge"]))
            color, linestyle = style[label]
            _draw_edge(ax, u, v, color, linestyle)
        for node, (x, y) in positions.items():
            ax.scatter(x, y, s=55, facecolor="white", edgecolor="#222222", linewidth=0.8, zorder=3)
            ax.text(x, y, node, ha="center", va="center", fontsize=6.5, zorder=4)
        ax.set_title(_display_scenario(scenario), pad=2)
    legend = [
        Line2D([0], [0], color=GREEN, lw=1.2, label="TP"),
        Line2D([0], [0], color=VERMILLION, lw=1.2, label="FP"),
        Line2D([0], [0], color=ORANGE, lw=1.2, linestyle=(0, (3, 2)), label="FN"),
    ]
    fig.legend(handles=legend, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.00))
    _save(fig, "Manuscript_Fig3_M4_Edge_Failure_Networks")


def make_figure4_endpoint_validation() -> None:
    endpoint = _read_csv(RESULT_DIR / "m4_topology_benchmark_summary.csv")
    sensitivity = _read_csv(RESULT_DIR / "m4_sparsity_sensitivity.csv")
    if endpoint.empty and sensitivity.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.75), gridspec_kw={"wspace": 0.42})
    ax = axes[0]
    if endpoint.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Endpoint validation output missing", ha="center", va="center")
    else:
        cols = [column for column in ("precision", "recall", "f1") if column in endpoint.columns]
        x = np.arange(len(endpoint))
        width = 0.22
        colors = [BLUE, GREEN, ORANGE]
        offsets = np.linspace(-width, width, len(cols))
        for offset, column, color in zip(offsets, cols, colors):
            values = pd.to_numeric(endpoint[column], errors="coerce").fillna(0)
            ax.bar(x + offset, values, width=width, color=color, label=column.upper())
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Score")
        ax.set_xticks(x)
        ax.set_xticklabels([_display_scenario(value) for value in endpoint["scenario"]], rotation=25, ha="right")
        ax.legend(frameon=False, loc="upper left")
        ax.set_title("Endpoint-derived MODPATH agreement")
    _panel_label(ax, "a")

    ax = axes[1]
    if sensitivity.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "Sparsity output missing", ha="center", va="center")
    else:
        x = pd.to_numeric(sensitivity["node_fraction"], errors="coerce")
        f1 = pd.to_numeric(sensitivity["mean_f1"], errors="coerce")
        f1_std = pd.to_numeric(sensitivity.get("std_f1", 0), errors="coerce").fillna(0)
        recall = pd.to_numeric(sensitivity["mean_recall"], errors="coerce")
        ax.errorbar(x, f1, yerr=f1_std, marker="o", color=BLUE, linewidth=1.0, markersize=3.5, capsize=2, label="F1")
        ax.plot(x, recall, marker="s", color=GREEN, linewidth=1.0, markersize=3.5, label="Recall")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Available node fraction")
        ax.set_ylabel("Score")
        ax.legend(frameon=False, loc="lower right")
        ax.set_title("Node-sparsity sensitivity")
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=5))
    _panel_label(ax, "b")
    _save(fig, "Manuscript_Fig4_M4_MODPATH_Endpoint_Validation")


def make_figure5_external_archive_validation() -> None:
    summary = _read_csv(RESULT_DIR / "external_modpath_archive_summary.csv")
    rank = _read_csv(RESULT_DIR / "external_modpath_travel_time_rank.csv")
    rank_summary = _read_csv(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv")
    structure = _read_csv(RESULT_DIR / "external_modpath_pathline_structure.csv")
    if summary.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.0), gridspec_kw={"wspace": 0.45, "hspace": 0.48})
    top = summary.iloc[0]

    ax = axes[0, 0]
    counts = [
        float(top.get("true_positive_edges", 0)),
        float(top.get("false_positive_edges", 0)),
        float(top.get("false_negative_edges", 0)),
    ]
    ax.bar(["TP", "FP", "FN"], counts, color=[GREEN, VERMILLION, ORANGE], width=0.62)
    ax.set_ylabel("Directed edges")
    ax.set_title("External MODPATH edge agreement")
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
    _panel_label(ax, "a")

    ax = axes[0, 1]
    if not rank.empty:
        x = pd.to_numeric(rank["endpoint_travel_time_mean"], errors="coerce")
        y = pd.to_numeric(rank["pathline_elapsed_time_median"], errors="coerce")
        valid = (x > 0) & (y > 0)
        ax.scatter(x[valid], y[valid], s=12, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.25)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Endpoint travel time mean")
        ax.set_ylabel("Pathline elapsed time median")
        rho = rank_summary.iloc[0].get("spearman_rho", np.nan) if not rank_summary.empty else np.nan
        ax.text(
            0.96,
            0.96,
            f"rho = {float(rho):.2f}",
            transform=ax.transAxes,
            va="top",
            ha="right",
            fontsize=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5},
        )
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Travel-time rank output missing", ha="center", va="center")
    ax.set_title("Travel-time rank diagnostic")
    _panel_label(ax, "b")

    ax = axes[1, 0]
    if not structure.empty:
        cells = pd.to_numeric(structure["n_compressed_cells"], errors="coerce").dropna()
        ax.hist(cells, bins=24, color=GRAY, edgecolor="white", linewidth=0.35)
        ax.axvline(cells.median(), color=VERMILLION, linewidth=1.0, label="Median")
        ax.set_xlabel("Compressed cells per pathline")
        ax.set_ylabel("Particles")
        ax.legend(frameon=False, loc="upper right")
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Pathline structure output missing", ha="center", va="center")
    ax.set_title("Full pathline structure")
    _panel_label(ax, "c")

    ax = axes[1, 1]
    metrics = [
        ("Edge F1", "edge_f1"),
        ("Direction", "direction_agreement_rate"),
        ("Source-receptor", "mean_source_receptor_overlap"),
        ("Endpoint projection", "endpoint_projection_preservation_rate"),
    ]
    values = [float(top.get(column, np.nan)) for _, column in metrics]
    ax.barh(np.arange(len(metrics)), values, color=[GREEN, SKY, BLUE, PURPLE], height=0.55)
    ax.set_yticks(np.arange(len(metrics)))
    ax.set_yticklabels([label for label, _ in metrics])
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Agreement")
    ax.set_title("Projected graph evidence")
    _panel_label(ax, "d")

    _save(fig, "Manuscript_Fig5_M4_External_MODPATH_Archive_Validation")


def _hull_xy(x: pd.Series, y: pd.Series) -> np.ndarray:
    from scipy.spatial import ConvexHull, QhullError

    points = np.column_stack([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")])
    points = points[np.isfinite(points).all(axis=1)]
    points = np.unique(points, axis=0)
    if len(points) < 3:
        return np.empty((0, 2))
    try:
        hull = ConvexHull(points)
    except QhullError:
        return np.empty((0, 2))
    return points[hull.vertices]


def make_figure6_capture_travel_time_validation() -> None:
    structure = _read_csv(RESULT_DIR / "external_modpath_pathline_structure.csv")
    capture = _read_csv(RESULT_DIR / "external_modpath_capture_envelope_overlap.csv")
    time = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time.csv")
    time_summary = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv")
    if structure.empty and capture.empty and time.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.1), gridspec_kw={"wspace": 0.45, "hspace": 0.52})

    ax = axes[0, 0]
    if not structure.empty:
        targets = list(structure["endpoint_target_node"].dropna().drop_duplicates())
        colors = [GREEN, BLUE, PURPLE, ORANGE, SKY]
        for idx, target in enumerate(targets[:5]):
            subset = structure[structure["endpoint_target_node"] == target]
            color = colors[idx % len(colors)]
            ax.scatter(
                pd.to_numeric(subset["endpoint_start_x"], errors="coerce"),
                pd.to_numeric(subset["endpoint_start_y"], errors="coerce"),
                s=5,
                color=color,
                alpha=0.22,
                linewidth=0,
            )
            hull = _hull_xy(subset["endpoint_start_x"], subset["endpoint_start_y"])
            if len(hull):
                closed = np.vstack([hull, hull[0]])
                ax.plot(closed[:, 0], closed[:, 1], color=color, linewidth=0.9, label=target.replace("cell_", ""))
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Start x")
        ax.set_ylabel("Start y")
        ax.legend(frameon=False, title="Target cell", loc="upper left", fontsize=5.5, title_fontsize=5.5)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Capture-envelope output missing", ha="center", va="center")
    ax.set_title("Capture-envelope point clouds")
    _panel_label(ax, "a")

    ax = axes[0, 1]
    if not capture.empty:
        labels = [str(value).replace("cell_", "") for value in capture["target_node"]]
        x = np.arange(len(capture))
        ax.bar(x, pd.to_numeric(capture["capture_envelope_iou"], errors="coerce"), color=GREEN, width=0.55)
        ax.plot(x, pd.to_numeric(capture["source_cell_jaccard"], errors="coerce"), color=BLUE, marker="o", linewidth=1.0, label="Source cells")
        ax.set_ylim(0, 1.05)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel("Overlap")
        ax.legend(frameon=False, loc="lower right")
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Capture overlap output missing", ha="center", va="center")
    ax.set_title("Capture-envelope overlap")
    _panel_label(ax, "b")

    ax = axes[1, 0]
    if not time.empty:
        x = pd.to_numeric(time["particle_endpoint_time_median"], errors="coerce")
        y = pd.to_numeric(time["hydrosheaf_modpath_weight_mean"], errors="coerce")
        valid = (x > 0) & (y > 0)
        ax.scatter(x[valid], y[valid], s=13, color=BLUE, alpha=0.76, edgecolor="white", linewidth=0.25)
        lo = float(min(x[valid].min(), y[valid].min()))
        hi = float(max(x[valid].max(), y[valid].max()))
        ax.plot([lo, hi], [lo, hi], color=GRAY, linewidth=0.9, linestyle=(0, (3, 2)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Endpoint particle time median")
        ax.set_ylabel("Hydrosheaf MODPATH weight")
        rho = time_summary.iloc[0].get("spearman_rho", np.nan) if not time_summary.empty else np.nan
        ax.text(0.05, 0.95, f"rho = {float(rho):.2f}", transform=ax.transAxes, va="top", fontsize=6)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Harmonized travel-time output missing", ha="center", va="center")
    ax.set_title("Harmonized travel-time agreement")
    _panel_label(ax, "c")

    ax = axes[1, 1]
    if not time.empty:
        x = pd.to_numeric(time["particle_endpoint_time_median"], errors="coerce")
        y = pd.to_numeric(time["hydrosheaf_modpath_weight_mean"], errors="coerce")
        valid = (x > 0) & (y > 0)
        diff = np.log10(y[valid]) - np.log10(x[valid])
        ax.hist(diff, bins=22, color=GRAY, edgecolor="white", linewidth=0.35)
        ax.axvline(0.0, color=VERMILLION, linewidth=1.0)
        ax.set_xlabel("log10(weight / endpoint time)")
        ax.set_ylabel("Edges")
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Harmonized travel-time output missing", ha="center", va="center")
    ax.set_title("Travel-time scale residual")
    _panel_label(ax, "d")

    _save(fig, "Manuscript_Fig6_M4_Capture_Travel_Time_Validation")


def make_edfigure1_prior_audit() -> None:
    priors = _read_csv(RESULT_DIR / "modpath_informed_priors.csv")
    if priors.empty:
        return
    fig, ax = plt.subplots(figsize=(3.5, 2.2))
    modes = priors["prior_mode"].astype(str)
    x = np.arange(len(priors))
    width = 0.24
    ax.bar(x - width, pd.to_numeric(priors["n_input_hydrosheaf_edges"], errors="coerce"), width, color=GRAY, label="Input Hydrosheaf")
    ax.bar(x, pd.to_numeric(priors["n_modpath_prior_edges"], errors="coerce"), width, color=ORANGE, label="MODPATH priors")
    ax.bar(x + width, pd.to_numeric(priors["n_output_edges"], errors="coerce"), width, color=BLUE, label="Output graph")
    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    ax.set_ylabel("Directed edges")
    ax.set_title("Prior-mode edge counts")
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        0.98,
        0.96,
        "Not independent\nvalidation",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color=VERMILLION,
        fontsize=6,
    )
    _save(fig, "Manuscript_EDFig1_M4_MODPATH_Prior_Mode_Audit")


def main() -> None:
    make_figure1_workflow()
    make_figure2_performance()
    make_figure3_edge_networks()
    make_figure4_endpoint_validation()
    make_figure5_external_archive_validation()
    make_figure6_capture_travel_time_validation()
    make_edfigure1_prior_audit()
    print(f"Wrote M4 manuscript-ready figures to {FIGURE_DIR}")


if __name__ == "__main__":
    sys.exit(main())
