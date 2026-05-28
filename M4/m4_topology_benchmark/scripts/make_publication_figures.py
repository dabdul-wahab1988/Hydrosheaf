"""Create Q1 manuscript figures matching figure_Table_guide.txt specification.

Produces 4 main-text figures + 7 supplementary figures from M4 benchmark results.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

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
PUBLIC_DIR = RESULT_DIR / "public_archives" / "savage"
FIGURE_DIR = BENCHMARK_ROOT / "figures" / "Manuscript_Ready"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Nature-style palette
BLUE = "#0072B2"
SKY = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#E69F00"
VERMILLION = "#D55E00"
PURPLE = "#CC79A7"
GRAY = "#666666"
LIGHT_GRAY = "#D9D9D9"
DARK = "#333333"
FN_DARK = "#8C510A"

SCENARIO_LABELS = {
    "spatial_only": "Spatial only",
    "head_gradient": "Head gradient",
    "head_depth": "Head depth",
    "hydrostratigraphic": "Hydrostrat.",
    "sparse_node": "Sparse node",
    "negative_random": "Random",
    "negative_wrong_direction": "Wrong dir.",
    "negative_shortcut": "Shortcut",
    "head_gradient_bayesian_hodge": "Hodge pruned",
    "real_head_projected_gradient": "Proj. gradient",
}

EVIDENCE_LEVELS = {
    # spatial_only uses no hydraulic data and achieves F1 = 0.0 → geometry-only control (level 0)
    "spatial_only": 0,
    "head_gradient": 3,
    "head_depth": 4,
    "hydrostratigraphic": 5,
    # sparse_node is a sensitivity analysis (node-subsampling sweep), NOT a scenario — excluded from here
    "negative_random": 0,
    "negative_wrong_direction": 0,
    "negative_shortcut": 0,
    "head_gradient_bayesian_hodge": 3,
    "real_head_projected_gradient": 3,
}

EVIDENCE_LABELS = {
    0: "Control\n(no skill)",
    3: "Head\ngradient",
    4: "Head\ndepth",
    5: "Hydrostrat.",
    6: "Prior-assisted\n(not independent)",
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.axisbelow": True,  # Keep grid lines behind plots
        "figure.dpi": 150,
    }
)


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _first_existing_column(df: pd.DataFrame, *columns: str) -> pd.Series:
    for column in columns:
        if column in df.columns:
            return df[column]
    return pd.Series(np.nan, index=df.index)


def _save(fig: plt.Figure, name: str) -> None:
    for suffix in (".png", ".svg"):
        out = FIGURE_DIR / f"{name}{suffix}"
        tmp = out.with_name(f"{out.stem}.tmp{out.suffix}")
        fig.savefig(tmp, dpi=600, bbox_inches="tight")
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
    # Position panel labels consistently at the top-left of the subplot cell to avoid overlaps
    try:
        fig = ax.figure
        ss = ax.get_subplotspec()
        if ss is not None:
            bbox = ss.get_position(fig)
            fig.text(
                bbox.x0 - 0.025, bbox.y1 + 0.018, label,
                transform=fig.transFigure,
                fontsize=10.0, fontweight="bold", va="bottom", ha="left",
                clip_on=False
            )
            return
    except Exception:
        pass

    # Fallback to ax.text if subplotspec is not available
    ax.text(
        -0.15, 1.04, label, transform=ax.transAxes,
        fontsize=10.0, fontweight="bold", va="bottom", ha="left",
        clip_on=False
    )



def _add_gridlines(ax: plt.Axes, x: bool = False, y: bool = True, which: str = "major") -> None:
    # Add subtle dotted gridlines — subtler for Q1 journal readability
    ax.grid(True, which=which, linestyle="--", alpha=0.4, linewidth=0.5, color="#E0E0E0")


def _display_scenario(value: str) -> str:
    return SCENARIO_LABELS.get(value, value.replace("_", " "))


# ---------------------------------------------------------------------------
# MAIN FIGURE 1: From MODPATH particles to falsifiable groundwater graphs
# ---------------------------------------------------------------------------

def make_figure1_workflow() -> None:
    """4-panel validation workflow: particles, reference graph, inference, metrics."""
    nodes_df = _read_csv(PUBLIC_DIR / "modpath_node_mapping.csv")
    ref_edges = _read_csv(PUBLIC_DIR / "modpath_reference_edges.csv")
    pathline = _read_csv(RESULT_DIR / "external_modpath_pathline_structure.csv")
    pathline_points = _read_pathline_sample(PUBLIC_DIR / "modpath_pathlines_standardised.csv")
    edge_class = _read_csv(RESULT_DIR / "edge_classification.csv")
    perf = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")

    fig = plt.figure(figsize=(4.5, 7.5))
    # fig.suptitle removed per Q1 guidelines
    fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.10, hspace=0.35)
    gs = fig.add_gridspec(2, 1)

    # --- Panel a: MODPATH particles ---
    ax = fig.add_subplot(gs[0])
    if not pathline_points.empty:
        for _, group in pathline_points.groupby("particle_id", sort=False):
            group = group.sort_values("sequence")
            ax.plot(
                pd.to_numeric(group["x"], errors="coerce"),
                pd.to_numeric(group["y"], errors="coerce"),
                color=GRAY,
                alpha=0.08,
                linewidth=0.20,
                zorder=1,
            )
    if not pathline.empty:
        n_show = min(800, len(pathline))
        sample = pathline.sample(n=n_show, random_state=42) if len(pathline) > n_show else pathline
        ax.scatter(
            pd.to_numeric(sample["endpoint_start_x"], errors="coerce"),
            pd.to_numeric(sample["endpoint_start_y"], errors="coerce"),
            s=5.0, color=BLUE, alpha=0.6, linewidth=0, label="Start points", zorder=2,
        )
        ax.scatter(
            pd.to_numeric(sample["endpoint_end_x"], errors="coerce"),
            pd.to_numeric(sample["endpoint_end_y"], errors="coerce"),
            s=5.0, color=VERMILLION, alpha=0.6, linewidth=0, label="End points", zorder=2,
        )
        n_parts = pathline["particle_id"].nunique()
        ax.set_title(f"MODPATH particles (n={n_parts})")
        ax.legend(
            frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20),
            ncol=2, markerscale=2, fontsize=6.5, handletextpad=0.4,
            columnspacing=1.4,
        )
    ax.set_aspect("auto")
    ax.set_xlabel("x (model coords)")
    ax.set_ylabel("y (model coords)")
    _add_gridlines(ax, x=True, y=True)
    _panel_label(ax, "a")

    # --- Panel b: MODPATH reference transition grid ---
    ax = fig.add_subplot(gs[1])
    if not nodes_df.empty and not ref_edges.empty:
        pos = _build_node_positions(nodes_df)
        cell_support = _transition_cell_support(ref_edges, weight_col="particle_count")
        _draw_grid_cells(ax, pos, cell_support, cmap_name="YlGn", alpha_multiplier=1.08)
        pruned = _draw_weighted_reference_edges(ax, ref_edges, pos, min_particle_count=15)
        _draw_cell_centers(ax, pos, s=5.0, color=GRAY, alpha=0.55, zorder=5)
        n_edges = len(ref_edges)
        ax.set_title(f"Reference transition grid ({len(pruned)}/{n_edges} edges)")
        handles = [
            Line2D([0], [0], color=GREEN, lw=0.85, alpha=0.65, label="Edge support\nby particle count"),
        ]
        ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=5.8,
                  handlelength=1.6, borderpad=0.2)
        _set_graph_extent(ax, pos)
    ax.set_aspect("auto")
    ax.axis("off")
    _panel_label(ax, "b")

    _save(fig, "Fig1_From_MODPATH_particles_to_falsifiable_groundwater_graphs")


def _build_node_positions(nodes_df: pd.DataFrame) -> Dict[str, Tuple[float, float]]:
    pos = {}
    for _, row in nodes_df.iterrows():
        node_id = str(row["node_id"])
        x = pd.to_numeric(row["x"], errors="coerce")
        y = pd.to_numeric(row["y"], errors="coerce")
        if pd.notna(x) and pd.notna(y):
            pos[node_id] = (float(x), float(y))
    return pos


def _read_pathline_sample(
    path: Path,
    max_particles: int = 160,
    max_points_per_particle: int = 80,
    random_state: int = 42,
) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    cols = ["particle_id", "x", "y", "sequence"]
    try:
        df = pd.read_csv(path, usecols=cols)
    except (ValueError, pd.errors.EmptyDataError):
        return pd.DataFrame()
    if df.empty:
        return df

    particle_ids = pd.Series(df["particle_id"].dropna().unique())
    if len(particle_ids) > max_particles:
        particle_ids = particle_ids.sample(n=max_particles, random_state=random_state)
    sampled = df[df["particle_id"].isin(set(particle_ids))].copy()

    rows = []
    for _, group in sampled.groupby("particle_id", sort=False):
        group = group.sort_values("sequence")
        if len(group) > max_points_per_particle:
            idx = np.linspace(0, len(group) - 1, max_points_per_particle).round().astype(int)
            group = group.iloc[np.unique(idx)]
        rows.append(group)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=cols)


def _infer_cell_size(pos: Dict[str, Tuple[float, float]]) -> Tuple[float, float]:
    if not pos:
        return 1.0, 1.0
    xs = sorted({xy[0] for xy in pos.values()})
    ys = sorted({xy[1] for xy in pos.values()})
    xdiffs = [b - a for a, b in zip(xs, xs[1:]) if b > a]
    ydiffs = [b - a for a, b in zip(ys, ys[1:]) if b > a]
    dx = float(np.median(xdiffs)) if xdiffs else 1.0
    dy = float(np.median(ydiffs)) if ydiffs else dx
    return dx * 0.72, dy * 0.72


def _transition_cell_support(edges_df: pd.DataFrame, weight_col: str = "particle_count") -> Dict[str, float]:
    support: Dict[str, float] = {}
    if edges_df.empty:
        return support
    for _, row in edges_df.iterrows():
        weight = row.get(weight_col, 1.0)
        weight = pd.to_numeric(weight, errors="coerce")
        if pd.isna(weight):
            weight = 1.0
        for col in ("u", "v"):
            node = str(row[col])
            support[node] = support.get(node, 0.0) + float(weight)
    return support


def _draw_grid_cells(
    ax,
    pos: Dict[str, Tuple[float, float]],
    support: Dict[str, float],
    cmap_name: str = "Greens",
    vmax: float | None = None,
    alpha_multiplier: float = 1.0,
    edge_alpha: float = 1.0,
    zorder: int = 0,
) -> None:
    if not pos:
        return
    dx, dy = _infer_cell_size(pos)
    values = [float(v) for v in support.values() if pd.notna(v)]
    vmax = float(vmax if vmax is not None else (max(values) if values else 1.0))
    vmax = max(vmax, 1.0)
    cmap = plt.get_cmap(cmap_name)

    for node, (x, y) in pos.items():
        value = float(support.get(node, 0.0))
        if value > 0:
            scaled = min(value / vmax, 1.0)
            face = cmap(0.14 + 0.78 * scaled)
            alpha = (0.30 + 0.40 * scaled) * alpha_multiplier
        else:
            face = "white"
            alpha = 0.20 * alpha_multiplier
        rect = Rectangle(
            (x - dx / 2, y - dy / 2),
            dx,
            dy,
            facecolor=face,
            edgecolor=LIGHT_GRAY,
            linewidth=0.22,
            alpha=min(alpha, 1.0),
            zorder=zorder,
        )
        ax.add_patch(rect)


def _draw_cell_centers(
    ax,
    pos: Dict[str, Tuple[float, float]],
    s: float = 5.0,
    color: str = GRAY,
    alpha: float = 0.55,
    zorder: int = 5,
) -> None:
    if not pos:
        return
    xs, ys = zip(*pos.values())
    ax.scatter(xs, ys, s=s, facecolor="white", edgecolor=color, linewidth=0.35,
               alpha=alpha, zorder=zorder)


def _attach_reference_support(edge_class_df: pd.DataFrame, ref_edges: pd.DataFrame) -> pd.DataFrame:
    df = edge_class_df.copy()
    if "edge" in df.columns:
        df = df.rename(columns={"edge": "edge_id"})
    if "particle_count" in ref_edges.columns:
        support = ref_edges[["edge_id", "particle_count"]].copy()
        support["particle_count"] = pd.to_numeric(support["particle_count"], errors="coerce")
        df = df.merge(support, on="edge_id", how="left")
    else:
        df["particle_count"] = np.nan
    df["_diagnostic_weight"] = 1.0
    ref_mask = df["classification"].isin(["TP", "FN"]) & df["particle_count"].notna()
    df.loc[ref_mask, "_diagnostic_weight"] = pd.to_numeric(
        df.loc[ref_mask, "particle_count"], errors="coerce"
    ).fillna(1.0)
    return df


def _diagnostic_cmap(label: str) -> str:
    return {"TP": "Greens", "FP": "Oranges", "FN": "PuOr"}.get(label, "Greys")


def _draw_graph_edges(ax, edges_df, u_col, v_col, pos, color, alpha=0.6, lw=0.4):
    for _, row in edges_df.iterrows():
        u = str(row[u_col])
        v = str(row[v_col])
        if u in pos and v in pos:
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle="-|>", color=color,
                                        lw=lw, alpha=alpha, mutation_scale=6))


def _draw_weighted_reference_edges(
    ax,
    edges_df: pd.DataFrame,
    pos: Dict[str, Tuple[float, float]],
    min_particle_count: int = 15,
) -> pd.DataFrame:
    df = edges_df.copy()
    if "particle_count" in df.columns:
        df["_particle_support"] = pd.to_numeric(df["particle_count"], errors="coerce").fillna(0)
    else:
        df["_particle_support"] = 1.0

    pruned = df[df["_particle_support"] >= min_particle_count].copy()
    if pruned.empty:
        pruned = df.nlargest(min(80, len(df)), "_particle_support").copy()

    s_min = float(pruned["_particle_support"].min())
    s_max = float(pruned["_particle_support"].max())
    denom = max(s_max - s_min, 1.0)

    for _, row in pruned.sort_values("_particle_support").iterrows():
        u = str(row["u"])
        v = str(row["v"])
        if u not in pos or v not in pos:
            continue
        support = float(row["_particle_support"])
        scaled = (support - s_min) / denom
        lw = 0.18 + 0.67 * scaled
        alpha = 0.18 + 0.35 * scaled
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="-|>",
                color=GREEN,
                lw=lw,
                alpha=alpha,
                mutation_scale=5.5 + 1.5 * scaled,
                shrinkA=1.5,
                shrinkB=1.5,
            ),
        )
    return pruned


def _draw_edge_subset(
    ax,
    edges_df: pd.DataFrame,
    pos,
    color,
    alpha=0.55,
    lw=0.45,
    linestyle="solid",
    mutation_scale=4.5,
    weight_col: str | None = None,
    zorder: int = 4,
):
    if weight_col and weight_col in edges_df.columns and not edges_df.empty:
        weights = pd.to_numeric(edges_df[weight_col], errors="coerce").fillna(1.0)
        w_min = float(weights.min())
        w_max = float(weights.max())
        denom = max(w_max - w_min, 1.0)
    else:
        weights = None
        w_min = 0.0
        denom = 1.0

    for _, row in edges_df.iterrows():
        u = str(row["u"])
        v = str(row["v"])
        if u in pos and v in pos:
            if weight_col and weight_col in row.index:
                weight = pd.to_numeric(row[weight_col], errors="coerce")
                scaled = 0.0 if pd.isna(weight) else (float(weight) - w_min) / denom
            else:
                scaled = 0.0
            edge_lw = lw * (0.75 + 0.65 * scaled)
            edge_alpha = min(alpha * (0.78 + 0.30 * scaled), 1.0)
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            ax.annotate(
                "",
                xy=(x2, y2),
                xytext=(x1, y1),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=color,
                    lw=edge_lw,
                    alpha=edge_alpha,
                    linestyle=linestyle,
                    mutation_scale=mutation_scale * (0.92 + 0.22 * scaled),
                    shrinkA=1.3,
                    shrinkB=1.3,
                ),
                zorder=zorder,
            )


def _set_graph_extent(ax, pos: Dict[str, Tuple[float, float]], pad_fraction: float = 0.06) -> None:
    if not pos:
        return
    xs, ys = zip(*pos.values())
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xpad = max((xmax - xmin) * pad_fraction, 1.0)
    ypad = max((ymax - ymin) * pad_fraction, 1.0)
    ax.set_xlim(xmin - xpad, xmax + xpad)
    ax.set_ylim(ymin - ypad, ymax + ypad)


def _draw_graph_nodes(ax, pos, s=6, color=GRAY):
    if not pos:
        return
    xs, ys = zip(*pos.values())
    ax.scatter(xs, ys, s=s, facecolor="white", edgecolor=color, linewidth=0.4, zorder=3)


def _draw_classified_edges(ax, edge_class_df, pos):
    styles = {
        "TP": (GREEN, 0.55, 0.5),
        "FP": (VERMILLION, 0.40, 0.5),
        "FN": (ORANGE, 0.45, 0.5),
    }
    for _, row in edge_class_df.iterrows():
        label = str(row["classification"])
        if label == "TN" or label not in styles:
            continue
        u = str(row["u"])
        v = str(row["v"])
        if u in pos and v in pos:
            color, alpha, lw = styles[label]
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle="-|>", color=color,
                                        lw=lw, alpha=alpha, mutation_scale=5))


def _draw_evidence_ladder(ax, perf, priors):
    """Draw evidence ladder as a structured horizontal bar diagram."""
    ax.set_xlim(0, 8)
    ax.set_ylim(-0.5, 11.5)
    ax.axis("off")

    # sparse_node excluded — it is a sensitivity analysis, not an inference scenario.
    # spatial_only is now a level-0 geometry-only control.
    independent_scenarios = [
        ("negative_random", "Negative control (random edges)"),
        ("negative_wrong_direction", "Negative control (wrong direction)"),
        ("negative_shortcut", "Negative control (shortcut edges)"),
        ("spatial_only", "Spatial proximity only (geometry control)"),
        ("head_gradient", "Head-gradient constrained"),
        ("head_gradient_bayesian_hodge", "Head-gradient (Hodge pruned)"),
        ("real_head_projected_gradient", "Proj. gradient constrained"),
        ("head_depth", "Head-depth constrained"),
        ("hydrostratigraphic", "Hydrostratigraphic constrained"),
    ]
    if not perf.empty:
        independent_scenarios = [(s, label) for s, label in independent_scenarios
                                 if s in perf["scenario"].values]

    level_colors = {
        0: LIGHT_GRAY,
        1: SKY, 2: SKY,
        3: BLUE, 4: BLUE, 5: BLUE,
        6: ORANGE, 7: PURPLE,
    }

    y = 0
    for scenario_key, label in independent_scenarios:
        level = EVIDENCE_LEVELS.get(scenario_key, 0)
        color = level_colors.get(level, LIGHT_GRAY)
        rect = Rectangle((0.3, y - 0.35), 4.2, 0.65, facecolor=color, edgecolor="white",
                          linewidth=0.5, alpha=0.75)
        ax.add_patch(rect)
        ax.text(0.5, y, label, va="center", fontsize=6.3, color=DARK)

        if not perf.empty:
            row = perf[perf["scenario"] == scenario_key]
            if not row.empty:
                f1 = float(row["f1"].iloc[0]) if pd.notna(row["f1"].iloc[0]) else 0
                ax.text(4.8, y, f"F1={f1:.3f}", va="center", fontsize=6.0, color=DARK,
                        fontweight="bold" if f1 > 0.5 else "normal")

        ax.text(5.8, y, f"Level {level}", va="center", fontsize=6.0, color=GRAY)
        y += 0.85

    # Divider
    ax.axhline(y + 0.1, xmin=0.35, xmax=0.95, color=GRAY, linewidth=0.6, linestyle=(0, (3, 2)))
    y += 0.5

    # Prior-assisted section
    ax.text(0.5, y, "Prior-assisted (not independent)", va="center", fontsize=7.0,
            fontweight="bold", color=VERMILLION)
    y += 0.7

    prior_modes = [
        ("modpath_prior_override", "Override mode"),
        ("modpath_prior_merge", "Merge mode"),
        ("modpath_prior_only", "Only-prior mode"),
    ]
    for scenario_key, label in prior_modes:
        rect = Rectangle((0.3, y - 0.3), 4.2, 0.55, facecolor=ORANGE, edgecolor="white",
                          linewidth=0.5, alpha=0.5)
        ax.add_patch(rect)
        ax.text(0.5, y, label, va="center", fontsize=6.3, color=DARK)
        ax.text(5.8, y, "Level 6", va="center", fontsize=6.0, color=GRAY)
        y += 0.65

    # Public-archive projection diagnostic
    rect = Rectangle((0.3, y - 0.3), 4.2, 0.55, facecolor=PURPLE, edgecolor="white",
                      linewidth=0.5, alpha=0.5)
    ax.add_patch(rect)
    ax.text(0.5, y, "Public archive projection (Savage)", va="center", fontsize=6.3, color="white")
    ax.text(5.8, y, "Diagnostic", va="center", fontsize=6.0, color=GRAY)

    ax.set_title("Evidence ladder")


# ---------------------------------------------------------------------------
# MAIN FIGURE 2: Independent topology performance across evidence levels
# ---------------------------------------------------------------------------

def make_figure2_performance() -> None:
    """4-panel 2×2: (a) evidence level & F1, (b) P/R/F1 bars,
    (c) node-sparsity sensitivity P/R/F1 merged, (d) scale-mismatch check."""
    perf = _read_csv(RESULT_DIR / "independent_graph_vs_modpath.csv")
    sparsity = _read_csv(PUBLIC_DIR / "node_sparsity_sensitivity.csv")
    if perf.empty:
        return

    ind = perf.copy()
    # sparse_node is a sensitivity analysis (node-subsampling sweep), not a scenario —
    # it lives exclusively in panels c and d. spatial_only is a geometry-only control (level 0).
    ordered = ["spatial_only", "head_gradient", "head_gradient_bayesian_hodge",
               "real_head_projected_gradient", "head_depth", "hydrostratigraphic",
               "negative_random", "negative_wrong_direction", "negative_shortcut"]
    ordered = [s for s in ordered if s in ind["scenario"].values]
    ind["_order"] = ind["scenario"].map({k: i for i, k in enumerate(ordered)})
    ind = ind[ind["scenario"].isin(ordered)].sort_values("_order")

    fig = plt.figure(figsize=(7.5, 6.2))
    # fig.suptitle removed per Q1 guidelines
    fig.subplots_adjust(left=0.10, right=0.96, top=0.92, bottom=0.20, hspace=0.52, wspace=0.38)
    gs = fig.add_gridspec(2, 2)

    # --- Panel a: Evidence level and F1 ---
    ax = fig.add_subplot(gs[0, 0])
    scenarios = ind["scenario"].tolist()
    levels = [EVIDENCE_LEVELS.get(s, 0) for s in scenarios]
    f1_vals = pd.to_numeric(ind["f1"], errors="coerce").fillna(0).values
    colors = [_evidence_color(lv) for lv in levels]
    y_pos = np.arange(len(scenarios))
    bars = ax.barh(y_pos, f1_vals, color=colors, height=0.55, edgecolor="white", linewidth=0.4)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([_display_scenario(s) for s in scenarios])
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("F1 score")
    ax.invert_yaxis()
    
    # Legend: spatial_only is now L0 (geometry-only control); no L1-2 tier remains
    legend_elements = [
        Line2D([0], [0], color=LIGHT_GRAY, lw=3, label="L0: Controls (incl. spatial-only)"),
        Line2D([0], [0], color=BLUE, lw=3, label="L3-5: Hydraulic evidence"),
    ]
    ax.legend(handles=legend_elements, frameon=False, loc="lower right", fontsize=5.5)
    ax.set_title("Evidence level and F1")
    ax.axvline(0.5, color=LIGHT_GRAY, linewidth=0.6, linestyle=(0, (3, 2)))
    _add_gridlines(ax, x=True, y=False)
    _panel_label(ax, "a")

    # --- Panel b: Precision, Recall, F1 grouped bars (no counts) ---
    ax = fig.add_subplot(gs[0, 1])
    x = np.arange(len(ind))
    w = 0.22
    prec = pd.to_numeric(ind["precision"], errors="coerce").fillna(0).values
    rec = pd.to_numeric(ind["recall"], errors="coerce").fillna(0).values
    f1 = f1_vals

    ax.bar(x - w, prec, w, color=SKY, label="Precision", edgecolor="white", linewidth=0.4)
    ax.bar(x, rec, w, color=BLUE, label="Recall", edgecolor="white", linewidth=0.4)
    ax.bar(x + w, f1, w, color=GREEN, label="F1", edgecolor="white", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([_display_scenario(s) for s in scenarios], rotation=35, ha="right")
    ax.set_ylim(0, 1.0)  # Bound strictly to 1.0
    ax.set_ylabel("Score")
    ax.set_title("Precision / Recall / F1 (grouped bars)")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "b")

    # --- Panel c: Node-sparsity sensitivity — Precision, Recall, F1 merged ---
    ax = fig.add_subplot(gs[1, 0])
    if not sparsity.empty:
        sparse_c = sparsity.sort_values("node_fraction").copy()
        xc = pd.to_numeric(sparse_c["node_fraction"], errors="coerce")
        for mc, sc, col, lbl in [
            ("mean_precision", "std_precision", SKY,   "Precision"),
            ("mean_recall",    "std_recall",    BLUE,  "Recall"),
            ("mean_f1",        "std_f1",        GREEN, "F1"),
        ]:
            yc   = pd.to_numeric(sparse_c[mc], errors="coerce")
            yerrc = 1.96 * pd.to_numeric(sparse_c[sc], errors="coerce").fillna(0)
            ax.errorbar(xc, yc, yerr=yerrc, marker="o", markersize=3.8,
                        linewidth=1.1, capsize=2.0, color=col, label=lbl)
        ax.set_xlim(0.05, 1.05)
        ax.set_ylim(0, 1.0)
        ax.set_xlabel("Node fraction retained")
        ax.set_ylabel("Score")
        ax.set_title("Node-sparsity sensitivity")
        ax.text(0.95, 0.05, "mean ± 95% CI", transform=ax.transAxes, fontsize=5.5,
                color=GRAY, ha="right", va="bottom")
        if "successful_trials" in sparse_c.columns and "planned_trials" in sparse_c.columns:
            n_ok_c  = pd.to_numeric(sparse_c["successful_trials"], errors="coerce")
            n_tot_c = pd.to_numeric(sparse_c["planned_trials"],    errors="coerce")
            ylo, yhi = ax.get_ylim()
            for xi, ok, tot in zip(xc, n_ok_c, n_tot_c):
                if pd.notna(ok) and pd.notna(tot):
                    ax.annotate(f"{int(ok)}/{int(tot)}", xy=(xi, ylo),
                                xytext=(xi, ylo - 0.07*(yhi-ylo)),
                                fontsize=5.5, color=GRAY, ha="center", va="top",
                                rotation=45, annotation_clip=False)
            ax.text(0.02, -0.13, "n = successful/planned trials",
                    transform=ax.transAxes, fontsize=5.5, color=GRAY, ha="left", va="top")
        _add_gridlines(ax, x=True, y=True)
    else:
        ax.text(0.5, 0.5, "No sparsity data", ha="center", va="center", transform=ax.transAxes)
    _panel_label(ax, "c")

    # --- Panel d: Scale-mismatch diagnostic ---
    ax = fig.add_subplot(gs[1, 1])
    ref_len = pd.to_numeric(ind["median_reference_length"], errors="coerce").fillna(0)
    inf_len = pd.to_numeric(ind["median_inferred_length"], errors="coerce").fillna(0)
    ax.bar(x - w / 2, ref_len, w, color=GRAY, label="Reference", edgecolor="white", linewidth=0.4)
    ax.bar(x + w / 2, inf_len, w, color=BLUE, label="Inferred", edgecolor="white", linewidth=0.4)
    scale = ind["scale_mismatch"].astype(str).str.lower().eq("true")
    for i, is_mm in enumerate(scale):
        if is_mm:
            ymax = max(float(ref_len.iloc[i]), float(inf_len.iloc[i]))
            ax.plot(i, ymax * 1.08, marker="v", color=VERMILLION, markersize=5, clip_on=False)
    ax.set_xticks(x)
    ax.set_xticklabels([_display_scenario(s) for s in scenarios], rotation=35, ha="right")
    ax.set_ylabel("Euclidean edge distance in model units", fontsize=7.5)
    ax.set_title("Scale-mismatch check (▼ = mismatch)")
    ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "d")

    # Two separate shared legends at the bottom, visually separated
    perf_handles = [
        Line2D([0], [0], color=SKY, lw=2.5, label="Precision"),
        Line2D([0], [0], color=BLUE, lw=2.5, label="Recall"),
        Line2D([0], [0], color=GREEN, lw=2.5, label="F1"),
    ]
    scale_handles = [
        Line2D([0], [0], color=GRAY, lw=2.5, label="Reference"),
        Line2D([0], [0], color=BLUE, lw=2.5, label="Inferred"),
    ]

    fig.text(0.10, 0.06, "Performance metrics:", fontsize=7.5, fontweight="bold", va="center")
    fig.legend(handles=perf_handles, frameon=False, ncol=3, loc="center left",
               bbox_to_anchor=(0.25, 0.06), fontsize=7)

    fig.text(0.60, 0.06, "Scale diagnostic:", fontsize=7.5, fontweight="bold", va="center")
    fig.legend(handles=scale_handles, frameon=False, ncol=2, loc="center left",
               bbox_to_anchor=(0.75, 0.06), fontsize=7)

    _save(fig, "Fig2_Independent_topology_performance_across_evidence_levels")


def _evidence_color(level: int) -> str:
    return {0: LIGHT_GRAY, 1: SKY, 2: SKY, 3: BLUE, 4: BLUE, 5: BLUE}.get(level, GRAY)


# ---------------------------------------------------------------------------
# MAIN FIGURE 3: Edge-level falsification networks
# ---------------------------------------------------------------------------

def _get_classified_edges(scenario: str, nodes_df: pd.DataFrame, ref_edges: pd.DataFrame) -> pd.DataFrame:
    # Build samples
    samples = []
    for _, row in nodes_df.iterrows():
        samples.append({
            "site_id": str(row["node_id"]),
            "lat": float(row["y"]),
            "lon": float(row["x"]),
            "elevation": float(row["z"]) if not pd.isna(row["z"]) else 0.0,
        })

    ref_edge_list = [(str(row["u"]), str(row["v"])) for _, row in ref_edges.iterrows()]
    ref_set = set(ref_edge_list)
    pos = _build_node_positions(nodes_df)

    from hydrosheaf.graph.build import infer_edges_from_coordinates

    if scenario == "spatial_only":
        obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=True)
        inf_edges = [(e.u, e.v) for e in obj]
    elif scenario == "head_gradient":
        obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
        inf_edges = [(e.u, e.v) for e in obj]
    elif scenario == "head_gradient_bayesian_hodge":
        obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
        inf_edges = [(e.u, e.v) for e in obj]
    elif scenario == "real_head_projected_gradient":
        obj = infer_edges_from_coordinates(samples, max_neighbors=2, allow_uphill=False)
        inf_edges = [(e.u, e.v) for e in obj]
    elif scenario == "head_depth":
        depth_samples = []
        for s in samples:
            s_copy = dict(s)
            elev = float(s_copy.get("elevation", 0))
            s_copy["depth_tier"] = "shallow" if elev > -10 else ("mid" if elev > -50 else "deep")
            depth_samples.append(s_copy)
        obj = infer_edges_from_coordinates(depth_samples, max_neighbors=3, allow_uphill=False)
        inf_edges = [(e.u, e.v) for e in obj]
    elif scenario == "hydrostratigraphic":
        strat_samples = []
        for s in samples:
            s_copy = dict(s)
            elev = float(s_copy.get("elevation", 0))
            s_copy["aquifer_unit"] = "shallow" if elev > -20 else "deep"
            strat_samples.append(s_copy)
        candidates_raw = infer_edges_from_coordinates(strat_samples, max_neighbors=5, allow_uphill=False)
        strat_obj = []
        for e in candidates_raw:
            u_sample = next((x for x in strat_samples if x["site_id"] == e.u), None)
            v_sample = next((x for x in strat_samples if x["site_id"] == e.v), None)
            if u_sample and v_sample:
                if u_sample.get("aquifer_unit") == v_sample.get("aquifer_unit"):
                    strat_obj.append(e)
        inf_edges = [(e.u, e.v) for e in strat_obj]
    elif scenario == "sparse_node":
        rng = np.random.default_rng(20260521)
        n_half = len(samples) // 2
        idxs = rng.choice(len(samples), size=n_half, replace=False)
        sub_samples = [samples[i] for i in idxs]
        sub_ids = {s["site_id"] for s in sub_samples}
        obj = infer_edges_from_coordinates(sub_samples, max_neighbors=2, allow_uphill=False)
        inf_edges = [(e.u, e.v) for e in obj]
        ref_set = {(u, v) for u, v in ref_edge_list if u in sub_ids and v in sub_ids}
    elif scenario == "negative_random":
        n_random = len(ref_edge_list)
        random_pairs = []
        node_list = list(pos.keys())
        rng_rand = np.random.default_rng(20260521)
        for _ in range(n_random):
            u = str(rng_rand.choice(node_list))
            v = str(rng_rand.choice(node_list))
            while u == v or (u, v) in random_pairs:
                v = str(rng_rand.choice(node_list))
            random_pairs.append((u, v))
        inf_edges = random_pairs
    elif scenario == "negative_wrong_direction":
        inf_edges = [(v, u) for (u, v) in ref_edge_list]
    elif scenario == "negative_shortcut":
        ref_graph = {}
        for u, v in ref_edge_list:
            ref_graph.setdefault(u, []).append(v)
        shortcuts = []
        for u, v in ref_edge_list:
            for w in ref_graph.get(v, []):
                if w != u and (u, w) not in ref_edge_list:
                    shortcuts.append((u, w))
        inf_edges = list(set(shortcuts))[:len(ref_edge_list)]
    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    inf_set = set(inf_edges)
    rows = []
    # True positives
    for u, v in sorted(ref_set & inf_set):
        rows.append({"edge_id": f"{u}->{v}", "u": u, "v": v, "classification": "TP"})
    # False negatives
    for u, v in sorted(ref_set - inf_set):
        rows.append({"edge_id": f"{u}->{v}", "u": u, "v": v, "classification": "FN"})
    # False positives
    for u, v in sorted(inf_set - ref_set):
        rows.append({"edge_id": f"{u}->{v}", "u": u, "v": v, "classification": "FP"})

    return pd.DataFrame(rows)


def _classified_metrics(classified: pd.DataFrame) -> Dict[str, float]:
    if classified.empty or "classification" not in classified.columns:
        return {"tp": 0, "fp": 0, "fn": 0, "precision": 0.0, "recall": 0.0, "f1": 0.0}
    tp = int((classified["classification"] == "TP").sum())
    fp = int((classified["classification"] == "FP").sum())
    fn = int((classified["classification"] == "FN").sum())
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1}


def _format_metric_title(label: str, metrics: Dict[str, float]) -> str:
    return (
        f"{label}\n"
        f"TP={int(metrics['tp'])}, FP={int(metrics['fp'])}, FN={int(metrics['fn'])}; "
        f"P={metrics['precision']:.2f}, R={metrics['recall']:.2f}, F1={metrics['f1']:.2f}"
    )


def _sample_edges(edges: pd.DataFrame, max_edges: int, random_state: int) -> pd.DataFrame:
    if len(edges) <= max_edges:
        return edges
    return edges.sample(n=max_edges, random_state=random_state)


def _draw_reference_context(
    ax,
    ref_edges: pd.DataFrame,
    pos: Dict[str, Tuple[float, float]],
    alpha: float = 0.12,
) -> None:
    if ref_edges.empty:
        return
    for _, row in ref_edges.iterrows():
        u = str(row["u"])
        v = str(row["v"])
        if u in pos and v in pos:
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            ax.plot([x1, x2], [y1, y2], color=LIGHT_GRAY, linewidth=0.22,
                    alpha=alpha, zorder=1)


def make_figure3_edge_networks() -> None:
    """Diagnostic controls for edge-level graph-inference failure modes."""
    nodes_df = _read_csv(PUBLIC_DIR / "modpath_node_mapping.csv")
    ref_edges = _read_csv(PUBLIC_DIR / "modpath_reference_edges.csv")

    if nodes_df.empty:
        return

    pos = _build_node_positions(nodes_df)

    fig = plt.figure(figsize=(7.2, 6.4))
    # fig.suptitle removed per Q1 guidelines
    fig.subplots_adjust(left=0.08, right=0.96, top=0.95, bottom=0.16, hspace=0.48, wspace=0.32)
    gs = fig.add_gridspec(2, 2)
    
    styles = [
        ("TP", GREEN, 0.58, 0.42, "solid", 4.5),
        ("FP", VERMILLION, 0.58, 0.42, "solid", 4.5),
        ("FN", FN_DARK, 0.86, 0.80, (0, (3, 2)), 6.0),
    ]

    # Panel a: Head gradient
    ax = fig.add_subplot(gs[0, 0])
    head = _attach_reference_support(_get_classified_edges("head_gradient", nodes_df, ref_edges), ref_edges)
    _draw_graph_nodes(ax, pos, s=1, color=LIGHT_GRAY)
    for klass, color, alpha, lw, ls, mut in styles:
        subset = head[head["classification"] == klass]
        _draw_edge_subset(
            ax, subset, pos, color=color, alpha=alpha,
            lw=lw, linestyle=ls, mutation_scale=mut
        )
    ax.set_title("Head gradient", fontsize=8.5)
    active_nodes = set(head["u"].astype(str)) | set(head["v"].astype(str))
    active_pos = {n: pos[n] for n in active_nodes if n in pos}
    _set_graph_extent(ax, active_pos, pad_fraction=0.02)
    metrics = _classified_metrics(head)
    ax.text(0.03, 0.93, f"P = {metrics['precision']:.2f}, R = {metrics['recall']:.2f}, F1 = {metrics['f1']:.2f}",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=LIGHT_GRAY, linewidth=0.5))
    ax.set_aspect("auto")
    ax.axis("off")
    _panel_label(ax, "a")

    # Panel b: Real head-projected gradient
    ax = fig.add_subplot(gs[0, 1])
    proj = _attach_reference_support(_get_classified_edges("real_head_projected_gradient", nodes_df, ref_edges), ref_edges)
    _draw_graph_nodes(ax, pos, s=1, color=LIGHT_GRAY)
    for klass, color, alpha, lw, ls, mut in styles:
        subset = proj[proj["classification"] == klass]
        _draw_edge_subset(
            ax, subset, pos, color=color, alpha=alpha,
            lw=lw, linestyle=ls, mutation_scale=mut
        )
    ax.set_title("Real head-projected gradient", fontsize=8.5)
    active_nodes = set(proj["u"].astype(str)) | set(proj["v"].astype(str))
    active_pos = {n: pos[n] for n in active_nodes if n in pos}
    _set_graph_extent(ax, active_pos, pad_fraction=0.02)
    metrics = _classified_metrics(proj)
    ax.text(0.03, 0.93, f"P = {metrics['precision']:.2f}, R = {metrics['recall']:.2f}, F1 = {metrics['f1']:.2f}",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=LIGHT_GRAY, linewidth=0.5))
    ax.set_aspect("auto")
    ax.axis("off")
    _panel_label(ax, "b")

    # Panel c: Spatial only (geometry-only control — no hydraulic data, F1 = 0.0)
    ax = fig.add_subplot(gs[1, 0])
    spatial = _attach_reference_support(_get_classified_edges("spatial_only", nodes_df, ref_edges), ref_edges)
    _draw_graph_nodes(ax, pos, s=1, color=LIGHT_GRAY)
    for klass, color, alpha, lw, ls, mut in styles:
        subset = spatial[spatial["classification"] == klass]
        _draw_edge_subset(
            ax, subset, pos, color=color, alpha=alpha,
            lw=lw, linestyle=ls, mutation_scale=mut
        )
    ax.set_title("Spatial only (geometry control)", fontsize=8.5)
    active_nodes_sp = set(spatial["u"].astype(str)) | set(spatial["v"].astype(str))
    active_pos_sp = {n: pos[n] for n in active_nodes_sp if n in pos}
    _set_graph_extent(ax, active_pos_sp if active_pos_sp else pos, pad_fraction=0.02)
    metrics = _classified_metrics(spatial)
    ax.text(0.03, 0.93, f"P = {metrics['precision']:.2f}, R = {metrics['recall']:.2f}, F1 = {metrics['f1']:.2f}",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=LIGHT_GRAY, linewidth=0.5))
    ax.set_aspect("auto")
    ax.axis("off")
    _panel_label(ax, "c")

    # Panel d: Random edges
    ax = fig.add_subplot(gs[1, 1])
    random_control = _attach_reference_support(_get_classified_edges("negative_random", nodes_df, ref_edges), ref_edges)
    _draw_reference_context(ax, ref_edges, pos, alpha=0.10)
    random_fp = _sample_edges(random_control[random_control["classification"] == "FP"], max_edges=20, random_state=5)
    random_fn = _sample_edges(random_control[random_control["classification"] == "FN"], max_edges=20, random_state=6)
    
    _draw_grid_cells(ax, pos, _transition_cell_support(random_control, weight_col="_diagnostic_weight"),
                     cmap_name="Oranges", vmax=8)
    _draw_edge_subset(ax, random_fp, pos, color=VERMILLION, alpha=0.18, lw=0.35, mutation_scale=4.2)
    _draw_edge_subset(ax, random_fn, pos, color=FN_DARK, alpha=0.30, lw=0.50,
                      linestyle=(0, (3, 2)), mutation_scale=5.0)
    ax.set_title("Random edges", fontsize=8.5)
    _set_graph_extent(ax, pos)
    ax.set_aspect("auto")
    ax.axis("off")
    _panel_label(ax, "d")

    legend_elements = [
        Line2D([0], [0], color=GREEN, lw=1.8, label="TP"),
        Line2D([0], [0], color=VERMILLION, lw=1.8, label="FP"),
        Line2D([0], [0], color=FN_DARK, lw=1.8, linestyle=(0, (3, 2)), label="FN"),
        Line2D([0], [0], color=LIGHT_GRAY, lw=1.0, label="Reference"),
    ]
    fig.legend(handles=legend_elements, frameon=False, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, 0.045), fontsize=7)

    _save(fig, "Fig3_Edge_level_falsification_networks")

# -------------------------------------------------------------------------------------------
# MAIN FIGURE 4: External archive and guarded diagnostic extensions
# ---------------------------------------------------------------------------

def make_figure4_external_archive() -> None:
    """4-panel synthesis: (a) Savage archive projection, (b) Source-start point-cloud envelopes,
    (c) Prior-transfer consistency, (d) Empty."""
    summary = _read_csv(RESULT_DIR / "external_modpath_archive_summary.csv")
    pathline = _read_csv(RESULT_DIR / "external_modpath_pathline_structure.csv")
    time_df = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time.csv")
    time_summary = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv")

    if summary.empty:
        return

    # 2x2 layout with 4th panel empty
    fig = plt.figure(figsize=(7.2, 6.4))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.93, bottom=0.12, hspace=0.45, wspace=0.42)

    top = summary[summary["validation_tier"] == "tier_1_savage"]
    if top.empty:
        top = summary.iloc[[0]]
    top_row = top.iloc[0]

    # --- Panel a: Savage archive projection ---
    ax = fig.add_subplot(2, 2, 1)
    metrics = [
        ("Particles", "n_particles", BLUE),
        ("Pathline pts", "n_pathline_points", SKY),
        ("Ref edges", "n_endpoint_edges", GREEN),
        ("TP", "true_positive_edges", GREEN),
        ("FP", "false_positive_edges", VERMILLION),
        ("FN", "false_negative_edges", ORANGE),
    ]
    labels = [m[0] for m in metrics]
    values = [float(top_row.get(m[1], 0)) for m in metrics]
    colors = [m[2] for m in metrics]
    # Set 0 to 1 for log scale
    values = [max(1, v) for v in values]
    ax.bar(labels, values, color=colors, width=0.55, edgecolor="white", linewidth=0.4)
    ax.set_yscale("log")
    ax.set_ylabel("Count (log scale)")
    ax.set_title("Savage archive projection")
    ax.tick_params(axis="x", rotation=35)
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "a")

    # --- Panel b: Source-start point-cloud envelopes ---
    ax_b = fig.add_subplot(2, 2, 2)
    if not pathline.empty:
        targets = list(pathline["endpoint_target_node"].dropna().drop_duplicates())
        colors_pts = [GREEN, BLUE, PURPLE, ORANGE, SKY]
        for idx, target in enumerate(targets[:5]):
            subset = pathline[pathline["endpoint_target_node"] == target]
            color = colors_pts[idx % len(colors_pts)]
            ax_b.scatter(
                pd.to_numeric(subset["endpoint_start_x"], errors="coerce"),
                pd.to_numeric(subset["endpoint_start_y"], errors="coerce"),
                s=2, color=color, alpha=0.18, linewidth=0,
            )
            hull = _convex_hull_xy(subset["endpoint_start_x"], subset["endpoint_start_y"])
            if len(hull):
                closed = np.vstack([hull, hull[0]])
                ax_b.plot(closed[:, 0], closed[:, 1], color=color, linewidth=0.8,
                        label=str(target).replace("cell_", ""))
        ax_b.set_aspect("equal", adjustable="box")
        ax_b.set_xlabel("Start x")
        ax_b.set_ylabel("Start y")
        
        # Move receptor legend below panel b as a horizontal row
        handles, leg_labels = ax_b.get_legend_handles_labels()
        ax_b.legend(handles, leg_labels, frameon=False, ncol=5, loc="upper center",
                   bbox_to_anchor=(0.5, -0.22), fontsize=6.5)
                   
    ax_b.set_title("Source-start point-cloud envelopes")
    _add_gridlines(ax_b, x=True, y=True)
    _panel_label(ax_b, "b")

    # --- Panel c: Prior-transfer consistency ---
    ax_c = fig.add_subplot(2, 2, 3)
    if not time_df.empty:
        x = pd.to_numeric(time_df["particle_endpoint_time_median"], errors="coerce")
        y = pd.to_numeric(time_df["hydrosheaf_modpath_weight_mean"], errors="coerce")
        valid = (x > 0) & (y > 0)
        ax_c.scatter(x[valid], y[valid], s=14, color=BLUE, alpha=0.7, edgecolor="white", linewidth=0.3)
        if valid.any():
            lo = float(min(x[valid].min(), y[valid].min()))
            hi = float(max(x[valid].max(), y[valid].max()))
            ax_c.plot([lo, hi], [lo, hi], color=GRAY, linewidth=0.8, linestyle=(0, (3, 2)))
            ax_c.set_xscale("log")
            ax_c.set_yscale("log")
            if not time_summary.empty:
                rho = time_summary.iloc[0].get("spearman_rho", np.nan)
                tau = time_summary.iloc[0].get("kendall_tau", np.nan)
                ax_c.text(0.95, 0.05, f"rho = {float(rho):.3f}\ntau = {float(tau):.3f}",
                         transform=ax_c.transAxes, va="bottom", ha="right", fontsize=9.0,
                         bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5})
        ax_c.set_xlabel("Endpoint particle time median (days)")
        ax_c.set_ylabel("MODPATH-derived Hydrosheaf edge weight")
    ax_c.set_title("Prior-transfer consistency")
    _add_gridlines(ax_c, x=True, y=True, which="both")
    _panel_label(ax_c, "c")

    # --- Panel d: Empty ---
    ax_d = fig.add_subplot(2, 2, 4)
    ax_d.axis("off")

    _save(fig, "Fig4_External_archive_validation_and_guarded_diagnostics")


def _convex_hull_xy(x: pd.Series, y: pd.Series) -> np.ndarray:
    from scipy.spatial import ConvexHull, QhullError
    points = np.column_stack([pd.to_numeric(x, errors="coerce"),
                               pd.to_numeric(y, errors="coerce")])
    points = points[np.isfinite(points).all(axis=1)]
    points = np.unique(points, axis=0)
    if len(points) < 3:
        return np.empty((0, 2))
    try:
        hull = ConvexHull(points)
    except QhullError:
        return np.empty((0, 2))
    return points[hull.vertices]


# ---------------------------------------------------------------------------
# SUPPLEMENTARY FIGURES
# ---------------------------------------------------------------------------

    pass


def make_supp_figure_s1_workflow_provenance() -> None:
    """S1: Full M4 workflow and file provenance diagram."""
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    # fig.suptitle removed per Q1 guidelines
    ax.set_xlim(0, 11)
    ax.set_ylim(0.2, 5.8)
    ax.axis("off")

    tiers = [
        ("Data", 0.5, [
            "MODPATH archive\n(Savage, G. Miami,\nLong Island)",
            "YAML configs\n(savage.yaml,\ngreat_miami.yaml)",
        ], BLUE),
        ("Ingestion", 3.2, [
            "run_public_archive_\npipeline.py",
            "phase2_savage_\npipeline.py",
            "inspect_public_\narchive.py",
        ], GREEN),
        ("Validation", 5.9, [
            "phase2b_independent_\nvalidation.py",
            "validate_independent_\ngraph_against_modpath()",
            "edge_confusion()",
        ], PURPLE),
        ("Output", 8.6, [
            "independent_graph_\nvs_modpath.csv",
            "edge_classification.csv",
            "external_modpath_*\nresult files",
        ], ORANGE),
    ]

    for label, x, items, color in tiers:
        ax.text(x + 0.7, 5.5, label, ha="center", fontsize=8, fontweight="bold", color=color)
        for i, item in enumerate(items):
            y = 4.3 - i * 1.35
            rect = Rectangle((x, y - 0.45), 1.4, 0.85, facecolor=color, edgecolor="white",
                              linewidth=0.6, alpha=0.25)
            ax.add_patch(rect)
            ax.text(x + 0.7, y, item, ha="center", va="center", fontsize=6.5, color=DARK)

        if x > 1:
            ax.annotate("", xy=(x, 4.0), xytext=(x - 0.3, 4.0),
                        arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=1.5, mutation_scale=10.0, facecolor=GRAY))

    ax.text(5.5, 0.4, "Scripts / Results / Figures / Tables",
            fontsize=7.5, color=GRAY, ha="center", fontweight="bold")

    _save(fig, "FigS1_Full_M4_workflow_and_file_provenance")


def make_supp_figure_s2_full_scenario_networks() -> None:
    """S2: Full scenario-by-scenario network plots."""
    nodes_df = _read_csv(PUBLIC_DIR / "modpath_node_mapping.csv")
    ref_edges = _read_csv(PUBLIC_DIR / "modpath_reference_edges.csv")

    if nodes_df.empty:
        return

    pos = _build_node_positions(nodes_df)
    # sparse_node is a sensitivity analysis, not a scenario — excluded from scenario panels.
    # 9 inference/control scenarios → 3×3 grid.
    scenarios_ordered = ["spatial_only", "head_gradient", "head_gradient_bayesian_hodge",
                         "real_head_projected_gradient", "head_depth",
                         "hydrostratigraphic",
                         "negative_random", "negative_wrong_direction", "negative_shortcut"]
    scenarios_ordered = [s for s in scenarios_ordered if s in EVIDENCE_LEVELS]

    n = len(scenarios_ordered)
    cols = 3
    rows = (n + cols - 1) // cols  # 9 scenarios → 3 rows × 3 cols
    fig, axes = plt.subplots(rows, cols, figsize=(7.5, 2.5 * rows))
    fig.subplots_adjust(left=0.04, right=0.96, top=0.92, bottom=0.12, hspace=0.50, wspace=0.18)
    axes = np.atleast_1d(axes).flatten()
    # Hide unused axes
    for extra in range(n, len(axes)):
        axes[extra].set_visible(False)

    # Calculate uniform limits across all panels to make it a clean diagnostic comparison
    xs, ys = zip(*pos.values())
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xpad = (xmax - xmin) * 0.04
    ypad = (ymax - ymin) * 0.04

    for idx, scenario in enumerate(scenarios_ordered):
        ax = axes[idx]
        subset = _get_classified_edges(scenario, nodes_df, ref_edges)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(xmin - xpad, xmax + xpad)
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ax.axis("off")

        sub_ids = None
        if scenario == "sparse_node":
            rng = np.random.default_rng(20260521)
            n_half = len(pos) // 2
            node_list = list(pos.keys())
            idxs = rng.choice(len(pos), size=n_half, replace=False)
            sub_ids = {node_list[i] for i in idxs}

        if ref_edges is not None and not ref_edges.empty:
            for _, row in ref_edges.iterrows():
                u = str(row["u"])
                v = str(row["v"])
                if scenario == "sparse_node" and sub_ids is not None:
                    if u not in sub_ids or v not in sub_ids:
                        continue
                if u in pos and v in pos:
                    x1, y1 = pos[u]
                    x2, y2 = pos[v]
                    ax.plot([x1, x2], [y1, y2], color=LIGHT_GRAY, linewidth=0.15,
                            alpha=0.25, zorder=0)

        styles = {"TP": (GREEN, 0.70, 0.65), "FP": (VERMILLION, 0.60, 0.65),
                  "FN": (ORANGE, 0.60, 0.65)}
        for _, row in subset.iterrows():
            label = str(row["classification"])
            if label == "TN" or label not in styles:
                continue
            u = str(row["u"])
            v = str(row["v"])
            if u in pos and v in pos:
                color, alpha, lw = styles[label]
                x1, y1 = pos[u]
                x2, y2 = pos[v]
                ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                            arrowprops=dict(arrowstyle="-|>", color=color,
                                            lw=lw, alpha=alpha, mutation_scale=5))

        if scenario == "sparse_node" and sub_ids is not None:
            active_pos = {k: v for k, v in pos.items() if k in sub_ids}
            inactive_pos = {k: v for k, v in pos.items() if k not in sub_ids}
            if inactive_pos:
                ixs, iys = zip(*inactive_pos.values())
                ax.scatter(ixs, iys, s=0.8, facecolor="none", edgecolor=LIGHT_GRAY, linewidth=0.15, zorder=2)
            if active_pos:
                axs_val, ays = zip(*active_pos.values())
                ax.scatter(axs_val, ays, s=1.5, facecolor="white", edgecolor=GRAY, linewidth=0.2, zorder=3)
        else:
            xs, ys = zip(*pos.values())
            ax.scatter(xs, ys, s=1.5, facecolor="white", edgecolor=GRAY, linewidth=0.2, zorder=3)

        tp = int((subset["classification"] == "TP").sum())
        fp = int((subset["classification"] == "FP").sum())
        fn = int((subset["classification"] == "FN").sum())
        
        # Clean title matching user recommendation format to avoid label/title overlaps
        ax.set_title(f"{chr(97 + idx)} {_display_scenario(scenario)}\nTP = {tp}, FP = {fp}, FN = {fn}",
                     fontsize=7.2, pad=4)

    for idx in range(n, len(axes)):
        axes[idx].axis("off")

    legend_elements = [
        Line2D([0], [0], color=GREEN, lw=1.5, label="TP"),
        Line2D([0], [0], color=VERMILLION, lw=1.5, label="FP"),
        Line2D([0], [0], color=ORANGE, lw=1.5, label="FN"),
    ]
    fig.legend(handles=legend_elements, frameon=False, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 0.01), fontsize=7)

    _save(fig, "FigS2_Full_scenario_by_scenario_network_diagnostics")


def make_supp_figure_s3_prior_audit() -> None:
    """S3: MODPATH prior-mode audit."""
    priors = _read_csv(RESULT_DIR / "modpath_informed_priors.csv")
    if priors.empty:
        return

    # Reduce figsize slightly to fit the legend nicely on the right
    fig, ax = plt.subplots(figsize=(4.0, 2.6))
    fig.subplots_adjust(left=0.18, right=0.72, top=0.88, bottom=0.18)
    modes = priors["prior_mode"].astype(str)
    x = np.arange(len(priors))
    w = 0.24
    input_edges = _first_existing_column(
        priors,
        "n_hydrosheaf_input_edges",
        "n_input_hydrosheaf_edges",
    )
    ax.bar(x - w, pd.to_numeric(input_edges, errors="coerce"),
           w, color=GRAY, label="Input Hydrosheaf")
    ax.bar(x, pd.to_numeric(priors["n_modpath_prior_edges"], errors="coerce"),
           w, color=ORANGE, label="MODPATH priors")
    ax.bar(x + w, pd.to_numeric(priors["n_output_edges"], errors="coerce"),
           w, color=BLUE, label="Output graph")
    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    # y-axis label smaller: fontsize=7.5
    ax.set_ylabel("Directed edges", fontsize=7.5)
    # Title size reduced by 30-40%: fontsize=7.5
    ax.set_title("MODPATH prior-mode audit", fontsize=7.5)
    
    # Place text annotation inside or below
    ax.text(0.5, -0.25, "Prior-assisted only; excluded from independent validation.",
            transform=ax.transAxes, ha="center", va="top", fontsize=6.0, color=GRAY)
    
    _add_gridlines(ax, x=False, y=True)
    # Move legend outside the plotting area
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=6.5)
    
    _save(fig, "FigS3_MODPATH_prior_mode_audit")


def make_supp_figure_s4_pathline_complexity() -> None:
    """S4: Full pathline sequence complexity."""
    structure = _read_csv(RESULT_DIR / "external_modpath_pathline_structure.csv")
    if structure.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2),
                             gridspec_kw={"wspace": 0.42, "hspace": 0.50})
    # fig.suptitle removed per Q1 guidelines
    fig.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.07, hspace=0.50, wspace=0.42)

    ax = axes[0, 0]
    cells = pd.to_numeric(structure["n_compressed_cells"], errors="coerce").dropna()
    ax.hist(cells, bins=30, color=GRAY, edgecolor="white", linewidth=0.4)
    ax.axvline(cells.median(), color=VERMILLION, linewidth=1.0, label=f"Median={cells.median():.0f}")
    ax.set_xlabel("Compressed cells per pathline")
    ax.set_ylabel("Particles")
    ax.legend(frameon=False, fontsize=7.5) # Increased legend fontsize
    ax.set_title("Compressed pathline cell count")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "a")

    ax = axes[0, 1]
    n_points = pd.to_numeric(structure["n_pathline_points"], errors="coerce").dropna()
    ax.hist(n_points, bins=30, color=SKY, edgecolor="white", linewidth=0.4)
    ax.axvline(n_points.median(), color=VERMILLION, linewidth=1.0, label=f"Median={n_points.median():.0f}")
    ax.set_xlabel("Pathline points")
    ax.set_ylabel("Particles")
    ax.legend(frameon=False, fontsize=7.5) # Increased legend fontsize
    ax.set_title("Pathline point count per particle")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "b")

    ax = axes[1, 0]
    edge_counts = structure.groupby("endpoint_edge_id").size()
    ax.hist(edge_counts, bins=25, color=BLUE, edgecolor="white", linewidth=0.4)
    ax.axvline(edge_counts.median(), color=VERMILLION, linewidth=1.0,
               label=f"Median={edge_counts.median():.0f}")
    ax.set_xlabel("Particles per edge")
    ax.set_ylabel("Edges")
    ax.legend(frameon=False, fontsize=7.5) # Increased legend fontsize
    ax.set_title("Particles per source-receptor edge")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "c")

    ax = axes[1, 1]
    transitions = pd.to_numeric(structure["n_cell_transitions"], errors="coerce").dropna()
    ax.hist(transitions, bins=30, color=PURPLE, edgecolor="white", linewidth=0.4)
    ax.axvline(transitions.median(), color=VERMILLION, linewidth=1.0,
               label=f"Median={transitions.median():.0f}")
    ax.set_xlabel("Cell transitions")
    ax.set_ylabel("Particles")
    ax.legend(frameon=False, fontsize=7.5) # Increased legend fontsize
    ax.set_title("Cell transitions per pathline")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "d")

    _save(fig, "FigS4_Full_pathline_sequence_complexity")


def make_supp_figure_s5_travel_time_audit() -> None:
    """S5: Travel-time rank limitation and harmonisation audit."""
    rank = _read_csv(RESULT_DIR / "external_modpath_travel_time_rank.csv")
    rank_summary = _read_csv(RESULT_DIR / "external_modpath_travel_time_rank_summary.csv")
    time_df = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time.csv")
    time_summary = _read_csv(RESULT_DIR / "external_modpath_harmonized_travel_time_summary.csv")

    if rank.empty and time_df.empty:
        return

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5),
                             gridspec_kw={"wspace": 0.42, "hspace": 0.52})
    # fig.suptitle removed per Q1 guidelines
    fig.subplots_adjust(left=0.10, right=0.95, top=0.95, bottom=0.07, hspace=0.52, wspace=0.42)

    # Panel a: Raw endpoint vs pathline time
    ax = axes[0, 0]
    if not rank.empty:
        x = pd.to_numeric(rank["endpoint_travel_time_mean"], errors="coerce")
        y = pd.to_numeric(rank["pathline_elapsed_time_median"], errors="coerce")
        valid = (x > 0) & (y > 0)
        ax.scatter(x[valid], y[valid], s=14, color=BLUE, alpha=0.7, edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Endpoint travel time mean (days)")
        ax.set_ylabel("Pathline elapsed time median (days)")
        if not rank_summary.empty:
            top = rank_summary.iloc[0]
            rho = top.get("spearman_rho", np.nan)
            tau = top.get("kendall_tau", np.nan)
            ax.text(0.95, 0.95, f"Raw\nρ = {float(rho):.3f}\nτ = {float(tau):.3f}",
                    transform=ax.transAxes, va="top", ha="right", fontsize=9.0,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    ax.set_title("Raw travel-time rank (discordant)")
    _add_gridlines(ax, x=True, y=True, which="both")
    _panel_label(ax, "a")

    # Panel b: Harmonised endpoint time vs edge weights
    ax = axes[0, 1]
    if not time_df.empty:
        x = pd.to_numeric(time_df["particle_endpoint_time_median"], errors="coerce")
        y = pd.to_numeric(time_df["hydrosheaf_modpath_weight_mean"], errors="coerce")
        valid = (x > 0) & (y > 0)
        ax.scatter(x[valid], y[valid], s=14, color=GREEN, alpha=0.7, edgecolor="white", linewidth=0.3)
        lo = float(min(x[valid].min(), y[valid].min()))
        hi = float(max(x[valid].max(), y[valid].max()))
        ax.plot([lo, hi], [lo, hi], color=GRAY, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Endpoint particle time median (days)")
        ax.set_ylabel("Hydrosheaf MODPATH-derived edge weight")
        if not time_summary.empty:
            top = time_summary.iloc[0]
            rho = top.get("spearman_rho", np.nan)
            tau = top.get("kendall_tau", np.nan)
            ax.text(0.95, 0.05, f"Harmonised\nρ = {float(rho):.3f}\nτ = {float(tau):.3f}",
                    transform=ax.transAxes, va="bottom", ha="right", fontsize=9.0,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    ax.set_title("Harmonised travel-time agreement")
    _add_gridlines(ax, x=True, y=True, which="both")
    _panel_label(ax, "b")

    # Panel c: Rank residual comparison
    ax = axes[1, 0]
    if not rank.empty and not time_df.empty:
        # Compare raw vs harmonised spearman rho
        raw_rho = float(rank_summary.iloc[0].get("spearman_rho", 0)) if not rank_summary.empty else 0
        harm_rho = float(time_summary.iloc[0].get("spearman_rho", 0)) if not time_summary.empty else 0
        ax.bar(["Raw time\n(discordant)", "Harmonised\n(concordant)"], [raw_rho, harm_rho],
               color=[VERMILLION, GREEN], width=0.45)
        ax.axhline(0, color=GRAY, linewidth=0.5)
        ax.set_ylabel("Spearman ρ")
        ax.set_title("Rank agreement: raw vs harmonised")
        for i, v in enumerate([raw_rho, harm_rho]):
            ax.text(i, v + 0.03 * np.sign(v), f"{v:.3f}", ha="center", fontsize=6, fontweight="bold")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "c")

    # Panel d: Harmonised residual histogram
    ax = axes[1, 1]
    if not time_df.empty:
        x = pd.to_numeric(time_df["particle_endpoint_time_median"], errors="coerce")
        y = pd.to_numeric(time_df["hydrosheaf_modpath_weight_mean"], errors="coerce")
        valid = (x > 0) & (y > 0)
        diff = np.log10(y[valid]) - np.log10(x[valid])
        ax.hist(diff, bins=25, color=GRAY, edgecolor="white", linewidth=0.4)
        ax.axvline(0.0, color=VERMILLION, linewidth=1.0)
        ax.set_xlabel("log10(weight / endpoint time)")
        ax.set_ylabel("Edges")
        mad = float(time_summary.iloc[0].get("median_abs_log10_difference", np.nan)) if not time_summary.empty else np.nan
        ax.text(0.95, 0.95, f"Median |Δ| = {mad:.4f}", transform=ax.transAxes,
                va="top", ha="right", fontsize=9.0,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    ax.set_title("Harmonised scale residual")
    _add_gridlines(ax, x=False, y=True)
    _panel_label(ax, "d")

    _save(fig, "FigS5_Travel_time_rank_limitation_and_harmonisation_audit")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print("Generating main-text figures...")
    make_figure1_workflow()
    make_figure2_performance()
    make_figure3_edge_networks()
    make_figure4_external_archive()

    print("Generating supplementary figures...")
    make_supp_figure_s1_workflow_provenance()
    make_supp_figure_s2_full_scenario_networks()
    make_supp_figure_s3_prior_audit()
    make_supp_figure_s4_pathline_complexity()
    make_supp_figure_s5_travel_time_audit()

    print(f"All figures written to {FIGURE_DIR}")


if __name__ == "__main__":
    sys.exit(main())
