"""Create data-driven M3 manuscript figures."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.nuclear.input_history import build_default_tritium_input
from hydrosheaf.nuclear.multi_tracer import build_atmospheric_tracer_input


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = BENCHMARK_ROOT / "figures" / "Manuscript_Ready"
RESULT_DIR = BENCHMARK_ROOT / "results"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    }
)


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _save(fig: plt.Figure, name: str) -> None:
    out = FIGURE_DIR / name
    tmp = out.with_suffix(out.suffix + ".tmp.png")
    fig.savefig(tmp, dpi=300, bbox_inches="tight")
    tmp.replace(out)
    plt.close(fig)


def _primary_results() -> pd.DataFrame:
    primary = _read_csv(RESULT_DIR / "m3_usgs_benchmark_results.csv")
    if not primary.empty:
        return primary
    design = _read_csv(RESULT_DIR / "m3_design_matrix_results.csv")
    if design.empty:
        return design
    if "screened_dgm_gases" in set(design["scenario_id"]):
        return design[design["scenario_id"] == "screened_dgm_gases"].copy()
    if "parity_reported_corrected" in set(design["scenario_id"]):
        return design[design["scenario_id"] == "parity_reported_corrected"].copy()
    return design.copy()


def plot_fig1_atmospheric_histories() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.4))
    histories = [
        ("3H", build_default_tritium_input(), "TU", "#c2410c"),
        ("SF6", build_atmospheric_tracer_input("SF6"), "pptv", "#1d4ed8"),
        ("CFC-12", build_atmospheric_tracer_input("CFC12"), "pptv", "#15803d"),
    ]
    for ax, (label, hist, unit, color) in zip(axes, histories):
        ax.plot(hist.years, hist.values, color=color, lw=2.0)
        ax.fill_between(hist.years, hist.values, color=color, alpha=0.12)
        ax.set_title(label, loc="left", fontweight="bold")
        ax.set_xlabel("Recharge year")
        ax.set_ylabel(unit)
        ax.grid(True, ls=":", alpha=0.4)
    axes[0].set_yscale("log")
    _save(fig, "Manuscript_Fig1_Atmospheric_Histories.png")


def plot_fig2_design_matrix_performance() -> None:
    summary = _read_csv(RESULT_DIR / "m3_design_matrix_summary.csv")
    if summary.empty:
        return
    summary = summary.sort_values("median_abs_log10_error")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.4, 1.0]})
    y = np.arange(len(summary))
    axes[0].barh(y, summary["median_abs_log10_error"], color="#2563eb")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(summary["scenario_id"])
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Median absolute log10 error")
    axes[0].set_title("a  Scenario accuracy", loc="left", fontweight="bold")

    axes[1].scatter(
        summary["within_factor_2"],
        summary["within_factor_10"],
        s=np.maximum(summary["metric_rows"], 1) * 10,
        color="#0f766e",
        alpha=0.75,
        edgecolor="white",
    )
    for _, row in summary.iterrows():
        axes[1].annotate(str(row["scenario_id"]).replace("_", "\n"), (row["within_factor_2"], row["within_factor_10"]), fontsize=7)
    axes[1].set_xlabel("Within factor 2")
    axes[1].set_ylabel("Within factor 10")
    axes[1].set_xlim(-0.03, 1.03)
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_title("b  Practical accuracy", loc="left", fontweight="bold")
    axes[1].grid(True, ls=":", alpha=0.4)
    _save(fig, "Manuscript_Fig2_Design_Matrix_Performance.png")


def plot_fig3_parity_plot() -> None:
    df = _primary_results().dropna(subset=["ref_age", "est_age_multi"])
    df = df[(df["ref_age"] > 0) & (df["est_age_multi"] > 0)].copy()
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, 5.8))
    classes = sorted(df["age_class"].dropna().unique())
    cmap = plt.get_cmap("tab10")
    for idx, age_class in enumerate(classes):
        subset = df[df["age_class"] == age_class]
        ax.scatter(
            subset["ref_age"],
            subset["est_age_multi"],
            s=42,
            alpha=0.78,
            color=cmap(idx),
            edgecolor="white",
            linewidth=0.5,
            label=str(age_class).replace("_", " "),
        )
    lims = [0.1, max(df["ref_age"].max(), df["est_age_multi"].max()) * 1.2]
    ax.plot(lims, lims, color="black", ls="--", lw=1.3)
    ax.fill_between(lims, np.array(lims) / 2.0, np.array(lims) * 2.0, color="#94a3b8", alpha=0.18)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("USGS reference age (yr)")
    ax.set_ylabel("Hydrosheaf estimated age (yr)")
    metric = df["log10_error"].dropna()
    ax.text(
        0.04,
        0.96,
        f"N = {len(df)}\nMedian |log10 error| = {metric.median():.3f}\nWithin factor 2 = {df['within_factor_2'].mean():.2f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.9},
    )
    ax.legend(frameon=False, loc="lower right")
    ax.grid(True, ls=":", alpha=0.4)
    _save(fig, "Manuscript_Fig3_USGS_Benchmark_Parity.png")


def plot_fig4_real_graph_benchmark() -> None:
    graph = _read_csv(RESULT_DIR / "m3_real_usgs_graph_benchmark.csv")
    if graph.empty:
        return
    graph = graph[graph["prior_strength"] != "none"].copy()
    graph = graph.sort_values(["control_type", "delta_rmse_graph_minus_single"])
    labels = graph["graph_family"] + "\n" + graph["prior_strength"]
    colors = np.where(graph["control_type"].eq("negative_control"), "#dc2626", "#2563eb")
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.axhline(0.0, color="black", lw=1.0)
    ax.bar(np.arange(len(graph)), graph["delta_rmse_graph_minus_single"], color=colors, alpha=0.82)
    ax.set_xticks(np.arange(len(graph)))
    ax.set_xticklabels(labels, rotation=55, ha="right")
    ax.set_ylabel("Delta RMSE, graph - single (log10 yr)")
    ax.set_title("Real-USGS graph-prior benchmark", loc="left", fontweight="bold")
    _save(fig, "Manuscript_Fig4_Real_USGS_Graph_Benchmark.png")


def plot_fig5_age_class_diagnostics() -> None:
    df = _primary_results().dropna(subset=["log10_error", "age_class"])
    if df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    age_classes = sorted(df["age_class"].unique())
    data = [df.loc[df["age_class"] == age_class, "log10_error"].dropna() for age_class in age_classes]
    axes[0].boxplot(data, tick_labels=age_classes, showfliers=False)
    axes[0].set_ylabel("Absolute log10 error")
    axes[0].set_title("a  Error by age class", loc="left", fontweight="bold")
    axes[0].tick_params(axis="x", rotation=35)

    diag = df.dropna(subset=["young_gas_proxy_coherence", "log10_error"])
    axes[1].scatter(
        diag["young_gas_proxy_coherence"],
        diag["log10_error"],
        c=np.log10(np.maximum(diag["ref_age"], 0.1)),
        cmap="viridis",
        s=42,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.4,
    )
    axes[1].set_xlabel("Young-tracer proxy coherence (log10-age SD)")
    axes[1].set_ylabel("Absolute log10 error")
    axes[1].set_title("b  Discordance diagnostic", loc="left", fontweight="bold")
    axes[1].grid(True, ls=":", alpha=0.4)
    _save(fig, "Manuscript_Fig5_Age_Class_Diagnostics.png")


def main() -> int:
    plot_fig1_atmospheric_histories()
    plot_fig2_design_matrix_performance()
    plot_fig3_parity_plot()
    plot_fig4_real_graph_benchmark()
    plot_fig5_age_class_diagnostics()
    print(f"Wrote M3 manuscript figures to {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
