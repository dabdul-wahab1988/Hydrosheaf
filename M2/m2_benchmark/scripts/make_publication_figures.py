from __future__ import annotations

import math
from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import yaml


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = BENCHMARK_ROOT / "figures"
RESULT_DIR = BENCHMARK_ROOT / "results"
EXTERNAL_DIR = BENCHMARK_ROOT / "external"


COLORS = {
    "young": "#3b82f6",
    "intermediate": "#22c55e",
    "old": "#f59e0b",
    "fossil": "#ef4444",
    "mixed": "#8b5cf6",
    "line": "#334155",
    "muted": "#64748b",
}


def _save(fig: plt.Figure, name: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )


def _wrap_labels(labels: list[str], width: int = 14) -> list[str]:
    return [fill(str(label), width=width) for label in labels]


def _short_node(label: str) -> str:
    replacements = {
        "Virtual_Endmember_PC1_High": "VE1H",
        "Virtual_Endmember_PC1_Low": "VE1L",
        "Virtual_Endmember_PC2_High": "VE2H",
        "Virtual_Endmember_PC2_Low": "VE2L",
    }
    out = str(label)
    for source, target in replacements.items():
        out = out.replace(source, target)
    if out.startswith("NGW-"):
        out = out.replace("NGW-", "")
    return out


def load_truth() -> dict:
    with (BENCHMARK_ROOT / "config" / "ground_truth.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def plot_compact_process_network() -> None:
    truth = load_truth()
    graph = nx.DiGraph()
    node_classes = {}
    for node_id, node in truth["nodes"].items():
        graph.add_node(node_id)
        node_classes[node_id] = node.get("age_class", "intermediate")
    for node_id, node in truth["generated_nodes"].items():
        graph.add_node(node_id)
        node_classes[node_id] = node.get("age_class", "intermediate")

    edge_labels = {}
    for edge in truth["generation_edges"]:
        graph.add_edge(edge["u"], edge["v"], process=edge["process"])
        edge_labels[(edge["u"], edge["v"])] = {
            "evaporation": "E",
            "carbonate_weathering": "C",
            "evaporite_weathering": "V",
            "gypsum_halite": "V",
            "lateral_mixing": "L",
            "dilution_recharge_pulse": "D",
            "nitrate_recharge": "N",
            "confined_exchange": "X",
            "regional_discharge": "R",
            "deep_discharge": "R",
        }.get(edge["process"], edge["process"][:1].upper())
    for edge in truth.get("lateral_truth_edges", []):
        graph.add_edge(edge["u"], edge["v"], process="lateral")
        edge_labels[(edge["u"], edge["v"])] = "Lateral"

    pos = {
        "MODERN": (0.02, 0.35),
        "RCH": (0.14, 0.52),
        "SH15": (0.28, 0.68),
        "SH30": (0.28, 0.36),
        "INT1": (0.43, 0.58),
        "INT2": (0.56, 0.52),
        "MIX20": (0.70, 0.56),
        "MIX40": (0.82, 0.52),
        "MIX60": (0.94, 0.47),
        "DILUTE": (1.07, 0.43),
        "DEEP": (1.20, 0.39),
        "DISCH": (1.34, 0.46),
        "LAT_SAL": (0.70, 0.83),
        "LAT_CARB": (0.94, 0.16),
    }

    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    ax.set_title("Figure 2. Locked synthetic process network", fontsize=14, fontweight="bold", pad=12)
    node_colors = [COLORS.get(node_classes.get(node, ""), COLORS["muted"]) for node in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=920, edgecolors="#1f2937", linewidths=1.3, ax=ax)
    nx.draw_networkx_labels(graph, pos, font_size=8.5, font_weight="bold", ax=ax)
    nx.draw_networkx_edges(
        graph,
        pos,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=14,
        width=1.8,
        edge_color=COLORS["line"],
        connectionstyle="arc3,rad=0.08",
        ax=ax,
    )

    label_pos = {
        ("RCH", "SH15"): (0.20, 0.61),
        ("RCH", "SH30"): (0.20, 0.43),
        ("SH15", "INT1"): (0.36, 0.63),
        ("INT1", "INT2"): (0.50, 0.55),
        ("INT2", "MIX20"): (0.63, 0.54),
        ("MIX20", "MIX40"): (0.76, 0.54),
        ("MIX40", "MIX60"): (0.88, 0.50),
        ("MIX60", "DILUTE"): (1.00, 0.45),
        ("DILUTE", "DEEP"): (1.14, 0.41),
        ("DEEP", "DISCH"): (1.27, 0.43),
        ("LAT_SAL", "MIX20"): (0.68, 0.70),
        ("LAT_SAL", "MIX40"): (0.79, 0.68),
        ("LAT_CARB", "MIX60"): (0.92, 0.30),
    }
    for edge, label in edge_labels.items():
        if edge in label_pos:
            ax.text(
                *label_pos[edge],
                label,
                fontsize=7.0,
                color="#1f2937",
                ha="center",
                va="center",
                bbox=dict(facecolor="white", edgecolor="#d1d5db", alpha=0.88, boxstyle="round,pad=0.18"),
            )

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=label, markerfacecolor=color, markeredgecolor="#1f2937", markersize=8)
        for label, color in [
            ("young", COLORS["young"]),
            ("intermediate", COLORS["intermediate"]),
            ("old", COLORS["old"]),
            ("mixed", COLORS["mixed"]),
            ("fossil", COLORS["fossil"]),
        ]
    ]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05), frameon=False, ncol=5, fontsize=8)
    ax.text(
        0.02,
        0.07,
        "Edge codes: E evaporation, C carbonate, V evaporite, L lateral, D dilution, X exchange, R discharge",
        fontsize=7.5,
        color="#374151",
        transform=ax.transAxes,
    )
    ax.set_xlim(-0.06, 1.42)
    ax.set_ylim(0.05, 0.93)
    ax.axis("off")
    _save(fig, "figure2_process_network.png")


def plot_external_validation_composite() -> None:
    age = pd.read_csv(EXTERNAL_DIR / "usgs_age" / "results" / "usgs_age_validation.csv")
    modpath = pd.read_csv(EXTERNAL_DIR / "modpath" / "results" / "modpath_topology_summary.csv").iloc[0]
    phreeqc = pd.read_csv(RESULT_DIR / "phreeqc_forward_validation.csv")
    dgmeta = pd.read_csv(EXTERNAL_DIR / "dgmeta" / "results" / "dgmeta_hydrosheaf_comparison.csv")

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))

    ax = axes[0, 0]
    plot_age = age[(age["reference_mean_age_years"] > 0) & (age["hydrosheaf_age_years"] > 0)].copy()
    counts = plot_age["supported_tracer_count"].clip(upper=5)
    sc = ax.scatter(
        plot_age["reference_mean_age_years"],
        plot_age["hydrosheaf_age_years"],
        c=counts,
        cmap="viridis",
        s=18,
        alpha=0.60,
        edgecolors="none",
    )
    lim_min = 0.1
    lim_max = max(plot_age["reference_mean_age_years"].max(), plot_age["hydrosheaf_age_years"].max()) * 1.15
    ax.plot([lim_min, lim_max], [lim_min, lim_max], color="#111827", lw=1.1)
    ax.plot([lim_min, lim_max], [lim_min * 10, lim_max * 10], color="#9ca3af", lw=0.8, ls="--")
    ax.plot([lim_min * 10, lim_max], [lim_min, lim_max / 10], color="#9ca3af", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)
    ax.set_title("USGS public tracer-age comparison")
    ax.set_xlabel("Published mean age (yr)")
    ax.set_ylabel("Hydrosheaf age (yr)")
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Supported tracers")
    _panel_label(ax, "A")

    ax = axes[0, 1]
    labels = ["TP", "FP", "FN"]
    values = [modpath["true_positive_edges"], modpath["false_positive_edges"], modpath["false_negative_edges"]]
    ax.bar(labels, values, color=["#2563eb", "#ef4444", "#f97316"])
    ax.set_title("MODPATH topology agreement")
    ax.set_ylabel("Directed edges")
    ax.text(
        0.96,
        0.88,
        f"F1={modpath['edge_f1']:.2f}\ndirection={modpath['direction_agreement_rate']:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="#d1d5db", boxstyle="round,pad=0.35"),
    )
    _panel_label(ax, "B")

    ax = axes[1, 0]
    ax.scatter(phreeqc["rmse_mmolL"], phreeqc["nse"], s=18, alpha=0.55, color="#0f766e", edgecolors="none")
    ax.axhline(0.5, color="#ef4444", lw=1.0, ls="--")
    ax.axhline(0.0, color="#9ca3af", lw=0.8)
    ax.set_title("Live PHREEQC forward diagnostics")
    ax.set_xlabel("RMSE (mmol/L)")
    ax.set_ylabel("NSE")
    ax.text(
        0.96,
        0.10,
        f"n={len(phreeqc)}\nmedian NSE={phreeqc['nse'].median():.2f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        bbox=dict(facecolor="white", edgecolor="#d1d5db", boxstyle="round,pad=0.35"),
    )
    _panel_label(ax, "C")

    ax = axes[1, 1]
    dg = dgmeta.dropna(subset=["reference_temperature_c", "hydrosheaf_temperature_c"]).copy()
    color_map = {"UA": "#2563eb", "CE": "#f59e0b"}
    for model, group in dg.groupby("hydrosheaf_model"):
        ax.scatter(
            group["reference_temperature_c"],
            group["hydrosheaf_temperature_c"],
            s=18,
            alpha=0.60,
            label=str(model),
            color=color_map.get(str(model), COLORS["muted"]),
            edgecolors="none",
        )
    lo = min(dg["reference_temperature_c"].min(), dg["hydrosheaf_temperature_c"].min()) - 1
    hi = max(dg["reference_temperature_c"].max(), dg["hydrosheaf_temperature_c"].max()) + 1
    ax.plot([lo, hi], [lo, hi], color="#111827", lw=1.0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_title("DGMETA recharge temperature")
    ax.set_xlabel("DGMETA reference (C)")
    ax.set_ylabel("Hydrosheaf estimate (C)")
    ax.legend(frameon=False, title="model", fontsize=8)
    _panel_label(ax, "D")

    fig.tight_layout()
    _save(fig, "figure4_external_validation_composite.png")


def plot_residence_time_supplement() -> None:
    ages = pd.read_csv(RESULT_DIR / "age_inference_validation.csv")
    consistency = pd.read_csv(RESULT_DIR / "age_network_consistency.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))

    ax = axes[0]
    sample = ages.groupby("site_id").median(numeric_only=True).reset_index()
    ax.scatter(sample["true_mrt_years"], sample["single_node_lpm_years"], label="single-node", color="#94a3b8", s=35)
    ax.scatter(sample["true_mrt_years"], sample["network_bayesian_years"], label="network", color="#2563eb", s=35)
    lim = [1, max(sample["true_mrt_years"].max(), sample["single_node_lpm_years"].max()) * 1.1]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Synthetic age agreement")
    ax.set_xlabel("True MRT (yr)")
    ax.set_ylabel("Estimated MRT (yr)")
    ax.legend(frameon=False, fontsize=8)
    _panel_label(ax, "A")

    ax = axes[1]
    labels = ["RCH", "DILUTE", "DEEP"]
    data = [ages.loc[ages["site_id"] == label, "network_bayesian_years"] for label in labels]
    ax.violinplot(data, positions=np.arange(len(labels)), showmeans=True)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_title("Network posteriors")
    ax.set_ylabel("MRT (yr)")
    _panel_label(ax, "B")

    ax = axes[2]
    grouped = consistency.groupby("edge_confidence_threshold").agg(
        consistency=("downstream_age_consistency_fraction", "median"),
        tau=("kendall_tau", "median"),
    )
    ax.plot(grouped.index, grouped["consistency"], "o-", label="age-order fraction", color="#2563eb")
    ax.plot(grouped.index, grouped["tau"], "s--", label="Kendall tau", color="#f59e0b")
    ax.set_ylim(-0.1, 1.05)
    ax.set_title("Downstream age consistency")
    ax.set_xlabel("Edge confidence threshold")
    ax.legend(frameon=False, fontsize=8)
    _panel_label(ax, "C")

    fig.tight_layout()
    _save(fig, "figure_s4_residence_time_network_update.png")


def plot_sensitivity_figure() -> None:
    unc = pd.read_csv(RESULT_DIR / "sensitivity_uncertainty_summary.csv")
    unc = unc.sort_values("iqr_contribution", ascending=True)
    unc["relative_percent"] = 100.0 * unc["iqr_contribution"] / unc["iqr_contribution"].sum()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.05, 1.0]})
    ax = axes[0]
    ax.barh(unc["factor"], unc["relative_percent"], color="#475569")
    ax.set_title("Relative uncertainty contribution")
    ax.set_xlabel("Share of summed contribution (%)")
    ax.set_xlim(0, max(unc["relative_percent"].max() * 1.12, 10))
    for y, value in enumerate(unc["relative_percent"]):
        ax.text(value + 1.0, y, f"{value:.1f}%", va="center", fontsize=8)
    _panel_label(ax, "A")

    ax = axes[1]
    ax.axis("off")
    ax.set_title("Uncertainty cascade")
    steps = [
        ("Tracer and\nchemistry errors", 0.16, 0.72),
        ("Age posterior\nand graph weights", 0.50, 0.72),
        ("Transport and\nreaction choice", 0.84, 0.72),
        ("PHREEQC and\nresidual filters", 0.50, 0.42),
        ("Final edge\nconfidence", 0.50, 0.14),
    ]
    for text, x, y in steps:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#334155", linewidth=1.1),
        )
    for start, end in [
        ((0.27, 0.72), (0.39, 0.72)),
        ((0.61, 0.72), (0.73, 0.72)),
        ((0.78, 0.63), (0.58, 0.48)),
        ((0.50, 0.32), (0.50, 0.23)),
    ]:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.4, color="#111827"))
    _panel_label(ax, "B")
    fig.tight_layout()
    _save(fig, "figure5_sensitivity_uncertainty.png")


def plot_sparse_algorithm() -> None:
    fig, ax = plt.subplots(figsize=(10.5, 2.7))
    ax.axis("off")
    steps = ["Transport fit", "Residual vector", "L1 reaction fit", "SI feasibility", "Objective score", "Ranked pathway"]
    xs = np.linspace(0.08, 0.92, len(steps))
    for x, step in zip(xs, steps):
        ax.text(
            x,
            0.52,
            step,
            ha="center",
            va="center",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.32", facecolor="#f8fafc", edgecolor="#0f766e", linewidth=1.2),
        )
    for x1, x2 in zip(xs[:-1], xs[1:]):
        ax.annotate("", xy=(x2 - 0.055, 0.52), xytext=(x1 + 0.055, 0.52), arrowprops=dict(arrowstyle="->", lw=1.3))
    ax.set_title("Figure S1. Sparse inverse hydrogeochemical fitting algorithm", fontsize=13, pad=10)
    _save(fig, "figure_s1_sparse_fitting_algorithm.png")


def plot_e6_nonlinear_validation() -> None:
    e6 = pd.read_csv(EXTERNAL_DIR / "e6_nonlinear" / "results" / "e6_nonlinear_validation.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    ax.hist(e6["nse"], bins=14, color="#2563eb", alpha=0.78)
    ax.axvline(e6["nse"].median(), color="#111827", lw=1.2, ls="--", label=f"median={e6['nse'].median():.2f}")
    ax.set_title("Nonlinear PHREEQC NSE")
    ax.set_xlabel("NSE")
    ax.set_ylabel("Edges")
    ax.legend(frameon=False, fontsize=8)
    _panel_label(ax, "A")

    ax = axes[1]
    minerals = [
        ("calcite", "#2563eb"),
        ("gypsum", "#f59e0b"),
        ("halite", "#0f766e"),
    ]
    all_true = []
    all_inf = []
    for mineral, color in minerals:
        x = e6[f"true_{mineral}"]
        y = e6[f"inf_{mineral}"]
        all_true.extend(x.tolist())
        all_inf.extend(y.tolist())
        ax.scatter(x, y, s=22, alpha=0.65, color=color, label=mineral)
    lim = max(max(map(abs, all_true)), max(map(abs, all_inf)), 0.01) * 1.15
    ax.plot([-lim, lim], [-lim, lim], color="#111827", lw=1.0)
    ax.axhline(0, color="#d1d5db", lw=0.8)
    ax.axvline(0, color="#d1d5db", lw=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_title("Reaction extent recovery")
    ax.set_xlabel("PHREEQC truth extent")
    ax.set_ylabel("Hydrosheaf inferred extent")
    ax.legend(frameon=False, fontsize=8)
    _panel_label(ax, "B")

    fig.tight_layout()
    _save(fig, "figure_s5_e6_nonlinear_validation.png")


def plot_dgmeta_diagnostics() -> None:
    dg = pd.read_csv(EXTERNAL_DIR / "dgmeta" / "results" / "dgmeta_hydrosheaf_comparison.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))

    ax = axes[0]
    temp = dg.dropna(subset=["reference_temperature_c", "hydrosheaf_temperature_c"])
    for model, group in temp.groupby("hydrosheaf_model"):
        ax.scatter(group["reference_temperature_c"], group["hydrosheaf_temperature_c"], s=18, alpha=0.60, label=str(model), edgecolors="none")
    lo = min(temp["reference_temperature_c"].min(), temp["hydrosheaf_temperature_c"].min()) - 1
    hi = max(temp["reference_temperature_c"].max(), temp["hydrosheaf_temperature_c"].max()) + 1
    ax.plot([lo, hi], [lo, hi], color="#111827", lw=1.0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_title("Recharge temperature")
    ax.set_xlabel("DGMETA reference (C)")
    ax.set_ylabel("Hydrosheaf (C)")
    ax.legend(frameon=False, fontsize=8)
    _panel_label(ax, "A")

    ax = axes[1]
    air = dg.dropna(subset=["reference_excess_air_cm3kg", "hydrosheaf_excess_air_cm3kg"])
    for model, group in air.groupby("hydrosheaf_model"):
        ax.scatter(group["reference_excess_air_cm3kg"], group["hydrosheaf_excess_air_cm3kg"], s=18, alpha=0.60, label=str(model), edgecolors="none")
    if not air.empty:
        lo = min(air["reference_excess_air_cm3kg"].min(), air["hydrosheaf_excess_air_cm3kg"].min()) - 0.5
        hi = max(air["reference_excess_air_cm3kg"].max(), air["hydrosheaf_excess_air_cm3kg"].max()) + 0.5
        ax.plot([lo, hi], [lo, hi], color="#111827", lw=1.0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    ax.set_title("Excess air")
    ax.set_xlabel("DGMETA reference (cm3/kg)")
    ax.set_ylabel("Hydrosheaf (cm3/kg)")
    _panel_label(ax, "B")

    ax = axes[2]
    resid = dg.dropna(subset=["temperature_error_c"]).copy()
    groups = [group["temperature_error_c"].abs() for _, group in resid.groupby("hydrosheaf_model")]
    labels = [str(model) for model, _ in resid.groupby("hydrosheaf_model")]
    ax.boxplot(groups, tick_labels=labels, showfliers=False)
    ax.set_title("Absolute temperature residuals")
    ax.set_ylabel("|error| (C)")
    _panel_label(ax, "C")

    fig.tight_layout()
    _save(fig, "figure_s6_dgmeta_diagnostics.png")


def plot_e1_residual_stratification() -> None:
    age = pd.read_csv(EXTERNAL_DIR / "usgs_age" / "results" / "usgs_age_validation.csv")
    age = age.dropna(subset=["log10_error"]).copy()
    age["age_class"] = pd.cut(
        age["reference_mean_age_years"],
        bins=[0, 70, 1000, math.inf],
        labels=["young (<70 yr)", "intermediate", "old (>1000 yr)"],
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4))

    ax = axes[0]
    count_labels = sorted(age["supported_tracer_count"].dropna().unique())
    data = [age.loc[age["supported_tracer_count"] == count, "log10_error"] for count in count_labels]
    ax.boxplot(data, tick_labels=[str(int(count)) for count in count_labels], showfliers=False)
    ax.axhline(0, color="#111827", lw=1.0)
    ax.set_title("Residual by tracer count")
    ax.set_xlabel("Supported tracer count")
    ax.set_ylabel("log10 age error")
    _panel_label(ax, "A")

    ax = axes[1]
    class_labels = ["young (<70 yr)", "intermediate", "old (>1000 yr)"]
    data = [age.loc[age["age_class"] == label, "log10_error"] for label in class_labels]
    ax.boxplot(data, tick_labels=_wrap_labels(class_labels, 12), showfliers=False)
    ax.axhline(0, color="#111827", lw=1.0)
    ax.set_title("Residual by age class")
    ax.set_ylabel("log10 age error")
    _panel_label(ax, "B")

    fig.tight_layout()
    _save(fig, "figure_s7_e1_tracer_age_residuals.png")


def plot_e4_pilot_network() -> None:
    northern_path = EXTERNAL_DIR / "northern_ghana" / "results" / "northern_ghana_edge_results.csv"
    public_path = EXTERNAL_DIR / "usgs_public_chem" / "results" / "public_chem_edge_results.csv"
    pilot_path = EXTERNAL_DIR / "pilot" / "results" / "pilot_edge_results.csv"

    if northern_path.exists():
        path = northern_path
        title_suffix = " (Corrected Northern Ghana)"
    elif public_path.exists():
        path = public_path
        title_suffix = " (USGS Public NAWQA)"
    elif pilot_path.exists():
        path = pilot_path
        title_suffix = " (Pilot Dataset)"
    else:
        print("Skipping Figure S8: No edge results found.")
        return

    pilot = pd.read_csv(path)
    top = pilot.nsmallest(15, "objective_score").copy()
    graph = nx.DiGraph()
    for _, row in top.iterrows():
        graph.add_edge(
            row["u"],
            row["v"],
            weight=float(row["chemistry_r2"]),
            score=float(row["objective_score"]),
            season=row.get("season", ""),
        )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1.2, 1.0]})
    ax = axes[0]
    pos = nx.spring_layout(graph, seed=11, k=1.05, iterations=150)
    node_colors = ["#e0f2fe" if str(node).startswith("Virtual") else "#dcfce7" for node in graph.nodes]
    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, edgecolors="#1f2937", node_size=620, ax=ax)
    if "season" in top.columns:
        season_colors = {"Wet": "#2563eb", "Dry": "#f59e0b"}
        edge_colors = [season_colors.get(graph.edges[edge].get("season"), "#334155") for edge in graph.edges]
        handles = [
            plt.Line2D([0], [0], color=color, lw=2.2, label=season)
            for season, color in season_colors.items()
        ]
        ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=8)
    else:
        edge_colors = "#334155"
    nx.draw_networkx_edges(graph, pos, arrows=True, arrowstyle="-|>", arrowsize=12, edge_color=edge_colors, width=1.2, ax=ax)
    labels = {node: _short_node(str(node)) for node in graph.nodes}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=6.2, ax=ax)
    ax.set_title(f"Top sparse-data edges{title_suffix}")
    ax.axis("off")
    _panel_label(ax, "A")

    ax = axes[1]
    top10 = top.head(10).sort_values("chemistry_r2")
    edge_labels = [f"{_short_node(row.u)} -> {_short_node(row.v)}" for row in top10.itertuples()]
    ax.barh(edge_labels, top10["chemistry_r2"], color="#0f766e")
    lower = max(0.0, min(0.90, float(top10["chemistry_r2"].min()) - 0.01))
    ax.set_xlim(lower, 1.0)
    ax.set_title("Chemistry fit for top edges")
    ax.set_xlabel("Chemistry R2")
    _panel_label(ax, "B")

    fig.tight_layout()
    _save(fig, "figure_s8_e4_sparse_pilot_network.png")


def main() -> None:
    plot_compact_process_network()
    plot_external_validation_composite()
    plot_residence_time_supplement()
    plot_sensitivity_figure()
    plot_sparse_algorithm()
    plot_e6_nonlinear_validation()
    plot_dgmeta_diagnostics()
    plot_e1_residual_stratification()
    plot_e4_pilot_network()
    print("Publication figures regenerated.")


if __name__ == "__main__":
    main()
