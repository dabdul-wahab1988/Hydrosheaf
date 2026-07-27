"""Create data-driven M3 manuscript figures."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
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
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 10,
        "figure.dpi": 150,
    }
)


LABEL_MAP = {
    # ── Scenario IDs ────────────────────────────────────────────────────
    "tracerlpm_strict_parity":             "Strict Parity",
    "tracerlpm_parity_hier_oldwater":      "Parity + Hierarchical Old-Water Prior",
    "tracerlpm_parity_agefractions":       "Parity + Age Fractions",
    "screened_dgm_gases":                  "Screened Young Gas",
    "parity_reported_corrected":           "Reported-Model Parity",
    "hydrosheaf_selection_corrected":      "Hydrosheaf Selection",
    "reported":                            "Reported-Model Parity",
    "selection":                           "Hydrosheaf Selection",
    # Ablation scenarios
    "ablation_no_he4":                     "Ablation (No $^4$He)",
    "ablation_raw_c14":                    "Ablation (Raw $^{14}$C)",
    "ablation_raw_gases":                  "Ablation (Raw Gases)",
    "oldwater_c14_ensemble":               "Old-Water $^{14}$C Ensemble",
    "oldwater_ensemble_he4_uncertainty":   "Old-Water Ensemble + $^4$He Uncertainty",
    "oldwater_he4_uncertainty":            "Old-Water $^4$He Uncertainty",
    "tracer_old_only":                     "Old Tracers Only",
    "tracer_young_only":                   "Young Tracers Only",
    # ── Age classes ─────────────────────────────────────────────────────
    "modern_le_50":                        "Modern \u226450 yr",
    "intermediate_50_1000":               "Intermediate 50\u20131,000 yr",
    "old_1000_30000":                      "Old 1,000\u201330,000 yr",
    "very_old_gt_30000":                   "Very old >30,000 yr",
    # ── Study Unit short codes ──────────────────────────────────────────
    "BNRF": "Basin & Range Basin-Fill",
    "BNRC": "Basin & Range Carbonate",
    "BISC": "Biscayne Aquifer",
    "CACB": "California Coastal Basin",
    "CMOR": "Cambrian-Ordovician",
    "CVAL": "Central Valley",
    "CLOW": "Coastal Lowlands",
    "COPL": "Colorado Plateaus",
    "CLPT": "Columbia Plateau Basalt",
    "EDTR": "Edwards-Trinity",
    "FLOR": "Floridan",
    "GLAC": "Glacial Aquifer",
    "HPAQ": "High Plains Aquifer",
    "METX": "Mississippi Embayment",
    "NACP": "Northern Atlantic Coastal Plain",
    "OZRK": "Ozark Plateaus",
    "PDBR": "Piedmont & Blue Ridge",
    "RIOG": "Rio Grande Aquifer",
    "SECP": "Southeastern Coastal Plain",
    "VRPD": "Valley & Ridge / Piedmont",
    # Study Unit full-name aliases (as they appear in the study_unit column)
    "Basin and Range basin-fill aquifers (BNRF)":                                  "Basin & Range Basin-Fill",
    "Basin and Range carbonate-rock aquifers (BNRC)":                              "Basin & Range Carbonate",
    "Biscayne aquifer (BISC)":                                                     "Biscayne Aquifer",
    "California Coastal Basin aquifers (CACB)":                                    "California Coastal Basin",
    "Cambrian-Ordovician aquifer system (CMOR)":                                   "Cambrian-Ordovician",
    "Central Valley aquifer system (CVAL)":                                        "Central Valley",
    "Coastal lowlands aquifer system (CLOW)":                                      "Coastal Lowlands",
    "Colorado Plateaus aquifers (COPL)":                                           "Colorado Plateaus",
    "Columbia Plateau basaltic-rock aquifers (CLPT)":                             "Columbia Plateau Basalt",
    "Edwards-Trinity aquifer system (EDTR)":                                       "Edwards-Trinity",
    "Floridan aquifer system (FLOR)":                                              "Floridan",
    "Glacial aquifer system (GLAC)":                                               "Glacial Aquifer",
    "High Plains aquifer (HPAQ)":                                                  "High Plains Aquifer",
    "Mississippi embayment-Texas coastal uplands aquifer system (METX)":          "Mississippi Embayment",
    "Northern Atlantic Coastal Plain aquifer system (NACP)":                      "N. Atlantic Coastal Plain",
    "Ozark Plateaus aquifer system (OZRK)":                                        "Ozark Plateaus",
    "Piedmont and Blue Ridge crystalline-rock aquifers (PDBR)":                   "Piedmont & Blue Ridge",
    "Rio Grande aquifer system (RIOG)":                                            "Rio Grande Aquifer",
    "Southeastern Coastal Plain aquifer system (SECP)":                           "Southeastern Coastal Plain",
    "Valley and Ridge and Piedmont and Blue Ridge carbonate-rock aquifers (VRPD)": "Valley & Ridge / Piedmont",
    # ── Graph families ──────────────────────────────────────────────────
    "depth_constrained":               "Depth-Constrained Graph",
    "hydraulic_proxy":                 "Hydraulic Proxy Graph",
    "hydraulic_proxy_constrained":     "Constrained Hydraulic Proxy",
    "randomized_control":              "Randomized Control Graph",
    "negative_control":                "Randomized Control Graph",
    "randomized_negative_control":     "Randomized Negative Control",
    "coordinate_global":               "Global Coordinate Graph",
    "study_unit_coordinate":           "Study-Unit Coordinate Graph",
    "parameter_smooth_aicc":           "Parameter Smoothing (AICc)",
    "wrong_direction_negative_control": "Wrong-Direction Control",
    # ── Prior strength ──────────────────────────────────────────────────
    "weak":   "Weak",
    "medium": "Medium",
    "strong": "Strong",
}


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _design_summary() -> pd.DataFrame:
    frames = []
    for candidate in (
        RESULT_DIR / "m3_tracerlpm_strict_parity_full_summary.csv",
        RESULT_DIR / "m3_tracerlpm_parity_hier_oldwater_full_summary.csv",
        RESULT_DIR / "m3_tracerlpm_parity_agefractions_full_summary.csv",
        RESULT_DIR / "m3_hydrosheaf_selection_corrected_full_summary.csv",
        RESULT_DIR / "m3_parity_reported_corrected_full_summary.csv",
        RESULT_DIR / "m3_design_matrix_summary.csv",
    ):
        df = _read_csv(candidate)
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    summary = pd.concat(frames, ignore_index=True)
    if "scenario_id" in summary.columns:
        summary = summary.drop_duplicates("scenario_id", keep="first")
    return summary


def _save(fig: plt.Figure, name: str) -> None:
    out = FIGURE_DIR / name
    tmp = out.with_suffix(out.suffix + ".tmp.png")
    fig.savefig(tmp, dpi=300, bbox_inches="tight")
    try:
        if out.exists():
            try:
                out.unlink()
            except Exception:
                pass
        tmp.replace(out)
    except Exception:
        time.sleep(1.0)
        try:
            tmp.replace(out)
        except Exception:
            print(f"Warning: Could not overwrite {name} due to lock. Saving to {name}.new.png instead.")
            alt = out.with_name(out.stem + ".new" + out.suffix)
            try:
                tmp.replace(alt)
            except Exception:
                pass
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except Exception:
                pass
        plt.close(fig)


def _primary_results(scenario_id: str | None = None) -> pd.DataFrame:
    """Load primary pointwise results, preferring strict parity scenarios."""
    candidates = (
        RESULT_DIR / "m3_tracerlpm_strict_parity_full.csv",
        RESULT_DIR / "m3_tracerlpm_parity_hier_oldwater_full.csv",
        RESULT_DIR / "m3_tracerlpm_parity_agefractions_full.csv",
        RESULT_DIR / "m3_design_matrix_results.csv",
        RESULT_DIR / "m3_phase4_screened_full_results.csv",
        RESULT_DIR / "m3_phase4_younggas_full_results.csv",
        RESULT_DIR / "m3_phase4_younggas_results.csv",
        RESULT_DIR / "m3_usgs_benchmark_results.csv",
    )
    for candidate in candidates:
        df = _read_csv(candidate)
        if df.empty:
            continue
        if scenario_id and "scenario_id" in df.columns and scenario_id in set(df["scenario_id"].dropna()):
            return df[df["scenario_id"] == scenario_id].copy()
        for preferred in ("tracerlpm_strict_parity", "tracerlpm_parity_hier_oldwater", "screened_dgm_gases", "parity_reported_corrected"):
            if "scenario_id" in df.columns and preferred in set(df["scenario_id"].dropna()):
                return df[df["scenario_id"] == preferred].copy()
        return df
    return pd.DataFrame()


def _intelligent_ticks(ax, axis="both", nbins=5, log=False):
    """Apply MaxNLocator for intelligent tick reduction on linear axes."""
    if log:
        return  # log axes auto-manage ticks
    loc = mticker.MaxNLocator(nbins=nbins, prune="both")
    if axis in ("x", "both"):
        ax.xaxis.set_major_locator(loc)
    if axis in ("y", "both"):
        ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=nbins, prune="both"))


# ---------------------------------------------------------------------------
# Figure 1: Atmospheric tracer input histories
# ---------------------------------------------------------------------------
def _benchmark_tritium_history():
    """Return the tritium forcing the benchmark actually used.

    Panel (a) previously plotted ``build_default_tritium_input()`` -- a stylised
    synthetic curve -- while the caption claimed GNIP/WISER provenance and the
    benchmark itself was forced with nearest-station WISER histories.  The panel
    therefore showed readers a curve that was neither GNIP data nor the function
    used in the analysis.  Plot a representative benchmark history instead, built
    through the same code path as the run (nearest WISER station, continued to
    the present so that post-record recharge years are not clamped to the
    station's last observed value).
    """
    from hydrosheaf.nuclear.tracer_inputs import (
        SiteInputContext,
        build_site_tracer_histories,
    )

    # Representative continental US public-supply site (mid-latitude, interior).
    ctx = SiteInputContext("representative", 2010.0, 40.0, -95.0)
    return build_site_tracer_histories(ctx)["3H"]


def plot_fig1_atmospheric_histories() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    histories = [
        ("a  Tritium ($^3$H)", _benchmark_tritium_history(), "Tritium (TU)", "#c2410c"),
        ("b  Sulfur Hexafluoride (SF$_6$)", build_atmospheric_tracer_input("SF6"), "SF6 (pptv)", "#1d4ed8"),
        ("c  CFC-12", build_atmospheric_tracer_input("CFC12"), "CFC-12 (pptv)", "#15803d"),
    ]
    flat_axes = axes.flatten()
    annotations = [
        # (year, label, text, xytext)
        [(1963, None, "Peak nuclear\ntesting ~1963", (1970, None))],   # 3H
        [(1980, None, "Industrial\nuptake\n~1980s", (1990, None))],     # SF6
        [(1989, None, "Montreal\nProtocol\n1989", (1995, None))],       # CFC-12
    ]
    for idx, (label, hist, unit, color) in enumerate(histories):
        ax = flat_axes[idx]
        ax.plot(hist.years, hist.values, color=color, lw=2.2)
        ax.fill_between(hist.years, hist.values, color=color, alpha=0.12)
        ax.set_title(label, loc="left", fontweight="bold", fontsize=13)
        ax.set_xlabel("Recharge year")
        ax.set_ylabel(unit)
        ax.grid(True, ls=":", alpha=0.4)
        if "$^3$H" in label or "3H" in label:
            ax.set_yscale("log")
            _intelligent_ticks(ax, axis="x", nbins=5, log=False)
        else:
            _intelligent_ticks(ax, axis="both", nbins=5)

        # Annotations
        for (yr, _, ann_text, (ann_yr, _)) in annotations[idx]:
            yvals = np.array(hist.values)
            xvals = np.array(hist.years)
            # Station histories are irregularly sampled, so match the nearest
            # available year rather than requiring an exact hit (an exact-match
            # test silently dropped the annotation on real WISER records).
            if xvals.size and np.min(np.abs(xvals - yr)) <= 2.0:
                yi = float(yvals[np.argmin(np.abs(xvals - yr))])
                ax.annotate(
                    ann_text,
                    xy=(yr, yi),
                    xytext=(ann_yr, yi * (3.0 if "$^3$H" in label or "3H" in label else 1.4)),
                    fontsize=9,
                    color="#374151",
                    arrowprops=dict(arrowstyle="->", color="#374151", lw=0.8),
                    ha="left",
                )

    # Hide the 4th panel
    flat_axes[3].set_visible(False)
    plt.tight_layout()
    _save(fig, "Manuscript_Fig1_Atmospheric_Histories.png")


# ---------------------------------------------------------------------------
# Figure 2: Design-matrix scenario performance
# ---------------------------------------------------------------------------
def plot_fig2_design_matrix_performance() -> None:
    summary = _design_summary()
    if summary.empty:
        return
    summary = summary.sort_values("median_abs_log10_error")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.2, 1.0]})
    y = np.arange(len(summary))
    cmap = plt.get_cmap("tab10")

    # Panel a: colored bars
    for idx, (_, row) in enumerate(summary.iterrows()):
        axes[0].barh(y[idx], row["median_abs_log10_error"], color=cmap(idx), height=0.6, edgecolor="white", lw=0.5)

    # Annotate best-performing scenario
    best_idx = int(np.argmin(summary["median_abs_log10_error"].values))
    best_val = summary["median_abs_log10_error"].values[best_idx]
    axes[0].annotate(
        f"Best: {best_val:.3f}",
        xy=(best_val, best_idx),
        xytext=(best_val + 0.02, best_idx - 0.5),
        fontsize=9,
        color="#0f766e",
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#0f766e", lw=0.8),
    )

    axes[0].set_yticks(y)
    scenario_labels = [LABEL_MAP.get(str(sid), str(sid).replace("_", " ")) for sid in summary["scenario_id"]]
    axes[0].set_yticklabels(scenario_labels, fontsize=11)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Median absolute log₁₀ error")
    axes[0].set_title("a  Scenario accuracy", loc="left", fontweight="bold")
    _intelligent_ticks(axes[0], axis="x", nbins=5)

    # Panel b: color-matched scatter; size = sample count
    scale = 1.25
    bubble_alpha = 0.45
    bubble_edge_width = 0.3
    for idx, (_, row) in enumerate(summary.iterrows()):
        axes[1].scatter(
            row["within_factor_2"],
            row["within_factor_10"],
            s=np.maximum(row["metric_rows"], 1) * scale,
            color=cmap(idx),
            alpha=bubble_alpha,
            edgecolor="black",
            linewidth=bubble_edge_width,
        )

    # Annotate the "ideal" corner
    axes[1].annotate(
        "Ideal\ncorner",
        xy=(1.0, 1.0),
        xytext=(0.72, 0.80),
        fontsize=9,
        color="#475569",
        arrowprops=dict(arrowstyle="->", color="#475569", lw=0.8),
    )

    axes[1].set_xlabel("Within factor 2 (fraction)")
    axes[1].set_ylabel("Within factor 10 (fraction)")
    axes[1].set_xlim(-0.05, 1.10)
    axes[1].set_ylim(-0.05, 1.10)
    axes[1].set_title("b  Practical accuracy  (bubble $\\propto$ N wells)", loc="left", fontweight="bold", fontsize=11)
    axes[1].grid(True, ls=":", alpha=0.4)
    _intelligent_ticks(axes[1], axis="both", nbins=5)

    # Size legend
    size_handles = []
    for n in [50, 200, 500]:
        size_handles.append(
            plt.scatter([], [], s=n * scale, color="#64748b", alpha=bubble_alpha,
                        edgecolor="black", linewidth=bubble_edge_width, label=f"N = {n}")
        )
    axes[1].legend(
        handles=size_handles,
        title="Sample count",
        title_fontsize=10,
        frameon=True,
        facecolor="white",
        edgecolor="#cbd5e1",
        loc="lower right",
        fontsize=10,
        scatterpoints=1,
    )
    plt.tight_layout()
    _save(fig, "Manuscript_Fig2_Design_Matrix_Performance.png")


def _mode_label_from_df(df: pd.DataFrame) -> str:
    """Return a human-readable mode label."""
    if "scenario_id" in df.columns and not df.empty:
        scenarios = set(df["scenario_id"].dropna().unique())
        if len(scenarios) == 1:
            sid = next(iter(scenarios))
            return LABEL_MAP.get(sid, sid)
    if "model_strategy" in df.columns and not df.empty:
        strategies = set(df["model_strategy"].dropna().unique())
        if len(strategies) == 1:
            strat = next(iter(strategies))
            return LABEL_MAP.get(strat, strat)
    return "benchmark"


# ---------------------------------------------------------------------------
# Figure 3: USGS benchmark parity plot
# ---------------------------------------------------------------------------
def plot_fig3_parity_plot() -> None:
    df = _primary_results().dropna(subset=["ref_age", "est_age_multi"])
    df = df[df["ref_age"] > 0].copy()
    if df.empty:
        return
    mode_label = _mode_label_from_df(df)
    df["plot_ref_age"] = np.maximum(df["ref_age"], 0.01)
    df["plot_est_age"] = np.maximum(df["est_age_multi"], 0.01)
    fig, ax = plt.subplots(figsize=(6.5, 6.2))
    classes = sorted(df["age_class"].dropna().unique())
    cmap = plt.get_cmap("tab10")
    for idx, age_class in enumerate(classes):
        subset = df[df["age_class"] == age_class]
        clean_ac = LABEL_MAP.get(age_class, str(age_class).replace("_", " ").title())
        ax.scatter(
            subset["plot_ref_age"],
            subset["plot_est_age"],
            s=35,
            alpha=0.72,
            color=cmap(idx),
            edgecolor="white",
            linewidth=0.4,
            label=clean_ac,
        )
    lims = [0.01, max(df["plot_ref_age"].max(), df["plot_est_age"].max(), 1e5) * 1.2]
    ax.plot(lims, lims, color="black", ls="--", lw=1.4, label="1:1 line")
    ax.fill_between(lims, np.array(lims) / 2.0, np.array(lims) * 2.0, color="#94a3b8", alpha=0.18, label="×2 envelope")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("USGS reference age (yr)")
    ax.set_ylabel("Hydrosheaf-estimated age (yr)")

    metric = df["log10_error"].dropna()
    within_raw = df.get("within_factor_2", pd.Series(np.nan, index=df.index))
    within_numeric = pd.to_numeric(within_raw, errors="coerce")
    within_bool = within_raw.astype(str).str.lower().map({"true": 1.0, "false": 0.0})
    within_factor_2 = within_numeric.where(within_numeric.notna(), within_bool).fillna(0.0).mean()

    # Stats box
    ax.text(
        0.04, 0.96,
        f"N = {len(df)}\nMedian |log₁₀ error| = {metric.median():.3f}\nWithin ×2 = {within_factor_2:.0%}\nScenario: {mode_label}",
        transform=ax.transAxes, va="top", ha="left", fontsize=10,
        bbox={"facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.9, "boxstyle": "round,pad=0.4"},
    )

    # Annotation for ×2 band
    ax.annotate(
        "×2 uncertainty\nenvelope",
        xy=(2, 4),
        xytext=(20, 0.15),
        fontsize=9,
        color="#475569",
        arrowprops=dict(arrowstyle="->", color="#475569", lw=0.8),
    )

    ax.legend(frameon=False, loc="lower right", fontsize=10)
    ax.grid(True, ls=":", alpha=0.4)
    _save(fig, "Manuscript_Fig3_USGS_Benchmark_Parity.png")


# ---------------------------------------------------------------------------
# Figure 4: Real-USGS graph-prior benchmark
# ---------------------------------------------------------------------------
def plot_fig4_real_graph_benchmark() -> None:
    graph = _read_csv(RESULT_DIR / "m3_real_usgs_graph_benchmark.csv")
    if graph.empty:
        return
    graph = graph[graph["prior_strength"] != "none"].copy()
    graph = graph.sort_values(["control_type", "delta_rmse_graph_minus_single"])

    clean_labels = []
    for _, row in graph.iterrows():
        fam = LABEL_MAP.get(row["graph_family"], row["graph_family"].replace("_", " ").title())
        strength = row["prior_strength"].title()
        clean_labels.append(f"{fam}\n({strength})")

    colors = np.where(graph["control_type"].eq("negative_control"), "#dc2626", "#0f766e")
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.axhline(0.0, color="black", lw=1.2, ls="-", zorder=3)
    ax.bar(np.arange(len(graph)), graph["delta_rmse_graph_minus_single"],
           color=colors, alpha=0.85, edgecolor="white", lw=0.5, zorder=2)
    ax.set_xticks(np.arange(len(graph)))
    ax.set_xticklabels(clean_labels, rotation=40, ha="right", fontsize=10)
    ax.set_ylabel("ΔRMSE, graph − single-well (log₁₀ yr)")
    ax.set_title("Real-USGS graph-prior benchmark", loc="left", fontweight="bold")
    ax.grid(True, ls=":", alpha=0.4, zorder=1)
    _intelligent_ticks(ax, axis="y", nbins=5)

    # Annotations placed INSIDE the axes using axes-fraction coordinates
    # so they never overlap with the rotated x-tick labels below the plot
    ax.text(
        0.02, 0.10,
        "Lower ΔRMSE = less degradation / closer to baseline",
        transform=ax.transAxes,
        fontsize=9, color="#0f766e", fontstyle="italic",
        va="bottom", ha="left",
    )
    ax.text(
        0.98, 0.92,
        "Graph degrades prediction →",
        transform=ax.transAxes,
        fontsize=9, color="#dc2626", fontstyle="italic",
        va="top", ha="right",
    )

    # Legend patches
    from matplotlib.patches import Patch
    legend_patches = [
        Patch(facecolor="#0f766e", alpha=0.85, label="Physical graph prior"),
        Patch(facecolor="#dc2626", alpha=0.85, label="Randomized control (negative)"),
    ]
    ax.legend(handles=legend_patches, frameon=True, facecolor="white",
              edgecolor="#cbd5e1", loc="upper left", fontsize=10)

    plt.tight_layout()
    _save(fig, "Manuscript_Fig4_Real_USGS_Graph_Benchmark.png")




# ---------------------------------------------------------------------------
# Figure 5: Age-class diagnostics
# ---------------------------------------------------------------------------
def plot_fig5_age_class_diagnostics() -> None:
    df = _primary_results().dropna(subset=["log10_error", "age_class"])
    if df.empty:
        return
    mode_label = _mode_label_from_df(df)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.0))
    age_classes = sorted(df["age_class"].unique())
    clean_ticks = [LABEL_MAP.get(ac, str(ac).replace("_", " ").title()) for ac in age_classes]
    data = [df.loc[df["age_class"] == age_class, "log10_error"].dropna() for age_class in age_classes]

    bp = axes[0].boxplot(data, tick_labels=clean_ticks, showfliers=False, patch_artist=True,
                         medianprops=dict(color="black", lw=1.5))
    cmap = plt.get_cmap("tab10")
    for patch, idx in zip(bp["boxes"], range(len(age_classes))):
        patch.set_facecolor(cmap(idx))
        patch.set_alpha(0.75)

    # Annotate worst age class
    medians = [d.median() for d in data]
    worst_idx = int(np.argmax(medians))
    axes[0].annotate(
        f"Highest median error\n({medians[worst_idx]:.2f})",
        xy=(worst_idx + 1, medians[worst_idx] + 0.05),
        xytext=(worst_idx + 1.45, medians[worst_idx] + 0.20),
        fontsize=9, color="#b45309",
        arrowprops=dict(arrowstyle="->", color="#b45309", lw=0.8),
    )

    axes[0].set_ylabel("Absolute log₁₀ error")
    axes[0].set_title("a  Absolute log₁₀ age error by reference-age class", loc="left", fontweight="bold")
    axes[0].tick_params(axis="x", rotation=30, labelsize=10)
    _intelligent_ticks(axes[0], axis="y", nbins=5)

    diag = df.dropna(subset=["young_gas_proxy_coherence", "log10_error"])
    sc = axes[1].scatter(
        diag["young_gas_proxy_coherence"],
        diag["log10_error"],
        c=np.log10(np.maximum(diag["ref_age"], 0.1)),
        cmap="viridis",
        s=35,
        alpha=0.72,
        edgecolor="white",
        linewidth=0.3,
    )
    cbar = plt.colorbar(sc, ax=axes[1], pad=0.02)
    cbar.set_label("log₁₀(reference age, yr)", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # Annotate high-discordance region
    axes[1].annotate(
        "High discordance\n(poor gas coherence)",
        xy=(diag["young_gas_proxy_coherence"].quantile(0.9), diag["log10_error"].median()),
        xytext=(diag["young_gas_proxy_coherence"].quantile(0.6), diag["log10_error"].quantile(0.75)),
        fontsize=9, color="#374151",
        arrowprops=dict(arrowstyle="->", color="#374151", lw=0.8),
    )

    axes[1].set_xlabel("Young-tracer discordance (log₁₀-age SD)")
    axes[1].set_ylabel("Absolute log₁₀ error")
    axes[1].set_title("b  Discordance diagnostic", loc="left", fontweight="bold")
    axes[1].grid(True, ls=":", alpha=0.4)
    _intelligent_ticks(axes[1], axis="both", nbins=5)
    fig.suptitle(f"Age-Class Diagnostics — {mode_label}", fontweight="bold", fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, "Manuscript_Fig5_Age_Class_Diagnostics.png")


# ---------------------------------------------------------------------------
# Figure 6: Cross-validation results
# ---------------------------------------------------------------------------
def plot_fig6_cross_validation_results() -> None:
    df_3h = _read_csv(RESULT_DIR / "m3_cv_benchmark_3H.csv")
    df_sf6 = _read_csv(RESULT_DIR / "m3_cv_benchmark_SF6.csv")
    df_c14 = _read_csv(RESULT_DIR / "m3_cv_benchmark_14C.csv")

    if df_3h.empty or df_sf6.empty or df_c14.empty:
        print("Warning: Cross-validation CSV files missing. Cannot plot Fig 6.")
        return

    def get_rmse(df, err_col):
        errs = df[err_col].dropna()
        return np.sqrt(np.mean(errs**2)) if len(errs) > 0 else 0.0

    rmse_data = {
        "3H":  [get_rmse(df_3h, c)  for c in ("err_single_abs", "err_phys_hyd_abs", "err_phys_dep_abs", "err_rand_abs")],
        "SF6": [get_rmse(df_sf6, c) for c in ("err_single_abs", "err_phys_hyd_abs", "err_phys_dep_abs", "err_rand_abs")],
        "14C": [get_rmse(df_c14, c) for c in ("err_single_abs", "err_phys_hyd_abs", "err_phys_dep_abs", "err_rand_abs")],
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    flat_axes = axes.flatten()
    labels = ["Single-Well\nBaseline", "Hydraulic-Proxy\nGraph", "Depth-Constrained\nGraph", "Randomized\nControl"]
    x = np.arange(len(labels))
    width = 0.35

    # Panel a: young-tracer error change relative to the single-well baseline.
    # Previously ³H (TU) and SF₆ (pptv) shared one panel on two y-axes with
    # different baselines, so the bar heights were not comparable -- the figure
    # note conceded as much.  Expressing both as percentage change from the
    # single-well baseline puts them on one honest, dimensionless axis.
    ax0 = flat_axes[0]

    def _pct_change(vals):
        base = vals[0]
        return [100.0 * (v - base) / base if base > 0 else np.nan for v in vals]

    pct_3h = _pct_change(rmse_data["3H"])
    pct_sf6 = _pct_change(rmse_data["SF6"])
    ax0.bar(x - width / 2, pct_3h, width, label="³H", color="#c2410c",
            alpha=0.85, edgecolor="white", lw=0.5)
    ax0.bar(x + width / 2, pct_sf6, width, label="SF₆", color="#1d4ed8",
            alpha=0.85, edgecolor="white", lw=0.5)
    ax0.axhline(0, color="#374151", lw=1.0)

    # Annotate BOTH the depth-constrained result and the randomised negative
    # control.  The randomised control also lowers ³H RMSE, which is the
    # decisive context for interpreting the depth-constrained gain; labelling
    # only the supporting bar would misrepresent the comparison.
    for idx, colour in ((2, "#c2410c"), (3, "#b91c1c")):
        val = pct_3h[idx]
        if not np.isfinite(val):
            continue
        ax0.annotate(
            f"{val:+.1f}%",
            xy=(x[idx] - width / 2, val),
            xytext=(x[idx] - width / 2, val - 4.0),
            fontsize=9.5, color=colour, fontweight="bold", ha="center", va="top",
        )
    # Describe what the controls actually did rather than asserting a fixed
    # conclusion: the sign of the randomised-control effect is the whole point of
    # the negative-control design and must be read from the data.
    rand_val, cand_val = pct_3h[3], pct_3h[2]
    if np.isfinite(rand_val) and np.isfinite(cand_val):
        if rand_val < 0 and cand_val < 0:
            note = "the randomised control also lowers ³H RMSE"
        elif rand_val > 0 and cand_val > 0:
            note = "candidate and randomised graphs both degrade ³H;\nthe control degrades far more"
        else:
            note = "candidate and randomised controls differ in sign"
        ax0.text(
            0.5, -0.34, note, transform=ax0.transAxes,
            fontsize=9, color="#b91c1c", ha="center", va="top",
        )

    ax0.set_ylabel("Change in withheld-tracer RMSE\nvs single-well baseline (%)", fontsize=11)
    ax0.set_xticks(x)
    ax0.set_xticklabels(labels, rotation=20, ha="right", fontsize=10)
    ax0.set_title("(a) Young-tracer cross-validation error", loc="left", fontweight="bold")
    ax0.legend(frameon=False, fontsize=10, loc="upper left")
    ax0.grid(True, ls=":", alpha=0.3)
    ax0.spines["top"].set_visible(False)
    ax0.spines["right"].set_visible(False)
    _intelligent_ticks(ax0, axis="y", nbins=5)

    # Panel b: ¹⁴C RMSE
    ax1 = flat_axes[1]
    bar_colors = ["#475569", "#0284c7", "#0f766e", "#dc2626"]
    bars = ax1.bar(x, rmse_data["14C"], color=bar_colors, alpha=0.85, edgecolor="white", lw=0.5)
    ax1.set_ylabel("¹⁴C RMSE (pmc)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=20, ha="right", fontsize=10)
    ax1.set_title("(b) Old-tracer cross-validation error", loc="left", fontweight="bold")
    ax1.grid(True, ls=":", alpha=0.3)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    _intelligent_ticks(ax1, axis="y", nbins=4)

    # Annotate the randomized control bar to highlight degradation
    rand_c14 = rmse_data["14C"][3]
    base_c14 = rmse_data["14C"][0]
    pct_degrade = 100 * (rand_c14 - base_c14) / base_c14 if base_c14 > 0 else 0
    ax1.annotate(
        f"+{pct_degrade:.1f}%\n(random noise)",
        xy=(x[3], rand_c14),
        xytext=(x[3] - 0.5, rand_c14 * 0.85),
        fontsize=9, color="#dc2626", fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#dc2626", lw=0.8),
    )

    # Panel c: ³H Observed vs Predicted scatter
    ax2 = flat_axes[2]
    valid = df_3h[df_3h["true_val"] > 0.1].dropna(subset=["pred_single", "pred_phys_dep"])
    ax2.scatter(valid["true_val"], valid["pred_single"], s=12, color="#dc2626", alpha=0.55, label="Single-well baseline", edgecolor="none")
    ax2.scatter(valid["true_val"], valid["pred_phys_dep"], s=12, color="#0f766e", alpha=0.65, label="Depth-constrained graph", edgecolor="none")

    lims = [0.1, max(valid["true_val"].max(), valid["pred_single"].max(), valid["pred_phys_dep"].max()) * 1.15]
    ax2.plot(lims, lims, color="black", ls="--", lw=1.2, label="1:1 line")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_xlabel("Observed ³H (TU)")
    ax2.set_ylabel("Predicted ³H (TU)")
    ax2.set_title("(c) ³H observed–predicted cross-validation", loc="left", fontweight="bold")
    ax2.legend(frameon=False, loc="lower right", fontsize=10)
    ax2.grid(True, ls=":", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Annotate 1:1 line
    ax2.annotate(
        "1:1 perfect\nprediction",
        xy=(6, 6),
        xytext=(0.6, 14),
        fontsize=9, color="#374151",
        arrowprops=dict(arrowstyle="->", color="#374151", lw=0.8),
    )

    # Clarification of low-end cluster
    ax2.annotate(
        "Plotting floor &\ndetection limit (0.1 TU)",
        xy=(0.115, 0.25),
        xytext=(0.45, 0.115),
        fontsize=8.5,
        color="#475569",
        arrowprops=dict(arrowstyle="->", color="#475569", lw=0.8),
    )

    # Panel d: is the ³H gain distributed, or tail-driven?
    # RMSE is dominated by a small number of degenerate young-water fits.  Showing
    # RMSE beside the median absolute error makes clear whether a graph prior
    # improves typical predictions or only suppresses outliers.  Reporting RMSE
    # alone would let a tail-only effect read as a general accuracy gain.
    ax3 = flat_axes[3]
    err_cols = ("err_single_abs", "err_phys_hyd_abs", "err_phys_dep_abs", "err_rand_abs")
    med_3h = [float(df_3h[c].dropna().median()) if df_3h[c].notna().any() else np.nan
              for c in err_cols]
    pct_med = _pct_change(med_3h)
    ax3.bar(x - width / 2, pct_3h, width, label="RMSE", color="#c2410c",
            alpha=0.85, edgecolor="white", lw=0.5)
    ax3.bar(x + width / 2, pct_med, width, label="Median absolute error",
            color="#a16207", alpha=0.85, edgecolor="white", lw=0.5)
    ax3.axhline(0, color="#374151", lw=1.0)
    ax3.set_ylabel("Change vs single-well baseline (%)", fontsize=11)
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=20, ha="right", fontsize=10)
    ax3.set_title("(d) ³H error structure: RMSE vs typical error", loc="left",
                  fontweight="bold")
    ax3.legend(frameon=False, fontsize=10, loc="upper left")
    ax3.grid(True, ls=":", alpha=0.3)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    _intelligent_ticks(ax3, axis="y", nbins=5)

    fig.text(
        0.1, 0.02,
        "Panels (a) and (d) report change relative to the single-well baseline, so ³H and SF₆ are directly comparable. "
        "Panel (d) separates RMSE (outlier-sensitive) from the median absolute error (typical case).",
        fontsize=10, style="italic",
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    _save(fig, "Manuscript_Fig6_Cross_Validation_Results.png")


# ---------------------------------------------------------------------------
# Figure 7: Residual Diagnostics
# ---------------------------------------------------------------------------
def plot_fig7_residual_diagnostics() -> None:
    """Plot boxplots of signed log10 residuals grouped by Age Class, Tracer Availability, and Old-Water Status.
    Layout: 2x2 grid with the bottom-right panel left empty.
    """
    df = _read_csv(RESULT_DIR / "m3_tracerlpm_strict_parity_full.csv")
    if df.empty:
        return

    # Filter for valid values
    df = df[df["ref_age"].notna() & df["est_age_multi"].notna()].copy()
    df = df[(df["ref_age"] > 0) & (df["est_age_multi"] > 0)]
    if df.empty:
        return

    # Calculate signed residual (log10 est - log10 ref)
    df["residual"] = np.log10(df["est_age_multi"]) - np.log10(df["ref_age"])

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    cmap = plt.get_cmap("tab10")

    ax_a = axes[0, 0]
    ax_b = axes[0, 1]
    ax_c = axes[1, 0]
    ax_empty = axes[1, 1]

    # Hide the empty bottom-right panel
    ax_empty.set_visible(False)

    # --- Panel (a) Age Class ---
    ac_order = ["modern_le_50", "intermediate_50_1000", "old_1000_30000", "very_old_gt_30000"]
    ac_present = [ac for ac in ac_order if ac in df["age_class"].values]
    ac_labels = [f"{LABEL_MAP.get(ac, ac.replace('_', ' ').title())}\n(N={len(df[df['age_class'] == ac])})" for ac in ac_present]
    data_ac = [df[df["age_class"] == ac]["residual"].dropna() for ac in ac_present]

    bp_ac = ax_a.boxplot(data_ac, tick_labels=ac_labels, showfliers=False, patch_artist=True,
                         medianprops=dict(color="black", lw=1.5))
    for patch, idx in zip(bp_ac["boxes"], range(len(ac_present))):
        patch.set_facecolor(cmap(idx))
        patch.set_alpha(0.75)

    ax_a.axhline(0.0, color="black", ls="--", lw=1.2, zorder=0)
    ax_a.set_ylabel("Signed log₁₀ error\n(estimated age − reference age)", fontsize=11)
    ax_a.set_title("a  Signed Residuals by Age Class", loc="left", fontweight="bold", fontsize=12)
    ax_a.tick_params(axis="x", rotation=15, labelsize=9)
    ax_a.grid(True, ls=":", alpha=0.3)

    # --- Panel (b) Tracer Availability ---
    def get_tracer_class(t_str):
        if not isinstance(t_str, str):
            return "Unknown"
        tracers = [t.strip() for t in t_str.split(",")]
        has_young = any("3H" in t or "3He" in t or "SF6" in t or "CFC" in t for t in tracers)
        has_old = any("14C" in t or "4He" in t for t in tracers)
        if has_young and has_old:
            return "Mixed (Young & Old)"
        elif has_old:
            return "Old Only ($^{14}$C/$^4$He)"
        elif has_young:
            return "Young Only ($^3$H/Gases)"
        return "Unknown"

    df["tracer_class"] = df["LPM_TracersMod"].apply(get_tracer_class)
    tc_order = ["Young Only ($^3$H/Gases)", "Old Only ($^{14}$C/$^4$He)", "Mixed (Young & Old)"]
    tc_present = [tc for tc in tc_order if tc in df["tracer_class"].values]
    tc_labels = [f"{tc}\n(N={len(df[df['tracer_class'] == tc])})" for tc in tc_present]
    data_tc = [df[df["tracer_class"] == tc]["residual"].dropna() for tc in tc_present]

    bp_tc = ax_b.boxplot(data_tc, tick_labels=tc_labels, showfliers=False, patch_artist=True,
                         medianprops=dict(color="black", lw=1.5))
    for patch, idx in zip(bp_tc["boxes"], range(len(tc_present))):
        patch.set_facecolor(cmap(idx + len(ac_present)))
        patch.set_alpha(0.75)

    ax_b.axhline(0.0, color="black", ls="--", lw=1.2, zorder=0)
    ax_b.set_ylabel("Signed log₁₀ error\n(estimated age − reference age)", fontsize=11)
    ax_b.set_title("b  Signed Residuals by Tracer Availability", loc="left", fontweight="bold", fontsize=12)
    ax_b.tick_params(axis="x", rotation=15, labelsize=9)
    ax_b.grid(True, ls=":", alpha=0.3)

    # --- Panel (c) Old-Water Diagnostic Status ---
    status_order = ["none", "single_tracer", "agreement", "tension", "conflict"]
    status_present = [st for st in status_order if st in df["old_groundwater_status"].values]
    status_labels = [f"{st.replace('_', ' ').title()}\n(N={len(df[df['old_groundwater_status'] == st])})" for st in status_present]
    data_st = [df[df["old_groundwater_status"] == st]["residual"].dropna() for st in status_present]

    bp_st = ax_c.boxplot(data_st, tick_labels=status_labels, showfliers=False, patch_artist=True,
                         medianprops=dict(color="black", lw=1.5))
    for patch, idx in zip(bp_st["boxes"], range(len(status_present))):
        patch.set_facecolor(cmap(idx + len(ac_present) + len(tc_present)))
        patch.set_alpha(0.75)

    ax_c.axhline(0.0, color="black", ls="--", lw=1.2, zorder=0)
    ax_c.set_ylabel("Signed log₁₀ error\n(estimated age − reference age)", fontsize=11)
    ax_c.set_title("c  Signed Residuals by Old-Water Status", loc="left", fontweight="bold", fontsize=12)
    ax_c.tick_params(axis="x", rotation=15, labelsize=9)
    ax_c.grid(True, ls=":", alpha=0.3)

    plt.tight_layout()
    _save(fig, "Manuscript_Fig7_Residual_Diagnostics.png")


# ---------------------------------------------------------------------------
# Supplementary Figure 1: Gas correction audit
# ---------------------------------------------------------------------------
def plot_supp1_gas_audit() -> None:
    df = _read_csv(RESULT_DIR / "m3_gas_correction_audit.csv")
    if df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # ── Panel a: Dissolved Gas Correction Model Distribution by Age Class ──
    # NOTE: For the USGS benchmark wells, the gas correction (UA model) did not
    # alter measured tracer concentrations (any_gas_value_changed = False for all
    # 80 wells), so parity_log10_error == raw_gas_log10_error. Instead, we show
    # the distribution of correction framework types applied per age class, which
    # provides methodological context for the gas-correction step.
    ac_order = ["modern_le_50", "intermediate_50_1000", "old_1000_30000", "very_old_gt_30000"]
    ac_present = [ac for ac in ac_order if ac in df["age_class"].values]
    ac_labels = [LABEL_MAP.get(ac, ac.replace("_", " ").title()) for ac in ac_present]

    dgm_counts = (
        df.groupby(["age_class", "dgm_name"])
        .size()
        .reset_index(name="count")
    )
    dgm_names = sorted(dgm_counts["dgm_name"].dropna().unique())
    cmap_dgm = plt.get_cmap("tab10")
    x = np.arange(len(ac_present))
    width = 0.8 / max(len(dgm_names), 1)
    for di, dgm in enumerate(dgm_names):
        heights = []
        for ac in ac_present:
            row = dgm_counts[(dgm_counts["age_class"] == ac) & (dgm_counts["dgm_name"] == dgm)]
            heights.append(int(row["count"].values[0]) if len(row) > 0 else 0)
        axes[0].bar(
            x + di * width - (len(dgm_names) - 1) * width / 2,
            heights, width * 0.9,
            color=cmap_dgm(di), alpha=0.85, edgecolor="white", lw=0.5,
            label=dgm if dgm else "Unknown",
        )

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(ac_labels, rotation=20, ha="right", fontsize=11)
    axes[0].set_ylabel("Number of wells")
    axes[0].set_title("a  Dissolved Gas Correction Model by Age Class", loc="left", fontweight="bold")
    axes[0].legend(frameon=True, facecolor="white", edgecolor="#cbd5e1",
                   loc="upper right", fontsize=9, title="Correction model", title_fontsize=9)
    axes[0].grid(True, ls=":", alpha=0.3)
    _intelligent_ticks(axes[0], axis="y", nbins=5)

    # Annotate: note corrections had no net tracer impact
    axes[0].text(
        0.02, 0.97,
        "UA corrections did not materially\nalter ³H–³He apparent ages\nfor these benchmark wells",
        transform=axes[0].transAxes, va="top", ha="left", fontsize=8.5,
        color="#475569",
        bbox={"facecolor": "#f8fafc", "edgecolor": "#cbd5e1", "alpha": 0.9, "boxstyle": "round,pad=0.3"},
    )

    # ── Panel b: ³H–³He apparent age scatter ────────────────────────────
    valid_ages = df.dropna(subset=["corrected_3h3he_apparent_age", "raw_3h3he_apparent_age"])
    if not valid_ages.empty:
        ax_max = max(valid_ages["raw_3h3he_apparent_age"].max(),
                     valid_ages["corrected_3h3he_apparent_age"].max(), 30) * 1.1
        cmap_ac = plt.get_cmap("tab10")
        for ci, ac in enumerate(ac_present):
            sub = valid_ages[valid_ages["age_class"] == ac]
            if sub.empty:
                continue
            axes[1].scatter(
                sub["raw_3h3he_apparent_age"],
                sub["corrected_3h3he_apparent_age"],
                s=25, color=cmap_ac(ci), alpha=0.75, edgecolor="white", lw=0.3,
                label=LABEL_MAP.get(ac, ac.replace("_", " ").title()),
            )
        axes[1].plot([0, ax_max], [0, ax_max], color="black", ls="--", lw=1.2, label="1:1 line")
        axes[1].set_xlim(0, ax_max)
        axes[1].set_ylim(0, ax_max)
        axes[1].legend(frameon=False, fontsize=9)
    else:
        axes[1].text(0.5, 0.5, "Insufficient ³H–³He data\nfor this benchmark",
                     transform=axes[1].transAxes, ha="center", va="center",
                     fontsize=11, color="#64748b")

    axes[1].set_xlabel("Uncorrected ³H–³He apparent age (yr)")
    axes[1].set_ylabel("Dissolved-gas corrected ³H–³He apparent age (yr)")
    axes[1].set_title("b  Corrected versus uncorrected ³H–³He apparent age", loc="left", fontweight="bold")
    axes[1].grid(True, ls=":", alpha=0.3)
    _intelligent_ticks(axes[1], axis="both", nbins=5)

    plt.tight_layout()
    _save(fig, "Suppl_Fig1_Gas_Correction_Audit.png")


# ---------------------------------------------------------------------------
# Supplementary Figure 2: Graph edge properties
# ---------------------------------------------------------------------------
def plot_supp2_graph_edges() -> None:
    df = _read_csv(RESULT_DIR / "m3_real_usgs_graph_edges.csv")
    if df.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    families = df["family"].unique()
    colors = {
        "depth_constrained":               "#0f766e",
        "hydraulic_proxy_constrained":     "#1d4ed8",
        "study_unit_coordinate":           "#0284c7",
        "coordinate_global":               "#475569",
        "wrong_direction_negative_control": "#b91c1c",
        "randomized_negative_control":     "#dc2626",
    }

    def is_control(fam: str) -> bool:
        return "control" in fam or "randomized" in fam

    # Panel a setup
    valid_dist = df["distance_km"].dropna()
    min_dist = max(valid_dist.min(), 0.1) if not valid_dist.empty else 0.1
    max_dist = valid_dist.max() if not valid_dist.empty else 100.0
    bins_a = np.logspace(np.log10(min_dist), np.log10(max_dist), 30)

    for fam in families:
        subset = df[df["family"] == fam]
        label  = LABEL_MAP.get(fam, fam.replace("_", " ").title())
        color  = colors.get(fam, "#64748b")
        if is_control(fam):
            axes[0].hist(subset["distance_km"], bins=bins_a, alpha=0.9, label=label,
                         color=color, density=True, histtype="step", lw=2, ls="--")
        else:
            axes[0].hist(subset["distance_km"], bins=bins_a, alpha=0.35, label=label,
                         color=color, density=True, histtype="stepfilled")

    axes[0].set_xscale("log")
    axes[0].set_xlabel("Spatial Edge Distance (km, log scale)")
    axes[0].set_ylabel("Density")
    axes[0].set_title("a  Edge Spatial Distance Distribution", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="upper right", fontsize=10)
    axes[0].grid(True, ls=":", alpha=0.3)
    _intelligent_ticks(axes[0], axis="both", nbins=5, log=True)

    for fam in families:
        subset = df[df["family"] == fam].dropna(subset=["depth_diff_m"])
        if subset.empty:
            continue
        label  = LABEL_MAP.get(fam, fam.replace("_", " ").title())
        color  = colors.get(fam, "#64748b")
        if is_control(fam):
            axes[1].hist(subset["depth_diff_m"], bins=30, alpha=0.9, label=label,
                         color=color, density=True, histtype="step", lw=2, ls="--")
        else:
            axes[1].hist(subset["depth_diff_m"], bins=30, alpha=0.35, label=label,
                         color=color, density=True, histtype="stepfilled")

    axes[1].axvline(0.0, color="black", ls="--", lw=1.2, zorder=3, label="Equal depth")
    axes[1].set_xlabel("Depth Difference between Nodes (m)\n(positive = downstream node deeper than upstream node)")
    axes[1].set_ylabel("Density")
    axes[1].set_title("b  Vertical (Depth) Gradient Distribution", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, loc="upper right", fontsize=10)
    axes[1].grid(True, ls=":", alpha=0.3)
    _intelligent_ticks(axes[1], axis="both", nbins=5)

    plt.tight_layout()
    _save(fig, "Suppl_Fig2_Graph_Edge_Properties.png")


# ---------------------------------------------------------------------------
# Supplementary Figure 3: CFC cross-validation
# ---------------------------------------------------------------------------
def plot_supp3_cfc_cv() -> None:
    df_cfc11 = _read_csv(RESULT_DIR / "m3_cv_benchmark_CFC11.csv")
    df_cfc12 = _read_csv(RESULT_DIR / "m3_cv_benchmark_CFC12.csv")
    if df_cfc11.empty or df_cfc12.empty:
        return

    def get_rmse(df, err_col):
        errs = df[err_col].dropna()
        return np.sqrt(np.mean(errs**2)) if len(errs) > 0 else 0.0

    rmse_cfc11 = [get_rmse(df_cfc11, c) for c in ("err_single_abs", "err_phys_hyd_abs", "err_phys_dep_abs", "err_rand_abs")]
    rmse_cfc12 = [get_rmse(df_cfc12, c) for c in ("err_single_abs", "err_phys_hyd_abs", "err_phys_dep_abs", "err_rand_abs")]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    labels = ["Single-Well\nBaseline", "Hydraulic-Proxy\nGraph", "Depth-Constrained\nGraph", "Randomized\nControl"]
    x = np.arange(len(labels))
    colors = ["#475569", "#1d4ed8", "#0f766e", "#dc2626"]

    # CFC-11 Panel
    axes[0].bar(x, rmse_cfc11, color=colors, alpha=0.85, edgecolor="white", lw=0.5)
    axes[0].set_ylabel("CFC-11 RMSE (pptv)")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, fontsize=10)
    axes[0].set_title("a  CFC-11 cross-validation error, N = 28", loc="left", fontweight="bold")
    axes[0].grid(True, ls=":", alpha=0.3)
    
    # Percentage-change annotations for CFC-11
    for i in range(1, 4):
        val = rmse_cfc11[i]
        pct = 100 * (val - rmse_cfc11[0]) / rmse_cfc11[0]
        sign = "+" if pct >= 0 else "−"
        label = f"{sign}{abs(pct):.1f}%"
        axes[0].text(x[i], val + 0.01 * max(rmse_cfc11), label, ha="center", va="bottom", fontsize=9, fontweight="bold", color="#1e293b")
    
    axes[0].set_ylim(0, max(rmse_cfc11) * 1.15)
    _intelligent_ticks(axes[0], axis="y", nbins=4)

    # CFC-12 Panel
    axes[1].bar(x, rmse_cfc12, color=colors, alpha=0.85, edgecolor="white", lw=0.5)
    axes[1].set_ylabel("CFC-12 RMSE (pptv)")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, fontsize=10)
    axes[1].set_title("b  CFC-12 cross-validation error, N = 16", loc="left", fontweight="bold")
    axes[1].grid(True, ls=":", alpha=0.3)

    # Percentage-change annotations for CFC-12
    for i in range(1, 4):
        val = rmse_cfc12[i]
        pct = 100 * (val - rmse_cfc12[0]) / rmse_cfc12[0]
        sign = "+" if pct >= 0 else "−"
        label = f"{sign}{abs(pct):.1f}%"
        axes[1].text(x[i], val + 0.01 * max(rmse_cfc12), label, ha="center", va="bottom", fontsize=9, fontweight="bold", color="#1e293b")

    axes[1].set_ylim(0, max(rmse_cfc12) * 1.15)
    _intelligent_ticks(axes[1], axis="y", nbins=4)

    plt.tight_layout()
    _save(fig, "Suppl_Fig3_CFC_CV_Performance.png")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> int:
    plot_fig1_atmospheric_histories()
    plot_fig2_design_matrix_performance()
    plot_fig3_parity_plot()
    plot_fig4_real_graph_benchmark()
    plot_fig5_age_class_diagnostics()
    plot_fig6_cross_validation_results()
    plot_fig7_residual_diagnostics()
    plot_supp1_gas_audit()
    plot_supp2_graph_edges()
    plot_supp3_cfc_cv()
    print(f"Wrote M3 manuscript figures to {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
