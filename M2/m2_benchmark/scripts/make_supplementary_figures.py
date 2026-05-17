from __future__ import annotations

from pathlib import Path
import uuid
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MaxNLocator

# --- Path Configuration ---
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DIR = BENCHMARK_ROOT / "external"
RESULT_DIR = BENCHMARK_ROOT / "results"
# Output directly to Manuscript folder
FIGURE_DIR = BENCHMARK_ROOT / "figures" / "Manuscript_Ready"

# --- Typography Standards (Matching Fig 5 style) ---
FONT_TITLE = 18
FONT_LABEL = 16
FONT_TICK = 16
FONT_LEGEND = 14
FONT_ANNOTATE = 11

COLOR_PRIMARY = "#3b82f6"  # Blue
COLOR_ACCENT = "#ef4444"   # Red
COLOR_SUCCESS = "#22c55e"  # Green
COLOR_WARNING = "#f59e0b"  # Gold
GRID_ALPHA = 0.1

def _save_supp(fig: plt.Figure, name: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    target = FIGURE_DIR / name
    temp = target.with_name(f"{target.stem}.{uuid.uuid4().hex}.tmp{target.suffix}")
    fig.savefig(temp, dpi=300, bbox_inches="tight", format=target.suffix.lstrip("."))
    temp.replace(target)
    plt.close(fig)

def plot_s1_age_parity() -> None:
    """Figure S1: Regional Residence-Time Validation (USGS Dataset)."""
    path = EXTERNAL_DIR / "usgs_age" / "results" / "usgs_age_validation.csv"
    if not path.exists():
        print("Skipping S1: USGS age data not found.")
        return

    df = pd.read_csv(path)
    df_clean = df.dropna(subset=["reference_mean_age_years", "hydrosheaf_age_years"]).copy()
    df_clean = df_clean[df_clean["reference_mean_age_years"] <= 50000]

    def get_category(row):
        tracers = str(row["supported_tracers"])
        if "14C" in tracers: return "14C Inclusive"
        if any(x in tracers for x in ["SF6", "CFC", "3H"]): return "Modern (3H/SF6/CFC)"
        return "Other"

    df_clean["Category"] = df_clean.apply(get_category, axis=1)

    log_ref = np.log10(np.maximum(df_clean["reference_mean_age_years"], 0.1))
    log_inf = np.log10(np.maximum(df_clean["hydrosheaf_age_years"], 0.1))

    r2_global = np.corrcoef(log_ref, log_inf)[0,1]**2
    mae_log = np.mean(np.abs(log_ref - log_inf))
    within_ci = (df_clean["inside_hydrosheaf_ci"] == True).mean() * 100

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(4, 4, left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_histx = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    colors = {"Modern (3H/SF6/CFC)": COLOR_PRIMARY, "14C Inclusive": COLOR_ACCENT, "Other": "#64748b"}
    for cat, color in colors.items():
        sub = df_clean[df_clean["Category"] == cat]
        ax_main.scatter(sub["reference_mean_age_years"], sub["hydrosheaf_age_years"],
                        c=color, alpha=0.5, s=30, label=cat, edgecolors='white', linewidth=0.3)

    lims = np.logspace(-1, 5, 100)
    ax_main.plot(lims, lims, color="#1f2937", ls="--", alpha=0.8, lw=1.5, label="1:1 Line")
    ax_main.fill_between(lims, lims/10, lims*10, color="#94a3b8", alpha=0.15, label="±1 Order of Mag")

    bins = np.logspace(-1, 5, 40)
    ax_histx.hist(df_clean["reference_mean_age_years"], bins=bins, color="#cbd5e1", edgecolor="white")
    ax_histy.hist(df_clean["hydrosheaf_age_years"], bins=bins, orientation='horizontal', color="#cbd5e1", edgecolor="white")

    ax_main.set_xscale("log")
    ax_main.set_yscale("log")
    ax_main.set_xlim(0.1, 100000)
    ax_main.set_ylim(0.1, 100000)
    ax_main.set_xlabel("Reference Tracer Age (Years)", fontsize=FONT_LABEL, fontweight="bold")
    ax_main.set_ylabel("Hydrosheaf Inferred Age (Years)", fontsize=FONT_LABEL, fontweight="bold")
    ax_main.tick_params(labelsize=FONT_TICK)

    stats_text = (f"Global $R^2 = {r2_global:.2f}$\n"
                  f"MAE = {mae_log:.2f} orders\n"
                  f"Within 95% CI = {within_ci:.1f}%")
    ax_main.text(0.05, 0.95, stats_text, transform=ax_main.transAxes,
                 verticalalignment='top', fontsize=FONT_ANNOTATE, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="#d1d5db", alpha=0.9))

    ax_main.legend(loc="lower right", fontsize=FONT_LEGEND, frameon=True)
    ax_main.grid(True, which="major", ls="-", alpha=GRID_ALPHA)

    ax_histx.axis("off")
    ax_histy.axis("off")
    ax_histx.set_title("Supplementary Figure S1: Regional Residence-Time Validation (USGS Dataset)",
                       fontsize=FONT_TITLE, fontweight="bold", pad=20)

    _save_supp(fig, "Manuscript_Supp_FigS1_Public_Age_Validation.png")

def plot_s2_geochemical_validation() -> None:
    """Figure S2: Geochemical Validation Suite (RMSE, NSE, and SI)."""
    path = BENCHMARK_ROOT / "results" / "phreeqc_forward_validation.csv"
    if not path.exists():
        print(f"Skipping S2: PHREEQC forward data not found at {path}")
        return

    df = pd.read_csv(path)
    # 2x2 Grid, one empty panel
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    ((ax1, ax2), (ax3, ax4)) = axes

    # (a) RMSE Histogram
    ax1.hist(df["rmse_mmolL"], bins=20, color=COLOR_PRIMARY, alpha=0.7, edgecolor="white")
    ax1.axvline(df["rmse_mmolL"].median(), color=COLOR_ACCENT, ls="--", lw=2,
                label=f"Median: {df['rmse_mmolL'].median():.3f}")
    ax1.set_title("(a) Geochemical Forward RMSE", fontsize=FONT_TITLE, fontweight="bold")
    ax1.set_xlabel("RMSE (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax1.set_ylabel("Frequency", fontsize=FONT_LABEL, fontweight="bold")
    ax1.tick_params(labelsize=FONT_TICK)
    ax1.legend(fontsize=FONT_LEGEND)

    # (b) Model Efficiency (NSE)
    df_clean = df.dropna(subset=['nse'])
    ax2.scatter(df_clean["rmse_mmolL"], df_clean["nse"], s=40, alpha=0.4, color="#0f766e")
    ax2.axhline(0.5, color=COLOR_ACCENT, linestyle="--", lw=1.5, label="Threshold (NSE=0.5)")
    ax2.set_title("(b) Model Efficiency (Nash-Sutcliffe)", fontsize=FONT_TITLE, fontweight="bold")
    ax2.set_xlabel("RMSE (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax2.set_ylabel("NSE", fontsize=FONT_LABEL, fontweight="bold")
    ax2.set_ylim(-0.5, 1.05)
    ax2.tick_params(labelsize=FONT_TICK)
    ax2.legend(loc='lower right', fontsize=FONT_LEGEND)

    # (c) Mineral Saturation Indices
    si_cols = [c for c in df.columns if c.startswith("si_")]
    if si_cols:
        labels = [c.replace("si_", "").capitalize() for c in si_cols]
        ax3.boxplot([df[c].dropna() for c in si_cols], patch_artist=True,
                    boxprops=dict(facecolor="#f1f5f9", color="#334155"),
                    medianprops=dict(color=COLOR_ACCENT))
        ax3.set_xticklabels(labels, rotation=45, fontsize=FONT_TICK)
        ax3.axhline(0, color="black", ls="-", lw=1)
        ax3.axhline(0.2, color="gray", ls=":", alpha=0.5)
        ax3.axhline(-0.2, color="gray", ls=":", alpha=0.5)
        ax3.set_title("(c) Mineral Saturation Indices", fontsize=FONT_TITLE, fontweight="bold")
        ax3.set_ylabel("Saturation Index (SI)", fontsize=FONT_LABEL, fontweight="bold")
        ax3.tick_params(axis='y', labelsize=FONT_TICK)

    # Hide the 4th panel
    ax4.axis("off")

    for ax in [ax1, ax2, ax3]:
        ax.grid(True, which="major", ls="-", alpha=GRID_ALPHA)

    plt.suptitle("Supplementary Figure S2: Geochemical & Thermodynamic Validation",
                 fontsize=FONT_TITLE + 2, fontweight="bold", y=0.98)
    fig.tight_layout()
    _save_supp(fig, "Manuscript_Supp_FigS2_Geochemical_Validation.png")

def plot_s4_field_performance() -> None:
    """Figure S3: Field Data Discovery Match (Ghana Site Residuals)."""
    path = RESULT_DIR / "field_discovery_results.csv"
    if not path.exists():
        print("Skipping S3: Field discovery results not found.")
        return

    df = pd.read_csv(path)
    df["site"] = df["edge_id"].apply(lambda x: "Manu" if "Manu" in x else "Talensi")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    sites = [("Manu", ax1, COLOR_PRIMARY, "Lower Anayari (Manu)"),
             ("Talensi", ax2, COLOR_WARNING, "Talensi Mining Area")]

    for site_name, ax, color, title in sites:
        data = df[df["site"] == site_name]["chemistry_r2"]
        ax.hist(data, bins=15, color=color, alpha=0.7, edgecolor="white")
        ax.axvline(data.median(), color=COLOR_ACCENT, ls="--", lw=2,
                   label=f"Median $R^2$: {data.median():.2f}")
        ax.set_title(f"Discovery Fit: {title}", fontsize=FONT_TITLE, fontweight="bold")
        ax.set_xlabel(r"Chemical Match ($R^2$)", fontsize=FONT_LABEL, fontweight="bold")
        ax.set_ylabel("Edge Count", fontsize=FONT_LABEL, fontweight="bold")
        ax.set_xlim(0, 1.05)
        ax.legend(loc="upper left", fontsize=FONT_LEGEND)
        ax.grid(True, which="major", ls="-", alpha=GRID_ALPHA)
        ax.tick_params(labelsize=FONT_TICK)

    fig.tight_layout()
    _save_supp(fig, "Manuscript_Supp_FigS3_Ghana_Field_Residuals.png")

def main() -> None:
    print("Generating M2 Supplementary Figures for Manuscript...")
    plot_s1_age_parity()
    plot_s2_geochemical_validation()
    plot_s4_field_performance()
    print("Done. Figures saved to 'figures/Manuscript_Ready/'.")

if __name__ == "__main__":
    main()
