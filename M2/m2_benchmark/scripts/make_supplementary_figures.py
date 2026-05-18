from __future__ import annotations

from pathlib import Path
import uuid
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import MaxNLocator

# --- Path Configuration ---
BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXTERNAL_DIR = BENCHMARK_ROOT / "external"
RESULT_DIR = BENCHMARK_ROOT / "results"
M3_RESULT_DIR = PROJECT_ROOT / "M3" / "m3_age_benchmark" / "results"
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


def _load_public_age_validation() -> tuple[pd.DataFrame, str]:
    m3_candidates = [
        (M3_RESULT_DIR / "m3_tracerlpm_strict_parity_full.csv", "M3 strict TracerLPM parity"),
        (M3_RESULT_DIR / "m3_design_matrix_results.csv", "M3 design-matrix USGS benchmark"),
        (M3_RESULT_DIR / "m3_phase4_screened_full_results.csv", "M3 full screened USGS benchmark"),
        (M3_RESULT_DIR / "m3_usgs_benchmark_results.csv", "M3 primary USGS benchmark"),
    ]
    for path, label in m3_candidates:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if {"ref_age", "est_age_multi"}.issubset(df.columns):
            preferred_scenario = None
            if "scenario_id" in df.columns:
                scenarios = set(df["scenario_id"].dropna())
                for pref in ("tracerlpm_strict_parity", "tracerlpm_parity_hier_oldwater", "screened_dgm_gases", "parity_reported_corrected"):
                    if pref in scenarios:
                        preferred_scenario = pref
                        break
                if preferred_scenario:
                    df = df[df["scenario_id"] == preferred_scenario].copy()
                    label = f"{label} ({preferred_scenario})"
            out = pd.DataFrame(
                {
                    "reference_mean_age_years": pd.to_numeric(df["ref_age"], errors="coerce"),
                    "hydrosheaf_age_years": pd.to_numeric(df["est_age_multi"], errors="coerce"),
                    "supported_tracers": df.get("tracer_mode", pd.Series("", index=df.index)).astype(str),
                    "log10_error": pd.to_numeric(df.get("log10_error", pd.Series(np.nan, index=df.index)), errors="coerce"),
                    "within_factor_2": df.get("within_factor_2", pd.Series(np.nan, index=df.index)),
                    "inside_hydrosheaf_ci": df.get("within_factor_10", pd.Series(np.nan, index=df.index)),
                }
            )
            out = out.dropna(subset=["reference_mean_age_years", "hydrosheaf_age_years"])
            if not out.empty:
                return out, label

    legacy_path = EXTERNAL_DIR / "usgs_age" / "results" / "usgs_age_validation.csv"
    if legacy_path.exists():
        return pd.read_csv(legacy_path), "M2 legacy E1 USGS age benchmark"
    return pd.DataFrame(), ""


def plot_s1_age_parity() -> None:
    """Figure S1: Regional Residence-Time Validation (USGS Dataset)."""
    df, source_label = _load_public_age_validation()
    if df.empty:
        print("Skipping S1: USGS age data not found.")
        return

    df_clean = df.dropna(subset=["reference_mean_age_years", "hydrosheaf_age_years"]).copy()
    df_clean["plot_reference_age_years"] = np.maximum(df_clean["reference_mean_age_years"], 0.01)
    df_clean["plot_hydrosheaf_age_years"] = np.maximum(df_clean["hydrosheaf_age_years"], 0.01)

    def get_category(row):
        tracers = str(row["supported_tracers"])
        if "14C" in tracers: return "14C Inclusive"
        if any(x in tracers for x in ["SF6", "CFC", "3H"]): return "Modern (3H/SF6/CFC)"
        return "Other"

    df_clean["Category"] = df_clean.apply(get_category, axis=1)

    log_ref = np.log10(df_clean["plot_reference_age_years"])
    log_inf = np.log10(df_clean["plot_hydrosheaf_age_years"])

    r2_global = np.corrcoef(log_ref, log_inf)[0,1]**2
    metric_log = pd.to_numeric(
        df_clean.get("log10_error", pd.Series(np.nan, index=df_clean.index)),
        errors="coerce",
    ).abs()
    mae_log = float(metric_log.mean()) if metric_log.notna().any() else np.mean(np.abs(log_ref - log_inf))
    within_ci = (df_clean["inside_hydrosheaf_ci"] == True).mean() * 100 if "inside_hydrosheaf_ci" in df_clean else np.nan
    within_raw = df_clean.get("within_factor_2", pd.Series(np.nan, index=df_clean.index))
    within_numeric = pd.to_numeric(within_raw, errors="coerce")
    within_bool = within_raw.astype(str).str.lower().map({"true": 1.0, "false": 0.0})
    within_factor_2 = within_numeric.where(within_numeric.notna(), within_bool).fillna(0.0).mean() * 100

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(4, 4, left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.05, hspace=0.05)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_histx = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_histy = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    colors = {"Modern (3H/SF6/CFC)": COLOR_PRIMARY, "14C Inclusive": COLOR_ACCENT, "Other": "#64748b"}
    for cat, color in colors.items():
        sub = df_clean[df_clean["Category"] == cat]
        ax_main.scatter(sub["plot_reference_age_years"], sub["plot_hydrosheaf_age_years"],
                        c=color, alpha=0.5, s=30, label=cat, edgecolors='white', linewidth=0.3)

    lim_max = max(
        1e5,
        float(df_clean["plot_reference_age_years"].max()),
        float(df_clean["plot_hydrosheaf_age_years"].max()),
    ) * 1.25
    lims = np.logspace(-2, np.log10(lim_max), 120)
    ax_main.plot(lims, lims, color="#1f2937", ls="--", alpha=0.8, lw=1.5, label="1:1 Line")
    ax_main.fill_between(lims, lims/10, lims*10, color="#94a3b8", alpha=0.15, label="±1 Order of Mag")

    bins = np.logspace(-2, np.log10(lim_max), 40)
    ax_histx.hist(df_clean["plot_reference_age_years"], bins=bins, color="#cbd5e1", edgecolor="white")
    ax_histy.hist(df_clean["plot_hydrosheaf_age_years"], bins=bins, orientation='horizontal', color="#cbd5e1", edgecolor="white")

    ax_main.set_xscale("log")
    ax_main.set_yscale("log")
    ax_main.set_xlim(0.01, lim_max)
    ax_main.set_ylim(0.01, lim_max)
    ax_main.set_xlabel("Reference Tracer Age (Years)", fontsize=FONT_LABEL, fontweight="bold")
    ax_main.set_ylabel("Hydrosheaf Inferred Age (Years)", fontsize=FONT_LABEL, fontweight="bold")
    ax_main.tick_params(labelsize=FONT_TICK)

    stats_text = (f"N = {len(df_clean)} ({metric_log.notna().sum()} finite log errors)\n"
                  f"Global $R^2 = {r2_global:.2f}$\n"
                  f"Mean |log10 err| = {mae_log:.2f}\n"
                  f"Within factor 2 = {within_factor_2:.1f}%")
    if np.isfinite(within_ci):
        stats_text += f"\nLegacy CI/proxy = {within_ci:.1f}%"
    ax_main.text(0.05, 0.95, stats_text, transform=ax_main.transAxes,
                 verticalalignment='top', fontsize=FONT_ANNOTATE, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="#d1d5db", alpha=0.9))
    ax_main.text(0.05, 0.04, source_label, transform=ax_main.transAxes,
                 fontsize=FONT_ANNOTATE, color="#64748b", ha="left", va="bottom")

    ax_main.legend(loc="lower right", fontsize=FONT_LEGEND, frameon=True)
    ax_main.grid(True, which="major", ls="-", alpha=GRID_ALPHA)

    ax_histx.axis("off")
    ax_histy.axis("off")
    ax_histx.set_title("Supplementary Figure S1: Public USGS Residence-Time Screening Check",
                       fontsize=FONT_TITLE, fontweight="bold", pad=20)

    _save_supp(fig, "Manuscript_Supp_FigS1_Public_Age_Validation.png")

def plot_s2_geochemical_validation() -> None:
    """Figure S2: Geochemical Validation Suite (RMSE, NSE, and SI)."""
    path = BENCHMARK_ROOT / "results" / "phreeqc_forward_validation.csv"
    if not path.exists():
        print(f"Skipping S2: PHREEQC forward data not found at {path}")
        return

    df = pd.read_csv(path)
    # 1x3 Grid for clean supplementary layout
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    ax1, ax2, ax3 = axes

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

    # (c) Percent Bias (PBIAS) Distribution
    pbias_clean = df["pbias_percent"].dropna()
    ax3.hist(pbias_clean, bins=20, color=COLOR_WARNING, alpha=0.7, edgecolor="white")
    ax3.axvline(pbias_clean.median(), color=COLOR_ACCENT, ls="--", lw=2,
                label=f"Median: {pbias_clean.median():.1f}%")
    ax3.axvline(0, color="black", ls="-", lw=1)
    ax3.set_title("(c) Percent Bias (PBIAS)", fontsize=FONT_TITLE, fontweight="bold")
    ax3.set_xlabel("PBIAS (%)", fontsize=FONT_LABEL, fontweight="bold")
    ax3.set_ylabel("Frequency", fontsize=FONT_LABEL, fontweight="bold")
    ax3.tick_params(labelsize=FONT_TICK)
    ax3.legend(fontsize=FONT_LEGEND)

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
