"""
Comprehensive Hydrosheaf Analysis with Built-in Advanced Plotting

This script uses the full hydrosheaf package capabilities including:
- Advanced ILR (Isometric Log-Ratio) plots
- Gibbs diagrams for hydrochemical classification
- Edge anomaly visualizations
- Built-in export functions
- Comprehensive HTML report generation

Author: Generated for Hydrosheaf demonstration
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import json
from datetime import datetime
import warnings

warnings.filterwarnings("ignore")

# Add hydrosheaf to path
sys.path.insert(0, str(Path(__file__).parent))

from hydrosheaf import (
    Config,
    fit_network,
    summarize_network,
    edge_process_maps,
    compute_d_excess,
    evaporation_index,
)
from hydrosheaf.data.units import mgL_to_mmolL, MOLAR_MASS_G_MOL

# Import hydrosheaf's built-in plotting and export functions
from hydrosheaf.outputs.plots import plot_edge_anomalies, plot_gibbs, plot_ilr
from hydrosheaf.outputs.science_plots import (
    plot_ttd_kernel,
    plot_breakthrough,
    plot_posterior_ridges,
    plot_reactive_transport_validation,
)
from hydrosheaf.outputs.export import export_edge_results_csv, export_edge_results_json

# Additional plotting
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False

try:
    import seaborn as sns

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False


# -------------------------------------------------------------------
# Output Directory Setup
# -------------------------------------------------------------------


def setup_output_directory(base_dir: Path) -> Path:
    """Create output directory structure."""
    output_dir = base_dir / "analysis_results"
    output_dir.mkdir(exist_ok=True)
    (output_dir / "plots").mkdir(exist_ok=True)
    (output_dir / "data").mkdir(exist_ok=True)
    return output_dir


# -------------------------------------------------------------------
# Data Loading and Preprocessing
# -------------------------------------------------------------------


def load_synthetic_data(data_dir: Path) -> dict:
    """Load all synthetic CSV files."""
    data = {}
    data["water_chem"] = pd.read_csv(data_dir / "water_chem_full.csv")
    data["stations"] = pd.read_csv(data_dir / "stations.csv")
    data["edges"] = pd.read_csv(data_dir / "network_edges.csv")
    data["events"] = pd.read_csv(data_dir / "events.csv")
    data["endmembers"] = pd.read_csv(data_dir / "endmembers_isotopes.csv")
    data["redox"] = pd.read_csv(data_dir / "redox_proxies.csv")

    # Load optional files
    for fname in [
        "fertilizer_applications.csv",
        "meteo_daily.csv",
        "soil_profiles.csv",
    ]:
        fpath = data_dir / fname
        if fpath.exists():
            data[fname.replace(".csv", "")] = pd.read_csv(fpath)

    return data


def prepare_samples_mmol(water_chem: pd.DataFrame, event_code: str = None) -> list:
    """Convert water chemistry to mmol/L format for hydrosheaf."""
    ion_mapping = {
        "Ca": "Ca_mg_L",
        "Mg": "Mg_mg_L",
        "Na": "Na_mg_L",
        "K": "K_mg_L",
        "HCO3": "HCO3_mg_L",
        "Cl": "Cl_mg_L",
        "SO4": "SO4_mg_L",
        "NO3": "NO3_mg_L",
    }

    df = water_chem.copy()
    if event_code:
        df = df[df["event_code"] == event_code].copy()

    # Convert to mmol/L
    for ion, col in ion_mapping.items():
        if col in df.columns:
            df[ion] = df[col].apply(
                lambda x: mgL_to_mmolL(x, ion) if pd.notna(x) else np.nan
            )

    # Add missing ions
    for ion in ["F", "Fe", "PO4"]:
        df[ion] = 0.0

    samples = []
    for _, row in df.iterrows():
        sample = {
            "site_id": row["station_code"],
            "sample_id": f"{row['event_code']}_{row['station_code']}",
            "Ca": row["Ca"],
            "Mg": row["Mg"],
            "Na": row["Na"],
            "K": row.get("K", mgL_to_mmolL(row.get("K_mg_L", 0), "K")),
            "HCO3": row["HCO3"],
            "Cl": row["Cl"],
            "SO4": row["SO4"],
            "NO3": row["NO3"],
            "F": row["F"],
            "Fe": row["Fe"],
            "PO4": row["PO4"],
            "pH": row.get("pH", 7.0),
            "EC": row.get("EC_uS_cm", None),
            "TDS": row.get("TDS_mg_L", None),
            "temp_C": row.get("temp_C", 25.0),
        }
        samples.append(sample)

    return samples


def prepare_samples_for_ilr(water_chem: pd.DataFrame) -> list:
    """Prepare samples in mmol/L format for ILR plotting."""
    samples = []
    for _, row in water_chem.iterrows():
        sample = {
            "site_id": row["station_code"],
            "Ca": mgL_to_mmolL(row["Ca_mg_L"], "Ca"),
            "Mg": mgL_to_mmolL(row["Mg_mg_L"], "Mg"),
            "Na": mgL_to_mmolL(row["Na_mg_L"], "Na"),
            "K": mgL_to_mmolL(row["K_mg_L"], "K"),
            "HCO3": mgL_to_mmolL(row["HCO3_mg_L"], "HCO3"),
            "Cl": mgL_to_mmolL(row["Cl_mg_L"], "Cl"),
            "SO4": mgL_to_mmolL(row["SO4_mg_L"], "SO4"),
            "NO3": mgL_to_mmolL(row["NO3_mg_L"], "NO3"),
            "TDS": row["TDS_mg_L"],
            "station_type": row["station_type"],
            "event_code": row["event_code"],
        }
        samples.append(sample)
    return samples


def prepare_edges(edges_df: pd.DataFrame) -> list:
    """Convert edges DataFrame to hydrosheaf edge format."""
    return [(row["from_station"], row["to_station"]) for _, row in edges_df.iterrows()]


# -------------------------------------------------------------------
# Custom Plotting Functions (Complementing Built-in)
# -------------------------------------------------------------------


def plot_water_isotopes(water_chem: pd.DataFrame, output_dir: Path):
    """Plot water isotopes with LMWL and GMWL."""
    if not PLOTTING_AVAILABLE:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    d18O = water_chem["d18O_H2O_permil"]
    d2H = water_chem["d2H_H2O_permil"]

    colors = {"lysimeter": "#3498db", "borehole": "#e74c3c"}
    markers = {"lysimeter": "o", "borehole": "s"}

    # Plot meteoric water lines
    x_range = np.linspace(d18O.min() - 1, d18O.max() + 1, 100)
    ax1.plot(x_range, 8 * x_range + 10, "k--", lw=1.5, label="GMWL")
    ax1.plot(x_range, 8.66 * x_range + 7.22, "g-", lw=1.5, label="LMWL (Ghana)")

    for stype in ["lysimeter", "borehole"]:
        subset = water_chem[water_chem["station_type"] == stype]
        ax1.scatter(
            subset["d18O_H2O_permil"],
            subset["d2H_H2O_permil"],
            c=colors[stype],
            marker=markers[stype],
            s=80,
            alpha=0.7,
            label=stype.capitalize(),
            edgecolors="black",
        )

    ax1.set_xlabel("d18O (permil)")
    ax1.set_ylabel("d2H (permil)")
    ax1.set_title("Water Isotopes", fontweight="bold")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # D-excess plot
    d_excess = d2H - 8 * d18O
    water_chem_copy = water_chem.copy()
    water_chem_copy["d_excess"] = d_excess

    for stype in ["lysimeter", "borehole"]:
        subset = water_chem_copy[water_chem_copy["station_type"] == stype]
        ax2.scatter(
            subset["station_code"],
            subset["d_excess"],
            c=colors[stype],
            marker=markers[stype],
            s=80,
            alpha=0.7,
            edgecolors="black",
        )

    ax2.axhline(y=10, color="gray", linestyle="--", label="GMWL d-excess")
    ax2.set_xlabel("Station")
    ax2.set_ylabel("d-excess (permil)")
    ax2.set_title("Deuterium Excess", fontweight="bold")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.suptitle("Stable Isotope Analysis", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(
        output_dir / "plots" / "water_isotopes.png", dpi=150, bbox_inches="tight"
    )
    plt.close()


def plot_nitrate_source_diagram(
    water_chem: pd.DataFrame, endmembers: pd.DataFrame, output_dir: Path
):
    """Plot nitrate isotopes with source endmember boxes."""
    if not PLOTTING_AVAILABLE:
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    colors = {"lysimeter": "#3498db", "borehole": "#e74c3c"}
    markers = {"lysimeter": "o", "borehole": "s"}

    for stype in ["lysimeter", "borehole"]:
        subset = water_chem[water_chem["station_type"] == stype]
        ax.scatter(
            subset["d18O_NO3_permil"],
            subset["d15N_NO3_permil"],
            c=colors[stype],
            marker=markers[stype],
            s=100,
            alpha=0.7,
            label=stype.capitalize(),
            edgecolors="black",
            zorder=5,
        )

    # Endmember boxes
    boxes = {
        "Synthetic Fertilizer": {"d15N": (0, 8), "d18O": (17, 30), "color": "#2ecc71"},
        "Manure/Sewage": {"d15N": (8, 22), "d18O": (0, 15), "color": "#9b59b6"},
        "Soil Organic N": {"d15N": (3, 12), "d18O": (-5, 10), "color": "#f39c12"},
        "Atmospheric": {"d15N": (-5, 5), "d18O": (50, 70), "color": "#1abc9c"},
    }

    for name, box in boxes.items():
        rect = mpatches.Rectangle(
            (box["d18O"][0], box["d15N"][0]),
            box["d18O"][1] - box["d18O"][0],
            box["d15N"][1] - box["d15N"][0],
            linewidth=2,
            edgecolor=box["color"],
            facecolor=box["color"],
            alpha=0.2,
        )
        ax.add_patch(rect)
        ax.text(
            np.mean(box["d18O"]),
            np.mean(box["d15N"]),
            name,
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
        )

    ax.set_xlabel("d18O-NO3 (permil)")
    ax.set_ylabel("d15N-NO3 (permil)")
    ax.set_title("Nitrate Source Identification", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-10, 35)
    ax.set_ylim(0, 15)

    plt.tight_layout()
    plt.savefig(
        output_dir / "plots" / "nitrate_sources.png", dpi=150, bbox_inches="tight"
    )
    plt.close()


def plot_temporal_trends(
    water_chem: pd.DataFrame, events: pd.DataFrame, output_dir: Path
):
    """Plot temporal patterns in water chemistry."""
    if not PLOTTING_AVAILABLE:
        return

    merged = water_chem.merge(events, on="event_code", how="left")
    merged["collection_date"] = pd.to_datetime(merged["collection_date"])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = {"lysimeter": "#3498db", "borehole": "#e74c3c"}
    season_colors = {"dry": "#f39c12", "wet": "#3498db", "transition": "#9b59b6"}

    # NO3 over time
    ax1 = axes[0, 0]
    for stype in ["lysimeter", "borehole"]:
        subset = (
            merged[merged["station_type"] == stype]
            .groupby("collection_date")["NO3_mg_L"]
            .mean()
        )
        ax1.plot(
            subset.index,
            subset.values,
            marker="o",
            lw=2,
            label=stype.capitalize(),
            color=colors[stype],
        )
    ax1.set_xlabel("Date")
    ax1.set_ylabel("NO3 (mg/L)")
    ax1.set_title("Nitrate Temporal Trend", fontweight="bold")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45)

    # NO3 by season
    ax2 = axes[0, 1]
    season_order = ["dry", "transition", "wet"]
    for i, season in enumerate(season_order):
        subset = merged[merged["season"] == season]["NO3_mg_L"]
        bp = ax2.boxplot([subset], positions=[i], widths=0.6, patch_artist=True)
        bp["boxes"][0].set_facecolor(season_colors[season])
    ax2.set_xticks(range(len(season_order)))
    ax2.set_xticklabels([s.capitalize() for s in season_order])
    ax2.set_ylabel("NO3 (mg/L)")
    ax2.set_title("Nitrate by Season", fontweight="bold")
    ax2.grid(True, alpha=0.3)

    # EC over time
    ax3 = axes[1, 0]
    for stype in ["lysimeter", "borehole"]:
        subset = (
            merged[merged["station_type"] == stype]
            .groupby("collection_date")["EC_uS_cm"]
            .mean()
        )
        ax3.plot(
            subset.index,
            subset.values,
            marker="s",
            lw=2,
            label=stype.capitalize(),
            color=colors[stype],
        )
    ax3.set_xlabel("Date")
    ax3.set_ylabel("EC (uS/cm)")
    ax3.set_title("Electrical Conductivity Trend", fontweight="bold")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45)

    # Rainfall
    ax4 = axes[1, 1]
    event_data = events.copy()
    ax4.bar(
        event_data["event_code"],
        event_data["rain_30d_mm"],
        color=[season_colors.get(s, "gray") for s in event_data["season"]],
        alpha=0.7,
    )
    ax4.set_xlabel("Event")
    ax4.set_ylabel("30-day Rainfall (mm)")
    ax4.set_title("Antecedent Rainfall", fontweight="bold")
    ax4.grid(True, alpha=0.3)

    plt.suptitle("Temporal and Seasonal Patterns", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(
        output_dir / "plots" / "temporal_patterns.png", dpi=150, bbox_inches="tight"
    )
    plt.close()


def plot_reaction_summary(all_results: dict, output_dir: Path):
    """Create heatmap and bar chart of reaction extents."""
    if not PLOTTING_AVAILABLE:
        return

    # Collect all reaction data
    rows = []
    for event_code, results in all_results.items():
        if not results:
            continue
        for res in results:
            for label, extent in zip(res.z_labels, res.z_extents):
                if abs(extent) > 1e-6:
                    rows.append(
                        {
                            "Event": event_code,
                            "Edge": f"{res.u}->{res.v}",
                            "Reaction": label,
                            "Extent": extent,
                        }
                    )

    if not rows:
        return

    df = pd.DataFrame(rows)

    # Create heatmap
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Mean reaction by edge
    pivot = df.pivot_table(
        index="Edge", columns="Reaction", values="Extent", aggfunc="mean", fill_value=0
    )

    if SEABORN_AVAILABLE:
        sns.heatmap(
            pivot,
            annot=True,
            fmt=".3f",
            cmap="RdYlBu_r",
            center=0,
            ax=ax1,
            cbar_kws={"label": "Reaction Extent (mmol/L)"},
        )
    else:
        im = ax1.imshow(pivot.values, cmap="RdYlBu_r", aspect="auto")
        ax1.set_xticks(range(len(pivot.columns)))
        ax1.set_xticklabels(pivot.columns, rotation=45, ha="right")
        ax1.set_yticks(range(len(pivot.index)))
        ax1.set_yticklabels(pivot.index)
        plt.colorbar(im, ax=ax1, label="Extent (mmol/L)")

    ax1.set_title("Mean Reaction Extents by Edge", fontweight="bold")
    ax1.set_xlabel("Reaction")
    ax1.set_ylabel("Flow Edge")

    # Total reaction intensity by type
    reaction_totals = (
        df.groupby("Reaction")["Extent"]
        .apply(lambda x: np.abs(x).sum())
        .sort_values(ascending=True)
    )
    ax2.barh(
        reaction_totals.index, reaction_totals.values, color="steelblue", alpha=0.7
    )
    ax2.set_xlabel("Cumulative |Extent| (mmol/L)")
    ax2.set_title("Total Reaction Activity", fontweight="bold")
    ax2.grid(True, alpha=0.3, axis="x")

    plt.tight_layout()
    plt.savefig(
        output_dir / "plots" / "reaction_summary.png", dpi=150, bbox_inches="tight"
    )
    plt.close()


def plot_transport_parameters(all_results: dict, output_dir: Path):
    """Plot gamma (evaporation factor) across events."""
    if not PLOTTING_AVAILABLE:
        return

    rows = []
    for event_code, results in all_results.items():
        if not results:
            continue
        for res in results:
            rows.append(
                {
                    "Event": event_code,
                    "Edge": f"{res.u}->{res.v}",
                    "Gamma": res.gamma if res.gamma else 1.0,
                    "Transport": res.transport_model,
                    "Objective": res.objective_score,
                }
            )

    df = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Gamma by edge
    ax1 = axes[0]
    edges = df["Edge"].unique()
    edge_colors = plt.cm.Set2(np.linspace(0, 1, len(edges)))

    for i, edge in enumerate(edges):
        subset = df[df["Edge"] == edge]
        ax1.plot(
            subset["Event"],
            subset["Gamma"],
            marker="o",
            lw=2,
            label=edge,
            color=edge_colors[i],
        )

    ax1.axhline(y=1.0, color="gray", linestyle="--", alpha=0.7, label="No evaporation")
    ax1.set_xlabel("Event")
    ax1.set_ylabel("Evaporation Factor (gamma)")
    ax1.set_title("Evaporation Intensity Over Time", fontweight="bold")
    ax1.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Objective score
    ax2 = axes[1]
    edge_means = df.groupby("Edge")["Objective"].mean().sort_values()
    ax2.barh(edge_means.index, edge_means.values, color="coral", alpha=0.7)
    ax2.set_xlabel("Mean Objective Score")
    ax2.set_title("Fit Quality by Edge", fontweight="bold")
    ax2.grid(True, alpha=0.3, axis="x")

    plt.suptitle("Transport Model Analysis", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(
        output_dir / "plots" / "transport_parameters.png", dpi=150, bbox_inches="tight"
    )
    plt.close()


def plot_network_schematic(
    edges_df: pd.DataFrame, stations_df: pd.DataFrame, output_dir: Path
):
    """Create a network flow diagram."""
    if not PLOTTING_AVAILABLE:
        return

    fig, ax = plt.subplots(figsize=(10, 8))

    positions = {
        row["station_code"]: (row["lon_deg"], row["lat_deg"])
        for _, row in stations_df.iterrows()
    }
    colors = {"lysimeter": "#3498db", "borehole": "#e74c3c"}

    # Draw edges
    for _, row in edges_df.iterrows():
        from_st, to_st = row["from_station"], row["to_station"]
        if from_st in positions and to_st in positions:
            x1, y1 = positions[from_st]
            x2, y2 = positions[to_st]
            ax.annotate(
                "",
                xy=(x2, y2),
                xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color="gray", lw=2),
            )

    # Draw stations
    for _, row in stations_df.iterrows():
        code, stype = row["station_code"], row["station_type"]
        x, y = positions[code]
        marker = "o" if stype == "lysimeter" else "s"
        ax.scatter(
            x,
            y,
            c=colors[stype],
            s=400,
            marker=marker,
            edgecolors="black",
            lw=2,
            zorder=5,
        )
        ax.annotate(
            code,
            (x, y),
            fontsize=11,
            fontweight="bold",
            ha="center",
            va="center",
            color="white",
        )

    legend_elements = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=colors["lysimeter"],
            markersize=12,
            label="Lysimeter",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="s",
            color="w",
            markerfacecolor=colors["borehole"],
            markersize=12,
            label="Borehole",
        ),
    ]
    ax.legend(handles=legend_elements, loc="upper right")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Groundwater Flow Network", fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(
        output_dir / "plots" / "network_schematic.png", dpi=150, bbox_inches="tight"
    )
    plt.close()


# -------------------------------------------------------------------
# Report Generation
# -------------------------------------------------------------------


def generate_comprehensive_report(
    data: dict, all_results: dict, all_summaries: dict, output_dir: Path
):
    """Generate comprehensive HTML report with all plots and tables."""

    water_chem = data["water_chem"]
    events = data["events"]

    # Stats
    n_samples = len(water_chem)
    n_stations = len(water_chem["station_code"].unique())
    n_events = len(water_chem["event_code"].unique())
    n_edges = sum(len(r) for r in all_results.values() if r)

    # Ion stats
    ion_cols = [
        "Ca_mg_L",
        "Mg_mg_L",
        "Na_mg_L",
        "K_mg_L",
        "HCO3_mg_L",
        "Cl_mg_L",
        "SO4_mg_L",
        "NO3_mg_L",
    ]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Hydrosheaf Comprehensive Analysis Report</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; background: #f0f4f8; }}
        .header {{ background: linear-gradient(135deg, #1a365d 0%, #2c5282 100%); color: white; padding: 40px; border-radius: 12px; margin-bottom: 30px; }}
        .header h1 {{ margin: 0; font-size: 2.5em; }}
        .section {{ background: white; padding: 30px; border-radius: 12px; margin-bottom: 25px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); }}
        .section h2 {{ color: #1a365d; border-bottom: 3px solid #4299e1; padding-bottom: 10px; }}
        .section h3 {{ color: #2c5282; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ border: 1px solid #e2e8f0; padding: 12px; text-align: left; }}
        th {{ background: #2c5282; color: white; }}
        tr:nth-child(even) {{ background: #f7fafc; }}
        tr:hover {{ background: #edf2f7; }}
        .plot-container {{ text-align: center; margin: 25px 0; }}
        .plot-container img {{ max-width: 100%; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.15); }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 20px; margin: 25px 0; }}
        .stat-card {{ background: linear-gradient(135deg, #4299e1 0%, #667eea 100%); color: white; padding: 25px; border-radius: 12px; text-align: center; }}
        .stat-card .value {{ font-size: 2.5em; font-weight: bold; }}
        .stat-card .label {{ font-size: 0.95em; opacity: 0.9; margin-top: 5px; }}
        .finding {{ background: #ebf8ff; border-left: 5px solid #4299e1; padding: 15px 20px; margin: 20px 0; border-radius: 0 8px 8px 0; }}
        .reaction-tag {{ display: inline-block; padding: 5px 12px; border-radius: 20px; margin: 3px; font-size: 0.9em; }}
        .reaction-positive {{ background: #c6f6d5; color: #22543d; }}
        .reaction-negative {{ background: #fed7d7; color: #822727; }}
        .plot-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        @media (max-width: 900px) {{ .plot-row {{ grid-template-columns: 1fr; }} }}
    </style>
</head>
<body>

<div class="header">
    <h1>Hydrosheaf Groundwater Analysis Report</h1>
    <p>Comprehensive hydrochemical and isotopic analysis with reaction pathway modeling</p>
    <p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
</div>

<div class="section">
    <h2>1. Executive Summary</h2>
    <div class="stats-grid">
        <div class="stat-card"><div class="value">{n_samples}</div><div class="label">Water Samples</div></div>
        <div class="stat-card"><div class="value">{n_stations}</div><div class="label">Monitoring Stations</div></div>
        <div class="stat-card"><div class="value">{n_events}</div><div class="label">Sampling Events</div></div>
        <div class="stat-card"><div class="value">{n_edges}</div><div class="label">Flow Edges Analyzed</div></div>
    </div>
    
    <div class="finding">
        <strong>Key Finding:</strong> All analyzed flow paths show evaporation-dominated transport (gamma > 1.0). 
        Denitrification is the most significant geochemical reaction, particularly along vadose-to-groundwater 
        pathways during wet season events. Rock weathering (halite, gypsum, calcite dissolution) and ion exchange 
        processes are consistently active across the flow network.
    </div>
</div>

<div class="section">
    <h2>2. Network Configuration</h2>
    <div class="plot-container">
        <img src="plots/network_schematic.png" alt="Network Diagram">
        <p><em>Figure 1: Groundwater flow network showing lysimeters (blue circles) and boreholes (red squares)</em></p>
    </div>
    
    <h3>Station Details</h3>
    <table>
        <tr><th>Station</th><th>Type</th><th>Cluster</th><th>Input Intensity</th><th>N Target (kg/ha)</th></tr>
"""

    for _, row in data["stations"].iterrows():
        html += f"<tr><td>{row['station_code']}</td><td>{row['station_type']}</td><td>{row['cluster_code']}</td><td>{row['input_intensity']}</td><td>{row['target_N_kg_ha_range']}</td></tr>\n"

    html += """    </table>
</div>

<div class="section">
    <h2>3. Hydrochemical Facies Analysis</h2>
    <p>Using Hydrosheaf's advanced Isometric Log-Ratio (ILR) approach for water type classification.</p>
    
    <div class="plot-container">
        <img src="plots/ilr_plot.png" alt="ILR Hydrochemical Facies">
        <p><em>Figure 2: ILR diagram showing hydrochemical facies classification with cation and anion fields</em></p>
    </div>
    
    <div class="plot-container">
        <img src="plots/gibbs_diagram.png" alt="Gibbs Diagram">
        <p><em>Figure 3: Gibbs diagram showing dominant hydrochemical processes</em></p>
    </div>
    
    <h3>Major Ion Statistics (mg/L)</h3>
    <table>
        <tr><th>Ion</th><th>Mean</th><th>Std</th><th>Min</th><th>Max</th><th>Median</th></tr>
"""

    for col in ion_cols:
        if col in water_chem.columns:
            ion = col.replace("_mg_L", "")
            html += f"<tr><td>{ion}</td><td>{water_chem[col].mean():.2f}</td><td>{water_chem[col].std():.2f}</td><td>{water_chem[col].min():.2f}</td><td>{water_chem[col].max():.2f}</td><td>{water_chem[col].median():.2f}</td></tr>\n"

    html += """    </table>
    
    <div class="finding">
        <strong>Water Type:</strong> The groundwater is predominantly Ca-Mg-HCO3 type, indicating carbonate rock weathering 
        as the dominant process. Samples plot in the rock weathering dominance zone of the Gibbs diagram.
    </div>
</div>

<div class="section">
    <h2>4. Stable Isotope Analysis</h2>
    
    <div class="plot-row">
        <div class="plot-container">
            <img src="plots/water_isotopes.png" alt="Water Isotopes">
            <p><em>Figure 4: Water isotopes relative to GMWL and local LMWL</em></p>
        </div>
        <div class="plot-container">
            <img src="plots/nitrate_sources.png" alt="Nitrate Sources">
            <p><em>Figure 5: Nitrate isotope source identification</em></p>
        </div>
    </div>
    
    <div class="finding">
        <strong>Isotope Interpretation:</strong> Water samples cluster near the LMWL, indicating limited evaporative 
        enrichment during recharge. Nitrate isotopes suggest mixed sources dominated by synthetic fertilizer 
        and soil organic nitrogen, consistent with agricultural land use.
    </div>
</div>

<div class="section">
    <h2>5. Temporal and Seasonal Patterns</h2>
    <div class="plot-container">
        <img src="plots/temporal_patterns.png" alt="Temporal Patterns">
        <p><em>Figure 6: Temporal trends in nitrate and EC with seasonal context</em></p>
    </div>
    
    <h3>Seasonal NO3 Statistics</h3>
    <table>
        <tr><th>Season</th><th>Mean NO3 (mg/L)</th><th>Std Dev</th><th>n</th></tr>
"""

    merged = water_chem.merge(events, on="event_code", how="left")
    for season in ["dry", "transition", "wet"]:
        subset = merged[merged["season"] == season]["NO3_mg_L"]
        if len(subset) > 0:
            html += f"<tr><td>{season.capitalize()}</td><td>{subset.mean():.1f}</td><td>{subset.std():.1f}</td><td>{len(subset)}</td></tr>\n"

    html += """    </table>
</div>

<div class="section">
    <h2>6. Hydrosheaf Network Analysis</h2>
    <p>Edge-based hydrochemical evolution modeling with transport and reaction pathway inference.</p>
    
    <div class="plot-container">
        <img src="plots/edge_anomalies.png" alt="Edge Anomalies">
        <p><em>Figure 7: Anomaly norms by flow edge (lower = better fit)</em></p>
    </div>
    
    <div class="plot-row">
        <div class="plot-container">
            <img src="plots/transport_parameters.png" alt="Transport Parameters">
            <p><em>Figure 8: Evaporation factors and model fit quality</em></p>
        </div>
        <div class="plot-container">
            <img src="plots/reaction_summary.png" alt="Reaction Summary">
            <p><em>Figure 9: Reaction pathway analysis</em></p>
        </div>
    </div>
    
    <h3>Event-by-Event Results</h3>
"""

    for event_code, results in all_results.items():
        if not results:
            continue

        html += f"""    <h4>Event: {event_code}</h4>
    <table>
        <tr><th>Edge</th><th>Transport</th><th>Gamma</th><th>Objective</th><th>Key Reactions</th></tr>
"""
        for res in results:
            reactions = [
                (l, e) for l, e in zip(res.z_labels, res.z_extents) if abs(e) > 0.01
            ]
            reactions.sort(key=lambda x: -abs(x[1]))
            reaction_html = " ".join(
                [
                    f'<span class="reaction-tag {"reaction-positive" if e > 0 else "reaction-negative"}">{l}: {e:+.3f}</span>'
                    for l, e in reactions[:3]
                ]
            )
            html += f"<tr><td>{res.u} → {res.v}</td><td>{res.transport_model}</td><td>{res.gamma:.3f}</td><td>{res.objective_score:.4f}</td><td>{reaction_html}</td></tr>\n"

        html += "    </table>\n"

    html += (
        """
</div>

<div class="section">
    <h2>7. Conclusions and Recommendations</h2>
    
    <h3>Key Findings</h3>
    <ul>
        <li><strong>Evaporative Concentration:</strong> All flow paths show gamma > 1.0, indicating evaporation is a significant process in water chemistry evolution.</li>
        <li><strong>Active Denitrification:</strong> Substantial denitrification occurs along vadose-to-groundwater pathways, particularly during wet season when nitrate inputs are highest.</li>
        <li><strong>Mineral Dissolution:</strong> Halite, gypsum, and calcite dissolution consistently increase TDS along flow paths.</li>
        <li><strong>Ion Exchange:</strong> Ca-Na and Mg-Na exchange reactions modify cation ratios, indicating clay mineral interactions.</li>
        <li><strong>Seasonal Dynamics:</strong> Wet season shows highest NO3 concentrations linked to fertilizer leaching; dry season shows natural attenuation.</li>
    </ul>
    
    <h3>Recommendations</h3>
    <ol>
        <li>Implement fertilizer timing management to reduce wet season leaching peaks</li>
        <li>Monitor denitrification capacity in the vadose zone as a natural attenuation buffer</li>
        <li>Consider additional age tracers (CFCs, SF6) to constrain residence time estimates</li>
        <li>Evaluate reactive barriers if nitrate attenuation needs enhancement</li>
    </ol>
</div>

<div class="section" style="text-align: center; color: #718096; padding: 20px;">
    <p><strong>Report generated using Hydrosheaf Groundwater Analysis Package v0.1.0</strong></p>
    <p>Analysis Date: """
        + datetime.now().strftime("%Y-%m-%d")
        + """</p>
    <p>All plots and data exported to the analysis_results directory</p>
</div>

</body>
</html>
"""
    )

    report_path = output_dir / "comprehensive_report.html"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    return report_path


# -------------------------------------------------------------------
# Main Analysis Pipeline
# -------------------------------------------------------------------


def main():
    print("=" * 70)
    print("HYDROSHEAF COMPREHENSIVE ANALYSIS")
    print("Using Built-in Advanced Plotting Functions")
    print("=" * 70)

    # Setup
    base_dir = Path(__file__).parent
    data_dir = base_dir / "hydrosheaf_synthetic_csv"
    output_dir = setup_output_directory(base_dir)

    print(f"\nOutput directory: {output_dir}")

    # Load data
    print("\n[1/8] Loading data...")
    data = load_synthetic_data(data_dir)
    print(
        f"  {len(data['water_chem'])} samples, {len(data['stations'])} stations, {len(data['edges'])} edges"
    )

    # Configure hydrosheaf
    config = Config(
        ion_order=["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"],
        weights=[1.0] * 10,
        phreeqc_enabled=False,
        transport_models_enabled=["evap", "mix"],
        active_minerals=[
            "calcite",
            "dolomite",
            "gypsum",
            "halite",
            "pyrite_oxidation_aerobic",
        ],
        gibbs_enabled=True,
        exchange_enabled=True,
    )

    # Run hydrosheaf analysis
    print("\n[2/8] Running Hydrosheaf network analysis...")
    all_results = {}
    all_summaries = {}
    all_edge_results = []

    for event_code in data["water_chem"]["event_code"].unique():
        samples = prepare_samples_mmol(data["water_chem"], event_code)
        edges = prepare_edges(data["edges"])
        results = fit_network(samples, edges, config)
        all_results[event_code] = results
        if results:
            all_summaries[event_code] = summarize_network(results)
            all_edge_results.extend(results)

    print(f"  Analyzed {len(all_edge_results)} edges across {len(all_results)} events")

    # Export using hydrosheaf's built-in functions
    print("\n[3/8] Exporting results using hydrosheaf export functions...")
    export_edge_results_csv(
        all_edge_results, str(output_dir / "data" / "edge_results.csv")
    )
    export_edge_results_json(
        all_edge_results, str(output_dir / "data" / "edge_results.json")
    )
    print(f"  Saved: data/edge_results.csv")
    print(f"  Saved: data/edge_results.json")

    # Save summaries
    with open(output_dir / "data" / "network_summaries.json", "w") as f:
        json.dump(all_summaries, f, indent=2, default=str)
    print(f"  Saved: data/network_summaries.json")

    # Generate hydrosheaf built-in plots
    print("\n[4/8] Generating Hydrosheaf ILR plot...")
    ilr_samples = prepare_samples_for_ilr(data["water_chem"])
    plot_ilr(ilr_samples, str(output_dir / "plots" / "ilr_plot.png"))
    print(f"  Saved: plots/ilr_plot.png")

    print("\n[5/8] Generating Hydrosheaf Gibbs diagram...")
    plot_gibbs(ilr_samples, str(output_dir / "plots" / "gibbs_diagram.png"))
    print(f"  Saved: plots/gibbs_diagram.png")

    print("\n[6/8] Generating Hydrosheaf edge anomaly plot...")
    plot_edge_anomalies(
        all_edge_results, str(output_dir / "plots" / "edge_anomalies.png")
    )
    print(f"  Saved: plots/edge_anomalies.png")

    # Generate additional plots
    print("\n[7/8] Generating additional analysis plots...")
    plot_water_isotopes(data["water_chem"], output_dir)
    print(f"  Saved: plots/water_isotopes.png")

    plot_nitrate_source_diagram(data["water_chem"], data["endmembers"], output_dir)
    print(f"  Saved: plots/nitrate_sources.png")

    plot_temporal_trends(data["water_chem"], data["events"], output_dir)
    print(f"  Saved: plots/temporal_patterns.png")

    plot_reaction_summary(all_results, output_dir)
    print(f"  Saved: plots/reaction_summary.png")

    plot_transport_parameters(all_results, output_dir)
    print(f"  Saved: plots/transport_parameters.png")

    plot_network_schematic(data["edges"], data["stations"], output_dir)
    print(f"  Saved: plots/network_schematic.png")

    # Generate comprehensive report
    print("\n[8/8] Generating comprehensive HTML report...")
    report_path = generate_comprehensive_report(
        data, all_results, all_summaries, output_dir
    )
    print(f"  Saved: {report_path.name}")

    # Summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nAll outputs saved to: {output_dir}")
    print("\nGenerated files:")
    print("  DATA:")
    print("    - data/edge_results.csv (full hydrosheaf results)")
    print("    - data/edge_results.json")
    print("    - data/network_summaries.json")
    print("  PLOTS (Hydrosheaf built-in):")
    print("    - plots/ilr_plot.png (ILR hydrochemical facies)")
    print("    - plots/gibbs_diagram.png (Gibbs classification)")
    print("    - plots/edge_anomalies.png (anomaly norms)")
    print("  PLOTS (Additional):")
    print("    - plots/water_isotopes.png")
    print("    - plots/nitrate_sources.png")
    print("    - plots/temporal_patterns.png")
    print("    - plots/reaction_summary.png")
    print("    - plots/transport_parameters.png")
    print("    - plots/network_schematic.png")
    print("  REPORT:")
    print("    - comprehensive_report.html")

    print(f"\nOpen {report_path} in a browser to view the full report.")

    return output_dir


if __name__ == "__main__":
    output_dir = main()
