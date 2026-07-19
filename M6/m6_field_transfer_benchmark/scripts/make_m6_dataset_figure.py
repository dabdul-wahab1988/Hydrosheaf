"""Create the M6 dataset/evidence-ceiling figure with a sourced Ghana basemap."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK = REPO_ROOT / "M6" / "m6_field_transfer_benchmark"
RESULTS = BENCHMARK / "results"
BOUNDARY = BENCHMARK / "data" / "ghana_adm0_natural_earth_110m.geojson"
OUT = BENCHMARK / "figures" / "r_publication"

BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#D55E00"
GRAY = "#6B6B6B"
PALE = "#F2F2F2"
DATASET_COLORS = {
    "northern_ghana": BLUE,
    "manu": GREEN,
    "talensi": ORANGE,
}
DATASET_LABELS = {
    "northern_ghana": "Northern Ghana",
    "manu": "Lower Anayari",
    "talensi": "Talensi",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.5,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "grid.color": "#E6E6E6",
            "grid.linewidth": 0.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(ax: plt.Axes, label: str, *, x: float = -0.12, y: float = 1.07) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )


def main() -> None:
    configure_style()
    coordinates = pd.read_csv(RESULTS / "m6_map_coordinates.csv")
    ladder = pd.read_csv(RESULTS / "m6_tier_ladder.csv")
    availability = pd.read_csv(RESULTS / "m6_variable_availability.csv")
    boundary = json.loads(BOUNDARY.read_text(encoding="utf-8"))
    outline = np.asarray(
        boundary["features"][0]["geometry"]["coordinates"][0], dtype=float
    )

    fig = plt.figure(figsize=(7.08, 5.35))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.35, 0.85], hspace=0.46, wspace=0.42)

    ax = fig.add_subplot(grid[0, 0])
    ax.add_patch(
        Polygon(
            outline,
            closed=True,
            facecolor=PALE,
            edgecolor=GRAY,
            linewidth=0.7,
            zorder=0,
        )
    )
    for dataset in ["northern_ghana", "manu", "talensi"]:
        subset = coordinates[coordinates["dataset"] == dataset]
        ax.scatter(
            subset["longitude"],
            subset["latitude"],
            s=10,
            color=DATASET_COLORS[dataset],
            alpha=0.72,
            edgecolors="white",
            linewidths=0.18,
            zorder=2,
        )
    ax.set_xlim(-3.45, 1.30)
    ax.set_ylim(4.55, 11.35)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title("Sampling locations (masked coordinates)")
    ax.grid(lw=0.4)
    ax.legend(
        handles=[
            Line2D(
                [],
                [],
                marker="o",
                ls="",
                color=DATASET_COLORS[key],
                label=DATASET_LABELS[key],
                markersize=4,
            )
            for key in ["northern_ghana", "manu", "talensi"]
        ],
        frameon=False,
        loc="lower left",
    )
    panel_label(ax, "a", x=-0.12, y=1.16)

    ax = fig.add_subplot(grid[0, 1])
    datasets = ["northern_ghana", "manu", "talensi"]
    tiers = [f"Tier {index}" for index in range(5)]
    attained = (
        ladder.pivot(index="dataset", columns="tier", values="attained")
        .loc[datasets, tiers]
        .to_numpy(float)
    )
    tier_cmap = LinearSegmentedColormap.from_list("tier", [PALE, BLUE])
    ax.imshow(attained, cmap=tier_cmap, vmin=0, vmax=1, aspect="auto")
    for row in range(attained.shape[0]):
        for column in range(attained.shape[1]):
            value = attained[row, column]
            ax.text(
                column,
                row,
                "✓" if value else "—",
                ha="center",
                va="center",
                color="white" if value else GRAY,
                fontsize=9,
                fontweight="bold",
            )
    ax.set_xticks(np.arange(5), [f"T{index}" for index in range(5)])
    ax.set_yticks(
        np.arange(3),
        [DATASET_LABELS[dataset] for dataset in datasets],
    )
    ax.set_title("Evidence-tier attainment")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.0,
        -0.18,
        "T4 = maximum M6 chemistry/metadata tier,\nnot complete age–head–screen evidence.",
        transform=ax.transAxes,
        fontsize=5.8,
        color=GRAY,
        va="top",
    )
    panel_label(ax, "b")

    ax = fig.add_subplot(grid[1, :])
    variable_order = [
        "Ca",
        "Mg",
        "Na",
        "K",
        "HCO3",
        "Cl",
        "SO4",
        "NO3",
        "F",
        "Fe",
        "d18O",
        "d2H",
        "Sr_mgL",
        "SiO2_mgL",
        "Calcite_SI",
        "Aquifer_Type",
    ]
    variable_labels = [
        "Ca",
        "Mg",
        "Na",
        "K",
        "HCO₃",
        "Cl",
        "SO₄",
        "NO₃",
        "F",
        "Fe",
        "δ¹⁸O",
        "δ²H",
        "Sr",
        "SiO₂",
        "SI",
        "Aquifer",
    ]
    present = (
        availability.pivot(index="dataset", columns="variable", values="frac_present")
        .loc[datasets, variable_order]
        .to_numpy(float)
    )
    availability_cmap = LinearSegmentedColormap.from_list(
        "availability", ["#F7F7F7", BLUE]
    )
    image = ax.imshow(present, cmap=availability_cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(
        np.arange(len(variable_order)), variable_labels, rotation=38, ha="right"
    )
    ax.set_yticks(
        np.arange(3),
        [DATASET_LABELS[dataset] for dataset in datasets],
    )
    ax.set_title("Variable availability")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(
        image, ax=ax, orientation="horizontal", pad=0.27, fraction=0.10
    )
    colorbar.set_ticks([0, 0.5, 1])
    colorbar.set_ticklabels(["0%", "50%", "100%"])
    colorbar.set_label("Fraction present")
    panel_label(ax, "c")

    fig.subplots_adjust(left=0.14, right=0.99, top=0.95, bottom=0.13)
    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / "figure1_dataset_tier_design"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=600)
    fig.savefig(
        stem.with_suffix(".tif"),
        bbox_inches="tight",
        dpi=300,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    print(f"M6 dataset/evidence figure -> {stem}")


if __name__ == "__main__":
    main()
