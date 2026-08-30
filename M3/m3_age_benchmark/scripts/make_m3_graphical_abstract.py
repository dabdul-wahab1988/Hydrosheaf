"""Create the accuracy-locked M3 graphical abstract from audited outputs."""

from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "M3" / "M3_geochemistry" / "Graphical Abstract.png"


def panel(ax, xy, width, height, title, lines, color):
    x, y = xy
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.6, edgecolor=color, facecolor="#ffffff", transform=ax.transAxes,
    ))
    ax.add_patch(FancyBboxPatch(
        (x, y + height - 0.075), width, 0.075,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=0, facecolor=color, transform=ax.transAxes,
    ))
    ax.text(x + 0.018, y + height - 0.038, title, color="white",
            fontsize=14, fontweight="bold", va="center", transform=ax.transAxes)
    top = y + height - 0.105
    bottom = y + 0.055
    step = (top - bottom) / max(1, len(lines) - 1)
    for idx, (lead, body) in enumerate(lines):
        yy = top - idx * step
        ax.text(x + 0.024, yy, lead, fontsize=10.2, fontweight="bold",
                color="#111827", va="top", transform=ax.transAxes)
        ax.text(x + 0.14, yy, fill(body, width=40), fontsize=9.5, color="#374151",
                va="top", transform=ax.transAxes, linespacing=1.25)


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(figsize=(16, 9), dpi=200)
    fig.patch.set_facecolor("#eef3f5")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.04, 0.955, "M3 — Conditional value of graph priors in groundwater-age inference",
            fontsize=23, fontweight="bold", color="#102a43", va="top", transform=ax.transAxes)
    ax.text(0.04, 0.905,
            "Audited public-data benchmark • model-derived reference agreement • leakage-guarded tracer withholding",
            fontsize=13.5, color="#486581", va="top", transform=ax.transAxes)

    panel(ax, (0.04, 0.50), 0.43, 0.34, "1  Evidence and reporting guard", [
        ("Input", "1,272 harmonised USGS public-supply rows"),
        ("Guard", "only identifiable, reportable fits are interpreted"),
        ("Support", "329 strict-configuration fits; 289 fraction-sensitivity fits"),
        ("Target", "published LPM products are model-derived references, not true ages"),
    ], "#236b8e")
    panel(ax, (0.53, 0.50), 0.43, 0.34, "2  Reference-workflow agreement", [
        ("Strict", "median |log₁₀ discrepancy| 0.0279; 87.5% within ×2"),
        ("Fractions", "0.0216; 91.7% within ×2 (sensitivity only)"),
        ("Reason", "reported fractions and ages share USGS LPM provenance"),
        ("Withdrawn", "hierarchical old-water prior; unsupported gas audit"),
    ], "#2a7f62")
    panel(ax, (0.04, 0.105), 0.43, 0.31, "3  Graph-age benchmark (N = 329)", [
        ("Best ΔRMSE", "−0.00135 log₁₀, but within-×2 fell 87.5% → 86.9%"),
        ("Decision", "no tested candidate graph gave a robust improvement"),
        ("Controls", "wrong direction +0.1008; randomised +0.655 log₁₀ RMSE"),
    ], "#b45f06")
    panel(ax, (0.53, 0.105), 0.43, 0.31, "4  Target-withheld prediction", [
        ("³H", "20.818 TU baseline; best graph 20.818 TU — negligible"),
        ("SF₆", "candidate graphs increased RMSE by 0.5–3.0%"),
        ("¹⁴C", "hydraulic result differed by −0.003%; not meaningful"),
        ("CFCs", "4/6 reportable fits and zero eligible graph edges: non-estimable"),
    ], "#7c3f83")
    for y in (0.67, 0.265):
        ax.add_patch(FancyArrowPatch((0.475, y), (0.525, y), arrowstyle="-|>",
                                    mutation_scale=18, linewidth=1.5, color="#829ab1",
                                    transform=ax.transAxes))
    ax.text(0.5, 0.054,
            "Conclusion: Hydrosheaf M3 is a controlled diagnostic and falsification benchmark.\n"
            "It is not field validation and does not establish true groundwater ages or universal graph benefit.",
            ha="center", va="center", fontsize=10.4, fontweight="bold", color="#102a43",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.42", facecolor="#d9eaf0", edgecolor="#236b8e", linewidth=1.4))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
