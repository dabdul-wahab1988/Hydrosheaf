"""Optional plotting helpers."""

from typing import List, Optional, Dict, Any, Union
import numpy as np

from ..inference.edge_fit import EdgeResult
from ..data.units import MOLAR_MASS_G_MOL, CHARGE_EQUIV


def _ensure_dataframe(data: Union[List[Dict[str, Any]], Any]) -> Any:
    """Ensure input is a pandas DataFrame."""
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError("pandas is required for plotting.") from exc

    if isinstance(data, pd.DataFrame):
        return data
    return pd.DataFrame(data)


def _calculate_meq(df: Any) -> Any:
    """Calculate meq/L for major ions from mmol/L or mg/L input.

    Assumes input columns are mmol/L (Hydrosheaf default).
    If columns like 'Ca_meq' already exist, uses them.
    """
    df = df.copy()

    # If mmol/L columns exist (standard Hydrosheaf), convert to meq/L
    # Charge = valence
    for ion, charge in CHARGE_EQUIV.items():
        col_mmol = ion  # e.g. 'Ca'
        col_meq = f"{ion}_meq"

        if col_meq not in df.columns:
            if col_mmol in df.columns:
                # mmol/L -> meq/L: value * charge
                df[col_meq] = df[col_mmol] * charge
            else:
                # Try mg/L fallback if mmol not found but known mg/L pattern exists
                col_mg = f"{ion}_mg_L"
                if col_mg in df.columns:
                    # mg/L -> mmol/L -> meq/L
                    mmol = df[col_mg] / MOLAR_MASS_G_MOL[ion]
                    df[col_meq] = mmol * charge
                else:
                    df[col_meq] = 0.0  # Default to 0 if missing

    # Also ensure TDS exists
    if "TDS" not in df.columns:
        # Approximate TDS if missing: sum of major ions in mg/L
        tds = 0.0
        for ion, mm in MOLAR_MASS_G_MOL.items():
            if ion in df.columns:
                tds += df[ion] * mm
        df["TDS"] = tds

    return df


def classify_water_type(
    row: Any, cation_threshold: float = 0.5, anion_threshold: float = 0.5
) -> str:
    """Classify water type based on meq/L concentrations."""
    # Cations
    ca = row.get("Ca_meq", 0)
    mg = row.get("Mg_meq", 0)
    na = row.get("Na_meq", 0)
    k = row.get("K_meq", 0)
    total_cations = ca + mg + na + k

    if total_cations == 0:
        return "Unknown"

    cations = []
    if ca / total_cations > cation_threshold:
        cations.append("Ca")
    elif ca / total_cations > 0.25:
        cations.append("Ca")

    if mg / total_cations > cation_threshold:
        cations.append("Mg")
    elif mg / total_cations > 0.25:
        cations.append("Mg")

    if na / total_cations > cation_threshold:
        cations.append("Na")
    elif na / total_cations > 0.25:
        cations.append("Na")

    if k / total_cations > cation_threshold:
        cations.append("K")
    elif k / total_cations > 0.25:
        cations.append("K")

    # Anions
    hco3 = row.get("HCO3_meq", 0)
    cl = row.get("Cl_meq", 0)
    so4 = row.get("SO4_meq", 0)
    total_anions = hco3 + cl + so4

    if total_anions == 0:
        return "Unknown"

    anions = []
    if hco3 / total_anions > anion_threshold:
        anions.append("HCO3")
    elif hco3 / total_anions > 0.25:
        anions.append("HCO3")

    if cl / total_anions > anion_threshold:
        anions.append("Cl")
    elif cl / total_anions > 0.25:
        anions.append("Cl")

    if so4 / total_anions > anion_threshold:
        anions.append("SO4")
    elif so4 / total_anions > 0.25:
        anions.append("SO4")

    if not cations or not anions:
        return "Mixed"

    return "-".join(cations + anions)


def plot_edge_anomalies(results: List[EdgeResult], path: Optional[str] = None) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required for plotting.") from exc

    edge_ids = [result.edge_id for result in results]
    values = [result.anomaly_norm for result in results]

    fig, ax = plt.subplots(figsize=(max(6.0, len(edge_ids) * 0.6), 4.0))
    ax.bar(edge_ids, values, color="#4C72B0")
    ax.set_xlabel("Edge")
    ax.set_ylabel("Anomaly norm")
    ax.set_title("Edge anomaly norms")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()

    if path:
        fig.savefig(path, dpi=150)
        plt.close(fig)
    else:
        plt.show()


def plot_gibbs(
    data: Union[List[Dict[str, Any]], Any], path: Optional[str] = None
) -> None:
    """Create Gibbs plots (TDS vs Na/(Na+Ca) and Cl/(Cl+HCO3))."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required for plotting.") from exc

    df = _ensure_dataframe(data)
    df_meq = _calculate_meq(df)

    tds = df_meq["TDS"]
    # Avoid division by zero
    na_ratio = df_meq["Na_meq"] / (df_meq["Na_meq"] + df_meq["Ca_meq"] + 1e-9)
    cl_ratio = df_meq["Cl_meq"] / (df_meq["Cl_meq"] + df_meq["HCO3_meq"] + 1e-9)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Gibbs Plot - Cations
    ax1.scatter(na_ratio, tds, c="steelblue", edgecolor="black", s=80, alpha=0.7)
    ax1.set_yscale("log")
    ax1.set_ylim(1, 100000)
    ax1.set_xlim(-0.05, 1.05)
    ax1.set_xlabel("Na$^{+}$ / (Na$^{+}$ + Ca$^{2+}$)", fontsize=12)
    ax1.set_ylabel("TDS (mg/L)", fontsize=12)
    ax1.set_title("Gibbs Plot - Cations", fontsize=12)
    ax1.grid(True, which="both", ls="--", alpha=0.5)

    # Annotations
    ax1.text(0.9, 20, "Precipitation\nDominance", fontsize=10, ha="center", va="center")
    ax1.text(
        0.25, 500, "Rock Weathering\nDominance", fontsize=10, ha="center", va="center"
    )
    ax1.text(0.9, 5000, "Evaporation\nDominance", fontsize=10, ha="center", va="top")
    ax1.text(
        0.02,
        0.98,
        "(a)",
        transform=ax1.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )

    # Gibbs Plot - Anions
    ax2.scatter(cl_ratio, tds, c="steelblue", edgecolor="black", s=80, alpha=0.7)
    ax2.set_yscale("log")
    ax2.set_ylim(1, 100000)
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_xlabel("Cl$^{-}$ / (Cl$^{-}$ + HCO$_3^{-}$)", fontsize=12)
    ax2.set_ylabel("TDS (mg/L)", fontsize=12)
    ax2.set_title("Gibbs Plot - Anions", fontsize=12)
    ax2.grid(True, which="both", ls="--", alpha=0.5)

    # Annotations
    ax2.text(0.9, 20, "Precipitation\nDominance", fontsize=10, ha="center", va="center")
    ax2.text(
        0.25, 500, "Rock Weathering\nDominance", fontsize=10, ha="center", va="center"
    )
    ax2.text(0.9, 5000, "Evaporation\nDominance", fontsize=10, ha="center", va="top")
    ax2.text(
        0.02,
        0.98,
        "(b)",
        transform=ax2.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
    )

    plt.tight_layout()
    if path:
        plt.savefig(path, dpi=300)
        plt.close(fig)
    else:
        plt.show()


def plot_ilr(
    data: Union[List[Dict[str, Any]], Any], path: Optional[str] = None
) -> None:
    """Create Isometric Log-Ratio (ILR) plot for hydrochemical facies."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("matplotlib is required for plotting.") from exc

    df = _ensure_dataframe(data)
    df_meq = _calculate_meq(df)

    # Classify water types for coloring
    df_meq["water_type"] = df_meq.apply(classify_water_type, axis=1)
    unique_types = df_meq["water_type"].unique()
    colors = plt.cm.get_cmap("Paired")(np.linspace(0, 1, len(unique_types)))
    color_map = dict(zip(unique_types, colors))

    # Calculate ILR coordinates
    # Protect against log(0)
    eps = 1e-9
    ca = df_meq["Ca_meq"] + eps
    mg = df_meq["Mg_meq"] + eps
    na_k = df_meq["Na_meq"] + df_meq["K_meq"] + eps
    cl = df_meq["Cl_meq"] + eps
    so4 = df_meq["SO4_meq"] + eps
    hco3 = df_meq["HCO3_meq"] + eps

    z1 = np.sqrt(2 / 3) * np.log(np.sqrt(ca * mg) / na_k)
    z2 = 1 / np.sqrt(2) * np.log(ca / mg)
    z3 = np.sqrt(2 / 3) * np.log(np.sqrt(cl * so4) / hco3)
    z4 = 1 / np.sqrt(2) * np.log(cl / so4)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(9, 7))
    point_colors = [color_map[wt] for wt in df_meq["water_type"]]

    # Helper for arrow annotation
    def add_arrow(ax, x, y, dx, dy, text, rotation=0):
        ax.arrow(
            x,
            y,
            dx,
            dy,
            head_width=0.2,
            head_length=0.3,
            fc="black",
            ec="black",
            length_includes_head=True,
        )
        ax.text(
            x + dx * 2.5,
            y + dy * 2.5,
            text,
            fontsize=8,
            ha="center",
            va="center",
            rotation=rotation,
        )

    # Upper Left: Ca/Mg vs Cl/SO4
    ax1.scatter(z2, z4, c=point_colors, edgecolors="black", s=50, linewidth=0.5)
    ax1.axhline(0, ls="--", color="gray", alpha=0.6)
    ax1.axvline(0, ls="--", color="gray", alpha=0.6)
    ax1.set_xlim(-10, 10)
    ax1.set_ylim(-10, 10)
    ax1.set_ylabel("[SO$_4^{2-}$|Cl$^-$]", fontsize=12)
    ax1.set_title("[Ca$^{2+}$|Mg$^{2+}$]", fontsize=12, pad=20)
    add_arrow(ax1, -0.5, 8.5, -1, 0, "Mg$^{2+}$ > Ca$^{2+}$")
    add_arrow(ax1, 0.5, 8.5, 1, 0, "Ca$^{2+}$ > Mg$^{2+}$")
    add_arrow(ax1, -8.5, -1, 0, -1, "SO$_4^{2-}$ > Cl$^-$", 90)
    add_arrow(ax1, -8.5, 1, 0, 1, "Cl$^-$ > SO$_4^{2-}$", 90)

    # Upper Right: Anions
    ax2.scatter(z3, z4, c=point_colors, edgecolors="black", s=50, linewidth=0.5)
    ax2.set_xlim(-10, 10)
    ax2.set_ylim(-10, 10)
    ax2.set_yticks([])
    ax2.set_title("[HCO$_3^-$|Cl$^-$, SO$_4^{2-}$]", fontsize=12, pad=20)
    ax2.set_ylabel("[SO$_4^{2-}$|Cl$^-$]", fontsize=12)
    ax2.yaxis.set_label_position("right")

    # Anion fields
    dum1 = np.full(23, 0.5)
    dum2 = np.logspace(-10, np.log10(0.5), 23)
    dum3 = 1 - dum1 - dum2
    # SO4 type
    ax2.plot(
        np.sqrt(2 / 3) * np.log(np.sqrt(dum1 * dum2) / dum3),
        np.sqrt(1 / 2) * np.log(dum2 / dum1),
        ls="--",
        color="gray",
        alpha=0.6,
    )
    ax2.text(1, -1.5, "SO$_4^{2-}$ type", fontsize=9)
    # HCO3 type
    ax2.plot(
        np.sqrt(2 / 3) * np.log(np.sqrt(dum3 * dum2) / dum1),
        np.sqrt(1 / 2) * np.log(dum3 / dum2),
        ls="--",
        color="gray",
        alpha=0.6,
    )
    ax2.text(-6, -4, "HCO$_3^{-}$ type", fontsize=9)
    # Cl type
    ax2.plot(
        np.sqrt(2 / 3) * np.log(np.sqrt(dum3 * dum1) / dum2),
        np.sqrt(1 / 2) * np.log(dum1 / dum3),
        ls="--",
        color="gray",
        alpha=0.6,
    )
    ax2.text(1, 3, "Cl$^{-}$ type", fontsize=9)

    # Lower Left: Cations
    ax3.scatter(z2, z1, c=point_colors, edgecolors="black", s=50, linewidth=0.5)
    ax3.set_xlim(-10, 10)
    ax3.set_ylim(-10, 10)
    ax3.set_xlabel("[Ca$^{2+}$|Mg$^{2+}$]", fontsize=12)
    ax3.set_ylabel("[Ca$^{2+}$, Mg$^{2+}$|Na$^+$ + K$^+$]", fontsize=12)

    # Cation fields (reusing dummy vars logic structure)
    # Mg type
    ax3.plot(
        np.sqrt(1 / 2) * np.log(dum2 / dum1),
        np.sqrt(2 / 3) * np.log(np.sqrt(dum1 * dum2) / dum3),
        ls="--",
        color="gray",
        alpha=0.6,
    )
    ax3.text(-4, 4, "Mg$^{2+}$ type", fontsize=9)
    # Na+K type
    ax3.plot(
        np.sqrt(1 / 2) * np.log(dum3 / dum2),
        np.sqrt(2 / 3) * np.log(np.sqrt(dum3 * dum2) / dum1),
        ls="--",
        color="gray",
        alpha=0.6,
    )
    ax3.text(2, -5, "Na$^+$+K$^+$ type", fontsize=9)
    # Ca type
    ax3.plot(
        np.sqrt(1 / 2) * np.log(dum1 / dum3),
        np.sqrt(2 / 3) * np.log(np.sqrt(dum3 * dum1) / dum2),
        ls="--",
        color="gray",
        alpha=0.6,
    )
    ax3.text(2, 0.5, "Ca$^{2+}$ type", fontsize=9)

    # Lower Right
    ax4.scatter(z3, z1, c=point_colors, edgecolors="black", s=50, linewidth=0.5)
    ax4.axhline(0, ls="--", color="gray", alpha=0.6)
    ax4.axvline(0, ls="--", color="gray", alpha=0.6)
    ax4.set_xlim(-10, 10)
    ax4.set_ylim(-10, 10)
    ax4.set_yticks([])
    ax4.set_xlabel("[HCO$_3^-$|Cl$^-$, SO$_4^{2-}$]", fontsize=12)
    ax4.set_ylabel("[Ca$^{2+}$, Mg$^{2+}$|Na$^+$ + K$^+$]", fontsize=12)
    ax4.yaxis.set_label_position("right")

    add_arrow(ax4, -1, 8.5, -1, 0, "HCO$_3^-$ > Cl$^-$, SO$_4^{2-}$")
    add_arrow(ax4, 1, 8.5, 1, 0, "Cl$^-$, SO$_4^{2-}$ > HCO$_3^-$")
    add_arrow(ax4, -8.5, -1, 0, -1, "Ca, Mg > Na+K", 90)
    add_arrow(ax4, -8.5, 1, 0, 1, "Ca, Mg > Na+K", 90)

    # Legend
    legend_elements = [
        plt.scatter([], [], c=color_map[wt], label=wt, edgecolors="black")
        for wt in unique_types
    ]
    fig.legend(
        handles=legend_elements,
        loc="center left",
        bbox_to_anchor=(0.85, 0.5),
        fontsize=10,
        frameon=True,
    )

    plt.subplots_adjust(wspace=0, hspace=0, right=0.75)

    if path:
        plt.savefig(path, bbox_inches="tight", dpi=300)
        plt.close(fig)
    else:
        plt.show()
