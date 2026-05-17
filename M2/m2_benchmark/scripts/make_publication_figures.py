from __future__ import annotations

import math
import uuid
from pathlib import Path
from textwrap import fill

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import yaml
from matplotlib.ticker import MaxNLocator


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
M6_BENCHMARK_ROOT = Path(__file__).resolve().parents[3] / "M6" / "m6_robustness_benchmark"
# Direct output to Manuscript folder
FIGURE_DIR = BENCHMARK_ROOT / "figures" / "Manuscript_Ready"
RESULT_DIR = BENCHMARK_ROOT / "results"
EXTERNAL_DIR = BENCHMARK_ROOT / "external"
M6_RESULT_DIR = M6_BENCHMARK_ROOT / "results"


# --- Typography Standards (Matching Fig 5 style) ---
FONT_TITLE = 18
FONT_LABEL = 16
FONT_TICK = 16
FONT_LEGEND = 14
FONT_ANNOTATE = 11

COLORS = {
    "young": "#3b82f6",
    "intermediate": "#22c55e",
    "old": "#f59e0b",
    "fossil": "#ef4444",
    "mixed": "#8b5cf6",
    "line": "#334155",
    "muted": "#64748b",
}


PROCESS_LABEL_MAP = {
    "calcite": "Calcite Dissolution",
    "calcite_open": "Calcite Dissolution (Open)",
    "calcite_closed": "Calcite Dissolution (Closed)",
    "dolomite": "Dolomite Dissolution",
    "dolomite_open": "Dolomite Dissolution (Open)",
    "dolomite_closed": "Dolomite Dissolution (Closed)",
    "gypsum": "Gypsum Dissolution",
    "albite": "Plagioclase Weathering",
    "anorthite": "Plagioclase Weathering",
    "k_feldspar": "K-Feldspar Weathering",
    "microcline": "K-Feldspar Weathering",
    "biotite": "Ferromagnesian Weathering",
    "chlorite": "Ferromagnesian Weathering",
    "enstatite": "Ferromagnesian Weathering",
    "diopside": "Ferromagnesian Weathering",
    "halite": "Halite Dissolution",
    "fluorite": "Fluorite Dissolution",
    "NO3src": "Nitrate Input",
    "denit": "Denitrification",
    "CaNa_exch": "Ca-Na Ion Exchange",
    "MgNa_exch": "Mg-Na Ion Exchange",
    "CaK_exch": "Ca-K Ion Exchange",
    "NaK_exch": "Na-K Ion Exchange",
    "pyrite_oxidation_aerobic": "Pyrite Oxidation (Aerobic)",
    "pyrite_oxidation_denit": "Pyrite Oxidation (Denit)",
    "pyrite_net": "Pyrite Oxidation (Net)",
    "albite_acid": "Plagioclase Weathering (Acid)",
    "anorthite_acid": "Plagioclase Weathering (Acid)",
    "chlorite_acid": "Ferromagnesian Weathering (Acid)",
    "k_feldspar_acid": "K-Feldspar Weathering (Acid)",
    "biotite_acid": "Ferromagnesian Weathering (Acid)",
    "sulfate_reduction": "Microbial Sulfate Reduction",
    "nitrification": "Agricultural Nitrification",
    "potash_fertilizer": "Potash Fertilizer (KCl)",
    "road_salt_calcium": "De-Icing Salt (CaCl2)",
    "road_salt_magnesium": "De-Icing Salt (MgCl2)",
    "iron_reduction": "Microbial Iron Reduction",
    "Conservative": "Evaporative Concentration",
}


def _save(fig: plt.Figure, name: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    target = FIGURE_DIR / name
    temp = target.with_name(f"{target.stem}.{uuid.uuid4().hex}.tmp{target.suffix}")
    fig.savefig(temp, dpi=300, bbox_inches="tight", format=target.suffix.lstrip("."))
    temp.replace(target)
    plt.close(fig)


def _panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=FONT_LABEL,
        fontweight="bold",
        va="top",
        ha="left",
    )


def _wrap_labels(labels: list[str], width: int = 14) -> list[str]:
    return [fill(str(label), width=width) for label in labels]


def load_truth() -> dict:
    with (BENCHMARK_ROOT / "config" / "ground_truth.yaml").open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def plot_manuscript_fig1_architecture() -> None:
    """Figure 1: Export the Mermaid architecture diagram."""
    mermaid_path = Path(__file__).resolve().parents[3] / "main_Figure1.txt"
    out_path = FIGURE_DIR / "Manuscript_Fig1_Architecture.mermaid"
    if mermaid_path.exists():
        FIGURE_DIR.mkdir(parents=True, exist_ok=True)
        with mermaid_path.open("r") as f:
            content = f.read()
        with out_path.open("w") as f:
            f.write(content)
        print(f"Exported Mermaid diagram to {out_path.name}")


def plot_manuscript_fig2_topology_validation() -> None:
    """Figure 2: Physical Prior & Topology Validation (MODPATH agreement)."""
    modpath_path = EXTERNAL_DIR / "modpath" / "results" / "modpath_topology_summary.csv"
    if not modpath_path.exists():
        print("Skipping Figure 2: MODPATH results not found.")
        return

    modpath = pd.read_csv(modpath_path).iloc[0]
    fig, ax = plt.subplots(figsize=(8, 7.5))

    labels = ["True\nPositive", "False\nPositive", "False\nNegative"]
    values = [modpath["true_positive_edges"], modpath["false_positive_edges"], modpath["false_negative_edges"]]
    bars = ax.bar(labels, values, color=["#2563eb", "#ef4444", "#f97316"], alpha=0.85, width=0.6)

    ax.set_title("Topology Validation vs. MODPATH Reference", fontsize=FONT_TITLE, fontweight="bold", pad=20)
    ax.set_ylabel("Count of Directed Edges", fontsize=FONT_LABEL, fontweight="bold")

    ax.tick_params(axis='both', which='major', labelsize=FONT_TICK)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))

    stats_text = (f"F1-Score: {modpath['edge_f1']:.2f}\n"
                  f"Precision: {modpath['edge_precision']:.2f}\n"
                  f"Recall: {modpath['edge_recall']:.2f}")
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, ha="right", va="top",
            fontsize=FONT_ANNOTATE, fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="#d1d5db", boxstyle="round,pad=0.5"))

    ax.text(0.5, -0.25, "Comparison of Hydrosheaf inferred topology against physical flow simulation.",
            transform=ax.transAxes, ha="center", va="center", fontsize=FONT_ANNOTATE, style='italic', color="#475569")

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{int(height)}', ha='center', va='bottom', fontsize=FONT_ANNOTATE, fontweight="bold")

    fig.tight_layout()
    _save(fig, "Manuscript_Fig2_Topology_Validation.png")


def plot_manuscript_fig3_synthetic_validation() -> None:
    """Figure 3: Comprehensive Synthetic Validation (The 'Master Proof').
    Consolidates transport recovery, reaction sparsity, and sensitivity analysis.
    """
    import seaborn as sns
    from matplotlib.gridspec import GridSpec

    edge_path = RESULT_DIR / "edge_fit_results.csv"
    truth = load_truth()

    if not edge_path.exists():
        print("Skipping Figure 3: Synthetic results not found.")
        return

    df = pd.read_csv(edge_path)

    # Filter for the 'baseline' scenario
    baseline = df[(df["scenario"] == "complete") & (df["topology_variant"] == "full")].copy()

    fig = plt.figure(figsize=(22, 17))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.2], hspace=0.45, wspace=0.35)

    # --- PANEL A: Transport Parameter Recovery (Violin-Box Hybrid) ---
    ax_a = fig.add_subplot(gs[0, 0])
    trans_data = []
    for edge in truth["generation_edges"]:
        eid = edge["edge_id"]
        sub = baseline[baseline["edge_id"] == eid]

        # Gamma
        true_gamma = 1.0 / edge.get("gamma_inv", 1.0)
        temp = sub[["gamma"]].copy().rename(columns={"gamma": "Value"})
        temp["Parameter"] = r"$\gamma$ (Evap)"
        temp["True"] = true_gamma
        trans_data.append(temp)

        # Mixing (if active)
        if edge.get("f", 0) > 0:
            temp_f = sub[["f"]].copy().rename(columns={"f": "Value"})
            temp_f["Parameter"] = r"$f$ (Mix)"
            temp_f["True"] = edge["f"]
            trans_data.append(temp_f)

    df_a = pd.concat(trans_data).reset_index(drop=True)
    sns.violinplot(data=df_a, x="Parameter", y="Value", ax=ax_a, inner=None, color="#f1f5f9", linewidth=1.5)
    sns.boxplot(data=df_a, x="Parameter", y="Value", ax=ax_a, width=0.15, boxprops={'zorder': 2}, color="white")
    sns.stripplot(data=df_a, x="Parameter", y="Value", ax=ax_a, size=3, alpha=0.3, color="#2563eb", jitter=0.2)

    # Plot true values as horizontal markers
    unique_params = df_a["Parameter"].unique()
    for i, param in enumerate(unique_params):
        true_val = df_a[df_a["Parameter"] == param]["True"].iloc[0]
        ax_a.plot([i - 0.2, i + 0.2], [true_val, true_val], color="red", ls="--", lw=2, label="Ground Truth" if i==0 else "")

    ax_a.set_title("A. Physical Transport Parameter Recovery", fontsize=FONT_TITLE, fontweight="bold", pad=20)
    ax_a.set_ylabel("Parameter Value / Fraction", fontsize=FONT_LABEL, fontweight="bold")
    ax_a.tick_params(labelsize=FONT_TICK)
    ax_a.legend(frameon=False, fontsize=FONT_LEGEND)

    # --- PANEL B: Reaction Extent & Sparsity (JointGrid Style) ---
    # We use a nested GridSpec for the JointGrid look within the main figure
    gs_b = gs[0, 1].subgridspec(4, 4, hspace=0.05, wspace=0.05)
    ax_b_main = fig.add_subplot(gs_b[1:, 0:3])
    ax_b_histx = fig.add_subplot(gs_b[0, 0:3], sharex=ax_b_main)
    ax_b_histy = fig.add_subplot(gs_b[1:, 3], sharey=ax_b_main)

    minerals = ["calcite", "gypsum", "halite"]
    plot_data_b = []
    for m in minerals:
        # Construct truth mapping for each edge explicitly (including false edges as 0)
        truth_map = {e["edge_id"]: e["reactions"].get(m, 0.0) for e in truth["generation_edges"]}
        temp = baseline[["edge_id", f"reaction_{m}"]].copy()
        temp["True"] = temp["edge_id"].apply(lambda x: truth_map.get(x, 0.0))
        temp = temp.rename(columns={f"reaction_{m}": "Inferred"})
        temp["Mineral"] = PROCESS_LABEL_MAP.get(m, m.capitalize())
        plot_data_b.append(temp[["Inferred", "True", "Mineral"]])

    df_b = pd.concat(plot_data_b).reset_index(drop=True)
    sns.scatterplot(data=df_b, x="True", y="Inferred", hue="Mineral", palette="viridis",
                    s=80, alpha=0.6, edgecolor="white", linewidth=0.5, ax=ax_b_main)

    lim = max(df_b["True"].max(), df_b["Inferred"].max(), 0.1) * 1.1
    ax_b_main.plot([-0.05, lim], [-0.05, lim], "k--", alpha=0.3, lw=1.5, zorder=0)

    # Marginals
    sns.kdeplot(data=df_b, x="True", ax=ax_b_histx, fill=True, alpha=0.3, hue="Mineral", legend=False)
    sns.kdeplot(data=df_b, y="Inferred", ax=ax_b_histy, fill=True, alpha=0.3, hue="Mineral", legend=False)

    ax_b_histx.axis("off")
    ax_b_histy.axis("off")
    ax_b_main.set_xlabel("True Extent (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax_b_main.set_ylabel("Inferred Extent (mmol/L)", fontsize=FONT_LABEL, fontweight="bold")
    ax_b_main.tick_params(labelsize=FONT_TICK)
    ax_b_main.legend(frameon=False, fontsize=FONT_LEGEND, loc="upper left")

    # Title on the top axis (histogram) to avoid overlap
    ax_b_histx.set_title("B. Reaction Recovery & Signal Sparsity", fontsize=FONT_TITLE, fontweight="bold", pad=20)

    # --- PANEL C: Missing Data Performance Decay ---
    ax_c = fig.add_subplot(gs[1, 0])
    # Group by scenario
    scen_map = {
        "complete": "1. All Data",
        "ion_incomplete": "2. Missing Major Ions",
        "head_absent": "3. No Head Data"
    }
    df["Scenario_Label"] = df["scenario"].map(scen_map)
    # Ensure order
    df["Scenario_Label"] = pd.Categorical(df["Scenario_Label"], categories=sorted(scen_map.values()), ordered=True)

    sns.boxplot(data=df, x="Scenario_Label", y="chemistry_r2", ax=ax_c, hue="Scenario_Label", palette="coolwarm", width=0.5, legend=False)
    ax_c.set_title("C. Model Robustness Stress-Test", fontsize=FONT_TITLE, fontweight="bold", pad=20)
    ax_c.set_ylabel(r"Discovery Accuracy ($R^2$)", fontsize=FONT_LABEL, fontweight="bold")
    ax_c.set_xlabel("Degradation Tier", fontsize=FONT_LABEL, fontweight="bold")
    ax_c.set_ylim(0, 1.05)
    ax_c.tick_params(labelsize=FONT_TICK)
    ax_c.grid(True, ls=":", alpha=0.3)

    # --- PANEL D: Isotopic Forensic Proof (Forensic Consistency) ---
    ax_d = fig.add_subplot(gs[1, 1])
    iso_path = RESULT_DIR / "res_isotopes.csv"
    if iso_path.exists():
        df_d = pd.read_csv(iso_path)
        # REVIEWER FIX: Plot the actual isotopic SHIFT (Delta d18O) instead of absolute states
        sns.regplot(data=df_d, x="true_shift", y="inf_shift", ax=ax_d, color="#db2777",
                    scatter_kws={'s': 80, 'alpha': 0.6, 'edgecolor': 'white'},
                    line_kws={'color': 'black', 'ls': '--', 'lw': 1.5})

        # 1:1 Line
        d_lim = [min(df_d["true_shift"].min(), df_d["inf_shift"].min()) - 0.2,
                 max(df_d["true_shift"].max(), df_d["inf_shift"].max()) + 0.2]
        ax_d.plot(d_lim, d_lim, "k:", alpha=0.3)

        ax_d.set_title("D. Isotopic Forensic Consistency Check", fontsize=FONT_TITLE, fontweight="bold", pad=20)
        ax_d.set_xlabel(r"True $\Delta \delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
        ax_d.set_ylabel(r"Inferred $\Delta \delta^{18}$O Shift (permil)", fontsize=FONT_LABEL, fontweight="bold")
        ax_d.tick_params(labelsize=FONT_TICK)

        # Add R2 annotation
        r2_iso = np.corrcoef(df_d["true_shift"], df_d["inf_shift"])[0,1]**2
        ax_d.text(0.05, 0.95, fr"$R^2 = {r2_iso:.3f}$", transform=ax_d.transAxes,
                  fontsize=FONT_LEGEND, fontweight="bold", verticalalignment='top',
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#d1d5db"))
    else:
        ax_d.text(0.5, 0.5, "Isotope results not found", ha="center", va="center")
        ax_d.set_title("D. Isotopic Forensic Consistency", fontsize=FONT_TITLE, fontweight="bold")

    plt.suptitle("Information-Theoretic Proof of Hydrosheaf Inversion Accuracy\n(Synthetic Benchmark Suite: 100 Realisations)",
                 fontsize=26, fontweight="bold", y=0.98)

    plt.tight_layout()
    _save(fig, "Manuscript_Fig3_Synthetic_Validation.png")


def plot_manuscript_fig4_ghana_network() -> None:
    """Figure 4: The 'Ghana Discovery' 4-Panel Spatial Process Network."""
    path = RESULT_DIR / "field_discovery_results.csv"
    if not path.exists():
        print("Skipping Figure 4: Field discovery results not found.")
        return

    field = pd.read_csv(path)

    manu_df = pd.read_csv(BENCHMARK_ROOT.parents[1] / "data" / "LowerAnayari" / "manu.csv")
    talensi_df = pd.read_csv(BENCHMARK_ROOT.parents[1] / "data" / "Talensi_MiningArea" / "talensi.csv")

    manu_pos = {str(row["Sample ID"]): (float(row["X coordinate"]), float(row["Y coordinate"])) for _, row in manu_df.iterrows() if pd.notna(row["X coordinate"])}
    manu_elev = {str(row["Sample ID"]): float(row["Elevation"]) for _, row in manu_df.iterrows() if pd.notna(row["Elevation"])}

    talensi_pos = {str(row["Code"]): (float(row["Longitude"]), float(row["Latitude"])) for _, row in talensi_df.iterrows() if pd.notna(row["Longitude"])}
    talensi_elev = {str(row["Code"]): float(row["Elevation"]) for _, row in talensi_df.iterrows() if pd.notna(row["Elevation"])}

    FAMILY_COLORS = {
        "Plagioclase Weathering": "#2563eb", "K-Feldspar Weathering": "#7c3aed", "Ferromagnesian Weathering": "#92400e",
        "Carbonate Dissolution": "#16a34a", "Evaporite Dissolution": "#0891b2", "Anthropogenic Input": "#dc2626",
        "Redox / Oxidation": "#db2777", "Evaporative Concentration": "#64748b",
    }

    def get_family(proc):
        p = proc.lower()
        if any(x in p for x in ["albite", "anorthite"]): return "Plagioclase Weathering"
        if any(x in p for x in ["k_feldspar", "microcline"]): return "K-Feldspar Weathering"
        if any(x in p for x in ["biotite", "chlorite", "enstatite", "diopside"]): return "Ferromagnesian Weathering"
        if any(x in p for x in ["calcite", "dolomite", "magnesite", "aragonite"]): return "Carbonate Dissolution"
        if any(x in p for x in ["gypsum", "halite", "sylvite", "fluorite", "apatite"]): return "Evaporite Dissolution"
        if any(x in p for x in ["no3", "potash", "nitrification", "road_salt"]): return "Anthropogenic Input"
        if any(x in p for x in ["denit", "pyrite", "sulfate_reduction", "iron_reduction"]): return "Redox / Oxidation"
        return "Evaporative Concentration"

    def get_dominant_info(row):
        extents = {}
        for k, v in row.items():
            if k.startswith("extent_") and pd.notna(v):
                lbl = k.replace("extent_", "")
                if lbl in ["pyrite_oxidation_aerobic", "pyrite_net"]:
                    extents["pyrite"] = extents.get("pyrite", 0) + abs(v)
                else:
                    extents[lbl] = abs(v)

        if not extents or max(extents.values()) < 0.02: return "Evaporative Concentration", 0.5
        dom_proc = max(extents, key=extents.get)
        return get_family(dom_proc), extents[dom_proc]

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    ((ax_m_all, ax_m_top), (ax_t_all, ax_t_top)) = axes

    import matplotlib.patheffects as path_effects

    # Load Edge-Level PSI Data
    psi_path = RESULT_DIR / "top_edges_psi.csv"
    edge_psi_dict = {}
    if psi_path.exists():
        psi_df = pd.read_csv(psi_path)
        # Helper to harmonize family names from CSV to the new figure legend
        fam_map = {
            "Plagioclase": "Plagioclase Weathering",
            "K-Feldspar": "K-Feldspar Weathering",
            "Ferromagnesian": "Ferromagnesian Weathering",
            "Carbonates": "Carbonate Dissolution",
            "Evaporites": "Evaporite Dissolution",
            "Anthropogenic": "Anthropogenic Input",
            "Redox": "Redox / Oxidation",
            "Conservative": "Evaporative Concentration"
        }
        for _, row in psi_df.iterrows():
            raw_fam = str(row["family"])
            clean_fam = fam_map.get(raw_fam, raw_fam)
            edge_psi_dict[str(row["edge_id"])] = (clean_fam, float(row["psi"]))

    def draw_site_panel(ax, edges_df, pos_dict, elev_dict, title, is_top=False):
        G = nx.DiGraph()
        for _, row in edges_df.iterrows():
            u, v, eid = str(row["u"]), str(row["v"]), str(row["edge_id"])
            if u in pos_dict and v in pos_dict:
                fam, mag = get_dominant_info(row)
                if is_top and eid in edge_psi_dict:
                    fam, psi = edge_psi_dict[eid]
                    # Filter: Only show edges with > 20% probability
                    if psi > 0.20:
                        G.add_edge(u, v, family=fam, magnitude=mag, psi=psi)
                elif not is_top:
                    G.add_edge(u, v, family=fam, magnitude=mag, psi=0.5)

        nodes = list(G.nodes())
        node_colors = [elev_dict.get(n, 0) for n in nodes]

        node_coll = nx.draw_networkx_nodes(
            G, pos_dict, nodelist=nodes, node_color=node_colors,
            cmap="terrain", edgecolors="#1f2937", node_size=60 if is_top else 80,
            ax=ax, alpha=0.9, linewidths=1.0
        )

        for u, v, d in G.edges(data=True):
            if not is_top:
                # Standardized thin width
                width = 1.5
                alpha_val = 0.8
            else:
                # Standardized thin width for PSI panels as well
                psi = d["psi"]
                width = 1.5
                alpha_val = 0.3 + (psi * 0.7)

            nx.draw_networkx_edges(
                G, pos_dict, edgelist=[(u,v)], arrows=True,
                arrowsize=10 if not is_top else int(10 + psi * 5),
                width=width, edge_color=FAMILY_COLORS[d["family"]],
                ax=ax, connectionstyle="arc3,rad=0.1", alpha=alpha_val
            )

            if is_top:
                # Add PSI text label to the edge
                mid_x = (pos_dict[u][0] + pos_dict[v][0]) / 2
                mid_y = (pos_dict[u][1] + pos_dict[v][1]) / 2
                # Fading: Use alpha=0.3 for low-probability labels
                label_alpha = 1.0 if d["psi"] >= 0.50 else 0.3
                ax.annotate(f"{d['psi']*100:.0f}%", (mid_x, mid_y),
                            color=FAMILY_COLORS[d["family"]], fontsize=8, fontweight="bold",
                            alpha=label_alpha,
                            path_effects=[path_effects.withStroke(linewidth=2, foreground="white", alpha=label_alpha)])

        if is_top:
            for node, (x, y) in pos_dict.items():
                if node in G:
                    ax.annotate(node, (x, y), xytext=(6, 6), textcoords='offset points',
                               fontsize=7.5, fontweight="bold",
                               path_effects=[path_effects.withStroke(linewidth=2, foreground="white")])

        ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold", pad=12)
        return node_coll

    m_edges = field[field["edge_id"].str.startswith("Manu")].copy()
    draw_site_panel(ax_m_all, m_edges, manu_pos, manu_elev, "(a) Lower Anayari: All Discovered Paths")
    draw_site_panel(ax_m_top, m_edges, manu_pos, manu_elev, "(b) Lower Anayari: Process-Identification Probability (PSI)", is_top=True)

    t_edges = field[field["edge_id"].str.startswith("Talensi")].copy()
    draw_site_panel(ax_t_all, t_edges, talensi_pos, talensi_elev, "(c) Talensi: All Discovered Paths")
    sm = draw_site_panel(ax_t_top, t_edges, talensi_pos, talensi_elev, "(d) Talensi: Process-Identification Probability (PSI)", is_top=True)

    for ax in axes.flat:
        ax.set_xlabel("Longitude / Easting", fontsize=FONT_LABEL, fontweight="bold")
        ax.set_ylabel("Latitude / Northing", fontsize=FONT_LABEL, fontweight="bold")
        ax.tick_params(axis='both', which='major', labelsize=FONT_TICK, left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.set_axis_on()
        ax.grid(True, ls=":", alpha=0.3)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=c, lw=4, label=k) for k, c in FAMILY_COLORS.items()]
    fig.legend(handles=legend_elements, loc="lower center", ncol=4,
                      title="Geochemical Process Families", bbox_to_anchor=(0.5, -0.05),
                      fontsize=FONT_LEGEND, title_fontsize=FONT_LABEL)

    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("Relative Elevation (m a.s.l.)", fontsize=FONT_LABEL, fontweight="bold")
    cbar.ax.tick_params(labelsize=FONT_TICK)

    plt.suptitle("Spatiotemporal Process Discovery Networks", fontsize=22, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0.02, 0.9, 0.94])
    _save(fig, "Manuscript_Fig4_Ghana_Process_Network.png")


def plot_manuscript_fig5_residence_time_validation() -> None:
    """Figure 5: Synthetic and public residence-time validation."""
    synthetic_path = RESULT_DIR / "age_inference_validation.csv"
    public_path = EXTERNAL_DIR / "usgs_age" / "results" / "usgs_age_validation.csv"
    if not synthetic_path.exists() and not public_path.exists():
        print("Skipping Figure 5: residence-time validation data not found.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(16, 7.5))
    ax_syn, ax_pub = axes

    if synthetic_path.exists():
        syn = pd.read_csv(synthetic_path).dropna(
            subset=["true_mrt_years", "network_bayesian_years"]
        )
        class_colors = {
            "young": COLORS["young"],
            "intermediate": COLORS["intermediate"],
            "old": COLORS["old"],
            "fossil": COLORS["fossil"],
            "mixed": COLORS["mixed"],
        }
        for age_class, color in class_colors.items():
            sub = syn[syn["age_class"] == age_class]
            if sub.empty:
                continue
            ax_syn.scatter(
                sub["true_mrt_years"],
                sub["network_bayesian_years"],
                s=70,
                alpha=0.72,
                edgecolor="white",
                linewidth=0.5,
                color=color,
                label=age_class.capitalize(),
            )
        lim_min = max(0.5, min(syn["true_mrt_years"].min(), syn["network_bayesian_years"].min()) * 0.6)
        lim_max = max(syn["true_mrt_years"].max(), syn["network_bayesian_years"].max()) * 1.6
        lims = np.logspace(np.log10(lim_min), np.log10(lim_max), 100)
        ax_syn.plot(lims, lims, color="#111827", ls="--", lw=1.5)
        ax_syn.fill_between(lims, lims / 3, lims * 3, color="#94a3b8", alpha=0.14)
        log_true = np.log10(np.maximum(syn["true_mrt_years"], 0.1))
        log_inf = np.log10(np.maximum(syn["network_bayesian_years"], 0.1))
        r2 = np.corrcoef(log_true, log_inf)[0, 1] ** 2
        mae = np.mean(np.abs(syn["network_bayesian_years"] - syn["true_mrt_years"]))
        ax_syn.text(
            0.05,
            0.95,
            f"$R^2$ = {r2:.2f}\nMAE = {mae:.1f} y",
            transform=ax_syn.transAxes,
            ha="left",
            va="top",
            fontsize=FONT_ANNOTATE,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d1d5db"),
        )
        ax_syn.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")
    else:
        ax_syn.text(0.5, 0.5, "Synthetic age validation not found", ha="center", va="center")

    if public_path.exists():
        pub = pd.read_csv(public_path).dropna(
            subset=["reference_mean_age_years", "hydrosheaf_age_years"]
        )
        pub = pub[pub["reference_mean_age_years"] <= 50000].copy()
        categories = {
            "Modern tracers": pub["supported_tracers"].astype(str).str.contains("3H|SF6|CFC", regex=True),
            "14C inclusive": pub["supported_tracers"].astype(str).str.contains("14C", regex=False),
        }
        plotted = pd.Series(False, index=pub.index)
        palette = {"Modern tracers": COLORS["young"], "14C inclusive": COLORS["fossil"], "Other": COLORS["muted"]}
        for label, mask in categories.items():
            sub = pub[mask & ~plotted]
            plotted.loc[sub.index] = True
            if not sub.empty:
                ax_pub.scatter(
                    sub["reference_mean_age_years"],
                    sub["hydrosheaf_age_years"],
                    s=32,
                    alpha=0.45,
                    edgecolor="white",
                    linewidth=0.25,
                    color=palette[label],
                    label=label,
                )
        other = pub[~plotted]
        if not other.empty:
            ax_pub.scatter(
                other["reference_mean_age_years"],
                other["hydrosheaf_age_years"],
                s=28,
                alpha=0.35,
                color=palette["Other"],
                label="Other",
            )
        lims = np.logspace(-1, 5, 100)
        ax_pub.plot(lims, lims, color="#111827", ls="--", lw=1.5)
        ax_pub.fill_between(lims, lims / 10, lims * 10, color="#94a3b8", alpha=0.14)
        log_ref = np.log10(np.maximum(pub["reference_mean_age_years"], 0.1))
        log_inf = np.log10(np.maximum(pub["hydrosheaf_age_years"], 0.1))
        r2 = np.corrcoef(log_ref, log_inf)[0, 1] ** 2
        mae_log = np.mean(np.abs(log_ref - log_inf))
        ax_pub.text(
            0.05,
            0.95,
            f"$R^2$ = {r2:.2f}\nMAE = {mae_log:.2f} log units",
            transform=ax_pub.transAxes,
            ha="left",
            va="top",
            fontsize=FONT_ANNOTATE,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d1d5db"),
        )
        ax_pub.legend(frameon=False, fontsize=FONT_LEGEND, loc="lower right")
    else:
        ax_pub.text(0.5, 0.5, "Public age validation not found", ha="center", va="center")

    for ax, title in [
        (ax_syn, "(a) Synthetic Network Bayesian Age Recovery"),
        (ax_pub, "(b) Public Tracer-Age Agreement"),
    ]:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Reference / True Age (years)", fontsize=FONT_LABEL, fontweight="bold")
        ax.set_ylabel("Hydrosheaf Age (years)", fontsize=FONT_LABEL, fontweight="bold")
        ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold", pad=14)
        ax.tick_params(labelsize=FONT_TICK)
        ax.grid(True, which="major", ls=":", alpha=0.25)

    plt.suptitle("Residence-Time Validation across Synthetic and Public Tracer Evidence", fontsize=22, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save(fig, "Manuscript_Fig5_Residence_Time_Validation.png")


def plot_manuscript_fig6_optimal_model_selection() -> None:
    """Figure 6: The Information-Theoretic 'Optimal Discovery' Plateau."""
    path = M6_RESULT_DIR / "m6_regularization_path.csv"
    if not path.exists():
        print("Skipping Figure 6: Regularization path data not found.")
        return

    df_all = pd.read_csv(path)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 13))

    def plot_panel(ax, site_name, label):
        df = df_all[df_all["site"] == site_name].copy()
        ax_right = ax.twinx()
        lns1 = ax.plot(df["lambda"], df["residual_norm"], "o-", color="#3b82f6", label="Residual Norm", markersize=5, linewidth=1.5, alpha=0.8)
        lns2 = ax_right.plot(df["lambda"], df["aicc"], "s-", color="#ef4444", label="AICc", markersize=5, linewidth=1.5, alpha=0.8)
        best_lambda = df.loc[df["aicc"].idxmin(), "lambda"]
        ax.axvline(best_lambda, color="#1f2937", ls="--", alpha=0.4)
        ax_right.annotate(fr"Optimal $\lambda$ = {best_lambda:.4f}", xy=(best_lambda, 0.5), xycoords=("data", "axes fraction"), ha="center", va="center", rotation=90, fontsize=10, fontweight="bold", bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#d1d5db", alpha=0.9))
        ax.set_xscale("log")
        ax.set_xlabel(r"Regularization Strength ($\lambda$)", fontsize=FONT_LABEL, fontweight="bold")
        ax.set_ylabel("Residual Norm", color="#3b82f6", fontsize=FONT_LABEL, fontweight="bold")
        ax_right.set_ylabel("AICc", color="#ef4444", fontsize=FONT_LABEL, fontweight="bold")
        ax.tick_params(axis='both', which='major', labelsize=FONT_TICK)
        ax_right.tick_params(axis='both', which='major', labelsize=FONT_TICK)
        ax.set_title(f"({label}) {'Lower Anayari' if site_name == 'Manu' else 'Talensi Mining Area'}", fontsize=FONT_TITLE, fontweight="bold", pad=12)
        return lns1 + lns2

    h = plot_panel(ax1, "Manu", "a")
    plot_panel(ax2, "Talensi", "b")
    fig.legend(h, ["Residual Norm", "AICc"], loc="lower center", ncol=2, frameon=True, fontsize=FONT_LEGEND, edgecolor="#d1d5db", bbox_to_anchor=(0.5, -0.02))
    plt.suptitle("Optimal Model Selection via AICc Minimum", fontsize=FONT_TITLE, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save(fig, "Manuscript_Fig6_Optimal_Model_Selection.png")


def plot_manuscript_fig7_psi_robustness_guarantee() -> None:
    """Figure 7: The Phase Stability Matrix (PSI) Heatmap."""
    path = M6_RESULT_DIR / "m6_phase_stability_index.csv"
    if not path.exists():
        print("Skipping Figure 7: PSI data not found.")
        return

    df = pd.read_csv(path)
    df["label"] = df["mineral"].map(lambda x: PROCESS_LABEL_MAP.get(x, x))
    pivot_df = df.pivot_table(index="label", columns="site", values="psi", aggfunc="max").fillna(0)
    pivot_df["total"] = pivot_df.sum(axis=1)
    pivot_df = pivot_df.sort_values("total", ascending=False).drop(columns="total")

    import seaborn as sns
    fig, ax = plt.subplots(figsize=(10, 11))
    sns.heatmap(pivot_df * 100, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': 'Identification Probability (%)'}, linewidths=.5, ax=ax, annot_kws={"weight": "bold", "size": 14})

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=FONT_TICK)
    cbar.set_label('Identification Probability (%)', fontsize=FONT_LABEL, fontweight="bold")

    ax.set_title("Phase Stability Matrix (PSI)\nDiscovery Robustness across Geologic Provinces", fontsize=FONT_TITLE, fontweight="bold", pad=20)
    ax.set_xticklabels(["Lower Anayari (Basement)", "Talensi (Mining/Agriculture)"], rotation=0, fontweight="bold", fontsize=FONT_LABEL)
    ax.tick_params(axis='y', labelsize=FONT_TICK)
    ax.set_ylabel("")
    ax.set_xlabel("")

    fig.tight_layout()
    _save(fig, "Manuscript_Fig7_PSI_Robustness_Guarantee.png")


def main() -> None:
    plot_manuscript_fig1_architecture()
    plot_manuscript_fig2_topology_validation()
    plot_manuscript_fig3_synthetic_validation()
    plot_manuscript_fig4_ghana_network()
    plot_manuscript_fig5_residence_time_validation()
    plot_manuscript_fig6_optimal_model_selection()
    plot_manuscript_fig7_psi_robustness_guarantee()
    print("M2 Manuscript-Ready figures generated in 'figures/Manuscript_Ready/'.")


if __name__ == "__main__":
    main()
