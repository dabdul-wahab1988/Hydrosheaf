"""Generate database-backed M5 publication figures.

This script reads only ``results/m5_results.duckdb``. It is intentionally
separate from the analysis workflow: all scientific calculations must already
exist in the results database before plotting starts.
"""
from __future__ import annotations

from pathlib import Path
import textwrap

import duckdb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
DATABASE_PATH = BENCHMARK_DIR / "results" / "m5_results.duckdb"
FIGURE_DIR = BENCHMARK_DIR / "figures" / "publication"
SUPP_DIR = FIGURE_DIR / "supplementary"

METHOD_LABELS = {
    "bounded_ls": "Bounded LS",
    "lasso": "Lasso",
    "elastic_net": "Elastic net",
    "thermo_elastic_net": "Thermo elastic net",
    "hydrosheaf_guarded": "Hydrosheaf Guarded",
    "hydrosheaf_core": "Hydrosheaf-Core",
}
TIER_LABELS = {
    "core": "Core",
    "plus_lite": "Plus-lite",
    "enhanced": "Enhanced",
    "available_plus_lite": "Available Plus-lite",
}
PALETTE = {
    "hydrosheaf_guarded": "#0f766e",
    "hydrosheaf_core": "#2563eb",
    "thermo_elastic_net": "#7c3aed",
    "elastic_net": "#64748b",
    "lasso": "#94a3b8",
    "bounded_ls": "#475569",
    "core": "#64748b",
    "plus_lite": "#0f766e",
    "enhanced": "#b45309",
    "available_plus_lite": "#0f766e",
}


def setup_theme() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 450,
            "font.family": "DejaVu Sans",
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "grid.linewidth": 0.35,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def connect() -> duckdb.DuckDBPyConnection:
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            f"{DATABASE_PATH} not found. Run export_m5_results_database.py first."
        )
    return duckdb.connect(str(DATABASE_PATH), read_only=True)


def tables(con: duckdb.DuckDBPyConnection) -> set[str]:
    return {row[0] for row in con.execute("SHOW TABLES").fetchall()}


def q(con: duckdb.DuckDBPyConnection, sql: str) -> pd.DataFrame:
    return con.execute(sql).df()


def save(fig: plt.Figure, name: str, supplementary: bool = False) -> None:
    directory = SUPP_DIR if supplementary else FIGURE_DIR
    directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(directory / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(directory / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def label_panel(ax: plt.Axes, label: str) -> None:
    ax.text(
        0.015,
        0.985,
        label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        ha="left",
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5),
        zorder=10,
    )


def clean_label(value: object) -> str:
    text = str(value)
    return METHOD_LABELS.get(text, TIER_LABELS.get(text, text.replace("_", " ").title()))


def figure1_reproducible_workflow(con: duckdb.DuckDBPyConnection) -> None:
    catalog = q(con, "SELECT * FROM m5_table_catalog")
    summary = q(con, "SELECT * FROM analysis_summary").iloc[0].to_dict()
    fig = plt.figure(figsize=(7.2, 4.7), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1.0])
    ax_flow = fig.add_subplot(gs[0, :])
    ax_flow.axis("off")
    nodes = [
        ("Analysis\nPython + PHREEQC", 0.06, "#e0f2fe"),
        ("Locked results\nCSV + JSON", 0.30, "#dcfce7"),
        ("DuckDB evidence\nm5_results.duckdb", 0.54, "#fef3c7"),
        ("Publication figures\nR/Python plotting", 0.78, "#ede9fe"),
    ]
    for text, x, color in nodes:
        box = FancyBboxPatch(
            (x, 0.32),
            0.18,
            0.38,
            boxstyle="round,pad=0.03,rounding_size=0.025",
            facecolor=color,
            edgecolor="#334155",
            linewidth=0.8,
            transform=ax_flow.transAxes,
        )
        ax_flow.add_patch(box)
        ax_flow.text(
            x + 0.09,
            0.51,
            text,
            ha="center",
            va="center",
            fontsize=9,
            transform=ax_flow.transAxes,
        )
    for x in [0.245, 0.485, 0.725]:
        ax_flow.annotate(
            "",
            xy=(x + 0.045, 0.51),
            xytext=(x, 0.51),
            arrowprops=dict(arrowstyle="->", color="#334155", lw=1.0),
            xycoords=ax_flow.transAxes,
            textcoords=ax_flow.transAxes,
        )
    ax_flow.set_title("M5 separates analysis, evidence storage and plotting", loc="left")

    ax_counts = fig.add_subplot(gs[1, 0])
    counts = pd.DataFrame(
        {
            "item": ["PHREEQC\nscenarios", "Inverse\nfits", "Field\npairs", "DB\ntables"],
            "count": [
                summary.get("n_phreeqc_scenarios", 0),
                summary.get("n_benchmark_fits", 0),
                summary.get("n_field_pairs", 0),
                len(catalog),
            ],
        }
    )
    sns.barplot(data=counts, x="item", y="count", color="#0f766e", ax=ax_counts)
    ax_counts.set_xlabel("")
    ax_counts.set_ylabel("Count")
    ax_counts.set_title("Evidence volume", loc="left")
    label_panel(ax_counts, "a")

    ax_rows = fig.add_subplot(gs[1, 1])
    top = catalog.sort_values("row_count", ascending=False).head(8).copy()
    top["table_name"] = top["table_name"].str.replace("_", "\n")
    sns.barplot(data=top, y="table_name", x="row_count", color="#2563eb", ax=ax_rows)
    ax_rows.set_xlabel("Rows")
    ax_rows.set_ylabel("")
    ax_rows.set_title("Largest database tables", loc="left")
    label_panel(ax_rows, "b")

    ax_kind = fig.add_subplot(gs[1, 2])
    kind = catalog.groupby("source_kind", as_index=False)["row_count"].sum()
    sns.barplot(data=kind, x="source_kind", y="row_count", color="#b45309", ax=ax_kind)
    ax_kind.set_xlabel("")
    ax_kind.set_ylabel("Rows")
    ax_kind.set_title("Stored result sources", loc="left")
    label_panel(ax_kind, "c")
    save(fig, "figure1_database_workflow")


def figure2_model_performance(con: duckdb.DuckDBPyConnection) -> None:
    df = q(
        con,
        """
        SELECT method, phase_f1, class_f1, false_discovery_rate,
               extent_rmse_mmolL, reconstruction_rmse_mmolL
        FROM benchmark_fits
        WHERE panel = 'full_11' AND noise_level = 0.03
        """,
    )
    order = [
        "bounded_ls",
        "lasso",
        "elastic_net",
        "thermo_elastic_net",
        "hydrosheaf_guarded",
        "hydrosheaf_core",
    ]
    df["method_label"] = pd.Categorical(
        df["method"].map(clean_label), [clean_label(x) for x in order], ordered=True
    )
    summary = (
        df.groupby("method_label", observed=True)
        .agg(
            phase_f1=("phase_f1", "mean"),
            class_f1=("class_f1", "mean"),
            fdr=("false_discovery_rate", "mean"),
            extent_rmse=("extent_rmse_mmolL", "mean"),
            recon_rmse=("reconstruction_rmse_mmolL", "mean"),
        )
        .reset_index()
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.8), constrained_layout=True)
    metrics = [
        ("class_f1", "Equivalence-class F1", "#0f766e", "a"),
        ("fdr", "False-discovery rate", "#b45309", "b"),
        ("extent_rmse", "Extent RMSE (mmol/L)", "#2563eb", "c"),
    ]
    for ax, (metric, title, color, label) in zip(axes, metrics):
        sns.barplot(data=summary, y="method_label", x=metric, color=color, ax=ax)
        ax.set_xlabel(title)
        ax.set_ylabel("")
        ax.set_title(title, loc="left")
        label_panel(ax, label)
    save(fig, "figure2_database_model_performance")


def figure3_equifinality_elri(con: duckdb.DuckDBPyConnection) -> None:
    classes = q(con, "SELECT * FROM equivalence_classes WHERE ambiguous = TRUE")
    elri = q(con, "SELECT * FROM tableS16_evidence_lifted_resolution")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.7), constrained_layout=True)

    class_counts = classes.assign(members_wrapped=classes["members"].str.replace(";", "\n"))
    sns.barplot(
        data=class_counts,
        y="members_wrapped",
        x="n_members",
        color="#64748b",
        ax=axes[0],
    )
    axes[0].set_xlabel("Members")
    axes[0].set_ylabel("")
    axes[0].set_title("Structurally ambiguous classes", loc="left")
    label_panel(axes[0], "a")

    plot = elri.copy()
    plot["tier"] = pd.Categorical(
        plot["data_tier"].map(clean_label),
        ["Core", "Plus-lite", "Enhanced"],
        ordered=True,
    )
    plot["members"] = plot["members"].str.replace(";", " / ")
    sns.barplot(
        data=plot,
        y="members",
        x="mean_elri",
        hue="tier",
        palette=["#64748b", "#0f766e", "#b45309"],
        ax=axes[1],
    )
    axes[1].set_xlabel("Mean ELRI")
    axes[1].set_ylabel("")
    axes[1].set_title("Evidence-lifted resolution", loc="left")
    axes[1].legend(title="", frameon=False)
    label_panel(axes[1], "b")

    calcite = plot[plot["members"].eq("anorthite / calcite")]
    sns.pointplot(
        data=calcite,
        x="tier",
        y="top_member_true_active_fraction_when_class_active",
        color="#0f766e",
        ax=axes[2],
    )
    axes[2].set_ylim(0, 1)
    axes[2].set_xlabel("")
    axes[2].set_ylabel("Top member correct\nwhen class active")
    axes[2].set_title("Calcite/anorthite remains conditional", loc="left")
    axes[2].tick_params(axis="x", rotation=20)
    label_panel(axes[2], "c")
    save(fig, "figure3_database_equifinality_elri")


def figure4_data_tiers(con: duckdb.DuckDBPyConnection) -> None:
    tiers = q(con, "SELECT * FROM data_tier_experiment")
    elri = q(con, "SELECT * FROM data_tier_evidence_lifted_resolution")
    tiers["tier"] = pd.Categorical(
        tiers["data_tier"].map(clean_label),
        ["Core", "Plus-lite", "Enhanced"],
        ordered=True,
    )
    elri["tier"] = pd.Categorical(
        elri["data_tier"].map(clean_label),
        ["Core", "Plus-lite", "Enhanced"],
        ordered=True,
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.6), constrained_layout=True)
    sns.boxplot(data=tiers, x="tier", y="class_f1", color="#bae6fd", ax=axes[0])
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Class F1")
    axes[0].set_title("Recovery by evidence tier", loc="left")
    label_panel(axes[0], "a")
    sns.boxplot(
        data=tiers,
        x="tier",
        y="false_discovery_rate",
        color="#fed7aa",
        ax=axes[1],
    )
    axes[1].set_xlabel("")
    axes[1].set_ylabel("False-discovery rate")
    axes[1].set_title("Overinterpretation control", loc="left")
    label_panel(axes[1], "b")
    sns.boxplot(
        data=elri,
        x="tier",
        y="evidence_lifted_resolution_index",
        color="#bbf7d0",
        ax=axes[2],
    )
    axes[2].set_xlabel("")
    axes[2].set_ylabel("ELRI")
    axes[2].set_title("Ambiguity lifting", loc="left")
    label_panel(axes[2], "c")
    for ax in axes:
        ax.tick_params(axis="x", rotation=20)
    save(fig, "figure4_database_data_tiers")


def figure5_thermo_phreeqc(con: duckdb.DuckDBPyConnection) -> None:
    thermo = q(con, "SELECT * FROM thermodynamic_threshold_sensitivity")
    baseline = q(con, "SELECT * FROM phreeqc_inverse_baseline")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.6), constrained_layout=True)
    sns.lineplot(
        data=thermo,
        x="si_threshold",
        y="class_f1",
        hue="archetype",
        marker="o",
        ax=axes[0],
    )
    axes[0].set_xlabel("SI threshold")
    axes[0].set_ylabel("Class F1")
    axes[0].set_title("Thermodynamic gate sensitivity", loc="left")
    axes[0].legend(title="", frameon=False, fontsize=7)
    label_panel(axes[0], "a")

    sns.histplot(
        data=baseline,
        x="models_found",
        bins=20,
        color="#2563eb",
        ax=axes[1],
    )
    axes[1].set_xlabel("Feasible PHREEQC models")
    axes[1].set_ylabel("Scenarios")
    axes[1].set_title("Conventional inverse multiplicity", loc="left")
    label_panel(axes[1], "b")

    sns.scatterplot(
        data=baseline,
        x="minimal_models_found",
        y="first_minimal_class_f1",
        hue="archetype",
        s=24,
        alpha=0.75,
        ax=axes[2],
    )
    axes[2].set_xlabel("Minimal models found")
    axes[2].set_ylabel("First-minimal class F1")
    axes[2].set_title("Multiplicity lowers confidence", loc="left")
    axes[2].legend(title="", frameon=False, fontsize=7)
    label_panel(axes[2], "c")
    save(fig, "figure5_database_thermo_phreeqc")


def figure6_field_transfer(con: duckdb.DuckDBPyConnection) -> None:
    external = q(con, "SELECT * FROM external_field_evidence_lifted_resolution")
    ghana = q(con, "SELECT * FROM ghana_field_pairs")
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.8), constrained_layout=True)
    external["tier"] = external["data_tier"].map(clean_label)
    external["members"] = external["members"].str.replace(";", " / ")
    sns.boxplot(
        data=external,
        x="dataset",
        y="evidence_lifted_resolution_index",
        hue="tier",
        ax=axes[0],
        palette=["#64748b", "#0f766e"],
    )
    axes[0].set_xlabel("")
    axes[0].set_ylabel("ELRI")
    axes[0].set_title("External field transfer", loc="left")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].legend(title="", frameon=False, fontsize=7)
    label_panel(axes[0], "a")

    class_summary = (
        external.groupby(["dataset", "members"], as_index=False)[
            "evidence_lifted_resolution_index"
        ]
        .median()
        .pivot(index="members", columns="dataset", values="evidence_lifted_resolution_index")
    )
    sns.heatmap(
        class_summary,
        annot=True,
        fmt=".2f",
        cmap="YlGnBu",
        cbar_kws={"label": "Median ELRI"},
        ax=axes[1],
    )
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    axes[1].set_title("Class-specific field ambiguity", loc="left")
    label_panel(axes[1], "b")

    counts = (
        ghana["resolution_class"]
        .value_counts(normalize=True)
        .rename_axis("resolution_class")
        .reset_index(name="fraction")
    )
    counts["resolution_class"] = counts["resolution_class"].str.replace("_", "\n")
    sns.barplot(data=counts, x="resolution_class", y="fraction", color="#b45309", ax=axes[2])
    axes[2].set_ylim(0, 1)
    axes[2].set_xlabel("")
    axes[2].set_ylabel("Fraction")
    axes[2].set_title("2025 Ghana field resolution", loc="left")
    label_panel(axes[2], "c")
    save(fig, "figure6_database_field_transfer")


def supplementary_figures(con: duckdb.DuckDBPyConnection) -> None:
    # S1: reaction stoichiometry heatmap
    stoich = q(con, "SELECT * FROM tableS1_reaction_stoichiometry")
    ion_cols = [col for col in stoich.columns if col in {"Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"}]
    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    matrix = stoich.set_index("Reaction")[ion_cols]
    sns.heatmap(matrix, cmap="vlag", center=0, linewidths=0.3, ax=ax)
    ax.set_title("Full M5 reaction stoichiometric dictionary", loc="left")
    save(fig, "figureS1_database_reaction_dictionary", True)

    # S2: structural diagnostics
    diag = q(con, "SELECT * FROM identifiability_diagnostics")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), constrained_layout=True)
    sns.barplot(data=diag, x="panel", y="rank", color="#0f766e", ax=axes[0])
    axes[0].tick_params(axis="x", rotation=35)
    axes[0].set_title("Panel rank", loc="left")
    sns.barplot(data=diag, x="panel", y="nullity", color="#b45309", ax=axes[1])
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].set_title("Reaction-space nullity", loc="left")
    save(fig, "figureS2_database_identifiability_diagnostics", True)

    # S3: method/panel/noise class F1
    fits = q(con, "SELECT * FROM benchmark_fits")
    pivot = fits.pivot_table(
        index="method",
        columns="panel",
        values="class_f1",
        aggfunc="mean",
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    sns.heatmap(pivot, cmap="YlGnBu", annot=True, fmt=".2f", ax=ax)
    ax.set_title("Mean equivalence-class F1 across panels", loc="left")
    save(fig, "figureS3_database_panel_method_heatmap", True)

    # S4: reaction confusion
    confusion = q(con, "SELECT * FROM tableS6_reaction_confusion_matrices")
    hyd = confusion[confusion["Method"].eq("hydrosheaf_core")].copy()
    if hyd.empty:
        hyd = confusion.copy()
    hyd["F1_proxy"] = 2 * hyd["Precision"] * hyd["Recall"] / (
        hyd["Precision"] + hyd["Recall"]
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    sns.barplot(
        data=hyd.sort_values("F1_proxy", ascending=False),
        x="F1_proxy",
        y="Reaction",
        color="#2563eb",
        ax=ax,
    )
    ax.set_xlabel("Reaction F1 proxy")
    ax.set_title("Hydrosheaf-Core reaction recovery profile", loc="left")
    save(fig, "figureS4_database_reaction_recovery", True)

    # S5: regularization paths
    paths = q(con, "SELECT * FROM regularization_paths")
    numeric_reactions = [
        col
        for col in paths.columns
        if col not in {"scenario_id", "archetype", "lambda_l1", "lambda_l2"}
        and pd.api.types.is_numeric_dtype(paths[col])
    ][:10]
    long = paths.melt(
        id_vars=["lambda_l1"],
        value_vars=numeric_reactions,
        var_name="reaction",
        value_name="extent",
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    sns.lineplot(data=long, x="lambda_l1", y="extent", hue="reaction", ax=ax, lw=0.9)
    ax.set_xscale("log")
    ax.set_title("Regularization path extents", loc="left")
    ax.legend(ncol=2, frameon=False, fontsize=6)
    save(fig, "figureS5_database_regularization_paths", True)

    # S6: measurement value
    nb = q(con, "SELECT * FROM next_best_measurement")
    ranking = (
        nb.groupby("candidate_measurement", as_index=False)["measurement_value_score"]
        .mean()
        .sort_values("measurement_value_score", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(6.6, 3.6), constrained_layout=True)
    sns.barplot(data=ranking, x="measurement_value_score", y="candidate_measurement", color="#0f766e", ax=ax)
    ax.set_xlabel("Mean measurement value score")
    ax.set_ylabel("")
    ax.set_title("Next-best measurement ranking", loc="left")
    save(fig, "figureS6_database_measurement_value", True)

    # S7: bootstrap support
    bootstrap = q(con, "SELECT * FROM bootstrap_support_stability")
    boot = (
        bootstrap.groupby("reaction", as_index=False)["support_frequency"]
        .mean()
        .sort_values("support_frequency", ascending=False)
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    sns.barplot(data=boot, x="support_frequency", y="reaction", color="#7c3aed", ax=ax)
    ax.set_xlabel("Mean bootstrap support frequency")
    ax.set_ylabel("")
    ax.set_title("Support stability across bootstrap perturbations", loc="left")
    save(fig, "figureS7_database_bootstrap_support", True)

    # S8: evidence gates
    evidence = q(con, "SELECT * FROM tableS12_hydrosheaf_core_evidence_gates")
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    sns.scatterplot(
        data=evidence,
        x="mean_evidence_score",
        y="median_penalty_scale",
        hue="family",
        size="true_active_fraction",
        ax=ax,
    )
    ax.set_title("Hydrosheaf-Core evidence gates", loc="left")
    ax.legend(frameon=False, fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left")
    save(fig, "figureS8_database_core_evidence", True)

    # S9: data-tier reaction evidence
    tier_evidence = q(con, "SELECT * FROM tableS15_data_tier_reaction_evidence")
    pivot = tier_evidence.pivot_table(
        index="reaction",
        columns="data_tier",
        values="mean_combined_evidence_score",
    )
    fig, ax = plt.subplots(figsize=(6.8, 5.0), constrained_layout=True)
    sns.heatmap(pivot, cmap="viridis", annot=True, fmt=".2f", ax=ax)
    ax.set_title("Combined reaction evidence by data tier", loc="left")
    save(fig, "figureS9_database_data_tier_evidence", True)

    # S10: external field ELRI
    external = q(con, "SELECT * FROM tableS17_external_field_evidence_lifted_resolution")
    external["members"] = external["members"].str.replace(";", " / ")
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    sns.barplot(
        data=external,
        y="members",
        x="median_elri",
        hue="dataset",
        ax=ax,
    )
    ax.set_xlabel("Median ELRI")
    ax.set_ylabel("")
    ax.set_title("External field-transfer ELRI by ambiguous class", loc="left")
    ax.legend(frameon=False, fontsize=7)
    save(fig, "figureS10_database_external_field_elri", True)


def main() -> None:
    setup_theme()
    with connect() as con:
        available = tables(con)
        required = {
            "benchmark_fits",
            "data_tier_experiment",
            "data_tier_evidence_lifted_resolution",
            "external_field_evidence_lifted_resolution",
            "m5_table_catalog",
        }
        missing = required - available
        if missing:
            raise RuntimeError(f"Database is missing required tables: {sorted(missing)}")
        figure1_reproducible_workflow(con)
        figure2_model_performance(con)
        figure3_equifinality_elri(con)
        figure4_data_tiers(con)
        figure5_thermo_phreeqc(con)
        figure6_field_transfer(con)
        supplementary_figures(con)
    print("Generated database-backed M5 publication figures.")


if __name__ == "__main__":
    main()
