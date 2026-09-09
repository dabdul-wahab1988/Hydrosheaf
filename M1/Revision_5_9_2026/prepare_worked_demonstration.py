from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[2]
ROOT = Path(__file__).resolve().parent
SOURCE = REPO / "M7" / "m7_nonuniqueness_benchmark"
RESULTS = SOURCE / "results" / "supporting_validation"
OUT = ROOT / "Worked_Demonstration"
FIGURES = ROOT / "Figures"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_value(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO, text=True, encoding="utf-8"
    ).strip()


def load_inputs() -> dict[str, object]:
    return {
        "manifest": json.loads((RESULTS / "manifest.json").read_text(encoding="utf-8")),
        "methods": pd.read_csv(RESULTS / "method_summary.csv"),
        "ages": pd.read_csv(RESULTS / "bayesian_age_diagnostics.csv"),
        "topology": pd.read_csv(RESULTS / "topology_posterior_diagnostics.csv"),
        "phreeqc": pd.read_csv(RESULTS / "phreeqc_constraint_audit.csv"),
        "reactions": pd.read_csv(RESULTS / "reaction_summary.csv"),
        "bootstrap": pd.read_csv(RESULTS / "age_incremental_bootstrap.csv"),
        "edges": pd.read_csv(RESULTS / "locked_test_edge_results.csv"),
    }


def validate(data: dict[str, object]) -> None:
    manifest = data["manifest"]
    methods = data["methods"]
    ages = data["ages"]
    topology = data["topology"]
    phreeqc = data["phreeqc"]
    edges = data["edges"]
    assert manifest["protocol_stage"] == "fresh_seed_confirmatory_after_initial_age_contrast_failure"
    assert manifest["external_generator"] == "official MODFLOW 6 + MODPATH 7 + independent nonlinear chemistry"
    assert manifest["locked_test_seeds"] == list(range(4101, 4113))
    assert len(ages) == 12 and bool(ages["converged"].all())
    assert len(topology) == 12 and bool(topology["converged"].all())
    assert int(edges.shape[0]) == 825 and int(edges["is_true_edge"].sum()) == 103
    assert int(phreeqc["n_samples"].sum()) == 144
    assert int(phreeqc["n_successful_samples"].sum()) == 144
    assert int(phreeqc["n_candidate_edge_fits"].sum()) == 825
    assert set(methods["method"]) == {
        "hydraulic_chemistry",
        "hydraulic_chemistry_age",
        "age_permuted_control",
    }


def metric_summary(data: dict[str, object]) -> dict[str, object]:
    manifest = data["manifest"]
    methods = data["methods"].set_index("method")
    ages = data["ages"]
    topology = data["topology"]
    phreeqc = data["phreeqc"]
    reactions = data["reactions"].set_index("true_process")
    bootstrap = data["bootstrap"]
    edges = data["edges"]

    def boot(contrast: str, metric: str) -> dict[str, float]:
        row = bootstrap[
            (bootstrap["contrast"] == contrast) & (bootstrap["metric"] == metric)
        ].iloc[0]
        return {
            "mean_difference": float(row["mean_difference"]),
            "ci95_low": float(row["ci95_low"]),
            "ci95_high": float(row["ci95_high"]),
        }

    all_reactions = reactions.loc["ALL"]
    summary = {
        "evidence_class": "controlled-synthetic integration benchmark",
        "claim_boundary": (
            "Model-conditioned synthetic truth tests execution and internal integration; "
            "it is not independent field validation or evidence of management-ready performance."
        ),
        "design": {
            "development_cases": int(manifest["n_development_cases"]),
            "locked_test_cases": int(manifest["n_locked_test_cases"]),
            "nodes_per_case": 12,
            "locked_test_seeds": manifest["locked_test_seeds"],
            "candidate_edges": int(len(edges)),
            "candidate_contained_true_edges": int(edges["is_true_edge"].sum()),
            "mean_candidate_recall": float(manifest["candidate_recall_test"]),
            "external_generator": manifest["external_generator"],
            "truth_blind_inference": True,
        },
        "age_inference": {
            "tracers": ["3H", "39Ar"],
            "model": "piston-flow model with exact-grid Bayesian quadrature",
            "draws_per_chain": 500,
            "chains": 4,
            "mean_mae_years": float(ages["mae_years"].mean()),
            "mean_bias_years": float(ages["bias_years"].mean()),
            "mean_interval95_coverage": float(ages["interval95_coverage"].mean()),
            "max_r_hat": float(ages["r_hat_max"].max()),
            "minimum_bulk_ess": float(ages["ess_bulk_min"].min()),
            "minimum_tail_ess": float(ages["ess_tail_min"].min()),
            "divergences": int(ages["divergences"].sum()),
            "converged_cases": int(ages["converged"].sum()),
        },
        "edge_scoring": {
            method: {
                key: float(methods.loc[method, key])
                for key in ("pr_auc", "roc_auc", "brier", "precision", "recall", "f1", "mcc")
            }
            for method in methods.index
        },
        "topology_posterior": {
            "retained_samples_per_chain": 2500,
            "burn_in_per_chain": 750,
            "chains": 4,
            "updates_per_retained_sample": int(manifest["topology_updates_per_sample"]),
            "mean_map_precision": float(topology["map_precision"].mean()),
            "mean_map_recall": float(topology["map_recall"].mean()),
            "mean_map_f1": float(topology["map_f1"].mean()),
            "max_edge_r_hat": float(topology["edge_r_hat_max"].max()),
            "minimum_edge_ess": float(topology["edge_ess_min"].min()),
            "converged_cases": int(topology["converged"].sum()),
        },
        "inverse_reaction_fitting": {
            "phreeqc_samples": int(phreeqc["n_samples"].sum()),
            "successful_phreeqc_samples": int(phreeqc["n_successful_samples"].sum()),
            "candidate_edge_fits": int(phreeqc["n_candidate_edge_fits"].sum()),
            "fits_with_active_direction_constraints": int(
                phreeqc["n_edges_with_active_direction_constraints"].sum()
            ),
            "fits_with_material_objective_change": int(
                phreeqc["n_edges_with_material_objective_change"].sum()
            ),
            "candidate_contained_truth_edges": int(all_reactions["n"]),
            "constrained_family_accuracy": float(all_reactions["constrained_family_accuracy"]),
            "unconstrained_family_accuracy": float(all_reactions["unconstrained_family_accuracy"]),
            "material_change_fraction_on_truth_edges": float(
                all_reactions["phreeqc_material_change_fraction"]
            ),
            "known_limitation": (
                "The constrained classifier recovered denitrification, sulfate reduction, "
                "silicate weathering and most iron reduction cases, but not either carbonate family."
            ),
        },
        "uncertainty": {
            "case_block_bootstrap_resamples": 10000,
            "resampling_unit": "independent MODFLOW case",
            "bootstrap_seed": 7781,
            "age_minus_hydraulic_chemistry_f1": boot(
                "hydraulic_chemistry_age_minus_hydraulic_chemistry", "f1"
            ),
            "age_minus_hydraulic_chemistry_pr_auc": boot(
                "hydraulic_chemistry_age_minus_hydraulic_chemistry", "pr_auc"
            ),
            "age_minus_permuted_f1": boot(
                "hydraulic_chemistry_age_minus_age_permuted_control", "f1"
            ),
            "age_minus_permuted_pr_auc": boot(
                "hydraulic_chemistry_age_minus_age_permuted_control", "pr_auc"
            ),
            "failed_runs": 0,
            "retention_rule": (
                "All pre-specified locked cases were retained; convergence and PHREEQC success "
                "were reported rather than used to remove cases."
            ),
        },
    }
    return summary


def write_tables(data: dict[str, object], summary: dict[str, object]) -> None:
    methods = data["methods"].copy()
    methods.insert(0, "workflow_component", "locked edge scoring")
    reactions = data["reactions"].copy()
    reactions.insert(0, "workflow_component", "inverse reaction family")
    pd.concat([methods, reactions], ignore_index=True, sort=False).to_csv(
        OUT / "worked_demonstration_metrics.csv", index=False
    )

    case_dir = RESULTS / "cases" / "locked_test_4101"
    diag = json.loads((case_dir / "diagnostics.json").read_text(encoding="utf-8"))
    case_edges = data["edges"][data["edges"]["seed"] == 4101].copy()
    case_edges["map_edge"] = case_edges["edge_id"].isin(diag["topology"]["map_edges"])
    case_edges[
        [
            "edge_id", "u", "v", "is_true_edge", "hydraulic_probability",
            "age_cost", "chemistry_objective_constrained", "fusion_probability",
            "dominant_reaction_family", "true_process", "map_edge",
        ]
    ].to_csv(OUT / "case_4101_edge_audit.csv", index=False)

    (OUT / "worked_demonstration_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def draw_uncertainty_figure() -> None:
    fig, ax = plt.subplots(figsize=(13.2, 7.4))
    ax.set_xlim(0, 13.2)
    ax.set_ylim(0, 7.4)
    ax.axis("off")
    ax.text(6.6, 7.08, "Uncertainty is sampled through the complete inference chain",
            ha="center", va="center", fontsize=17, weight="bold", color="#17365D")
    columns = [
        (0.25, 2.25, "Evidence source", "#D9EAF7"),
        (2.65, 5.85, "Quantification in each draw", "#E8F2E3"),
        (6.10, 8.65, "Propagated state", "#FBE5D6"),
        (8.90, 12.95, "Reported diagnostic", "#E4DFEC"),
    ]
    for x0, x1, label, colour in columns:
        ax.add_patch(FancyBboxPatch((x0, 6.35), x1-x0, 0.48,
                                   boxstyle="round,pad=0.03", fc=colour, ec="#666666"))
        ax.text((x0+x1)/2, 6.59, label, ha="center", va="center", fontsize=11, weight="bold")
    rows = [
        ("Tracer measurements", "Analytical-error and\ndetection-limit model", "Tracer realisation", "Age/RTD distribution; fit residuals"),
        ("Recharge forcing", "Alternative input histories\nand age-model classes", "Age-model realisation", "Age-class probability; sensitivity"),
        ("Hydraulic connectivity", "Candidate-edge ledger; stochastic\ngraph or particle-track ensemble", "Directed graph realisation", "Edge-inclusion probability; topology CI"),
        ("Reaction system", "Phase-list alternatives; analytical\nerror; saturation/redox filters", "Constrained reaction fit", "Viable solutions; objective change"),
        ("Model structure", "Pre-specified model alternatives;\ncase-block bootstrap", "Integrated diagnosis", "Cross-method consistency; 95% CI"),
    ]
    ys = [5.55, 4.48, 3.41, 2.34, 1.27]
    colors = ["#D9EAF7", "#E8F2E3", "#FBE5D6", "#FFF2CC"]
    for y, row in zip(ys, rows):
        spans = [(0.25, 2.25), (2.65, 5.85), (6.10, 8.65), (8.90, 12.95)]
        for idx, ((x0, x1), text) in enumerate(zip(spans, row)):
            ax.add_patch(FancyBboxPatch((x0, y-0.36), x1-x0, 0.72,
                                       boxstyle="round,pad=0.04", fc=colors[idx], ec="#7F7F7F"))
            ax.text((x0+x1)/2, y, text, ha="center", va="center", fontsize=9.2)
        for x0, x1 in [(2.25, 2.65), (5.85, 6.10), (8.65, 8.90)]:
            ax.annotate("", xy=(x1-0.03, y), xytext=(x0+0.03, y),
                        arrowprops=dict(arrowstyle="-|>", color="#4472C4", lw=1.5))
    ax.text(6.6, 0.48,
            "Report draws, seeds, convergence, failures and retention rules. Independent field evidence remains a separate validation gate.",
            ha="center", va="center", fontsize=10.2, color="#7F6000",
            bbox=dict(boxstyle="round,pad=0.35", fc="#FFF2CC", ec="#BF9000"))
    fig.tight_layout()
    fig.savefig(FIGURES / "Figure_4_Fully_Revised.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "Figure_4_Fully_Revised.svg", bbox_inches="tight")
    plt.close(fig)


def draw_worked_demo_figure(data: dict[str, object], summary: dict[str, object]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.0))
    fig.suptitle("Controlled-synthetic worked demonstration on 12 locked aquifer twins",
                 fontsize=16, weight="bold", color="#17365D", y=0.985)

    # A: representative graph from locked seed 4101.
    ax = axes[0, 0]
    case_dir = RESULTS / "cases" / "locked_test_4101"
    obs = pd.read_csv(case_dir / "blind_observations.csv")
    truth = pd.read_csv(case_dir / "heldout_truth.csv")
    diag = json.loads((case_dir / "diagnostics.json").read_text(encoding="utf-8"))
    xy = {r.site_id: (r.x_m, r.y_m) for r in obs.itertuples()}
    truth_ids = set(truth["edge_id"])
    for edge_id in diag["topology"]["map_edges"]:
        u, v = edge_id.split("->")
        color = "#2E8B57" if edge_id in truth_ids else "#C44E52"
        ax.annotate("", xy=xy[v], xytext=xy[u],
                    arrowprops=dict(arrowstyle="->", color=color, alpha=0.75, lw=1.4))
    ax.scatter(obs["x_m"], obs["y_m"], s=48, c="#4472C4", ec="white", zorder=3)
    for r in obs.itertuples():
        ax.text(r.x_m, r.y_m+55, r.site_id.split("_")[-1], fontsize=7, ha="center")
    ax.set_title("A  Posterior MAP graph, seed 4101\n(green=true positive; red=false positive)", loc="left", fontsize=11, weight="bold")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.grid(alpha=0.2)

    # B: locked edge discrimination.
    ax = axes[0, 1]
    methods = data["methods"].set_index("method")
    labels = ["Hydraulic +\nchemistry", "+ age gate", "Age-permuted\ncontrol"]
    keys = ["hydraulic_chemistry", "hydraulic_chemistry_age", "age_permuted_control"]
    x = np.arange(3)
    width = 0.34
    ax.bar(x-width/2, [methods.loc[k, "pr_auc"] for k in keys], width, label="PR-AUC", color="#4472C4")
    ax.bar(x+width/2, [methods.loc[k, "f1"] for k in keys], width, label="F1", color="#ED7D31")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 0.62)
    ax.set_ylabel("Locked-test score")
    ax.set_title("B  Edge scoring", loc="left", fontsize=11, weight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)

    # C: age inference by independent test aquifer.
    ax = axes[1, 0]
    ages = data["ages"]
    ax.plot(ages["seed"], ages["mae_years"], marker="o", color="#4472C4", label="MAE (years)")
    ax2 = ax.twinx()
    ax2.plot(ages["seed"], ages["interval95_coverage"], marker="s", color="#70AD47", label="95% coverage")
    ax.axhline(ages["mae_years"].mean(), color="#4472C4", ls="--", lw=1)
    ax2.axhline(ages["interval95_coverage"].mean(), color="#70AD47", ls="--", lw=1)
    ax.set_xlabel("Locked test seed")
    ax.set_ylabel("Age MAE (years)", color="#4472C4")
    ax2.set_ylabel("Interval coverage", color="#548235")
    ax2.set_ylim(0.68, 1.04)
    ax.set_title("C  Bayesian 3H + 39Ar age inference", loc="left", fontsize=11, weight="bold")
    ax.grid(alpha=0.2)

    # D: reaction performance and uncertainty guardrail.
    ax = axes[1, 1]
    all_rxn = data["reactions"].set_index("true_process").loc["ALL"]
    ax.bar([0, 1], [all_rxn["unconstrained_family_accuracy"], all_rxn["constrained_family_accuracy"]],
           color=["#A5A5A5", "#70AD47"], width=0.62)
    ax.set_xticks([0, 1], ["Unconstrained", "PHREEQC-\nconstrained"])
    ax.set_ylim(0, 0.75)
    ax.set_ylabel("Reaction-family accuracy")
    ax.set_title("D  Inverse reaction fits and propagated uncertainty", loc="left", fontsize=11, weight="bold")
    ax.grid(axis="y", alpha=0.2)
    ax.text(0.5, 0.70,
            "103 candidate-contained truth edges\n825/825 fits constrained; 681 materially changed\n10,000 case-block bootstraps; zero failed runs\nAge gate: ΔF1 = 0.000; ΔPR-AUC = −0.0017",
            transform=ax.transAxes, ha="center", va="top", fontsize=9.2,
            bbox=dict(boxstyle="round,pad=0.4", fc="#FFF2CC", ec="#BF9000"))
    ax.text(0.5, 0.07, "Controlled synthetic evidence — not field validation",
            transform=ax.transAxes, ha="center", fontsize=9.5, color="#C00000", weight="bold")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIGURES / "Figure_7_Worked_Demonstration.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIGURES / "Figure_7_Worked_Demonstration.svg", bbox_inches="tight")
    plt.close(fig)


def write_manifest(summary: dict[str, object]) -> None:
    source_files = [
        RESULTS / "manifest.json",
        RESULTS / "method_summary.csv",
        RESULTS / "bayesian_age_diagnostics.csv",
        RESULTS / "topology_posterior_diagnostics.csv",
        RESULTS / "phreeqc_constraint_audit.csv",
        RESULTS / "reaction_summary.csv",
        RESULTS / "age_incremental_bootstrap.csv",
        RESULTS / "locked_test_edge_results.csv",
        SOURCE / "scripts" / "run_supporting_validation.py",
        SOURCE / "scripts" / "strong_inference.py",
    ]
    manifest = {
        "package_id": "M1-WORKED-DEMONSTRATION-20260906",
        "created_utc_date": date.today().isoformat(),
        "purpose": "Publication extract from the locked M7.2 controlled-synthetic integration benchmark.",
        "repository_commit": git_value("rev-parse", "HEAD"),
        "repository_status": git_value("status", "--short", "--branch"),
        "source_root": str(SOURCE.relative_to(REPO)).replace("\\", "/"),
        "source_files": [
            {
                "path": str(path.relative_to(REPO)).replace("\\", "/"),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in source_files
        ],
        "validation": {
            "input_assertions": "PASS",
            "n_locked_test_cases": summary["design"]["locked_test_cases"],
            "failed_runs": summary["uncertainty"]["failed_runs"],
            "claim_guardrail": summary["claim_boundary"],
        },
        "generated_outputs": [
            "worked_demonstration_summary.json",
            "worked_demonstration_metrics.csv",
            "case_4101_edge_audit.csv",
            "../Figures/Figure_4_Fully_Revised.png",
            "../Figures/Figure_4_Fully_Revised.svg",
            "../Figures/Figure_7_Worked_Demonstration.png",
            "../Figures/Figure_7_Worked_Demonstration.svg",
        ],
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    readme = """# M1 controlled-synthetic worked demonstration

This directory is a publication extract from the locked M7.2 supporting-validation run. The six development aquifer twins (seeds 2101–2106) were used to freeze the edge-fusion models; twelve independent test aquifer twins (seeds 4101–4112) were then evaluated without exposing held-out truth to inference. The external generator uses official MODFLOW 6 and MODPATH 7 plus independent nonlinear chemistry.

The package demonstrates executable integration of 3H/39Ar age inference, candidate-graph construction, edge scoring, PHREEQC-constrained reaction fitting, topology-posterior sampling and case-block bootstrap uncertainty. It does not constitute field validation. In particular, the locked aggregate comparison found no improvement from the age gate over hydraulic-plus-chemistry scoring, and carbonate-family reaction attribution remained unresolved.

Reproduce the extract from the repository root with:

    .\\.venv\\Scripts\\python.exe M1\\Revision_5_9_2026\\prepare_worked_demonstration.py

The script validates the locked file dimensions and convergence flags before rewriting any derived output. `analysis_manifest.json` records source paths, SHA-256 hashes, the repository commit and the working-tree status.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    validate(data)
    summary = metric_summary(data)
    write_tables(data, summary)
    draw_uncertainty_figure()
    draw_worked_demo_figure(data, summary)
    write_manifest(summary)
    print(json.dumps({
        "status": "PASS",
        "locked_cases": summary["design"]["locked_test_cases"],
        "candidate_edges": summary["design"]["candidate_edges"],
        "outputs": 7,
    }))


if __name__ == "__main__":
    main()
