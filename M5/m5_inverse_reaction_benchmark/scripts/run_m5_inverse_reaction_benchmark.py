"""Run M5 sparse inverse-reaction validation benchmarks."""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf.config import Config
from hydrosheaf.validation import validate_sparse_inverse_reaction_model


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
TABLES_DIR = BENCHMARK_DIR / "tables"
DOCS_DIR = BENCHMARK_DIR / "docs"


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_violations_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["reaction", "extent", "lower_bound", "upper_bound", "constraint"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def flatten_l1_rows(rows):
    out = []
    for row in rows:
        out.append(
            {
                "lambda_l1": row["lambda_l1"],
                "transport_residual_norm": row["transport_residual_norm"],
                "reaction_residual_norm": row["reaction_residual_norm"],
                "residual_reduction_fraction": row["residual_reduction_fraction"],
                "l1_norm": row["l1_norm"],
                "n_selected_reactions": row["n_selected_reactions"],
                "selected_reactions": ";".join(item["reaction"] for item in row["selected_reactions"]),
                "converged": row["converged"],
            }
        )
    return out


def flatten_missing_rows(rows):
    out = []
    for row in rows:
        out.append(
            {
                "missing_ions": ";".join(row["missing_ions"]),
                "lambda_l1": row["lambda_l1"],
                "reaction_residual_norm": row["reaction_residual_norm"],
                "l1_norm": row["l1_norm"],
                "n_selected_reactions": row["n_selected_reactions"],
                "selected_reactions": ";".join(item["reaction"] for item in row["selected_reactions"]),
                "converged": row["converged"],
            }
        )
    return out


def main() -> None:
    for directory in (RESULTS_DIR, TABLES_DIR, DOCS_DIR):
        directory.mkdir(parents=True, exist_ok=True)

    config = Config(
        ion_order=["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3"],
        weights=[1.0, 1.0, 0.4, 1.0, 0.4, 1.0, 1.0],
        conservative_weights=[1.0, 1.0, 0.4, 1.0, 0.4, 1.0, 1.0],
    )
    reaction_labels = ["calcite", "dolomite", "gypsum", "halite", "denit"]
    reaction_matrix = [
        [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.8, 0.0, 0.0, -1.0],
    ]
    upstream = [1.0, 0.5, 2.0, 3.0, 2.0, 1.0, 0.6]
    post_transport = [1.0, 0.5, 2.0, 3.0, 2.0, 1.0, 0.6]
    downstream = [2.2, 0.5, 2.15, 4.15, 2.15, 2.05, 0.15]
    phreeqc_bounds = {
        "lb": [0.0, 0.0, 0.0, 0.0, 0.0],
        "ub": [float("inf"), 0.0, float("inf"), float("inf"), float("inf")],
        "constraints_active": {
            "calcite": "dissolution_only",
            "dolomite": "precipitation_only",
            "gypsum": "dissolution_only",
            "halite": "dissolution_only",
            "denit": "dissolution_only",
        },
    }

    report = validate_sparse_inverse_reaction_model(
        upstream=upstream,
        downstream=downstream,
        post_transport=post_transport,
        reaction_matrix=reaction_matrix,
        reaction_labels=reaction_labels,
        config=config,
        lambda_grid=[0.0, 0.02, 0.1, 0.5, 1.0],
        missing_ion_sets=[["SO4"], ["NO3"], ["HCO3"], ["Ca", "SO4"]],
        phreeqc_bounds=phreeqc_bounds,
    )

    write_csv(RESULTS_DIR / "l1_penalty_sensitivity.csv", flatten_l1_rows(report["l1_penalty_sensitivity"]))
    write_csv(RESULTS_DIR / "missing_ion_sensitivity.csv", flatten_missing_rows(report["missing_ion_sensitivity"]))
    write_violations_csv(
        RESULTS_DIR / "thermodynamic_bound_violations.csv",
        report["thermodynamic_bound_violations"],
    )
    summary = [
        {
            "model_framing": report["model_framing"],
            "best_lambda_l1": report["best_fit"]["lambda_l1"],
            "best_reaction_residual_norm": report["best_fit"]["reaction_residual_norm"],
            "best_l1_norm": report["best_fit"]["l1_norm"],
            "best_n_selected_reactions": report["best_fit"]["n_selected_reactions"],
            "selected_reactions": ";".join(item["reaction"] for item in report["best_fit"]["selected_reactions"]),
            "thermodynamic_bound_violation": report["flags"]["thermodynamic_bound_violation"],
            "missing_ion_sensitive": report["flags"]["missing_ion_sensitive"],
            "l1_penalty_sensitive": report["flags"]["l1_penalty_sensitive"],
            "linearized_dictionary_limit": report["flags"]["linearized_dictionary_limit"],
            "not_fully_coupled_phreeqc_inverse_solver": report[
                "not_a_fully_coupled_nonlinear_phreeqc_inverse_solver"
            ],
        }
    ]
    write_csv(TABLES_DIR / "table1_inverse_reaction_validation_summary.csv", summary)
    lines = [
        "# M5 Inverse Reaction Benchmark Results",
        "",
        f"- L1 penalty rows tested: {len(report['l1_penalty_sensitivity'])}.",
        f"- Missing-ion stress tests: {len(report['missing_ion_sensitivity'])}.",
        f"- Thermodynamic bound violations in best fit: {len(report['thermodynamic_bound_violations'])}.",
        f"- Best selected reactions: {summary[0]['selected_reactions']}.",
        "",
        "Interpretation guardrail:",
        "",
        report["claim_guardrail"],
        "",
        "This benchmark supports sparse linear inverse reaction fitting with PHREEQC",
        "screening and forward-validation diagnostics. It does not validate Hydrosheaf",
        "as a fully coupled nonlinear PHREEQC inverse solver.",
        "",
    ]
    (DOCS_DIR / "m5_results_summary.md").write_text("\n".join(lines), encoding="utf-8")
    print("Wrote M5 sparse inverse-reaction benchmark outputs.")


if __name__ == "__main__":
    main()
