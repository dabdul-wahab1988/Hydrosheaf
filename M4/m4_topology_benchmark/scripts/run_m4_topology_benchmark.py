"""Run controlled M4 topology-validation benchmarks."""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf.validation import (
    apply_modpath_informed_graph_priors,
    validate_independent_graph_against_modpath,
)


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
TABLES_DIR = BENCHMARK_DIR / "tables"
DOCS_DIR = BENCHMARK_DIR / "docs"

Edge = Tuple[str, str]


def scenarios() -> List[Dict[str, object]]:
    reference = [("R", "A"), ("A", "B"), ("B", "C"), ("C", "D")]
    candidates = [
        ("R", "A"),
        ("R", "B"),
        ("A", "B"),
        ("A", "C"),
        ("B", "C"),
        ("B", "D"),
        ("C", "D"),
        ("A", "D"),
    ]
    lengths = {
        "R->A": 100.0,
        "R->B": 340.0,
        "A->B": 120.0,
        "A->C": 360.0,
        "B->C": 130.0,
        "B->D": 400.0,
        "C->D": 140.0,
        "A->D": 600.0,
    }
    return [
        {
            "scenario": "well_resolved_chain",
            "purpose": "Independent graph recovers MODPATH chain with no errors.",
            "reference": reference,
            "inferred": reference,
            "candidates": candidates,
            "lengths": lengths,
        },
        {
            "scenario": "false_positive_shortcut",
            "purpose": "Hydrosheaf adds a plausible but non-MODPATH shortcut edge.",
            "reference": reference,
            "inferred": [("R", "A"), ("A", "B"), ("A", "C"), ("B", "C"), ("C", "D")],
            "candidates": candidates,
            "lengths": lengths,
        },
        {
            "scenario": "false_negative_missing_link",
            "purpose": "Hydrosheaf misses one MODPATH edge and preserves a partial chain.",
            "reference": reference,
            "inferred": [("R", "A"), ("A", "B"), ("C", "D")],
            "candidates": candidates,
            "lengths": lengths,
        },
        {
            "scenario": "scale_mismatch_shortcuts",
            "purpose": "Hydrosheaf infers coarse shortcut edges instead of local MODPATH links.",
            "reference": reference,
            "inferred": [("R", "B"), ("A", "C"), ("B", "D")],
            "candidates": candidates,
            "lengths": lengths,
        },
    ]


def edge_key(edge: Edge) -> str:
    return f"{edge[0]}->{edge[1]}"


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def independent_rows() -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    metric_rows: List[Dict[str, object]] = []
    edge_rows: List[Dict[str, object]] = []
    for scenario in scenarios():
        report = validate_independent_graph_against_modpath(
            scenario["inferred"],  # type: ignore[arg-type]
            scenario["reference"],  # type: ignore[arg-type]
            candidate_edges=scenario["candidates"],  # type: ignore[arg-type]
            edge_lengths=scenario["lengths"],  # type: ignore[arg-type]
        )
        metrics = report["metrics"]
        scale = report["scale_mismatch"]
        metric_rows.append(
            {
                "scenario": scenario["scenario"],
                "validation_mode": report["validation_mode"],
                "purpose": scenario["purpose"],
                "precision": metrics["precision"],
                "recall": metrics["recall"],
                "f1": metrics["f1"],
                "false_positive_rate": metrics["false_positive_rate"],
                "false_negative_rate": metrics["false_negative_rate"],
                "tp": metrics["tp"],
                "fp": metrics["fp"],
                "fn": metrics["fn"],
                "tn": metrics["tn"],
                "scale_mismatch": scale["scale_mismatch"],
                "median_reference_length": scale["median_reference_length"],
                "median_inferred_length": scale["median_inferred_length"],
            }
        )

        classes = [
            ("TP", metrics["true_positives"]),
            ("FP", metrics["false_positives"]),
            ("FN", metrics["false_negatives"]),
            ("TN", metrics["true_negatives"]),
        ]
        for label, edges in classes:
            for edge in edges:
                edge_rows.append(
                    {
                        "scenario": scenario["scenario"],
                        "edge": edge_key(edge),
                        "classification": label,
                    }
                )
    return metric_rows, edge_rows


def prior_rows() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    reference = scenarios()[0]["reference"]
    base_graph = [("R", "A"), ("A", "C")]
    for mode in ("override", "merge", "only"):
        report = apply_modpath_informed_graph_priors(
            hydrosheaf_edges=base_graph,
            modpath_edges=reference,  # type: ignore[arg-type]
            mode=mode,
            default_p_uv=0.95,
        )
        rows.append(
            {
                "validation_mode": report["validation_mode"],
                "prior_mode": mode,
                "n_input_hydrosheaf_edges": report["n_input_hydrosheaf_edges"],
                "n_modpath_prior_edges": report["n_modpath_prior_edges"],
                "n_output_edges": report["n_output_edges"],
                "not_independent_validation": report["not_independent_validation"],
                "output_edges": ";".join(report["output_edges"]),
            }
        )
    return rows


def write_summary(metric_rows: Sequence[Mapping[str, object]], prior_mode_rows: Sequence[Mapping[str, object]]) -> None:
    perfect = [row for row in metric_rows if float(row["f1"]) == 1.0]
    with_errors = [row for row in metric_rows if float(row["fp"]) > 0 or float(row["fn"]) > 0]
    scale_mismatch = [row for row in metric_rows if bool(row["scale_mismatch"])]
    lines = [
        "# M4 Topology Benchmark Results",
        "",
        f"- Independent graph scenarios tested: {len(metric_rows)}.",
        f"- Perfect reproduction scenarios: {len(perfect)}.",
        f"- Scenarios with false positives or false negatives: {len(with_errors)}.",
        f"- Scenarios with scale mismatch: {len(scale_mismatch)}.",
        f"- MODPATH-informed prior modes tested separately: {len(prior_mode_rows)}.",
        "",
        "Interpretation guardrail:",
        "",
        "Only the independent graph-inference rows evaluate Hydrosheaf topology skill",
        "against MODPATH advective connectivity. The MODPATH-informed prior rows show",
        "how MODPATH connectivity can enter Hydrosheaf as prior information and must",
        "not be reported as independent validation.",
        "",
    ]
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    (DOCS_DIR / "m4_results_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    for directory in (RESULTS_DIR, TABLES_DIR, DOCS_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    metric_rows, edge_rows = independent_rows()
    priors = prior_rows()
    write_csv(RESULTS_DIR / "independent_graph_vs_modpath.csv", metric_rows)
    write_csv(RESULTS_DIR / "edge_classification.csv", edge_rows)
    write_csv(RESULTS_DIR / "modpath_informed_priors.csv", priors)
    write_csv(TABLES_DIR / "table1_topology_validation_summary.csv", metric_rows)
    write_summary(metric_rows, priors)
    print(f"Wrote {len(metric_rows)} independent validation rows and {len(priors)} prior-mode rows.")


if __name__ == "__main__":
    main()
