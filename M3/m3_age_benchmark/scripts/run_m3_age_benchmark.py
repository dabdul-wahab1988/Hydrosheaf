"""Run M3 graph-regularized groundwater-age benchmarks.

The benchmark is deliberately synthetic and controlled. It tests when a graph
prior improves age inference, when it degrades inference, and whether randomized
negative-control graphs produce misleading improvements.
"""
from __future__ import annotations

import csv
import math
import random
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hydrosheaf.nuclear import audit_graph_age_coherence
from hydrosheaf.validation import regression_metrics


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
TABLES_DIR = BENCHMARK_DIR / "tables"
FIGURES_DIR = BENCHMARK_DIR / "figures"
DOCS_DIR = BENCHMARK_DIR / "docs"

PRIOR_STRENGTHS = {
    "none": 0.0,
    "weak": 0.25,
    "medium": 0.55,
    "strong": 0.85,
}


def rmse(true_ages: Mapping[str, float], estimated: Mapping[str, float]) -> float:
    metrics = regression_metrics(
        [true_ages[node] for node in true_ages],
        [estimated[node] for node in true_ages],
        log10=True,
    )
    return float(metrics["rmse"])


def count_violations(edges: Iterable[Tuple[str, str]], ages: Mapping[str, float]) -> int:
    audit = audit_graph_age_coherence(edges, {node: {"age_years": age} for node, age in ages.items()})
    return int(audit["n_violations"])


def graph_regularize(
    single_ages: Mapping[str, float],
    edges: Sequence[Tuple[str, str]],
    strength: float,
    *,
    iterations: int = 8,
    min_increment_years: float = 0.0,
) -> Dict[str, float]:
    """Apply a monotonic directed-edge prior to apparent ages."""

    ages = {node: float(age) for node, age in single_ages.items()}
    if strength <= 0:
        return ages
    alpha = min(max(float(strength), 0.0), 1.0)
    for _ in range(iterations):
        for upstream, downstream in edges:
            if upstream not in ages or downstream not in ages:
                continue
            required = ages[upstream] + min_increment_years
            if ages[downstream] < required:
                ages[downstream] = (1.0 - alpha) * ages[downstream] + alpha * required
    return ages


def randomized_edges(nodes: Sequence[str], edge_count: int, seed: int) -> List[Tuple[str, str]]:
    rng = random.Random(seed)
    candidates = [(a, b) for a in nodes for b in nodes if a != b]
    rng.shuffle(candidates)
    return candidates[:edge_count]


def scenarios() -> List[Dict[str, object]]:
    return [
        {
            "scenario": "helpful_graph",
            "purpose": "Correct graph fixes a downstream-younger-than-upstream tracer artefact.",
            "true_ages": {"A": 8.0, "B": 35.0, "C": 120.0, "D": 420.0},
            "single_ages": {"A": 30.0, "B": 14.0, "C": 105.0, "D": 390.0},
            "edges": [("A", "B"), ("B", "C"), ("C", "D")],
        },
        {
            "scenario": "degrading_graph",
            "purpose": "Wrong graph direction forces a physically false monotonic pattern.",
            "true_ages": {"A": 8.0, "B": 35.0, "C": 120.0, "D": 420.0},
            "single_ages": {"A": 9.0, "B": 31.0, "C": 130.0, "D": 405.0},
            "edges": [("D", "C"), ("C", "B"), ("B", "A")],
        },
        {
            "scenario": "mixed_recharge_ambiguity",
            "purpose": "A local young recharge node should not be over-smoothed by a regional graph.",
            "true_ages": {"A": 12.0, "B": 55.0, "C": 18.0, "D": 210.0},
            "single_ages": {"A": 11.0, "B": 60.0, "C": 20.0, "D": 240.0},
            "edges": [("A", "B"), ("B", "C"), ("C", "D")],
        },
    ]


def benchmark_rows() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for scenario in scenarios():
        name = str(scenario["scenario"])
        purpose = str(scenario["purpose"])
        true_ages = dict(scenario["true_ages"])  # type: ignore[arg-type]
        single_ages = dict(scenario["single_ages"])  # type: ignore[arg-type]
        edges = list(scenario["edges"])  # type: ignore[arg-type]
        single_rmse = rmse(true_ages, single_ages)
        before_violations = count_violations(edges, single_ages)

        graph_sets = [("scenario_graph", edges)]
        graph_sets.append(
            (
                "randomized_negative_control",
                randomized_edges(list(true_ages), len(edges), seed=100 + len(rows)),
            )
        )

        for graph_label, graph_edges in graph_sets:
            for prior_label, strength in PRIOR_STRENGTHS.items():
                graph_ages = graph_regularize(single_ages, graph_edges, strength)
                graph_rmse = rmse(true_ages, graph_ages)
                after_violations = count_violations(graph_edges, graph_ages)
                rows.append(
                    {
                        "scenario": name,
                        "purpose": purpose,
                        "graph_label": graph_label,
                        "prior_strength": prior_label,
                        "prior_weight": strength,
                        "rmse_single_log10": single_rmse,
                        "rmse_graph_log10": graph_rmse,
                        "delta_rmse_graph_minus_single": graph_rmse - single_rmse,
                        "improved": graph_rmse < single_rmse,
                        "n_violations_before": before_violations,
                        "n_violations_after": after_violations,
                        "edges": ";".join(f"{u}->{v}" for u, v in graph_edges),
                    }
                )
    return rows


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows: Sequence[Mapping[str, object]]) -> None:
    summary_rows = [
        row
        for row in rows
        if row["graph_label"] == "scenario_graph" and row["prior_strength"] != "none"
    ]
    negative_rows = [
        row
        for row in rows
        if row["graph_label"] == "randomized_negative_control" and row["prior_strength"] != "none"
    ]
    out_rows = []
    for row in summary_rows:
        out_rows.append(
            {
                "scenario": row["scenario"],
                "prior_strength": row["prior_strength"],
                "rmse_single_log10": row["rmse_single_log10"],
                "rmse_graph_log10": row["rmse_graph_log10"],
                "delta_rmse_graph_minus_single": row["delta_rmse_graph_minus_single"],
                "improved": row["improved"],
                "n_violations_after": row["n_violations_after"],
            }
        )
    write_csv(TABLES_DIR / "table1_m3_benchmark_summary.csv", out_rows)

    improved = [row for row in summary_rows if bool(row["improved"])]
    degraded = [row for row in summary_rows if not bool(row["improved"])]
    neg_improved = [row for row in negative_rows if bool(row["improved"])]
    lines = [
        "# M3 Age Benchmark Results",
        "",
        f"- Scenario-graph rows tested: {len(summary_rows)}.",
        f"- Rows where graph regularization improved RMSE: {len(improved)}.",
        f"- Rows where graph regularization degraded or failed to improve RMSE: {len(degraded)}.",
        f"- Randomized negative-control rows with apparent improvement: {len(neg_improved)}.",
        "",
        "Interpretation guardrail:",
        "",
        "Hydrosheaf should claim graph regularization improves age inference only for the",
        "specific tracer/graph/prior settings where RMSE decreases and randomized graphs",
        "do not produce comparable improvement. Degrading scenarios must be reported as",
        "evidence that graph priors can harm inference when topology, recharge structure,",
        "or tracer contamination assumptions are wrong.",
        "",
    ]
    (DOCS_DIR / "m3_results_summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_figure(rows: Sequence[Mapping[str, object]]) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    rows = [
        row
        for row in rows
        if row["graph_label"] == "scenario_graph" and row["prior_strength"] != "none"
    ]
    labels = [f"{row['scenario']}\n{row['prior_strength']}" for row in rows]
    values = [float(row["delta_rmse_graph_minus_single"]) for row in rows]
    colors = ["#2f855a" if value < 0 else "#c53030" for value in values]
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.axhline(0.0, color="#444444", linewidth=1)
    ax.bar(range(len(values)), values, color=colors)
    ax.set_ylabel("Delta RMSE, graph - single (log10 years)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_title("M3 graph regularization benchmark")
    fig.tight_layout()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES_DIR / "figure1_graph_regularization_delta_rmse.png", dpi=200)
    plt.close(fig)


def main() -> None:
    for directory in (RESULTS_DIR, TABLES_DIR, FIGURES_DIR, DOCS_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    rows = benchmark_rows()
    write_csv(RESULTS_DIR / "graph_regularization_scenarios.csv", rows)
    write_summary(rows)
    write_figure(rows)
    print(f"Wrote {len(rows)} benchmark rows to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
