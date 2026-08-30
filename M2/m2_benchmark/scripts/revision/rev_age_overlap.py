"""Revision analysis 4 (CAGEO-D-26-00847): wide age-range handling in the age-ordering constraint.

Answers Reviewer 1 Major 12: how the framework handles wide, overlapping
posterior age ranges in the age-ordering criterion tau_j >= tau_i - eps.

We reproduce the benchmark's posterior model (evaluate_age_inference) and run
the actual audit_graph_age_coherence implementation over the true directed
edges, then quantify:
  - fraction of edges violating the ordering criterion,
  - among violations, the share whose posterior intervals OVERLAP
    (retained, flagged "not resolved at stated uncertainty"),
  - the share of SEVERE violations (non-overlapping reversals beyond the
    log10 severity threshold), which receive age-coherence failure codes.

Output: M2/m2_benchmark/results/revision/age_overlap_stats.csv
"""
from __future__ import annotations

import sys
import math
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
BENCHMARK_ROOT = PROJECT_ROOT / "M2" / "m2_benchmark"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(BENCHMARK_ROOT / "scripts"))

from run_m2_benchmark import build_noiseless_nodes, load_truth, make_config  # noqa: E402
from hydrosheaf.nuclear.age_coherence import audit_graph_age_coherence  # noqa: E402

OUT_DIR = BENCHMARK_ROOT / "results" / "revision"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    truth = load_truth(BENCHMARK_ROOT / "config" / "ground_truth.yaml")
    config = make_config(truth)
    nodes = build_noiseless_nodes(truth, config)
    seed = int(truth["benchmark"]["random_seed"])

    ordered_edges = list(truth["generation_edges"]) + list(truth.get("lateral_truth_edges", []))
    edge_pairs = [(str(e["u"]), str(e["v"])) for e in ordered_edges]

    stats = []
    for realisation in range(100):
        rng = np.random.default_rng(seed + 3000 + realisation)
        records = {}
        for node_id, node in nodes.items():
            true_age = float(node["mrt_years"])
            age_class = str(node["age_class"])
            sigma = 0.25 if age_class == "young" else 0.45
            if age_class in {"old", "fossil"}:
                sigma = 0.65
            if age_class == "mixed":
                sigma = 1.10
            network = true_age * math.exp(rng.normal(0.0, sigma * 0.48))
            width = 3.92 * sigma * 0.48 * true_age
            records[str(node_id)] = {
                "age_years": float(network),
                "ci_low_years": float(max(0.0, network - width / 2.0)),
                "ci_high_years": float(network + width / 2.0),
                "flags": {},
            }
        report = audit_graph_age_coherence(edge_pairs, records)
        n_checked = int(report["n_checked"])
        n_violations = int(report["n_violations"])
        n_severe = int(report["n_severe_violations"])
        stats.append(
            {
                "realisation": realisation,
                "n_checked": n_checked,
                "n_violations": n_violations,
                "n_severe_violations": n_severe,
                "violation_fraction": float(report["violation_fraction"]),
                "severe_fraction": (n_severe / n_checked) if n_checked else float("nan"),
            }
        )
    out = pd.DataFrame(stats)
    summary = pd.DataFrame(
        [
            {"metric": "n_edges_checked", "value": int(out["n_checked"].median())},
            {"metric": "violation_fraction_mean", "value": float(out["violation_fraction"].mean())},
            {"metric": "severe_violation_fraction_mean", "value": float(out["severe_fraction"].mean())},
            {"metric": "violations_with_overlapping_intervals_fraction", "value": float(
                1.0 - out["severe_fraction"].mean() / max(out["violation_fraction"].mean(), 1e-12)
            )},
            {"metric": "age_order_consistency_index", "value": float(1.0 - out["violation_fraction"].mean())},
        ]
    )
    summary.to_csv(OUT_DIR / "age_overlap_stats.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote {OUT_DIR / 'age_overlap_stats.csv'}")


if __name__ == "__main__":
    main()
