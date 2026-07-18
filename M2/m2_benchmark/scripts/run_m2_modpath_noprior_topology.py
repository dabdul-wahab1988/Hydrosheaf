"""M2 no-prior (heuristic-only) topology validation against the Savage MODPATH reference.

PROVENANCE CORRECTION SCRIPT.

The submitted M2 manuscript originally reported a no-prior topology result of
215 candidate edges / 168 TP / 47 FP / 6 FN (precision 0.78, recall 0.96,
F1 0.86) in Section 3.4. A full provenance audit (2026-07-10) found that this
figure exists in no result file, is producible by no code path in this
repository or its git history, and exceeds the information limit of the
benchmark (the outlet-convergent Savage structure caps honest recovery at
147 TP). It is treated as an unreconciled drafting error.

This script computes the CORRECT, reproducible no-prior result by rerunning
the exact independent head-gradient inference used by the dedicated M4
topology benchmark (phase2b_independent_validation.py) on the same
standardized Savage inputs:

  - nodes:  M4/.../results/public_archives/savage/modpath_node_mapping.csv
  - edges:  M4/.../results/public_archives/savage/modpath_reference_edges.csv
    (both derived from the public USGS archive, DOI 10.5066/F7J102FK)

Expected output (must match M4/results/independent_graph_vs_modpath.csv):

  no-prior head-gradient: 302 candidates, TP=147, FP=155, FN=27
                          precision=0.487, recall=0.845, F1=0.618

The prior-assisted row (TP=174, FP=0, FN=0, F1=1.00) is carried over from
M2's E2 report unchanged; it validates MODPATH prior INGESTION fidelity only,
not independent inference.

Writes: M2/m2_benchmark/results/modpath_noprior_topology.csv
        M2/m2_benchmark/external/modpath/results/e2c_noprior_topology_report.md
"""
from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.graph.build import infer_edges_from_coordinates
from hydrosheaf.validation import validate_independent_graph_against_modpath

SAVAGE_RESULTS = (
    PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results"
    / "public_archives" / "savage"
)
M2_RESULTS = PROJECT_ROOT / "M2" / "m2_benchmark" / "results"
M2_MODPATH_DOCS = (
    PROJECT_ROOT / "M2" / "m2_benchmark" / "external" / "modpath" / "results"
)


def nodes_to_samples(nodes: pd.DataFrame) -> list[dict]:
    """Grid-cell nodes -> Hydrosheaf sample dicts (same mapping as M4 phase2b)."""
    samples = []
    for _, row in nodes.iterrows():
        samples.append({
            "site_id": str(row["node_id"]),
            "lat": float(row["y"]),
            "lon": float(row["x"]),
            "elevation": float(row["z"]) if not pd.isna(row["z"]) else 0.0,
        })
    return samples


def main() -> None:
    nodes = pd.read_csv(SAVAGE_RESULTS / "modpath_node_mapping.csv")
    ref = pd.read_csv(SAVAGE_RESULTS / "modpath_reference_edges.csv")
    ref_edges = [(str(u), str(v)) for u, v in zip(ref["u"], ref["v"]) if u != v]

    samples = nodes_to_samples(nodes)
    node_ids = sorted(str(n) for n in nodes["node_id"])
    candidate_universe = [(u, v) for u in node_ids for v in node_ids if u != v]

    # Independent no-prior inference: elevation-as-head, downhill only, k=2.
    # Identical call to M4 phase2b_independent_validation head_gradient scenario.
    inferred_obj = infer_edges_from_coordinates(
        samples, max_neighbors=2, allow_uphill=False
    )
    inferred = [(e.u, e.v) for e in inferred_obj]

    report = validate_independent_graph_against_modpath(
        inferred, ref_edges, candidate_edges=candidate_universe
    )
    m = report["metrics"]

    rows = [
        {
            "mode": "no_prior_head_gradient",
            "independent_validation": True,
            "n_reference_edges": m["n_reference_edges"],
            "n_inferred_edges": m["n_inferred_edges"],
            "tp": m["tp"], "fp": m["fp"], "fn": m["fn"],
            "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
            "interpretation": (
                "Independent heuristic inference (elevation-as-head, downhill, "
                "k=2). High recall with substantial overconnection: FP exceeds "
                "TP, so inferred edges are screening-level hypotheses."
            ),
        },
        {
            "mode": "prior_assisted_ingestion",
            "independent_validation": False,
            "n_reference_edges": 174,
            "n_inferred_edges": 174,
            "tp": 174, "fp": 0, "fn": 0,
            "precision": 1.0, "recall": 1.0, "f1": 1.0,
            "interpretation": (
                "MODPATH physics-prior ingestion fidelity only (E2 report); "
                "not independent topology inference."
            ),
        },
    ]
    out = pd.DataFrame(rows)
    M2_RESULTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(M2_RESULTS / "modpath_noprior_topology.csv", index=False)

    # Cross-check against the M4 audited result.
    m4_csv = (
        PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results"
        / "independent_graph_vs_modpath.csv"
    )
    match_note = "M4 cross-check not found"
    if m4_csv.exists():
        m4 = pd.read_csv(m4_csv)
        hg = m4[m4["scenario"] == "head_gradient"].iloc[0]
        agree = (
            int(hg["tp"]) == m["tp"] and int(hg["fp"]) == m["fp"]
            and int(hg["fn"]) == m["fn"]
        )
        match_note = (
            "MATCHES M4 head_gradient exactly" if agree
            else f"MISMATCH vs M4 (M4: tp={hg['tp']}, fp={hg['fp']}, fn={hg['fn']})"
        )

    stamp = datetime.now(timezone.utc).isoformat()
    M2_MODPATH_DOCS.mkdir(parents=True, exist_ok=True)
    md = f"""# E2c No-Prior Topology Validation Report (Provenance Correction)

Run timestamp UTC: {stamp}

Source: USGS Savage MODFLOW/MODPATH archive, DOI `10.5066/F7J102FK`
(standardized node mapping and 174-edge directed reference shared with the
M4 topology benchmark).

## Why this report exists

The submitted M2 manuscript (Section 3.4) originally reported a no-prior
topology result of 215 candidates / 168 TP / 47 FP / 6 FN (F1 = 0.86). A
provenance audit found no computational source for those numbers in this
repository or its history. This report supersedes them with the reproducible
independent result below, identical in method to the dedicated M4 benchmark.

## Result

| mode | independent | candidates | TP | FP | FN | precision | recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| no-prior head-gradient | yes | {m['n_inferred_edges']} | {m['tp']} | {m['fp']} | {m['fn']} | {m['precision']:.3f} | {m['recall']:.3f} | {m['f1']:.3f} |
| prior-assisted ingestion | no (prior transfer) | 174 | 174 | 0 | 0 | 1.000 | 1.000 | 1.000 |

Cross-check vs M4 `independent_graph_vs_modpath.csv` (head_gradient): **{match_note}**.

## Interpretation guardrails

- The F1 = 1.00 row validates MODPATH physics-prior ingestion fidelity only.
- The no-prior result is high recall with substantial overconnection
  (FP > TP): independently inferred edges are screening-level candidate
  hypotheses requiring corroboration (chemistry residuals, tracer evidence,
  or MODPATH agreement) before process interpretation.
- The recall ceiling (147/174) is an information limit of the benchmark:
  TP and FP candidates converge on the same outlet cells, and the
  distinction depends on internal model permeability/pumping that
  coordinates and head proxies cannot carry.
"""
    (M2_MODPATH_DOCS / "e2c_noprior_topology_report.md").write_text(
        md, encoding="utf-8"
    )

    print(out.to_string(index=False))
    print(f"\nCross-check: {match_note}")
    print(f"Wrote {M2_RESULTS / 'modpath_noprior_topology.csv'}")
    print(f"Wrote {M2_MODPATH_DOCS / 'e2c_noprior_topology_report.md'}")


if __name__ == "__main__":
    main()
