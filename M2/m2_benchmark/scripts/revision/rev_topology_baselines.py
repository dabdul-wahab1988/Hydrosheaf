"""Revision analysis 2 (CAGEO-D-26-00847): topology baseline comparisons.

Answers Reviewer 1 Major 11 and Reviewer 2 Major 2: compare the no-prior
topology inference against simple baseline graph-construction rules using the
same Savage MODPATH reference (174 directed edges, DOI 10.5066/F7J102FK):

  1. head-gradient (elevation-as-head, downhill, k=2)  -- the canonical no-prior result
  2. elevation-drop (downhill from strictly higher node)
  3. distance-based k-nearest (k=2, symmetric)
  4. conservative-tracer ordering (Cl ascending along flow; screening-level)

Output: M2/m2_benchmark/results/revision/topology_baselines.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.graph.build import infer_edges_from_coordinates  # noqa: E402
from hydrosheaf.validation import validate_independent_graph_against_modpath  # noqa: E402

SAVAGE_RESULTS = PROJECT_ROOT / "M4" / "m4_topology_benchmark" / "results" / "public_archives" / "savage"
OUT_DIR = PROJECT_ROOT / "M2" / "m2_benchmark" / "results" / "revision"


def nodes_to_samples(nodes: pd.DataFrame) -> list[dict]:
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
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nodes = pd.read_csv(SAVAGE_RESULTS / "modpath_node_mapping.csv")
    ref = pd.read_csv(SAVAGE_RESULTS / "modpath_reference_edges.csv")
    ref_edges = [(str(u), str(v)) for u, v in zip(ref["u"], ref["v"]) if u != v]
    node_ids = sorted(str(n) for n in nodes["node_id"])
    candidate_universe = [(u, v) for u, v in [(a, b) for a in node_ids for b in node_ids] if u != v]
    ref_set = set(ref_edges)

    def score(inferred: list[tuple[str, str]]) -> dict:
        report = validate_independent_graph_against_modpath(inferred, ref_edges, candidate_edges=candidate_universe)
        m = report["metrics"]
        return {
            "n_reference_edges": m["n_reference_edges"],
            "n_inferred_edges": m["n_inferred_edges"],
            "tp": m["tp"], "fp": m["fp"], "fn": m["fn"],
            "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
        }

    rows = []

    # 1. head-gradient (canonical no-prior): elevation-as-head, downhill, k=2
    inferred = [(e.u, e.v) for e in infer_edges_from_coordinates(nodes_to_samples(nodes), max_neighbors=2, allow_uphill=False)]
    rows.append({"baseline": "head_gradient_downhill_k2", "independent": True, **score(inferred)})

    # Coordinate/elevation conventions follow the framework's inference path:
    # x/y are projected metres; elevation uses the archive's z with NaN imputed
    # to 0.0 (nodes_to_samples), which is what the canonical construction sees.
    pos = {str(r["node_id"]): (float(r["x"]), float(r["y"]), float(r["z"]) if not pd.isna(r["z"]) else 0.0) for _, r in nodes.iterrows()}

    # 2. distance k-nearest (k=2, no direction filter): proximity-only baseline
    def euclid(a, b):
        return float(np.hypot(a[0] - b[0], a[1] - b[1]))

    dist_edges = []
    for u in node_ids:
        d = sorted((euclid(pos[u], pos[v]), v) for v in node_ids if v != u)
        for _, v in d[:2]:
            dist_edges.append((u, v))
    rows.append({"baseline": "distance_knn_k2", "independent": True, **score(dist_edges)})

    # 3. elevation-drop over all downhill pairs (NaN z imputed to 0.0,
    #    identical convention to the framework's downhill filter)
    elev_down = [(u, v) for u in node_ids for v in node_ids if u != v and pos[u][2] > pos[v][2]]
    rows.append({"baseline": "elevation_drop_all_pairs", "independent": True, **score(elev_down)})

    # 4. conservative-tracer ordering (Cl) -- screening-level: no field Cl data in
    #    the Savage archive, so this baseline is reported as not computable here
    #    and the framework's conservative-tracer evidence is documented in the SI.
    rows.append({
        "baseline": "conservative_tracer_ordering",
        "independent": True,
        "n_reference_edges": len(ref_edges),
        "n_inferred_edges": np.nan,
        "tp": np.nan, "fp": np.nan, "fn": np.nan,
        "precision": np.nan, "recall": np.nan, "f1": np.nan,
    })

    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "topology_baselines.csv", index=False)
    print(out.to_string(index=False))
    print(f"Wrote {OUT_DIR / 'topology_baselines.csv'}")


if __name__ == "__main__":
    main()
