"""Revision analysis 5 (CAGEO-D-26-00847): candidate-edge construction scaling.

Answers Reviewer 2 Minor 2: report runtime and scalability of candidate-edge
construction (potential O(n^2) growth unless pruned).

We time the framework's coordinate-based edge inference
(infer_edges_from_coordinates, max_neighbors=3) on synthetic node sets of
increasing size and record wall time and candidate-edge counts.

Output: M2/m2_benchmark/results/revision/runtime_scaling.csv
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PROJECT_ROOT))

from hydrosheaf.graph.build import infer_edges_from_coordinates  # noqa: E402

OUT_DIR = PROJECT_ROOT / "M2" / "m2_benchmark" / "results" / "revision"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    rows = []
    for n in [10, 20, 40, 80, 160, 320]:
        lats = rng.uniform(10.0, 11.0, n)
        lons = rng.uniform(-2.0, -0.5, n)
        zs = rng.uniform(150.0, 300.0, n)
        samples = [
            {"site_id": f"N{i}", "lat": float(lats[i]), "lon": float(lons[i]), "elevation": float(zs[i])}
            for i in range(n)
        ]
        start = time.perf_counter()
        edges = infer_edges_from_coordinates(samples, max_neighbors=3, allow_uphill=False)
        elapsed = time.perf_counter() - start
        rows.append(
            {
                "n_nodes": n,
                "n_candidate_edges": len(edges),
                "wall_seconds": round(elapsed, 4),
                "edges_per_node": round(len(edges) / n, 2),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / "runtime_scaling.csv", index=False)
    print(out.to_string(index=False))
    print(f"Wrote {OUT_DIR / 'runtime_scaling.csv'}")


if __name__ == "__main__":
    main()
