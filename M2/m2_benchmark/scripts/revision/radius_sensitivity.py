"""Graph-stage sensitivity to the configured edge search radius.

This is a topology-only audit for the revision response. It reruns the
documented probabilistic builder and the sheaf-inspired selection stage at
3, 5, 7.5 and 10 km for both Ghana sites. Chemistry fitting is deliberately
not rerun here, so the output is a sensitivity diagnostic for the retained
graph rather than a new field-performance estimate.

Output: M2/m2_benchmark/results/revision/radius_sensitivity.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[4]
BENCHMARK_ROOT = PROJECT_ROOT / "M2" / "m2_benchmark"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(BENCHMARK_ROOT / "scripts"))

from run_m2_field_documented_pipeline import (  # noqa: E402
    _config,
    _manu_samples,
    _talensi_samples,
)
from hydrosheaf.inference.network_fit import infer_edges  # noqa: E402
from hydrosheaf.sheaf.topology_refine import refine_edges_with_sheaf  # noqa: E402


RADII_KM = (3.0, 5.0, 7.5, 10.0)
OUT_PATH = BENCHMARK_ROOT / "results" / "revision" / "radius_sensitivity.csv"


def _topology(edges: list[object]) -> set[tuple[str, str]]:
    return {(str(edge.u), str(edge.v)) for edge in edges}


def main() -> None:
    rows: list[dict[str, object]] = []
    for site, samples in (("Manu", _manu_samples()), ("Talensi", _talensi_samples())):
        candidates_by_radius: dict[float, list[object]] = {}
        retained_by_radius: dict[float, list[object]] = {}
        topology_by_radius: dict[float, set[tuple[str, str]]] = {}
        for radius in RADII_KM:
            config = _config(site, samples, documented=True)
            config.edge_radius_km = radius
            candidates = infer_edges(samples, method="probabilistic", config=config)
            retained = refine_edges_with_sheaf(samples, candidates, config)
            candidates_by_radius[radius] = candidates
            retained_by_radius[radius] = retained
            topology_by_radius[radius] = _topology(retained)

        baseline = topology_by_radius[5.0]
        for radius in RADII_KM:
            topology = topology_by_radius[radius]
            union = topology | baseline
            jaccard = len(topology & baseline) / len(union) if union else 1.0
            retained = retained_by_radius[radius]
            rows.append(
                {
                    "site": site,
                    "radius_km": radius,
                    "n_candidates": len(candidates_by_radius[radius]),
                    "n_retained": len(retained),
                    "n_primary": sum(
                        getattr(edge, "type", "") == "primary" for edge in retained
                    ),
                    "n_lateral": sum(
                        getattr(edge, "type", "") == "lateral" for edge in retained
                    ),
                    "jaccard_vs_5km": round(jaccard, 3),
                }
            )

    output = pd.DataFrame(rows)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUT_PATH, index=False)
    print(output.to_string(index=False))
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
