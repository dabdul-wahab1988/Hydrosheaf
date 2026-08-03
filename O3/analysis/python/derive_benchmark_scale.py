"""TAB-4: benchmark scale and computational cost, all three components.

Counts are read directly from each component's own manifests/summaries; where
a total is a sum across already-published per-file counts (e.g. M3's five
tracer-specific cross-validation files), the sum is disclosed as such.

Run:  .venv/Scripts/python.exe O3/analysis/python/derive_benchmark_scale.py
"""

from __future__ import annotations

import json

import pandas as pd

from _common import M3, M4, M5, write


def m3_scale() -> dict:
    design = pd.read_csv(M3 / "results/m3_design_matrix_summary.csv")
    n_scenarios = int(design["scenario_id"].nunique())
    n_rows_total = int(design["total_rows"].sum())
    cv_files = ["3H", "14C", "CFC11", "CFC12", "SF6"]
    n_cv_rows = 0
    for tracer in cv_files:
        p = M3 / f"results/m3_cv_benchmark_{tracer}.csv"
        if p.exists():
            n_cv_rows += len(pd.read_csv(p, usecols=[0]))
    return dict(
        component="Age / residence time",
        external_reference="USGS national groundwater-age release "
                           "(TracerLPM Table 4, reported configuration)",
        n_scenarios=n_scenarios,
        n_fits=n_rows_total,
        n_cross_validation_rows=n_cv_rows,
        median_runtime_per_fit_ms=pd.NA,
        runtime_recorded=False,
        note="Design matrix: 10 scenarios x 1,272 rows each; leakage-guarded "
            "tracer-withholding cross-validation adds 5 tracer-specific runs")


def m4_scale() -> dict:
    ind = pd.read_csv(M4 / "results/independent_graph_vs_modpath.csv")
    n_independent = int(len(ind))
    archives = 0
    n_edges = 0
    for tier in ["tier_1_savage", "tier_2_great_miami", "tier_3_long_island"]:
        p = M4 / f"results/{tier}_archive_summary.csv"
        if p.exists():
            row = pd.read_csv(p).iloc[0]
            if pd.notna(row.get("n_particles")) and row.get("n_particles", 0) > 0:
                archives += 1
                n_edges += int(row.get("n_endpoint_edges") or 0)
    return dict(
        component="Topology",
        external_reference="MODPATH particle-tracking connectivity "
                           "(3 public MODFLOW/MODPATH archives)",
        n_scenarios=n_independent,
        n_fits=n_edges,
        n_cross_validation_rows=pd.NA,
        median_runtime_per_fit_ms=pd.NA,
        runtime_recorded=False,
        note=f"{n_independent} independent graph-inference scenarios; "
            f"{archives} of 3 archives active (Long Island is a documented "
            "fallback stub with no active validation rows); edge counts are "
            "endpoint-derived reference edges, not per-fit runtime")


def m5_scale() -> dict:
    summary = json.loads((M5 / "results/analysis_summary.json").read_text())
    t5 = pd.read_csv(M5 / "tables/tableS5_complete_model_metrics.csv")
    runtime = t5["Median runtime (ms)"].median()
    return dict(
        component="Reaction",
        external_reference="240-scenario live-PHREEQC factorial synthetic "
                           "benchmark (carbonate/crystalline/evaporitic/mixed)",
        n_scenarios=int(summary["n_phreeqc_scenarios"]),
        n_fits=int(summary["n_benchmark_fits"]),
        n_cross_validation_rows=pd.NA,
        median_runtime_per_fit_ms=float(runtime),
        runtime_recorded=True,
        note="21,600 factorial inverse fits = 240 scenarios x 5 comparator "
            "methods x noise/panel/archetype combinations; only component "
            "with per-fit runtime recorded across all comparator methods")


def main() -> None:
    df = pd.DataFrame([m3_scale(), m4_scale(), m5_scale()])
    write(df, "benchmark_scale.csv")


if __name__ == "__main__":
    main()
