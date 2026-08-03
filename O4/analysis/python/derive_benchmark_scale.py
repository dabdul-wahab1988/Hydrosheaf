"""TAB-4: design/case/bootstrap scale and field-/archive-transfer scope for
M6, M7 and M8. Counts are read directly from each component's own locked
DECISIONS.md / README.md / result-file row counts; no experiment is re-run to
produce this table (DECISIONS.md D2).

Run:  .venv/Scripts/python.exe O4/analysis/python/derive_benchmark_scale.py
"""

from __future__ import annotations

import pandas as pd

from _common import M6, M7, write

ROWS = [
    dict(
        component="M6 robustness", design_unit="wells x evidence tiers",
        primary_scale="160 wells x 5 tiers = 800 tier-well rows",
        replicate_or_bootstrap_count="not applicable (deterministic per-well classification, seed 1234)",
        external_reference="independent charge-balance QC; synthetic extended model with known truth (ED2)",
        field_or_archive_transfer="Northern Ghana (160 wells), Talensi (63 samples), Lower Anayari (41 samples)",
        source="results/m6_tier_ablation.csv; results/m6_dataset_readiness.csv",
    ),
    dict(
        component="M7 identifiability (M7.3)", design_unit="independent MODFLOW6/MODPATH7 cases",
        primary_scale="6 development cases + 12 untouched locked-test cases",
        replicate_or_bootstrap_count="50,000 age-importance particles per case/tracer regime; 64 reaction bootstraps per case; 10,000 case-block bootstrap replicates",
        external_reference="fresh independent MODFLOW 6/MODPATH 7 generator (imports no HydroSheaf code)",
        field_or_archive_transfer="Northern Ghana workbook audited for component-diagnostic readiness only, not topology/age validation",
        source="docs/m7_3_results.md",
    ),
    dict(
        component="M7 identifiability (M7.4 sheaf-vs-graph)", design_unit="independent scalar-section cases",
        primary_scale="64 held-out cases across 4 scenario strata",
        replicate_or_bootstrap_count="10,000 paired case-block bootstrap replicates",
        external_reference="independent sheaf/graph generator (imports no HydroSheaf code)",
        field_or_archive_transfer="none (controlled synthetic only)",
        source="results/RUN-M7-SHEAF-VS-GRAPH-20260729-01",
    ),
    dict(
        component="M7 identifiability (M7.5 robust hybrid)", design_unit="independent scalar-section cases",
        primary_scale="64 development cases (8401-8464) + 128 locked-test cases (8501-8628)",
        replicate_or_bootstrap_count="10,000 paired case-block bootstrap replicates",
        external_reference="independent sheaf/graph generator (imports no HydroSheaf code)",
        field_or_archive_transfer="none (controlled synthetic only)",
        source="results/RUN-M7-ROBUST-HYBRID-SHEAF-20260729-01",
    ),
    dict(
        component="M8 calibration (matched-model transport)", design_unit="calibrations (16 fixed designs + 1 OED experiment)",
        primary_scale="4,000 fixed-design + 4,500 optimal-design = 8,500 calibrations",
        replicate_or_bootstrap_count="250 paired replicates per design/strategy; 2,000-resample percentile bootstrap",
        external_reference="known synthetic truth (dispersivity 2.0 m, decay 0.005 1/d)",
        field_or_archive_transfer="none (controlled synthetic only)",
        source="manuscript/artifacts/m8_transport_parameter_summary.csv",
    ),
    dict(
        component="M8 calibration (independent-model robustness)", design_unit="calibrations against independently generated numerical truth",
        primary_scale="1,000 locked test calibrations (250 per strategy x 4 strategies)",
        replicate_or_bootstrap_count="80 development replicates (oracle selection); 5,000 paired bootstrap resamples",
        external_reference="implicit finite-volume/upwind numerical solver sharing no code with the calibration forward model",
        field_or_archive_transfer="none (controlled synthetic only)",
        source="results/RUN-M8-INDEPENDENT-20260728-01/strategy_summary.csv",
    ),
    dict(
        component="M8 calibration (kinetic structural control)", design_unit="PHREEQC kinetic designs",
        primary_scale="3 designs (single-time, multi-time, multi-time + surface area)",
        replicate_or_bootstrap_count="not applicable (deterministic local Fisher-information analysis)",
        external_reference="known synthetic truth (k=1e-10 mol/m^2/s, A=0.1 m^2/L)",
        field_or_archive_transfer="none (controlled synthetic only)",
        source="manuscript/artifacts/m8_kinetics_structural_summary.csv",
    ),
    dict(
        component="M8 calibration (frontier active learning)", design_unit="independent MODFLOW6/MODPATH7 aquifer cases",
        primary_scale="24 untouched locked-test cases, up to 5 sequential actions each",
        replicate_or_bootstrap_count="5,000 paired case-bootstrap resamples",
        external_reference="independent MODFLOW 6/MODPATH 7 heterogeneous-aquifer generator with nonlinear synthetic geochemistry",
        field_or_archive_transfer="none (controlled synthetic only)",
        source="provenance/runs/RUN-M8-FRONTIER-AL-20260728-01/strategy_summary.csv",
    ),
]


def main() -> None:
    write(pd.DataFrame(ROWS), "benchmark_scale.csv")


if __name__ == "__main__":
    main()
