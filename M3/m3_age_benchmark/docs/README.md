# M3 documentation status

The authoritative current claim boundary is defined by:

- `../DECISIONS.md`;
- `../results/m3_manuscript_analysis_manifest.json`;
- `../results/m3_design_matrix_manifest.json`;
- `../results/m3_real_usgs_graph_benchmark_manifest.json`;
- the five `m3_cv_benchmark_*_manifest.json` files; and
- the `../tables/Manuscript_Ready` and `../figures/Manuscript_Ready` outputs.

QA files dated before the 2026-07-28 accuracy lock are historical diagnostics.
Their counts and metrics may predate the reference-free reportability guard,
corrected tracer-withholding leakage control, current input forcing, or scenario
withdrawals. They must not be quoted as current manuscript results unless their
values are reproduced in the authoritative manifests and current tables.

In particular:

- the hierarchical old-water prior is withdrawn;
- the age-fraction result is a shared-provenance sensitivity, not independent
  validation;
- the MRVA release is not a like-for-like graph-age replication because the
  required reported-LPM table is unavailable;
- the gas-correction effect is non-estimable in the audited paired fields; and
- no candidate graph gives a robust reference-agreement or target-withheld
  predictive improvement in the current public-data benchmark.
