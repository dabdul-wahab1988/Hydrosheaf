| Benchmark | External reference | Scenarios | Fits | Cross-validation rows | Median runtime per fit | Runtime recorded |
|---|---|---:|---:|---:|---:|:---:|
| Age / residence time | USGS national groundwater-age release (TracerLPM Table 4) | 13 design-matrix scenarios | 16,536 (10 scenarios x 1,272 rows, plus 3 ablations) | 375 (leakage-guarded tracer-withholding CV, 5 tracers) | not recorded | no |
| Topology | MODPATH particle-tracking connectivity (3 public MODFLOW/MODPATH archives) | 12 independent graph-inference rows | 242 endpoint-derived reference edges (Savage + Great Miami; Long Island is a documented fallback stub) | not applicable | not recorded | no |
| Reaction | 240-scenario live-PHREEQC factorial synthetic benchmark | 240 PHREEQC scenarios | 21,600 factorial inverse fits (240 x 5 comparator methods x noise/panel/archetype grid) | not applicable | 20.8 ms (median across methods) | yes |

Notes. Counts are read directly from each benchmark's own manifests
(`m3_design_matrix_summary.csv`, `independent_graph_vs_modpath.csv` plus the
three `tier_*_archive_summary.csv` files, `analysis_summary.json`); no
scenario was re-run to produce this table. The reaction benchmark is the
only one to record per-fit runtime across every comparator method, which
this paper reports as a genuine difference in benchmark design rather than
normalising it away.
