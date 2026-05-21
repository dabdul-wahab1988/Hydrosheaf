# M4 Topology Benchmark Manifest

Generated outputs:

- `results/independent_graph_vs_modpath.csv`: scenario metrics for independent
  Hydrosheaf graph inference against MODPATH reference connectivity.
- `results/edge_classification.csv`: TP/FP/FN/TN edge-level classifications.
- `results/modpath_informed_priors.csv`: MODPATH-informed prior-mode output;
  these rows are explicitly not independent validation.
- `results/external_modpath_archive_summary.csv`: M4 ingestion of the M2 USGS
  Savage MODFLOW/MODPATH external archive validation.
- `results/external_modpath_edge_agreement.csv`: edge-level endpoint/pathline
  agreement from the public MODPATH archive.
- `results/external_modpath_pathline_structure.csv`: particle-level raw pathline
  sequence structure, endpoint projection checks, and elapsed-time summaries.
- `results/external_modpath_capture_envelope_overlap.csv`: receptor-level
  point-cloud capture-envelope overlap from endpoint release coordinates and
  pathline first-point coordinates.
- `results/external_modpath_capture_envelope_summary.csv`: summary of
  capture-envelope IoU and source-cell Jaccard overlap.
- `results/external_modpath_travel_time_rank.csv`: edge-level endpoint-time versus
  pathline-elapsed-time rank diagnostic.
- `results/external_modpath_harmonized_travel_time.csv`: endpoint total-time
  versus MODPATH-derived Hydrosheaf edge-weight comparison.
- `results/external_modpath_harmonized_travel_time_summary.csv`: harmonised
  travel-time rank and scale agreement summary.
- `results/external_modpath_bootstrap_ci.csv`: bootstrap intervals for particle
  and edge diagnostics.
- `tables/table1_topology_validation_summary.csv`: manuscript-ready summary of
  precision, recall, F1, false-positive rate, and false-negative rate.
- `tables/Manuscript_Ready/Manuscript_Table1_M4_Benchmark_Design.csv`: analysis
  contract separating implemented, independent, prior-informed, and future
  evidence blocks.
- `tables/Manuscript_Ready/Manuscript_Table2_Independent_Topology_Metrics.csv`:
  Q1-facing independent topology metrics with claim guardrails.
- `tables/Manuscript_Ready/Manuscript_Table3_Edge_Failure_Diagnostics.csv`:
  true-positive, false-positive, false-negative, and scale-mismatch diagnostics.
- `tables/Manuscript_Ready/Manuscript_Table4_MODPATH_Prior_Mode_Audit.csv`:
  prior-mode audit showing that MODPATH-informed rows are not independent
  validation.
- `tables/Manuscript_Ready/Manuscript_Table5_External_MODPATH_Archive_Validation.csv`:
  public USGS Savage archive validation summary.
- `tables/Manuscript_Ready/Manuscript_Table6_External_Pathline_Time_Diagnostics.csv`:
  full pathline sequence, endpoint projection, bootstrap, and travel-time rank
  diagnostics.
- `tables/Manuscript_Ready/Manuscript_Table7_Capture_Travel_Time_Validation.csv`:
  point-cloud capture-envelope and harmonised travel-time validation.
- `tables/Manuscript_Ready/Manuscript_Table8_MODPATH_Endpoint_Scenarios.csv`:
  endpoint-derived MODPATH validation summary when generated.
- `tables/Manuscript_Ready/Manuscript_Table9_Sparsity_Sensitivity.csv`:
  controlled node-sparsity sensitivity when generated.
- `tables/Manuscript_Ready/Manuscript_Table10_Reproducibility.csv`: Python and
  package version record for the M4 analysis scripts.
- `figures/Manuscript_Ready/Manuscript_Fig1_M4_Benchmark_Workflow.png` and
  `.svg`: workflow schematic separating validation from prior use.
- `figures/Manuscript_Ready/Manuscript_Fig2_M4_Independent_Topology_Performance.png`
  and `.svg`: multi-panel independent topology performance figure.
- `figures/Manuscript_Ready/Manuscript_Fig3_M4_Edge_Failure_Networks.png` and
  `.svg`: edge-level failure diagnostic network panels.
- `figures/Manuscript_Ready/Manuscript_Fig4_M4_MODPATH_Endpoint_Validation.png`
  and `.svg`: endpoint-derived MODPATH and sparsity-sensitivity figure when
  generated.
- `figures/Manuscript_Ready/Manuscript_Fig5_M4_External_MODPATH_Archive_Validation.png`
  and `.svg`: public USGS Savage archive topology, pathline, and time-rank
  diagnostic figure.
- `figures/Manuscript_Ready/Manuscript_Fig6_M4_Capture_Travel_Time_Validation.png`
  and `.svg`: point-cloud capture-envelope overlap and harmonised endpoint-time
  agreement with MODPATH-derived Hydrosheaf edge weights.
- `figures/Manuscript_Ready/Manuscript_EDFig1_M4_MODPATH_Prior_Mode_Audit.png`
  and `.svg`: extended-data prior-mode audit.
- `docs/m4_results_summary.md`: narrative guardrails for M4 claims.
- `docs/02_figures.md`: figure evidence register.
- `docs/03_tables.md`: table evidence register.

Claim rule:

M4 can claim reduced-order topology reproduction only for independent graph
inference rows. MODPATH-informed priors should be described as a prior-assisted
Hydrosheaf mode, not as independent evidence of topology inference skill.

Pending evidence:

- model-supplied polygon/raster capture-zone overlap;
- independent travel-time prediction outside MODPATH-derived prior transfer;
- processed Great Miami, Long Island, and coastal external-validation archives.
