# Render every M3.1 figure and map. Run from the repository root:
#   Rscript M3.1/analysis/r/make_all_figures.R
# Requires manuscript/artifacts/data/*.csv (M3.1/analysis/python/build_exports.py)
# and, for FIG-1, M3.1/analysis/r/data/study_area_boundaries.gpkg
# (M3.1/analysis/r/fetch_boundaries.R, one-time, needs network).

scripts <- c(
  "fig01_site_map.R",
  "fig02_atmospheric_histories.R",
  "fig03_design_matrix_performance.R",
  "fig04_strict_parity_scatter.R",
  "fig05_graph_benchmark.R",
  "fig06_tracer_withholding.R",
  "fig07_infeasibility_audit.R",
  "figS2_edge_geometry.R",
  "figS3_cfc_withholding.R",
  "figS4_network_dating_demo.R"
)

for (s in scripts) {
  message("=== ", s, " ===")
  source(file.path("M3.1", "analysis", "r", s), local = new.env())
}
message("All M3.1 figures written to M3.1/manuscript/artifacts/figures/")
