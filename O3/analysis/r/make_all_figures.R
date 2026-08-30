# Regenerate every O3 figure from the Python exports. Run from the
# repository root.

figs <- c("fig01_pipelines.R", "fig02_benchmark_scale.R", "fig03_headline_pattern.R",
          "fig04_calibration_gap.R", "fig05_within_component.R", "fig06_field_transfer.R")

for (f in figs) {
  message("=== ", f, " ===")
  source(file.path("O3", "analysis", "r", f))
}
message("All O3 figures regenerated.")
