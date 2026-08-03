# One-command O4 figure build. Run from the repository root.
scripts <- c(
  "O4/analysis/r/fig01_taxonomy.R",
  "O4/analysis/r/fig02_central_divergence.R",
  "O4/analysis/r/fig03_robustness_gradient.R",
  "O4/analysis/r/fig04_integration_value.R",
  "O4/analysis/r/fig05_calibration_gap.R",
  "O4/analysis/r/fig06_scale_and_scope.R"
)
for (s in scripts) {
  message("=== ", s, " ===")
  source(s)
}
message("All O4 figures written to O4/manuscript/artifacts/figures/")
