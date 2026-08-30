# Regenerate every manuscript figure from the read-only Python exports.
# Run from the repository root:
#   Rscript M2.3/analysis/r/make_all_figures.R

scripts <- c(
  "fig01_architecture.R",
  "fig02_field_setting.R",
  "fig03_data_availability.R",
  "fig04_synthetic_recovery.R",
  "fig05_programme_gates.R",
  "fig06_external_comparison.R",
  "fig07_field_application.R"
)

for (s in scripts) {
  message("== ", s)
  source(file.path("M2.3", "analysis", "r", s), local = new.env())
}
message("all figures regenerated")
