# Figure 3: active design-matrix reference agreement and reportability.
# Metrics use only identifiable reportable fits; archived withdrawn scenarios
# are excluded. Lower discrepancy does not imply field-validated true-age
# accuracy. Carried over from M3 Figure 2.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

scenario_labels <- c(
  tracerlpm_strict_parity              = "Strict reported-configuration\nparity",
  tracerlpm_parity_agefractions        = "Reported-output-constrained\nfraction sensitivity",
  hydrosheaf_selection_corrected       = "Hydrosheaf model\nselection",
  parity_reported_corrected            = "Reported-model\nparity",
  screened_dgm_gases                   = "Screened young-gas\ncorrection",
  oldwater_he4_uncertainty             = "Old-water 4He\nuncertainty",
  oldwater_c14_ensemble                = "Old-water 14C\nensemble",
  oldwater_ensemble_he4_uncertainty    = "Old-water ensemble\n4He uncertainty",
  ablation_no_he4                      = "Ablation:\nno 4He",
  ablation_raw_c14                     = "Ablation:\nraw 14C",
  tracer_young_only                    = "Young tracers\nonly"
)

perf <- read_export("Manuscript_Table2_Design_Matrix_Performance.csv") |>
  filter(scenario_id %in% names(scenario_labels), !is.na(median_abs_log10_error)) |>
  mutate(scenario_label = factor(scenario_labels[scenario_id],
                                 levels = rev(unname(scenario_labels[order(-median_abs_log10_error)]))))

p_error <- ggplot(perf, aes(x = median_abs_log10_error, y = scenario_label)) +
  geom_segment(aes(x = 0, xend = median_abs_log10_error,
                   y = scenario_label, yend = scenario_label),
               colour = GRID, linewidth = 2.2) +
  geom_point(aes(size = identifiable_rows), colour = PAL[1]) +
  scale_size_continuous(name = "Reportable N", range = c(1.4, 5)) +
  labs(tag = "(a)", title = "Median absolute log10 reference discrepancy",
       x = "Median |log10(estimated / reference)|", y = NULL) +
  theme_hs()

p_within2 <- ggplot(perf, aes(x = within_factor_2, y = scenario_label)) +
  geom_segment(aes(x = 0, xend = within_factor_2,
                   y = scenario_label, yend = scenario_label),
               colour = GRID, linewidth = 2.2) +
  geom_point(aes(size = identifiable_rows), colour = PAL[2]) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_size_continuous(name = "Reportable N", range = c(1.4, 5), guide = "none") +
  labs(tag = "(b)", title = "Share within a factor of two",
       x = "Within-factor-2 agreement", y = NULL) +
  theme_hs() +
  theme(axis.text.y = element_blank())

p <- p_error + p_within2 +
  plot_layout(widths = c(1.35, 1), guides = "collect") +
  plot_annotation(
    title = "Active design-matrix reference agreement and reportability",
    subtitle = str_wrap("Metrics use only identifiable reportable fits; withdrawn scenarios are excluded. Lower discrepancy does not imply field-validated true-age accuracy.", 100),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 8, colour = INK_MUTED))
  ) &
  theme(legend.position = "bottom")

save_fig(p, "FIG-3_design_matrix_performance", width_mm = 190, height_mm = 130)
