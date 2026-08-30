# Supplementary Figure S3: CFC-11 and CFC-12 target-withheld RMSE after the
# reportability guard. No reportable CFC nodes formed eligible graph edges,
# so identical bars indicate a non-estimable graph effect, not equivalence.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

variant_labels <- c(
  rmse_single_baseline    = "Single-node\nbaseline",
  rmse_hydraulic_graph    = "Hydraulic\ngraph",
  rmse_depth_graph        = "Depth\ngraph",
  rmse_randomised_control = "Randomised\ncontrol"
)

cv <- read_export("tracer_withholding_summary.csv") |>
  filter(tracer %in% c("CFC11", "CFC12"))

long <- cv |>
  tidyr::pivot_longer(cols = starts_with("rmse_"), names_to = "variant", values_to = "rmse") |>
  mutate(
    variant_label = factor(variant_labels[variant], levels = unname(variant_labels)),
    tracer_label = sprintf("%s (N = %d/%d)", tracer, reportable_rows, eligible_rows)
  )

p <- ggplot(long, aes(x = variant_label, y = rmse)) +
  geom_col(width = 0.6, fill = PAL[4]) +
  geom_text(aes(label = scales::number(rmse, accuracy = 0.001)),
            vjust = -0.4, size = 2.2, colour = INK) +
  facet_wrap(~tracer_label, scales = "free_y") +
  expand_limits(y = 0) +
  labs(
    title = "CFC-11 and CFC-12 target-withheld RMSE",
    subtitle = str_wrap("No reportable CFC node formed an eligible graph edge, so identical values across variants are a non-estimable comparison, not evidence of equivalence.", 100),
    x = NULL, y = "RMSE (pptv)"
  ) +
  theme_hs()

save_fig(p, "FIG-S3_cfc_withholding", width_mm = 140, height_mm = 90)
