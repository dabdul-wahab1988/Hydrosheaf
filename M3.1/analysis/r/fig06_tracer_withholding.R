# Figure 6: leakage-guarded target-tracer withholding. Each target is removed
# before every single-node and graph-node refit, and only reportable fits are
# scored in native units. Candidate graphs show no meaningful gain. Carried
# over from M3 Figure 6.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

# CFC11 and CFC12 carry zero eligible graph edges (non-estimable graph
# comparison; see Supplementary Figure S3), so the main figure keeps the
# three tracers where a graph-vs-baseline comparison is actually possible.
cv <- read_export("tracer_withholding_summary.csv") |>
  filter(tracer %in% c("3H", "SF6", "14C"))

variant_labels <- c(
  rmse_single_baseline    = "Single-node\nbaseline",
  rmse_hydraulic_graph    = "Hydraulic\ngraph",
  rmse_depth_graph        = "Depth\ngraph",
  rmse_randomised_control = "Randomised\ncontrol"
)

long <- cv |>
  tidyr::pivot_longer(cols = starts_with("rmse_"), names_to = "variant", values_to = "rmse") |>
  mutate(
    variant_label = factor(variant_labels[variant], levels = unname(variant_labels)),
    tracer_label = sprintf("%s (N = %d/%d)", tracer, reportable_rows, eligible_rows)
  )
tracer_order <- cv |> arrange(tracer) |> pull(tracer)
long$tracer_label <- factor(long$tracer_label, levels = unique(long$tracer_label[order(match(long$tracer, tracer_order))]))

p <- ggplot(long, aes(x = variant_label, y = rmse, fill = variant_label == "Single-node\nbaseline")) +
  geom_col(width = 0.62) +
  geom_text(aes(label = scales::number(rmse, accuracy = 0.001)),
            vjust = -0.4, size = 2.1, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = "#5F6368", `FALSE` = PAL[1]), guide = "none") +
  facet_wrap(~tracer_label, scales = "free_y", nrow = 1) +
  expand_limits(y = 0) +
  labs(
    title = "Leakage-guarded target-tracer withholding",
    subtitle = str_wrap("The withheld tracer is excluded from every fit and from graph construction; only reportable fits are scored, in native units. No candidate graph shows a meaningful gain over the single-node baseline.", 110),
    x = NULL, y = "RMSE (native tracer unit)"
  ) +
  theme_hs() +
  theme(axis.text.x = element_text(size = 6), strip.text = element_text(size = 7))

save_fig(p, "FIG-6_tracer_withholding", width_mm = 190, height_mm = 100)
