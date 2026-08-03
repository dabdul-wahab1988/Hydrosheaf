# Figure 6: comparisons against evidence generated outside this study.
#
# Panel (a) contrasts held-out agreement with published lumped-parameter age
# outputs against the calibrated emulation of the same outputs. The distinction
# is the point of the panel: calibrating to a reference measures emulation, not
# independent agreement (DECISION D4).
# Panel (b) reports the no-prior directed-topology comparison, whose recall far
# exceeds its precision.

source(file.path("M2.3", "analysis", "r", "_theme.R"))

parity <- read_export("public_age_parity_points.csv")
psum   <- read_export("public_age_parity_summary.csv")
topo   <- read_export("topology_comparison.csv")

# ---- (a) published-age parity ----------------------------------------------
long <- bind_rows(
  parity |> transmute(reference = log10_reported_age, estimate = log10_est_age,
                      mode = "Held-out uncalibrated strict parity"),
  parity |> transmute(reference = log10_reported_age,
                      estimate = log10_calibrated_age,
                      mode = "Held-out calibrated emulation")
) |>
  mutate(mode = factor(mode, levels = c("Held-out uncalibrated strict parity",
                                        "Held-out calibrated emulation")))

lab <- psum |>
  mutate(mode = factor(mode, levels = levels(long$mode)),
         txt = sprintf("R2 = %.3f\nmedian |error| = %.3f\nwithin factor 2 = %.0f%%\nn = %d",
                       log10_r2, median_abs_log10_error,
                       100 * within_factor_2, n))

lim <- range(c(long$reference, long$estimate), na.rm = TRUE)

p_parity <- ggplot(long, aes(reference, estimate)) +
  geom_abline(slope = 1, intercept = 0, colour = INK, linetype = "dashed",
              linewidth = 0.4) +
  geom_abline(slope = 1, intercept = c(-log10(2), log10(2)), colour = INK_MUTED,
              linetype = "dotted", linewidth = 0.3) +
  geom_point(aes(colour = mode), size = 0.7, alpha = 0.35, stroke = 0) +
  geom_text(data = lab, aes(x = lim[1], y = lim[2], label = txt),
            hjust = 0, vjust = 1, size = 2.0, colour = INK,
            inherit.aes = FALSE) +
  scale_colour_manual(values = c(PAL[2], PAL[1]), guide = "none") +
  facet_wrap(~ mode, nrow = 1, labeller = labeller(mode = label_wrap_gen(28))) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = "Agreement with published lumped-parameter age outputs",
       subtitle = paste("Dotted lines bound a factor of two. The calibrated panel",
                        "is fitted to the reference and measures emulation, not",
                        "independent agreement."),
       x = expression("Published age, log"[10] * " years"),
       y = expression("Inferred age, log"[10] * " years")) +
  theme_hs() +
  theme(plot.subtitle = element_text(size = 6.5))

# ---- (b) directed topology --------------------------------------------------
tp <- topo |>
  mutate(mode_label = recode(mode,
           no_prior_head_gradient = "No prior\n(independent inference)",
           prior_assisted_ingestion = "Prior-assisted\n(ingestion fidelity only)")) |>
  select(mode_label, independent_validation,
         Precision = precision_recomputed, Recall = recall_recomputed,
         `F1` = f1_recomputed) |>
  pivot_longer(c(Precision, Recall, `F1`), names_to = "metric",
               values_to = "value") |>
  mutate(metric = factor(metric, levels = c("Precision", "Recall", "F1")))

p_topo <- ggplot(tp, aes(value, metric, fill = independent_validation)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.72), hjust = -0.18,
            size = 2.1, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = PAL[1], `FALSE` = "#B9C0C7"),
                    labels = c(`TRUE` = "Independent", `FALSE` = "Not independent"),
                    name = NULL) +
  facet_wrap(~ mode_label, nrow = 1) +
  scale_x_continuous(limits = c(0, 1.18), breaks = c(0, 0.5, 1),
                     expand = c(0, 0)) +
  labs(title = "Directed-topology comparison against a particle-tracking reference",
       subtitle = paste("The independent mode recovers most reference edges but",
                        "proposes more false than true ones"),
       x = "Score", y = NULL) +
  theme_hs() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 6.5))

fig <- p_parity / p_topo +
  plot_layout(heights = c(1.35, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "FIG-6_external_comparison", width_mm = 190, height_mm = 145)
