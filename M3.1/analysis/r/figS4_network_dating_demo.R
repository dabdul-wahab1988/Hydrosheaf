# Supplementary Figure S4: controlled-synthetic network-dating ambiguity and
# recovery diagnostic. Uses known simulated ages to test implementation
# behaviour; it is not field validation and does not justify a general
# graph-improvement claim. Adapted from
# M3/m3_age_benchmark/scripts/plot_m3_network_dating_demo.R onto the M3.1
# theme (self-contained; no longer depends on the M6 theme file).

source(file.path("M3.1", "analysis", "r", "_theme.R"))

demo <- read_export("m3_network_dating_demo.csv")
summ <- read_export("m3_network_dating_demo_summary.csv")

al <- demo |>
  filter(ambiguous == 1) |>
  select(true_age, age_single, age_network) |>
  pivot_longer(-true_age, names_to = "method", values_to = "est") |>
  mutate(method = recode(method, age_single = "Single-node",
                         age_network = "Network (flow-ordered)"))

p_a <- ggplot(al, aes(true_age, est, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = INK_MUTED) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_colour_manual(values = c("Single-node" = PAL[2],
                                "Network (flow-ordered)" = PAL[1]), name = NULL) +
  labs(tag = "(a)", subtitle = "Synthetic ambiguity: age recovery",
       x = "True age (yr)", y = "Estimated age (yr)") +
  theme_hs()

bl <- summ |>
  select(subset, Single = wf2_single, Network = wf2_network) |>
  pivot_longer(-subset, names_to = "method", values_to = "wf2") |>
  mutate(subset = factor(subset, levels = c("unambiguous", "ambiguous", "all"),
                         labels = c("Unambiguous", "Ambiguous", "All")))

p_b <- ggplot(bl, aes(subset, wf2, fill = method)) +
  geom_col(position = position_dodge(0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * wf2)),
            position = position_dodge(0.72), vjust = -0.4, size = 2.4) +
  scale_fill_manual(values = c("Single" = PAL[2], "Network" = PAL[1]), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08),
                     expand = expansion(c(0, 0.02))) +
  labs(tag = "(b)", subtitle = "Synthetic within-factor-2 agreement",
       x = NULL, y = "Within factor 2 of true age") +
  theme_hs()

p <- (p_a | p_b) +
  plot_annotation(
    caption = str_wrap("Controlled-synthetic capability demonstration: flow-ordering resolves selected tritium bomb-peak aliases. This is not field validation and does not establish general graph benefit.", 110),
    theme = theme(plot.caption = element_text(size = 7, colour = INK_MUTED, hjust = 0))
  )

save_fig(p, "FIG-S4_network_dating_demo", width_mm = 175, height_mm = 92)
