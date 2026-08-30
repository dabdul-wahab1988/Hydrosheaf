# Figure 3: the central cross-component result. Capture-type versus
# correctness-type metrics, independent/uncalibrated evaluation, all three
# components on their own native scale (paired bars, not one merged axis;
# see O3/analysis/python/derive_headline_metrics.py docstring for why).

source(file.path("O3", "analysis", "r", "_theme.R"))

df <- read_export("headline_metrics.csv") |>
  filter(axis %in% c("capture", "correctness")) |>
  mutate(
    axis = factor(recode(axis, capture = "Capture\n(fraction detected)",
                          correctness = "Correctness\n(fraction correct)"),
                  levels = c("Capture\n(fraction detected)",
                             "Correctness\n(fraction correct)")),
    metric_label = recode(metric,
      recall = "Recall (edges)", precision = "Precision (edges)",
      reportability_rate = "Reportability rate", within_factor_2 = "Within factor of 2",
      recall_pooled = "Recall (phases, pooled)", precision_macro = "Precision (phases, macro)")
  )

p <- ggplot(df, aes(x = axis, y = value, fill = component)) +
  geom_col(width = 0.62, position = position_dodge2(preserve = "single")) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.35, colour = INK,
                position = position_dodge2(width = 0.62, padding = 0.55)) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -1.3, size = 2.5,
            colour = INK, position = position_dodge2(width = 0.62, preserve = "single")) +
  geom_text(aes(label = metric_label, y = -0.05), size = 1.9, colour = INK_MUTED,
            angle = 90, hjust = 1, position = position_dodge2(width = 0.62, preserve = "single")) +
  facet_wrap(~component, nrow = 1) +
  scale_fill_manual(values = PAL_COMPONENT, guide = "none") +
  scale_y_continuous(limits = c(-0.55, 1.05), breaks = seq(0, 1, 0.25),
                      expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Recall/sensitivity exceeds precision/specificity in all three layers",
       subtitle = "Independent, prior-free, uncalibrated evaluation; error bars are 95% Wilson score intervals (Methods, Table 2)",
       x = NULL, y = "Value (0-1 scale)") +
  theme_o3(base_size = 9) +
  theme(strip.text = element_text(size = 8), panel.spacing = unit(4, "mm"))

save_fig(p, "FIG-2", width_mm = 190, height_mm = 110)
