# Figure 4: strict reported-configuration emulation against model-derived
# USGS reference ages, N = 329 reportable subset. Evaluates reference-workflow
# agreement, not independently known groundwater ages (DECISION: never
# describe this axis as true age). Carried over from M3 Figure 3.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

scatter <- read_export("strict_parity_scatter.csv")

age_breaks <- c(1, 10, 100, 1000, 10000, 100000)
age_class_labels <- c(
  "modern_le_50"       = "Modern (≤50 yr)",
  "intermediate_50_1000" = "Intermediate (50–1,000 yr)",
  "old_1000_30000"     = "Old (1,000–30,000 yr)",
  "very_old_gt_30000"  = "Very old (>30,000 yr)"
)
scatter <- scatter |>
  mutate(age_class_label = factor(age_class_labels[age_class],
                                  levels = unname(age_class_labels)))

p <- ggplot(scatter, aes(x = ref_age, y = est_age_multi)) +
  geom_abline(slope = 1, intercept = 0, colour = INK_MUTED, linewidth = 0.3) +
  geom_abline(slope = 1, intercept = log10(2), colour = INK_MUTED,
              linewidth = 0.25, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -log10(2), colour = INK_MUTED,
              linewidth = 0.25, linetype = "dashed") +
  geom_point(aes(colour = age_class_label), size = 1.1, alpha = 0.75, stroke = 0) +
  scale_x_log10(breaks = age_breaks, labels = scales::comma) +
  scale_y_log10(breaks = age_breaks, labels = scales::comma) +
  scale_colour_manual(values = PAL, name = "Age class") +
  coord_fixed() +
  labs(
    title = "Strict reported-configuration emulation vs. USGS reference ages",
    subtitle = str_wrap("N = 329 reportable fits; dashed lines mark a factor of two. Reference ages are model-derived, not observed truth.", 78),
    x = "USGS reference age (years, log scale)",
    y = "Hydrosheaf emulated age (years, log scale)"
  ) +
  theme_hs() +
  theme(legend.position = "bottom")

save_fig(p, "FIG-4_strict_parity_scatter", width_mm = 165, height_mm = 165)
