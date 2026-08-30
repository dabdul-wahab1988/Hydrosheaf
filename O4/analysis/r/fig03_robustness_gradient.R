# FIG-3: M6 tier-ablation cliff -- mean MRS (flat, internal) vs percent
# non-identifiable (collapsing, external), five tiers, with the
# evidence-gate-off negative control.
source("O4/analysis/r/_theme.R")

rg <- read_export("robustness_gradient.csv") %>%
  filter(condition == "field_dataset_tier_ablation") %>%
  arrange(tier_order) %>%
  mutate(axis_x = factor(axis_x, levels = unique(axis_x)))

p_top <- ggplot(rg, aes(x = axis_x, y = internal_signal_value, group = 1)) +
  geom_line(colour = PAL_SIGNAL["internal"], linewidth = 0.6) +
  geom_point(colour = PAL_SIGNAL["internal"], size = 2) +
  geom_text(aes(label = sprintf("%.1f", internal_signal_value)), vjust = -1, size = 2.4) +
  coord_cartesian(ylim = c(60, 75)) +
  labs(y = "Mean MRS\n(internal, fit-quality)", x = NULL,
       title = "Internal fit quality stays flat while identifiability collapses",
       subtitle = "M6 Northern Ghana, 160 wells, 5-level evidence-tier ablation") +
  theme_o4()

p_bot <- ggplot(rg, aes(x = axis_x, y = external_signal_value, group = 1)) +
  geom_line(colour = PAL_SIGNAL["external"], linewidth = 0.6) +
  geom_point(colour = PAL_SIGNAL["external"], size = 2) +
  geom_text(aes(label = sprintf("%.1f%%", external_signal_value)), vjust = -1, size = 2.4) +
  coord_cartesian(ylim = c(-2, 65)) +
  labs(y = "% wells non-identifiable\n(external, verified)", x = "Evidence tier removed (right to left)",
       caption = "Evidence-gate-off negative control returns 0% non-identifiable at every tier (not plotted; see Table 3), confirming the collapse is the framework's conservative prior, not a classifier artefact.") +
  theme_o4()

p <- (p_top / p_bot)
save_fig(p, "FIG-3", width_mm = 140, height_mm = 130)
