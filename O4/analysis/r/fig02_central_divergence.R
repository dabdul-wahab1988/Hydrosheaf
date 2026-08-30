# FIG-2: central result. Paired internal confidence signal vs external truth
# signal, one panel per component, on each component's own native scale, in
# increasing order of stress.
source("O4/analysis/r/_theme.R")

rg <- read_export("robustness_gradient.csv") %>% filter(condition == "field_dataset_tier_ablation") %>% arrange(tier_order)
rg_long <- bind_rows(
  rg %>% transmute(axis_x, signal = "internal", value = internal_signal_value, order = tier_order),
  rg %>% transmute(axis_x, signal = "external", value = external_signal_value, order = tier_order)
) %>% mutate(axis_x = factor(axis_x, levels = unique(axis_x[order(order)])))

p_m6 <- ggplot(rg_long, aes(x = axis_x, y = value, colour = signal, group = signal)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
  scale_colour_manual(values = c(internal = unname(PAL_SIGNAL["internal"]), external = unname(PAL_SIGNAL["external"])),
                       labels = c(internal = "internal: mean MRS", external = "external: % non-identifiable")) +
  labs(title = "M6 robustness", x = NULL, y = "Value (native scale)", colour = NULL) +
  theme_o4() + theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1))

iv <- read_export("integration_value.csv") %>%
  mutate(axis_x = factor(axis_x, levels = axis_x))
iv_long <- bind_rows(
  iv %>% transmute(axis_x, signal = "internal", value = internal_signal_value),
  iv %>% transmute(axis_x, signal = "external", value = external_signal_value)
)
p_m7 <- ggplot(iv_long, aes(x = axis_x, y = value, fill = signal)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  scale_fill_manual(values = c(internal = unname(PAL_SIGNAL["internal"]), external = unname(PAL_SIGNAL["external"])),
                     labels = c(internal = "internal: entropy change", external = "external: PR-AUC change")) +
  labs(title = "M7 identifiability", x = NULL, y = "Change (native scale)", fill = NULL) +
  theme_o4() + theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1))

me <- read_export("calibration_matched_extremes.csv")
me_long <- bind_rows(
  me %>% transmute(axis_x, signal = "internal", value = internal_signal_2_value),
  me %>% transmute(axis_x, signal = "external", value = external_signal_value)
) %>% mutate(axis_x = factor(axis_x, levels = unique(axis_x)))
p_m8 <- ggplot(me_long, aes(x = axis_x, y = value, fill = signal)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(internal = unname(PAL_SIGNAL["internal"]), external = unname(PAL_SIGNAL["external"])),
                     labels = c(internal = "internal: 95% coverage", external = "external: median |log10 error|")) +
  labs(title = "M8 calibration", x = NULL, y = "Value (native scale)", fill = NULL,
       caption = "Success rate is 1.0 in every M8 design (not shown); coverage stays within a narrow band while recovery error spans more than an order of magnitude at constant success.") +
  theme_o4() + theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1, size = 6))

p <- (p_m6 | p_m7 | p_m8) +
  plot_annotation(title = "Internal confidence signals and external truth signals diverge under stress in all three layers",
                   subtitle = "Each panel uses its own native metric and scale; panels are not numerically comparable to one another") &
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8, colour = INK_MUTED))

save_fig(p, "FIG-2", width_mm = 190, height_mm = 105)
