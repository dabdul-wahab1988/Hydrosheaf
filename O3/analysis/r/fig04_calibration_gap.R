# Figure 4: screening/independent versus calibrated, emulated, or
# prior-informed evaluation. The two calibration exercises (age, reaction)
# are not equally rigorous; panel subtitles carry the caveat rather than
# letting the numbers speak alone (Outline.md central argument, point 4).

source(file.path("O3", "analysis", "r", "_theme.R"))

df <- read_export("calibration_gap.csv")

age <- df |> filter(component == "Age / residence time", metric == "log10_r2") |>
  mutate(condition = recode(condition, independent = "Independent\n(held-out, uncalibrated)",
                            calibrated_emulation = "Calibrated\n(emulation, non-independent)"))
p1 <- ggplot(age, aes(x = condition, y = value, fill = condition)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = c("Independent\n(held-out, uncalibrated)" = "#0072B2",
                               "Calibrated\n(emulation, non-independent)" = "#9AA0A6"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Age: log10 R2", subtitle = "Calibration fit on the same\nheld-out folds it is scored against",
       x = NULL, y = "log10 R2") +
  theme_o3(base_size = 9)

rx <- df |> filter(component == "Reaction") |>
  mutate(condition = recode(condition,
    held_out_archetype_transfer = "Held-out archetype\n(genuine transfer)",
    chance_level_four_class = "Chance level\n(4 classes)"))
p2 <- ggplot(rx, aes(x = condition, y = value, fill = condition)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.4, size = 2.6) +
  scale_fill_manual(values = c("Held-out archetype\n(genuine transfer)" = "#009E73",
                               "Chance level\n(4 classes)" = "#9AA0A6"), guide = "none") +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Reaction: MRS 4-class accuracy",
       subtitle = "Trained on 3 archetypes, tested on\na 4th, untouched archetype",
       x = NULL, y = "Accuracy") +
  theme_o3(base_size = 9)

topo <- df |> filter(component == "Topology")
topo <- topo |>
  mutate(label = recode(condition, independent = "Independent",
                        prior_informed_override = "Prior: override",
                        prior_informed_merge = "Prior: merge",
                        prior_informed_only = "Prior: only"),
         label = factor(label, levels = c("Independent", "Prior: override",
                                          "Prior: merge", "Prior: only")),
         panel = ifelse(metric == "f1", "F1 (independent only)", "Output-graph edges"))
p3 <- ggplot(topo, aes(x = label, y = value, fill = is_independent)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(metric == "f1", sprintf("%.3f", value),
                               scales::comma(value))), vjust = -0.4, size = 2.4) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#9AA0A6"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Topology: F1 vs prior-informed\nedge count",
       subtitle = "Different metrics; separate axes",
       x = NULL, y = "Value") +
  theme_o3(base_size = 8.5) +
  theme(axis.text.x = element_text(size = 6, angle = 25, hjust = 1),
        strip.text = element_text(size = 6.3))

if (has_patchwork) {
  p <- p1 + p2 + p3 + plot_layout(nrow = 1)
} else {
  p <- p1
}

save_fig(p, "FIG-3", width_mm = 190, height_mm = 105)
