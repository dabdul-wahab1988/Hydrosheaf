# Figure 5: the locked controlled-synthetic programme gates.
#
# Each panel plots the observed value against its predeclared threshold. A
# single common axis is used within each panel; no panel carries two scales.

source(file.path("M2.3", "analysis", "r", "_theme.R"))

metrics <- read_export("locked_gate_metrics.csv")
thresholds <- read_export("locked_gate_thresholds.csv")
contrasts <- read_export("prospective_policy_contrasts.csv")

# ---- (a) age and reaction: observed value against its acceptance floor -------
gate_spec <- tibble::tribble(
  ~component, ~metric,                          ~thr_key,               ~direction,
  "age",      "coverage_including_abstention",  "target_coverage",      "at least",
  "age",      "relative_width",                 "max_relative_width",   "at most",
  "age",      "acceptance_rate",                "minimum_acceptance_rate", "at least",
  "reaction", "coverage",                       "minimum_coverage",     "at least",
  "reaction", "max_classwise_ece",              "max_classwise_ece",    "at most",
  "reaction", "false_commitment_rate",          "max_false_commitment", "at most",
  "reaction", "selective_risk",                 "max_selective_risk",   "at most",
  "kinetics", "predictive_rmse",                "max_predictive_rmse",  "at most",
  "kinetics", "interval_coverage",              "minimum_interval_coverage", "at least",
  "kinetics", "identifiability_rate",           "minimum_identifiability_rate", "at least"
)

gates <- gate_spec |>
  left_join(metrics |> select(component, metric, label, value),
            by = c("component", "metric")) |>
  left_join(thresholds |> rename(thr_key = requirement),
            by = c("component", "thr_key")) |>
  mutate(component = recode(component, age = "Age", reaction = "Reaction family",
                            kinetics = "Kinetics"),
         component = factor(component, levels = c("Age", "Reaction family",
                                                  "Kinetics")),
         label = factor(label, levels = rev(unique(label))),
         passes = ifelse(direction == "at least", value >= threshold,
                         value <= threshold))

p_gates <- ggplot(gates, aes(y = label)) +
  geom_segment(aes(x = threshold, xend = value, yend = label),
               colour = "#C8CDD2", linewidth = 1.4,
               lineend = "round") +
  geom_point(aes(x = threshold), colour = INK_MUTED, shape = 124, size = 2.4) +
  geom_point(aes(x = value, colour = passes), size = 2.0) +
  geom_text(aes(x = value, label = sprintf("%.3g", value)),
            vjust = -1.1, size = 2.0, colour = INK) +
  scale_colour_manual(values = c(`TRUE` = PAL[1], `FALSE` = PAL[2]),
                      labels = c(`TRUE` = "Meets threshold",
                                 `FALSE` = "Below threshold"), name = NULL) +
  facet_grid(component ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.14))) +
  labs(title = "Component gates vs thresholds",
       subtitle = "Tick = threshold; dot = observed",
       x = "Metric value", y = NULL) +
  theme_hs() +
  theme(legend.position = "bottom", strip.placement = "outside",
        panel.grid.major.y = element_blank())

# ---- (b) prospective decision utility ---------------------------------------
util <- metrics |>
  filter(component == "integrated",
         metric %in% c("hydrosheaf_mean_utility_per_cost",
                       "random_mean_utility_per_cost",
                       "strongest_specialist_mean_utility_per_cost")) |>
  mutate(policy = recode(metric,
                         hydrosheaf_mean_utility_per_cost = "HydroSheaf",
                         random_mean_utility_per_cost = "Random",
                         strongest_specialist_mean_utility_per_cost = "Strongest specialist"),
         policy = stats::reorder(policy, value))

p_util <- ggplot(util, aes(value, policy, fill = policy == "HydroSheaf")) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, colour = INK_MUTED, linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", value),
                hjust = ifelse(value < 0, 1.18, -0.18)),
            size = 2.2, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = PAL[1], `FALSE` = "#B9C0C7"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.22))) +
  labs(title = "Prospective decision utility",
       subtitle = "Six locked truth-blind cases",
       x = "Mean utility per cost", y = NULL) +
  theme_hs() +
  theme(panel.grid.major.y = element_blank())

# ---- (c) paired contrasts ---------------------------------------------------
# Only the contrasts in which HydroSheaf is the left policy are reported; the
# random-versus-specialist contrast is not a claim of this study.
ctr <- contrasts |>
  filter(left_policy == "hydrosheaf") |>
  mutate(label = recode(contrast,
                        hydrosheaf_vs_random = "vs random",
                        hydrosheaf_vs_specialist = "vs strongest specialist"))

p_ctr <- ggplot(ctr, aes(mean_delta, label)) +
  geom_vline(xintercept = 0, colour = INK, linetype = "dashed",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.10,
                 colour = PAL[1], linewidth = 0.5) +
  geom_point(size = 2.0, colour = PAL[1]) +
  geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]", mean_delta, ci_low, ci_high)),
            vjust = -1.2, size = 2.0, colour = INK) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(title = "Paired utility differences",
       subtitle = paste("95% within-run descriptive intervals,",
                        "not population bounds"),
       x = "Difference in utility per cost", y = NULL) +
  theme_hs() +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 6.5))

fig <- p_gates | (p_util / p_ctr)
fig <- fig + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "FIG-5_programme_gates", width_mm = 190, height_mm = 118)
