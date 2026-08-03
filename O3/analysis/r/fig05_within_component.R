# Figure 5: within-component identifiability detail, small multiples on each
# component's own retained-row taxonomy (Table 1), one consistent metric per
# panel so bar heights are directly comparable within a panel.

source(file.path("O3", "analysis", "r", "_theme.R"))

tax <- read_export("taxonomy.csv")

keep_metric <- c(
  "within-factor-2 agreement (reportable fits)",
  "F1 (edge classification)",
  "Phase F1 (mean, macro across scenarios)"
)

df <- tax |>
  filter(headline_metric %in% keep_metric, !is.na(headline_value)) |>
  mutate(
    scenario = str_replace_all(scenario, "_", " "),
    scenario = str_trunc(scenario, 26),
    row_type = factor(row_type, levels = c(
      "independent", "non_independent", "informed_control", "negative_control",
      "sensitivity_diagnostic", "prior_informed", "conventional_baseline"))
  )

row_type_pal <- c(
  independent = "#0072B2", non_independent = "#9AA0A6", informed_control = "#E69F00",
  negative_control = "#B0B0B0", sensitivity_diagnostic = "#CC79A7",
  prior_informed = "#9AA0A6", conventional_baseline = "#D55E00"
)

p <- ggplot(df, aes(x = reorder(scenario, headline_value), y = headline_value,
                     fill = row_type)) +
  geom_col(width = 0.7) +
  coord_flip() +
  facet_wrap(~component, scales = "free", ncol = 1) +
  scale_fill_manual(values = row_type_pal, name = "Row type",
                    guide = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Within-component detail: retained benchmark rows by type",
       x = NULL, y = "Component's own headline metric (Table 1)") +
  theme_o3(base_size = 8.3) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))

save_fig(p, "FIG-4", width_mm = 150, height_mm = 200)
