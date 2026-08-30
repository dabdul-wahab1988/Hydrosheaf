# Figure 6: field- and archive-transfer scope across the three components.
# Two different score types (ELRI for reaction, edge F1 for topology) are
# faceted on free y-scales rather than forced onto one axis; age has no
# field-transfer row and is called out in the subtitle rather than silently
# dropped.

source(file.path("O3", "analysis", "r", "_theme.R"))

df <- read_export("field_transfer.csv") |>
  filter(!is.na(score_value))

p <- ggplot(df, aes(x = reorder(dataset, score_value), y = score_value,
                     fill = component)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("n=%d", n_edges)), hjust = -0.12, size = 2.4,
            colour = INK_MUTED) +
  coord_flip() +
  facet_wrap(~score_name, scales = "free", ncol = 1) +
  scale_fill_manual(values = PAL_COMPONENT, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(title = "Field- and archive-transfer scope is uneven across components",
       subtitle = "No field-transfer benchmark exists for the age layer; see Table 4",
       x = NULL, y = "Score (component-specific; see caption)") +
  theme_o3(base_size = 9) +
  theme(legend.position = "bottom")

save_fig(p, "FIG-6", width_mm = 150, height_mm = 130)
