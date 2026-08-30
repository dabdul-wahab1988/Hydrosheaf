# Supplementary Figure S2: spatial-distance and vertical-difference
# distributions for graph edges in the N = 329 strict reportable-node
# benchmark. Diagnoses hypothesised edge geometry, not verified flow
# connections.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

family_labels <- c(
  coordinate_global               = "Coordinate (global)",
  study_unit_coordinate            = "Coordinate (study unit)",
  depth_constrained                = "Depth constrained",
  hydraulic_proxy_constrained      = "Hydraulic proxy",
  parameter_smooth_aicc            = "Parameter smoothing",
  wrong_direction_negative_control = "Wrong-direction control",
  randomized_negative_control      = "Randomised control"
)

edges <- read_export("edge_geometry.csv") |>
  filter(family %in% names(family_labels)) |>
  mutate(family_label = factor(family_labels[family], levels = unname(family_labels)))

p_dist <- ggplot(edges, aes(x = family_label, y = distance_km)) +
  geom_boxplot(width = 0.5, outlier.size = 0.6, outlier.alpha = 0.4,
               colour = INK, fill = PAL[1], alpha = 0.25, linewidth = 0.3) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                     labels = scales::comma,
                     breaks = c(0, 1, 10, 100, 1000)) +
  coord_flip() +
  labs(tag = "(a)", title = "Edge distance", x = NULL, y = "Distance (km, log scale)") +
  theme_hs()

p_depth <- ggplot(edges, aes(x = family_label, y = depth_diff_m)) +
  geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.3) +
  geom_boxplot(width = 0.5, outlier.size = 0.6, outlier.alpha = 0.4,
               colour = INK, fill = PAL[2], alpha = 0.25, linewidth = 0.3) +
  coord_flip() +
  labs(tag = "(b)", title = "Downstream minus upstream depth", x = NULL, y = "Depth difference (m)") +
  theme_hs() +
  theme(axis.text.y = element_blank())

p <- p_dist + p_depth +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    title = "Graph-edge geometry by family, N = 329 reportable nodes",
    subtitle = "Hypothesised connectivity, not verified flow paths.",
    theme = theme(plot.title = element_text(size = 9.5, face = "bold"),
                  plot.subtitle = element_text(size = 7.5, colour = INK_MUTED))
  )

save_fig(p, "FIG-S2_edge_geometry", width_mm = 175, height_mm = 90)
