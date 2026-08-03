# Figure 5: graph-age effects on the N = 329 strict reportable nodes. The
# decision rule needs both axes jointly: a robust improvement requires lower
# RMSE (left of the vertical zero line) AND unchanged-or-better within-
# factor-2 agreement (at or above the horizontal zero line). Only the
# hydraulic-proxy family sits near the origin; every other family, and both
# negative controls most severely, falls in the "worse on both" quadrant.
# Carried over from M3 Figure 4, redesigned as a single joint-decision plot
# because a small-multiples RMSE-only view hid the numerically tiny
# hydraulic-proxy effect that is the paper's central near-miss.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

family_labels <- c(
  coordinate_global               = "Coordinate (global)",
  study_unit_coordinate            = "Coordinate (study unit)",
  depth_constrained                = "Depth constrained",
  hydraulic_proxy_constrained      = "Hydraulic proxy",
  parameter_smooth_aicc            = "Parameter smoothing",
  wrong_direction_negative_control = "Wrong-direction (negative control)",
  randomized_negative_control      = "Randomised (negative control)"
)
family_type <- c(
  coordinate_global = "Candidate", study_unit_coordinate = "Candidate",
  depth_constrained = "Candidate", hydraulic_proxy_constrained = "Candidate",
  parameter_smooth_aicc = "Candidate",
  wrong_direction_negative_control = "Negative control",
  randomized_negative_control = "Negative control"
)
strength_levels <- c("none", "weak", "medium", "strong")

graph <- read_export("Manuscript_Table4_Real_USGS_Graph_Benchmark.csv") |>
  filter(graph_family %in% names(family_labels), prior_strength != "none") |>
  mutate(
    family_label = factor(family_labels[graph_family], levels = unname(family_labels)),
    type = family_type[graph_family],
    delta_within_factor_2 = within_factor_2_graph - within_factor_2_single,
    prior_strength = factor(prior_strength, levels = strength_levels)
  )

p <- ggplot(graph, aes(x = delta_rmse_graph_minus_single, y = delta_within_factor_2,
                       colour = family_label, shape = type)) +
  geom_vline(xintercept = 0, colour = "#009E73", linewidth = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "#009E73", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.9, alpha = 0.85) +
  scale_colour_manual(
    values = setNames(c(PAL[1:5], "#8A9199", "#5F6368"), unname(family_labels)),
    name = NULL, guide = guide_legend(ncol = 2)
  ) +
  scale_shape_manual(values = c(Candidate = 16, `Negative control` = 4), name = NULL) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Graph-age effects on the N = 329 strict reportable nodes",
    subtitle = str_wrap("Each point is one family at one non-zero prior strength. A robust improvement needs lower RMSE and unchanged-or-better within-factor-2 agreement jointly (left of and above the dashed lines); only the hydraulic-proxy family approaches the origin, and even there within-factor-2 agreement still falls.", 105),
    x = "Delta log10 RMSE (graph minus single-node; negative is better)",
    y = "Delta within-factor-2 agreement (graph minus single-node; positive is better)"
  ) +
  theme_hs() +
  theme(legend.position = "right", legend.text = element_text(size = 6.4))

save_fig(p, "FIG-5_graph_benchmark", width_mm = 190, height_mm = 118)
