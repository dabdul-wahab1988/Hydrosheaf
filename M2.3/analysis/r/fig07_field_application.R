# Figure 7: the field application and its ceiling.
#
# Candidate edges are coordinate-neighbour constructions closed against the same
# chemistry used to fit them. The panels therefore report in-sample closure and
# the number of retained alternatives, never physical connectivity.

source(file.path("M2.3", "analysis", "r", "_theme.R"))

pts  <- read_export("field_edge_chemistry_points.csv")
fsum <- read_export("field_edge_chemistry_summary.csv")

pts <- pts |> filter(site != "unassigned")
site_levels <- c("Lower Anayari", "Talensi")
pts <- pts |> mutate(site = factor(site, levels = site_levels))

lab <- fsum |>
  filter(site %in% site_levels) |>
  mutate(site = factor(site, levels = site_levels),
         txt = sprintf("n = %d, median R2 = %.3f", n_candidate_edges,
                       median_chemistry_r2))
n_negative <- sum(pts$chemistry_r2 < 0)

# ---- (a) in-sample chemistry closure ---------------------------------------
p_r2 <- ggplot(pts, aes(chemistry_r2, site, colour = site)) +
  geom_jitter(height = 0.16, size = 0.85, alpha = 0.5, stroke = 0) +
  geom_boxplot(width = 0.30, outlier.shape = NA, fill = NA, colour = INK,
               linewidth = 0.32) +
  # A zero reference makes the negative closures readable as worse than the mean.
  geom_vline(xintercept = 0, colour = INK_MUTED, linewidth = 0.3,
             linetype = "dotted") +
  geom_text(data = lab, aes(x = -0.55, y = site, label = txt), hjust = 0,
            vjust = -1.5, size = 2.0, colour = INK, inherit.aes = FALSE) +
  scale_colour_manual(values = c("Lower Anayari" = PAL_DATASET[["Lower Anayari"]],
                                 "Talensi" = PAL_DATASET[["Talensi"]]),
                      guide = "none") +
  scale_x_continuous(limits = c(-0.6, 1.02),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = "In-sample chemistry closure",
       subtitle = sprintf("%d of %d site-assigned edges close worse than the mean",
                          n_negative, nrow(pts)),
       x = expression("Chemistry R"^2), y = NULL) +
  theme_hs() +
  theme(panel.grid.major.y = element_blank())

# ---- (b) retained reaction alternatives ------------------------------------
retained <- pts |> count(site, n_reactions_retained)

p_alt <- ggplot(retained, aes(n_reactions_retained, n, fill = site)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1),
           width = 0.8) +
  scale_fill_manual(values = c("Lower Anayari" = PAL_DATASET[["Lower Anayari"]],
                               "Talensi" = PAL_DATASET[["Talensi"]]),
                    name = NULL) +
  scale_x_continuous(breaks = scales::breaks_width(1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Retained reaction alternatives",
       subtitle = "Terms above 0.05 mmol/L per candidate edge",
       x = "Reaction terms retained", y = "Candidate edges") +
  theme_hs() +
  theme(legend.position = "bottom")

# ---- (c) measurement value --------------------------------------------------
mv <- read_export("field_measurement_value.csv") |>
  mutate(tier_label = factor(tier_label, levels = tier_label))

p_mv <- ggplot(mv, aes(pct_non_identifiable, tier_label)) +
  geom_col(fill = PAL[1], width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", pct_non_identifiable)),
            hjust = -0.18, size = 2.1, colour = INK) +
  scale_x_continuous(limits = c(0, 72), expand = c(0, 0)) +
  labs(title = "Measurement value",
       subtitle = sprintf(paste("Wells that remain non-identifiable as",
                                "determinands are added (n = %d)"),
                          mv$n_wells[1]),
       x = "Non-identifiable wells (%)", y = NULL) +
  theme_hs() +
  theme(panel.grid.major.y = element_blank())

fig <- p_r2 / p_alt / p_mv +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "FIG-7_field_application", width_mm = 140, height_mm = 168)
