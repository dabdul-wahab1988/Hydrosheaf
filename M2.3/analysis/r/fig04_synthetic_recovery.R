# Figure 4: what the 100-realisation synthetic benchmark recovers, and what it
# does not. The panels are deliberately ordered from the best-supported target
# (transport parameters, network-constrained age) to the least-supported one
# (point reaction extents).

source(file.path("M2.3", "analysis", "r", "_theme.R"))

transport <- read_export("synthetic_transport_points.csv")
reaction  <- read_export("synthetic_reaction_points.csv")
age_pts   <- read_export("synthetic_age_points.csv")
recov     <- read_export("synthetic_recovery_summary.csv")
age_cmp   <- read_export("synthetic_age_method_comparison.csv")

# ---- (a) transport parameter recovery ---------------------------------------
t_stat <- recov |> filter(quantity == "Transport parameter")
lim_t <- range(c(transport$true_value, transport$recovered_value), na.rm = TRUE)

transport <- transport |>
  mutate(process = str_to_sentence(str_replace_all(process, "_", " ")))

p_trans <- ggplot(transport, aes(true_value, recovered_value)) +
  geom_abline(slope = 1, intercept = 0, colour = INK, linetype = "dashed",
              linewidth = 0.4) +
  geom_point(aes(colour = process), size = 0.85, alpha = 0.5, stroke = 0) +
  scale_colour_manual(values = PAL, name = NULL) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 1.6,
                                                             alpha = 1))) +
  coord_equal(xlim = lim_t, ylim = lim_t) +
  annotate("text", x = lim_t[1], y = lim_t[2], hjust = 0, vjust = 1,
           size = 2.0, colour = INK,
           label = sprintf("median |error| = %.3f\nmodel correct = %.0f%%\nn = %d of %d",
                           t_stat$median_absolute_error,
                           100 * t_stat$model_selection_accuracy,
                           t_stat$n_rows, t_stat$n_rows + t_stat$n_not_recovered)) +
  labs(title = "Transport parameters",
       subtitle = "Mixing and evaporation fractions",
       x = "True value", y = "Recovered value") +
  theme_hs() +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.8))

# ---- (b) age: single-node versus network-constrained -------------------------
age_long <- age_pts |>
  select(true_mrt_years, `Single-node lumped model` = single_node_lpm_years,
         `Network-constrained Bayesian` = network_bayesian_years) |>
  pivot_longer(-true_mrt_years, names_to = "method", values_to = "estimate") |>
  filter(is.finite(estimate), is.finite(true_mrt_years),
         estimate > 0, true_mrt_years > 0)

age_note <- age_cmp |>
  arrange(desc(method)) |>
  mutate(txt = sprintf("%s: R2 %.3f, MAE %.0f yr",
                       ifelse(startsWith(method, "Network"), "Network",
                              "Single-node"), log10_r2, mae_years)) |>
  pull(txt) |> paste(collapse = "\n")

# The single-node cloud is drawn first so the network layer is not hidden.
p_age <- ggplot(age_long, aes(true_mrt_years, estimate, colour = method)) +
  geom_abline(slope = 1, intercept = 0, colour = INK, linetype = "dashed",
              linewidth = 0.4) +
  geom_point(data = ~ filter(.x, method == "Single-node lumped model"),
             size = 0.7, alpha = 0.30, stroke = 0) +
  geom_point(data = ~ filter(.x, method == "Network-constrained Bayesian"),
             size = 0.7, alpha = 0.30, stroke = 0) +
  scale_x_log10(labels = label_number(big.mark = ",")) +
  scale_y_log10(labels = label_number(big.mark = ",")) +
  scale_colour_manual(values = c("Single-node lumped model" = PAL[2],
                                 "Network-constrained Bayesian" = PAL[1]),
                      name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6,
                                                             alpha = 1))) +
  annotate("text", x = max(age_long$true_mrt_years),
           y = min(age_long$estimate), hjust = 1, vjust = 0, size = 2.0,
           colour = INK, label = age_note) +
  coord_equal() +
  labs(title = "Residence time",
       subtitle = "Network vs single-node",
       x = "True mean residence time (years)",
       y = "Inferred age (years)") +
  theme_hs() +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.8))

# ---- (c) reaction extents: the non-identifiable target ----------------------
r_stat <- recov |> filter(startsWith(quantity, "Reaction"))
active <- reaction |> filter(term_status == "active")
inactive <- reaction |> filter(term_status == "inactive")
lim_r <- range(c(reaction$true_extent_mmolL, reaction$recovered_extent_mmolL),
               na.rm = TRUE)

p_react <- ggplot() +
  geom_abline(slope = 1, intercept = 0, colour = INK, linetype = "dashed",
              linewidth = 0.4) +
  geom_point(data = inactive,
             aes(x = true_extent_mmolL, y = recovered_extent_mmolL),
             colour = "#9AA0A6", size = 0.55, alpha = 0.20, stroke = 0) +
  geom_point(data = active,
             aes(x = true_extent_mmolL, y = recovered_extent_mmolL),
             colour = PAL[2], size = 0.75, alpha = 0.40, stroke = 0) +
  geom_hline(yintercept = c(-0.05, 0.05), colour = PAL[1], linewidth = 0.3,
             linetype = "dotted") +
  coord_equal(xlim = lim_r, ylim = lim_r) +
  annotate("text", x = lim_r[1], y = lim_r[2], hjust = 0, vjust = 1,
           size = 2.0, colour = INK,
           label = sprintf("active: R2 = %.3f, MAE = %.2f mmol/L\ninactive > 0.05 mmol/L: %.1f%%",
                           r_stat$r_squared, r_stat$mean_absolute_error,
                           100 * r_stat$inactive_term_leakage_fraction)) +
  labs(title = "Reaction extents",
       subtitle = "Grey marks truth-zero terms",
       x = "True extent (mmol/L)", y = "Recovered extent (mmol/L)") +
  theme_hs() +
  theme(plot.margin = margin(5, 6, 26, 5))

fig <- p_trans + p_age + p_react +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "FIG-4_synthetic_recovery", width_mm = 190, height_mm = 100)
