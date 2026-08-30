# Figure 2: dataset and benchmark-scale overview.

source(file.path("O3", "analysis", "r", "_theme.R"))

scale_df <- read_export("benchmark_scale.csv") |>
  select(component, n_scenarios, n_fits) |>
  pivot_longer(c(n_scenarios, n_fits), names_to = "quantity", values_to = "n") |>
  mutate(quantity = recode(quantity, n_scenarios = "Scenarios", n_fits = "Fits"),
         quantity = factor(quantity, levels = c("Scenarios", "Fits")))

p1 <- ggplot(scale_df, aes(x = component, y = n, fill = component)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.35, size = 2.4,
            colour = INK) +
  facet_wrap(~quantity, scales = "free_y") +
  scale_fill_manual(values = PAL_COMPONENT, guide = "none") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Benchmark scale differs by roughly two orders of magnitude",
       x = NULL, y = "Count (log-free axis per panel)") +
  theme_o3(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

runtime_df <- read_export("benchmark_scale.csv") |>
  filter(runtime_recorded) |>
  select(component, median_runtime_per_fit_ms)

p2 <- ggplot(runtime_df, aes(x = component, y = median_runtime_per_fit_ms,
                              fill = component)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.1f ms", median_runtime_per_fit_ms)),
            vjust = -0.4, size = 2.4, colour = INK) +
  scale_fill_manual(values = PAL_COMPONENT, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Only one component records per-fit runtime",
       x = NULL, y = "Median runtime per fit (ms)") +
  theme_o3(base_size = 9)

if (has_patchwork) {
  p <- p1 / p2 + plot_layout(heights = c(1.6, 1))
} else {
  p <- p1
}

save_fig(p, "FIG-5", width_mm = 190, height_mm = 150)
