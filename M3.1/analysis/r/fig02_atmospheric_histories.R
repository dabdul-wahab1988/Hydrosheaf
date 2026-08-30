# Figure 2: atmospheric input histories for tritium, SF6 and CFC-12, showing
# why young tracers must be interpreted jointly (complementary and distinct
# ambiguity structure). Carried over from M3 Figure 1, replotted in R from a
# regenerated export.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

hist_data <- read_export("atmospheric_input_histories.csv")

panel_labels <- c(
  "3H"    = "(a) Tritium (3H)",
  "SF6"   = "(b) Sulfur hexafluoride (SF6)",
  "CFC12" = "(c) Chlorofluorocarbon-12 (CFC-12)"
)
hist_data <- hist_data |>
  mutate(panel = factor(panel_labels[tracer], levels = unname(panel_labels)))

p <- ggplot(hist_data, aes(x = recharge_year, y = value)) +
  geom_line(colour = PAL[1], linewidth = 0.45) +
  facet_wrap(~panel, scales = "free_y", ncol = 1) +
  labs(
    title = "Atmospheric input histories of transient young-water tracers",
    subtitle = "Complementary timing and shape give the tracer suite joint identifying power",
    x = "Recharge year", y = "Input concentration (tracer-specific unit)"
  ) +
  theme_hs()

save_fig(p, "FIG-2_atmospheric_histories", width_mm = 165, height_mm = 175)
