# Figure 2: the measured field datasets - location, analytical quality,
# hydrochemical composition and stable-isotope character.
#
# Only primary measured panels are plotted. The reconstructed Northern Ghana
# seasonal panel is excluded throughout (DECISION D2).

source(file.path("M2.3", "analysis", "r", "_theme.R"))
source(file.path("M2.3", "analysis", "r", "_map.R"))

samples <- read_export("field_samples_derived.csv") |>
  filter(primary_panel)
mwl     <- read_export("field_meteoric_water_lines.csv")
facies  <- read_export("field_facies_counts.csv")

ds_levels <- c("Northern Ghana", "Lower Anayari", "Talensi")
samples <- samples |> mutate(dataset = factor(dataset, levels = ds_levels))

# ---- (a) study-area map -----------------------------------------------------
p_map <- build_study_map(samples, PAL_DATASET, ds_levels)

# ---- (b) charge-balance error ----------------------------------------------
cbe <- samples |> filter(is.finite(cbe_percent))
tier_counts <- cbe |> count(dataset, quality_tier)

p_cbe <- ggplot(cbe, aes(x = dataset, y = cbe_percent)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -5, ymax = 5,
           fill = "#0072B2", alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5, ymax = 10,
           fill = "#E69F00", alpha = 0.10) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -10, ymax = -5,
           fill = "#E69F00", alpha = 0.10) +
  geom_hline(yintercept = 0, colour = INK_MUTED, linewidth = 0.25) +
  geom_jitter(aes(colour = quality_tier), width = 0.16, size = 0.75,
              alpha = 0.7, stroke = 0) +
  geom_boxplot(width = 0.32, outlier.shape = NA, fill = NA,
               colour = INK, linewidth = 0.32) +
  scale_colour_manual(values = PAL_TIER, name = NULL,
                      limits = c("quantitative", "screening", "exploratory")) +
  scale_x_discrete(labels = function(x) str_wrap(x, 9)) +
  labs(tag = "(b)", title = "Charge-balance error",
       subtitle = "Bands mark the +/-5% and +/-10% acceptance tiers",
       x = NULL, y = "Charge-balance error (%)") +
  theme_hs() +
  theme(legend.position = "bottom")

# ---- (c) hydrochemical facies ----------------------------------------------
fac <- facies |>
  mutate(dataset = factor(dataset, levels = ds_levels),
         facies = stats::reorder(facies, pct, FUN = sum))

p_fac <- ggplot(fac, aes(x = pct, y = facies, fill = dataset)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.12),
           width = 0.78) +
  scale_fill_manual(values = PAL_DATASET, name = NULL) +
  scale_x_continuous(labels = label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(tag = "(c)", title = "Dominant-ion facies",
       subtitle = "Share of samples in each dataset",
       x = "Samples", y = NULL) +
  theme_hs() +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank())

# ---- (d) stable isotopes ----------------------------------------------------
iso <- samples |> filter(is.finite(d18O), is.finite(d2H))
rng <- range(iso$d18O)
gmwl <- tibble(d18O = rng, d2H = 8 * rng + 10)
# Each local line is drawn only across the range its own samples span, so no
# regression is extrapolated into a region where that dataset has no data.
iso_rng <- iso |> group_by(dataset) |>
  summarise(lo = min(d18O), hi = max(d18O), .groups = "drop")
lmwl <- mwl |>
  mutate(dataset = factor(dataset, levels = ds_levels)) |>
  left_join(iso_rng, by = "dataset") |>
  group_by(dataset) |>
  group_modify(~ tibble(d18O = c(.x$lo, .x$hi),
                        d2H = .x$intercept + .x$slope * c(.x$lo, .x$hi))) |>
  ungroup()

iso_lab <- mwl |>
  mutate(label = sprintf("%s: d2H = %.2f d18O %+.1f (R2 = %.2f)",
                         dataset, slope, intercept, r_squared)) |>
  pull(label) |> paste(collapse = "\n")

p_iso <- ggplot(iso, aes(d18O, d2H)) +
  geom_line(data = gmwl, aes(d18O, d2H), colour = INK, linewidth = 0.45,
            linetype = "dashed", inherit.aes = FALSE) +
  geom_point(aes(colour = dataset), size = 1.0, alpha = 0.75, stroke = 0) +
  geom_line(data = lmwl, aes(d18O, d2H, colour = dataset), linewidth = 0.55,
            inherit.aes = FALSE) +
  annotate("text", x = max(rng), y = 8 * max(rng) + 10, label = "GMWL",
           hjust = 1.1, vjust = 1.4, size = 2.4, colour = INK) +
  scale_colour_manual(values = PAL_DATASET, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  labs(tag = "(d)", title = "Stable isotope composition",
       subtitle = "Local regressions vs the global line (GMWL)",
       x = expression(delta^18 * "O (per mil VSMOW)"),
       y = expression(delta^2 * "H (per mil VSMOW)"),
       caption = iso_lab) +
  theme_hs() +
  theme(legend.position = "none")

fig <- (p_map | p_cbe) / (p_fac | p_iso) &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "FIG-2_field_setting", width_mm = 190, height_mm = 165)
