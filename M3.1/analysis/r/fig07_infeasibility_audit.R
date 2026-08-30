# Figure 7 (new in M3.1): independent set-valued compatibility audit and the
# tracer-level local-infeasibility finding. Exploratory/development-stage
# diagnostic (configs/identified_ttd_protocol.yaml: status "development",
# claim_authority "implementation_only") -- shown as a secondary, honestly
# bounded result, not as confirmatory validation. Seven candidate
# explanations were tested for the infeasibility signature and none was
# supported (docs/m3_cfc_reconciliation_step1_20260731.md); no cause is
# attributed here.

source(file.path("M3.1", "analysis", "r", "_theme.R"))

pairwise <- read_export("infeasibility_pairwise.csv")
conditioning <- read_export("infeasibility_by_conditioning_size.csv")

tracer_order <- c("14C", "3H", "3H/3He", "SF6", "CFC11", "CFC113", "CFC12")
pairwise <- pairwise |>
  mutate(
    tracer_a = factor(tracer_a, levels = tracer_order),
    tracer_b = factor(tracer_b, levels = tracer_order)
  )
# Symmetrise for a full square heatmap (the source reports each unordered
# pair once).
pairwise_sym <- bind_rows(
  pairwise,
  pairwise |> rename(tracer_a = tracer_b, tracer_b = tracer_a)
)

p_heat <- ggplot(pairwise_sym, aes(x = tracer_a, y = tracer_b, fill = rate)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = scales::percent(rate, accuracy = 1)), size = 2.1, colour = INK) +
  scale_fill_gradient(low = SEQ_LOW, high = SEQ_HIGH, labels = scales::percent,
                      name = "Pairwise\ninfeasibility rate", na.value = "white") +
  coord_fixed() +
  labs(tag = "(a)", title = "Pairwise tracer infeasibility",
       x = NULL, y = NULL) +
  theme_hs() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        axis.text.y = element_text(size = 6.5),
        panel.grid.major = element_blank())

p_cond <- ggplot(conditioning, aes(x = factor(n_conditioning_tracers), y = rate, group = 1)) +
  geom_line(colour = PAL[1], linewidth = 0.4) +
  geom_point(colour = PAL[1], size = 1.6) +
  geom_text(aes(label = sprintf("%d/%d", infeasible, total)), vjust = -0.9, size = 2.1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.85)) +
  labs(tag = "(b)", title = "Infeasibility rises with conditioning-set size",
       x = "Number of conditioning tracers", y = "Local infeasibility rate") +
  theme_hs()

p <- p_heat + p_cond +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Exploratory set-valued compatibility audit (development-stage)",
    subtitle = str_wrap("An independent interval-bound diagnostic, not the scalar benchmark above. 27.85% of eligible folds admit no feasible transit-time distribution; seven candidate explanations were tested and none accounts for it (Supplementary Information S6-S7). No cause is attributed.", 120),
    theme = theme(plot.title = element_text(size = 10, face = "bold"),
                  plot.subtitle = element_text(size = 7.5, colour = INK_MUTED))
  )

save_fig(p, "FIG-7_infeasibility_audit", width_mm = 190, height_mm = 105)
