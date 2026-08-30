# FIG-5: M8 coverage-vs-recovery divergence under matched vs independent
# numerical truth, and the kinetic rank-one structural result.
source("O4/analysis/r/_theme.R")

mf <- read_export("calibration_model_form_shift.csv") %>%
  mutate(model_form = factor(model_form,
           levels = c("matched_analytical", "independent_numerical", "independent_numerical_baseline"),
           labels = c("Matched model\n(50d added)", "Independent truth\n(50d added)", "Independent truth\n(no new measurement)")))

p1 <- ggplot(mf, aes(x = model_form, y = internal_signal_value, fill = parameter)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "22", colour = "grey40", linewidth = 0.3) +
  annotate("text", x = 0.6, y = 0.97, label = "nominal 0.95", size = 2.2, colour = "grey40", hjust = 0) +
  scale_fill_manual(values = c(dispersivity = "#0072B2", decay = "#D55E00")) +
  labs(y = "Linearised 95% coverage\n(internal signal)", x = NULL, fill = NULL,
       title = "The same observation that helps one parameter degrades the other,\nand coverage collapses exactly where the point estimate improves",
       subtitle = "M8 optimal-design 50 d observation, matched analytical model vs independent numerical truth") +
  theme_o4() + theme(legend.position = "top")

p2 <- ggplot(mf, aes(x = model_form, y = external_signal_value, fill = parameter)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(dispersivity = "#0072B2", decay = "#D55E00"), guide = "none") +
  labs(y = "Median abs. log10 parameter error\n(external signal, vs known truth)", x = NULL,
       caption = "Dispersivity error falls under independent truth (0.826->0.167) while its own coverage collapses (0.788->0.02); decay error rises (0.137->0.154) while its coverage also collapses (0.64->0.004).") +
  theme_o4()

p <- (p1 / p2)
save_fig(p, "FIG-5", width_mm = 150, height_mm = 130)
