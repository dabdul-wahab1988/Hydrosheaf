# FIG-4: M7 entropy reduction (internal) vs predictive-skill change
# (external, PR-AUC), native evidence combination and three adverse
# controls, with 95% bootstrap intervals.
source("O4/analysis/r/_theme.R")

iv <- read_export("integration_value.csv") %>%
  mutate(is_adverse = str_detect(axis_x, "adverse"),
         axis_x = factor(axis_x, levels = rev(axis_x)))

p1 <- ggplot(iv, aes(x = axis_x, y = internal_signal_value, colour = is_adverse)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_pointrange(aes(ymin = internal_ci_low, ymax = internal_ci_high), fatten = 2.4) +
  coord_flip() +
  scale_colour_manual(values = c(`FALSE` = PAL_COMPONENT["M7 identifiability"], `TRUE` = "#8B0000"), guide = "none") +
  labs(y = "Mean edge-entropy change (nats)\ninternal signal", x = NULL,
       title = "Entropy falls under every condition, including the adverse controls",
       subtitle = "M7.3 native evidence increments and permuted-stream adverse control, 95% case-block bootstrap CI") +
  theme_o4()

p2 <- ggplot(iv, aes(x = axis_x, y = external_signal_value, colour = is_adverse)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_pointrange(aes(ymin = external_ci_low, ymax = external_ci_high), fatten = 2.4) +
  coord_flip() +
  scale_colour_manual(values = c(`FALSE` = PAL_COMPONENT["M7 identifiability"], `TRUE` = "#8B0000"), guide = "none") +
  labs(y = "PR-AUC change vs independent MODFLOW6/MODPATH7 truth\nexternal signal", x = NULL,
       caption = "Red = adverse control (age permuted); blue = native evidence increment. Entropy reduction and predictive-skill change diverge in sign for age (both native and permuted) and, more sharply, for the joint-misspecification control.") +
  theme_o4()

p <- (p1 / p2)
save_fig(p, "FIG-4", width_mm = 160, height_mm = 130)
