# FIG-1: three stress-test pipelines and the common internal-vs-external
# signal taxonomy (schematic; content sourced from taxonomy.csv row counts).
source("O4/analysis/r/_theme.R")

tax <- read_export("taxonomy.csv")

box <- tibble::tibble(
  component = c("M6 robustness", "M7 identifiability", "M8 calibration"),
  x = c(1, 2, 3),
  n_experiments = sapply(component, function(cc) sum(tax$component == substr(cc, 1, 2))),
  stress = c("Data limitation\n(evidence-tier ablation)",
             "Evidence misspecification\n(adverse controls)",
             "Model-form shift\n(independent numerical truth)"),
  internal = c("Mean MRS\n(fit quality)", "Posterior-entropy\nreduction", "Nominal 95%\ninterval coverage"),
  external = c("True identifiability\nclass", "PR-AUC / Brier / log-loss\nvs MODFLOW6/MODPATH7", "Realised parameter-\nrecovery error")
)

p <- ggplot(box, aes(x = x)) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 2.55, ymax = 3.15,
                fill = component), colour = "white", linewidth = 0.4, alpha = 0.92) +
  geom_text(aes(y = 2.85, label = component), colour = "white", fontface = "bold", size = 3.0) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 1.65, ymax = 2.35),
            fill = "grey96", colour = "grey70", linewidth = 0.3) +
  geom_text(aes(y = 2.0, label = stress), size = 2.5, lineheight = 0.95) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.95, ymax = 1.45),
            fill = "grey92", colour = PAL_SIGNAL["internal"], linewidth = 0.4) +
  geom_text(aes(y = 1.2, label = internal), size = 2.4, lineheight = 0.9) +
  geom_rect(aes(xmin = x - 0.42, xmax = x + 0.42, ymin = 0.15, ymax = 0.65),
            fill = "grey92", colour = PAL_SIGNAL["external"], linewidth = 0.4) +
  geom_text(aes(y = 0.4, label = external), size = 2.4, lineheight = 0.9) +
  geom_segment(aes(x = x, xend = x, y = 1.65, yend = 1.46), arrow = arrow(length = unit(1.5, "mm")), colour = "grey50") +
  geom_segment(aes(x = x, xend = x, y = 0.95, yend = 0.66), arrow = arrow(length = unit(1.5, "mm")), colour = "grey50") +
  annotate("segment", x = 0.5, xend = 3.5, y = 1.5, yend = 1.5, linetype = "22", colour = "grey60", linewidth = 0.3) +
  annotate("text", x = 0.56, y = 1.55, label = "internal signal (this paper's question)",
           hjust = 0, size = 2.2, colour = INK_MUTED, fontface = "italic") +
  annotate("text", x = 0.56, y = 0.10, label = "external / ground-truth signal", hjust = 0,
           size = 2.2, colour = INK_MUTED, fontface = "italic") +
  scale_fill_manual(values = PAL_COMPONENT, guide = "none") +
  coord_cartesian(xlim = c(0.45, 3.55), ylim = c(0, 3.25), expand = FALSE) +
  labs(title = "Three independently locked stress tests, one shared taxonomy",
       subtitle = "Each layer pairs a signal computed without reference to ground truth against a signal verified against it") +
  theme_o4() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

save_fig(p, "FIG-1", width_mm = 190, height_mm = 95)
