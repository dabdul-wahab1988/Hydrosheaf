## Supplementary map for a representative locked MODFLOW/MODPATH realization.
## This is a synthetic model-domain map in metres, not a geographic site map.

args <- commandArgs(FALSE)
HERE <- tryCatch(dirname(normalizePath(sub("--file=", "",
        grep("--file=", args, value = TRUE)))), error = function(e) getwd())
if (length(HERE) == 0 || is.na(HERE)) HERE <- getwd()
LOCAL_LIB <- normalizePath(file.path(HERE, "..", "..", "..", ".r-lib"),
                           mustWork = FALSE)
if (dir.exists(LOCAL_LIB)) .libPaths(c(LOCAL_LIB, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
})

BENCH <- normalizePath(file.path(HERE, ".."))
CASE <- file.path(BENCH, "results", "supporting_validation", "cases",
                  "locked_test_4101")
OUT <- file.path(BENCH, "figures", "supporting_validation")

required <- c("blind_observations.csv", "heldout_truth.csv",
              "modpath_pathline_truth.csv")
missing <- required[!file.exists(file.path(CASE, required))]
if (length(missing)) stop("Missing locked case files: ", paste(missing, collapse = ", "))

nodes <- read_csv(file.path(CASE, "blind_observations.csv"), show_col_types = FALSE) |>
  mutate(path = sub("^MF4101_(P[0-9]+)_.*$", "\\1", site_id),
         milestone = as.integer(sub("^.*_M([0-9]+)$", "\\1", site_id)))
truth <- read_csv(file.path(CASE, "heldout_truth.csv"), show_col_types = FALSE)
pathlines <- read_csv(file.path(CASE, "modpath_pathline_truth.csv"),
                      show_col_types = FALSE) |>
  arrange(particle, milestone)

edges <- truth |>
  left_join(nodes |> select(site_id, x_m, y_m), by = c("u" = "site_id")) |>
  rename(x = x_m, y = y_m) |>
  left_join(nodes |> select(site_id, x_m, y_m), by = c("v" = "site_id")) |>
  rename(xend = x_m, yend = y_m)

starts <- nodes |> filter(milestone == min(milestone))

p <- ggplot() +
  annotate("rect", xmin = 0, xmax = 2700, ymin = 0, ymax = 1500,
           fill = "#F6F8F9", colour = "#4D4D4D", linewidth = 0.5) +
  geom_vline(xintercept = seq(300, 2400, 300), colour = "white", linewidth = 0.35) +
  geom_hline(yintercept = seq(300, 1200, 300), colour = "white", linewidth = 0.35) +
  geom_path(data = pathlines,
            aes(x_m, y_m, group = particle),
            colour = "#6A6A6A", linewidth = 1.15, alpha = 0.32,
            lineend = "round") +
  geom_segment(data = edges,
               aes(x = x, y = y, xend = xend, yend = yend),
               colour = "#252525", linewidth = 0.55,
               arrow = arrow(length = grid::unit(5, "pt"), type = "closed"),
               lineend = "round") +
  geom_point(data = nodes,
             aes(x_m, y_m, fill = hydraulic_head, shape = path),
             size = 4.1, colour = "white", stroke = 0.65) +
  geom_point(data = nodes,
             aes(x_m, y_m, shape = path),
             size = 4.1, colour = "#303030", stroke = 0.22, fill = NA) +
  geom_label(data = starts,
             aes(x_m, y_m, label = path), nudge_x = -85, nudge_y = 65,
             size = 3.0, fontface = "bold", linewidth = 0.2,
             label.padding = grid::unit(0.10, "lines"), fill = alpha("white", 0.88)) +
  annotate("segment", x = 120, xend = 620, y = 95, yend = 95,
           linewidth = 1.4, colour = "#252525") +
  annotate("segment", x = c(120, 620), xend = c(120, 620),
           y = 75, yend = 115, linewidth = 0.45) +
  annotate("text", x = c(120, 620), y = 145,
           label = c("0", "500 m"), size = 3.0, hjust = c(0.5, 0.5)) +
  scale_fill_gradientn(colours = c("#440154", "#2A788E", "#7AD151", "#FDE725"),
                       name = "Hydraulic head (m)", guide = guide_colourbar(
                         title.position = "top", barwidth = grid::unit(42, "mm"),
                         barheight = grid::unit(3.5, "mm"))) +
  scale_shape_manual(values = c(P0 = 21, P1 = 22, P2 = 24),
                     name = "Flow path") +
  scale_x_continuous(breaks = seq(0, 2500, 500), labels = label_number(scale = 0.001),
                     limits = c(0, 2700), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(0, 1500, 500), labels = label_number(scale = 0.001),
                     limits = c(0, 1500), expand = expansion(mult = 0)) +
  coord_equal(clip = "on") +
  labs(title = "Locked synthetic flow domain",
       subtitle = "MODFLOW/MODPATH truth network, representative locked realization 4101",
       x = "Model x (km)", y = "Model y (km)",
       caption = paste0(
         "Arrows show the nine held-out directed flow edges; nodes are coloured by simulated hydraulic head.\n",
         "Coordinates are synthetic model space and are not georeferenced.")) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "#404040", size = 10.5,
                                     margin = margin(b = 7)),
        plot.caption = element_text(colour = "#4D4D4D", size = 8.2, hjust = 0,
                                    margin = margin(t = 7)),
        axis.title = element_text(colour = "#202020"),
        axis.text = element_text(colour = "#202020"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 9),
        plot.margin = margin(8, 10, 6, 8),
        plot.background = element_rect(fill = "white", colour = NA)) +
  guides(shape = guide_legend(title.position = "top", override.aes = list(size = 3.8)))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(OUT, "figure_s1_model_domain_map")
ggsave(paste0(stem, ".pdf"), p, width = 7.08, height = 4.7,
       units = "in", device = grDevices::cairo_pdf, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(paste0(stem, ".svg"), p, width = 7.08, height = 4.7,
         units = "in", device = svglite::svglite, bg = "white")
}
ggsave(paste0(stem, ".png"), p, width = 7.08, height = 4.7,
       units = "in", dpi = 600, bg = "white")
ggsave(paste0(stem, ".tif"), p, width = 7.08, height = 4.7,
       units = "in", dpi = 300, compression = "lzw", bg = "white")

message("M7 synthetic-domain map -> ", stem)
