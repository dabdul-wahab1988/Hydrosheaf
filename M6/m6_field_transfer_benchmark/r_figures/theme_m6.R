## Shared Nature-style theme and palettes for the M6 field-transfer figures.
suppressPackageStartupMessages({
  library(ggplot2)
})

## Resolution-class palette (ordered, colour-blind safe)
m6_resolution_levels <- c("non_identifiable", "partially_identifiable",
                          "equivalence_class", "identifiable")
m6_resolution_labels <- c(
  non_identifiable = "Non-identifiable",
  partially_identifiable = "Partially identifiable",
  equivalence_class = "Equivalence class",
  identifiable = "Identifiable"
)
m6_palette <- c(
  ## resolution classes
  non_identifiable = "#B2182B",
  partially_identifiable = "#EF8A62",
  equivalence_class = "#67A9CF",
  identifiable = "#2166AC",
  ## process families
  silicate = "#3B6EA8",
  carbonate = "#008C7A",
  evaporite = "#C77C2B",
  ion_exchange = "#8C6BB1",
  anthropogenic = "#D6604D",
  redox = "#7F7F7F",
  trace_mineral = "#B58900",
  none = "#CCCCCC",
  ## datasets
  northern_ghana = "#2166AC",
  talensi = "#C77C2B",
  manu = "#008C7A",
  ## edge sets
  provided_graph = "#2166AC",
  chemistry_knn = "#008C7A",
  geographic_nearest = "#C77C2B",
  random_perturbed = "#B2182B"
)

m6_fill_scale <- function(...) scale_fill_manual(values = m6_palette, ..., drop = FALSE)
m6_colour_scale <- function(...) scale_colour_manual(values = m6_palette, ..., drop = FALSE)

theme_m6 <- function(base_size = 13, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = grid::unit(2.5, "pt"),
      axis.text = element_text(colour = "black", size = base_size - 1.5),
      axis.title = element_text(colour = "black", size = base_size),
      plot.title = element_text(face = "bold", colour = "black",
                                size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(colour = "#2A2A2A", size = base_size - 0.5,
                                   hjust = 0, margin = margin(b = 4)),
      plot.tag = element_text(face = "bold", size = base_size + 5, colour = "black"),
      plot.tag.position = c(0.006, 0.99),
      plot.caption = element_text(colour = "#4D4D4D", size = base_size - 2.5, hjust = 0),
      plot.margin = margin(8, 10, 6, 6),
      legend.title = element_text(size = base_size - 1, colour = "black"),
      legend.text = element_text(size = base_size - 1.5, colour = "black"),
      legend.key.height = grid::unit(12, "pt"),
      legend.key.width = grid::unit(16, "pt"),
      legend.position = "bottom",
      legend.box.spacing = grid::unit(3, "pt"),
      legend.margin = margin(2, 4, 2, 4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1, colour = "black"),
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "#ECECEC"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

## Common patchwork assembly: collect legends into one aligned bottom strip,
## add bold a/b tags + shared caption.
m6_tag <- function(p, caption = NULL, collect = TRUE) {
  if (collect) p <- p + patchwork::plot_layout(guides = "collect")
  p + patchwork::plot_annotation(
    tag_levels = "a",
    caption = caption,
    theme = theme(plot.tag = element_text(face = "bold", size = 18),
                  plot.caption = element_text(size = 10.5, colour = "#4D4D4D", hjust = 0),
                  legend.position = "bottom", legend.box = "vertical",
                  legend.box.just = "center", legend.spacing.y = grid::unit(1, "pt")))
}

## Save a figure at 300 dpi in png/pdf/tif (Nature portfolio formats).
m6_save <- function(plot, path_noext, width, height) {
  dir.create(dirname(path_noext), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(path_noext, ".png"), plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  ggsave(paste0(path_noext, ".pdf"), plot, width = width, height = height,
         units = "in", bg = "white")
  ok <- try(ggsave(paste0(path_noext, ".tif"), plot, width = width, height = height,
                   units = "in", dpi = 300, bg = "white",
                   compression = "lzw"), silent = TRUE)
  invisible(plot)
}

m6_reslab <- function(x) factor(m6_resolution_labels[x],
                                levels = unname(m6_resolution_labels))
