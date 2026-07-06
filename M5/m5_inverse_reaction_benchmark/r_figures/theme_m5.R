suppressPackageStartupMessages({
  library(ggplot2)
})

m5_palette <- c(
  bounded_ls = "#737373",
  lasso = "#A6A6A6",
  elastic_net = "#525252",
  thermo_elastic_net = "#8C6BB1",
  hydrosheaf_core = "#3B6EA8",
  hydrosheaf_guarded = "#008C7A",
  phreeqc_inverse = "#111111",
  core = "#737373",
  plus_lite = "#008C7A",
  enhanced = "#C77C2B",
  available_plus_lite = "#008C7A",
  carbonate = "#008C7A",
  crystalline = "#3B6EA8",
  evaporitic = "#C77C2B",
  mixed = "#8C6BB1",
  Carbonate = "#008C7A",
  Crystalline = "#3B6EA8",
  Evaporitic = "#C77C2B",
  Mixed = "#8C6BB1"
)

m5_fill_scale <- function(...) {
  scale_fill_manual(values = m5_palette, ..., drop = FALSE)
}

m5_colour_scale <- function(...) {
  scale_colour_manual(values = m5_palette, ..., drop = FALSE)
}

theme_m5 <- function(base_size = 7.2, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = grid::unit(1.6, "pt"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 0.5),
      plot.title = element_text(
        face = "bold", colour = "black", size = base_size + 1.1, hjust = 0
      ),
      plot.subtitle = element_text(colour = "#333333", size = base_size - 0.4, hjust = 0),
      plot.caption = element_text(colour = "#4D4D4D", size = base_size - 0.5, hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.2, colour = "black"),
      legend.key.height = grid::unit(8, "pt"),
      legend.key.width = grid::unit(12, "pt"),
      legend.position = "bottom",
      legend.box.spacing = grid::unit(1, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size, colour = "black"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "#E5E5E5"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(4, 5, 4, 4, "pt"),
      plot.tag = element_text(face = "bold", size = base_size + 1, colour = "black")
    )
}

theme_m5_heatmap <- function(base_size = 7.2, base_family = "Arial") {
  theme_m5(base_size, base_family) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right"
    )
}

m5_patchwork_tags <- function() {
  patchwork::plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = element_text(face = "bold", size = 9, colour = "black")
    )
}

tag_panel <- function(plot, tag) {
  plot +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(face = "bold", size = 8.5, colour = "black"),
      plot.tag.position = c(0, 1)
    )
}

save_m5 <- function(plot, path, width = 183, height = 120, dpi = 600) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!grepl("\\.png$", path, ignore.case = TRUE)) {
    stop("save_m5 expects a .png path; .pdf and .tif are written beside it.")
  }
  ggplot2::ggsave(
    path, plot, width = width, height = height, units = "mm",
    device = ragg::agg_png, dpi = dpi, bg = "white"
  )
  ggplot2::ggsave(
    sub("\\.png$", ".tif", path, ignore.case = TRUE), plot,
    width = width, height = height, units = "mm",
    device = ragg::agg_tiff, dpi = dpi, bg = "white", compression = "lzw"
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(
    sub("\\.png$", ".pdf", path, ignore.case = TRUE), plot,
    width = width, height = height, units = "mm",
    device = pdf_device, bg = "white"
  )
}
