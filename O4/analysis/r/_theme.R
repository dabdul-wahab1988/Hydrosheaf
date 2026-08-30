# Shared figure theme, palette and IO helpers for the O4 manuscript.
#
# Authority: R consumes read-only tidy CSV exports written by the Python
# layer under O4/manuscript/artifacts/data/. No figure recomputes a reported
# statistic (DECISION D4). Palette and theme follow the O3/M2.3 house style
# for visual consistency across the three companion papers.

if (dir.exists(".r-lib")) .libPaths(c(".r-lib", .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

DATA_DIR <- file.path("O4", "manuscript", "artifacts", "data")
FIG_DIR  <- file.path("O4", "manuscript", "artifacts", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Fixed component-to-colour assignment, stable across every figure.
PAL_COMPONENT <- c(
  "M6 robustness"       = "#0072B2",
  "M7 identifiability"  = "#D55E00",
  "M8 calibration"      = "#009E73"
)
PAL_SIGNAL <- c("internal" = "#5F6368", "external" = "#CC3311")
INK        <- "#1A1A1A"
INK_MUTED  <- "#5F6368"
GRID       <- "#E4E6E8"

theme_o4 <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      text             = element_text(colour = INK),
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      hjust = 0, margin = margin(b = 3)),
      plot.subtitle    = element_text(size = base_size - 0.5, colour = INK_MUTED,
                                      hjust = 0, margin = margin(b = 6)),
      plot.caption     = element_text(size = base_size - 1.5, colour = INK_MUTED,
                                      hjust = 0),
      plot.tag         = element_text(size = base_size + 2, face = "bold"),
      axis.title       = element_text(size = base_size - 0.5, colour = INK_MUTED),
      axis.text        = element_text(size = base_size - 1, colour = INK_MUTED),
      panel.grid.major = element_line(colour = GRID, linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.text       = element_text(size = base_size - 0.5, face = "bold",
                                      colour = INK, hjust = 0),
      strip.background = element_blank(),
      legend.title     = element_text(size = base_size - 1, colour = INK_MUTED),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.32, "cm"),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.15, "cm"),
      plot.margin      = margin(5, 6, 4, 5)
    )
}

read_export <- function(name) {
  path <- file.path(DATA_DIR, name)
  if (!file.exists(path)) stop("missing export: ", path, call. = FALSE)
  suppressMessages(readr::read_csv(path, show_col_types = FALSE,
                                   progress = FALSE))
}

# Journal single column is 90 mm, double column 190 mm.
save_fig <- function(plot, id, width_mm, height_mm) {
  for (dev in c("pdf", "png")) {
    ggsave(
      filename = file.path(FIG_DIR, paste0(id, ".", dev)),
      plot = plot, width = width_mm, height = height_mm, units = "mm",
      dpi = 600, device = dev, bg = "white"
    )
  }
  message(sprintf("wrote %s (%.0f x %.0f mm)", id, width_mm, height_mm))
  invisible(plot)
}
