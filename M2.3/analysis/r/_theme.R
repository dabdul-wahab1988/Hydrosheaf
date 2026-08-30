# Shared figure theme, palette and IO helpers for the M2.3 manuscript.
#
# Authority: R consumes read-only tidy CSV exports written by the Python layer
# under M2.3/manuscript/artifacts/data/. No figure recomputes a reported
# statistic; figures render values that Python already derived (DECISION D6).
#
# Palette: Okabe-Ito subset, ordered by contrast against a white print surface.
# Validated with the dataviz six-check validator (light mode, categorical):
# lightness band PASS, chroma floor PASS, CVD separation PASS (worst adjacent
# pair dE 17.9 deutan, above the 12 target). The contrast check WARNs for the
# three lightest hues, which obliges visible labels; every figure below carries
# either a legend or direct labels, so identity is never colour-alone.

# Project-local library carries sf, ggspatial and geobounds.
if (dir.exists(".r-lib")) .libPaths(c(".r-lib", .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
})

DATA_DIR <- file.path("M2.3", "manuscript", "artifacts", "data")
FIG_DIR  <- file.path("M2.3", "manuscript", "artifacts", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Fixed hue order. Never cycled; a category beyond slot 6 folds into "Other".
PAL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")

# Stable entity-to-colour assignment, so a filtered panel never repaints.
PAL_DATASET <- c(
  "Northern Ghana" = "#0072B2",
  "Lower Anayari"  = "#D55E00",
  "Talensi"        = "#009E73"
)
PAL_TIER <- c(
  "quantitative" = "#0072B2",
  "screening"    = "#E69F00",
  "exploratory"  = "#D55E00",
  "unclassified" = "#9AA0A6"
)
# Sequential ramp: one hue, light to dark, for magnitude only.
SEQ_LOW  <- "#DCE9F5"
SEQ_HIGH <- "#00436B"

INK        <- "#1A1A1A"
INK_MUTED  <- "#5F6368"
GRID       <- "#E4E6E8"

theme_hs <- function(base_size = 9) {
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
