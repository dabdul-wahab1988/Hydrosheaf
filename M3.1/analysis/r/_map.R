# Study-area map for Figure 1: USGS national public-supply aquifer sites used
# in the M3.1 benchmark.
#
# Reads the locally cached GeoPackage written by fetch_boundaries.R, so this
# runs offline. Boundaries: geoBoundaries (William & Mary), CC BY 4.0.
# Unlike the Ghana map in M2.3/analysis/r/_map.R, the benchmark spans the
# continental United States rather than a sub-national region, so no locator
# inset is used; state outlines alone orient the reader.

suppressPackageStartupMessages({
  library(sf)
  library(ggspatial)
})

BOUNDARY_GPKG <- file.path("M3.1", "analysis", "r", "data",
                           "study_area_boundaries.gpkg")

build_study_map <- function(sites, pal_aqgroup) {
  usa    <- st_read(BOUNDARY_GPKG, layer = "usa", quiet = TRUE)
  states <- st_read(BOUNDARY_GPKG, layer = "states", quiet = TRUE)

  pts <- sites |>
    dplyr::filter(is.finite(lon), is.finite(lat)) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  # Contiguous United States only; the benchmark carries no Alaska, Hawaii or
  # territory sites (lon in [-124.3, -67.9], lat in [25.4, 48.9]), so framing
  # tighter than the full country avoids a mostly-empty canvas.
  bb <- st_bbox(pts)
  pad_x <- 1.2
  pad_y <- 1.0
  xlim <- unname(c(bb[["xmin"]] - pad_x, bb[["xmax"]] + pad_x))
  ylim <- unname(c(bb[["ymin"]] - pad_y, bb[["ymax"]] + pad_y))

  aq_levels <- names(pal_aqgroup)
  n_lab <- sites |>
    dplyr::count(AqGroup) |>
    dplyr::mutate(label = sprintf("%s (n = %d)", AqGroup, n))
  lab_lookup <- n_lab$label[match(aq_levels, n_lab$AqGroup)]

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = states, fill = "#F7F8F9", colour = "#C8CDD2",
                     linewidth = 0.25) +
    ggplot2::geom_sf(data = usa, fill = NA, colour = "#6B7280", linewidth = 0.5) +
    ggplot2::geom_sf(data = pts, ggplot2::aes(colour = AqGroup,
                                              shape = reportable_strict_parity),
                     size = 0.85, alpha = 0.80, stroke = 0.3) +
    ggplot2::scale_colour_manual(values = pal_aqgroup, breaks = aq_levels,
                                 labels = lab_lookup, name = "Aquifer lithology group") +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 4),
      labels = c(`TRUE` = "In the N = 329 reportable subset",
                 `FALSE` = "Not reportable (excluded upstream)"),
      name = NULL
    ) +
    annotation_scale(location = "bl", width_hint = 0.18, height = unit(0.10, "cm"),
                     text_cex = 0.55, line_width = 0.4,
                     pad_x = unit(0.15, "cm"), pad_y = unit(0.15, "cm")) +
    annotation_north_arrow(location = "br", which_north = "true",
                           height = unit(0.6, "cm"), width = unit(0.46, "cm"),
                           pad_x = unit(0.16, "cm"), pad_y = unit(0.18, "cm"),
                           style = north_arrow_minimal(text_size = 6)) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(
      title = "USGS public-supply sites used in the M3.1 benchmark",
      x = NULL, y = NULL
    ) +
    theme_hs() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 6.2),
      legend.title = ggplot2::element_text(size = 6.6),
      legend.key.size = ggplot2::unit(0.30, "cm"),
      panel.background = ggplot2::element_rect(fill = "#FBFCFD", colour = NA),
      panel.border = ggplot2::element_rect(fill = NA, colour = "#B9C0C7", linewidth = 0.35),
      panel.grid.major = ggplot2::element_line(colour = "#EDEFF1", linewidth = 0.2),
      axis.text = ggplot2::element_text(size = 5.6)
    )
}
