# Study-area map for Figure 2(a).
#
# Reads the locally cached GeoPackage written by fetch_boundaries.R, so this
# runs offline. Boundaries: geoBoundaries (William & Mary), CC BY 4.0.

suppressPackageStartupMessages({
  library(sf)
  library(ggspatial)
})

BOUNDARY_GPKG <- file.path("M2.3", "analysis", "r", "data",
                           "study_area_boundaries.gpkg")

# Regions represented in the sampled datasets.
STUDY_REGIONS <- c("Northern Region", "North East Region",
                   "Upper East Region", "Upper West Region",
                   "Savannah Region", "Oti Region")

build_study_map <- function(samples, pal_dataset, ds_levels) {
  ghana      <- st_read(BOUNDARY_GPKG, layer = "ghana", quiet = TRUE)
  regions    <- st_read(BOUNDARY_GPKG, layer = "regions", quiet = TRUE)
  neighbours <- st_read(BOUNDARY_GPKG, layer = "neighbours", quiet = TRUE)

  regions$in_study <- regions$shapeName %in% STUDY_REGIONS

  pts <- samples |>
    dplyr::filter(is.finite(Longitude), is.finite(Latitude)) |>
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326,
                 remove = FALSE)

  # Frame the main panel on the sampled area with a small margin.
  bb <- st_bbox(pts)
  pad_x <- 0.35
  pad_y <- 0.30
  # unname: st_bbox returns a named vector and the names propagate into
  # st_bbox() below, which then fails to build the inset frame polygon.
  xlim <- unname(c(bb[["xmin"]] - pad_x, bb[["xmax"]] + pad_x))
  ylim <- unname(c(bb[["ymin"]] - pad_y, bb[["ymax"]] + pad_y))

  n_lab <- samples |>
    dplyr::count(dataset) |>
    dplyr::mutate(label = sprintf("%s (n = %d)", dataset, n))
  lab_lookup <- n_lab$label[match(ds_levels, n_lab$dataset)]

  # Region name anchors, placed at the polygon point-on-surface.
  reg_lab <- regions[regions$in_study, ]
  reg_ctr <- suppressWarnings(st_point_on_surface(st_geometry(reg_lab)))
  reg_lab <- cbind(
    sf::st_drop_geometry(reg_lab),
    as.data.frame(st_coordinates(reg_ctr))
  )
  reg_lab$short <- sub(" Region$", "", reg_lab$shapeName)
  reg_lab <- reg_lab[reg_lab$X > xlim[1] & reg_lab$X < xlim[2] &
                     reg_lab$Y > ylim[1] & reg_lab$Y < ylim[2], ]
  # The Upper East centroid coincides with the Talensi cluster, so that one
  # label is shifted west into open ground.
  ue <- reg_lab$short == "Upper East"
  reg_lab$X[ue] <- reg_lab$X[ue] - 0.42
  reg_lab$Y[ue] <- reg_lab$Y[ue] + 0.12

  main <- ggplot() +
    geom_sf(data = neighbours, fill = "#F2F3F4", colour = "#C8CDD2",
            linewidth = 0.25) +
    geom_sf(data = regions, fill = "#FFFFFF", colour = "#C8CDD2",
            linewidth = 0.25) +
    geom_sf(data = regions[regions$in_study, ], fill = "#F7F4EC",
            colour = "#A9B0B7", linewidth = 0.35) +
    geom_sf(data = ghana, fill = NA, colour = "#6B7280", linewidth = 0.55) +
    geom_text(data = reg_lab, aes(X, Y, label = short), size = 1.9,
              colour = "#8A9199", fontface = "italic") +
    geom_sf(data = pts, aes(colour = dataset), size = 1.05, alpha = 0.85,
            stroke = 0) +
    scale_colour_manual(values = pal_dataset, breaks = ds_levels,
                        labels = lab_lookup, name = NULL) +
    annotation_scale(location = "br", width_hint = 0.24, height = unit(0.10, "cm"),
                     text_cex = 0.42, line_width = 0.4,
                     pad_x = unit(0.12, "cm"), pad_y = unit(0.12, "cm")) +
    annotation_north_arrow(location = "tr", which_north = "true",
                           height = unit(0.52, "cm"), width = unit(0.40, "cm"),
                           pad_x = unit(0.14, "cm"), pad_y = unit(0.16, "cm"),
                           style = north_arrow_minimal(text_size = 5)) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    # Tag is set here rather than by plot_annotation, because the inset makes
    # this panel a nested patchwork that automatic tagging would skip.
    labs(title = "Sampling locations, northern Ghana",
         tag = "(a)", x = NULL, y = NULL) +
    theme_hs() +
    # Legend sits top-left, the emptiest corner of the sampled area; the inset
    # goes bottom-right so that neither covers the Talensi cluster.
    theme(legend.position = c(0.015, 0.985),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.85),
                                           colour = "#E4E6E8", linewidth = 0.2),
          legend.margin = margin(2, 3, 2, 3),
          legend.text = element_text(size = 5.4),
          legend.key.size = unit(0.26, "cm"),
          panel.background = element_rect(fill = "#FBFCFD", colour = NA),
          panel.border = element_rect(fill = NA, colour = "#B9C0C7",
                                      linewidth = 0.35),
          panel.grid.major = element_line(colour = "#EDEFF1", linewidth = 0.2),
          axis.text = element_text(size = 4.8))

  # Inset: Ghana within its immediate neighbours, with the study frame marked.
  frame <- st_as_sfc(st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                               ymin = ylim[1], ymax = ylim[2]), crs = 4326))
  inset <- ggplot() +
    geom_sf(data = neighbours, fill = "#F2F3F4", colour = "#C8CDD2",
            linewidth = 0.2) +
    geom_sf(data = ghana, fill = "#E3EAF2", colour = "#6B7280",
            linewidth = 0.35) +
    geom_sf(data = frame, fill = NA, colour = "#D55E00", linewidth = 0.5) +
    annotate("text", x = -1.05, y = 6.6, label = "Ghana", size = 1.9,
             colour = "#374151") +
    coord_sf(xlim = c(-6.2, 3.0), ylim = c(4.2, 15.2), expand = FALSE) +
    theme_void() +
    theme(panel.background = element_rect(fill = "white", colour = "#B9C0C7",
                                          linewidth = 0.3),
          plot.tag = element_blank(),
          plot.margin = margin(0, 0, 0, 0))

  # Bottom-left is the only corner carrying neither sample points nor a region
  # label. ignore_tag keeps patchwork from spending a panel letter on the inset.
  main + patchwork::inset_element(inset, left = 0.012, bottom = 0.012,
                                  right = 0.245, top = 0.40,
                                  align_to = "panel", ignore_tag = TRUE)
}
