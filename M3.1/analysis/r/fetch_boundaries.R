# One-time fetch of administrative boundaries for the M3.1 study-area map.
#
# M3 benchmarks the USGS national public-supply aquifer release (continental
# United States), not the Ghanaian field datasets used elsewhere in this
# research programme, so this fetches US boundaries rather than reusing
# M2.3/analysis/r/data/study_area_boundaries.gpkg.
#
# Boundaries come from the geoBoundaries open database (William & Mary),
# CC BY 4.0. This script is the only part of the M3.1 figure pipeline that
# touches the network; it writes a local GeoPackage so every figure script
# afterwards runs offline and reproducibly.
#
# Run from the repository root:
#   Rscript M3.1/analysis/r/fetch_boundaries.R

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(geobounds)
  library(sf)
})

OUT_DIR <- file.path("M3.1", "analysis", "r", "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
GPKG <- file.path(OUT_DIR, "study_area_boundaries.gpkg")

gb_set_cache_dir(".r-lib/geobounds-cache")

message("fetching USA national boundary")
usa <- gb_get_adm0("USA")

message("fetching USA first-level administrative regions (states)")
states <- gb_get_adm1("USA")

message("fetching neighbouring national boundaries for context")
neighbours <- do.call(rbind, lapply(c("CAN", "MEX"), function(iso) {
  x <- gb_get_adm0(iso)
  x$iso <- iso
  x[, intersect(names(x), c("shapeName", "iso", "geometry"))]
}))

keep <- function(x, cols) x[, intersect(names(x), c(cols, "geometry"))]

st_write(keep(usa, "shapeName"), GPKG, layer = "usa",
         delete_dsn = TRUE, quiet = TRUE)
st_write(keep(states, "shapeName"), GPKG, layer = "states",
         append = FALSE, quiet = TRUE)
st_write(neighbours, GPKG, layer = "neighbours",
         append = FALSE, quiet = TRUE)

message(sprintf("wrote %s (%.1f KB)", GPKG, file.size(GPKG) / 1024))
message("layers: ", paste(st_layers(GPKG)$name, collapse = ", "))
message("states: ", nrow(states))
