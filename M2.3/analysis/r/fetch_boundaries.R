# One-time fetch of administrative boundaries for the study-area map.
#
# Boundaries come from the geoBoundaries open database (William & Mary), which
# is CC BY 4.0. This script is the only part of the pipeline that touches the
# network. It writes a local GeoPackage so that every figure script afterwards
# runs offline and reproducibly.
#
# Run from the repository root:
#   Rscript M2.3/analysis/r/fetch_boundaries.R

.libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(geobounds)
  library(sf)
})

OUT_DIR <- file.path("M2.3", "analysis", "r", "data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
GPKG <- file.path(OUT_DIR, "study_area_boundaries.gpkg")

gb_set_cache_dir(".r-lib/geobounds-cache")

message("fetching Ghana national boundary")
ghana <- gb_get_adm0("GHA")

message("fetching Ghana first-level administrative regions")
regions <- gb_get_adm1("GHA")

message("fetching neighbouring national boundaries for context")
neighbours <- do.call(rbind, lapply(c("BFA", "CIV", "TGO"), function(iso) {
  x <- gb_get_adm0(iso)
  x$iso <- iso
  x[, intersect(names(x), c("shapeName", "iso", "geometry"))]
}))

keep <- function(x, cols) x[, intersect(names(x), c(cols, "geometry"))]

st_write(keep(ghana, "shapeName"), GPKG, layer = "ghana",
         delete_dsn = TRUE, quiet = TRUE)
st_write(keep(regions, "shapeName"), GPKG, layer = "regions",
         append = FALSE, quiet = TRUE)
st_write(neighbours, GPKG, layer = "neighbours",
         append = FALSE, quiet = TRUE)

message(sprintf("wrote %s (%.1f KB)", GPKG, file.size(GPKG) / 1024))
message("layers: ", paste(st_layers(GPKG)$name, collapse = ", "))
message("regions: ", nrow(regions))
print(sort(regions$shapeName))
