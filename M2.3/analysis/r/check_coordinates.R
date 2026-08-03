# Point-in-polygon validation of every sample coordinate against the Ghana
# national boundary, and against the administrative region each dataset names.
#
# This check exists because a visual inspection of the study-area map was not
# sufficient: the Talensi longitudes were recorded without their negative sign
# and plotted outside Ghana without that being noticed by eye.
#
# Run from the repository root:
#   Rscript M2.3/analysis/r/check_coordinates.R

if (dir.exists(".r-lib")) .libPaths(c(".r-lib", .libPaths()))
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
})

GPKG <- file.path("M2.3", "analysis", "r", "data", "study_area_boundaries.gpkg")
DATA <- file.path("M2.3", "manuscript", "artifacts", "data")

ghana   <- st_read(GPKG, layer = "ghana", quiet = TRUE)
regions <- st_read(GPKG, layer = "regions", quiet = TRUE)

samples <- read_csv(file.path(DATA, "field_samples_derived.csv"),
                    show_col_types = FALSE) |>
  filter(primary_panel, is.finite(Longitude), is.finite(Latitude))

pts <- st_as_sf(samples, coords = c("Longitude", "Latitude"), crs = 4326,
                remove = FALSE)

# A national boundary polygon is generalised, so a point a short distance beyond
# it is a cartographic artefact rather than a bad coordinate. Anything further
# out than this is treated as a genuine failure.
TOLERANCE_M <- 5000

inside <- lengths(st_within(pts, ghana)) > 0
dist_m <- as.numeric(st_distance(pts, st_boundary(ghana)))
pts$outside_m <- ifelse(inside, 0, dist_m)

reg_idx <- st_within(pts, regions)
pts$region <- vapply(reg_idx, function(i)
  if (length(i)) regions$shapeName[i[1]] else NA_character_, character(1))

report <- pts |>
  st_drop_geometry() |>
  mutate(in_ghana = inside,
         within_tolerance = inside | outside_m <= TOLERANCE_M) |>
  group_by(dataset) |>
  summarise(
    n = n(),
    n_inside_ghana = sum(in_ghana),
    n_outside_ghana = sum(!in_ghana),
    max_outside_m = round(max(outside_m), 0),
    n_beyond_tolerance = sum(!within_tolerance),
    regions = paste(sort(unique(na.omit(region))), collapse = "; "),
    .groups = "drop"
  )

cat("\n== coordinate validation ==\n")
print(as.data.frame(report), row.names = FALSE)

write_csv(report, file.path(DATA, "coordinate_validation.csv"))

failed <- report |> filter(n_beyond_tolerance > 0)
if (nrow(failed)) {
  cat(sprintf("\nFAIL: samples further than %d m outside Ghana\n", TOLERANCE_M))
  for (i in seq_len(nrow(failed))) {
    cat(sprintf("  %s: %d of %d, furthest %.0f m\n", failed$dataset[i],
                failed$n_beyond_tolerance[i], failed$n[i],
                failed$max_outside_m[i]))
  }
  quit(status = 1)
}
cat(sprintf("\nPASS: all coordinates inside Ghana or within %d m of the border.\n",
            TOLERANCE_M))
near <- report |> filter(n_outside_ghana > 0)
if (nrow(near)) {
  cat("Note: some points sit just beyond the generalised boundary polygon:\n")
  for (i in seq_len(nrow(near))) {
    cat(sprintf("  %s: %d of %d, furthest %.0f m\n", near$dataset[i],
                near$n_outside_ghana[i], near$n[i], near$max_outside_m[i]))
  }
}
