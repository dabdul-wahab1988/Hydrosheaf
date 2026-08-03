# Figure 1 (new in M3.1): geographic distribution of the USGS national
# public-supply sites used in the benchmark, by aquifer lithology group, with
# the N = 329 strict reportable subset marked. M3 had no map figure.

source(file.path("M3.1", "analysis", "r", "_theme.R"))
source(file.path("M3.1", "analysis", "r", "_map.R"))

sites <- read_export("site_locations.csv")

p <- build_study_map(sites, PAL_AQGROUP)

save_fig(p, "FIG-1_site_map", width_mm = 190, height_mm = 115)
