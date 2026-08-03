# Figure 1: three independent benchmark pipelines and the common taxonomy
# they are harmonized into.
#
# Layout declared as data (tribble), not hand-placed, consistent with
# M2.3/analysis/r/fig01_architecture.R.

source(file.path("O3", "analysis", "r", "_theme.R"))

BOX_W <- 2.55
BOX_H <- 0.62

nodes <- tibble::tribble(
  ~id,        ~x,    ~y,   ~label,                                              ~lane,
  "m3ref",    0,     4.6,  "USGS national\ngroundwater-age release",            "M3",
  "m3bench",  0,     3.4,  "TracerLPM strict-parity\nemulation + graph falsification", "M3",
  "m3out",    0,     2.2,  "Reportability-guarded\nheld-out agreement",         "M3",
  "m4ref",    3.1,   4.6,  "MODPATH particle tracking\n(3 public archives)",    "M4",
  "m4bench",  3.1,   3.4,  "Independent graph inference\n+ negative controls",  "M4",
  "m4out",    3.1,   2.2,  "Edge precision / recall\n(no MODPATH access)",      "M4",
  "m5ref",    6.2,   4.6,  "240-scenario live-PHREEQC\nfactorial benchmark",    "M5",
  "m5bench",  6.2,   3.4,  "Sparse inversion vs\nconventional PHREEQC baseline","M5",
  "m5out",    6.2,   2.2,  "Phase F1 / false-discovery\nrate, all comparators", "M5",
  "tax",      3.1,   0.6,  "Common evidence taxonomy\n(independent / negative control /\nprior-informed / calibrated-emulation)", "O3",
  "pattern",  3.1,  -0.8,  "Recall exceeds precision under\nindependent, uncalibrated evaluation\nin all three layers",           "O3"
)

edges <- tibble::tribble(
  ~from,     ~to,
  "m3ref",   "m3bench",
  "m3bench", "m3out",
  "m4ref",   "m4bench",
  "m4bench", "m4out",
  "m5ref",   "m5bench",
  "m5bench", "m5out",
  "m3out",   "tax",
  "m4out",   "tax",
  "m5out",   "tax",
  "tax",     "pattern"
)

lane_fill <- c("M3" = "#DCE9F5", "M4" = "#FBE3D3", "M5" = "#D8F0E8", "O3" = "#EDEDED")

edge_df <- edges |>
  dplyr::left_join(nodes, by = c("from" = "id")) |>
  dplyr::rename(x0 = x, y0 = y) |>
  dplyr::left_join(nodes |> dplyr::select(id, x, y), by = c("to" = "id")) |>
  dplyr::rename(x1 = x, y1 = y) |>
  dplyr::mutate(
    y0 = y0 - BOX_H / 2 - 0.04,
    y1 = y1 + BOX_H / 2 + 0.04
  )

p <- ggplot() +
  geom_segment(data = edge_df, aes(x = x0, y = y0, xend = x1, yend = y1),
               colour = INK_MUTED, linewidth = 0.35,
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_tile(data = nodes, aes(x = x, y = y, fill = lane),
            width = BOX_W, height = BOX_H, colour = INK, linewidth = 0.3) +
  geom_text(data = nodes, aes(x = x, y = y, label = label),
            size = 2.35, colour = INK, lineheight = 0.88) +
  scale_fill_manual(values = lane_fill, guide = "none") +
  coord_fixed(clip = "off") +
  labs(title = "Three independently designed benchmarks, one evidence taxonomy",
       subtitle = "Each pipeline was locked on its own timetable against its own external reference; O3 adds only the bottom two layers") +
  theme_o3(base_size = 9) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

save_fig(p, "FIG-1", width_mm = 190, height_mm = 130)
