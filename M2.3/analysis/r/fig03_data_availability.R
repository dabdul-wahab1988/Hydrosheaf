# Figure 3: what each dataset measures, and the evidence ceiling that follows.
#
# Panel (a) is a presence matrix over the canonical variable schema. Because
# almost every carried variable is fully populated, the informative contrast is
# carried / partly carried / not carried rather than a continuous gradient.
# Panel (b) states which claims the available variables can and cannot support.

source(file.path("M2.3", "analysis", "r", "_theme.R"))

inv <- read_export("field_variable_inventory.csv")
ds_levels <- c("Northern Ghana", "Lower Anayari", "Talensi")
grp_levels <- c("physical setting", "field", "major ion", "minor ion",
                "trace/corroborating", "stable isotope")
grp_short <- c("physical setting" = "Setting", "field" = "Field",
               "major ion" = "Major ions", "minor ion" = "Minor",
               "trace/corroborating" = "Trace", "stable isotope" = "Isotopes")

comp <- inv |>
  filter(variable_group %in% grp_levels) |>
  mutate(
    dataset = factor(dataset, levels = rev(ds_levels)),
    variable_group = factor(grp_short[variable_group],
                            levels = unname(grp_short[grp_levels])),
    status = case_when(
      !present            ~ "Not carried",
      pct_complete >= 100 ~ "Complete",
      TRUE                ~ "Partly complete"),
    status = factor(status, levels = c("Complete", "Partly complete",
                                       "Not carried")),
    cell_label = ifelse(present & pct_complete < 100,
                        sprintf("%.0f", pct_complete), "")) |>
  arrange(variable_group, variable) |>
  mutate(label = factor(label, levels = unique(label)))

STATUS_FILL <- c("Complete" = "#0072B2", "Partly complete" = "#E69F00",
                 "Not carried" = "#EDEFF1")

p_comp <- ggplot(comp, aes(x = label, y = dataset, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.9) +
  geom_text(aes(label = cell_label), size = 2.0, colour = INK) +
  scale_fill_manual(values = STATUS_FILL, name = NULL, drop = FALSE) +
  facet_grid(~ variable_group, scales = "free_x", space = "free_x") +
  labs(title = "Variable availability across the measured datasets",
       subtitle = paste("Numbers give percentage completeness where a carried",
                        "variable is incomplete"),
       x = NULL, y = NULL) +
  theme_hs() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 5.6),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        panel.spacing = unit(0.2, "lines"),
        plot.subtitle = element_text(size = 7),
        plot.margin = margin(5, 14, 4, 5))

# ---- (b) evidence ceiling ---------------------------------------------------
ceiling_tbl <- tibble::tribble(
  ~requirement,                      ~status,      ~supports,
  "Major-ion chemistry",             "available",  "Chemistry and screening diagnostics",
  "Stable isotopes",                 "available",  "Source and mixing plausibility only",
  "Charge-balance screening",        "available",  "Analytical quality tiering",
  "Fluoride",                        "partial",    "Absent from Talensi",
  "Strontium and silica",            "partial",    "Northern Ghana only",
  "Redox potential",                 "partial",    "Talensi only, 86% complete",
  "Environmental age tracers",       "absent",     "No field age claim",
  "Screen intervals",                "absent",     "No screen-resolved exchange",
  "Repeated hydraulic heads",        "absent",     "No head-series validation",
  "Independent connectivity truth",  "absent",     "No directed-path validation",
  "Reaction-mechanism truth",        "absent",     "No unique-mechanism claim"
) |>
  mutate(requirement = factor(requirement, levels = rev(requirement)),
         status = factor(status, levels = c("available", "partial", "absent")))

p_ceiling <- ggplot(ceiling_tbl, aes(x = 1, y = requirement, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.9, width = 0.5) +
  geom_text(aes(x = 1.30, label = supports), hjust = 0, size = 2.15,
            colour = INK_MUTED) +
  scale_fill_manual(values = c(available = "#0072B2", partial = "#E69F00",
                               absent = "#D55E00"), name = NULL,
                    labels = c("Available", "Partial", "Absent")) +
  scale_x_continuous(limits = c(0.72, 2.35), expand = c(0, 0)) +
  labs(title = "Evidence ceiling for the field application",
       x = NULL, y = NULL) +
  theme_hs() +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top", legend.justification = "left")

fig <- p_comp / p_ceiling +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_fig(fig, "FIG-3_data_availability", width_mm = 190, height_mm = 152)
