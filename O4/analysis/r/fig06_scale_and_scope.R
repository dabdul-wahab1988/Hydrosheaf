# FIG-6: benchmark scale and field-/archive-transfer scope across all three
# components.
source("O4/analysis/r/_theme.R")

bs <- read_export("benchmark_scale.csv") %>%
  mutate(has_field_transfer = !str_detect(field_or_archive_transfer, "^none"),
         component_group = case_when(
           str_detect(component, "^M6") ~ "M6 robustness",
           str_detect(component, "^M7") ~ "M7 identifiability",
           TRUE ~ "M8 calibration"
         ))

p <- ggplot(bs, aes(x = reorder(component, has_field_transfer), y = 1, fill = component_group, alpha = has_field_transfer)) +
  geom_col(width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = PAL_COMPONENT, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25),
                      labels = c(`TRUE` = "field/archive transfer", `FALSE` = "controlled synthetic only")) +
  labs(x = NULL, y = NULL, alpha = NULL,
       title = "Field- and archive-transfer scope is uneven across the three layers",
       subtitle = "M6 transfers to three real Ghanaian datasets; M7 uses fresh independent MODFLOW6/MODPATH7 truth throughout; M8 has no field-transfer component at all",
       caption = "Bar presence indicates the design/replicate scale reported in Table 4; shading indicates whether any field or archive dataset enters that specific benchmark row.") +
  theme_o4() +
  theme(axis.text.y = element_text(size = 6.2), legend.position = "top",
        axis.text.x = element_blank(), panel.grid.major.y = element_blank())

save_fig(p, "FIG-6", width_mm = 170, height_mm = 100)
