## M6 field-transfer benchmark — six Nature-style main figures (2-column layouts).
## Usage: Rscript plot_m6_publication_figures.R
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(forcats); library(stringr); library(scales)
})

HERE <- tryCatch(dirname(normalizePath(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(HERE) == 0 || is.na(HERE)) HERE <- getwd()
source(file.path(HERE, "theme_m6.R"))
BENCH <- normalizePath(file.path(HERE, ".."))
RES <- file.path(BENCH, "results")
OUT <- file.path(BENCH, "figures", "r_publication")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) suppressMessages(read_csv(file.path(RES, f), show_col_types = FALSE))

AQ <- c("Middle Voltaian sedimentary aquifer"="Voltaian sed.",
        "Birimian fractured basement aquifer"="Birimian bsmt.",
        "Granitoid fractured basement aquifer"="Granitoid bsmt.",
        "Regolith/alluvial shallow aquifer"="Regolith/alluv.")
AQCOL <- c("Voltaian sed."="#3B6EA8","Birimian bsmt."="#008C7A",
           "Granitoid bsmt."="#C77C2B","Regolith/alluv."="#8C6BB1")
aqf <- function(x) factor(unname(AQ[x]), levels = unname(AQ))
famlab <- function() labs(fill = "Process family")
proc_labels <- c(silicate="Silicate", carbonate="Carbonate",
                 evaporite="Evapoconc.", ion_exchange="Cation exch.",
                 anthropogenic="Nitrate", redox="Redox", trace_mineral="Fluoride",
                 none="Unresolved")

## ===================================================== FIGURE 1: design ======
fig1 <- function() {
  coords <- rd("m6_map_coordinates.csv"); ladder <- rd("m6_tier_ladder.csv")
  avail <- rd("m6_variable_availability.csv")
  dscol <- c("Northern Ghana"="#2166AC","Lower Anayari"="#008C7A","Talensi"="#C77C2B")
  coords <- coords |> mutate(ds = recode(dataset, northern_ghana="Northern Ghana",
                                         manu="Lower Anayari", talensi="Talensi"))
  a <- ggplot(coords, aes(longitude, latitude, colour = ds)) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_colour_manual(values = dscol, name = NULL) +
    scale_x_continuous(breaks = pretty_breaks(4)) +
    scale_y_continuous(breaks = pretty_breaks(4)) +
    labs(subtitle = "Sampling locations", x = "Longitude (°E)", y = "Latitude (°N)") +
    theme_m6() + theme(legend.position = c(0.02, 0.02),
                       legend.justification = c(0, 0),
                       legend.background = element_rect(fill = alpha("white", 0.7), colour = NA))

  ladder <- ladder |> mutate(dataset = factor(dataset,
              levels = c("talensi","manu","northern_ghana"),
              labels = c("Talensi","Lower Anayari","Northern Ghana")),
              tier = fct_reorder(tier, tier_index))
  b <- ggplot(ladder, aes(tier, dataset, fill = factor(attained))) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = ifelse(attained == 1, "✓", "–"),
                  colour = factor(attained)), size = 4.6, fontface = "bold") +
    scale_fill_manual(values = c("0"="#F0F0F0","1"="#2166AC"), guide = "none") +
    scale_colour_manual(values = c("0"="#BBBBBB","1"="white"), guide = "none") +
    scale_x_discrete(labels = function(x) str_replace(x, "Tier ", "T")) +
    labs(subtitle = "Evidence-tier attainment", x = NULL, y = NULL) +
    theme_m6() + theme(panel.grid = element_blank(), axis.line = element_blank(),
                       axis.ticks = element_blank())

  varlev <- c("Ca","Mg","Na","K","HCO3","Cl","SO4","NO3","F","Fe",
              "d18O","d2H","Sr_mgL","SiO2_mgL","Calcite_SI","Aquifer_Type")
  varnice <- c(Sr_mgL="Sr", SiO2_mgL="SiO2", Calcite_SI="SI", Aquifer_Type="Aquifer",
               d18O="δ18O", d2H="δ2H")
  avail <- avail |> mutate(
    variable = factor(variable, levels = varlev,
                      labels = ifelse(varlev %in% names(varnice), varnice[varlev], varlev)),
    dataset = factor(dataset, levels = c("northern_ghana","manu","talensi"),
                     labels = c("N. Ghana","L. Anayari","Talensi")))
  c <- ggplot(avail, aes(variable, fct_rev(dataset), fill = frac_present)) +
    geom_tile(colour = "white", linewidth = 0.7) +
    scale_fill_gradient(low = "#F7F7F7", high = "#08519C", limits = c(0,1),
                        name = "Fraction present", labels = percent,
                        guide = guide_colourbar(barheight = 0.5, barwidth = 8)) +
    labs(subtitle = "Variable availability across datasets", x = NULL, y = NULL) +
    theme_m6() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5),
                       panel.grid = element_blank(), axis.line = element_blank(),
                       axis.ticks = element_blank())

  p <- m6_tag((a | b) / c + plot_layout(heights = c(1.15, 0.85)), collect = FALSE)
  m6_save(p, file.path(OUT, "figure1_dataset_tier_design"), 8.8, 7.4)
}

## ================================================= FIGURE 2: workflow ========
fig2 <- function() {
  nb <- function(x,y,lab,fill,w=2.0,h=1.0) data.frame(x=x,y=y,w=w,h=h,lab=lab,fill=fill)
  stageband <- data.frame(
    x = c(1.15, 3.55, 6.05, 8.35), y = 4.15,
    lab = c("Frozen M5 assets","Ghana field data","Transfer inference","Robustness"),
    fill = c("#DDE6F0","#E1EFE9","#F4EBDB","#F7E2DD"))
  nodes <- do.call(rbind, list(
    nb(1.15, 3.2, "M5 synthetic\nPHREEQC benchmark", "#DDE6F0"),
    nb(1.15, 1.9, "Frozen inverse model\n+ MRS calibration", "#DDE6F0"),
    nb(3.55, 3.2, "3 Ghana datasets\nharmonise → mmol/L", "#E1EFE9"),
    nb(3.55, 1.9, "CBE screening +\nevidence tiers 0–4", "#E1EFE9"),
    nb(6.05, 3.2, "Transferred inverse\nreaction-class inference", "#F4EBDB"),
    nb(6.05, 1.9, "Identifiability class\n+ evidence gates", "#F4EBDB"),
    nb(8.35, 3.2, "Robustness &\nuncertainty diagnostics", "#F7E2DD"),
    nb(8.35, 1.9, "Limitation map\n(what survives)", "#F7E2DD")))
  vseg <- data.frame(x = c(1.15,3.55,6.05,8.35), y = 2.7, xend = c(1.15,3.55,6.05,8.35), yend = 2.4)
  hseg <- data.frame(x = c(2.15,4.55,7.05), y = c(2.55,2.55,2.55),
                     xend = c(2.55,5.05,7.35), yend = c(2.55,2.55,2.55))
  ggplot() +
    geom_tile(data = stageband, aes(x, y, width = 2.05, height = 0.42, fill = fill)) +
    geom_text(data = stageband, aes(x, y, label = lab), size = 3.5, fontface = "bold",
              colour = "#333333") +
    geom_tile(data = nodes, aes(x, y, width = w, height = h, fill = fill),
              colour = "#5A5A5A", linewidth = 0.4) +
    geom_text(data = nodes, aes(x, y, label = lab), size = 3.4, lineheight = 0.95) +
    geom_segment(data = vseg, aes(x=x,y=y,xend=xend,yend=yend),
                 arrow = arrow(length = unit(4,"pt"), type = "closed"),
                 linewidth = 0.35, colour = "#5A5A5A") +
    geom_segment(data = hseg, aes(x=x,y=y,xend=xend,yend=yend),
                 arrow = arrow(length = unit(4,"pt"), type = "closed"),
                 linewidth = 0.35, colour = "#5A5A5A") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.05, 9.45), ylim = c(1.25, 4.55)) +
    labs(title = "Hydrosheaf field-transfer workflow: M5 benchmark → M6 Ghana field") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(8,8,8,8)) -> p
  m6_save(p, file.path(OUT, "figure2_workflow"), 8.8, 3.7)
}

## ============================================ FIGURE 3: NG stability =========
fig3 <- function() {
  pairs <- rd("m6_ng_field_pairs.csv") |> mutate(aq = aqf(aquifer))
  cls <- rd("m6_ng_class_support.csv") |> mutate(aq = aqf(aquifer),
           resolution_class = factor(resolution_class, levels = m6_resolution_levels))
  fam <- rd("m6_ng_family_by_aquifer.csv") |> mutate(aq = aqf(aquifer))
  nlab <- pairs |> count(aq) |> mutate(lab = paste0("n=", n))

  a <- ggplot(cls, aes(aq, n, fill = resolution_class)) +
    geom_col(position = "fill", width = 0.74) +
    m6_fill_scale(labels = m6_resolution_labels, name = NULL,
                  guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Identifiability class by aquifer", x = NULL, y = "Share of wells") +
    theme_m6() + theme(axis.text.x = element_text(angle = 18, hjust = 1))

  b <- ggplot(fam, aes(aq, n, fill = dominant_family)) +
    geom_col(position = "fill", width = 0.74) +
    m6_fill_scale(name = NULL, labels = proc_labels, guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Dominant seasonal process (Cl-corrected)", x = NULL,
         y = "Share of wells") +
    theme_m6() + theme(axis.text.x = element_text(angle = 18, hjust = 1))

  c <- ggplot(pairs, aes(aq, mrs, fill = aq)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.32) +
    geom_text(data = nlab, aes(aq, 82, label = lab), inherit.aes = FALSE, size = 3.2,
              colour = "#666666") +
    scale_fill_manual(values = AQCOL, guide = "none") +
    coord_cartesian(ylim = c(55, 84)) +
    labs(subtitle = "Mechanism Resolution Score", x = NULL, y = "MRS") +
    theme_m6() + theme(axis.text.x = element_text(angle = 18, hjust = 1))

  d <- ggplot(pairs, aes(aq, support_stability, fill = aq)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.32) +
    scale_fill_manual(values = AQCOL, guide = "none") +
    coord_cartesian(ylim = c(0.8, 1.0)) +
    labs(subtitle = "Bootstrap support stability", x = NULL, y = "Mean Jaccard") +
    theme_m6() + theme(axis.text.x = element_text(angle = 18, hjust = 1))

  p <- m6_tag((a | b) / (c | d),
    caption = "All aquifers report uniformly partially identifiable at full information (Tier 4); Hydrosheaf does not over-resolve.")
  m6_save(p, file.path(OUT, "figure3_ng_stability"), 8.8, 8.0)
}

## ============================================ FIGURE 4: tier ablation ========
fig4 <- function() {
  abl <- rd("m6_tier_ablation.csv"); tr <- rd("m6_tier_ablation_transitions.csv")
  tlab <- c(tier0_majors="T0\nmajors", tier1_isotopes="T1\n+iso",
            tier2_fluoride="T2\n+F", tier3_sr_sio2="T3\n+Sr/SiO2",
            tier4_full_metadata="T4\nfull")
  ordv <- names(tlab)
  abl <- abl |> mutate(tier = factor(tier, levels = ordv, labels = tlab),
                       resolution_class = factor(resolution_class, levels = m6_resolution_levels))
  tr <- tr |> mutate(tier = factor(tier, levels = ordv, labels = tlab),
                     xi = as.integer(factor(tier, levels = tlab)))

  a <- ggplot(abl, aes(tier, fill = resolution_class)) +
    geom_bar(position = "fill", width = 0.78) +
    m6_fill_scale(labels = m6_resolution_labels, name = NULL,
                  guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Identifiability vs evidence tier", x = NULL, y = "Share of wells") +
    theme_m6()

  b <- ggplot(tr, aes(xi)) +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
             fill = "#B2182B", alpha = 0.07) +
    annotate("segment", x = 3.5, xend = 3.0, y = 0.30, yend = 0.10,
             arrow = arrow(length = unit(4,"pt"), type="closed"),
             colour = "#B2182B", linewidth = 0.4) +
    annotate("text", x = 3.5, y = 0.36, label = "Sr/SiO2 removed:\n52% lose identifiability",
             size = 3.3, colour = "#B2182B", hjust = 0.5, lineheight = 0.95) +
    geom_line(aes(y = mean_mrs/100, colour = "Mean MRS / 100"), linewidth = 0.9) +
    geom_point(aes(y = mean_mrs/100, colour = "Mean MRS / 100"), size = 2) +
    geom_line(aes(y = frac_non_identifiable, colour = "Non-identifiable"), linewidth = 0.9) +
    geom_point(aes(y = frac_non_identifiable, colour = "Non-identifiable"), size = 2) +
    scale_colour_manual(values = c("Non-identifiable"="#B2182B","Mean MRS / 100"="#2166AC"),
                        name = NULL) +
    scale_x_continuous(breaks = 1:5, labels = tlab) +
    scale_y_continuous(labels = percent, limits = c(0, 0.75),
                       sec.axis = sec_axis(~.*100, name = "Mean MRS")) +
    labs(subtitle = "Fit quality stays high while resolution collapses",
         x = NULL, y = "Fraction non-identifiable") + theme_m6()

  c <- ggplot(tr, aes(tier, frac_family_changed_vs_tier4)) +
    geom_col(fill = "#C77C2B", width = 0.66) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.05))) +
    labs(subtitle = "Process-label flips vs Tier 4", x = NULL,
         y = "Share of wells reassigned") + theme_m6()

  d <- ggplot(tr, aes(xi, mean_stability)) +
    geom_line(colour = "#008C7A", linewidth = 0.9) +
    geom_point(colour = "#008C7A", size = 2) +
    scale_x_continuous(breaks = 1:5, labels = tlab) +
    coord_cartesian(ylim = c(0.8, 1.0)) +
    labs(subtitle = "Support stability also stays high",
         x = NULL, y = "Mean Jaccard stability") + theme_m6()

  p <- m6_tag((a | b) / (c | d),
    caption = "Removing Sr/SiO2 collapses identifiability while MRS and support stability stay high.")
  m6_save(p, file.path(OUT, "figure4_tier_ablation"), 8.8, 8.0)
}

## ============================================ FIGURE 5: external =============
fig5 <- function() {
  summ <- rd("m6_external_summary.csv")
  extc <- rd("m6_external_transfer.csv") |>
    mutate(resolution_class = factor(resolution_class, levels = m6_resolution_levels),
           ds = recode(dataset, manu = "Lower Anayari\n(Tier 2)", talensi = "Talensi\n(Tier 1)"))
  summ <- summ |> mutate(is_ref = grepl("Ghana", dataset),
                         dataset = recode(dataset,
                           "manu (Tier2)" = "Lower Anayari (T2)",
                           "talensi (Tier1)" = "Talensi (T1)",
                           "N.Ghana ref (Tier1)" = "N. Ghana ref (T1)",
                           "N.Ghana ref (Tier2)" = "N. Ghana ref (T2)"),
                         dataset = fct_reorder(dataset, mean_mrs))

  a <- ggplot(summ, aes(mean_mrs, dataset, fill = is_ref)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = round(mean_mrs, 0)), hjust = 1.3, size = 3.3, colour = "white") +
    scale_fill_manual(values = c("FALSE"="#C77C2B","TRUE"="#2166AC"),
                      labels = c("External","N. Ghana ref"), name = NULL) +
    labs(subtitle = "Mean MRS vs matched-tier reference",
         x = "Mean MRS", y = NULL) + theme_m6()

  d <- ggplot(summ, aes(frac_non_identifiable, dataset, fill = is_ref)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = percent(frac_non_identifiable, accuracy = 1)),
              hjust = -0.15, size = 3.3, colour = "#333333") +
    scale_fill_manual(values = c("FALSE"="#C77C2B","TRUE"="#2166AC"),
                      labels = c("External","N. Ghana ref"), name = NULL) +
    scale_x_continuous(labels = percent, limits = c(0, 1.12),
                       breaks = c(0,0.25,0.5,0.75,1)) +
    labs(subtitle = "Non-identifiable share: transfer cost",
         x = "Fraction non-identifiable", y = NULL) + theme_m6()

  b <- ggplot(extc, aes(ds, fill = resolution_class)) +
    geom_bar(position = "fill", width = 0.6) +
    m6_fill_scale(labels = m6_resolution_labels, name = NULL,
                  guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Identifiability distribution", x = NULL, y = "Share of edges") +
    theme_m6()

  c <- ggplot(extc, aes(ds, fill = dominant_family)) +
    geom_bar(position = "fill", width = 0.6) +
    m6_fill_scale(name = NULL, labels = proc_labels, guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Dominant process", x = NULL, y = "Share of edges") + theme_m6()

  p <- m6_tag((a | b) / (d | c),
    caption = "Lower Anayari's silicate-dominant waters are 96% non-identifiable without Sr/SiO2, despite good charge balance.")
  m6_save(p, file.path(OUT, "figure5_external_transfer"), 8.8, 8.0)
}

## ============================================ FIGURE 6: limitation map =======
fig6 <- function() {
  lim <- rd("m6_limitation_map.csv"); cvb <- rd("m6_conservative_vs_bestfit.csv")
  rank_res <- c(non_identifiable=1, partially_identifiable=2, equivalence_class=3, identifiable=4)
  lim2 <- lim |> mutate(w = frac * rank_res[resolution_class]) |>
    group_by(dataset, dominant_family) |>
    summarise(score = sum(w), n = sum(n), .groups = "drop") |>
    mutate(dataset = factor(dataset, levels = c("talensi","manu","northern_ghana"),
                            labels = c("Talensi","L. Anayari","N. Ghana")),
           family = factor(unname(proc_labels[dominant_family]),
                           levels = unname(proc_labels)))
  a <- ggplot(lim2, aes(dataset, fct_rev(family), fill = score)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = n), size = 3.6, colour = "#222222") +
    scale_fill_gradientn(colours = c("#B2182B","#EF8A62","#67A9CF","#2166AC"),
                         limits = c(1,4), name = "Resolution",
                         breaks = 1:4, labels = c("non","partial","class","ident."),
                         guide = guide_colourbar(barheight = 0.5, barwidth = 9,
                                                 title.position = "top")) +
    labs(subtitle = "Process-class identifiability (cell = n wells)",
         x = NULL, y = NULL) + theme_m6() +
    theme(panel.grid = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank())

  cvb2 <- cvb |> count(bestfit_n_reactions, conservative_n_admitted)
  b <- ggplot(cvb2, aes(bestfit_n_reactions, conservative_n_admitted, size = n)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#999999") +
    geom_point(alpha = 0.55, colour = "#2166AC") +
    scale_size_continuous(range = c(1.5, 6), name = "Wells") +
    annotate("text", x = 4, y = 9.4, label = "conservative admits\nmore alternatives",
             size = 3.2, colour = "#666666", hjust = 0, lineheight = 0.95) +
    labs(subtitle = "Conservative vs single best-fit reporting",
         x = "Best-fit reactions committed", y = "Class members admitted") + theme_m6()

  p <- m6_tag((a | b) + plot_layout(widths = c(1.05, 1)))
  m6_save(p, file.path(OUT, "figure6_limitation_map"), 9.0, 4.6)
}

main <- function() {
  for (f in list(fig1, fig2, fig3, fig4, fig5, fig6)) {
    ok <- try(f(), silent = FALSE)
    if (inherits(ok, "try-error")) cat("  [warn] a figure failed\n") else cat("figure done\n")
  }
  cat("All main figures written to", OUT, "\n")
}
main()
