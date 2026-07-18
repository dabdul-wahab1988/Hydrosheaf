## M6 supplementary figures (S1–S10), Nature style — enlarged fonts, collected
## legends, and multi-panel composition where it adds value.
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
OUT <- file.path(BENCH, "figures", "r_publication", "supplementary")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) suppressMessages(read_csv(file.path(RES, f), show_col_types = FALSE))

dscol <- c("Talensi"="#C77C2B","Lower Anayari"="#008C7A","Northern Ghana"="#2166AC")
dslab <- function(x) factor(recode(x, talensi="Talensi", manu="Lower Anayari",
                                   northern_ghana="Northern Ghana"),
                            levels = c("Talensi","Lower Anayari","Northern Ghana"))
proc_labels <- c(silicate="Silicate", carbonate="Carbonate", evaporite="Evapoconc.",
                 ion_exchange="Cation exch.", anthropogenic="Nitrate", redox="Redox",
                 trace_mineral="Fluoride", none="Unresolved")
AQ <- c("Middle Voltaian sedimentary aquifer"="Voltaian sed.",
        "Birimian fractured basement aquifer"="Birimian bsmt.",
        "Granitoid fractured basement aquifer"="Granitoid bsmt.",
        "Regolith/alluvial shallow aquifer"="Regolith/alluv.")
AQCOL <- c("Voltaian sed."="#3B6EA8","Birimian bsmt."="#008C7A",
           "Granitoid bsmt."="#C77C2B","Regolith/alluv."="#8C6BB1")
aqf <- function(x) factor(unname(AQ[x]), levels = unname(AQ))
save_s <- function(p, name, w, h) m6_save(p, file.path(OUT, name), w, h)

## S1 — Hydrochemical context (2x2: Gibbs x2, cation-anion, TDS by dataset) ----
s1 <- function() {
  h <- rd("m6_hydrochem_context.csv") |> mutate(ds = dslab(dataset))
  base <- function(g) g + scale_colour_manual(values = dscol, name = NULL) + theme_m6()
  a <- base(ggplot(h, aes(na_ratio, tds_mgL, colour = ds)) +
    geom_point(size = 1.5, alpha = 0.6) + scale_y_log10() +
    labs(subtitle = "Gibbs: cation ratio", x = "Na / (Na + Ca)", y = "TDS (mg/L)"))
  b <- base(ggplot(h, aes(cl_ratio, tds_mgL, colour = ds)) +
    geom_point(size = 1.5, alpha = 0.6) + scale_y_log10() +
    labs(subtitle = "Gibbs: anion ratio", x = "Cl / (Cl + HCO3)", y = "TDS (mg/L)"))
  c <- base(ggplot(h, aes(nak_pct, hco3_pct, colour = ds)) +
    geom_point(size = 1.5, alpha = 0.6) +
    labs(subtitle = "Major-ion composition", x = "Na+K (% cations)",
         y = "HCO3 (% anions)"))
  d <- ggplot(h, aes(ds, tds_mgL, fill = ds)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.32) +
    scale_fill_manual(values = dscol, guide = "none") + scale_y_log10() +
    labs(subtitle = "Mineralisation by dataset", x = NULL, y = "TDS (mg/L)") +
    theme_m6() + theme(axis.text.x = element_text(angle = 15, hjust = 1))
  save_s(m6_tag((a | b) / (c | d)), "figureS1_hydrochem_context", 8.8, 8.0)
}

## S2 — Charge balance (distribution + QC-class composition) ------------------
s2 <- function() {
  cbe <- rd("m6_cbe_distribution.csv") |> mutate(ds = dslab(dataset))
  a <- ggplot(cbe, aes(ds, cbe_pct, fill = ds)) +
    geom_hline(yintercept = c(-5, 5), linetype = "dashed", colour = "#999999") +
    geom_hline(yintercept = c(-10, 10), linetype = "dotted", colour = "#BBBBBB") +
    geom_violin(alpha = 0.5, linewidth = 0.3) +
    geom_boxplot(width = 0.15, outlier.size = 0.4, linewidth = 0.3) +
    scale_fill_manual(values = dscol, guide = "none") +
    coord_cartesian(ylim = c(-40, 20)) +
    labs(subtitle = "Charge-balance error (dashed ±5%, dotted ±10%)",
         x = NULL, y = "CBE (%)") +
    theme_m6() + theme(axis.text.x = element_text(angle = 12, hjust = 1))
  cls <- cbe |> mutate(cbe_class = factor(cbe_class,
              levels = c("quantitative","screening","exploratory","unusable"))) |>
    count(ds, cbe_class)
  b <- ggplot(cls, aes(ds, n, fill = cbe_class)) +
    geom_col(position = "fill", width = 0.7) +
    scale_fill_manual(values = c(quantitative="#2166AC", screening="#67A9CF",
                       exploratory="#EF8A62", unusable="#B2182B"), name = NULL) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Quality-class composition", x = NULL, y = "Share of samples") +
    theme_m6() + theme(axis.text.x = element_text(angle = 12, hjust = 1))
  save_s(m6_tag(a | b), "figureS2_charge_balance", 9.2, 4.6)
}

## S3 — Seasonal change (evapoconcentration + reactive signal by aquifer) -----
s3 <- function() {
  p <- rd("m6_ng_field_pairs.csv") |> mutate(aq = aqf(aquifer))
  a <- ggplot(p, aes(aq, evapo_factor, fill = aq)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "#999999") +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.32) +
    scale_fill_manual(values = AQCOL, guide = "none") +
    labs(subtitle = "Dry/wet evapoconcentration (Cl-based)", x = NULL,
         y = "Cl_dry / Cl_wet") +
    theme_m6() + theme(axis.text.x = element_text(angle = 18, hjust = 1))
  b <- ggplot(p, aes(aq, reactive_rmse, fill = aq)) +
    geom_boxplot(width = 0.6, outlier.size = 0.4, linewidth = 0.32) +
    scale_fill_manual(values = AQCOL, guide = "none") +
    labs(subtitle = "Reactive residual after transport correction", x = NULL,
         y = "Reactive RMSE (mmol/L)") +
    theme_m6() + theme(axis.text.x = element_text(angle = 18, hjust = 1))
  save_s(m6_tag(a | b, collect = FALSE), "figureS3_seasonal_change", 9.2, 4.6)
}

## S4 — Alternative edge networks (composition + divergence metrics) ----------
s4 <- function() {
  f <- rd("m6_edge_family_distribution.csv")
  es <- c(provided_graph="Provided graph", chemistry_knn="Chemistry kNN",
          geographic_nearest="Geographic", random_perturbed="Random")
  f <- f |> mutate(edge_set = factor(unname(es[edge_set]), levels = unname(es)))
  a <- ggplot(f, aes(edge_set, frac, fill = dominant_family)) +
    geom_col(width = 0.72) + m6_fill_scale(name = NULL, labels = proc_labels,
                                           guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Inferred process-network composition", x = NULL,
         y = "Share of edges") +
    theme_m6() + theme(axis.text.x = element_text(angle = 15, hjust = 1))
  s <- rd("m6_edge_network_summary.csv") |>
    mutate(edge_set = factor(unname(es[edge_set]), levels = unname(es))) |>
    select(edge_set, `TVD vs provided` = tvd_vs_provided,
           `Mean stability` = mean_stability) |>
    pivot_longer(-edge_set)
  b <- ggplot(s, aes(edge_set, value, fill = name)) +
    geom_col(position = position_dodge(0.7), width = 0.65) +
    scale_fill_manual(values = c("TVD vs provided"="#C77C2B","Mean stability"="#2166AC"),
                      name = NULL) +
    labs(subtitle = "Divergence and stability vs provided graph", x = NULL, y = NULL) +
    theme_m6() + theme(axis.text.x = element_text(angle = 15, hjust = 1))
  save_s(m6_tag(a | b), "figureS4_edge_networks", 9.6, 4.8)
}

## S5 — Tier x process composition heatmap (single) ---------------------------
s5 <- function() {
  abl <- rd("m6_tier_ablation.csv")
  tlab <- c(tier0_majors="T0", tier1_isotopes="T1", tier2_fluoride="T2",
            tier3_sr_sio2="T3", tier4_full_metadata="T4")
  d <- abl |> mutate(tier = factor(tlab[tier], levels = unname(tlab)),
                     family = factor(unname(proc_labels[dominant_family]),
                                     levels = unname(proc_labels))) |>
    count(tier, family) |> group_by(tier) |> mutate(frac = n / sum(n)) |> ungroup()
  p <- ggplot(d, aes(tier, fct_rev(family), fill = frac)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = n), size = 3.4) +
    scale_fill_gradient(low = "#F7F7F7", high = "#08519C", labels = percent,
                        name = "Share of wells",
                        guide = guide_colourbar(barwidth = 10, barheight = 0.5,
                                                title.position = "top")) +
    labs(title = "Figure S5. Dominant-process composition across evidence tiers",
         x = NULL, y = NULL) +
    theme_m6() + theme(panel.grid = element_blank(), axis.line = element_blank(),
                       axis.ticks = element_blank())
  save_s(p, "figureS5_tier_family_heatmap", 8.0, 5.2)
}

## S6 — Diagnostic distributions (support stability + MRS histograms) ---------
s6 <- function() {
  p <- rd("m6_ng_field_pairs.csv")
  a <- ggplot(p, aes(support_stability)) +
    geom_histogram(bins = 24, fill = "#2166AC", colour = "white", linewidth = 0.2) +
    labs(subtitle = "Bootstrap support stability", x = "Mean Jaccard (5% noise)",
         y = "Wells") + theme_m6()
  b <- ggplot(p, aes(mrs)) +
    geom_histogram(bins = 24, fill = "#008C7A", colour = "white", linewidth = 0.2) +
    labs(subtitle = "Mechanism Resolution Score", x = "MRS", y = "Wells") + theme_m6()
  save_s(m6_tag(a | b, collect = FALSE), "figureS6_diagnostic_distributions", 9.2, 4.4)
}

## S7 — Reactive residual vs evapoconcentration (single, coloured by process) -
s7 <- function() {
  p <- rd("m6_ng_field_pairs.csv")
  g <- ggplot(p, aes(evapo_factor, reactive_rmse, colour = dominant_family)) +
    geom_point(size = 1.8, alpha = 0.75) +
    m6_colour_scale(name = NULL, labels = proc_labels, guide = guide_legend(nrow = 1)) +
    labs(title = "Figure S7. Reactive residual versus conservative evapoconcentration",
         x = "Evapoconcentration factor (Cl_dry / Cl_wet)",
         y = "Reactive RMSE (mmol/L)") + theme_m6()
  save_s(g, "figureS7_reactive_vs_evapo", 8.4, 5.4)
}

## S8 — MRS vs support stability diagnostic (single) --------------------------
s8 <- function() {
  p <- rd("m6_ng_field_pairs.csv")
  g <- ggplot(p, aes(support_stability, mrs, colour = dominant_family)) +
    geom_point(size = 1.8, alpha = 0.75) +
    m6_colour_scale(name = NULL, labels = proc_labels, guide = guide_legend(nrow = 1)) +
    labs(title = "Figure S8. Mechanism Resolution Score versus support stability",
         x = "Support stability (mean Jaccard)", y = "MRS") + theme_m6()
  save_s(g, "figureS8_mrs_diagnostic", 8.4, 5.4)
}

## S9 — External transfer detail (faceted by dataset) ------------------------
s9 <- function() {
  e <- rd("m6_external_transfer.csv") |>
    mutate(ds = recode(dataset, manu = "Lower Anayari (Tier 2)",
                       talensi = "Talensi (Tier 1)"))
  g <- ggplot(e, aes(support_stability, mrs, colour = dominant_family)) +
    geom_point(size = 2, alpha = 0.8) + facet_wrap(~ds) +
    m6_colour_scale(name = NULL, labels = proc_labels, guide = guide_legend(nrow = 1)) +
    labs(title = "Figure S9. External-dataset transfer: resolution versus stability",
         x = "Support stability (mean Jaccard)", y = "MRS") + theme_m6()
  save_s(g, "figureS9_external_detail", 9.0, 5.0)
}

## S10 — Identifiability composition per process family and dataset ----------
s10 <- function() {
  lim <- rd("m6_limitation_map.csv") |>
    mutate(resolution_class = factor(resolution_class, levels = m6_resolution_levels),
           family = factor(unname(proc_labels[dominant_family]),
                           levels = unname(proc_labels)))
  g <- ggplot(lim, aes(family, frac, fill = resolution_class)) +
    geom_col(width = 0.72) + facet_wrap(~dslab(dataset)) +
    m6_fill_scale(labels = m6_resolution_labels, name = NULL,
                  guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(title = "Figure S10. Identifiability composition per process family and dataset",
         x = NULL, y = "Share") +
    theme_m6() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_s(g, "figureS10_limitation_detail", 10.5, 5.2)
}

for (f in list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10)) {
  ok <- try(f(), silent = FALSE)
  if (inherits(ok, "try-error")) cat("  [warn] a supplementary figure failed\n")
  else cat("supp figure done\n")
}
cat("Supplementary figures written to", OUT, "\n")
