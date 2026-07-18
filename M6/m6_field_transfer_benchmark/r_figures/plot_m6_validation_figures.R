## M6 validation & reviewer-response figures (Fig 7 synthetic validation,
## Fig 8 circularity + transport sensitivity). Nature style, 2-column.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(forcats); library(stringr); library(scales)
})
HERE <- tryCatch(dirname(normalizePath(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(HERE) == 0 || is.na(HERE)) HERE <- getwd()
source(file.path(HERE, "theme_m6.R"))
BENCH <- normalizePath(file.path(HERE, ".."))
RES <- file.path(BENCH, "results"); OUT <- file.path(BENCH, "figures", "r_publication")
rd <- function(f) suppressMessages(read_csv(file.path(RES, f), show_col_types = FALSE))
proc_labels <- c(silicate="Silicate", carbonate="Carbonate", evaporite="Evapoconc.",
                 ion_exchange="Cation exch.", anthropogenic="Nitrate", redox="Redox",
                 trace_mineral="Fluoride", none="Unresolved")
PAN <- c(T0_majors="T0\nmajors", T2_plus_F="T2\n+F", T3_plus_Sr_SiO2="T3\n+Sr/SiO2")
panf <- function(x) factor(PAN[x], levels = unname(PAN))

## ---- FIGURE 7: synthetic-field validation (ground truth) -------------------
fig7 <- function() {
  rec <- rd("m6_synthetic_recovery_by_tier.csv")
  val <- rd("m6_synthetic_validation.csv")
  st  <- rd("m6_synthetic_structural.csv")
  dist <- rd("m6_synthetic_class_distribution.csv")

  arch <- rec |> mutate(archetype = str_to_title(archetype), panel = panf(panel))
  a <- ggplot(arch, aes(archetype, dominant_accuracy, fill = panel)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    scale_fill_manual(values = c("T0\nmajors"="#EF8A62","T2\n+F"="#FDDBC7",
                       "T3\n+Sr/SiO2"="#2166AC"), name = NULL) +
    scale_y_continuous(labels = percent, limits = c(0, 1.02),
                       expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Dominant-process accuracy vs known truth", x = NULL,
         y = "Accuracy") +
    theme_m6() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

  lev <- val |> group_by(panel) |>
    summarise(`Exact mineral` = mean(phase_f1), `Process family` = mean(family_f1),
              .groups = "drop") |>
    pivot_longer(-panel) |> mutate(panel = panf(panel))
  b <- ggplot(lev, aes(panel, value, colour = name, group = name)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
    scale_colour_manual(values = c("Exact mineral"="#B2182B","Process family"="#2166AC"),
                        name = NULL) +
    scale_y_continuous(labels = percent, limits = c(0.4, 0.85)) +
    labs(subtitle = "Recovery: family resolvable, exact mineral not", x = NULL,
         y = "Mean F1 vs truth") + theme_m6()

  st2 <- st |> mutate(panel = panf(panel))
  sc <- max(st2$rank)
  c <- ggplot(st2, aes(panel)) +
    geom_col(aes(y = rank, fill = "Matrix rank"), width = 0.55) +
    geom_line(aes(y = silicate_coherence * sc, colour = "Silicate coherence", group = 1),
              linewidth = 0.9) +
    geom_point(aes(y = silicate_coherence * sc, colour = "Silicate coherence"), size = 2.4) +
    scale_fill_manual(values = c("Matrix rank"="#BFD3E6"), name = NULL) +
    scale_colour_manual(values = c("Silicate coherence"="#B2182B"), name = NULL) +
    scale_y_continuous(name = "Matrix rank",
                       sec.axis = sec_axis(~./sc, name = "Silicate coherence")) +
    labs(subtitle = "Structural identifiability (gate-independent)", x = NULL) +
    theme_m6()

  dd <- dist |> mutate(panel = panf(panel),
        true_class = factor(true_class, levels = m6_resolution_levels))
  d <- ggplot(dd, aes(panel, n, fill = true_class)) +
    geom_col(position = "fill", width = 0.7) +
    m6_fill_scale(labels = m6_resolution_labels, name = NULL,
                  guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "True identifiability class spans full range", x = NULL,
         y = "Share of wells") + theme_m6()

  p <- m6_tag((a | b) / (c | d),
    caption = "Sr/SiO2 improves process-family recovery and structural rank against known truth; exact-mineral recovery stays limited at every tier.")
  m6_save(p, file.path(OUT, "figure7_synthetic_validation"), 8.8, 8.2)
}

## ---- FIGURE 8: circularity resolution + transport sensitivity --------------
fig8 <- function() {
  gs <- rd("m6_field_gate_structural.csv")
  ts <- rd("m6_transport_sensitivity.csv")
  tl <- c(tier0_majors="T0\nmajors", tier1_isotopes="T1\n+iso", tier2_fluoride="T2\n+F",
          tier3_sr_sio2="T3\n+Sr/SiO2", tier4_full_metadata="T4\nfull")
  gs <- gs |> mutate(tier = factor(tl[tier], levels = unname(tl)),
                     xi = as.integer(tier))
  a <- ggplot(gs, aes(xi)) +
    geom_ribbon(aes(ymin = gated_ci_lo, ymax = gated_ci_hi, fill = "Gated 95% CI"),
                alpha = 0.2) +
    geom_line(aes(y = gated_frac_non_identifiable, colour = "Gated (with prior)"),
              linewidth = 0.9) +
    geom_point(aes(y = gated_frac_non_identifiable, colour = "Gated (with prior)"), size = 2.3) +
    geom_line(aes(y = ungated_frac_non_identifiable, colour = "Ungated (classifier only)"),
              linewidth = 0.9, linetype = "22") +
    geom_point(aes(y = ungated_frac_non_identifiable, colour = "Ungated (classifier only)"),
               size = 2.3) +
    scale_colour_manual(values = c("Gated (with prior)"="#B2182B",
                        "Ungated (classifier only)"="#2166AC"), name = NULL) +
    scale_fill_manual(values = c("Gated 95% CI"="#B2182B"), name = NULL) +
    scale_x_continuous(breaks = 1:5, labels = unname(tl)) +
    scale_y_continuous(labels = percent, limits = c(0, 0.75)) +
    labs(subtitle = "Field non-identifiable share: gate vs classifier",
         x = NULL, y = "Fraction non-identifiable") + theme_m6()

  ts2 <- ts |> mutate(treatment = factor(treatment, levels = c("none","cl","recharge"),
                      labels = c("No correction","Cl-conservative","Recharge endmember")))
  b <- ggplot(ts2, aes(treatment, frac, fill = dominant_family)) +
    geom_col(width = 0.68) +
    m6_fill_scale(name = NULL, labels = proc_labels, guide = guide_legend(nrow = 1)) +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Process attribution vs transport treatment", x = NULL,
         y = "Share of wells") +
    theme_m6() + theme(axis.text.x = element_text(angle = 12, hjust = 1))

  p <- m6_tag(a | b,
    caption = "The field collapse is the conservative prior (ungated: all stay partially identifiable); attribution is stable across principled corrections.")
  m6_save(p, file.path(OUT, "figure8_circularity_sensitivity"), 9.6, 5.2)
}

ok7 <- try(fig7(), silent = FALSE); cat(if (inherits(ok7,"try-error")) "fig7 FAIL\n" else "fig7 done\n")
ok8 <- try(fig8(), silent = FALSE); cat(if (inherits(ok8,"try-error")) "fig8 pending (data not ready)\n" else "fig8 done\n")
