## M3 network-dating ambiguity demonstration — Nature-style figure.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr); library(patchwork)
})
HERE <- tryCatch(dirname(normalizePath(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(HERE) == 0 || is.na(HERE)) HERE <- getwd()
THEME <- normalizePath(file.path(HERE, "..", "..", "..", "M6",
         "m6_field_transfer_benchmark", "r_figures", "theme_m6.R"))
source(THEME)
BENCH <- normalizePath(file.path(HERE, ".."))
RES <- file.path(BENCH, "results")
OUT <- file.path(BENCH, "figures", "Manuscript_Ready")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) suppressMessages(read_csv(file.path(RES, f), show_col_types = FALSE))

demo <- rd("m3_network_dating_demo.csv")
summ <- rd("m3_network_dating_demo_summary.csv")

## a: age recovery — single (aliased) vs network (ordering-resolved), ambiguous nodes
al <- demo |> filter(ambiguous == 1) |>
  select(true_age, age_single, age_network) |>
  pivot_longer(-true_age, names_to = "method", values_to = "est") |>
  mutate(method = recode(method, age_single = "Single-node",
                         age_network = "Network (flow-ordered)"))
a <- ggplot(al, aes(true_age, est, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#999999") +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_colour_manual(values = c("Single-node" = "#B2182B",
                      "Network (flow-ordered)" = "#2166AC"), name = NULL) +
  labs(subtitle = "Ambiguous-zone age recovery",
       x = "True age (yr)", y = "Estimated age (yr)") + theme_m6()

## b: within-factor-2 by subset, single vs network
bl <- summ |> select(subset, Single = wf2_single, Network = wf2_network) |>
  pivot_longer(-subset, names_to = "method", values_to = "wf2") |>
  mutate(subset = factor(subset, levels = c("unambiguous", "ambiguous", "all"),
                         labels = c("Unambiguous", "Ambiguous", "All")))
b <- ggplot(bl, aes(subset, wf2, fill = method)) +
  geom_col(position = position_dodge(0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * wf2)),
            position = position_dodge(0.72), vjust = -0.4, size = 3.0) +
  scale_fill_manual(values = c("Single" = "#B2182B", "Network" = "#2166AC"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08),
                     expand = expansion(c(0, 0.02))) +
  labs(subtitle = "Accuracy gain is confined to the ambiguous regime",
       x = NULL, y = "Within factor 2 of true age") + theme_m6()

p <- m6_tag(a | b,
  caption = "Network-enhanced dating adds no accuracy where single-node already works, but resolves tritium bomb-peak ambiguity via flow-ordering.")
m6_save(p, file.path(OUT, "Suppl_Fig4_Network_Dating_Ambiguity"), 9.2, 4.8)
cat("M3 network-dating figure done\n")
