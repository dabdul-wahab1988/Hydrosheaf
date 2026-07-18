## M7 coupled-integration benchmark — Nature-style keystone figure.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(forcats); library(scales)
})
HERE <- tryCatch(dirname(normalizePath(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) getwd())
if (length(HERE) == 0 || is.na(HERE)) HERE <- getwd()
THEME <- normalizePath(file.path(HERE, "..", "..", "..", "M6",
         "m6_field_transfer_benchmark", "r_figures", "theme_m6.R"))
source(THEME)
BENCH <- normalizePath(file.path(HERE, ".."))
RES <- file.path(BENCH, "results"); OUT <- file.path(BENCH, "figures", "r_publication")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) suppressMessages(read_csv(file.path(RES, f), show_col_types = FALSE))
STREAM <- c(geometry_only="Geometry", chemistry_only="Chemistry", age_only="Age", joint="Joint")
SCOL <- c(Geometry="#C77C2B", Chemistry="#008C7A", Age="#8C6BB1", Joint="#2166AC")

fig <- function() {
  coh <- rd("m7_age_coherence_demo.csv")
  gain <- rd("m7_integration_gain.csv") |>
    mutate(stream = factor(unname(STREAM[stream]), levels = unname(STREAM)))
  trap <- rd("m7_trap_rejection.csv")

  ## a: age <-> topology coherence — upstream vs downstream tritium age.
  ## True edges lie above the 1:1 line (downstream older); trap A (age reversal)
  ## falls below it and is flagged by the age-coherence audit.
  coh <- coh |> mutate(edge = recode(label, "true" = "True edge",
                                     "trapA" = "Trap A (age reversal)"))
  a <- ggplot(coh, aes(upstream_age, downstream_age, colour = edge)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#999999") +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("True edge"="#2166AC",
                        "Trap A (age reversal)"="#B2182B"), name = NULL) +
    labs(subtitle = "Age-topology: downstream older",
         x = "Upstream age (yr)", y = "Downstream age (yr)") + theme_m6()

  ## b: false-connection (trap) acceptance by stream — integration gain
  b <- ggplot(gain, aes(stream, trap_accept_rate, fill = stream)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = percent(trap_accept_rate, accuracy = 1)), vjust = -0.4,
              size = 3.2) +
    scale_fill_manual(values = SCOL, guide = "none") +
    scale_y_continuous(labels = percent, expand = expansion(c(0, 0.08))) +
    labs(subtitle = "False connections accepted",
         x = NULL, y = "Trap edges accepted") + theme_m6()

  ## c: which stream rejects which trap type (complementarity)
  tl <- trap |> select(trap_type, Geometry = geometry_rejects,
        Chemistry = chemistry_rejects, Age = age_rejects, Joint = joint_rejects) |>
    pivot_longer(-trap_type, names_to = "stream", values_to = "reject") |>
    mutate(stream = factor(stream, levels = c("Geometry","Chemistry","Age","Joint")),
           trap_type = recode(trap_type, trapA = "Trap A\n(age reversal)",
                     trapB = "Trap B\n(chem. impossible)", spurious = "Spurious\n(distant)"))
  c <- ggplot(tl, aes(stream, fct_rev(trap_type), fill = reject)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = percent(reject, accuracy = 1)), size = 3.0) +
    scale_fill_gradient(low = "#F7F7F7", high = "#2166AC", limits = c(0, 1),
                        labels = percent, name = "Rejected",
                        guide = guide_colourbar(barwidth = 9, barheight = 0.5,
                                                title.position = "top")) +
    labs(subtitle = "Complementarity by trap type",
         x = NULL, y = NULL) + theme_m6() +
    theme(panel.grid = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank())

  ## d: precision / recall by stream
  pr <- gain |> select(stream, Precision = precision, Recall = recall) |>
    pivot_longer(-stream, names_to = "metric", values_to = "value")
  d <- ggplot(pr, aes(stream, value, fill = metric)) +
    geom_col(position = position_dodge(0.72), width = 0.66) +
    scale_fill_manual(values = c("Precision"="#2166AC","Recall"="#EF8A62"), name = NULL) +
    scale_y_continuous(labels = percent, limits = c(0, 1.02),
                       expand = expansion(c(0, 0.02))) +
    labs(subtitle = "Precision & recall by stream", x = NULL, y = NULL) +
    theme_m6()

  p <- m6_tag((a | b) / (c | d),
    caption = "Controlled coupled twin: integrating age + topology + chemistry rejects false connections that single streams admit.")
  m6_save(p, file.path(OUT, "figure_m7_integration"), 8.8, 8.2)
}
ok <- try(fig(), silent = FALSE)
cat(if (inherits(ok, "try-error")) "M7 figure FAIL\n" else "M7 figure done\n")
