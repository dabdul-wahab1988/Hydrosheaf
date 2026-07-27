suppressPackageStartupMessages({
  library(DBI)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(viridis)
})

source(file.path("r_figures", "theme_m5.R"))

db_path <- file.path("results", "m5_results.duckdb")
out_dir <- file.path("figures", "r_publication")
supp_dir <- file.path(out_dir, "supplementary")

method_order <- c(
  "bounded_ls",
  "lasso",
  "elastic_net",
  "thermo_elastic_net",
  "hydrosheaf_core",
  "hydrosheaf_guarded"
)
method_names <- c(
  bounded_ls = "Bounded least squares",
  lasso = "Lasso",
  elastic_net = "Elastic net",
  thermo_elastic_net = "Thermodynamic elastic net",
  hydrosheaf_core = "Hydrosheaf-Core",
  hydrosheaf_guarded = "Hydrosheaf-Guarded",
  phreeqc_inverse = "PHREEQC inverse"
)
tier_names <- c(
  core = "Core",
  plus_lite = "Plus-lite",
  enhanced = "Enhanced",
  available_plus_lite = "Available Plus-lite"
)
tier_order <- c("Core", "Plus-lite", "Enhanced", "Available Plus-lite")
ion_order <- c("Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4")
panel_names <- c(
  full_11 = "Full 11-ion",
  major_8 = "Major 8-ion",
  major_7 = "Major 7-ion",
  core_5 = "Core 5-ion",
  no_alkalinity = "No alkalinity",
  no_redox_trace = "No redox/trace",
  sparse_5 = "Sparse 5-ion",
  cation_anion = "Cation-anion"
)
dataset_names <- c(
  "NorthernGhana.xlsx" = "Northern Ghana\nworkbook",
  Talensi = "Talensi\nmining area",
  LowerAnayari = "Lower Anayari"
)
dataset_short_names <- c(
  "NorthernGhana.xlsx" = "Northern\nGhana",
  Talensi = "Talensi",
  LowerAnayari = "Lower\nAnayari"
)
archetype_names <- c(
  carbonate = "Carbonate",
  crystalline = "Crystalline",
  evaporitic = "Evaporitic",
  mixed = "Mixed"
)
table_names <- c(
  benchmark_fits = "Sparse inverse fits",
  reaction_recovery = "Reaction recovery",
  heldout_ion_results = "Held-out ion tests",
  hydrosheaf_core_evidence = "Hydrosheaf-Core evidence",
  hydrosheaf_core_evidence_lifted_resolution = "Hydrosheaf-Core ELRI",
  phreeqc_inverse_baseline_models = "PHREEQC inverse models",
  data_tier_reaction_evidence = "Tier reaction evidence",
  data_tier_evidence_lifted_resolution = "Tier ELRI",
  external_field_evidence_lifted_resolution = "External field ELRI",
  external_field_reaction_evidence = "External field evidence",
  regularization_paths = "Regularisation paths",
  bootstrap_support_stability = "Bootstrap support"
)
setup_names <- c(
  strict_5pct = "Strict 5%",
  relaxed_20pct = "Relaxed 20%"
)
reaction_names <- c(
  CaNa_exch = "Ca-Na exchange",
  MgNa_exch = "Mg-Na exchange",
  NO3src = "Nitrate source",
  NaCa_exch = "Na-Ca exchange",
  NaMg_exch = "Na-Mg exchange",
  albite = "Albite",
  anorthite = "Anorthite",
  apatite = "Apatite",
  calcite = "Calcite",
  denit = "Denitrification",
  dolomite = "Dolomite",
  fluorite = "Fluorite",
  gypsum = "Gypsum",
  halite = "Halite",
  k_feldspar = "K-feldspar",
  pyrite_oxidation = "Pyrite oxidation"
)
# The canonical raw field workbook carries Region/District but no independent
# aquifer-type, geology-group, or lithology classification for these wells
# (see NORTHERN_GHANA_WORKBOOK in run_m5_inverse_reaction_benchmark.py), so
# panels are faceted by Region rather than aquifer type.
region_names <- c(
  "Northern Region" = "Northern",
  "North East Region" = "North East",
  "Upper East Region" = "Upper East",
  "Upper West Region" = "Upper West"
)

recode_named <- function(x, map, levels = NULL) {
  y <- unname(map[as.character(x)])
  missing <- is.na(y)
  if (any(missing)) {
    y[missing] <- as.character(x)[missing]
  }
  y <- as.character(y)
  if (is.null(levels)) {
    y
  } else {
    factor(y, levels = levels)
  }
}

clean_method <- function(x, include_phreeqc = FALSE) {
  levels <- unname(method_names[method_order])
  if (include_phreeqc) {
    levels <- c(levels, method_names[["phreeqc_inverse"]])
  }
  recode_named(x, method_names, levels)
}

clean_tier <- function(x) {
  recode_named(x, tier_names, tier_order)
}

clean_setup <- function(x) {
  recode_named(x, setup_names, unname(setup_names))
}

clean_panel <- function(x) {
  y <- unname(panel_names[as.character(x)])
  missing <- is.na(y)
  if (any(missing)) {
    y[missing] <- gsub("_", " ", as.character(x)[missing])
  }
  as.character(y)
}

clean_dataset <- function(x) {
  y <- unname(dataset_names[as.character(x)])
  missing <- is.na(y)
  if (any(missing)) {
    y[missing] <- as.character(x)[missing]
  }
  as.character(y)
}

clean_dataset_short <- function(x) {
  y <- unname(dataset_short_names[as.character(x)])
  missing <- is.na(y)
  if (any(missing)) {
    y[missing] <- clean_dataset(x)[missing]
  }
  as.character(y)
}

clean_archetype <- function(x) {
  y <- unname(archetype_names[as.character(x)])
  missing <- is.na(y)
  if (any(missing)) {
    y[missing] <- tools::toTitleCase(gsub("_", " ", as.character(x)[missing]))
  }
  as.character(y)
}

clean_table <- function(x) {
  y <- unname(table_names[as.character(x)])
  missing <- is.na(y)
  if (any(missing)) {
    y[missing] <- tools::toTitleCase(gsub("_", " ", as.character(x)[missing]))
  }
  as.character(y)
}

clean_reaction <- function(x) {
  recode_named(x, reaction_names)
}

clean_region <- function(x) {
  recode_named(x, region_names)
}

clean_class_label <- function(x) {
  vapply(
    as.character(x),
    function(value) {
      normalised <- gsub("\\s*/\\s*", ";", value)
      parts <- trimws(unlist(strsplit(normalised, ";", fixed = TRUE)))
      paste(clean_reaction(parts), collapse = " / ")
    },
    character(1)
  )
}

source_map <- c(
  analysis_summary = "results/analysis_summary.json",
  benchmark_fits = "results/benchmark_fits.csv",
  bootstrap_support_stability = "results/bootstrap_support_stability.csv",
  data_tier_experiment = "results/data_tier_experiment.csv",
  data_tier_evidence_lifted_resolution = "results/data_tier_evidence_lifted_resolution.csv",
  data_tier_reaction_evidence = "results/data_tier_reaction_evidence.csv",
  equivalence_classes = "results/equivalence_classes.csv",
  external_field_evidence_lifted_resolution = "results/external_field_evidence_lifted_resolution.csv",
  ghana_field_pairs = "results/ghana_field_pairs.csv",
  ghana_field_hydrosheaf_core_evidence = "results/ghana_field_hydrosheaf_core_evidence.csv",
  heldout_ion_results = "results/heldout_ion_results.csv",
  hydrosheaf_core_evidence = "results/hydrosheaf_core_evidence.csv",
  identifiability_diagnostics = "results/identifiability_diagnostics.csv",
  m5_table_catalog = "results/m5_results_database_catalog.csv",
  next_best_measurement = "results/next_best_measurement.csv",
  phreeqc_inverse_baseline = "results/phreeqc_inverse_baseline.csv",
  reaction_recovery = "results/reaction_recovery.csv",
  regularization_paths = "results/regularization_paths.csv",
  tableS1_reaction_stoichiometry = "tables/tableS1_reaction_stoichiometry.csv",
  tableS12_hydrosheaf_core_evidence_gates = "tables/tableS12_hydrosheaf_core_evidence_gates.csv",
  tableS15_data_tier_reaction_evidence = "tables/tableS15_data_tier_reaction_evidence.csv",
  tableS16_evidence_lifted_resolution = "tables/tableS16_evidence_lifted_resolution.csv",
  tableS17_external_field_evidence_lifted_resolution = "tables/tableS17_external_field_evidence_lifted_resolution.csv",
  thermodynamic_threshold_sensitivity = "results/thermodynamic_threshold_sensitivity.csv"
)

jsonlite_available <- requireNamespace("jsonlite", quietly = TRUE)
duckdb_available <- requireNamespace("duckdb", quietly = TRUE)

con <- NULL
if (duckdb_available && file.exists(db_path)) {
  con <- tryCatch(
    DBI::dbConnect(duckdb::duckdb(), db_path, read_only = TRUE),
    error = function(err) {
      message("Could not open DuckDB database from R: ", conditionMessage(err))
      NULL
    }
  )
}

if (!is.null(con)) {
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tbl_db <- function(name) DBI::dbReadTable(con, name)
  message("Reading M5 plotting data from ", db_path)
} else {
  message("R duckdb package is unavailable; reading generated CSV/JSON mirrors.")
  tbl_db <- function(name) {
    path <- source_map[[name]]
    if (is.null(path) || !file.exists(path)) {
      stop("No source mapping for table: ", name)
    }
    if (grepl("\\.json$", path)) {
      if (!jsonlite_available) {
        stop("jsonlite is required for JSON fallback source: ", path)
      }
      data <- jsonlite::fromJSON(path, simplifyVector = FALSE)
      scalars <- data[vapply(data, function(value) {
        !is.list(value) && length(value) == 1
      }, logical(1))]
      return(as.data.frame(scalars, stringsAsFactors = FALSE))
    }
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

safe_div <- function(num, den) ifelse(den > 0, num / den, NA_real_)

as_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.numeric(x)) {
    return(x != 0)
  }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

summarise_recovery <- function(recovery, method_value = "hydrosheaf_core") {
  recovery |>
    mutate(
      true_active_flag = as_bool(true_active),
      recovered_active_flag = as_bool(recovered_active)
    ) |>
    filter(method == method_value, panel == "full_11", noise_level == 0.03) |>
    group_by(reaction, family) |>
    summarise(
      tp = sum(true_active_flag & recovered_active_flag, na.rm = TRUE),
      fp = sum(!true_active_flag & recovered_active_flag, na.rm = TRUE),
      fn = sum(true_active_flag & !recovered_active_flag, na.rm = TRUE),
      precision = safe_div(tp, tp + fp),
      recall = safe_div(tp, tp + fn),
      f1 = safe_div(2 * precision * recall, precision + recall),
      .groups = "drop"
    ) |>
    arrange(desc(f1))
}

stoich_long <- function(stoich, reactions = NULL) {
  ion_cols <- intersect(ion_order, names(stoich))
  if (!is.null(reactions)) {
    stoich <- stoich[stoich$Reaction %in% reactions, , drop = FALSE]
  }
  do.call(
    rbind,
    lapply(
      ion_cols,
      function(ion) {
        data.frame(
          Reaction = clean_reaction(stoich$Reaction),
          Family = stoich$Family,
          ion = ion,
          coefficient = stoich[[ion]],
          stringsAsFactors = FALSE
        )
      }
    )
  ) |>
    mutate(ion = factor(ion, levels = ion_order))
}

heatmap_scale <- scale_fill_gradient2(
  low = "#3B6EA8", mid = "white", high = "#C77C2B", midpoint = 0,
  name = "Coefficient"
)

figure1_database_design <- function() {
  catalog <- tbl_db("m5_table_catalog")
  summary <- tbl_db("analysis_summary")
  fits <- tbl_db("benchmark_fits")

  n_value <- function(frame, col) {
    if (col %in% names(frame)) as.numeric(frame[[col]][1]) else NA_real_
  }
  evidence_counts <- data.frame(
    metric = factor(
      c("PHREEQC truth cases", "Sparse inverse fits", "Field pairs", "Stored result rows"),
      levels = rev(c("PHREEQC truth cases", "Sparse inverse fits", "Field pairs", "Stored result rows"))
    ),
    value = c(
      n_value(summary, "n_phreeqc_scenarios"),
      n_value(summary, "n_benchmark_fits"),
      n_value(summary, "n_field_pairs"),
      sum(catalog$row_count, na.rm = TRUE)
    )
  )
  p1a <- ggplot(evidence_counts, aes(value, metric)) +
    geom_col(width = 0.65, fill = "#008C7A") +
    geom_text(aes(label = comma(value)), hjust = -0.08, size = 2.15) +
    scale_x_log10(labels = comma, expand = expansion(mult = c(0, 0.22))) +
    labs(
      title = "Benchmark evidence scale",
      subtitle = "Current run; labels show counts",
      x = "Count (log10 scale)",
      y = NULL
    ) +
    theme_m5()

  top_tables <- catalog |>
    arrange(desc(row_count)) |>
    slice_head(n = 8) |>
    mutate(table_label = factor(clean_table(table_name), levels = rev(clean_table(table_name))))
  p1b <- ggplot(top_tables, aes(row_count, table_label)) +
    geom_col(width = 0.68, fill = "#3B6EA8") +
    scale_x_continuous(labels = comma) +
    labs(title = "Largest evidence tables", subtitle = "Rows imported into DuckDB", x = "Rows", y = NULL) +
    theme_m5()

  panel_perf <- fits |>
    group_by(panel, method) |>
    summarise(class_f1 = mean(class_f1, na.rm = TRUE), .groups = "drop") |>
    mutate(method_label = clean_method(method), panel_label = clean_panel(panel))
  p1c <- ggplot(panel_perf, aes(panel_label, method_label, fill = class_f1)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_viridis(option = "C", limits = c(0, 1), name = "Class F1") +
    labs(
      title = "Recovery depends on measured chemistry",
      subtitle = "Mean equivalence-class F1 by method and ion panel",
      x = "Measured chemistry panel",
      y = NULL
    ) +
    theme_m5_heatmap() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  save_m5(
    (tag_panel(p1a, "a") | tag_panel(p1b, "b")) /
      tag_panel(p1c, "c") +
      plot_layout(heights = c(1, 1.08)),
    file.path(out_dir, "figure1_r_database_design.png"),
    183,
    118
  )
}

figure2_model_performance <- function() {
  fits <- tbl_db("benchmark_fits") |>
    filter(panel == "full_11", noise_level == 0.03) |>
    mutate(method_label = clean_method(method))
  phreeqc <- tbl_db("phreeqc_inverse_baseline")

  p2a <- ggplot(fits, aes(class_f1, method_label, fill = method)) +
    geom_boxplot(width = 0.62, outlier.alpha = 0.18, linewidth = 0.25) +
    stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 1.5, stroke = 0.25) +
    m5_fill_scale(guide = "none") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Class-level recovery",
      subtitle = "Full 11-ion panel, 3% noise; white dot = mean",
      x = "Equivalence-class F1 (higher is better)",
      y = "Inverse method"
    ) +
    theme_m5()

  p2b <- ggplot(fits, aes(false_discovery_rate, method_label, fill = method)) +
    geom_boxplot(width = 0.62, outlier.alpha = 0.18, linewidth = 0.25) +
    stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 1.5, stroke = 0.25) +
    m5_fill_scale(guide = "none") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "False support selections",
      subtitle = "Same scenarios as panel a",
      x = "False-discovery rate (lower is better)",
      y = NULL
    ) +
    theme_m5()

  method_summary <- fits |>
    group_by(method) |>
    summarise(
      class_f1 = mean(class_f1, na.rm = TRUE),
      fdr = mean(false_discovery_rate, na.rm = TRUE),
      extent_rmse = mean(extent_rmse_mmolL, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      method_label = as.character(clean_method(method, include_phreeqc = TRUE)),
      display_method = method
    )
  phreeqc_summary <- phreeqc |>
    summarise(
      class_f1 = mean(first_minimal_class_f1, na.rm = TRUE),
      fdr = mean(first_minimal_class_false_discovery_rate, na.rm = TRUE),
      extent_rmse = NA_real_,
      .groups = "drop"
    ) |>
    mutate(method_label = "PHREEQC inverse", display_method = "phreeqc_inverse")
  pareto <- bind_rows(method_summary, phreeqc_summary)

  p2c <- ggplot(pareto, aes(fdr, class_f1, colour = display_method)) +
    geom_point(size = 2.2) +
    ggrepel::geom_text_repel(
      aes(label = method_label),
      size = 2.05,
      min.segment.length = 0,
      segment.size = 0.2,
      box.padding = 0.16,
      seed = 7,
      show.legend = FALSE
    ) +
    m5_colour_scale(guide = "none") +
    coord_cartesian(xlim = c(0.10, 0.85), ylim = c(0.30, 0.68)) +
    labs(
      title = "Performance frontier",
      subtitle = "PHREEQC = first-minimal inverse baseline",
      x = "Mean false-discovery rate",
      y = "Mean equivalence-class F1"
    ) +
    theme_m5()

  save_m5(
    (tag_panel(p2a, "a") | tag_panel(p2b, "b")) /
      tag_panel(p2c, "c") +
      plot_layout(heights = c(1, 0.86)),
          file.path(out_dir, "figure2_r_model_performance.png"), 183, 125)
}

figure3_equifinality_elri <- function() {
  stoich <- tbl_db("tableS1_reaction_stoichiometry")
  ambiguous <- tbl_db("equivalence_classes") |>
    mutate(ambiguous_flag = tolower(as.character(ambiguous)) %in% c("true", "1")) |>
    filter(ambiguous_flag) |>
    mutate(members_label = factor(clean_class_label(members), levels = rev(clean_class_label(members))))
  elri <- tbl_db("tableS16_evidence_lifted_resolution") |>
    mutate(tier = clean_tier(data_tier), members_label = clean_class_label(members))

  ambiguous_reactions <- unique(unlist(strsplit(paste(ambiguous$members, collapse = ";"), ";")))
  s_long <- stoich_long(stoich, ambiguous_reactions)
  s_long$Reaction <- factor(s_long$Reaction, levels = clean_reaction(ambiguous_reactions))

  p3a <- ggplot(ambiguous, aes(n_members, members_label)) +
    geom_col(width = 0.55, fill = "#737373") +
    scale_x_continuous(breaks = c(0, 1, 2), limits = c(0, 2.2)) +
    labs(
      title = "Exact non-unique classes",
      subtitle = "Sparse major-ion mass balance cannot separate these members",
      x = "Number of indistinguishable reaction members",
      y = NULL
    ) +
    theme_m5()

  p3b <- ggplot(s_long, aes(ion, Reaction, fill = coefficient)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    heatmap_scale +
    labs(title = "Sparse-ion stoichiometry", x = "Measured ion", y = "Reaction member") +
    theme_m5_heatmap() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p3c <- ggplot(elri, aes(mean_elri, members_label, fill = tier)) +
    geom_col(position = position_dodge(width = 0.74), width = 0.62) +
    scale_fill_manual(values = c("Core" = "#737373", "Plus-lite" = "#008C7A", "Enhanced" = "#C77C2B")) +
    labs(
      title = "Conditional evidence separation",
      x = "Mean ELRI (0 = unresolved, 1 = resolved)",
      y = NULL,
      fill = "Measurement tier"
    ) +
    theme_m5()

  calcite <- elri |>
    filter(members_label == "Anorthite / Calcite")
  p3d <- ggplot(calcite, aes(tier, top_member_true_active_fraction_when_class_active, group = 1)) +
    geom_line(linewidth = 0.35, colour = "#008C7A") +
    geom_point(size = 1.8, colour = "#008C7A") +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(
      title = "Calcite/anorthite discrimination",
      subtitle = "Correct only when diagnostic evidence supports a top member",
      x = "Measurement tier",
      y = "Top member correct when class active"
    ) +
    theme_m5() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  save_m5(
    (tag_panel(p3a, "a") | tag_panel(p3b, "b")) /
      (tag_panel(p3c, "c") | tag_panel(p3d, "d")),
          file.path(out_dir, "figure3_r_equifinality_elri.png"), 183, 125)
}

figure4_data_tiers <- function() {
  tiers <- tbl_db("data_tier_experiment") |>
    mutate(tier = clean_tier(data_tier))
  elri <- tbl_db("data_tier_evidence_lifted_resolution") |>
    mutate(tier = clean_tier(data_tier))
  tier_evidence <- tbl_db("tableS15_data_tier_reaction_evidence") |>
    mutate(
      tier = clean_tier(data_tier),
      reaction_label = factor(clean_reaction(reaction), levels = clean_reaction(unique(reaction)))
    )

  p4a <- ggplot(tiers, aes(tier, class_f1, fill = data_tier)) +
    geom_boxplot(width = 0.58, outlier.alpha = 0.15, linewidth = 0.25) +
    stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 1.5, stroke = 0.25) +
    scale_fill_manual(values = c(core = "#737373", plus_lite = "#008C7A", enhanced = "#C77C2B"), guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Reaction-class recovery", x = "Measurement tier", y = "Class F1 (higher is better)") +
    theme_m5()

  p4b <- ggplot(tiers, aes(tier, false_discovery_rate, fill = data_tier)) +
    geom_boxplot(width = 0.58, outlier.alpha = 0.15, linewidth = 0.25) +
    stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 1.5, stroke = 0.25) +
    scale_fill_manual(values = c(core = "#737373", plus_lite = "#008C7A", enhanced = "#C77C2B"), guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Overinterpretation control", x = "Measurement tier", y = "FDR (lower is better)") +
    theme_m5()

  p4c <- ggplot(elri, aes(tier, evidence_lifted_resolution_index, fill = data_tier)) +
    geom_boxplot(width = 0.58, outlier.alpha = 0.15, linewidth = 0.25) +
    scale_fill_manual(values = c(core = "#737373", plus_lite = "#008C7A", enhanced = "#C77C2B"), guide = "none") +
    labs(title = "Within-class ambiguity lifting", x = "Measurement tier", y = "ELRI") +
    theme_m5()

  p4d <- ggplot(tier_evidence, aes(tier, reaction_label, fill = mean_combined_evidence_score)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis(option = "C", limits = c(0, 1), name = "Evidence") +
    labs(title = "Reaction-level evidence gates", x = "Measurement tier", y = "Reaction") +
    theme_m5_heatmap()

  save_m5(
    (tag_panel(p4a, "a") | tag_panel(p4b, "b")) /
      (tag_panel(p4c, "c") | tag_panel(p4d, "d")) +
      plot_layout(widths = c(1, 1.18)),
          file.path(out_dir, "figure4_r_data_tiers.png"), 183, 125)
}

figure5_phreeqc_thermo <- function() {
  thermo <- tbl_db("thermodynamic_threshold_sensitivity")
  phreeqc <- tbl_db("phreeqc_inverse_baseline") |>
    mutate(
      setup = factor(phreeqc_inverse_setup, levels = c("strict_5pct", "relaxed_20pct")),
      setup_label = clean_setup(phreeqc_inverse_setup)
    )

  p5a <- thermo |>
    group_by(si_threshold, archetype) |>
    summarise(class_f1 = mean(class_f1, na.rm = TRUE), fdr = mean(false_discovery_rate, na.rm = TRUE), .groups = "drop") |>
    ggplot(aes(si_threshold, class_f1, colour = archetype)) +
    geom_line(linewidth = 0.35) +
    geom_point(size = 1.3) +
    m5_colour_scale(labels = archetype_names) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Thermodynamic screening sensitivity",
      x = "SI exclusion threshold",
      y = "Mean class F1"
    ) +
    theme_m5()

  p5b <- ggplot(phreeqc, aes(models_found + 1, setup_label, fill = setup)) +
    geom_boxplot(width = 0.58, outlier.alpha = 0.20, linewidth = 0.25) +
    scale_x_log10(labels = comma) +
    scale_fill_manual(values = c(strict_5pct = "#737373", relaxed_20pct = "#C77C2B"), guide = "none") +
    labs(
      title = "PHREEQC inverse-model multiplicity",
      x = "Feasible inverse models + 1 (log10 scale)",
      y = "Uncertainty setting"
    ) +
    theme_m5()

  p5c <- ggplot(phreeqc, aes(models_found + 1, first_minimal_class_f1, colour = archetype)) +
    geom_point(alpha = 0.58, size = 1.15) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(limits = c(0, 1)) +
    m5_colour_scale(labels = archetype_names) +
    labs(
      title = "First-minimal model uncertainty",
      x = "Feasible inverse models + 1",
      y = "First-minimal class F1"
    ) +
    theme_m5()

  p5d_data <- phreeqc |>
    group_by(archetype, setup) |>
    summarise(
      success = mean(as_bool(phreeqc_inverse_success), na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(is.finite(success)) |>
    mutate(archetype_label = clean_archetype(archetype))
  p5d <- ggplot(p5d_data, aes(archetype_label, success, fill = setup)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.60) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c(strict_5pct = "#737373", relaxed_20pct = "#C77C2B"),
                      labels = c("Strict 5%", "Relaxed 20%")) +
    labs(title = "Feasibility requires fallback", x = "Scenario family", y = "Successful inverse runs") +
    theme_m5() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  save_m5(
    (tag_panel(p5a, "a") | tag_panel(p5b, "b")) /
      (tag_panel(p5c, "c") | tag_panel(p5d, "d")),
          file.path(out_dir, "figure5_r_phreeqc_thermo.png"), 183, 125)
}

figure6_field_transfer <- function() {
  ghana <- tbl_db("ghana_field_pairs") |>
    mutate(region_label = clean_region(region))
  external <- tbl_db("external_field_evidence_lifted_resolution") |>
    mutate(
      tier = clean_tier(data_tier),
      members_label = clean_class_label(members),
      dataset_label = clean_dataset(dataset)
    )
  external_summary <- tbl_db("tableS17_external_field_evidence_lifted_resolution") |>
    mutate(
      tier = clean_tier(data_tier),
      members_label = clean_class_label(members),
      dataset_label = clean_dataset(dataset),
      dataset_short_label = clean_dataset_short(dataset)
    )

  p6a <- ghana |>
    count(resolution_class) |>
    mutate(frac = n / sum(n), resolution_label = gsub("_", "\n", resolution_class)) |>
    ggplot(aes(resolution_label, frac)) +
    geom_col(width = 0.62, fill = "#3B6EA8") +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(title = "Northern Ghana transfer audit", x = "Resolution class", y = "Fraction of wet-dry pairs") +
    theme_m5()

  p6b <- ggplot(external, aes(dataset_label, evidence_lifted_resolution_index, fill = tier)) +
    geom_boxplot(width = 0.62, outlier.alpha = 0.15, linewidth = 0.25) +
    scale_fill_manual(values = c("Core" = "#737373", "Available Plus-lite" = "#008C7A")) +
    labs(title = "External sparse-field transfer", x = "Field dataset", y = "ELRI") +
    theme_m5() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p6c <- ggplot(
    external_summary,
    aes(dataset_short_label, members_label, fill = median_elri)
  ) +
    geom_tile(colour = "white", linewidth = 0.22) +
    scale_fill_viridis(option = "C", limits = c(0, 0.5), name = "Median ELRI") +
    facet_wrap(~tier, nrow = 1) +
    labs(title = "Class-specific field ambiguity", x = "Field dataset", y = "Ambiguous class") +
    theme_m5_heatmap() +
    theme(axis.text.x = element_text(size = 6.4))

  p6d <- ggplot(
    ghana,
    aes(tds_delta_consistency_score, mean_hydrosheaf_core_evidence_score, colour = region_label)
  ) +
    geom_point(alpha = 0.72, size = 1.25) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(
      values = c(
        "Northern" = "#40BFC1",
        "North East" = "#F8766D",
        "Upper East" = "#7CAE00",
        "Upper West" = "#C77CFF"
      )
    ) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 1.9))) +
    labs(
      title = "Plausibility, not reaction truth",
      x = "TDS consistency score",
      y = "Mean core evidence score",
      colour = "Region"
    ) +
    theme_m5() +
    theme(
      legend.position = c(0.03, 0.05),
      legend.justification = c(0, 0),
      legend.background = element_rect(fill = scales::alpha("white", 0.82), colour = NA),
      legend.key.height = grid::unit(3.4, "mm"),
      legend.key.width = grid::unit(3.4, "mm"),
      legend.text = element_text(size = 5.8),
      legend.title = element_text(size = 6.1)
    )

  save_m5(
    (tag_panel(p6a, "a") | tag_panel(p6b, "b")) /
      (tag_panel(p6c, "c") | tag_panel(p6d, "d")) +
      plot_layout(heights = c(1, 1.12)),
          file.path(out_dir, "figure6_r_field_transfer.png"), 183, 135)
}

supplementary_figures <- function() {
  stoich <- tbl_db("tableS1_reaction_stoichiometry")
  diag <- tbl_db("identifiability_diagnostics")
  fits <- tbl_db("benchmark_fits")
  recovery <- tbl_db("reaction_recovery")
  nb <- tbl_db("next_best_measurement")
  paths <- tbl_db("regularization_paths")
  bootstrap <- tbl_db("bootstrap_support_stability")
  evidence <- tbl_db("tableS12_hydrosheaf_core_evidence_gates")
  tier_evidence <- tbl_db("tableS15_data_tier_reaction_evidence") |>
    mutate(
      tier = clean_tier(data_tier),
      reaction_label = factor(clean_reaction(reaction), levels = clean_reaction(unique(reaction)))
    )
  external <- tbl_db("tableS17_external_field_evidence_lifted_resolution") |>
    mutate(
      tier = clean_tier(data_tier),
      members_label = clean_class_label(members),
      dataset_label = clean_dataset(dataset)
    )
  phreeqc <- tbl_db("phreeqc_inverse_baseline")

  s1 <- ggplot(stoich_long(stoich), aes(ion, Reaction, fill = coefficient)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    heatmap_scale +
    labs(title = "Complete M5 reaction stoichiometric dictionary", x = NULL, y = NULL) +
    theme_m5_heatmap() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_m5(s1, file.path(supp_dir, "figureS1_r_reaction_dictionary.png"), 145, 150)

  diag_long <- rbind(
    data.frame(panel = diag$panel, metric = "Rank", value = diag$rank),
    data.frame(panel = diag$panel, metric = "Nullity", value = diag$nullity),
    data.frame(panel = diag$panel, metric = "Mutual coherence", value = diag$mutual_coherence)
  )
  s2 <- ggplot(diag_long, aes(panel, value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.68), width = 0.58) +
    scale_fill_manual(values = c("Rank" = "#008C7A", "Nullity" = "#C77C2B", "Mutual coherence" = "#3B6EA8")) +
    labs(title = "Structural identifiability diagnostics", x = NULL, y = "Value") +
    theme_m5() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_m5(s2, file.path(supp_dir, "figureS2_r_identifiability_diagnostics.png"), 145, 85)

  s3_dat <- fits |>
    group_by(method, panel) |>
    summarise(class_f1 = mean(class_f1, na.rm = TRUE), .groups = "drop") |>
    mutate(method_label = clean_method(method), panel_label = clean_panel(panel))
  s3 <- ggplot(s3_dat, aes(panel_label, method_label, fill = class_f1)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_viridis(option = "C", limits = c(0, 1), name = "Class F1") +
    labs(title = "Method performance across measured-ion panels", x = NULL, y = NULL) +
    theme_m5_heatmap() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_m5(s3, file.path(supp_dir, "figureS3_r_panel_method_heatmap.png"), 145, 95)

  rr <- summarise_recovery(recovery) |>
    filter(is.finite(f1)) |>
    mutate(reaction_label = clean_reaction(reaction))
  s4 <- ggplot(rr, aes(f1, reorder(reaction_label, f1), fill = family)) +
    geom_col(width = 0.68) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = "Hydrosheaf-Core reaction recovery profile", x = "Reaction F1", y = NULL) +
    theme_m5()
  save_m5(s4, file.path(supp_dir, "figureS4_r_reaction_recovery.png"), 145, 105)

  s5_dat <- nb |>
    group_by(candidate_measurement) |>
    summarise(
      score = mean(measurement_value_score, na.rm = TRUE),
      gain = mean(realised_class_f1_gain, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(score))
  s5 <- ggplot(s5_dat, aes(score, reorder(candidate_measurement, score), fill = gain)) +
    geom_col(width = 0.68) +
    scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C77C2B", midpoint = 0, name = "F1 gain") +
    labs(title = "Next-best sparse measurement ranking", x = "Measurement value score", y = NULL) +
    theme_m5()
  save_m5(s5, file.path(supp_dir, "figureS5_r_measurement_value.png"), 145, 90)

  reaction_cols <- setdiff(names(paths), c("scenario_id", "archetype", "lambda_l1", "lambda_l2", "rmse_mmolL", "n_selected"))
  reaction_cols <- reaction_cols[seq_len(min(length(reaction_cols), 10))]
  path_long <- do.call(
    rbind,
    lapply(reaction_cols, function(reaction) {
      data.frame(lambda_l1 = paths$lambda_l1, reaction = reaction, extent = abs(paths[[reaction]]))
    })
  ) |>
    group_by(lambda_l1, reaction) |>
    summarise(mean_extent = mean(extent, na.rm = TRUE), .groups = "drop") |>
    mutate(reaction_label = clean_reaction(reaction))
  s6 <- ggplot(path_long, aes(lambda_l1, mean_extent, colour = reaction_label)) +
    geom_line(linewidth = 0.28) +
    scale_x_log10() +
    labs(title = "Regularisation path stability", x = expression(lambda[1]), y = "Mean |extent|") +
    theme_m5()
  save_m5(s6, file.path(supp_dir, "figureS6_r_regularization_paths.png"), 145, 90)

  s7 <- bootstrap |>
    group_by(reaction) |>
    summarise(support_frequency = mean(support_frequency, na.rm = TRUE), .groups = "drop") |>
    mutate(reaction_label = clean_reaction(reaction)) |>
    ggplot(aes(support_frequency, reorder(reaction_label, support_frequency))) +
    geom_col(width = 0.68, fill = "#8C6BB1") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = "Bootstrap support stability", x = "Mean support frequency", y = NULL) +
    theme_m5()
  save_m5(s7, file.path(supp_dir, "figureS7_r_bootstrap_support.png"), 145, 100)

  s8 <- ggplot(
    evidence,
    aes(mean_evidence_score, median_penalty_scale, colour = family, size = true_active_fraction)
  ) +
    geom_point(alpha = 0.82) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = "Hydrosheaf-Core evidence gates", x = "Mean evidence score", y = "Median penalty scale") +
    theme_m5()
  save_m5(s8, file.path(supp_dir, "figureS8_r_core_evidence_gates.png"), 145, 95)

  s9 <- ggplot(tier_evidence, aes(tier, reaction_label, fill = mean_combined_evidence_score)) +
    geom_tile(colour = "white", linewidth = 0.18) +
    scale_fill_viridis(option = "C", limits = c(0, 1), name = "Evidence") +
    labs(title = "Data-tier reaction evidence", x = NULL, y = NULL) +
    theme_m5_heatmap()
  save_m5(s9, file.path(supp_dir, "figureS9_r_data_tier_evidence.png"), 145, 115)

  s10 <- ggplot(external, aes(median_elri, members_label, fill = dataset_label)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    labs(title = "External field-transfer ELRI", x = "Median ELRI", y = NULL) +
    theme_m5()
  save_m5(s10, file.path(supp_dir, "figureS10_r_external_field_elri.png"), 145, 95)

  s11 <- phreeqc |>
    group_by(archetype, phreeqc_inverse_setup) |>
    summarise(
      models = mean(models_found, na.rm = TRUE),
      class_f1 = mean(first_minimal_class_f1, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(setup_label = clean_setup(phreeqc_inverse_setup)) |>
    ggplot(aes(models + 1, class_f1, colour = archetype)) +
    geom_point(size = 1.9) +
    facet_wrap(~setup_label) +
    scale_x_log10(labels = comma) +
    scale_y_continuous(limits = c(0, 1)) +
    m5_colour_scale() +
    labs(title = "PHREEQC inverse baseline by archetype", x = "Mean feasible models + 1", y = "Class F1") +
    theme_m5()
  save_m5(s11, file.path(supp_dir, "figureS11_r_phreeqc_archetypes.png"), 145, 85)
}

figure1_database_design()
figure2_model_performance()
figure3_equifinality_elri()
figure4_data_tiers()
figure5_phreeqc_thermo()
figure6_field_transfer()
supplementary_figures()

message("Generated Nature-style R M5 figures in ", out_dir)
