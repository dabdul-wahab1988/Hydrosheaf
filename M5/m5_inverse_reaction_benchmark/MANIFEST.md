# M5 Complete Output Manifest

## Analysis entry points

- `scripts/run_m5_all.py`: one-command analysis, tables, and figures.
- `scripts/run_m5_inverse_reaction_benchmark.py`: live-PHREEQC benchmark,
  factorial inversion, MRS, measurement value, bootstrap, and Ghana transfer.
- `scripts/make_m5_publication_tables.py`: one main and 17 supplementary tables.
- `scripts/make_m5_publication_figures.py`: six main and up to 18 supplementary figures.
- `scripts/export_m5_results_database.py`: single DuckDB evidence database and
  table catalog from all generated CSV/JSON result artefacts.
- `scripts/make_m5_database_figures.py`: database-backed Python fallback/QC
  publication figures.
- `r_figures/plot_m5_publication_figures.R`: Nature-style R manuscript figures
  from the current M5 results database or its generated CSV/JSON mirrors.
- `scripts/export_m5_artifact_manifest.py`: run-specific manifest for figures,
  tables, evidence maps, and plotting code.

## Core results

- `results/phreeqc_ground_truth.csv`
- `results/benchmark_fits.csv`
- `results/reaction_recovery.csv`
- `results/heldout_ion_results.csv`
- `results/equivalence_classes.csv`
- `results/identifiability_diagnostics.csv`
- `results/hyperparameter_selection.csv`
- `results/regularization_paths.csv`
- `results/mechanism_resolution_scores.csv`
- `results/mrs_calibration_model.json`
- `results/next_best_measurement.csv`
- `results/bootstrap_support_stability.csv`
- `results/thermodynamic_threshold_sensitivity.csv`
- `results/phreeqc_inverse_baseline.csv`
- `results/phreeqc_inverse_baseline_models.csv`
- `results/hydrosheaf_core_evidence.csv`
- `results/hydrosheaf_core_evidence_lifted_resolution.csv`
- `results/data_tier_experiment.csv`
- `results/data_tier_reaction_evidence.csv`
- `results/data_tier_evidence_lifted_resolution.csv`
- `results/data_tier_optional_diagnostics.csv`
- `results/ghana_field_pairs.csv`
- `results/ghana_field_reaction_extents.csv`
- `results/ghana_field_heldout_ions.csv`
- `results/ghana_field_class_support.csv`
- `results/ghana_field_hydrosheaf_core_evidence.csv`
- `results/ghana_evidence_lifted_resolution.csv`
- `results/external_field_transfer_pairs.csv`
- `results/external_field_reaction_evidence.csv`
- `results/external_field_evidence_lifted_resolution.csv`
- `results/m5_results.duckdb`
- `results/m5_results_database_catalog.csv`
- `results/analysis_summary.json`

## Displays

- Intended manuscript figures: six R Nature-style main figures in
  `figures/r_publication/figure1_r_*.{png,tif,pdf}` through
  `figures/r_publication/figure6_r_*.{png,tif,pdf}`.
- R supplementary figure set in
  `figures/r_publication/supplementary/figureS1_r_*.{png,tif,pdf}` through
  `figures/r_publication/supplementary/figureS11_r_*.{png,tif,pdf}`.
- Python database-backed fallback/QC figures in
  `figures/publication/figure1_database_*.{png,pdf}` through
  `figures/publication/figure6_database_*.{png,pdf}` and
  `figures/publication/supplementary/figureS1_database_*.{png,pdf}` through
  `figures/publication/supplementary/figureS10_database_*.{png,pdf}`.
- Legacy six main figures in `figures/figure1_*.{png,pdf}` through
  `figures/figure6_*.{png,pdf}`.
- Up to eighteen supplementary figures in
  `figures/supplementary/figureS1_*.{png,pdf}` through
  `figures/supplementary/figureS18_*.{png,pdf}`.
- `tables/table1_comparative_inverse_performance.csv`.
- `tables/tableS1_*.csv` through `tables/tableS17_*.csv`.

## Reproducibility and evidence maps

- `phreeqc_inputs/m5_factorial_benchmark.pqi`
- `phreeqc_inputs/m5_factorial_benchmark.out`
- `phreeqc_inputs/m5_factorial_selected.tsv`
- `phreeqc_inputs/inverse_baseline/M5S*_strict_5pct.{pqi,out}`
- `phreeqc_inputs/inverse_baseline/M5S*_relaxed_20pct.{pqi,out}`
- `docs/m5_results_summary.md`
- `docs/02_figures.md`
- `docs/03_tables.md`
- `docs/m5_artifact_manifest.csv`
- `docs/m5_artifact_manifest.md`

## Claim rule

M5 evaluates Hydrosheaf as an identifiability-aware sparse linear inverse
reaction model with PHREEQC-generated truth, thermodynamic screening, and
predictive diagnostics. It does not validate a fully coupled nonlinear
PHREEQC inverse or reactive-transport solver.
