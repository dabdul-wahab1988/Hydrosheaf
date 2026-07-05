# M5 Complete Output Manifest

## Analysis entry points

- `scripts/run_m5_all.py`: one-command analysis, tables, and figures.
- `scripts/run_m5_inverse_reaction_benchmark.py`: live-PHREEQC benchmark,
  factorial inversion, MRS, measurement value, bootstrap, and Ghana transfer.
- `scripts/make_m5_publication_tables.py`: one main and 16 supplementary tables.
- `scripts/make_m5_publication_figures.py`: six main and up to 18 supplementary figures.

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
- `results/analysis_summary.json`

## Displays

- Six main figures in `figures/figure1_*.{png,pdf}` through
  `figures/figure6_*.{png,pdf}`.
- Up to eighteen supplementary figures in
  `figures/supplementary/figureS1_*.{png,pdf}` through
  `figures/supplementary/figureS18_*.{png,pdf}`.
- `tables/table1_comparative_inverse_performance.csv`.
- `tables/tableS1_*.csv` through `tables/tableS16_*.csv`.

## Reproducibility and evidence maps

- `phreeqc_inputs/m5_factorial_benchmark.pqi`
- `phreeqc_inputs/m5_factorial_benchmark.out`
- `phreeqc_inputs/m5_factorial_selected.tsv`
- `phreeqc_inputs/inverse_baseline/M5S*_strict_5pct.{pqi,out}`
- `phreeqc_inputs/inverse_baseline/M5S*_relaxed_20pct.{pqi,out}`
- `docs/m5_results_summary.md`
- `docs/02_figures.md`
- `docs/03_tables.md`

## Claim rule

M5 evaluates Hydrosheaf as an identifiability-aware sparse linear inverse
reaction model with PHREEQC-generated truth, thermodynamic screening, and
predictive diagnostics. It does not validate a fully coupled nonlinear
PHREEQC inverse or reactive-transport solver.
