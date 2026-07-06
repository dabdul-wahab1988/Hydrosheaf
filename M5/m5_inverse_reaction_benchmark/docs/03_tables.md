# M5 Table Evidence Map

| Display | File | Primary evidence | Use |
|---|---|---|---|
| Table 1 | `table1_comparative_inverse_performance.csv` | `benchmark_fits.csv`; `phreeqc_inverse_baseline.csv` | Main comparator performance at 3% noise with the full ion panel plus the conventional PHREEQC inverse baseline. |
| Table S1 | `tableS1_reaction_stoichiometry.csv` | `scripts/m5_common.py` | Reaction vectors, families, signs, and SI mappings. |
| Table S2 | `tableS2_reaction_equivalence_classes.csv` | `equivalence_classes.csv` | Exact signed equivalence classes. |
| Table S3 | `tableS3_phreeqc_scenario_parameters.csv` | `phreeqc_ground_truth.csv` | Scenario-level PHREEQC design and quality control. |
| Table S4 | `tableS4_hyperparameter_grid.csv` | `hyperparameter_selection.csv` | Predeclared calibration grid and selected penalties. |
| Table S5 | `tableS5_complete_model_metrics.csv` | `benchmark_fits.csv` | Full method, noise, panel, and archetype results. |
| Table S6 | `tableS6_reaction_confusion_matrices.csv` | `reaction_recovery.csv` | Reaction-specific support confusion matrices. |
| Table S7 | `tableS7_missing_ion_and_measurement_value.csv` | `next_best_measurement.csv` | Ion-specific held-out error and realised measurement value. |
| Table S8 | `tableS8_thermodynamic_bounds.csv` | `thermodynamic_threshold_sensitivity.csv` | SI-gate sensitivity and residual ambiguity. |
| Table S9 | `tableS9_software_environment.csv` | `analysis_summary.json` | Reproducibility metadata. |
| Table S10 | `tableS10_northern_ghana_summary.csv` | `ghana_field_pairs.csv` | Aquifer-stratified chemistry-only transfer summary. |
| Table S11 | `tableS11_phreeqc_inverse_baseline.csv` | `phreeqc_inverse_baseline.csv`; `phreeqc_inverse_baseline_models.csv` | Strict 5% and relaxed 20% PHREEQC inverse-model feasibility, multiplicity, and recovery by archetype. |
| Table S12 | `tableS12_hydrosheaf_core_evidence_gates.csv` | `hydrosheaf_core_evidence.csv` | Reaction-wise sparse-data evidence scores, penalty scales, and synthetic recovery checks for Hydrosheaf-Core. |
| Table S13 | `tableS13_ghana_hydrosheaf_core_evidence.csv` | `ghana_field_hydrosheaf_core_evidence.csv` | Field evidence-gate support frequencies by aquifer and reaction. |
| Table S14 | `tableS14_data_tier_experiment.csv` | `data_tier_experiment.csv` | Core, Plus-lite, and Enhanced measurement-tier recovery under controlled synthetic optional diagnostics. |
| Table S15 | `tableS15_data_tier_reaction_evidence.csv` | `data_tier_reaction_evidence.csv` | Reaction-level optional diagnostic evidence and penalty scales by data tier. |
| Table S16 | `tableS16_evidence_lifted_resolution.csv` | `data_tier_evidence_lifted_resolution.csv` | Entropy-based evidence-lifted resolution index for ambiguous reaction classes across Core, Plus-lite, and Enhanced tiers. |
| Table S17 | `tableS17_external_field_evidence_lifted_resolution.csv` | `external_field_evidence_lifted_resolution.csv` | Field-transfer ELRI audit for `NorthernGhana.xlsx`, Talensi, and Lower Anayari chemistry. |

## Results Database

The complete M5 results database is `results/m5_results.duckdb`. Its table
catalog is `results/m5_results_database_catalog.csv`, and the run-specific
figure/table/code manifest is written to `docs/m5_artifact_manifest.csv` and
`docs/m5_artifact_manifest.md`.
