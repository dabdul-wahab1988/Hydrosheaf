# M5 Figure Evidence Map

## Intended R manuscript figures

The R figures in `figures/r_publication/` are the intended Nature-style
manuscript display layer. They are generated from `results/m5_results.duckdb`
when the R `duckdb` package is available, with a verified fallback to the
generated CSV/JSON mirrors from the same run.

| Display | File stem | Primary evidence | Defensible interpretation |
|---|---|---|---|
| Figure 1 | `figure1_r_database_design` | `analysis_summary.json`; `m5_results_database_catalog.csv`; `benchmark_fits.csv` | Documents the reproducible M5 evidence architecture and database scale. |
| Figure 2 | `figure2_r_model_performance` | `benchmark_fits.csv`; `phreeqc_inverse_baseline.csv` | Hydrosheaf Guarded is best on class recovery, but gains are moderate and must be interpreted under equifinality. |
| Figure 3 | `figure3_r_equifinality_elri` | `equivalence_classes.csv`; `tableS1_reaction_stoichiometry.csv`; `data_tier_evidence_lifted_resolution.csv` | Sparse chemistry has exact non-uniqueness; evidence lifts ambiguity conditionally rather than creating new mass-balance uniqueness. |
| Figure 4 | `figure4_r_data_tiers` | `data_tier_experiment.csv`; `data_tier_evidence_lifted_resolution.csv`; `tableS15_data_tier_reaction_evidence.csv` | Plus-lite and enhanced diagnostics reduce overinterpretation and improve conditional identifiability. |
| Figure 5 | `figure5_r_phreeqc_thermo` | `thermodynamic_threshold_sensitivity.csv`; `phreeqc_inverse_baseline.csv` | Thermodynamic gates and conventional PHREEQC inverse modelling expose feasible multiplicity rather than eliminating it. |
| Figure 6 | `figure6_r_field_transfer` | `ghana_field_pairs.csv`; `external_field_evidence_lifted_resolution.csv`; `tableS17_external_field_evidence_lifted_resolution.csv` | Ghana and external field datasets are plausibility/transfer audits, not reaction-truth validation. |

R Supplementary Figures S1-S11 report the reaction dictionary, structural
diagnostics, panel-by-method heatmap, Hydrosheaf-Core reaction recovery,
measurement-value ranking, regularisation paths, bootstrap support, evidence
gates, data-tier evidence, external-field ELRI, and PHREEQC archetype behaviour.

## Python fallback and QC figures

| Display | File stem | Primary evidence | Defensible interpretation |
|---|---|---|---|
| Figure 1 | `figure1_identifiability_workflow` | `phreeqc_ground_truth.csv`; `benchmark_fits.csv`; `heldout_ion_results.csv`; `ghana_field_pairs.csv` | Documents benchmark design and evidence scale. |
| Figure 2 | `figure2_fit_vs_mechanism` | `benchmark_fits.csv`; `reaction_recovery.csv` | Low residuals frequently coexist with wrong exact-phase support. |
| Figure 3 | `figure3_regularization_and_mrs` | `regularization_paths.csv`; `identifiability_diagnostics.csv`; `mechanism_resolution_scores.csv` | Regularisation and structural rank affect support; MRS discrimination is moderate on the held-out archetype. |
| Figure 4 | `figure4_measurement_value` | `next_best_measurement.csv`; `benchmark_fits.csv` | Held-out prediction identifies measurements with realised ambiguity-reduction value. |
| Figure 5 | `figure5_thermodynamic_screening` | `benchmark_fits.csv`; `thermodynamic_threshold_sensitivity.csv` | Thermodynamic bounds remove incompatible directions but do not resolve feasible equifinality. |
| Figure 6 | `figure6_northern_ghana` | `ghana_field_pairs.csv`; `ghana_field_heldout_ions.csv`; `ghana_field_class_support.csv`; `ghana_field_hydrosheaf_core_evidence.csv` | Northern Ghana is a chemistry-only seasonal transfer showing partial identifiability, not reaction truth. |

Supplementary Figures S1-S18 report the full reaction dictionary, structural
diagnostics, PHREEQC quality control, noise sensitivity, phase bias,
regularisation paths, MRS transfer, ion ablation, measurement rankings,
thermodynamic thresholds, bootstrap stability, runtime, Ghana quality control,
Ghana class ensembles, Hydrosheaf-Core evidence gates, Ghana field evidence,
and Core/Plus-lite/Enhanced data-tier performance.
