# M5 Figure Evidence Map

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
