# M2 Hydrosheaf Benchmark Results Summary

Realisations: 100

This package was generated from `config/ground_truth.yaml` and is isolated from prior Hydrosheaf result folders.

## Table 4 Snapshot

| benchmark                           | data_source                                                                                 | target_variable                                                          | performance_metric                                                                             | expected_evidence                                                                                        | key_reference                                                                       | m2_status                                          |
|:------------------------------------|:--------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------|:---------------------------------------------------|
| Synthetic aquifer benchmark         | locked M2 ground truth                                                                      | transport and reaction extents                                           | median transport absolute error=0.058; median reaction error=0.112                             | direct recovery of known processes                                                                       | SyntheticDataGuide.docx                                                             | completed                                          |
| Public tracer-age validation        | USGS public-supply aquifer groundwater-age data release, DOI 10.5066/P9W7T0DN               | published TracerLPM mean age and young/Holocene/Pleistocene fractions    | M3 canonical n=1272; identifiability-gated n=356; median |log10 error|=0.022; within factor 2=0.92         | screening-level agreement with independent public LPM/TracerLPM-style age estimates                      | Jurgens et al. USGS data release 10.5066/P9W7T0DN                                   | completed by M3 identifiability-gated public-age benchmark |
| MODFLOW/MODPATH topology validation | USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5 archive, DOI 10.5066/F7J102FK | directed endpoint edges                                                   | prior-assisted TP=174; FP=0; FN=0; F1=1.00; no-prior F1=0.62                                  | prior result is ingestion fidelity; no-prior result tests endpoint-connectivity recovery, not pathline geometry or travel time | Harte USGS data release 10.5066/F7J102FK                                            | completed topology-only comparison                 |
| Live PHREEQC forward validation     | USGS PHREEQC version 3 examples and databases, DOI 10.3133/tm6A43                           | major-ion concentration evolution and saturation-index residuals         | current proxy feasible fraction=0.88; planned: live PHREEQC RMSE, NSE/equivalent, SI residuals | inferred reaction pathways remain forward-feasible in PHREEQC                                            | Parkhurst and Appelo 2013 USGS TM 6-A43                                             | proxy completed; live external pending             |
| Data-limited pilot scenario         | Lower Anayari (Manu, 41 source samples) and Talensi (63 source samples) field pilot datasets  | end-to-end generated-edge and reaction outputs under sparse field inputs | n_edges=258; median chemistry R²=0.703; generated graph, no independent process truth          | workflow remains interpretable when optional tracers, source graph edges, or PHREEQC inputs are absent; sampling-density metadata were unavailable | data/FieldData/LowerAnayari/manu.csv; data/FieldData/Talensi_MiningArea/talensi.csv | completed field-hydrochemistry demonstration       |

## Interpretation Guardrails

- Public USGS age results are screening-level residence-time evidence; their reported log10 R² is the residual-based 1 - SS_res/SS_tot statistic, not squared Pearson correlation and not full TracerLPM family equivalence.
- MODPATH evidence validates endpoint/pathline topology conversion only, not geochemical process inference.
- PHREEQC evidence remains a proxy unless a live PHREEQC executable/backend is configured and rerun.
- Ghana field results are a generated-graph field-hydrochemistry demonstration without independent process-truth labels.

## Main Outputs

- `results/transport_recovery.csv`
- `results/reaction_recovery.csv`
- `results/missing_data_sensitivity.csv`
- `results/topology_robustness.csv`
- `results/age_inference_validation.csv`
- `results/phreeqc_forward_validation.csv`
- `tables/table1_module_architecture.csv` through `tables/table5_method_comparison.csv`
- manuscript-ready figures are generated by `scripts/make_publication_figures.py` and `scripts/make_supplementary_figures.py`; M2 Figure 6 is the field-filtering/configuration figure and does not consume M6 regularization output

## PHREEQC Note

The forward-validation table uses a deterministic PHREEQC-compatible mass-balance proxy because the benchmark is generated from locked saturation fields. If a PHREEQC executable is configured later, replace `linear_mass_balance_phreeqc_proxy` with full PHREEQC kinetic simulations while keeping the same input and output schema.
