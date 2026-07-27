# Supplementary Information: Figures and Tables

Field transfer of Hydrosheaf under data scarcity: robustness, uncertainty and limits across Ghanaian aquifers

This document contains the supplementary figures and supplementary tables referenced in the main text and in the Supplementary Methods. Extended Data Figures 1-3 are reported with the main manuscript, not repeated here. All supplementary figures use the shared Nature-portfolio theme described in the Methods; all supplementary tables are also deposited as machine-readable CSV alongside this document.

## Supplementary Figures

### Supplementary Figure S1

![Supplementary Figure S1. Hydrochemical context. (a) Gibbs cation ratio. (b) Gibbs anion ratio. (c) Major-ion composition. (d) Mineralisation (TDS) by dataset.](figures/r_publication/supplementary/figureS1_hydrochem_context.png)

**Supplementary Figure S1.** Hydrochemical context. (a) Gibbs cation ratio. (b) Gibbs anion ratio. (c) Major-ion composition. (d) Mineralisation (TDS) by dataset.

*Technical description:* 2x2 panel figure computed from harmonised major-ion concentrations for all three datasets.

*Interpretation:* Places the three datasets within established Gibbs water-rock-interaction classification; supports, but does not substitute for, the class-level reaction inference reported in the main text.

### Supplementary Figure S2

![Supplementary Figure S2. Charge balance. (a) Charge-balance-error distribution with +-5%/+-10% guides. (b) Quality-class composition by dataset.](figures/r_publication/supplementary/figureS2_charge_balance.png)

**Supplementary Figure S2.** Charge balance. (a) Charge-balance-error distribution with +-5%/+-10% guides. (b) Quality-class composition by dataset.

*Technical description:* Two-panel figure computed from the independently recomputed charge-balance error for all three datasets.

*Interpretation:* Visualises the quantitative/screening/exploratory split underlying Table 1 and identifies Talensi's charge-balance limitation.

### Supplementary Figure S3

![Supplementary Figure S3. Seasonal signal. (a) Dry/wet evapoconcentration factor by region. (b) Reactive residual after transport correction by region.](figures/r_publication/supplementary/figureS3_seasonal_change.png)

**Supplementary Figure S3.** Seasonal signal. (a) Dry/wet evapoconcentration factor by region. (b) Reactive residual after transport correction by region.

*Technical description:* Two-panel figure computed from the Cl-conservative transport correction applied to the 160 Northern Ghana well pairs.

*Interpretation:* Shows the magnitude of the seasonal transport signal that is removed before reaction inference, by administrative region (no independent aquifer-type classification exists for these boreholes).

### Supplementary Figure S4

![Supplementary Figure S4. Alternative edge networks. (a) Inferred process-network composition per edge set. (b) Divergence (total variation distance) and stability versus the chemistry-kNN reference set.](figures/r_publication/supplementary/figureS4_edge_networks.png)

**Supplementary Figure S4.** Alternative edge networks. (a) Inferred process-network composition per edge set. (b) Divergence (total variation distance) and stability versus the chemistry-kNN reference set.

*Technical description:* Two-panel figure comparing three Hydrosheaf-generated edge sets (chemistry-kNN, geographic-nearest, random/perturbed) against each other. An earlier revision also compared against a fourth, "provided graph" edge set imported from a retired antecedent-study workbook; it has been removed rather than substituted (Methods, DECISIONS.md).

*Interpretation:* Shows that network-level process attribution is sensitive to the assumed connectivity even though point-level identifiability is not (Figure 4c, Extended Data context).

### Supplementary Figure S5

![Supplementary Figure S5. Dominant-process composition across evidence tiers.](figures/r_publication/supplementary/figureS5_tier_family_heatmap.png)

**Supplementary Figure S5.** Dominant-process composition across evidence tiers.

*Technical description:* Single heatmap of dominant reaction family versus evidence tier for the 160 Northern Ghana wells.

*Interpretation:* Provides the full family-composition detail underlying the constant 35.6% family-flip fraction reported in the Results.

### Supplementary Figure S6

![Supplementary Figure S6. Diagnostic distributions. (a) Bootstrap support stability. (b) Mechanism Resolution Score.](figures/r_publication/supplementary/figureS6_diagnostic_distributions.png)

**Supplementary Figure S6.** Diagnostic distributions. (a) Bootstrap support stability. (b) Mechanism Resolution Score.

*Technical description:* Two-panel distributional figure across all Northern Ghana wells and tiers.

*Interpretation:* Documents that support stability and Mechanism Resolution Score remained high and narrow even where identifiability class collapsed under tier ablation.

### Supplementary Figure S7

![Supplementary Figure S7. Reactive residual versus conservative evapoconcentration, coloured by dominant process.](figures/r_publication/supplementary/figureS7_reactive_vs_evapo.png)

**Supplementary Figure S7.** Reactive residual versus conservative evapoconcentration, coloured by dominant process.

*Technical description:* Single scatter figure relating the fitted reactive residual magnitude to the estimated evapoconcentration factor.

*Interpretation:* Shows that reactive and conservative (transport) signals are separable across the dominant process classes recovered.

### Supplementary Figure S8

![Supplementary Figure S8. Mechanism Resolution Score versus support stability, coloured by dominant process.](figures/r_publication/supplementary/figureS8_mrs_diagnostic.png)

**Supplementary Figure S8.** Mechanism Resolution Score versus support stability, coloured by dominant process.

*Technical description:* Single scatter figure relating the two principal reliability diagnostics across all Northern Ghana wells.

*Interpretation:* Illustrates that Mechanism Resolution Score and support stability are correlated but not redundant diagnostics.

### Supplementary Figure S9

![Supplementary Figure S9. External-dataset transfer detail, faceted by dataset (Talensi, Lower Anayari).](figures/r_publication/supplementary/figureS9_external_detail.png)

**Supplementary Figure S9.** External-dataset transfer detail, faceted by dataset (Talensi, Lower Anayari).

*Technical description:* Per-edge diagnostic detail figure underlying Table 4 and Extended Data Figure 1.

*Interpretation:* Provides edge-level detail supporting the external-transfer summary statistics reported in the main text.

### Supplementary Figure S10

![Supplementary Figure S10. Identifiability composition per process family and dataset, faceted by dataset.](figures/r_publication/supplementary/figureS10_limitation_detail.png)

**Supplementary Figure S10.** Identifiability composition per process family and dataset, faceted by dataset.

*Technical description:* Detailed faceted heatmap underlying the summary limitation map in Figure 6a.

*Interpretation:* Provides the full per-family, per-dataset breakdown supporting the main-text limitation map.

### Supplementary Figure S11

![Supplementary Figure S11. Competing no-flow explanation sensitivity: fractions of evaluated edges exceeding the 0.5 and 0.8 HydroSheaf screening thresholds by dataset and candidate-edge construction.](figures/r_publication/supplementary/figureS11_null_sensitivity.png)

**Supplementary Figure S11.** Competing no-flow explanation sensitivity: fractions of evaluated edges exceeding the 0.5 and 0.8 HydroSheaf screening thresholds by dataset and candidate-edge construction.

*Technical description:* Bar-chart figure summarising the null-model sensitivity screen described in Supplementary Methods Section S12.

*Interpretation:* These are sensitivity flags, not calibrated field probabilities; high flag rates further restrict interpretation of inferred edges as physical flow paths.

## Supplementary Tables

### Supplementary Table S2

**Supplementary Table S2.** Charge-balance quality classes by dataset.

| dataset | exploratory | quantitative | screening |
| --- | --- | --- | --- |
| manu | 0 | 36 | 5 |
| northern_ghana | 7 | 294 | 19 |
| talensi | 58 | 0 | 5 |

*Technical description:* Full per-dataset counts underlying the quantitative/screening/exploratory summary in Table 1.

*Interpretation:* Confirms Talensi's charge-balance limitation quantitatively at full resolution.

### Supplementary Table S3

**Supplementary Table S3.** Northern Ghana region summary: wells, mean Mechanism Resolution Score, mean stability, fraction partially identifiable, and top dominant family, by administrative region. No independent aquifer-type classification exists for these boreholes (Methods, DECISIONS.md); an earlier revision reported this table by aquifer type and included a fraction-concordant-with-prior-labels column derived from a retired antecedent-study workbook, which has been removed rather than substituted.

| region | n_wells | mean_mrs | mean_stability | frac_partial | top_family |
| --- | --- | --- | --- | --- | --- |
| North East Region | 40 | 70.768 | 0.965 | 1.0 | silicate |
| Northern Region | 40 | 71.219 | 0.974 | 1.0 | silicate |
| Upper East Region | 40 | 71.096 | 0.976 | 1.0 | silicate |
| Upper West Region | 40 | 70.585 | 0.966 | 1.0 | silicate |

*Technical description:* Aggregated from the 160-well Tier-4 seasonal transfer experiment.

*Interpretation:* Provides the per-region detail underlying the seasonal/region-stratified transfer results reported in the main text.

### Supplementary Table S4

**Supplementary Table S4.** Full tier-ablation identifiability counts (non-identifiable versus partially identifiable) at every evidence tier.

| tier | non_identifiable | partially_identifiable |
| --- | --- | --- |
| tier0_majors | 96 | 64 |
| tier1_isotopes | 83 | 77 |
| tier2_fluoride | 83 | 77 |
| tier3_sr_sio2 | 1 | 159 |
| tier4_full_metadata | 0 | 160 |

*Technical description:* Raw well counts underlying the fractional results reported in Table 3.

*Interpretation:* Gives the complete count-level detail behind the pivotal tier-ablation result.

### Supplementary Table S5

**Supplementary Table S5.** Edge-set sensitivity summary: edges, total variation distance versus the chemistry-kNN reference set, mean Mechanism Resolution Score, mean stability, and fraction partially identifiable, for all three Hydrosheaf-generated edge sets. An earlier revision also compared against a fourth, "provided graph" edge set imported from a retired antecedent-study workbook, which has been removed rather than substituted (Methods, DECISIONS.md).

| edge_set | n_edges | tvd_vs_chemistry_knn | mean_mrs | mean_stability | frac_partial |
| --- | --- | --- | --- | --- | --- |
| chemistry_knn | 140 | 0.0 | 73.339 | 0.972 | 0.993 |
| geographic_nearest | 140 | 0.121 | 72.644 | 0.965 | 1.0 |
| random_perturbed | 140 | 0.050 | 72.756 | 0.967 | 1.0 |

*Technical description:* Computed from the edge-uncertainty experiment across the chemistry-kNN, geographic-nearest and random/perturbed edge sets.

*Interpretation:* Supports the Results claim that network composition is edge-set-sensitive while point-level identifiability is edge-invariant.

### Supplementary Table S6

**Supplementary Table S6.** External-dataset transfer summary by dominant reaction family for Talensi and Lower Anayari. The full per-edge output (source, target, dominant family, resolution class, MRS, support stability, reactive RMSE, and support count for all 129 Talensi and 85 Lower Anayari evaluated edges) underlying this summary is deposited as `results/m6_external_transfer.csv` alongside this document rather than reproduced row-by-row here.

| dataset | native_tier | n_edges | dominant family (n) |
| --- | --- | --- | --- |
| talensi | tier1_isotopes | 129 | silicate (39), redox (37), anthropogenic (37), evaporite (6), ion_exchange (5), carbonate (5) |
| manu | tier2_fluoride | 85 | silicate (80), carbonate (2), redox (2), ion_exchange (1) |

*Technical description:* Aggregated from the full per-edge output table underlying the external-transfer summary in Table 4; full per-edge detail is in `results/m6_external_transfer.csv`.

*Interpretation:* Confirms that Lower Anayari's edges are overwhelmingly silicate-dominant (80 of 85), the basis for the strontium/silica next-best-measurement recommendation reported in the main text; Talensi's dominant-family distribution is more heterogeneous, consistent with its additional charge-balance limitation.

### Supplementary Table S7

**Supplementary Table S7.** Uncertainty-metric definitions used throughout the study.

| Metric | Definition |
| --- | --- |
| Support stability | Mean Jaccard of selected reaction support across 5% analytical-noise bootstrap resamples (0-1, higher = more stable) |
| MRS | Mechanism Resolution Score (0-100) from frozen M5 calibration |
| Held-out ion RMSE | Leave-one-ion predictive error of the concentration change |
| Reactive RMSE | Residual after Cl-conservative transport correction and inversion |
| Identifiability class | non/partially/equivalence-class/identifiable, transferred classifier + evidence-lift gate |
| Evidence corroboration | Whether the dominant process family is supported by an available tracer (isotopes/Sr/SiO2/F/SI) |

*Technical description:* Reference table of the six principal uncertainty and reliability metrics and their definitions.

*Interpretation:* Provides a single reference point for interpreting the uncertainty metrics reported in every experiment.

### Supplementary Table S8

**Supplementary Table S8.** Software and computational environment.

| Component | Version / setting |
| --- | --- |
| Python | 3.14.6 |
| numpy | 2.4.6 |
| pandas | 2.3.3 |
| Inverse solver | M5 fit_inverse (FISTA, column-normalised, SI-bounded) |
| MRS calibration | frozen M5 mrs_calibration_model.json (not re-fit) |
| Seed | 1234 |

*Technical description:* Records the Python, numpy and pandas versions, the frozen inverse solver and Mechanism Resolution Score calibration, and the fixed random seed.

*Interpretation:* Supports full reproducibility of every reported result.

### Supplementary Table S9

**Supplementary Table S9.** Competing no-flow explanation sensitivity by dataset and candidate-edge construction. Northern Ghana carries no `null_common_lithology` flags because no independent aquifer-type, geology-group or lithology classification exists for these boreholes (Methods, DECISIONS.md); Talensi and Lower Anayari (`manu`) retain a fixed, site-level lithology label and so still contribute to that flag. An earlier revision also scored a fourth Northern Ghana edge set, "provided_graph", imported from a retired antecedent-study workbook, which has been removed rather than substituted.

| dataset | edge_set | n_edges | mean_null_score | median_null_score | fraction_gt_0_5 | fraction_gt_0_8 | fraction_retained_at_0_5 | top_null_flags |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| manu | chemistry_knn | 95 | 0.973 | 1.0 | 1.0 | 0.968 | 0.0 | null_chemistry_similar:95;null_common_lithology:95;null_isotope_proximity:51;null_shared_recharge:46;null_common_anthropogenic:35;null_spatial_autocorr:12 |
| northern_ghana | chemistry_knn | 140 | 0.992 | 1.0 | 1.0 | 1.0 | 0.0 | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:92;null_shared_recharge:82 |
| northern_ghana | geographic_nearest | 140 | 0.961 | 1.0 | 1.0 | 0.957 | 0.0 | null_chemistry_similar:140;null_common_anthropogenic:140;null_shared_recharge:90;null_isotope_proximity:88;null_spatial_autocorr:11 |
| northern_ghana | random_perturbed | 140 | 0.928 | 0.958 | 0.993 | 0.921 | 0.007 | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:83;null_shared_recharge:79 |
| talensi | chemistry_knn | 135 | 0.899 | 1.0 | 0.956 | 0.815 | 0.044 | null_common_lithology:135;null_chemistry_similar:130;null_common_anthropogenic:114;null_shared_recharge:79;null_isotope_proximity:70;null_spatial_autocorr:3 |

*Technical description:* Mean and median null scores, and fractions exceeding the 0.5 and 0.8 screening thresholds, for each dataset/edge-set combination.

*Interpretation:* Screening flags for competing non-flow explanations of edge similarity; not calibrated field probabilities.

