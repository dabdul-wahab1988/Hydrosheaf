# M6 Null-Model Sensitivity Extension

This is a separate edge-level sensitivity analysis. The locked M6 field-transfer
outputs are unchanged. Null scores quantify whether chemistry, recharge, metadata,
or common-source similarity creates a competing no-flow explanation.

## Outputs

- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M6/m6_field_transfer_benchmark/results/m6_null_edge_scores.csv`: edge-level scores and flags.
- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M6/m6_field_transfer_benchmark/results/m6_null_sensitivity_summary.csv`: dataset/edge-set summaries.
- `C:/Users/DicksonAbdul-Wahab/Desktop/NeutroProject/Groundwater/Hydrosheaf/M6/m6_field_transfer_benchmark/results/m6_null_model_run.json`: configuration and source hashes.

Scores >0.5 are reported as a competing no-flow explanation; scores >0.8 are
reported as a dominant no-flow explanation. These are screening thresholds from
Hydrosheaf's evidence-classification logic, not calibrated field probabilities.

## Summary

| dataset        | edge_set           |   n_edges |   mean_null_score |   median_null_score |   fraction_gt_0_5 |   fraction_gt_0_8 |   fraction_retained_at_0_5 | top_null_flags                                                                                                                                                         |
|:---------------|:-------------------|----------:|------------------:|--------------------:|------------------:|------------------:|---------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| manu           | chemistry_knn      |        95 |          0.973052 |            1        |          1        |          0.968421 |                 0          | null_chemistry_similar:95;null_common_lithology:95;null_isotope_proximity:51;null_shared_recharge:46;null_common_anthropogenic:35;null_spatial_autocorr:12             |
| northern_ghana | chemistry_knn      |       140 |          0.996733 |            1        |          1        |          1        |                 0          | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:92;null_shared_recharge:82;null_common_lithology:70;null_common_lithology_explicit:24  |
| northern_ghana | geographic_nearest |       140 |          0.964243 |            1        |          1        |          0.95     |                 0          | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:100;null_shared_recharge:90;null_common_lithology:55;null_common_lithology_explicit:13 |
| northern_ghana | provided_graph     |       140 |          0.97048  |            1        |          1        |          0.95     |                 0          | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:98;null_shared_recharge:90;null_common_lithology:74;null_common_lithology_explicit:22  |
| northern_ghana | random_perturbed   |       140 |          0.942395 |            0.996161 |          0.992857 |          0.928571 |                 0.00714286 | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:81;null_shared_recharge:79;null_common_lithology:40;null_common_lithology_explicit:10  |
| talensi        | chemistry_knn      |       135 |          0.898755 |            1        |          0.955556 |          0.814815 |                 0.0444444  | null_common_lithology:135;null_chemistry_similar:130;null_common_anthropogenic:114;null_shared_recharge:79;null_isotope_proximity:70;null_spatial_autocorr:3           |

## Guardrail

This extension does not establish field process truth. It tests sensitivity to
non-flow explanations and should be interpreted alongside the locked M6 evidence
gate, edge-set sensitivity, and identifiability results.
