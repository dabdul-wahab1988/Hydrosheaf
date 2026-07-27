# M6 Null-Model Sensitivity Extension

This is a separate edge-level sensitivity analysis. The locked M6 field-transfer
outputs are unchanged. Null scores quantify whether chemistry, recharge, metadata,
or common-source similarity creates a competing no-flow explanation.

## Outputs

- `M6/m6_field_transfer_benchmark/results/m6_null_edge_scores.csv`: edge-level scores and flags.
- `M6/m6_field_transfer_benchmark/results/m6_null_sensitivity_summary.csv`: dataset/edge-set summaries.
- `M6/m6_field_transfer_benchmark/results/m6_null_model_run.json`: configuration and source hashes.

Scores >0.5 are reported as a competing no-flow explanation; scores >0.8 are
reported as a dominant no-flow explanation. These are screening thresholds from
Hydrosheaf's evidence-classification logic, not calibrated field probabilities.

## Summary

| dataset | edge_set | n_edges | mean_null_score | median_null_score | fraction_gt_0_5 | fraction_gt_0_8 | fraction_retained_at_0_5 | top_null_flags |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| manu | chemistry_knn | 95 | 0.973052376639718 | 1.0 | 1.0 | 0.968421052631579 | 0.0 | null_chemistry_similar:95;null_common_lithology:95;null_isotope_proximity:51;null_shared_recharge:46;null_common_anthropogenic:35;null_spatial_autocorr:12 |
| northern_ghana | chemistry_knn | 140 | 0.9916133281807203 | 1.0 | 1.0 | 1.0 | 0.0 | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:92;null_shared_recharge:82 |
| northern_ghana | geographic_nearest | 140 | 0.9610743892003798 | 1.0 | 1.0 | 0.9571428571428572 | 0.0 | null_chemistry_similar:140;null_common_anthropogenic:140;null_shared_recharge:90;null_isotope_proximity:88;null_spatial_autocorr:11 |
| northern_ghana | random_perturbed | 140 | 0.9279353910126653 | 0.9583766043529961 | 0.9928571428571429 | 0.9214285714285714 | 0.007142857142857143 | null_chemistry_similar:140;null_common_anthropogenic:140;null_isotope_proximity:83;null_shared_recharge:79 |
| talensi | chemistry_knn | 135 | 0.8987553350534447 | 1.0 | 0.9555555555555556 | 0.8148148148148148 | 0.044444444444444446 | null_common_lithology:135;null_chemistry_similar:130;null_common_anthropogenic:114;null_shared_recharge:79;null_isotope_proximity:70;null_spatial_autocorr:3 |

## Guardrail

This extension does not establish field process truth. It tests sensitivity to
non-flow explanations and should be interpreted alongside the locked M6 evidence
gate, edge-set sensitivity, and identifiability results.
