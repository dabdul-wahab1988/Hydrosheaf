# M3 Phase-2 Design Matrix QA — superseded pre-reportability snapshot

> **Do not use these values for the manuscript.** This 2026-07-27 snapshot
> predates the reference-free identifiability/reportability guard and the final
> correctness lock. The authoritative design-matrix outputs are
> `m3_design_matrix_results.csv`, `m3_design_matrix_summary.csv`, and
> `m3_design_matrix_manifest.json`, regenerated on 2026-07-28.

- Generated: 2026-07-27T07:11:46.326550+00:00
- Config: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `M3\m3_age_benchmark\results\m3_design_matrix_results.csv`
- Output rows: 17808
- Unique scenarios: 14
- Row limit: full

## Scenario Metrics

```text
                      scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
                  ablation_no_he4        1272         1216              1216                0.190502    0.827294         0.578947          0.857730                    0  0.575369
                 ablation_raw_c14        1272         1247              1247                0.203143    0.865882         0.583801          0.842021                  276  0.564126
               ablation_raw_gases        1272         1249              1249                0.221153    0.898020         0.568455          0.840673                  278  0.530650
   hydrosheaf_selection_corrected        1272          481               481                0.440724    0.952317         0.384615          0.758836                  127  0.271992
            oldwater_c14_ensemble        1272         1248              1248                0.379409    1.225036         0.443109          0.750801                  278  0.127215
oldwater_ensemble_he4_uncertainty        1272         1265              1265                0.339840    1.193346         0.464032          0.766798                  295  0.166188
         oldwater_he4_uncertainty        1272         1266              1266                0.190502    0.833364         0.589258          0.859400                  295  0.593074
        parity_reported_corrected        1272         1249              1249                0.221153    0.898020         0.568455          0.840673                  278  0.530650
               screened_dgm_gases        1272         1249              1249                0.221153    0.898020         0.568455          0.840673                  278  0.530650
                  tracer_old_only        1272         1131              1131                0.767507    1.615145         0.346596          0.541998                  735 -0.530869
                tracer_young_only        1272         1169              1169                1.296522    1.682321         0.209581          0.425150                    0 -0.632471
    tracerlpm_parity_agefractions        1272         1249              1249                0.166134    0.798768         0.630104          0.852682                  278  0.628278
   tracerlpm_parity_hier_oldwater        1272         1250              1250                0.362071    1.205931         0.458400          0.751200                  278  0.155479
          tracerlpm_strict_parity        1272         1251              1251                0.205746    0.894837         0.573141          0.842526                  278  0.534661
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
                  ablation_no_he4         1199                  0.000000               -0.047313           0.120934                    30                  22
                 ablation_raw_c14         1247                  0.000000               -0.011236           0.172414                    28                  10
               ablation_raw_gases         1249                  0.000000                0.000000           0.000000                     0                   0
   hydrosheaf_selection_corrected          477                  0.011691               -0.005248           0.448637                    55                  65
            oldwater_c14_ensemble         1248                  0.000000                0.246874           0.087340                     5                 162
oldwater_ensemble_he4_uncertainty         1248                  0.000000                0.211503           0.186699                    38                 173
         oldwater_he4_uncertainty         1249                  0.000000               -0.046478           0.112890                    35                  11
               screened_dgm_gases         1249                  0.000000                0.000000           0.000000                     0                   0
                  tracer_old_only         1128                  0.105965                0.624124           0.138298                    29                 272
                tracer_young_only         1153                  0.555240                0.820567           0.261058                    91                 481
    tracerlpm_parity_agefractions         1245                  0.000000               -0.061175           0.471486                   108                  33
   tracerlpm_parity_hier_oldwater         1248                  0.000000                0.233116           0.121795                    15                 154
          tracerlpm_strict_parity         1249                  0.000000               -0.003016           0.013611                     6                   1
```
