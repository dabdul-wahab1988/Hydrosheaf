# M3 Phase-2 Design Matrix QA

- Generated: 2026-05-21T16:22:47.190236+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_design_matrix_results.csv`
- Output rows: 17808
- Unique scenarios: 14
- Row limit: full

## Scenario Metrics

```text
                      scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
                  ablation_no_he4        1272         1221              1221                0.190602    0.780718         0.576577          0.853399                    0  0.623836
                 ablation_raw_c14        1272         1252              1252                0.209045    0.800052         0.575080          0.842652                  276  0.629836
               ablation_raw_gases        1272         1254              1254                0.209559    0.824424         0.574163          0.842105                  278  0.606508
   hydrosheaf_selection_corrected        1272          493               493                0.487361    0.980279         0.385396          0.699797                  129  0.231146
            oldwater_c14_ensemble        1272         1253              1253                0.381483    1.166627         0.442139          0.751796                  278  0.212623
oldwater_ensemble_he4_uncertainty        1272         1270              1270                0.353535    1.132377         0.458268          0.767717                  295  0.253186
         oldwater_he4_uncertainty        1272         1271              1271                0.190429    0.753343         0.590873          0.859953                  295  0.669230
        parity_reported_corrected        1272         1254              1254                0.209559    0.824424         0.574163          0.842105                  278  0.606508
               screened_dgm_gases        1272         1254              1254                0.209559    0.824424         0.574163          0.842105                  278  0.606508
                  tracer_old_only        1272         1131              1131                0.765657    1.593109         0.351017          0.542882                  735 -0.489382
                tracer_young_only        1272         1199              1199                1.239441    1.611130         0.192661          0.425354                    0 -0.499535
    tracerlpm_parity_agefractions        1272         1249              1249                0.169311    0.736201         0.612490          0.851882                  278  0.684033
   tracerlpm_parity_hier_oldwater        1272         1253              1253                0.359651    1.154766         0.456504          0.751796                  278  0.228552
          tracerlpm_strict_parity        1272         1254              1254                0.207998    0.823234         0.575758          0.842105                  278  0.607643
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
                  ablation_no_he4         1204                  0.000000               -0.029908           0.116279                    29                  31
                 ablation_raw_c14         1252                  0.000000               -0.002394           0.181310                    17                  17
               ablation_raw_gases         1254                  0.000000                0.000000           0.000000                     0                   0
   hydrosheaf_selection_corrected          493                  0.010162                0.078355           0.464503                    56                  71
            oldwater_c14_ensemble         1253                  0.000000                0.247222           0.078212                     3                 169
oldwater_ensemble_he4_uncertainty         1253                  0.000000                0.213994           0.182761                    34                 184
         oldwater_he4_uncertainty         1254                  0.000000               -0.043882           0.116427                    34                  15
               screened_dgm_gases         1254                  0.000000                0.000000           0.000000                     0                   0
                  tracer_old_only         1131                  0.109854                0.659022           0.131742                    21                 282
                tracer_young_only         1184                  0.510509                0.804213           0.251689                    85                 514
    tracerlpm_parity_agefractions         1247                  0.000000               -0.048056           0.483561                    86                  41
   tracerlpm_parity_hier_oldwater         1253                  0.000000                0.235250           0.113328                    13                 161
          tracerlpm_strict_parity         1254                  0.000000               -0.001561           0.013557                     4                   2
```
