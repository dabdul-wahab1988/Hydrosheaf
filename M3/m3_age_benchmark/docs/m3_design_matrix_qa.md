# M3 Phase-2 Design Matrix QA

- Generated: 2026-05-17T11:13:31.923911+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_design_matrix_results.csv`
- Output rows: 880
- Unique scenarios: 11
- Row limit: 80

## Scenario Metrics

```text
                      scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows
                  ablation_no_he4          80           80                80                0.254461    0.819289         0.525000          0.750000                    0
                 ablation_raw_c14          80           80                80                0.333308    0.871348         0.400000          0.675000                   38
               ablation_raw_gases          80           80                80                0.333308    0.871414         0.400000          0.675000                   38
   hydrosheaf_selection_corrected          80           80                80                0.516213    1.343242         0.312500          0.650000                   38
            oldwater_c14_ensemble          80           80                80                0.338094    0.898741         0.400000          0.687500                   38
oldwater_ensemble_he4_uncertainty          80           80                80                0.313141    0.849778         0.475000          0.737500                   38
         oldwater_he4_uncertainty          80           80                80                0.313141    0.814172         0.475000          0.737500                   38
        parity_reported_corrected          80           80                80                0.333308    0.871414         0.400000          0.675000                   38
               screened_dgm_gases          80           80                80                0.333308    0.871414         0.400000          0.675000                   38
                  tracer_old_only          80           79                79                0.933986    1.307564         0.329114          0.493671                   38
                tracer_young_only          80           79                79                1.393501    1.971569         0.126582          0.405063                    0
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
                  ablation_no_he4           74                  0.000000               -0.038148           0.216216                     5                   0
                 ablation_raw_c14           74                  0.000000               -0.001013           0.027027                     0                   0
               ablation_raw_gases           74                  0.000000                0.000000           0.000000                     0                   0
   hydrosheaf_selection_corrected           74                  0.116824                0.330508           0.148649                     1                   9
            oldwater_c14_ensemble           74                  0.000000                0.024174           0.013514                     0                   1
oldwater_ensemble_he4_uncertainty           74                  0.000000               -0.016361           0.202703                     2                   0
         oldwater_he4_uncertainty           74                  0.000000               -0.044615           0.229730                     2                   0
               screened_dgm_gases           74                  0.000000                0.000000           0.000000                     0                   0
                  tracer_old_only           72                  0.000000                0.335447           0.041667                     0                   6
                tracer_young_only           74                  0.123472                0.927734           0.405405                     5                  27
```
