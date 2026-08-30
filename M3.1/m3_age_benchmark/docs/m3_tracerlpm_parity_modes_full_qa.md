# M3 Phase-2 Design Matrix QA

> **Historical QA snapshot — superseded 2026-07-28.** Do not use these values
> as current manuscript metrics.

- Generated: 2026-05-18T21:09:35.870676+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_tracerlpm_parity_modes_full.csv`
- Output rows: 5088
- Unique scenarios: 4
- Row limit: full

## Scenario Metrics

```text
                   scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
     parity_reported_corrected        1272         1254              1254                0.221215    0.842500         0.569378          0.838915                  278  0.589063
 tracerlpm_parity_agefractions        1272         1249              1249                0.167356    0.740121         0.614091          0.859087                  278  0.680659
tracerlpm_parity_hier_oldwater        1272         1253              1253                0.362514    1.164932         0.456504          0.750200                  278  0.214909
       tracerlpm_strict_parity        1272         1254              1254                0.220658    0.841340         0.572568          0.838118                  278  0.590195
```

## Paired Effects Versus `parity_reported_corrected`

```text
                   scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
 tracerlpm_parity_agefractions         1247                       0.0               -0.057170           0.491580                    89                  36
tracerlpm_parity_hier_oldwater         1253                       0.0                0.233552           0.122107                    17                 159
       tracerlpm_strict_parity         1254                       0.0               -0.001554           0.013557                     6                   2
```
