# M3 Phase-2 Design Matrix QA

- Generated: 2026-05-15T07:53:25.595180+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `M3\m3_age_benchmark\results\m3_phase4_younggas_results.csv`
- Output rows: 240
- Unique scenarios: 3
- Row limit: 80

## Scenario Metrics

```text
              scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows
       ablation_raw_gases          80           80                80                0.115949    0.470300           0.7375            0.9500                   21
parity_reported_corrected          80           80                80                0.166393    0.561093           0.6875            0.9250                   21
       screened_dgm_gases          80           80                80                0.115949    0.459719           0.7625            0.9375                   21
```

## Paired Effects Versus `parity_reported_corrected`

```text
       scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
ablation_raw_gases           80                       0.0               -0.073904             0.2250                     6                   2
screened_dgm_gases           80                       0.0               -0.082759             0.1875                     6                   0
```
