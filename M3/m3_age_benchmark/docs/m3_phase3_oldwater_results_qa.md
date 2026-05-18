# M3 Phase-2 Design Matrix QA

- Generated: 2026-05-15T06:20:11.040903+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `M3\m3_age_benchmark\results\m3_phase3_oldwater_results.csv`
- Output rows: 320
- Unique scenarios: 4
- Row limit: 80

## Scenario Metrics

```text
                      scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows
            oldwater_c14_ensemble          80           80                80                0.169934    0.618167           0.6750            0.9125                   21
oldwater_ensemble_he4_uncertainty          80           80                80                0.169934    0.618167           0.6750            0.9125                   21
         oldwater_he4_uncertainty          80           80                80                0.166393    0.561093           0.6875            0.9250                   21
        parity_reported_corrected          80           80                80                0.166393    0.561093           0.6875            0.9250                   21
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
            oldwater_c14_ensemble           80                       0.0                0.035245              0.025                     1                   2
oldwater_ensemble_he4_uncertainty           80                       0.0                0.035245              0.025                     1                   2
         oldwater_he4_uncertainty           80                       0.0                0.000000              0.000                     0                   0
```
