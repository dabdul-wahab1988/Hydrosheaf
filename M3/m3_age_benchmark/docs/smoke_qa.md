# M3 Phase-2 Design Matrix QA

- Generated: 2026-07-27T12:13:52.223477+00:00
- Config: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `C:\Users\THINKP~1\AppData\Local\Temp\claude\C--Users-ThinkPad-P1-G4-Documents-July-2026-NeutroProject-Groundwater-Hydrosheaf\8c2a7981-2fc7-4157-8ce8-003af41dcafa\scratchpad\smoke.csv`
- Output rows: 80
- Unique scenarios: 2
- Row limit: 40

## Scenario Metrics

```text
              scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
parity_reported_corrected          40           38                38                0.346022    1.249037         0.447368          0.789474                    4   0.15899
  tracerlpm_strict_parity          40           38                38                0.346022    1.249037         0.447368          0.789474                    4   0.15899
```

## Paired Effects Versus `parity_reported_corrected`

```text
            scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
tracerlpm_strict_parity           38                       0.0                     0.0                0.0                     0                   0
```
