# M3 Phase-2 Design Matrix QA

> **Historical QA snapshot — superseded 2026-07-28.** This targeted run is not
> current manuscript evidence.

- Generated: 2026-05-18T07:26:56.165197+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `M3\m3_age_benchmark\results\m3_selection_paircheck_fixed.csv`
- Output rows: 160
- Unique scenarios: 2
- Row limit: 80

## Scenario Metrics

```text
                   scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows
hydrosheaf_selection_corrected          80           80                80                0.754107    1.599396           0.1500            0.5875                   38
     parity_reported_corrected          80           80                80                0.355212    0.939030           0.3875            0.6625                   38
```

## Paired Effects Versus `parity_reported_corrected`

```text
                   scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
hydrosheaf_selection_corrected           75                  0.371318                0.532578           0.146667                     4                  23
```
