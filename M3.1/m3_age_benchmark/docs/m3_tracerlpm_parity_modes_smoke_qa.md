# M3 Phase-2 Design Matrix QA

> **Historical smoke-test snapshot — superseded 2026-07-28.** Do not use its
> values as current manuscript evidence.

- Generated: 2026-05-18T18:21:55.224054+00:00
- Config: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `M3\m3_age_benchmark\results\m3_tracerlpm_parity_modes_smoke.csv`
- Output rows: 240
- Unique scenarios: 3
- Row limit: 80

## Scenario Metrics

```text
                   scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
 tracerlpm_parity_agefractions          80           79                79                0.122069    0.625030         0.746835          0.936709                   20  0.815752
tracerlpm_parity_hier_oldwater          80           79                79                0.148508    0.639449         0.759494          0.924051                   20  0.807153
       tracerlpm_strict_parity          80           79                79                0.150025    0.636106         0.759494          0.924051                   20  0.809164
```

## Paired Effects Versus `parity_reported_corrected`

- No paired effect rows.
