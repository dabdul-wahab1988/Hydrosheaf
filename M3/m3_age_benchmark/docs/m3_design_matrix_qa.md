# M3 Phase-2 Design Matrix QA

- Generated: 2026-07-27T13:04:48.841706+00:00
- Config: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_design_matrix_results.csv`
- Output rows: 23030
- Unique scenarios: 14
- Row limit: full

## Scenario Metrics

```text
                      scenario_id  total_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
                  ablation_no_he4        1645         1577              1577                0.287742    1.028120         0.506024          0.799620                    0  0.351606
                 ablation_raw_c14        1645         1608              1608                0.284681    1.049055         0.511194          0.788557                  276  0.370504
               ablation_raw_gases        1645         1609              1609                0.274908    1.040379         0.520199          0.812927                  278  0.379897
   hydrosheaf_selection_corrected        1645          666               666                0.469818    0.912343         0.366366          0.768769                  127  0.240949
            oldwater_c14_ensemble        1645         1609              1609                0.447158    1.297217         0.402113          0.717837                  278  0.037528
oldwater_ensemble_he4_uncertainty        1645         1626              1626                0.411245    1.273281         0.418819          0.730627                  295  0.069368
         oldwater_he4_uncertainty        1645         1627              1627                0.277442    1.026365         0.516288          0.802704                  295  0.394941
        parity_reported_corrected        1645         1610              1610                0.301612    1.069606         0.499379          0.787578                  278  0.345245
               screened_dgm_gases        1645         1610              1610                0.301612    1.069606         0.499379          0.787578                  278  0.345245
                  tracer_old_only        1645         1198              1198                0.831789    1.601935         0.338063          0.533389                  735 -0.488090
                tracer_young_only        1645         1533              1533                0.986894    1.619777         0.228311          0.506849                    0 -0.516814
    tracerlpm_parity_agefractions        1645         1610              1610                0.239417    1.006682         0.547205          0.796894                  278  0.419479
   tracerlpm_parity_hier_oldwater        1645         1611              1611                0.431364    1.283164         0.414029          0.718187                  278  0.058656
          tracerlpm_strict_parity        1645         1612              1612                0.291562    1.067338         0.503102          0.789082                  278  0.348288
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
                  ablation_no_he4         1562                  0.000000               -0.036317           0.093470                    30                  22
                 ablation_raw_c14         1610                  0.000000               -0.008702           0.134161                    28                  10
               ablation_raw_gases         1611                  0.000000               -0.034543           0.101179                    41                   8
   hydrosheaf_selection_corrected          662                  0.000000               -0.108527           0.495468                    93                  91
            oldwater_c14_ensemble         1611                  0.000000                0.191247           0.068281                     5                 162
oldwater_ensemble_he4_uncertainty         1611                  0.000000                0.163846           0.145251                    38                 173
         oldwater_he4_uncertainty         1612                  0.000000               -0.036012           0.088089                    35                  11
               screened_dgm_gases         1612                  0.000000                0.000000           0.000620                     0                   0
                  tracer_old_only         1196                  0.089709                0.592834           0.162207                    37                 280
                tracer_young_only         1516                  0.016059                0.579897           0.299472                   109                 488
    tracerlpm_parity_agefractions         1608                  0.000000               -0.047365           0.365672                   108                  33
   tracerlpm_parity_hier_oldwater         1611                  0.000000                0.180589           0.094972                    15                 154
          tracerlpm_strict_parity         1612                  0.000000               -0.002337           0.011166                     6                   1
```
