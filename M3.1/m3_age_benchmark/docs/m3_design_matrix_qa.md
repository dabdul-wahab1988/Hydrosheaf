# M3 Phase-2 Design Matrix QA

- Generated: 2026-08-03T08:47:30.019320+00:00
- Config: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3.1\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3.1\m3_age_benchmark\results\m3_design_matrix_results.csv`
- Output rows: 16536
- Unique scenarios: 13
- Row limit: full

## Scenario Metrics

```text
                      scenario_id  total_rows  supported_rows  identifiable_rows  evaluated_rows  nonidentifiable_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
                  ablation_no_he4        1272            1272                 65            1272                  1207         65.0              65.0                0.131041    0.398084         0.692308          1.000000                  0.0 -2.997147
                 ablation_raw_c14        1272            1272                 61            1272                  1211         61.0              61.0                0.146889    0.418957         0.672131          1.000000                  0.0 -3.195868
               ablation_raw_gases        1272               0                 47            1272                  1225          NaN               NaN                     NaN         NaN              NaN               NaN                  NaN       NaN
   hydrosheaf_selection_corrected        1272            1272                309            1272                   963        309.0             309.0                0.130446    0.609790         0.673139          0.915858                 87.0  0.763909
            oldwater_c14_ensemble        1272            1272                 59            1272                  1213         59.0              59.0                0.146889    0.413759         0.677966          1.000000                  0.0 -4.243558
oldwater_ensemble_he4_uncertainty        1272            1272                 59            1272                  1213         59.0              59.0                0.146889    0.413759         0.677966          1.000000                  0.0 -4.243558
         oldwater_he4_uncertainty        1272            1272                 66            1272                  1206         66.0              66.0                0.136067    0.395457         0.696970          1.000000                  1.0 -2.617950
        parity_reported_corrected        1272            1272                 66            1272                  1206         66.0              66.0                0.136067    0.395457         0.696970          1.000000                  1.0 -2.617955
               screened_dgm_gases        1272            1272                 66            1272                  1206         66.0              66.0                0.136067    0.395457         0.696970          1.000000                  1.0 -2.617955
                  tracer_old_only        1272            1272                  0            1272                  1272          NaN               NaN                     NaN         NaN              NaN               NaN                  NaN       NaN
                tracer_young_only        1272            1272                 67            1272                  1205         65.0              65.0                0.192187    0.687181         0.538462          0.923077                  0.0 -0.099892
    tracerlpm_parity_agefractions        1272            1272                289            1272                   983        289.0             289.0                0.021613    0.196357         0.916955          0.996540                166.0  0.969780
          tracerlpm_strict_parity        1272            1272                329            1272                   943        329.0             329.0                0.027932    0.276882         0.875380          0.987842                177.0  0.937147
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
                  ablation_no_he4           65              0.000000e+00            0.000000e+00           0.000000                     0                   0
                 ablation_raw_c14           61              0.000000e+00            5.369846e-03           0.065574                     0                   0
   hydrosheaf_selection_corrected           64             -1.732147e-03           -4.456361e-02           0.531250                    11                   5
            oldwater_c14_ensemble           59              0.000000e+00           -3.481711e-04           0.169492                     0                   0
oldwater_ensemble_he4_uncertainty           59              0.000000e+00           -3.481711e-04           0.169492                     0                   0
         oldwater_he4_uncertainty           66              0.000000e+00           -7.874929e-07           0.015152                     0                   0
               screened_dgm_gases           66              0.000000e+00            0.000000e+00           0.000000                     0                   0
                tracer_young_only           38              5.444561e-08           -5.370691e-03           0.447368                     0                   0
    tracerlpm_parity_agefractions           40              0.000000e+00           -1.198396e-02           0.350000                     0                   1
          tracerlpm_strict_parity           44              0.000000e+00            2.291264e-03           0.181818                     0                   1
```
