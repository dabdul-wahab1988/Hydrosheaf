# M3 correctness rerun — audit-only pre-reportability metrics

> **Not a manuscript performance table.** The values below include raw scalar
> fits before the final reference-free reportability guard. They are retained to
> diagnose catastrophic branches and support the withdrawal of the hierarchical
> old-water prior. Publication metrics must be taken from the canonical
> `m3_design_matrix_*` outputs and `Manuscript_Ready` tables.

- Generated: 2026-07-28T02:53:16.180622+00:00
- Config: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\configs\design_matrix.yaml`
- Output: `M3\m3_age_benchmark\results\m3_correctness_rerun_full.csv`
- Output rows: 17808
- Unique scenarios: 14
- Row limit: full

## Scenario Metrics

```text
                      scenario_id  total_rows  supported_rows  metric_rows  finite_estimates  median_abs_log10_error  log10_rmse  within_factor_2  within_factor_10  calibrated_he4_rows  log10_r2
                  ablation_no_he4        1272            1272       1215.0            1215.0                0.206210    0.845809         0.592593          0.850206                  0.0  0.555818
                 ablation_raw_c14        1272            1272       1266.0            1266.0                0.133680    0.730755         0.664297          0.876777                295.0  0.687205
               ablation_raw_gases        1272               0          NaN               NaN                     NaN         NaN              NaN               NaN                  NaN       NaN
   hydrosheaf_selection_corrected        1272            1272       1254.0            1254.0                0.162733    0.763094         0.645933          0.871611                295.0  0.654424
            oldwater_c14_ensemble        1272            1272       1266.0            1266.0                0.261364    1.185134         0.530016          0.793049                295.0  0.178402
oldwater_ensemble_he4_uncertainty        1272            1272       1266.0            1266.0                0.294349    1.221662         0.503949          0.783570                295.0  0.126976
         oldwater_he4_uncertainty        1272            1272       1266.0            1266.0                0.173371    0.777659         0.638231          0.875197                295.0  0.645763
        parity_reported_corrected        1272            1272       1266.0            1266.0                0.142236    0.749335         0.662717          0.878357                295.0  0.671096
               screened_dgm_gases        1272            1272       1266.0            1266.0                0.142236    0.749335         0.662717          0.878357                295.0  0.671096
                  tracer_old_only        1272            1272       1184.0            1184.0                0.369749    1.253055         0.462838          0.647804                788.0  0.073116
                tracer_young_only        1272            1272       1174.0            1174.0                1.270030    1.638522         0.207836          0.427598                  0.0 -0.544965
    tracerlpm_parity_agefractions        1272            1272       1267.0            1267.0                0.093643    0.653677         0.710339          0.891871                295.0  0.749534
   tracerlpm_parity_hier_oldwater        1272            1272       1256.0            1256.0                0.362430    1.309757         0.464968          0.751592                283.0  0.004418
          tracerlpm_strict_parity        1272            1272       1268.0            1268.0                0.141210    0.745886         0.666404          0.879338                295.0  0.674602
```

## Paired Effects Versus `parity_reported_corrected`

```text
                      scenario_id  paired_rows  median_delta_log10_error  mean_delta_log10_error  improved_fraction  gained_factor_2_rows  lost_factor_2_rows
                  ablation_no_he4         1215                  0.000000                0.081610           0.041975                     3                  93
                 ablation_raw_c14         1266                  0.000000               -0.004627           0.181675                    20                  18
   hydrosheaf_selection_corrected         1248                  0.000555                0.030543           0.401442                   104                 136
            oldwater_c14_ensemble         1265                  0.000000                0.256969           0.091700                     4                 173
oldwater_ensemble_he4_uncertainty         1265                  0.000000                0.295462           0.119368                     4                 206
         oldwater_he4_uncertainty         1266                  0.000000                0.025434           0.035545                     0                  31
               screened_dgm_gases         1266                  0.000000                0.000000           0.000000                     0                   0
                  tracer_old_only         1179                  0.001641                0.416436           0.192536                    42                 274
                tracer_young_only         1172                  0.697097                0.904503           0.192833                    52                 573
    tracerlpm_parity_agefractions         1264                  0.000000               -0.059233           0.481804                    84                  26
   tracerlpm_parity_hier_oldwater         1253                  0.000000                0.378102           0.150040                    14                 262
          tracerlpm_strict_parity         1266                  0.000000               -0.002728           0.015008                     5                   1
```
