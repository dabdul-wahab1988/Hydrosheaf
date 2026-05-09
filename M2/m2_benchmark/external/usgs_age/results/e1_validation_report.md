# E1 USGS Public Tracer-Age Validation Report

Run timestamp UTC: 2026-05-04T16:15:49.120297+00:00

## Status

Completed from the local USGS release files supplied for DOI `10.5066/P9W7T0DN`.

## Source Tables Used

- Age/reference table: `Table_2_Ages.txt`
- Tracer table: `Table_3_Tracers.txt`
- Site table: `Table_1_Sites.txt`

## Outputs

- Validation rows: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\usgs_age\results\usgs_age_validation.csv`
- Summary metrics: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\external\usgs_age\results\usgs_age_validation_summary.csv`
- Figure 4A: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M2\m2_benchmark\figures\figure4a_public_tracer_age_agreement.png`

## Sample Counts By Hydrosheaf Tracer Method

- `3H+3H/3He+14C+SF6`: 268 samples
- `14C`: 228 samples
- `3H+14C+SF6`: 64 samples
- `3H+14C`: 47 samples
- `3H+3H/3He+14C`: 22 samples
- `3H+3H/3He+SF6`: 17 samples
- `3H+3H/3He+14C+CFC11+CFC12+CFC113`: 10 samples
- `14C+SF6`: 6 samples
- `SF6`: 5 samples
- `3H+SF6`: 5 samples
- `3H+3H/3He+14C+CFC11+CFC12`: 5 samples
- `3H`: 4 samples
- `3H+3H/3He+14C+CFC11`: 3 samples
- `3H+14C+CFC11+CFC12+CFC113`: 3 samples
- `14C+CFC11`: 1 samples
- `3H+3H/3He`: 1 samples

## Metrics

| status    | group                            |   n_samples |   log10_age_rmse |   median_log10_bias |   median_abs_log10_error |   ci_coverage_fraction | age_table        | tracer_table        | note                                                                                                                       |
|:----------|:---------------------------------|------------:|-----------------:|--------------------:|-------------------------:|-----------------------:|:-----------------|:--------------------|:---------------------------------------------------------------------------------------------------------------------------|
| completed | all                              |         689 |         1.23003  |          -0.140846  |                0.30656   |               0.801161 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 14C                              |         228 |         0.396124 |          -0.0495744 |                0.0599249 |               0.657895 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 14C+CFC11                        |           1 |         0.952621 |          -0.952621  |                0.952621  |               1        | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 14C+SF6                          |           6 |         0.784833 |          -0.756571  |                0.756571  |               1        | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H                               |           4 |         0.590998 |          -0.505869  |                0.512872  |               0.25     | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+14C                           |          47 |         2.32681  |          -2.16802   |                2.16802   |               0.829787 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+14C+CFC11+CFC12+CFC113        |           3 |         0.216078 |          -0.0909438 |                0.10991   |               1        | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+14C+SF6                       |          64 |         1.07139  |          -0.393641  |                0.651001  |               0.96875  | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He                        |           1 |         0.332253 |          -0.332253  |                0.332253  |               1        | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He+14C                    |          22 |         1.1228   |          -0.327463  |                0.403742  |               0.818182 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He+14C+CFC11              |           3 |         1.40088  |          -0.20796   |                1.05795   |               0.666667 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He+14C+CFC11+CFC12        |           5 |         2.05999  |          -2.20947   |                2.20947   |               0.6      | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He+14C+CFC11+CFC12+CFC113 |          10 |         1.50713  |          -0.38755   |                0.38755   |               0.9      | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He+14C+SF6                |         268 |         1.4343   |          -0.317231  |                0.514099  |               0.914179 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+3H/3He+SF6                    |          17 |         0.676641 |          -0.44183   |                0.44183   |               0.470588 | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | 3H+SF6                           |           5 |         0.762972 |          -0.814485  |                0.814485  |               0.4      | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |
| completed | SF6                              |           5 |         2.56561  |          -3.14947   |                3.14947   |               0.4      | Table_2_Ages.txt | Table_3_Tracers.txt | Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available. |

## Interpretation Guardrail

This E1 run is a public-data agreement test for Hydrosheaf's upgraded single-node multi-tracer age inference. The implementation now uses 3H/3He closed-system ingrowth and SF6/CFC atmospheric-equivalent screening ages in addition to 3H and 14C. The USGS reference ages were fitted with the full TracerLPM workflow and site-specific corrections, so remaining disagreement should be interpreted as evidence for future work on local recharge histories, dissolved-gas corrections, and mixture-model fitting.
