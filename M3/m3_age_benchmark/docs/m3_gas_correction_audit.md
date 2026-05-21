# M3 Gas-Correction Audit

- Generated: 2026-05-21T16:22:59.114834+00:00
- Audit output: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_gas_correction_audit.csv`
- Paired rows: 1272
- Rows where raw gases improved log10 error: 0
- Rows where raw gases degraded log10 error: 0
- Rows with any changed gas value: 0

## Summary By Age Class

```text
           age_class   n  median_delta_log10_error  mean_delta_log10_error  raw_improved_fraction  raw_degraded_fraction  raw_gained_factor_2  raw_lost_factor_2  any_gas_value_changed_fraction  median_changed_gas_fields
intermediate_50_1000 252                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
        modern_le_50 302                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
      old_1000_30000 517                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
   very_old_gt_30000 201                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
```

## Largest Raw-Gas Improvements

```text
    site_id            age_class  delta_log10_error_raw_minus_corrected  parity_est_age_multi  raw_gas_est_age_multi  ref_age changed_gas_fields dissolved_gas_correction
RGAQPAS1-49       old_1000_30000                                    0.0           4204.915178            4204.915178   2990.0                                     dgm_sf6
RGAQPAS1-48    very_old_gt_30000                                    0.0          50000.000000           50000.000000  47000.0                                     dgm_sf6
RGAQPAS1-45       old_1000_30000                                    0.0              0.250000               0.250000   2030.0                                     dgm_sf6
RGAQPAS1-41       old_1000_30000                                    0.0          17827.859405           17827.859405  22300.0                                     dgm_sf6
RGAQPAS1-37    very_old_gt_30000                                    0.0          29360.903933           29360.903933  30900.0                                     dgm_sf6
RGAQPAS1-32       old_1000_30000                                    0.0          10414.486120           10414.486120  13900.0                                     dgm_sf6
RGAQPAS1-31       old_1000_30000                                    0.0           6031.428848            6031.428848   5000.0                                     dgm_sf6
RGAQPAS1-30 intermediate_50_1000                                    0.0              0.250000               0.250000    346.0                                     dgm_sf6
RGAQPAS1-29       old_1000_30000                                    0.0          15440.656077           15440.656077  19800.0                                     dgm_sf6
RGAQPAS1-26 intermediate_50_1000                                    0.0              6.326359               6.326359    432.0                                     dgm_sf6
```

## Largest Raw-Gas Degradations

```text
    site_id         age_class  delta_log10_error_raw_minus_corrected  parity_est_age_multi  raw_gas_est_age_multi  ref_age changed_gas_fields dissolved_gas_correction
RGAQPAS1-59 very_old_gt_30000                                    0.0          47285.596925           47285.596925  30500.0                                     dgm_sf6
BISCPAS1-01      modern_le_50                                    0.0             11.079666              11.079666     25.6                                     dgm_sf6
BISCPAS1-02      modern_le_50                                    0.0             16.676512              16.676512     19.3                                     dgm_sf6
BISCPAS1-03      modern_le_50                                    0.0             16.852728              16.852728     21.4                                     dgm_sf6
BISCPAS1-05      modern_le_50                                    0.0             14.346055              14.346055     10.1                                     dgm_sf6
BISCPAS1-06      modern_le_50                                    0.0             63.279850              63.279850     24.0                                     dgm_sf6
RGAQPAS1-20    old_1000_30000                                    0.0           3914.457733            3914.457733   3700.0                                     dgm_sf6
RGAQPAS1-14    old_1000_30000                                    0.0          50000.000000           50000.000000  28100.0                                     dgm_sf6
RGAQPAS1-15    old_1000_30000                                    0.0          26956.666676           26956.666676  19400.0                                     dgm_sf6
RGAQPAS1-13    old_1000_30000                                    0.0              0.250000               0.250000   2130.0                                     dgm_sf6
```

## Interpretation Guardrail

A negative delta means the raw-gas ablation was closer to the USGS reference age than the corrected-gas parity run for the same site. This does not prove raw gases are generally superior; it identifies rows where the dissolved-gas correction pathway, tracer-mode masking, or Hydrosheaf atmospheric-equivalent assumptions need targeted review.

Rows with changed `tritium;he3_trit` should be inspected through `raw_3h3he_apparent_age`, `corrected_3h3he_apparent_age`, and `delta_3h3he_apparent_age_raw_minus_corrected`, because small concentration changes can produce large apparent-age changes when tritium is low.
