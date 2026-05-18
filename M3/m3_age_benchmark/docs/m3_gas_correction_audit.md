# M3 Gas-Correction Audit

- Generated: 2026-05-17T08:10:12.628447+00:00
- Audit output: `C:\Users\DicksonAbdul-Wahab\Desktop\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_gas_correction_audit.csv`
- Paired rows: 80
- Rows where raw gases improved log10 error: 0
- Rows where raw gases degraded log10 error: 0
- Rows with any changed gas value: 0

## Summary By Age Class

```text
           age_class  n  median_delta_log10_error  mean_delta_log10_error  raw_improved_fraction  raw_degraded_fraction  raw_gained_factor_2  raw_lost_factor_2  any_gas_value_changed_fraction  median_changed_gas_fields
intermediate_50_1000 20                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
        modern_le_50 20                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
      old_1000_30000 20                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
   very_old_gt_30000 20                       0.0                     0.0                    0.0                    0.0                    0                  0                             0.0                        0.0
```

## Largest Raw-Gas Improvements

```text
    site_id            age_class  delta_log10_error_raw_minus_corrected  parity_est_age_multi  raw_gas_est_age_multi  ref_age changed_gas_fields dissolved_gas_correction
BISCPAS1-04 intermediate_50_1000                                    0.0            262.817614             262.817614    154.0                                     dgm_sf6
BISCPAS1-25 intermediate_50_1000                                    0.0              0.250000               0.250000     52.6                                     dgm_sf6
BNRCPAS1-03 intermediate_50_1000                                    0.0            205.882353             205.882353     85.4                                     dgm_sf6
BNRCPAS1-09 intermediate_50_1000                                    0.0           2259.228660            2259.228660    102.0                                     dgm_sf6
BNRCPAS1-11 intermediate_50_1000                                    0.0            118.263269             118.263269     76.9                                     dgm_sf6
BNRCPAS1-12 intermediate_50_1000                                    0.0             54.529274              54.529274     65.3                                     dgm_sf6
BNRCPAS1-19 intermediate_50_1000                                    0.0             18.573415              18.573415     66.4                                     dgm_sf6
BNRFPAS1-02 intermediate_50_1000                                    0.0           5589.518805            5589.518805    747.0                                     dgm_sf6
BNRFPAS1-04 intermediate_50_1000                                    0.0             26.595271              26.595271     60.8                                     dgm_sf6
BNRFPAS1-07 intermediate_50_1000                                    0.0             38.081765              38.081765     52.5                                     dgm_sf6
```

## Largest Raw-Gas Degradations

```text
    site_id            age_class  delta_log10_error_raw_minus_corrected  parity_est_age_multi  raw_gas_est_age_multi  ref_age changed_gas_fields dissolved_gas_correction
BISCPAS1-04 intermediate_50_1000                                    0.0            262.817614             262.817614    154.0                                     dgm_sf6
BISCPAS1-25 intermediate_50_1000                                    0.0              0.250000               0.250000     52.6                                     dgm_sf6
BNRCPAS1-03 intermediate_50_1000                                    0.0            205.882353             205.882353     85.4                                     dgm_sf6
BNRCPAS1-09 intermediate_50_1000                                    0.0           2259.228660            2259.228660    102.0                                     dgm_sf6
BNRCPAS1-11 intermediate_50_1000                                    0.0            118.263269             118.263269     76.9                                     dgm_sf6
BNRCPAS1-12 intermediate_50_1000                                    0.0             54.529274              54.529274     65.3                                     dgm_sf6
BNRCPAS1-19 intermediate_50_1000                                    0.0             18.573415              18.573415     66.4                                     dgm_sf6
BNRFPAS1-02 intermediate_50_1000                                    0.0           5589.518805            5589.518805    747.0                                     dgm_sf6
BNRFPAS1-04 intermediate_50_1000                                    0.0             26.595271              26.595271     60.8                                     dgm_sf6
BNRFPAS1-07 intermediate_50_1000                                    0.0             38.081765              38.081765     52.5                                     dgm_sf6
```

## Interpretation Guardrail

A negative delta means the raw-gas ablation was closer to the USGS reference age than the corrected-gas parity run for the same site. This does not prove raw gases are generally superior; it identifies rows where the dissolved-gas correction pathway, tracer-mode masking, or Hydrosheaf atmospheric-equivalent assumptions need targeted review.

Rows with changed `tritium;he3_trit` should be inspected through `raw_3h3he_apparent_age`, `corrected_3h3he_apparent_age`, and `delta_3h3he_apparent_age_raw_minus_corrected`, because small concentration changes can produce large apparent-age changes when tritium is low.
