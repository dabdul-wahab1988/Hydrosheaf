# M3 Gas-Correction Audit

- Generated: 2026-07-27T13:05:33.397864+00:00
- Audit output: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_gas_correction_audit.csv`
- Paired rows: 1651
- Rows where raw gases improved log10 error: 167
- Rows where raw gases degraded log10 error: 83
- Rows with any changed gas value: 0

## Summary By Age Class

```text
           age_class   n  median_delta_log10_error  mean_delta_log10_error  raw_improved_fraction  raw_degraded_fraction  raw_gained_factor_2  raw_lost_factor_2  any_gas_value_changed_fraction  median_changed_gas_fields
intermediate_50_1000 369                       0.0               -0.038495               0.100271                0.04878                    9                  4                             0.0                        0.0
        modern_le_50 532                       0.0               -0.084860               0.244361                0.12218                   32                  4                             0.0                        0.0
      old_1000_30000 545                       0.0                0.000000               0.000000                0.00000                    0                  0                             0.0                        0.0
   very_old_gt_30000 205                       0.0                0.000000               0.000000                0.00000                    0                  0                             0.0                        0.0
```

## Largest Raw-Gas Improvements

```text
      site_id            age_class  delta_log10_error_raw_minus_corrected  parity_est_age_multi  raw_gas_est_age_multi  ref_age changed_gas_fields dissolved_gas_correction
      WSAC-15         modern_le_50                              -2.357894           1074.558912               4.713420      3.0                                     dgm_sf6
      WSAC-06         modern_le_50                              -2.192151            620.842861               3.988690      1.0                                     dgm_sf6
      ESAC-32         modern_le_50                              -1.817904            857.620710              13.043401     13.0                                     dgm_sf6
   MADCHOW-17         modern_le_50                              -1.800911            928.230447              14.680625     13.0                                     dgm_sf6
  S3-MACK-M05         modern_le_50                              -1.680614           1141.201111              18.522100     21.0                                     dgm_sf6
    CE-QPC-09         modern_le_50                              -1.533866            541.276979              14.211139     15.0                                     dgm_sf6
S4-TUSK-KAW09 intermediate_50_1000                              -1.528527           3880.884951             114.922032     51.0                                     dgm_sf6
  S3-MACK-K29         modern_le_50                              -1.494589           1074.558912              19.647407     26.0                                     dgm_sf6
   MADCHOW-29         modern_le_50                              -1.449680            471.908089              15.277873     16.0                                     dgm_sf6
      NSAC-02         modern_le_50                              -1.448533            575.842697              15.804277     18.0                                     dgm_sf6
```

## Largest Raw-Gas Degradations

```text
    site_id    age_class  delta_log10_error_raw_minus_corrected  parity_est_age_multi  raw_gas_est_age_multi  ref_age changed_gas_fields dissolved_gas_correction
    KERN-36 modern_le_50                               2.787217             19.121667           11715.000006     19.0                                     dgm_sf6
    DM-U-16 modern_le_50                               2.191318             17.537499            2724.490931      9.0                                     dgm_sf6
S3-MACK-K02 modern_le_50                               2.073229             23.968058           11378.333268     48.0                                     dgm_sf6
    RICE-14 modern_le_50                               2.036455             14.711758               0.010000      4.0                                     dgm_sf6
    DM-U-14 modern_le_50                               1.833719            182.627868           12453.333349     40.0                                     dgm_sf6
S3-MACK-K21 modern_le_50                               1.629060             22.479315            2880.091797     39.0                                     dgm_sf6
    KERN-24 modern_le_50                               1.510545            250.000000               0.010000      9.0                                     dgm_sf6
   MERMW-01 modern_le_50                               1.464398             13.729736               0.250000     10.0                                     dgm_sf6
    RICE-19 modern_le_50                               1.434658             14.702864               0.010000      2.0                                     dgm_sf6
  TRLKFP-01 modern_le_50                               1.029216             13.463134               0.250000      6.0                                     dgm_sf6
```

## Interpretation Guardrail

A negative delta means the raw-gas ablation was closer to the USGS reference age than the corrected-gas parity run for the same site. This does not prove raw gases are generally superior; it identifies rows where the dissolved-gas correction pathway, tracer-mode masking, or Hydrosheaf atmospheric-equivalent assumptions need targeted review.

Rows with changed `tritium;he3_trit` should be inspected through `raw_3h3he_apparent_age`, `corrected_3h3he_apparent_age`, and `delta_3h3he_apparent_age_raw_minus_corrected`, because small concentration changes can produce large apparent-age changes when tritium is low.
