# M3 Gas-Correction Audit

- Generated: 2026-07-28T06:41:32.894596+00:00
- Audit output: `C:\Users\ThinkPad P1 G4\Documents\July_2026\NeutroProject\Groundwater\Hydrosheaf\M3\m3_age_benchmark\results\m3_gas_correction_audit.csv`
- Paired rows: 0
- Rows where raw gases improved log10 error: 0
- Rows where raw gases degraded log10 error: 0
- Rows with any changed gas value: 0

## Summary By Age Class

```text
No rows.
```

## Largest Raw-Gas Improvements

```text
Empty DataFrame
Columns: [site_id, age_class, delta_log10_error_raw_minus_corrected, parity_est_age_multi, raw_gas_est_age_multi, ref_age, changed_gas_fields, dissolved_gas_correction]
Index: []
```

## Largest Raw-Gas Degradations

```text
Empty DataFrame
Columns: [site_id, age_class, delta_log10_error_raw_minus_corrected, parity_est_age_multi, raw_gas_est_age_multi, ref_age, changed_gas_fields, dissolved_gas_correction]
Index: []
```

## Interpretation Guardrail

A negative delta means the raw-gas ablation was closer to the USGS reference age than the corrected-gas parity run for the same site. This does not prove raw gases are generally superior; it identifies rows where the dissolved-gas correction pathway, tracer-mode masking, or Hydrosheaf atmospheric-equivalent assumptions need targeted review.

Rows with changed `tritium;he3_trit` should be inspected through `raw_3h3he_apparent_age`, `corrected_3h3he_apparent_age`, and `delta_3h3he_apparent_age_raw_minus_corrected`, because small concentration changes can produce large apparent-age changes when tritium is low.
