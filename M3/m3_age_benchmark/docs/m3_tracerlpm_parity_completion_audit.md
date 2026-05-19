# M3 TracerLPM-Parity Completion Audit

Generated: 2026-05-18

## Implementation Status

- Strict TracerLPM parity is implemented with reported LPM family, reported tracer mask, WISER nearest tritium histories, site-adjusted gas histories, UZ travel-time handling, and preserved USGS source columns.
- Gas observations now carry censored or contaminated-mixture likelihoods where appropriate.
- Hierarchical old-water priors are implemented with shrinkage/fallback behavior, but the full benchmark shows this mode is not currently suitable as the main validation result.
- Age-fraction constraints are implemented and produce the strongest independent full benchmark result.
- USGS-calibrated benchmark emulation is implemented with leave-study-unit-out folds and is labelled separately from independent validation.

## Full Benchmark Metrics

| Mode | N finite | Median abs log10 error | log10 RMSE | log10 R2 | Within factor 2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Reported-model parity | 1254 | 0.2212 | 0.8425 | 0.5891 | 0.5694 |
| Strict TracerLPM parity | 1254 | 0.2207 | 0.8413 | 0.5902 | 0.5726 |
| Strict parity + age fractions | 1249 | 0.1674 | 0.7401 | 0.6807 | 0.6141 |
| Strict parity + hierarchical old-water priors | 1253 | 0.3625 | 1.1649 | 0.2149 | 0.4565 |
| USGS-calibrated benchmark emulation | 675 | 0.1718 | 0.2752 | 0.9618 | 0.7689 |

## Recommendation

Use `tracerlpm_parity_agefractions` as the strongest independent parity benchmark if the manuscript can justify using USGS age fractions as reported TTD-shape constraints. Use `tracerlpm_strict_parity` as the conservative Fig. 5 mode if the figure must only compare scalar reported ages. Do not use hierarchical old-water priors as the headline result until the priors are recalibrated by study unit/aquifer with stronger geochemical support.

Label `m3_usgs_calibrated_parity.csv` strictly as USGS-calibrated benchmark emulation, not independent validation.
