# M3 Young-Tracer Diagnostic and Decision

Date: 2026-05-15

## Context

The previous M3 gas-correction audit showed that raw gas tracers occasionally resulted in better multi-tracer log10 errors than the USGS dissolved-gas-model (DGM) corrected concentrations (Table 6 `TRC_*` fields). To understand this, we conducted a targeted diagnostic comparing the simple Piston Flow Model (PFM) apparent ages of raw versus corrected `3H/3He` and `SF6` against the USGS reference ages.

## Diagnostic Results

We ran a PFM apparent age script (`diagnostic_young_tracers.py`) over the USGS national dataset, looking strictly at samples where both raw and corrected concentrations were available and differed.

**3H/3He Apparent Age:**
*   **Rows where raw and corrected ages differ:** 387
*   **Raw 3H/3He error better than corrected:** 212 out of 387 (54.8%)
*   **Median log10 error Raw:** 0.430
*   **Median log10 error Corrected:** 0.436

**SF6 Apparent Age:**
*   **Rows where raw and corrected ages differ:** 432
*   **Raw SF6 error better than corrected:** 222 out of 432 (51.4%)
*   **Median log10 error Raw:** 2.394
*   **Median log10 error Corrected:** 2.402
*   *(Note: The high log10 error for SF6 reflects that SF6 is a young tracer being compared to a dataset with many old/mixed reference ages, but the relative performance holds).*

## Interpretation

The diagnostic confirms the signal seen in the broader multi-tracer audit: **USGS DGM-corrected concentrations do not uniformly improve single-tracer PFM apparent age accuracy relative to raw measured concentrations.** In fact, in roughly half the cases where they differ, the raw concentration yields a PFM age closer to the final USGS reference age.

Why does this happen?
1.  **DGM Uncertainty:** Dissolved-gas correction depends on Ne, Ar, N2, temperature, and specific excess-air models (e.g., CE, UA). Errors or assumptions in these parameters can occasionally "over-correct" a sample.
2.  **Complex Mixing:** The USGS reference ages are often derived from complex lumped parameter models (e.g., Exponential Mixing, Dispersion) fitting multiple tracers simultaneously, not from single-tracer PFM. The "best fit" in a complex model space might not map linearly to PFM apparent age.
3.  **Contamination/Degradation:** Corrected values might reflect physical processes (like degassing or contamination) that make them *less* representative of pure advective age, but *more* physically accurate to the sample state.

## Decision: How to Treat Table 6 `TRC_*` Values

Given the M3 Experimental Design Plan's goal of splitting the benchmark into "USGS-parity mode" and "Hydrosheaf-selection mode", we decide the following:

1.  **Treat as Reported Model Inputs for Parity:** The `TRC_*` values in USGS Table 6 (e.g., `TRC_3H_TU`, `TRC_SF6_pptv`) are the *exact* inputs USGS used in their TracerLPM models to derive the reference ages. Therefore, for the **USGS-parity benchmark** (Phase 1), they MUST be used. We cannot replicate USGS estimates without using their chosen inputs.
2.  **Treat as Diagnostic-Only for Hydrosheaf Capability:** For the **Hydrosheaf-selection benchmark** (Phase 4 ablation experiments), the `TRC_*` values should not be treated as absolute ground-truth "corrected concentrations." Instead, the benchmark must explicitly run ablation scenarios (Raw vs. USGS-Corrected) and ultimately support Hydrosheaf reconstructing its own DGM corrections from primary noble gas data.
3.  **Integration into the Design Matrix:** The `m3_design_matrix.yaml` must include `gas_correction: [raw, usgs_reported]` as a primary experimental factor to quantify this effect across all age classes.

This decision aligns directly with the "Critical Flaws" identified in the M3 plan, explicitly quantifying the effect of dissolved-gas correction rather than assuming it is universally beneficial.
