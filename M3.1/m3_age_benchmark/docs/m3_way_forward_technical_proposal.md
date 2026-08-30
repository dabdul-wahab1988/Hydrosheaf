# Hydrosheaf USGS Benchmarking: Way Forward Technical Proposal

> **Historical proposal — superseded 2026-07-28.** Proposed mechanisms and
> expected improvements below are not achieved results. Use the locked current
> evidence and decisions instead.

- **Date:** 2026-05-15
- **Subject:** Roadmap to achieve strong improvement in M3 canonical age benchmarking.

## 1. Executive Summary

Hydrosheaf's refreshed benchmarking against the USGS national dataset showed that while young-gas screening is useful for identifying contamination (especially SF6), the canonical full age-grid-35 run remains only a screening-level public-age result: 1,272 eligible rows, median absolute log10 error 0.383, log10 RMSE 1.065, and 41.0% within a factor of 2. The remaining weakness is now dominated by structural limitations in tracer weighting, dissolved-gas handling, and old-groundwater uncertainty rather than by a missing full run.

The "Way Forward" focuses on moving from binary masking to probabilistic reliability weighting and bringing the dissolved-gas physical model inside Hydrosheaf to reduce dependence on external corrections.

---

## 2. Phase 1: Performance Baseline Correction (Immediate)

The investigation revealed that earlier full-screened runs used a coarse age grid (`age_steps = 8`) which made those outputs unsuitable as canonical manuscript evidence. The corrected full screened run has now been executed at `age_steps = 35`.

- **Current baseline:** 1,272 eligible rows, median absolute log10 error 0.383, log10 RMSE 1.065, within factor 2 = 0.410.
- **Goal:** Improve the canonical baseline through reliability weighting and stronger old-water handling.
- **Verification:** Confirm improvements against this refreshed full-run baseline, not against the older coarse-grid files.

---

## 3. Phase 2: Soft Tracer Reliability Weighting (Structural)

Current behavior uses a "hard" choice between raw and USGS-corrected values, or binary masking of contaminated tracers. This is brittle for borderline samples.

- **Redesign:** Implement a `TracerReliabilityScorer` that calculates a reliability index ($R \in [0, 1]$) for each tracer.
- **Reliability Factors:**
    - **Historical Plausibility:** How far above the maximum historical atmospheric concentration is the sample?
    - **Proxy Coherence:** Does the SF6 age agree with the 3H/3He age? (calculated during pre-processing).
    - **Uncertainty Inflation:** Instead of masking a tracer ($R=0$), inflate its $\sigma$ by $1/R$.
- **Objective Function Update:**
    $$ \chi^2 = \sum_{i} \left( \frac{O_i - P_i}{\sigma_i \cdot \text{Inflation}_i} \right)^2 $$
- **Benefit:** Allows Hydrosheaf to "listen less" to noisy or suspicious tracers without discarding them entirely, improving performance in the "intermediate" age class.

---

## 4. Phase 3: Internal Dissolved-Gas Optimization (Process Model)

Hydrosheaf currently treats USGS gas corrections as static inputs. This ignores the uncertainty in the recharge temperature (T), pressure (P), and excess air (A) estimation.

- **Development:** Integrate `hydrosheaf.nuclear.dissolved_gas` into the main `joint_lpm` fitting loop.
- **Action:** 
    - Use raw noble gas concentrations (Ne, Ar, Kr, Xe) to solve for T, P, and A.
    - Calculate atmospheric equivalents (SF6_atm, CFC_atm) with propagated uncertainty.
    - If noble gases are missing, use a regional prior for recharge conditions.
- **Benefit:** Makes Hydrosheaf a self-contained physical model. Reduces "external correction bias" where USGS errors are inherited by Hydrosheaf.

---

## 5. Phase 4: Region-Aware Old-Groundwater Calibration (Old Water)

Most remaining error in the national benchmark resides in the "Old" (>1,000y) and "Very Old" (>30,000y) classes, dominated by 14C and 4He.

- **Improvement:** 
    - **4He Rate Grouping:** Instead of a single national prior for 4He accumulation, use a lookup based on StudyUnit or Aquifer Type (consolidating data from the benchmark itself).
    - **Probabilistic 14C Ensemble:** Instead of "selecting" one 14C correction, use a weighted average of candidates (Netpath, Phreeqc, Fontes-Garnier) based on geochemical coherence.
- **Benefit:** Addresses the largest part of the remaining error budget.

---

## 6. Implementation Timeline

| Milestone | Task | Priority |
| :--- | :--- | :--- |
| **M3.1** | Canonical Full Rerun (35 steps) | High |
| **M3.2** | Soft-Weighting Prototype (80-row test) | High |
| **M3.3** | Internal Dissolved Gas Integration | Medium |
| **M3.4** | 4He Regionalization | Medium |
| **M3.5** | Final M3 Manuscript Benchmark | High |

---

## 7. Conclusion

The way forward for Hydrosheaf is to move from being a **benchmark harness** to a **comprehensive probabilistic engine**. By addressing age-grid resolution, soft-weighting tracer reliability, and internalizing the dissolved-gas physics, Hydrosheaf will provide a more robust and accurate representation of groundwater age across the full national spectrum.
