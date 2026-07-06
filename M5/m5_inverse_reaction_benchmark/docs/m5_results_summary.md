# M5 Complete Analysis Summary

- Live PHREEQC scenarios: 240.
- Factorial inverse fits: 21600.
- Northern Ghana quantitative wet-dry pairs: 138.
- Maximum PHREEQC-to-stoichiometric generation RMSE: 6.562e-04 mmol/L.
- Primary Hydrosheaf Guarded 3% noise full-panel phase F1 / equivalence-class F1: 0.563 / 0.609.
- Primary Hydrosheaf Guarded false-discovery rate: 0.361; Hydrosheaf-Core evidence-gated comparator class F1 / false-discovery rate: 0.606 / 0.383; legacy thermodynamic elastic-net false-discovery rate: 0.465.
- Hydrosheaf-Core mean reaction-evidence score: 0.622; fraction of high-evidence synthetic reaction rows that were truly active: 0.230.
- Data-tier experiment at 3% noise: Core 0.606 class F1, 0.383 FDR; Plus-lite 0.610 class F1, 0.364 FDR; Enhanced 0.614 class F1, 0.361 FDR.
- Evidence-lifted resolution of ambiguous reaction classes: Core 0.225 mean ELRI, 0.551 conditionally preferred/resolved; Plus-lite 0.243 mean ELRI, 0.628 conditionally preferred/resolved; Enhanced 0.232 mean ELRI, 0.628 conditionally preferred/resolved. ELRI quantifies within-class evidence separation, not new mass-balance uniqueness.
- Conventional PHREEQC inverse-model baseline success fraction (strict 5% plus relaxed 20% fallback): 0.996; strict-only success fraction: 0.458; first-minimal phase F1 / equivalence-class F1: 0.510 / 0.512.
- PHREEQC inverse-model multiplicity: 185.78 feasible models and 6.79 minimal models per scenario on average.
- Fraction of lowest-residual fits with phase F1 below 0.80: 0.550.
- Mixed-archetype held-out MRS classification accuracy: 0.489.
- Ghana median Hydrosheaf-Core evidence score / TDS consistency score: 0.681 / 0.941; pairs with optional SiO2/Sr/isotope support: 138.
- External field ELRI transfer: NorthernGhana.xlsx: 160 edges, 0.072 median ELRI; Talensi: 85 edges, 0.072 median ELRI; LowerAnayari: 49 edges, 0.010 median ELRI. These are field plausibility audits, not reaction-truth validation.

## Claim Guardrail

Hydrosheaf is evaluated as an identifiability-aware sparse linear inverse reaction model with PHREEQC-generated truth and thermodynamic screening, not as a fully coupled nonlinear PHREEQC inverse solver.

The Northern Ghana component is a chemistry-only seasonal transfer demonstration without independent flow-path, groundwater-age, or reaction-truth validation.
