# M6 Data-Readiness Audit

Independent audit of the three real Ghanaian datasets used in M6. Charge-balance error
(CBE) is recomputed from harmonised mmol/L ions (not taken from the source files);
Northern Ghana's independent CBE reproduces the workbook's own `Data_Class` split
(294 quantitative / 19 screening / 7 exploratory), validating the harmonisation.

## Summary

| Dataset | n samples | n sites | Seasonal | Native tier | Quantitative (|CBE|≤5%) | Screening (5–10%) | Exploratory (>10%) | Median |CBE| |
|---|---:|---:|:--:|---|---:|---:|---:|---:|
| Northern Ghana (Mendeley) | 320 | 160 | yes (wet/dry, Mar & Aug 2025) | Tier 4 | 294 | 19 | 7 | 1.5% |
| Lower Anayari (manu) | 41 | 41 | no | Tier 2 | 36 | 5 | 0 | 3.1% |
| Talensi (mining area) | 63 | 63 | no | Tier 1 | 0 | 5 | 58 | 29.9% |

## Interpretation

- **Northern Ghana** is the only dataset at the maximum M6 chemistry/metadata
  tier: it carries major ions, F, Sr, SiO₂, δ18O, δ2H, saturation indices,
  aquifer/geology/lithology metadata, QC classes and candidate graph edges. It
  supports the M6 reaction-component diagnostic workflow, not a fully observed
  age–head–screen integration.
- **Lower Anayari** is a clean Tier-2 external set (majors + F + isotopes; no Sr/SiO₂,
  no season, no metadata). 36/41 samples pass quantitative CBE — suitable for external
  reaction-class transfer, but silicate/carbonate corroboration is limited without Sr/SiO₂.
- **Talensi** is a genuinely stressed Tier-1 external set. It has **no F, Sr, SiO₂**,
  no season and no reported QC, and **58/63 samples exceed ±10% CBE** (median |CBE| ≈ 30%).
  The imbalance reflects Na–HCO₃ waters with a large anion excess and no reported cation
  completeness. Talensi is therefore a **screening-only** transfer target and a documented
  failure-mode case: it demonstrates where sparse, low-QC data limit defensible inference.

## Claim strength by dataset

- Northern Ghana: reaction-component field workflow; class-level reaction
  inference with evidence gates; conservative identifiability reporting.
- Lower Anayari: external reaction-class transfer; reduced corroboration; wider uncertainty.
- Talensi: screening only; charge-balance-limited; conclusions reported as low-confidence
  with explicit failure-mode flags.

## Honesty boundary

No dataset provides independent reaction, flow-path or age truth. The Mendeley inferred
labels (`Dominant_Process`, `Aquifer_Evolution_Label`, `Graph_Edges`) are used only as
concordance references. M6 reports transferability and robustness, not field validation.
