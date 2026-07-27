# M6 Data-Readiness Audit

Independent audit of the three real Ghanaian datasets used in M6, sourced strictly
from `data/FieldData/` (`NorthenGhana/NorthernGhana.xlsx`, `Talensi_MiningArea/talensi.csv`,
`LowerAnayari/manu.csv`; see `DECISIONS.md`). Charge-balance error (CBE) is recomputed
from harmonised mmol/L ions (not taken from any source file); no external quality-class
field exists in the canonical raw data to check this classification against.

## Summary

| Dataset | n samples | n sites | Seasonal | Native tier | Quantitative (|CBE|≤5%) | Screening (5–10%) | Exploratory (>10%) | Median |CBE| |
|---|---:|---:|:--:|---|---:|---:|---:|---:|
| Northern Ghana (canonical) | 320 | 160 | yes (wet/dry, Mar & Aug 2025) | Tier 4 | 294 | 19 | 7 | 1.5% |
| Lower Anayari (manu) | 41 | 41 | no | Tier 2 | 36 | 5 | 0 | 3.1% |
| Talensi (mining area) | 63 | 63 | no | Tier 1 | 0 | 5 | 58 | 29.9% |

## Interpretation

- **Northern Ghana** is the only dataset at the maximum M6 chemistry evidence
  tier: it carries major ions, F, Sr, SiO₂, δ18O, δ2H and Hydrosheaf-computed
  PHREEQC saturation indices. No independent aquifer-type, geology, lithology,
  land-use, QC-class or graph-edge data exists for these boreholes; stratified
  reporting elsewhere in M6 uses administrative region/district instead. It
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

No dataset provides independent reaction, flow-path or age truth. An earlier revision of
this workflow additionally read a different, antecedent study's own inferred labels for
the Northern Ghana boreholes (`Dominant_Process`, `Aquifer_Evolution_Label`, `Graph_Edges`,
from `Aquifers_Dataset_Mendeley.xlsx`); those fields, and the workbook that supplied them,
have been removed from this study entirely rather than retained as a concordance reference
(`DECISIONS.md`). M6 reports transferability and robustness, not field validation.
