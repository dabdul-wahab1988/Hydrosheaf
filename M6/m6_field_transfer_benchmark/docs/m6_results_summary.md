# M6 Results Summary (locked)

Deterministic run (seed 1234). All numbers below are reproducible via
`python scripts/run_m6_q1.py`. Northern Ghana chemistry comes only from the
canonical raw workbook (`data/FieldData/NorthenGhana/NorthernGhana.xlsx`); no
aquifer/geology/lithology metadata, prior process labels, or externally
provided graph edges are used (see `DECISIONS.md`).

Under revised Objective 6, `Tier 4` is the maximum M6 chemistry tier (majors,
isotopes, fluoride, Sr/SiO2, and independently PHREEQC-derived saturation
indices), not a fully observed age-head-screen field integration and not a
metadata/graph-edge tier. The objective-level interpretation is in
`docs/objective6_data_limited_synthesis.md`.

## E1 — Dataset readiness
- Independent charge-balance error classifies Northern Ghana as **294
  quantitative / 19 screening / 7 exploratory** of 320 records (median |CBE|
  1.5%). This is Hydrosheaf's own calculation from the raw ion panel; no
  external quality-class label exists in the canonical data to check it
  against.
- Lower Anayari: 36/41 quantitative, median |CBE| 3.1% (Tier 2).
- Talensi: **0 quantitative, 58/63 exploratory**, median |CBE| ≈ 30% (Tier 1) — a genuine
  low-QC external set (Na–HCO₃ anion excess, no reported cation completeness).

## E2 — Northern Ghana seasonal × region transfer (Tier 4)
- **All 160 wells classified partially identifiable** — reproducing the M5 transferred-
  calibration finding and confirming conservative, class-level (not unique-phase) reporting.
- Mean MRS ≈ 70.6–71.2; mean bootstrap support stability ≈ 0.965–0.976 across all four
  administrative regions (Northern, North East, Upper East, Upper West; no
  independent aquifer-type classification exists for these boreholes).
- Dominant Cl-corrected seasonal process: silicate weathering (114/160), redox (18),
  ion exchange (16), anthropogenic (5), carbonate (4), evaporite (3).
- No prior process-label concordance is reported: the only such reference available
  before this fix was a different, antecedent study's own inferred labels, not
  independent ground truth, and it is not part of the canonical field data (see
  `DECISIONS.md`).

## E3 — Tier ablation (H1 confirmed, sharply)
| Tier | non-identifiable | partial | mean MRS | family-flip vs T4 |
|---|---:|---:|---:|---:|
| T4 (full) | 0% | 100% | 71.0 | 0% |
| T3 (+Sr/SiO₂) | 0.6% | 99.4% | 70.7 | 35.6% |
| T2 (+F) | 51.9% | 48.1% | 70.7 | 35.6% |
| T1 (+isotopes) | 51.9% | 48.1% | 68.6 | 35.6% |
| T0 (majors) | 60.0% | 40.0% | 68.4 | 35.6% |

- **Removing Sr/SiO₂ (T3→T2) collapses 52% of wells to non-identifiable**, whereas removing
  the saturation-index lifter (T4→T3) barely changes anything (0→0.6%). Sr/SiO₂ is the
  load-bearing diagnostic for resolving silicate weathering — confirming H1.
- Mean MRS stays ≈ 68–71 throughout: **numerical fit quality does not reveal the loss of
  identifiability** (the central claim). Dominant-family labels flip for roughly a third
  of wells below Tier 4.

## E4 — Edge / flow-path uncertainty (H3 confirmed)
- Three Hydrosheaf-generated edge sets are compared against each other (chemistry-similarity
  kNN, geographic-nearest, random-perturbed); an earlier revision also compared against a
  fourth "provided graph" edge set imported from a retired antecedent-study workbook, which
  has been removed rather than substituted (see `DECISIONS.md`). Process-*network*
  composition still shifts materially with the assumed edge set (total-variation distance vs
  the chemistry-kNN reference: geographic-nearest 0.12, random-perturbed 0.05).
- But per-edge identifiability class is edge-invariant: 99.3–100% partially identifiable,
  mean MRS 72.6–73.3 across all three edge sets. Connectivity assumptions affect
  network-level attribution, not point-reaction identifiability.

## E5 — External sparse transfer (H4 confirmed, direction reversed for Talensi)
Talensi and Lower Anayari measure complete pH, temperature and major ions, so PHREEQC
saturation indices are computed for them exactly as for Northern Ghana (not withheld on the
assumption that they require Sr/SiO2 — see `DECISIONS.md`) and used as an evidence-lift-gate
lifter for both external datasets throughout.

| Dataset | non-identifiable | partial | mean MRS |
|---|---:|---:|---:|
| Talensi (Tier 1) | 37.2% | 62.8% | 68.8 |
| N. Ghana ref (Tier 1) | 54.2% | 45.8% | 70.9 |
| Lower Anayari (Tier 2) | **95.3%** | 4.7% | 69.9 |
| N. Ghana ref (Tier 2) | 53.3% | 46.7% | 73.1 |
- Talensi is now materially *lower* than its matched Northern Ghana reference (37.2% vs
  54.2%), reversing the pre-SI-fix "essentially the same rate" reading: Talensi's edges are
  rarely carbonate-dominant, so the newly available saturation-index lifter was rarely the
  limiting evidence for it, while it was the limiting evidence for a larger share of the NG
  reference edges.
- Lower Anayari's edges are silicate-dominant; silicate corroboration requires Sr/SiO2
  specifically (saturation indices do not corroborate it) → **95.3% non-identifiable**. The
  actionable message is unchanged: **Sr/SiO₂ is the next-best measurement** for Lower
  Anayari.

## E6 — Limitation map & reporting style
- Identifiable-tier by process class: silicate remains **non-identifiable in both external
  datasets** but partially identifiable in Northern Ghana (Sr/SiO2 still absent externally).
  Carbonate is now partially identifiable in all three datasets, since the small number of
  carbonate-dominant external edges are corroborated by the now-available saturation
  indices. Majors-corroborable processes (nitrate, redox, cation exchange) remain partially
  identifiable everywhere.
- Conservative class-level reporting admits **7.1** reaction alternatives on average versus
  **5.4** committed by single best-fit — quantifying the honesty premium of evidence-gated
  reporting.

## Validation and robustness diagnostics

See `docs/m6_robustness_diagnostics.md`.
- **Circularity resolved (Extended Data Figure ED3a).** With the evidence gate OFF, the field
  non-identifiable fraction is **0% at every tier** — the tier "collapse" is the
  conservative prior, not the classifier. Its direction is validated independently:
  extended-model structural rank rises 8→11 and silicate coherence falls 0.707→0.500
  when Sr/SiO2 are added.
- **Synthetic-field validation with known truth (Extended Data Figure ED2).** Sr/SiO2 improve
  ground-truth dominant-process recovery where theory predicts (cation-exchange
  0.45→0.66; carbonate 1.00; silicate ~0.98), while exact-mineral F1 stays ~0.49–0.74 at
  every tier — validating conservative class/family reporting. True identifiability class
  spans the full range (identifiable 35% / partial 30% / non 23% / equivalence 12%),
  demonstrating the framework discriminates.
- **Transport sensitivity (Extended Data Figure ED3b).** No correction is degenerate (93% evaporite);
  Cl-conservative and recharge-endmember corrections agree on 64% of dominant families and
  both give a silicate-dominated network — correction is necessary, conclusion is robust.
- **Talensi diagnosed** as an unmeasured-cation-pool imbalance (2.1 meq deficit, HCO3 73%
  of anions); screening-only.

## Headline
Hydrosheaf transfers to real Ghanaian data as a **conservative, evidence-gated** framework:
it reports all defensible field cases as partially identifiable rather than over-committing,
it exposes Sr/SiO₂ as the pivotal diagnostic for silicate-weathering resolution, and it shows
that good numerical fit and stable support persist even where identifiability has collapsed —
so fit quality alone must not be trusted in data-scarce aquifers.
