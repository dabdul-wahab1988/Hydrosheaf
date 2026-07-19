# M6 Results Summary (locked)

Deterministic run (seed 1234), 866 s. All numbers below are reproducible via
`python scripts/run_m6_all.py`. Inferred Mendeley labels are concordance references,
never validation targets.

Under revised Objective 6, `Tier 4` is the maximum M6 chemistry/metadata tier,
not a fully observed age–head–screen field integration. The objective-level
interpretation is in `docs/objective6_data_limited_synthesis.md`.

## E1 — Dataset readiness
- Independent charge-balance error reproduces the Northern Ghana workbook `Data_Class`
  split exactly: **294 quantitative / 19 screening / 7 exploratory** (median |CBE| 1.5%).
- Lower Anayari: 36/41 quantitative, median |CBE| 3.1% (Tier 2).
- Talensi: **0 quantitative, 58/63 exploratory**, median |CBE| ≈ 30% (Tier 1) — a genuine
  low-QC external set (Na–HCO₃ anion excess, no reported cation completeness).

## E2 — Northern Ghana seasonal × aquifer transfer (Tier 4)
- **All 160 wells classified partially identifiable** — reproducing the M5 transferred-
  calibration finding and confirming conservative, class-level (not unique-phase) reporting.
- Mean MRS ≈ 70–72; mean bootstrap support stability ≈ 0.96–0.98 across all four aquifers.
- Dominant Cl-corrected seasonal process: silicate weathering (106/160), redox (23),
  ion exchange (18), anthropogenic (6), carbonate (4), evaporite (3).
- Concordance with inferred prior `Dominant_Process` labels = **0.42** (regolith/alluvial
  0.61, granitoid 0.33, Voltaian 0.31) — moderate, confirming the inferred labels are not
  ground truth and that Hydrosheaf's reactive-seasonal signal diverges from static facies.

## E3 — Tier ablation (H1 confirmed, sharply)
| Tier | non-identifiable | partial | mean MRS | family-flip vs T4 |
|---|---:|---:|---:|---:|
| T4 (full) | 0% | 100% | 70.9 | 0% |
| T3 (+Sr/SiO₂) | 0.6% | 99.4% | 70.7 | 33% |
| T2 (+F) | 51.9% | 48.1% | 70.7 | 33% |
| T1 (+isotopes) | 51.9% | 48.1% | 68.6 | 33% |
| T0 (majors) | 60.0% | 40.0% | 68.4 | 33% |

- **Removing Sr/SiO₂ (T3→T2) collapses 52% of wells to non-identifiable**, whereas removing
  metadata/SI/edges (T4→T3) barely changes anything (0→0.6%). Sr/SiO₂ is the load-bearing
  diagnostic for resolving silicate weathering — confirming H1.
- Mean MRS stays ≈ 68–71 throughout: **numerical fit quality does not reveal the loss of
  identifiability** (the central claim). Process labels flip for a third of wells below Tier 4.

## E4 — Edge / flow-path uncertainty (H3 confirmed)
- Process-*network* composition shifts materially with the assumed edge set (total-variation
  distance vs provided graph: chemistry-kNN 0.23, random 0.17, geographic 0.09).
- But per-edge identifiability class is edge-invariant: 98.6–100% partially identifiable,
  mean MRS 71.9–73.3 across all four edge sets. Connectivity assumptions affect
  network-level attribution, not point-reaction identifiability.

## E5 — External sparse transfer (H4 confirmed)
| Dataset | non-identifiable | partial | mean MRS |
|---|---:|---:|---:|
| Talensi (Tier 1) | 36% | 64% | 69.4 |
| N. Ghana ref (Tier 1) | 27% | 73% | 71.1 |
| Lower Anayari (Tier 2) | **96.5%** | 3.5% | 70.5 |
| N. Ghana ref (Tier 2) | 27% | 73% | 73.0 |
- Lower Anayari's edges are 81/85 silicate-dominant; without Sr/SiO₂ the silicate signal
  cannot be corroborated → **96.5% non-identifiable**. The actionable message: **Sr/SiO₂ is
  the next-best measurement** for Lower Anayari. Talensi degrades through both missing
  tracers and severe charge imbalance.

## E6 — Limitation map & reporting style
- Identifiable-tier by process class: silicate & carbonate are **non-identifiable in both
  external datasets** but partially identifiable in Northern Ghana; majors-corroborable
  processes (nitrate, redox, cation exchange) remain partially identifiable everywhere.
- Conservative class-level reporting admits **7.3** reaction alternatives on average versus
  **5.8** committed by single best-fit — quantifying the honesty premium of evidence-gated
  reporting.

## Validation & reviewer-response (E7–E8; see docs/m6_reviewer_response.md)
- **Circularity resolved (Figure 8a).** With the evidence gate OFF, the field
  non-identifiable fraction is **0% at every tier** — the tier "collapse" is the
  conservative prior, not the classifier. Its direction is validated independently:
  extended-model structural rank rises 8→11 and silicate coherence falls 0.707→0.500
  when Sr/SiO2 are added.
- **Synthetic-field validation with known truth (Figure 7).** Sr/SiO2 improve
  ground-truth dominant-process recovery where theory predicts (cation-exchange
  0.45→0.66; carbonate 1.00; silicate ~0.98), while exact-mineral F1 stays ~0.49–0.74 at
  every tier — validating conservative class/family reporting. True identifiability class
  spans the full range (identifiable 35% / partial 30% / non 23% / equivalence 12%),
  demonstrating the framework discriminates.
- **Transport sensitivity (Figure 8b).** No correction is degenerate (93% evaporite);
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
