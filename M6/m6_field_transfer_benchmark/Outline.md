# M6 — Q1 Nature-Portfolio Manuscript Outline

Field transfer and robustness of Hydrosheaf under data scarcity across Ghanaian aquifers.

Development status: drafted 9 July 2026. Grounded in the audited real datasets
(`data/NorthenGhana/Aquifers_Dataset_Mendeley.xlsx`,
`data/Talensi_MiningArea/talensi.csv`, `data/LowerAnayari/manu.csv`). Replaces the
hard-coded `M6/m6_robustness_benchmark` scaffold.

---

## Recommended Journal Positioning

Primary target: **Communications Earth & Environment** (Nature Portfolio), Article
(~5,000 words main text). Positions M6 as the **robustness-and-transfer pillar** of
Hydrosheaf: it does not add a new inference method; it stress-tests whether the M5
identifiability-aware workflow survives real, uneven, data-scarce African aquifer
datasets.

Stretch target after independent field truth is secured: **Nature Water**.

Fallback (fully rigorous, applied): *Water Research* or *Journal of Hydrology*.

## Working Title

**Field transfer of Hydrosheaf under data scarcity: robustness, uncertainty and limits across Ghanaian aquifers**

Alternatives:

1. **How robust is hydrogeochemical process inference in data-limited aquifers? A Hydrosheaf transfer study across Ghanaian groundwater datasets**
2. **Which groundwater process interpretations survive sparse data? An evidence-gated transfer benchmark**
3. **Transferability, not validation: conservative reaction-class inference across three Ghanaian aquifer datasets**

## Role Within the Hydrosheaf PhD

M6 answers the fourth question in the Hydrosheaf inference chain (see M5 Outline):

1. M4 — Where can groundwater move? (candidate directed edges)
2. M3 — How long has groundwater moved? (residence-time distributions)
3. M5 — What happened chemically along a supported edge? (reaction classes, extents, identifiability)
4. **M6 — How stable is the resulting interpretation when the workflow is moved onto real, data-limited field datasets?**

M5 established what Hydrosheaf can recover under a controlled live-PHREEQC inverse
benchmark. M6 asks whether the workflow remains **useful, stable and honest** when
metadata, tracers, seasonality, connectivity and independent process truth are uneven.

## Central Claim and Claim Boundary

**Central claim.** Hydrosheaf provides a structured, evidence-gated way to state
*which* process interpretations survive sparse field data, *which* collapse, and
*where* additional measurements are scientifically worth collecting — without
overclaiming mechanistic truth.

**What M6 claims:** transferability, stability under seasonal variation and missing
diagnostics, uncertainty control, robustness to assumed connectivity, and transparent
limitation mapping.

**What M6 does NOT claim:** independent field validation of true reactions, true flow
paths, or true residence times. The Mendeley `Dominant_Process`, `Aquifer_Evolution_Label`
and `Graph_Edges` are *inferred/generated labels*, not independent ground truth, and are
used only as internal consistency references, never as validation targets. This boundary
is stated in the Abstract, Methods §2.6, and Discussion §4.4, and is the load-bearing
honesty of the paper.

## Research Questions

1. **RQ1 (Readiness).** What can and cannot be inferred from each Ghanaian dataset given its available evidence tier (major ions → +isotopes → +F → +Sr/SiO₂ → +metadata/graph)?
2. **RQ2 (Field transfer).** When the full M5 workflow is applied to 320 wet/dry Northern Ghana records, how stable are reaction-class support, evidence gates and identifiability class across season and aquifer type?
3. **RQ3 (Diagnostic robustness).** How much do process labels, reaction-class confidence and identifiability class change when Sr/SiO₂, isotopes, aquifer metadata or graph edges are removed?
4. **RQ4 (Connectivity robustness).** Do Hydrosheaf conclusions depend too strongly on assumed edges? How do graph-based, chemistry-similarity, geographic-nearest and random/perturbed edge sets change the inferred process network?
5. **RQ5 (External transfer).** When the workflow is transferred to genuinely sparse external datasets (Talensi n=63, Lower Anayari n=41), what transfers, how much does uncertainty increase, and what are the failure modes?
6. **RQ6 (Limits).** Across all datasets, which process classes are identifiable, weakly identifiable, or non-identifiable under field data scarcity?

## Testable Hypotheses

- **H1.** Removing Sr/SiO₂ or isotopes raises the fraction of samples classified *non-identifiable* and shifts silicate-vs-carbonate weathering attribution more than removing graph edges does.
- **H2.** Reaction-class support (equivalence-class level) is more stable across season and diagnostic removal than exact-phase support.
- **H3.** Inferred process labels change materially under edge-set perturbation only for edge-derived diagnostics (mixing/flowpath), not for point-sample reaction-class inference — isolating where connectivity assumptions matter.
- **H4.** External sparse transfer (Talensi, Manu) yields systematically lower identifiability and wider uncertainty than Northern Ghana at the same evidence tier, quantifying the cost of missing F/Sr/SiO₂ and metadata.
- **H5.** Charge-balance-error screening (Data_Class tiers) changes conclusions: exploratory (CBE>10%) samples inflate apparent process diversity relative to the 294 quantitative records.
- **H6.** The evidence-gated identifiability hierarchy produces fewer, more conservative, more reproducible process claims than a conventional single best-fit hydrochemical interpretation on the same data.

---

## Datasets (audited, real)

| Dataset | n samples | Evidence available | Best-supported tier | M6 role |
|---|---:|---|---|---|
| Northern Ghana (Mendeley workbook; == `NorthernGhana.xlsx` chemistry, use Mendeley as it adds metadata) | 320 (160 wells × wet/dry, Mar & Aug 2025) | majors, NO₃, F, Sr, SiO₂, δ18O, δ2H, SI (calcite/dolomite/gypsum/halite), CBE, SAR, facies, aquifer type (4), geology (4), lithology (12), land use, QC, 397 graph edges, feature matrix, evolution labels | **Tier 4** (full integrated) | Primary field-transfer + robustness testbed |
| Talensi (mining area) | 63 | majors, NO₃, Fe, Eh, δ18O, δ2H, coords; **no F, Sr, SiO₂, no season, no aquifer meta, no explicit QC** | **Tier 1** (majors + isotopes) | External sparse transfer #1 |
| Lower Anayari (manu) | 41 | majors, NO₃, F, Fe, δ18O, δ2H, coords; **no Sr, SiO₂, no season, no aquifer meta, no explicit QC** | **Tier 2** (majors + isotopes + F) | External sparse transfer #2 |

Northern Ghana `Data_Class`: 294 *Quantitative inverse modelling* (|CBE| ≤ 5%), 19
*Screening* (5–10%), 7 *Exploratory* (>10%). Primary analyses use the 294; the remainder
enters sensitivity reporting only.

## Evidence-Tier Ladder (grounded in real availability)

- **Tier 0** — major ions only (Ca, Mg, Na, K, HCO₃, Cl, SO₄, NO₃). Available: all three.
- **Tier 1** — Tier 0 + stable isotopes (δ18O, δ2H). Available: all three. *(Talensi caps here.)*
- **Tier 2** — Tier 1 + F. Available: Northern Ghana, Manu. *(Manu caps here.)*
- **Tier 3** — Tier 2 + Sr + SiO₂. Available: Northern Ghana only.
- **Tier 4** — Tier 3 + aquifer/geology metadata + saturation indices + QC class + graph edges. Available: Northern Ghana only.

The ladder is the experimental spine: Experiment 3 walks Northern Ghana *down* the ladder
(ablation), and Experiment 5 tests external datasets that are natively *stuck* at Tier 1–2.

---

## Article Architecture and Word Budget

Target **5,000 words**, Introduction → Conclusion (excludes Abstract, references, legends,
availability statements, SI).

| Section | Target words |
|---|---:|
| 1. Introduction | 850 |
| 2. Study design and datasets | 800 |
| 3. Hydrosheaf transfer workflow (Methods) | 1,100 |
| 4. Results | 1,750 |
| 5. Discussion | 700 |
| 6. Conclusion | 200 |
| **Total** | **5,000** |

Abstract: 180–200 words, unstructured, no references. Written only after analyses are locked.
Sequence: (1) data scarcity in African aquifer hydrogeochemistry; (2) the transfer/robustness
question; (3) the three-dataset, five-experiment design; (4) principal quantitative results
(tier-loss, seasonal/aquifer stability, edge sensitivity, external-transfer uncertainty);
(5) significance — a conservative, evidence-gated framework that says what survives sparse data.

---

# 1. Introduction — 850 words

## 1.1 Data scarcity in African aquifer hydrogeochemistry — 160 words
Groundwater is the dominant potable source across semi-arid West Africa; process
understanding (weathering, ion exchange, evapoconcentration, geogenic fluoride,
anthropogenic nitrate) underpins vulnerability and treatment decisions. Yet field
datasets are chronically incomplete: uneven tracers, single-season sampling, absent
lithology, unstated QC.

## 1.2 Limits of conventional inverse hydrogeochemistry under scarcity — 170 words
PHREEQC inverse modelling, saturation indices, Gibbs/Piper interpretation and endmember
mixing all assume adequate analytical coverage; on sparse data they either fail or return
a single confident-looking narrative. Equifinality is worst exactly where data are thinnest.

## 1.3 Why field transferability matters — 150 words
A method validated on a synthetic benchmark (M5) is only useful if it survives real,
uneven data. The scientific question is not "is Hydrosheaf right in Ghana?" but "does the
identifiability-aware workflow degrade gracefully, honestly and reproducibly as evidence
is removed?"

## 1.4 The gap M6 addresses — 140 words
No existing workflow reports, on real African aquifer data, a transparent map of which
process interpretations are stable, which collapse under missing diagnostics, and where a
single extra measurement would most reduce ambiguity.

## 1.5 Aim and objectives — 150 words
Aim: *to evaluate how Hydrosheaf transfers to Ghanaian data-limited aquifer datasets by
quantifying the stability of process inference, evidence gates, reaction-class interpretation
and uncertainty under seasonal variation, missing diagnostics, uncertain flow-path structure
and external-dataset transfer.* Five objectives = the five experiments.

## 1.6 Significance — 80 words
Scientific: replaces a single deterministic story with a measurable survival hierarchy.
Practical: identifies the highest-value next measurement per aquifer. Framework: supplies a
conservative honesty layer to Hydrosheaf without overclaiming validation.

---

# 2. Study Design and Datasets — 800 words

## 2.1 Ghanaian aquifer settings (250 words)
Voltaian sedimentary, Birimian and granitoid fractured basement, and regolith/alluvial
shallow aquifers; aridification context (SPEI, rainfall anomaly from the Climate sheet).
Talensi = artisanal-mining-affected basement; Lower Anayari = shallow regolith with geogenic F.

## 2.2 Variable availability and evidence tiers (300 words)
Present the tier ladder as a design, not an afterthought. Table 2 (variable availability
matrix). Explain the March/August 2025 paired design for Northern Ghana and the
single-campaign nature of the external sets.

## 2.3 Why these datasets represent realistic data-limited conditions (250 words)
Missingness, QC unevenness, inferred vs independent labels. State the honesty boundary here.

---

# 3. Hydrosheaf Transfer Workflow (Methods) — 1,100 words

## 3.1 Data harmonisation and unit conversion (150 words)
mg/L → meq/L, ion-name mapping across three schemas, coordinate handling, δ notation.
Reuse M5/M2 field-data contract. All rules in Supplementary Methods.

## 3.2 Charge-balance screening and QC tiers (130 words)
Recompute CBE independently; reproduce Data_Class tiers; define the 294/19/7 split and its
use.

## 3.3 Evidence tiers and the inference unit (150 words)
Define Tier 0–4. The inference unit is a directed edge where connectivity exists; a point
sample (or wet→dry within-well pair) where it does not. State chemistry-only vs integrated
explicitly per dataset.

## 3.4 Candidate-edge construction (180 words)
Four edge-generation strategies for Experiment 4: (a) provided `Graph_Edges`, (b)
chemistry-similarity kNN, (c) geographic-nearest, (d) random/perturbed null. Within-well
wet→dry pairs as a fifth, metadata-free connectivity.

## 3.5 Inverse reaction-class inference and evidence gates (200 words)
Reuse the M5 sparse linear inverse reaction model + PHREEQC screening + identifiability
diagnostics (rank, nullity, coherence, null-space) and the calibrated Mechanism Resolution
Score / four-class identifiability labels. Evidence gates = indicator-ion and thermodynamic
feasibility gates from M5.

## 3.6 Robustness, uncertainty and transfer diagnostics (180 words)
Tier-ablation protocol; edge-set perturbation; bootstrap support stability; seasonal
concordance; per-class identifiability shift; next-best-measurement. Uncertainty metrics
defined here, full spec in SI.

## 3.7 Reproducibility (110 words)
One-command `run_m6_all.py`; environment lock; R figure scripts; deposited derived CSVs;
model-card limitations. Mirror M5 package layout.

---

# 4. Results — 1,750 words

## 4.1 Dataset readiness and missingness (240 words) → **Figure 1**, **Table 1**, **Table 2**
Per-dataset tier attainment, missingness heatmap, CBE distributions, claim-strength matrix.
Verify chemistry equivalence of `NorthernGhana.xlsx` and Mendeley workbook.

## 4.2 Northern Ghana seasonal and aquifer-stratified transfer (330 words) → **Figure 3**
Reaction-class support and identifiability class by aquifer type × season. Wet→dry
concordance. Which processes (silicate weathering, cation exchange, nitrate loading,
evapoconcentration) are stably supported per aquifer. Report as class-level, conservative.

## 4.3 Robustness to removing diagnostics and metadata (330 words) → **Figure 4**, **Table 3**
Walk Northern Ghana down the tier ladder. Quantify identifiability-class migration and
process-label flips as Sr/SiO₂, isotopes, aquifer metadata, graph edges are removed
(tests H1, H2). Full-vs-reduced information table.

## 4.4 Edge and flow-path uncertainty (300 words) → **Figure**(S) + panel in Fig 3/6
Compare inferred process networks under the four edge sets (tests H3). Show which
conclusions are connectivity-invariant (point-sample reaction class) vs
connectivity-sensitive (mixing/flowpath diagnostics).

## 4.5 External transfer to Talensi and Lower Anayari (320 words) → **Figure 5**, **Table 4**
Apply workflow at native Tier 1 / Tier 2. Report transfer success, uncertainty increase vs
Northern Ghana at matched tier, and failure modes (e.g. F-dependent fluoride mobilisation
unresolvable in Talensi).

## 4.6 Limitation map and process-class identifiability (230 words) → **Figure 6**
Synthesise all experiments into an identifiable / weakly identifiable / non-identifiable
map per process class per dataset. Compare conservative evidence-gated reporting vs
conventional single best-fit (tests H6).

---

# 5. Discussion — 700 words

## 5.1 What Hydrosheaf adds beyond conventional chemistry interpretation (200 words)
## 5.2 Where it reduces overinterpretation / equifinality it cannot resolve (200 words)
## 5.3 Monitoring implications: next-best measurement per aquifer (150 words)
## 5.4 Interpretation and honesty boundary of the field application (150 words)
Field transfer of a chemistry audit — not validation of connectivity, age or reaction truth.

# 6. Conclusion — 200 words
Hydrosheaf is useful as a conservative, evidence-gated field-transfer framework, but not a
substitute for independent flow/reaction validation. Answer each RQ directly.

---

# Display-Item Plan

## Main Figures (Nature style, built in R)
1. **Ghana dataset map + evidence-tier design.** Panel a: map of the three study areas with sample counts and aquifer types. Panel b: evidence-tier ladder (Tier 0–4) with per-dataset attainment. Panel c: variable-availability / missingness heatmap.
2. **Hydrosheaf field-transfer workflow.** Schematic from the M5 synthetic benchmark to the M6 Ghana field application, showing harmonisation → CBE screening → tier gating → edge construction → inverse reaction-class inference → identifiability/evidence gates → robustness diagnostics.
3. **Northern Ghana wet/dry process-inference stability by aquifer type.** Reaction-class support and identifiability class across 4 aquifer types × 2 seasons; wet→dry concordance.
4. **Robustness loss under missing diagnostics and metadata.** Identifiability-class migration and process-label flips as the tier ladder is descended (Tier 4 → Tier 0).
5. **External transfer to Talensi and Lower Anayari.** Transfer success, identifiability distribution, and uncertainty increase vs Northern Ghana at matched tier; failure-mode annotations.
6. **Limitation map.** Identifiable / weakly identifiable / non-identifiable process classes per dataset; conservative evidence-gated vs conventional best-fit comparison.

## Main Tables
1. **Dataset readiness and claim strength** (per dataset: n, tier, QC, what can/cannot be claimed).
2. **Variable availability by dataset and tier.**
3. **Hydrosheaf outputs under full vs reduced information** (Northern Ghana tier ablation).
4. **External transfer performance and uncertainty summary** (Talensi, Manu).

## Supplementary Figures
- S1. Piper / Gibbs / key ion-ratio plots (all datasets).
- S2. Charge-balance-error distributions and Data_Class screening.
- S3. Aquifer-specific wet/dry chemistry changes.
- S4. Alternative edge networks (graph vs chemistry vs geographic vs random).
- S5. Tier-ablation sensitivity heatmaps (per process class).
- S6. Bootstrap support-frequency stability.
- S7. Seasonal concordance diagnostics.
- S8. Next-best-measurement rankings per aquifer.
- S9. External-dataset identifiability detail (Talensi, Manu).
- S10. Saturation-index distributions vs inferred processes.

## Supplementary Tables
- S1. Harmonisation and unit-conversion rules; ion-name crosswalk across the three schemas.
- S2. Full missingness and QC tables.
- S3. Aquifer / season summary statistics.
- S4. Complete tier-ablation results (all classes, all tiers).
- S5. Edge-sensitivity results (all four edge sets).
- S6. External-dataset transfer outputs (per sample).
- S7. Uncertainty-metric definitions and per-dataset values.
- S8. Software and computational environment.

## Supplementary Methods
Harmonisation rules, unit conversion, independent CBE calculation, evidence-tier
definitions, graph-edge generation algorithms, robustness/stress-test design, uncertainty
metrics, and the honesty-boundary protocol for inferred vs independent labels.

---

# Claim Guardrails

Use:
> Hydrosheaf is an identifiability-aware sparse linear inverse reaction model with
> PHREEQC-derived thermodynamic screening; M6 evaluates its field transferability and
> robustness under data scarcity.

Do not use:
> M6 validates the true reactions, flow paths, or groundwater ages of Ghanaian aquifers.

Do not equate:
- an inferred `Dominant_Process`/`Evolution_Label` with ground truth;
- a provided `Graph_Edge` with a validated flow path;
- transfer success with mechanistic validation;
- low charge-balance error with correct process attribution.

# Data & Analysis Status — LOCKED (9 July 2026)
The executable M6 package is built and run at `M6/m6_field_transfer_benchmark/`, mirroring
`M5/m5_inverse_reaction_benchmark/`. Deterministic (seed 1234), 866 s.

## Implemented evidence
- Three real Ghanaian datasets harmonised (320 Northern Ghana wet/dry + metadata/SI/edges;
  63 Talensi; 41 Lower Anayari). Independent CBE reproduces the workbook Data_Class split
  (294/19/7), validating harmonisation.
- Reuses the frozen M5 inverse solver and the frozen M5 MRS calibration as the *transferred*
  identifiability classifier (never re-fit). Cl-conservative transport correction isolates
  reactive residuals.
- Six experiments executed → 19 results CSVs; 4 main + 8 supplementary tables (CSV + MD);
  6 main + 10 supplementary Nature-style figures (PNG/PDF/TIF, 300 dpi).

## Principal locked results
- **E2:** all 160 Northern Ghana wells partially identifiable at Tier 4 (conservative,
  reproduces M5 transfer); mean MRS ≈70–72; prior-label concordance 0.42 (labels ≠ truth).
- **E3 (H1):** removing Sr/SiO₂ (Tier 3→2) sends 52% of wells to non-identifiable while mean
  MRS stays ≈70 — fit quality masks identifiability loss. Sr/SiO₂ is the pivotal diagnostic.
- **E4 (H3):** process-network composition shifts with edge choice (TVD ≤0.23) but per-edge
  identifiability is edge-invariant (~99% partial).
- **E5 (H4):** external transfer degrades — Talensi 36% non-identifiable, Lower Anayari 96.5%
  (silicate un-corroborable without Sr/SiO₂); Sr/SiO₂ = next-best measurement.
- **E6 (H6):** silicate/carbonate non-identifiable in both external sets, resolvable in
  Northern Ghana; conservative reporting admits 7.3 vs 5.8 best-fit reactions.

See `docs/m6_results_summary.md`, `docs/02_figures.md`, `docs/03_tables.md`,
`docs/m6_data_readiness_audit.md`, `docs/m6_artifact_manifest.md`.

Old `M6/m6_robustness_benchmark/` retained as reference scaffold only.
