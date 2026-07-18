# M6 Figure Descriptions

Each entry gives the artifact, the objective/RQ it serves, its locked data source, the
caption, an objective technical description, and a deductive scientific description.
Style: Nature portfolio, 300 dpi, PNG/PDF/TIF.

## Figure 1 — Dataset and evidence-tier design
- **Objective:** RQ1 (readiness). **Source:** `m6_map_coordinates.csv`, `m6_tier_ladder.csv`, `m6_variable_availability.csv`.
- **Caption:** Three Ghanaian groundwater datasets and their evidence tiers. (a) Sample
  locations by dataset. (b) Evidence-tier attainment (Tier 0 majors → Tier 4 full metadata).
  (c) Variable-availability heatmap.
- **Technical:** (a) 464 sample coordinates coloured by dataset; (b) 3×5 tile grid of tier
  attainment; (c) fraction-present of 16 variables across three datasets.
- **Scientific:** Establishes that only Northern Ghana reaches Tier 4; Lower Anayari caps at
  Tier 2 (no Sr/SiO₂) and Talensi at Tier 1 (no F/Sr/SiO₂), fixing the evidence ceiling for
  every downstream inference.
- **Figure style compliance:** in-panel tags a–c, ≤6 major ticks, white background, 300 dpi.

## Figure 2 — Hydrosheaf field-transfer workflow
- **Objective:** framing. **Source:** schematic.
- **Caption:** Transfer of the frozen M5 inverse model and MRS calibration onto three Ghana
  datasets through harmonisation, CBE screening, evidence tiers, candidate edges, transferred
  inference, evidence gates, robustness diagnostics and a limitation map.
- **Technical:** node–arrow flow diagram; colour blocks group the M5, data, inference and
  diagnostic stages.
- **Scientific:** Makes explicit that nothing is re-fit on field data — M6 is a transfer and
  stress test, not a re-calibration.

## Figure 3 — Northern Ghana seasonal/aquifer stability
- **Objective:** RQ2. **Source:** `m6_ng_class_support.csv`, `m6_ng_family_by_aquifer.csv`, `m6_ng_field_pairs.csv`.
- **Caption:** (a) Identifiability class by aquifer; (b) dominant Cl-corrected seasonal
  process by aquifer; (c) Mechanism Resolution Score by aquifer.
- **Technical:** stacked proportion bars (a, b) and boxplots (c) over 160 wells, 4 aquifers.
- **Scientific:** All aquifers report uniformly partially identifiable at full information —
  Hydrosheaf does not over-resolve — while dominant process and MRS vary modestly by aquifer,
  with silicate weathering prevailing across basement and sedimentary settings.

## Figure 4 — Robustness loss under diagnostic ablation
- **Objective:** RQ3 / H1. **Source:** `m6_tier_ablation.csv`, `m6_tier_ablation_transitions.csv`.
- **Caption:** (a) Identifiability distribution vs evidence tier; (b) fraction non-identifiable
  and mean MRS vs tier; (c) process-label flips vs Tier 4.
- **Technical:** proportion bars, dual-axis line plot, and a bar of family-change fraction.
- **Scientific:** The pivotal result — removing Sr/SiO₂ (Tier 3→2) sends 52% of wells to
  non-identifiable while mean MRS stays ≈70, demonstrating that numerical fit quality masks
  the loss of mechanistic resolution.

## Figure 5 — External sparse transfer
- **Objective:** RQ5 / H4. **Source:** `m6_external_transfer.csv`, `m6_external_summary.csv`.
- **Caption:** (a) Mean MRS, external vs matched-tier Northern Ghana reference;
  (b) identifiability distribution; (c) dominant process, per external dataset.
- **Technical:** horizontal bars and stacked proportion bars over 129 (Talensi) and 85
  (Lower Anayari) chemistry-kNN edges plus 240 matched-tier reference edges.
- **Scientific:** External transfer sharply increases non-identifiability — 96.5% for Lower
  Anayari, whose silicate-dominant waters cannot be resolved without Sr/SiO₂ — quantifying the
  cost of missing diagnostics and identifying the next-best measurement.

## Figure 6 — Limitation map
- **Objective:** RQ6 / H6. **Source:** `m6_limitation_map.csv`, `m6_conservative_vs_bestfit.csv`.
- **Caption:** (a) Process-class identifiability by dataset; (b) conservative class-level vs
  single best-fit reaction commitment.
- **Technical:** identifiability-scored heatmap (process × dataset) and a jittered
  best-fit-vs-conservative scatter with a 1:1 reference.
- **Scientific:** Silicate and carbonate weathering are non-identifiable in both external
  datasets but resolvable in Northern Ghana, whereas majors-corroborable processes are
  resolvable everywhere; conservative reporting consistently admits more alternatives than a
  single best fit, i.e. avoids over-commitment.

## Supplementary figures
All supplementary figures use the same Nature-style theme (enlarged fonts, collected/
aligned bottom legends) and are composed as multi-panel figures where more than one view
adds value.

- **S1** Hydrochemical context — 2×2: (a) Gibbs cation ratio, (b) Gibbs anion ratio,
  (c) major-ion composition, (d) mineralisation (TDS) by dataset.
- **S2** Charge balance — 2 panels: (a) CBE distribution (±5/±10% guides), (b) quality-class
  composition by dataset.
- **S3** Seasonal signal — 2 panels: (a) dry/wet evapoconcentration factor by aquifer,
  (b) reactive residual after transport correction by aquifer.
- **S4** Alternative edge networks — 2 panels: (a) inferred process-network composition per
  edge set, (b) divergence (TVD) and stability versus the provided graph.
- **S5** Dominant-process composition across evidence tiers (single heatmap).
- **S6** Diagnostic distributions — 2 panels: (a) bootstrap support stability, (b) MRS.
- **S7** Reactive residual versus conservative evapoconcentration (single, coloured by process).
- **S8** MRS versus support stability (single, coloured by process).
- **S9** External-dataset transfer detail — faceted by dataset (Talensi, Lower Anayari).
- **S10** Identifiability composition per process family and dataset — faceted by dataset.
