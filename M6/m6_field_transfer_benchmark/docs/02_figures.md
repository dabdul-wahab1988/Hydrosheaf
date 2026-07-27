# M6 Figure Descriptions

Each entry gives the artifact, the objective/RQ it serves, its locked data source, the
caption, an objective technical description, and a deductive scientific description.
Style: Nature portfolio, 300 dpi, PNG/PDF/TIF.

## Figure 1 — Dataset and evidence-tier design
- **Objective:** RQ1 (readiness). **Source:** `m6_map_coordinates.csv`,
  `m6_tier_ladder.csv`, `m6_variable_availability.csv`, and the unsimplified
  2021 geoBoundaries GHA ADM1 Shapefile (16 regions).
- **Caption:** Three Ghanaian groundwater datasets and their evidence tiers. (a) Sample
  locations by dataset. (b) Evidence-tier attainment (Tier 0 majors → Tier 4 + saturation
  indices). (c) Variable-availability heatmap.
- **Technical:** (a) supplied masked coordinates plotted on a sourced Ghana
  outline and regional boundaries, with no site-scale spatial interpretation;
  Talensi positive degrees-west longitudes are converted to signed WGS84 and
  out-of-bound Northern Ghana masks are constrained 5 km inside the national
  boundary for display only; (b) 3×5 tile grid of tier attainment; (c)
  fraction-present of 16 variables across three datasets.
- **Scientific:** Establishes that only Northern Ghana reaches Tier 4; Lower Anayari caps at
  Tier 2 (no Sr/SiO₂) and Talensi at Tier 1 (no F/Sr/SiO₂), fixing the evidence ceiling for
  every downstream inference. Tier 4 is the maximum M6 chemistry tier only; it does not
  imply that environmental age tracers, screen intervals, repeated heads, independent flow
  truth, or any aquifer/geology/lithology metadata are available.
- **Figure style compliance:** in-panel tags a–c, ≤6 major ticks, white background, 300 dpi.

## Figure 2 — Hydrosheaf field-transfer workflow
- **Objective:** framing. **Source:** schematic.
- **Caption:** Transfer of the frozen M5 inverse model and MRS calibration onto three Ghana
  datasets through harmonisation, CBE screening, evidence tiers, candidate edges, transferred
  inference, evidence gates, robustness diagnostics and a limitation map.
- **Technical:** node–arrow flow diagram; colour blocks group the M5, data, inference and
  diagnostic stages.
- **Scientific:** Makes explicit which frozen components are transferred and which field
  diagnostics are computed. M6 is a chemistry/metadata transfer and stress test, not a
  validation of residence time, exact directed topology or a field digital twin.

## Figure 3 — Northern Ghana seasonal/region stability
- **Objective:** RQ2. **Source:** `m6_ng_class_support.csv`, `m6_ng_family_by_aquifer.csv`, `m6_ng_field_pairs.csv`.
- **Caption:** (a) Identifiability class by region; (b) dominant Cl-corrected seasonal
  process by region; (c) Mechanism Resolution Score by region.
- **Technical:** stacked proportion bars (a, b) and boxplots (c) over 160 wells, 4
  administrative regions (no independent aquifer-type classification exists for these
  boreholes).
- **Scientific:** All regions report uniformly partially identifiable at full information —
  Hydrosheaf does not over-resolve — while dominant process and MRS vary modestly by region,
  with silicate weathering prevailing across all four regions.

## Figure 4 — Robustness loss under diagnostic ablation
- **Objective:** RQ3 / H1. **Source:** `m6_tier_ablation.csv`, `m6_tier_ablation_transitions.csv`.
- **Caption:** (a) Identifiability distribution vs evidence tier; (b) fraction non-identifiable
  and mean MRS vs tier; (c) process-label flips vs Tier 4.
- **Technical:** proportion bars, dual-axis line plot, and a bar of family-change fraction.
- **Scientific:** The pivotal result — removing Sr/SiO₂ (Tier 3→2) sends 52% of wells to
  non-identifiable while mean MRS stays ≈70, demonstrating that numerical fit quality masks
  the loss of mechanistic resolution.

## Figure 5 — Within-campaign seasonal hold-forward performance
- **Objective:** revised Objective 6. **Source:** the locked supporting-validation
  `field_prequential_summary.csv`, `field_prequential_predictions.csv` and
  `field_prequential_audit.json`.
- **Caption:** One-step seasonal hold-forward performance for persistence,
  expanding-mean and graph-ridge predictors, including paired well-block
  contrasts, ion-level errors and empirical interval coverage.
- **Technical:** four panels report overall MAE, paired bootstrap differences,
  ion-level descriptive errors and 90% predictive-interval coverage.
- **Scientific:** Graph ridge improves on persistence, but its incremental
  advantage over the expanding-mean baseline is small and its paired interval
  crosses zero. This is truth-free within-campaign chemistry prediction, not
  independent-basin validation of residence time, topology or reactions.

## Figure 6 — Limitation map
- **Objective:** RQ6 / H6. **Source:** `m6_limitation_map.csv`, `m6_conservative_vs_bestfit.csv`.
- **Caption:** (a) Process-class identifiability by dataset; (b) conservative class-level vs
  single best-fit reaction commitment.
- **Technical:** identifiability-scored heatmap (process × dataset) and a jittered
  best-fit-vs-conservative scatter with a 1:1 reference.
- **Scientific:** Silicate- and carbonate-related class assignments are non-identifiable in
  both external datasets and more constrained in Northern Ghana at the M6 class/equivalence
  level. This is not evidence of unique reaction-family recovery. Conservative reporting
  consistently admits more alternatives than a single best fit and therefore avoids
  over-commitment.

## Extended Data figures

- **Extended Data Figure 1 — External sparse transfer.** Mean MRS,
  identifiability and dominant-process distributions for Talensi and Lower
  Anayari against matched-tier Northern Ghana references. It quantifies the
  cost of missing diagnostics without treating prior labels as field truth.
- **Extended Data Figure 2 — Synthetic validation.** Positive-control recovery checks for the transferred
  component diagnostics. It supports implementation validity but is not independent field
  truth.
- **Extended Data Figure 3 — Circularity sensitivity.** Sensitivity of field interpretation to inputs that
  can enter both candidate construction and scoring. It is retained as an explicit
  anti-circularity diagnostic.

## Supplementary figures
All supplementary figures use the same Nature-style theme (enlarged fonts, collected/
aligned bottom legends) and are composed as multi-panel figures where more than one view
adds value.

- **S1** Hydrochemical context — 2×2: (a) Gibbs cation ratio, (b) Gibbs anion ratio,
  (c) major-ion composition, (d) mineralisation (TDS) by dataset.
- **S2** Charge balance — 2 panels: (a) CBE distribution (±5/±10% guides), (b) quality-class
  composition by dataset.
- **S3** Seasonal signal — 2 panels: (a) dry/wet evapoconcentration factor by region,
  (b) reactive residual after transport correction by region (no independent aquifer-type
  classification exists for these boreholes).
- **S4** Alternative edge networks — 2 panels: (a) inferred process-network composition per
  edge set, (b) divergence (TVD) and stability versus the chemistry-kNN reference set (three
  Hydrosheaf-generated edge sets compared; an earlier revision also compared against a
  fourth, imported "provided graph" edge set, since removed — `DECISIONS.md`).
- **S5** Dominant-process composition across evidence tiers (single heatmap).
- **S6** Diagnostic distributions — 2 panels: (a) bootstrap support stability, (b) MRS.
- **S7** Reactive residual versus conservative evapoconcentration (single, coloured by process).
- **S8** MRS versus support stability (single, coloured by process).
- **S9** External-dataset transfer detail — faceted by dataset (Talensi, Lower Anayari).
- **S10** Identifiability composition per process family and dataset — faceted by dataset.
- **S11** Competing no-flow explanation sensitivity — fractions of evaluated
  edges exceeding the 0.5 and 0.8 HydroSheaf screening thresholds by dataset
  and candidate-edge construction. These are sensitivity flags, not calibrated
  field probabilities.
