# M5 Data Readiness Audit for the Proposed Nature-Portfolio Manuscript

Audit date: 15 June 2026

> **Update (2026-07-26): the "richer workbook" (`Aquifers_Dataset_Mendeley.xlsx`)
> has been retired from the M5 pipeline.** It was found to be a derived
> product of a different, antecedent study ("Graph-inverted Ghanaian aquifers
> under aridification"; see `data/FieldData/NorthenGhana/SI.pdf`) that
> classified these same boreholes and pre-computed its own connectivity
> graph, not part of this project's own field data, and it is not present in
> `data/FieldData/`. All references below to the "richer workbook", its
> 294-record quantitative subset, aquifer/geology/lithology metadata, and its
> 397 graph edges describe a source that M5 no longer reads. M5 now uses
> only the raw `data/FieldData/NorthenGhana/NorthernGhana.xlsx` Dry/Wet
> sheets (all 160 boreholes, 320 records, no aquifer-type/lithology/SI
> metadata, no per-record sampling-date column). The accompanying
> `data/FieldData/NorthenGhana/SI.pdf` (Supplementary Method 1 and
> Supplementary Table 19) documents the sampling design for these same 160
> boreholes: a dry-season campaign of 5-24 March 2025 and a wet-season
> campaign of 5-24 August 2025, one record per borehole per season. See
> `M5/m5_inverse_reaction_benchmark/DECISIONS.md`.

## Executive decision

The repository contains enough data and code to build a strong
identifiability-centred M5 paper for a journal such as *Communications Earth &
Environment*, provided that the present M5 smoke test is replaced by a larger
live-PHREEQC factorial benchmark.

The repository does not yet support a *Nature Water*-level claim of jointly
validated field connectivity, residence time, and reaction mechanism. The
Northern Ghana data can support a chemistry-only field demonstration, not
ground-truth validation of reactions along known flow paths.

## Readiness matrix

| Proposed analysis | Evidence currently available | Readiness | Defensible use now | Required addition |
|---|---|---:|---|---|
| Stoichiometric null space and reaction-equivalence classes | The M2 dictionary has 11 ions and 12 reactions, rank 8, nullity 4, infinite condition number, and mutual coherence approximately 1. Exact dependencies include opposite ion-exchange directions. | Ready | Establish that low residual and sparse support do not imply a unique mechanism. | Implement SVD/null-space reporting and group indistinguishable reactions into equivalence classes. |
| Hydrogeochemical identifiability phase diagram | M2 contains 100 realisations, 10 true edges, four topology variants, four missing-data scenarios, and 7,000 reaction-recovery rows. | Partly ready | Produce a preliminary phase diagram for the existing locked aquifer design. | Generate a factorial set of distinct water chemistries, reaction assemblages, noise levels, missing-ion panels, and reaction sparsities using live PHREEQC. |
| Mechanism Resolution Score | M2 supplies known reaction labels and errors; M6 supplies 40 regularisation-path rows, 34 phase-stability rows, and six leave-one-mineral diagnostics. | Pilot ready | Prototype and internally calibrate a score combining rank, coherence, support stability, held-out error, and thermodynamic feasibility. | Expand synthetic archetypes and reserve entire archetypes as an untouched test set before claiming general calibration. |
| Predictive falsification by held-out ions | M2 and the Ghana dataset contain sufficiently complete major-ion vectors. | Ready after code change | Test whether inferred reactions predict ions excluded from fitting. | Replace the current zero-weight missing-ion comparison with true held-out prediction error and predeclare train/test ion panels. |
| Next-best chemical measurement | Complete Ghana records include major ions, water isotopes, Sr, and SiO2; the synthetic package also includes redox and nitrate-isotope variables. | Retrospective simulation ready | Hide measured variables, rank the next measurement by expected ambiguity reduction, then test against the known hidden value. | Add a reaction-specific expected-information-gain routine. Existing active learning is mainly topology oriented. Prospective field validation remains optional for CE&E but important for Nature Water. |
| Edge-conditioned chemical confidence | M2 has topology truth, topology variants, reaction truth, and known ages. M4 independently benchmarks topology against MODPATH archives. | Synthetic ready; field incomplete | Quantify how reaction confidence changes when topology probability is uncertain in synthetic experiments. | Preserve `edge_confidence` in Ghana result exports and obtain independent connectivity labels before making field-validation claims. M4 and Ghana are not paired observations. |
| Residence-time-conditioned reaction plausibility | M2 has true synthetic mean residence times. M3 has 1,272 public USGS age-model records. | Synthetic/module validation only | Demonstrate the concept synthetically and show that Hydrosheaf can ingest age uncertainty. | Obtain age tracers for the same Ghana wells or use an external dataset containing both chemistry and age tracers. M3 sites cannot be joined to Ghana wells. |
| Field hydrochemical demonstration | Northern Ghana provides 160 wells and 320 wet/dry records; 294 pass the quantitative charge-balance class. The richer workbook contains four aquifer types and 12 lithologies. | Date correction completed; provenance citation pending | Show field support stability, alternative pathways, held-out-ion skill, and chemistry-only mechanism classes. | Cite the final public data source and state that graph edges and reaction truth are not independently observed. |
| External field transfer test | Lower Anayari has 41 records, Talensi has 63, and another georeferenced chemistry workbook has 41 records. | Potentially useful | Use one or more datasets as locked external chemistry-transfer cases. | Document provenance, units, sampling design, lithology, and independent publication or repository citation. None currently supplies reaction truth or groundwater ages. |
| Live PHREEQC forward validation | M5 contains one small PHREEQC-generated example. M2 contains 1,000 forward-validation rows. | Not manuscript ready | Retain only as development evidence. | Correct M5 species/unit handling and run live PHREEQC. M2 is explicitly labelled `linear_mass_balance_phreeqc_proxy`, not live PHREEQC. |

## Evidence inventory

### M5 prototype

- One deterministic upstream/downstream pair.
- Seven ions, five reactions, five L1 penalties, and four missing-ion tests.
- The matrix is full rank in this small example, so it is poorly designed to
  demonstrate the broader dictionary's non-identifiability.
- The two current PHREEQC result files are inconsistent: one reports RMSE
  0.237 with two selected phases, while another reports RMSE 553.912 with four
  phases.
- Missing-ion tests set ion weights to zero but do not evaluate prediction of
  the excluded ions.

Conclusion: retain M5 as a software smoke test only.

### M2 synthetic benchmark

- 100 realisations of one locked aquifer design.
- 9,800 edge-fit rows and 7,000 reaction-recovery rows.
- 2,100 rows contain non-zero true reaction extents, while 1,729 rows are
  marked as false-positive recovery events.
- Four topology variants: full, sparse, dense, and reversed.
- Four data scenarios: complete, ion-incomplete, tracer-absent, and
  head-absent.
- Known transport processes, reaction extents, topology, and mean residence
  times.
- The PHREEQC forward table is a deterministic mass-balance proxy.

Conclusion: this is the main immediately usable training and pilot-validation
dataset, but it represents only one hydrogeochemical design.

### M6 robustness outputs

- 40 regularisation-path rows.
- 34 mineral phase-stability estimates.
- Six leave-one-mineral structural tests.

Conclusion: these outputs can seed the Mechanism Resolution Score, but their
sample size is too small for independent calibration.

### Northern Ghana field data

- 160 boreholes and 320 seasonal chemistry records.
- Complete major ions, d18O, d2H, Sr, and SiO2 in the corrected workbook.
- 294 records meet the quantitative inverse-modelling charge-balance class.
- The richer workbook includes aquifer, geology, lithology, saturation
  indices, and 397 processed graph edges.
- The corrected workbook supplies no source graph, tracer age, MODPATH truth,
  or reaction truth. Hydrosheaf generated 1,019 candidate edges.
- All 1,019 current Ghana edge fits have `phreeqc_ok = False`,
  `rt_validated = False`, no populated uncertainty method, and no populated
  exported edge confidence.
- Hydraulic residence-time estimates exist for 1,015 edges, but their median
  is approximately 1.74 million days and the maximum is approximately
  3.14 billion days. They are unvalidated hydraulic heuristics and must not be
  presented as measured groundwater ages.

Conclusion: use as a field hydrochemistry demonstration only.

### M3 and M4 supporting modules

- M3 contains 1,272 USGS age-model records across 20 study units and six
  aquifer groups. These validate the age module separately, not Ghana ages.
- The M3 manifest reports 5,155 graph-edge rows, but the expected
  `m3_real_usgs_graph_edges.csv` and
  `m3_real_usgs_graph_benchmark.csv` files are currently absent from the
  results directory.
- M4 contains independent topology experiments and external MODPATH archives.
  The strongest independent head-gradient scenario has F1 approximately
  0.618; MODPATH archive conversion is topology-only evidence.

Conclusion: M3 and M4 strengthen the modular Hydrosheaf thesis, but cannot be
merged with Ghana chemistry as if they were measurements from the same wells.

## Data-integrity issues and resolutions

1. Corrected on 15 June 2026: all 320 `Sampling_Date` values in
   `Aquifers_Dataset_Mendeley.xlsx` now record the March and August 2025
   campaigns.
2. Corrected on 15 June 2026: Supplementary Table 19 in `SI.pdf` now lists the
   exact dry- and wet-season dates in March and August 2025, consistent with
   the date ranges and explanatory prose.
3. The corrected `NorthernGhana.xlsx` does not supply sampling dates, and its
   source manifest states that the final public DOI or URL is not embedded.
4. The 397 graph edges in the richer Ghana workbook are processed analytical
   outputs, not independently observed flow paths.
5. High chemistry R2 on generated edges is in-sample consistency, not proof of
   physical connectivity or correct reactions.
6. Current field outputs do not contain live PHREEQC or tracer-age validation.

The remaining unresolved items, particularly provenance and evidence
independence, must be addressed in the repository record and manuscript before
submission.

## Recommended analysis programme

### Stage 1: executable with existing data

1. Implement null-space, singular-value, condition-number, and mutual-coherence
   diagnostics for every reaction dictionary and ion panel.
2. Define reaction-equivalence classes and score recovery at both exact-phase
   and equivalence-class levels.
3. Add true held-out-ion prediction and compare it with residual-only model
   selection.
4. Prototype the Mechanism Resolution Score using M2 truth and M6 stability
   outputs.
5. Run retrospective next-best-measurement experiments by masking observed
   Ghana and M2 variables.

### Stage 2: required new computational data

Build a live-PHREEQC factorial benchmark with:

- multiple hydrochemical archetypes rather than one aquifer;
- carbonate, evaporite, silicate, exchange, redox, and mixed assemblages;
- identifiable and deliberately non-identifiable dictionaries;
- several noise levels and missing-ion panels;
- alternative regularisation and comparator methods;
- complete scenario-level ground truth;
- locked calibration and external synthetic test partitions.

This is new simulated data, not a new field campaign.

### Stage 3: field evidence for CE&E

After confirming the final public source and provenance citation, use Northern
Ghana for:

- chemistry-only reaction equivalence and support-stability mapping;
- held-out-ion predictive validation;
- aquifer-stratified transfer analysis;
- comparison of conventional best-fit and identifiability-aware reporting.

Use Lower Anayari, Talensi, or the 41-record geology workbook as a locked
external transfer dataset only after their provenance and units are documented.

### Stage 4: additional evidence for a Nature Water stretch target

Collect matched observations from the same wells or independently supported
flow paths:

- tritium, carbon-14, SF6, or CFC age tracers;
- repeated hydraulic heads and screened intervals;
- lithology and mineralogical abundance;
- redox indicators and dissolved gases where relevant;
- prospective measurements selected by the next-best-measurement algorithm.

Without matched data of this kind, the paper should not claim joint field
validation of topology, age, and reaction mechanism.

## Final assessment

The most defensible high-novelty M5 paper is not “Hydrosheaf fits groundwater
chemistry.” It is:

> Hydrosheaf exposes and quantifies when groundwater reaction mechanisms are
> observationally distinguishable, predictively falsifiable, or confined to
> an equivalence class.

Most of the analytical foundation for that paper already exists. The critical
missing component is a large, corrected, live-PHREEQC benchmark with diverse
ground-truth regimes. New groundwater-age and flow-path measurements are only
essential if the target is raised from a strong Nature-Portfolio methods paper
to a jointly validated *Nature Water* field-mechanism paper.
