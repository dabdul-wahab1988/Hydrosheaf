# M2 External Validation Plan

This plan records what the M2 outline and table/figure guide require beyond the
locked synthetic benchmark already generated in this directory.

## What the M2 Documents Require

The M2 outline requires external validation in three Results subsections:

- Section 3.3: nuclear-tracer and residence-time validation.
- Section 3.4: MODFLOW/MODPATH topology validation.
- Section 3.5: PHREEQC-constrained hydrogeochemical validation.

`TablesFigures.docx` makes this operational in Table 4 and Figure 4:

- Table 4 must include public tracer-age datasets, MODFLOW/MODPATH topology
  comparison, and PHREEQC forward kinetic validation.
- Figure 4A must compare Hydrosheaf residence-time estimates with published
  TracerLPM/LUMPY-style ages from public datasets.
- Figure S2 should include a worked MODPATH-to-graph conversion example from a
  public model archive.
- Figure S3 should show live PHREEQC forward-validation diagnostics, not only a
  mass-balance proxy.

## External Validation Tiers To Add

### Tier E1: Public Groundwater-Age Validation

Primary dataset:

- USGS data release: `Data for Distribution of Groundwater Age in Aquifers Used
  for Public Supply in the Continental United States, 2004-2017`, DOI
  `10.5066/P9W7T0DN`.
- Access URL:
  `https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004`

Why this is the right dataset:

- It includes environmental tracer measurements and LPM results for 1,279 sites
  across 21 principal aquifers.
- It includes tritium, tritiogenic helium-3, sulfur hexafluoride, carbon-14,
  radiogenic helium-4, and TracerLPM age-distribution outputs.

Hydrosheaf task:

- Ingest the seven USGS tables.
- Map site metadata, tracer measurements, and reference LPM outputs into
  Hydrosheaf's residence-time validation schema.
- Use the canonical M3 age-fraction-constrained parity benchmark for samples
  with enough tracer and reported TTD-fraction data:

  `python M3\m3_age_benchmark\scripts\run_m3_design_matrix.py --full --age-steps 90 --scenario tracerlpm_parity_agefractions --output M3\m3_age_benchmark\results\m3_tracerlpm_parity_agefractions_full.csv`

- Compare Hydrosheaf age estimates against published LPM mean age and young /
  Holocene / Pleistocene fractions.

Required outputs:

- `external/usgs_age/input/usgs_age_sites.csv`
- `external/usgs_age/results/usgs_age_validation.csv`
- `external/usgs_age/results/usgs_age_validation_summary.csv`
- `M3/m3_age_benchmark/results/m3_tracerlpm_parity_agefractions_full.csv`
- `figures/Manuscript_Ready/Manuscript_Fig5_Residence_Time_Validation.png`
- `figures/Manuscript_Ready/Manuscript_Supp_FigS1_Public_Age_Validation.png`

Minimum metrics:

- log10 age RMSE.
- median bias in log10 years.
- fraction of estimates inside the published uncertainty range where available.
- stratified metrics by tracer type and aquifer group.

### Tier E2: MODFLOW/MODPATH Topology Validation

Primary candidate dataset:

- USGS data release: `MODFLOW-2005, MODPATH, and MOC3D used for groundwater flow
  simulation, pathlines analysis, and solute transport in the crystalline-rock
  aquifer in the vicinity of the Savage Municipal Water-Supply Well Superfund
  Site, Milford, New Hampshire`, DOI `10.5066/F7J102FK`.
- Access URL:
  `https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute`

Why this is a good first candidate:

- It is explicitly a MODFLOW-2005 and MODPATH5 archive.
- The aquifer is fractured/crystalline, which is closer to a data-limited,
  pathway-uncertain setting than an ideal basin model.
- The release contains input and output files, so it can support a reproducible
  pathline/endpoints-to-graph workflow.

Hydrosheaf task:

- Download the model archive and identify endpoint/pathline files.
- Use `hydrosheaf.physics.modpath` or a small adapter to parse MODPATH endpoints.
- Define graph nodes from wells, source zones, discharge cells, or receptor
  groups in the model archive.
- Convert particle paths into directed edge priors with edge probability and
  travel-time summaries.
- Compare Hydrosheaf-inferred graph edges against MODPATH-derived edges.

Required outputs:

- `external/modpath/input/modpath_endpoints.*`
- `external/modpath/results/modpath_graph_priors.csv`
- `external/modpath/results/modpath_topology_agreement.csv`
- `figures/Manuscript_Ready/Manuscript_Fig2_Topology_Validation.png`

Minimum metrics:

- edge true-positive, false-positive, and false-negative counts.
- direction agreement rate.
- source-receptor path overlap.
- travel-time consistency where MODPATH travel times are available.

### Tier E3: Live PHREEQC Forward Validation

Primary source:

- USGS PHREEQC Version 3 distribution, database files, and examples.
- USGS PHREEQC Techniques and Methods 6-A43 examples.
- Access URLs:
  `https://www.usgs.gov/software/phreeqc-version-3`
  `https://pubs.usgs.gov/tm/06/a43/`

Hydrosheaf task:

- Configure a live PHREEQC executable or robust `phreeqpython` backend.
- Replace the current `linear_mass_balance_phreeqc_proxy` rows with real PHREEQC
  speciation/kinetic runs.
- Run forward simulations from upstream chemistry plus Hydrosheaf-inferred
  reactions and compare against downstream chemistry.

Required outputs:

- `external/phreeqc/results/phreeqc_forward_validation.csv`
- `external/phreeqc/results/phreeqc_failure_log.csv`
- `figures/Manuscript_Ready/Manuscript_Supp_FigS2_Geochemical_Validation.png`

Minimum metrics:

- major-ion RMSE.
- saturation-index residuals.
- NSE or equivalent goodness-of-fit.
- feasibility pass/fail flag by pathway.

### Tier E4c: Northern Ghana Field-Hydrochemistry Demonstration

Selected dataset:

- `data/FieldData/NorthenGhana/NorthernGhana.xlsx`, with optional supplementary methods
  in `data/FieldData/NorthenGhana/SI.pdf`.

Hydrosheaf task:

- Treat this as field-hydrochemistry and data-limited workflow evidence, not as
  replacement for E1 tracer-age, E2 MODPATH topology, or E3 live PHREEQC truth.
- Use charge-balance screening to select quantitative wet/dry samples.
- Because the corrected workbook has no graph-edge sheet, use Hydrosheaf's
  probabilistic graph inference mechanism to build directed graph priors from
  coordinates and the `Elevation_m - Static_Water_Level_m` head proxy.
- Harmonize major ions from mg/L to mmol/L and run Hydrosheaf edge fitting by
  season.
- Report what Hydrosheaf can infer from public/citable field chemistry, while
  retaining the guardrail that no independent process-truth graph is available.

Required outputs:

- `external/northern_ghana/results/northern_ghana_samples.csv`
- `external/northern_ghana/results/northern_ghana_edge_results.csv`
- `external/northern_ghana/results/northern_ghana_validation_summary.csv`
- `external/northern_ghana/results/e4c_northern_ghana_report.md`
- `figures/Manuscript_Ready/Manuscript_Fig4_Ghana_Process_Network.png`
- `figures/Manuscript_Ready/Manuscript_Supp_FigS3_Ghana_Field_Residuals.png`

## Correction To Current M2 Package

The generated synthetic benchmark is valid for Section 3.2 and part of Table 4.
The public-age claim should now be sourced from the M3 age-fraction-constrained
USGS parity benchmark, while the MODPATH result should be described as
topology-only agreement. PHREEQC remains a proxy unless a live PHREEQC backend
is configured and rerun. The Ghana field network remains a field-hydrochemistry
demonstration because no independent process-truth graph is available.
