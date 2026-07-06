# M5 Experimental Design

M5 benchmarks Hydrosheaf as an identifiability-aware sparse inverse
hydrogeochemical framework. It does not claim fully unique reaction truth from
major ions alone.

## Experiment 1: Controlled Truth

Purpose: test reaction and equivalence-class recovery against PHREEQC-generated
synthetic truth.

Evidence:

- `results/phreeqc_ground_truth.csv`
- `results/benchmark_fits.csv`
- `results/reaction_recovery.csv`
- `results/phreeqc_inverse_baseline.csv`
- `tables/table1_comparative_inverse_performance.csv`

Methods compared: bounded least squares, lasso, elastic net, thermodynamic
elastic net, Hydrosheaf Guarded, Hydrosheaf-Core, and conventional PHREEQC
inverse modelling.

## Experiment 2: Identifiability

Purpose: quantify what can and cannot be resolved under sparse major-ion
chemistry.

Evidence:

- `results/identifiability_diagnostics.csv`
- `results/equivalence_classes.csv`
- `results/mechanism_resolution_scores.csv`
- `results/phreeqc_inverse_baseline_models.csv`

Main interpretation: low residuals do not imply unique mechanisms; rank,
nullity, equivalence classes, and PHREEQC inverse-model multiplicity define the
safe interpretive boundary.

## Experiment 3: Data Tiers

Purpose: estimate the marginal value of additional diagnostic measurements.

Evidence:

- `results/data_tier_experiment.csv`
- `results/data_tier_reaction_evidence.csv`
- `results/data_tier_evidence_lifted_resolution.csv`
- `results/data_tier_optional_diagnostics.csv`
- `tables/tableS14_data_tier_experiment.csv`
- `tables/tableS15_data_tier_reaction_evidence.csv`
- `tables/tableS16_evidence_lifted_resolution.csv`

Tiers:

- `core`: major ions plus sparse Hydrosheaf evidence gates.
- `plus_lite`: Core plus controlled synthetic SiO2, Sr, and water-isotope
  diagnostics.
- `enhanced`: Plus-lite plus controlled synthetic Br, DO, DOC, sulphate-isotope,
  and nitrate-isotope diagnostics.

The optional diagnostics are generated from controlled reaction truth with
measurement noise. They are a measurement-design experiment, not field-measured
validation.

Evidence-lifted resolution index (ELRI) is reported for ambiguous equivalence
classes. ELRI is entropy based: zero means evidence does not separate class
members, and values approaching one mean one member has dominant evidence
support. It measures conditional evidence separation, not new stoichiometric
uniqueness.

## Experiment 4: Field Transfer

Purpose: test whether the workflow remains useful in a real data-sparse region.

Evidence:

- `results/ghana_field_pairs.csv`
- `results/ghana_field_reaction_extents.csv`
- `results/ghana_field_hydrosheaf_core_evidence.csv`
- `results/external_field_transfer_pairs.csv`
- `results/external_field_evidence_lifted_resolution.csv`
- `tables/tableS10_northern_ghana_summary.csv`
- `tables/tableS13_ghana_hydrosheaf_core_evidence.csv`
- `tables/tableS17_external_field_evidence_lifted_resolution.csv`

Main interpretation: the Ghana component demonstrates plausible sparse-data
transfer and evidence auditing. It is not independent reaction-truth validation.
The external field-transfer extension repeats the ELRI audit on
`NorthernGhana.xlsx`, Talensi, and Lower Anayari chemistry using wet-dry or
nearest-neighbour field edges.
