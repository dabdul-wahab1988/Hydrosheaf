# M6 Table Descriptions

## Main tables
- **Table 1 — Dataset readiness and claim strength** (`table1_dataset_readiness`). Per
  dataset: n samples/sites, seasonality, native tier, quantitative/screening/exploratory CBE
  counts, median |CBE|, and the defensible claim strength. Serves RQ1.
- **Table 2 — Variable availability by dataset** (`table2_variable_availability`). Present/
  partial/absent status of 16 variables, showing the F/Sr/SiO₂/metadata gaps that cap the
  external datasets. Serves RQ1.
- **Table 3 — Outputs under full vs reduced information** (`table3_tier_ablation`). Northern
  Ghana identifiability, class/family flips, MRS and stability across Tiers 0–4. Serves RQ3/H1.
- **Table 4 — External transfer performance** (`table4_external_transfer`). MRS, stability and
  identifiability shares for Talensi and Lower Anayari against matched-tier Northern Ghana
  references. Serves RQ5/H4.
- **Table 5 — Field prequential benchmark** (`table5_field_prequential`). Overall
  one-step seasonal hold-forward MAE, RMSE and empirical 90% coverage for persistence,
  expanding-mean and graph-ridge predictors, with the paired well-block bootstrap contrast
  reported in the accompanying audit. This supports chemistry-prediction diagnostics only;
  it does not validate residence time or field topology.

## Supplementary tables
- **S2** charge-balance quality classes by dataset. **S3** Northern Ghana aquifer summary.
  **S4** full tier-ablation identifiability counts. **S5** edge-set sensitivity summary.
  **S6** per-edge external outputs. **S7** uncertainty-metric definitions. **S8** software and
  computational environment.

All tables are emitted as CSV and Markdown by `scripts/make_m6_tables.py`.
