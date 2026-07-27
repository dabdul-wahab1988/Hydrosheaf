# External Validation Workspace

This directory is reserved for the M2 external validation tiers required by
`M2outline.docx` and `TablesFigures.docx`.

## Directory Map

- `usgs_age/`: public tracer-age validation using the USGS public-supply aquifer
  groundwater-age release, DOI `10.5066/P9W7T0DN`.
- `modpath/`: MODFLOW/MODPATH topology validation using a public USGS model
  archive, starting with DOI `10.5066/F7J102FK`.
- `phreeqc/`: live PHREEQC forward validation using PHREEQC v3 examples and
  databases, DOI `10.3133/tm6A43`.
- `northern_ghana/`: field-hydrochemistry and data-limited workflow
  demonstration using `data/FieldData/NorthenGhana/NorthernGhana.xlsx`.
- `pilot/`: superseded local pilot workspace retained for provenance only; do
  not use `manu.xlsx` for current M2 evidence.

## Rule For Reporting

Synthetic benchmark outputs belong to Section 3.2. External validation claims
should only be made after the corresponding result files in this directory are
generated from the public or pilot sources listed in
`../docs/external_validation_plan.md`.
