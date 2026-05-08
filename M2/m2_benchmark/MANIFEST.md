# M2 Benchmark Manifest

This manifest describes the repository-curated M2 benchmark package. Bulky raw
public inputs are intentionally excluded from Git and can be restored from the
DOIs, URLs, README files, and source manifests recorded under `external/`.

## Committed Core Benchmark Assets

- `config/ground_truth.yaml`
- `data/baseline_samples_realisation_000.csv`
- `data/ground_truth_edges.csv`
- `data/ground_truth_nodes.csv`
- `data/hydrosheaf_reaction_dictionary.csv`
- `data/realisations/*.csv`
- `results/*.csv`
- `results/e3_phreeqc_forward_report.md`
- `tables/*.csv`
- `figures/*.png`
- `docs/02_figures.md`
- `docs/03_tables.md`
- `docs/external_validation_plan.md`
- `docs/m2_results_summary.md`
- `scripts/*.py`

## Committed External-Validation Evidence

- `external/usgs_age/results/*.csv`
- `external/usgs_age/results/*.md`
- `external/usgs_age/results/*manifest.json`
- `external/modpath/results/*.csv`
- `external/modpath/results/*.md`
- `external/modpath/results/*manifest.json`
- `external/phreeqc/input/README.md`
- `external/phreeqc/results/README.md`
- `external/dgmeta/results/*.csv`
- `external/dgmeta/results/*.md`
- `external/dgmeta/results/*manifest.json`
- `external/northern_ghana/results/*.csv`
- `external/northern_ghana/results/*.md`
- `external/northern_ghana/results/*manifest.json`
- `external/usgs_public_chem/results/*.csv`
- `external/usgs_public_chem/results/*.md`
- `external/usgs_public_chem/results/*manifest.json`

## Excluded Raw Inputs

The following are not committed because they are large, public, local, or
reproducible from source:

- downloaded USGS zip archives;
- extracted USGS public-source tables;
- MODPATH model archives and raw endpoint/pathline files;
- DGMETA workbook macros;
- superseded local pilot input files;
- M2 draft Word documents and duplicate top-level source downloads.

## Claim Guardrail

M2 can claim an integrated Hydrosheaf validation package only when captions and
text retain the tier-specific limitations in `docs/m2_results_summary.md`.
External validation rows should not be collapsed into a single claim of full
age, topology, and reactive-transport equivalence.
