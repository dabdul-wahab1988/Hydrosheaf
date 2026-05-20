# Examples

## Python Examples

Run the bundled small network demo:

```bash
python hydrosheaf/examples/demo_small_network.py
```

Run the research objectives demo (requires `hydrosheaf_synthetic_csv/` at repo root and optional PyMC/FloPy extras):

```bash
python examples/demo_research_objectives.py
```

## CLI Example

Run on the bundled CSV files:

```bash
hydrosheaf \
  --samples hydrosheaf/examples/sample_data.csv \
  --edges hydrosheaf/examples/edges.csv \
  --output results.json \
  --format json
```

## Quick QA

Run a single test for sanity:

```bash
python -m pytest tests/test_edge_fit.py::EdgeFitTests::test_edge_fit_evap_and_halite
```

Run the fast regression suite:

```bash
python -m pytest tests/test_regression_small.py
```

Run PHREEQC-specific tests when PHREEQC is installed:

```bash
python -m pytest tests/test_phreeqc*.py
```

Run uncertainty tests (can be slower):

```bash
python -m pytest tests/test_accuracy_uncertainty.py
```

Run temporal regression tests:

```bash
python -m pytest tests/test_accuracy_temporal.py
```

Run 3D flow network accuracy tests:

```bash
python -m pytest tests/test_accuracy_3d.py
```

Run reactive transport validation tests:

```bash
python -m pytest tests/test_accuracy_reactive.py
```

## Data Units Notes

- CLI input units are configurable via `--unit` (default `mmol/L`).
- Research objectives demo expects `*_mg_L` columns in the synthetic CSVs.

## Optional Dependencies

- PyMC + ArviZ enable Bayesian mixing and uncertainty workflows.
- FloPy + MODFLOW/MT3DMS executables enable numerical transport coupling.
- PHREEQC is required for thermodynamic constraints and reactive transport validation.

## Example Outputs

- JSON/CSV outputs come from `--format` and `--output` flags.
- Temporal module can write an extra `--temporal-output` JSON file.
- Vadose module can emit priors or nitrate breakthrough summaries.

## Troubleshooting

- If running examples on Windows, ensure required executables are on `PATH`.
- If PyMC is missing, Bayesian mixing and uncertainty steps are skipped.
- If FloPy is missing, transport coupling uses analytical fallback.

## Edge Inference from Coordinates

If your samples include `lat`, `lon`, and optional `elevation`, you can infer edges:

```bash
hydrosheaf \
  --samples samples.csv \
  --infer-edges \
  --max-neighbors 1 \
  --output results.json
```

Use `--head-key hydraulic_head` when hydraulic head data is available.

## Probabilistic Edge Inference

Use head/DTW when available and fall back to topography with uncertainty:

```bash
hydrosheaf \
  --samples samples.csv \
  --infer-edges \
  --infer-edges-method probabilistic \
  --edge-p-min 0.75 \
  --edge-radius-km 5 \
  --edge-max-neighbors 3 \
  --edge-head-key head_meas \
  --edge-dtw-key dtw \
  --edge-elevation-key elevation \
  --output results.json
```

## EC/TDS Calibration

Fit linear EC/TDS models from the dataset:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --calibrate-ec-tds
```

## Detection Limits

Handle detection limits with half‑limit substitution:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --detection-policy half
```

## Missing Ions

Impute missing ions with zeros (use cautiously):

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --missing-policy impute_zero
```

## PHREEQC Constraints

Enable PHREEQC (if installed) and set a tau threshold:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --phreeqc-enabled --si-threshold 0.2
```

## Endmembers JSON

Load a JSON file of endmembers:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --endmembers-json endmembers.json
```

## Isotope Integration

Fit a local meteoric water line and enable isotope penalties:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --isotope-enabled --fit-lmwl --isotope-weight 1.0 --isotope-d-excess-weight 0.5
```

To use the Gibrilla et al., 2022 LMWL directly:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --isotope-enabled --lmwl-a 8.66 --lmwl-b 7.22 --isotope-weight 1.0 --isotope-d-excess-weight 0.5
```

## Nitrate Source Discrimination

Enable manure vs. fertilizer discrimination using the `--nitrate-source-enabled` flag. This works best when Potassium (K), Chloride (Cl), and Nitrate (NO3) data are present.

### Basic Usage

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --nitrate-source-enabled
```

### With Custom Minimum Concentration Threshold

By default, nitrate source discrimination is only performed on samples with NO3 ≥ 10 mg/L. You can adjust this threshold:

```bash
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --nitrate-source-enabled --nitrate-source-min-conc 5.0
```

This is useful for:
- **Background exclusion**: Avoid attempting source discrimination on ambient/background groundwater
- **Detection limit alignment**: Match your laboratory's reliable quantification limits
- **Focus on contaminated areas**: Concentrate analysis on samples with significant nitrate enrichment

The output will contain:

- `nitrate_source_p_manure`: Probability of manure source (0.0 - 1.0), or `null` if below threshold.
- `nitrate_source_evidence`: Top reasons (e.g., `NO3/Cl_high_manure`, `denitrif_strong`).
- `nitrate_source_reason`: Explanation if discrimination was not performed (e.g., "Low Nitrate (Background < 10 mg/L)").

## Vadose Priors and Transport Coupling

Generate vadose-derived physics priors:

```bash
hydrosheaf \
  --samples samples.csv \
  --edges edges.csv \
  --vadose-enabled \
  --vadose-forcing forcing.csv \
  --vadose-profile profile.yaml \
  --vadose-links links.csv \
  --vadose-priors-output vadose_priors.csv \
  --output results.json
```

Predict nitrate breakthrough to the water table:

```bash
hydrosheaf \
  --vadose-enabled \
  --vadose-forcing forcing.csv \
  --vadose-profile profile.yaml \
  --vadose-links links.csv \
  --vadose-no3-loading no3_loading.csv \
  --vadose-no3-output vadose_breakthrough.csv \
  --vadose-no3-summary-output vadose_summary.json \
  --output results.json
```

## Temporal Module

```bash
hydrosheaf \
  --temporal-data timeseries.csv \
  --temporal-enabled \
  --residence-method cross_correlation \
  --temporal-output temporal_results.json \
  --output results.json
```

## Uncertainty Quantification

```bash
hydrosheaf \
  --samples samples.csv \
  --edges edges.csv \
  --uncertainty bootstrap \
  --bootstrap-n 1000 \
  --output results.json
```

## Reactive Transport Validation

```bash
hydrosheaf \
  --samples samples.csv \
  --edges edges.csv \
  --validate-forward \
  --rt-residence-time 30.0 \
  --rt-rmse-threshold 1.0 \
  --output results.json
```

## 3D Flow Networks

```bash
hydrosheaf \
  --samples wells_3d.csv \
  --output results_3d.json \
  --3d \
  --z-key screen_depth \
  --layer-file layers.yaml \
  --anisotropy 0.1
```

## Optional Hooks

- `--report-only results.json` renders the interpretation report without re-running.
- `--endmember` and `--endmember-id` add mixing endmembers from CLI.
- `--physics-priors` and `--modpath-endpoints` apply travel-time priors.

## Performance Notes

- Bayesian MCMC and reactive transport can take minutes per edge.
- Temporal interpolation with spline requires SciPy.
- PHREEQC subprocess mode is slower but more compatible with custom databases.

## Support Data

- Sample CSVs: `hydrosheaf/examples/sample_data.csv`, `hydrosheaf/examples/edges.csv`.
- Synthetic demo data: `hydrosheaf_synthetic_csv/` at repo root.


## Next Steps

- See `docs/INPUTS_REFERENCE.md` for expected input columns.
- See `docs/TECHNICAL_REFERENCE.md` for modeling details.
- See `docs/UNCERTAINTY_IMPLEMENTATION.md` for runtime considerations.
