# Examples

## CLI Example

Run on the bundled CSV files:

```
hydrosheaf \
  --samples hydrosheaf/examples/sample_data.csv \
  --edges hydrosheaf/examples/edges.csv \
  --output results.json \
  --format json
```

## Edge Inference from Coordinates

If your samples include `lat`, `lon`, and optional `elevation`, you can infer edges:

```
hydrosheaf \
  --samples samples.csv \
  --infer-edges \
  --max-neighbors 1 \
  --output results.json
```

Use `--head-key hydraulic_head` when hydraulic head data is available.

## Probabilistic Edge Inference

Use head/DTW when available and fall back to topography with uncertainty:

```
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

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --calibrate-ec-tds
```

## Detection Limits

Handle detection limits with half‑limit substitution:

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --detection-policy half
```

## Missing Ions

Impute missing ions with zeros (use cautiously):

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --missing-policy impute_zero
```

## PHREEQC Constraints

Enable PHREEQC (if installed) and set a tau threshold:

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --phreeqc-enabled --si-threshold 0.2
```

## Endmembers JSON

Load a JSON file of endmembers:

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --endmembers-json endmembers.json
```

## Isotope Integration

Fit a local meteoric water line and enable isotope penalties:

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --isotope-enabled --fit-lmwl --isotope-weight 1.0 --isotope-d-excess-weight 0.5
```

To use the Gibrilla et al., 2022 LMWL directly:

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --isotope-enabled --lmwl-a 8.66 --lmwl-b 7.22 --isotope-weight 1.0 --isotope-d-excess-weight 0.5
```

## Nitrate Source Discrimination

Enable manure vs. fertilizer discrimination using the `--nitrate-source-enabled` flag. This works best when Potassium (K), Chloride (Cl), and Nitrate (NO3) data are present.

### Basic Usage

```
hydrosheaf --samples samples.csv --edges edges.csv --output results.json --nitrate-source-enabled
```

### With Custom Minimum Concentration Threshold

By default, nitrate source discrimination is only performed on samples with NO3 ≥ 10 mg/L. You can adjust this threshold:

```
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
