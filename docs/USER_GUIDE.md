# Hydrosheaf User Guide

## 1. Quick Start

### Basic Run

Run the model on your prepared dataset using default minerals:

```powershell
python -m hydrosheaf.cli --samples data_prepared.csv --output results.json --infer-edges --phreeqc-enabled
```

### Advanced Run (Expert Mode)

Enable isotopes, auto-calibration, and interpretation report:

```powershell
python -m hydrosheaf.cli \
  --samples data_prepared.csv \
  --output results_final.json \
  --infer-edges \
  --phreeqc-enabled \
  --isotope-enabled \
  --auto-lmwl \
  --iso-consistency-weight 2.0 \
  --interpret
```

---

## 2. CLI Options Reference

### Data Inputs

* `--samples`: Path to CSV file containing ion concentrations (Site ID, Ca, Mg, Na... etc).
* `--unit`: Input units. Options: `mmol/L` (default), `mg/L`, `meq/L`.

### Geochemical Settings

* `--minerals`: Comma-separated list of allowed minerals.
  * *Example*: `"calcite,dolomite,gypsum,halite,albite,pyrite_oxidation_aerobic"`
* `--list-minerals`: Print all valid minerals in the library and exit.
* `--phreeqc-enabled`: Check saturation indices to constrain dissolution/precipitation.

### Isotope Settings

* `--isotope-enabled`: Turn on isotope physics.
* `--auto-lmwl`: Automatically fit the Local Meteoric Water Line from your data.
* `--iso-consistency-weight`: Weight for the penalty connecting isotope enrichment to Chloride increase. Set to `2.0` for strict physics.

### Interpretation & Output

* `--interpret`: Print a text-based geochemical report after the run finishes.
* `--report-only <file.json>`: Generate a report from an existing results file without re-running the model.

---

## 3. Interpreting the Report

The interpretation report provides high-level insights:

### Transport Physics

* **Evap-Dominant**: Water is concentrating via evaporation.
* **Mix-Dominant**: Water is mixing with a saline end-member (e.g., seawater).

### Isotope Consistency

* **Enrichment Slope**: ~3.5-6.0 indicates strong evaporation. ~8.0 indicates rain/recharge.
* **Consistency Penalty**: High values warn that chemical salinity is increasing faster than physical evaporation can explain (Process: **Rock Interaction**).

### Redox State

* **Oxic**: High Nitrate, Low Iron. Pyrite burns to Sulfate.
* **Reducing**: Low Nitrate, High Iron. Pyrite is stable; Aerobic oxidation is disabled.
