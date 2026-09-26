# Hydrosheaf

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18157915.svg)](https://doi.org/10.5281/zenodo.20339942)

**Hydrosheaf: A Graph-Sheaf Framework for Assumption-Audited Groundwater Inference across Hydrogeochemistry, Tracer Ages, and Flow-Topology Uncertainty**

## Overview

Hydrosheaf is a Python-based framework for groundwater inference across hydrogeochemistry, environmental tracers, and flow-network topology. It combines inverse geochemical modeling, nuclear-tracer age estimation, graph-based aquifer connectivity, uncertainty quantification, and benchmark workflows against public groundwater-age and MODFLOW/MODPATH reference datasets.

The current public source branch is intentionally source-only: it contains the Python package, runtime configuration, tests, and reproducibility/benchmark runners needed to inspect and exercise the implementation. Manuscripts, raw and derived research data, generated figures, and historical output trees remain local and are not treated as package inputs or public validation evidence.

The framework integrates:

- **Weighted Least Squares Optimization** for transport model selection.
- **Sparse LASSO Regression** with coordinate descent for parsimonious reaction fitting.
- **Thermodynamic Constraints** utilizing PHREEQC for saturation index calculations.
- **Field-Data Contracts and Unit Registries** for manifest-backed datasets, chemistry normalization, and configured trace-species/speciation inputs such as boron, silica, and arsenic.
- **Isotope and Nuclear-Tracer Hydrogeology** for process validation and groundwater-age inference.
- **Graph and Sheaf-Based Network Inference** for probabilistic flow connectivity and topology refinement.
- **Null-Model Screening** for ruling out non-connectivity explanations (shared lithology, endmember similarity, spatial proximity).
- **Bayesian Topology Posterior** for quantified edge-inclusion probabilities and uncertainty via Metropolis-Hastings sampling.
- **Sheaf Cohomology Diagnostics** for detecting global flow-consistency obstructions (cycles where chemistry constraints cannot be simultaneously satisfied).
- **Optimal Transport and Causal Discovery** for reaction-aware chemistry-plausibility screens and guarded causal direction support.
- **Active Learning** for recommending which wells to measure next based on topology uncertainty and validation gaps.
- **Certified Measurement Design** for bounded tracer and conditional hydrochemical models, with explicit ambiguity witnesses, sequential updates, cost accounting, and abstention when declared evidence is inadequate.
- **MODFLOW/MODPATH Benchmarking** for testing reduced-order graph topology against reference particle-tracking outputs.
- **History-aware TTD inversion** for stable-water-isotope and atmospheric-tracer response matrices conditioned on dated source histories, with explicit abstention reasons for missing or invalid histories.
- **Truth-blind topology-v2** for all-pairs candidate generation, soft physical evidence, optional calibration, tri-state `PRESENT`/`ABSENT`/`ABSTAIN` decisions, and bootstrap case metrics.
- **Fair M4 benchmark modes** that separate sparse geometry, archive-informed, and path-aware evidence while keeping MODPATH reference edges out of inference.
- **Reproducible execution** through deterministic seeds/epochs, file hashes, tree comparison, isolated output staging, and recipe-level replay checks.

The implementation distinguishes software contracts and controlled benchmark evidence from independent field or system-level validation. Inference routines fail closed or return `ABSTAIN` when required evidence is missing, contradictory, or outside the declared scope.

## Features

- **Transport Modeling**: Distinguishes between evaporation (Rayleigh distillation-like) and end-member mixing.
- **Reaction Path Modeling**: Identifies mineral dissolution/precipitation, redox reactions (denitrification), and ion exchange.
- **Sparsity**: Uses L1 regularization to find the simplest chemical explanation for observed data.
- **Physical Consistency**: Enforces thermodynamic bounds (e.g., minerals cannot precipitate from undersaturated solutions).
- **Network Inference**: Infers flow direction probabilities from hydraulic head data with uncertainty.
- **Nitrate Source Discrimination**: Distinguishes between manure and fertilizer sources using a hybrid approach: Dual Isotope Bayesian Mixing ($\delta^{15}\text{N}, \delta^{18}\text{O}$), with optional boron evidence, prioritized over CoDA-based hydrochemical statistics.
- **Reactive Transport**: Validates inverse results against kinetic rate laws (Arrhenius) and Damköhler numbers.
- **3D Flow Networks**: Analyzes layered aquifer systems with vertical anisotropy, geophysical observations, bedrock barriers, transit-time effects, and topographic priors.
- **Temporal Dynamics**: Resolves time-variant signals with cross-correlation residence times and seasonal decomposition.
- **Conditional Time-History TTDs**: Uses supplied dated source histories without silently extrapolating stable-isotope inputs; results are explicitly conditional on those histories.
- **Uncertainty Quantification**: Provides rigorous confidence intervals via Bayesian MCMC (NUTS) and Bias-Corrected Bootstrap (BCa).
- **Assumption Auditing**: Null-model screening rules out alternative explanations before accepting flow connectivity, with an evidence ladder from FALSIFIED to VALIDATED per edge.
- **Topology-v2 Uncertainty**: Truth-blind all-pairs candidate universes, calibrated scores, frozen threshold policies, missingness-aware feature rows, and tri-state decisions.
- **Fair Topology Benchmarking**: M4-A/M4-B/M4-C workflows using projected Savage coordinates, public MODFLOW CBC context, and path-aware physical evidence.
- **Global Consistency Checks**: Sheaf cohomology detects cycles where chemistry constraints cannot be simultaneously satisfied, computing obstruction energy and per-edge leverage scores.
- **Active Learning**: Recommends which wells to measure next based on variant disagreement, posterior uncertainty, geophysical ambiguity, missing tracer/boron evidence, and validation gaps.
- **Conditional Hydrochemical Design**: Uses declared reaction stoichiometry, finite extent bounds, concentration intervals, and site-specific measurement responses to select measurements under a robust minimax criterion. Its certificates apply only to that stated linear model.
- **Reproducibility Contracts**: Deterministic environment setup, SHA-256 manifests, exact tree comparisons, and isolated rerun verification.
- **Contract-Tested Source**: Public tests and configuration exercise package imports, validation gates, benchmark rules, and abstention behavior.
- **Evidence-Aware Optional Modules**: Input validation can resolve available topology, hydraulic, geophysical, and tracer modules from the signals actually present instead of silently fabricating missing evidence.

## Installation

### Core Package (CLI)

The current package metadata declares version `0.7.0` and requires Python
`3.10` or newer.

```bash
git clone https://github.com/dabdul-wahab1988/Hydrosheaf.git
cd Hydrosheaf
pip install .
```

Optional 3D visualization and VTK export support can be installed with:

```bash
pip install ".[viz3d]"
```

### External Dependencies

The core package and its Python validation workflows run from the declared Python
dependencies. External simulator workflows are separate adapters and require the
corresponding executables and model files to be installed and configured by the
user:

1. **Core inverse modeling** works within Python and does not require PEST++,
   MODFLOW, or MT3DMS executables.
2. **PEST++/MODFLOW/MT3DMS workflows** are optional and depend on the external
   binaries, templates, and input files required by the selected calibration or
   transport model.
3. The PEST++ adapter can attempt to download its configured PEST++ release
   (currently defaulting to version `5.2.25`) when the executable is missing;
   this requires network access and writes the runtime binary under the local
   `bin/<version>/` directory. MODFLOW and MT3DMS executables are not managed by
   this path. For reproducible runs, preinstall or pin external tools and record
   their versions in the run provenance.

## Quick Start

### Core Python Workflows

After `pip install .`, the public package provides the core Python workflows
for:

- Core inverse geochemical modeling (transport + reactions)
- PHREEQC thermodynamic constraints
- Isotope analysis and forensics
- Nitrate source discrimination (Bayesian MCMC)
- Network inference and topology refinement
- Assumption auditing (null models + evidence ladder)
- Topology uncertainty quantification (Bayesian posterior)
- Sheaf cohomology global-consistency diagnostics
- Optimal transport and causal direction screening
- 3D flow-network inference
- Temporal dynamics and residence time estimation
- History-aware TTD inference with conditional source histories and explicit abstention
- Truth-blind topology-v2 scoring and associated benchmark utilities
- Uncertainty quantification (Bootstrap, MCMC)
- Active learning measurement recommendations

Install `.[viz3d]` only when you need PyVista/VTK-based 3D plotting or VTK file export. Benchmark and replay commands are designed to write to isolated run directories and do not overwrite preserved historical output trees.

### Command-Line Entry Points

The package installs two CLI entry points:

```bash
hydrosheaf --help
hydrosheaf-cal --help
```

The public benchmark and reproducibility runners can be inspected directly from
the repository root:

```bash
python scripts/run_topology_v2_benchmark.py --help
python M4/m4_topology_benchmark/scripts/run_m4_fair_topology_benchmark.py --help
python scripts/reproduce_outputs.py --help
```

These runners are deliberately conservative. The topology-v2 runner writes to
an isolated `.codex_work/topology-v2` directory by default; the fair M4 runner
keeps MODPATH reference edges in the evaluator rather than the inference path;
and `reproduce_outputs.py` stages reruns separately from historical `outputs`
trees and can compare hashes in strict mode. A complete run still requires the
corresponding public/archive inputs to be available locally.

### Python API Usage

The recommended entry point is `fit_network_pipeline()` from `hydrosheaf.api`:

```python
from hydrosheaf.api import fit_network_pipeline
from hydrosheaf.config import Config

# Load your data
samples = {...}  # List of sample dictionaries
edges = [...]    # List of edge definitions
config = Config()

# Run the pipeline
results, diagnostics = fit_network_pipeline(
    samples=samples,
    edges=edges,
    config=config,
    auto_disable_missing=True
)

# Access results
for edge_result in results:
    print(f"Edge {edge_result.u} → {edge_result.v}:")
    print(f"  Transport model: {edge_result.transport_model}")
    print(f"  Reactions: {edge_result.z_labels}")
    print(f"  Reaction extents: {edge_result.z_extents}")
```

### History-Aware TTD API

The history-aware TTD wrapper is conditional on the dated source histories
provided by the caller. Stable-isotope rows are not silently extrapolated, and
invalid, incomplete, contradictory, or unsupported inputs return an explicit
`ABSTAIN` result with machine-readable reason codes.

```python
from hydrosheaf import fit_history_ttd
from hydrosheaf.nuclear.history_ttd import HistoryTracerObservation

observations = [
    HistoryTracerObservation(tracer="d18O", value=-8.5, sigma=0.2),
    HistoryTracerObservation(tracer="3H", value=4.1, sigma=0.5),
]

result = fit_history_ttd(
    observations=observations,
    sample_year=2020,
    age_grid_years=[0.0, 1.0, 5.0, 10.0, 20.0],
    source_histories={"d18O": dated_d18O_history},
)

if result.status == "ABSTAIN":
    print(result.abstention_reasons)
else:
    print(result.age_grid_years, result.g)
```

This is a model-conditioned inference result, not independent field validation
of groundwater ages, flow paths, or source-history reconstruction.

### Optional External Calibration Tooling

PEST++ source compilation is outside the Python package. If a calibration
workflow needs PEST++, the adapter can attempt to fetch its configured release
when no local executable is available, or you can install/build it separately
and provide the executable and model configuration to the adapter. MODFLOW and
MT3DMS remain separately managed external tools. None of this is required for
the core inverse-modeling, TTD, topology-v2, or reproducibility workflows.

### Configuration

Hydrosheaf uses a central `Config` dataclass. Key settings:

```python
from hydrosheaf.config import Config

config = Config(
    # Chemistry
    ion_order=["Ca", "Mg", "Na", "HCO3", "Cl", "SO4", "NO3", "F", "Fe", "PO4"],
    charge_balance_limit=0.1,
    
    # Inference
    lambda_l1=0.01,  # LASSO penalty
    transport_models_enabled=["evap", "mix"],
    
    # Thermodynamics
    phreeqc_enabled=True,
    si_threshold_tau=0.2,
    
    # Isotopes & Nitrate
    isotope_enabled=True,
    nitrate_source_enabled=True,
    
    # Sheaf topology & evidence auditing
    sheaf_cohomology_enabled=False,
    topology_posterior_enabled=False,
    assumption_calibration_enabled=False,
    evidence_ladder_enabled=False,
    
    # Optimal transport & causal screening
    ot_enabled=False,
    causal_discovery_enabled=False,
    
    # Temporal
    residence_time_method="bayesian_lag",
    
    # 3D Flow
    network_3d_enabled=False,
    vertical_anisotropy=0.1,
)
```

For the installed dependency and CLI contracts, see [pyproject.toml](pyproject.toml),
[hydrosheaf/cli.py](hydrosheaf/cli.py), and
[hydrosheaf/calibration/cli.py](hydrosheaf/calibration/cli.py).

### Data Input Format

Samples should be provided as a list of dictionaries with keys:

```python
samples = [
    {
        "site_id": "W001",
        "sample_id": "W001_2023-01",
        "Ca": 80.0,      # mg/L (converted internally to mmol/L)
        "Mg": 20.0,
        "Na": 15.0,
        "K": 2.0,
        "HCO3": 250.0,
        "Cl": 30.0,
        "SO4": 50.0,
        "NO3": 5.0,
        "F": 1.0,
        "Fe": 0.1,
        "PO4": 0.05,
        "pH": 7.5,
        "d18o": -8.5,    # Optional: δ18O (‰)
        "d2h": -60.0,    # Optional: δ2H (‰)
        "latitude": 40.5,
        "longitude": -74.0,
        "elevation": 50.0,
        "head": 45.0,    # Optional: hydraulic head (m)
    },
    # ... more samples
]
```

Edges should specify connectivity:

```python
edges = [
    {"u": "W001", "v": "W002", "distance_km": 2.5, "delta_h": 1.5},
    {"u": "W002", "v": "W003", "distance_km": 3.0, "delta_h": 2.0},
]
```



## Documentation and Public Source Map

The public branch keeps the executable source and tests close to the contracts
they implement:

- **[Public API](hydrosheaf/api.py)**: Pipeline, temporal-edge, network-prior, and history-aware TTD entry points.
- **[Field data](hydrosheaf/data/field.py)** and **[input validation](hydrosheaf/data/validation.py)**: Dataset manifests, optional-signal resolution, and fail-closed module gating.
- **[Chemical units and registry](hydrosheaf/data/units.py)**: Explicit species, unit-conversion, and trace-speciation contracts.
- **[History-aware TTD](hydrosheaf/nuclear/history_ttd.py)**: Conditional source-history inversion and abstention reason codes.
- **[Active certified design](hydrosheaf/acmd.py)** and **[polyhedral solver](hydrosheaf/acmd_polytope.py)**: Sequential decisions, bounded ambiguity, and auditable certificates.
- **[Hydrochemical design](hydrosheaf/reactive_transport/chem_acmd.py)**: Conditional reaction and mixing constraints with declared measurement responses.
- **[Assumption diagnostics](hydrosheaf/validation/assumption_diagnostics.py)**: Evidence-limited checks for selected model assumptions.
- **[Topology-v2](hydrosheaf/validation/topology_v2.py)**: Truth-blind candidate generation, calibration, thresholds, and metrics.
- **[M4 fair inputs](hydrosheaf/validation/m4_fair.py)**: Projected Savage-frame and archive-informed observation construction.
- **[M4 path-aware evidence](hydrosheaf/validation/m4_path_aware.py)**: Truth-blind head-gradient/CBC path features.
- **[PEST++ adapter](hydrosheaf/calibration/pestpp/runner.py)**: External executable resolution and versioned calibration-run support.
- **[Reproducibility core](hydrosheaf/reproducibility/core.py)**: Deterministic environment, hashes, manifests, and tree comparison.
- **[M4 benchmark runner](M4/m4_topology_benchmark/scripts/run_m4_fair_topology_benchmark.py)**: Isolated fair-benchmark orchestration.
- **[Development workflow](scripts/reproduce_outputs.py)**: Isolated recipe reruns and strict output comparison.

Manuscripts, private research reports, raw/derived field datasets, generated
figures, and historical output trees are intentionally not linked here because
they are not part of the public source release.

## Troubleshooting & Common Questions

### Q: Do I need to compile anything to use Hydrosheaf?
**A:** No compilation is required for the core Python package. Run `pip install .`.
Compilation or installation of an external simulator is only relevant to a
specific PEST++/MODFLOW/MT3DMS adapter workflow.

### Q: What if I do not have PEST++ installed?
**A:** Core Hydrosheaf inverse modeling, history-aware TTD, topology-v2, and
reproducibility workflows do not require PEST++. Install a compatible PEST++
release separately, or allow the PEST++ adapter to attempt its configured
download, only when using an adapter that explicitly needs it. MODFLOW and
MT3DMS binaries are not auto-managed.

### Q: Will it work on Linux/Mac?
**A:** The Python package is cross-platform:

- Core Hydrosheaf source works on Windows/Linux/macOS.
- External simulator binaries must be installed separately for the relevant
  platform.
- Most users do not need external simulator binaries.

### Q: Can I use Hydrosheaf without PHREEQC?
**A: Yes!** While the PHREEQC library is installed by default, you can disable thermodynamic constraints in your configuration if you don't need them.

### Q: I got an error about missing dependencies
**A:** Ensure you have installed the package using `pip install .`. For PyVista/VTK-based 3D plotting or VTK export, install the optional visualization extra with `pip install ".[viz3d]"`.


---

## Reference Data & External Citations

This project utilizes reference data and documentation from the following official sources:

*   **USGS Groundwater Age Distribution Data:**
    *   **Title:** Data for distribution of groundwater age in aquifers used for public supply, United States
    *   **Source of data:** Reference data
    *   **Repository name:** Other (USGS Data Release)
    *   **DOI:** [https://doi.org/10.5066/P9W7T0DN](https://doi.org/10.5066/P9W7T0DN)
*   **USGS Savage Well Site Model Archive (MODFLOW/MODPATH):**
    *   **Title:** MODFLOW-2005, MODPATH, and MOC3D model archive for the Savage Municipal Water-Supply Well site
    *   **Source of data:** Reference data
    *   **Repository name:** Other (USGS Data Release)
    *   **DOI:** [https://doi.org/10.5066/F7J102FK](https://doi.org/10.5066/F7J102FK)
*   **PHREEQC Documentation & Examples:**
    *   **Title:** PHREEQC version 3 examples and documentation
    *   **Source of data:** Reference data
    *   **Repository name:** Other (USGS Techniques and Methods)
    *   **DOI:** [https://doi.org/10.3133/tm6A43](https://doi.org/10.3133/tm6A43)

---

## Authors

**Dickson Abdul-Wahab**,
**Ebenezer Aquisman Asare**,
**Abdul Rashid Dickson**

## License

[MIT License](LICENSE) (See LICENSE file for details)
