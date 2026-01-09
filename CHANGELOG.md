# Changelog

## 0.3.0 (2026-01-09)

### New Features

- **Dual Isotope Nitrate Apportionment**: Integrated a Bayesian mixing model using $\delta^{15}\text{N}$ and $\delta^{18}\text{O}_{\text{NO}_3}$ to rigorously distinguish manure/sewage from fertilizer sources. Priorities isotopic evidence over hydrochemical proxies when available.
- **Comparison with Commercial Software**: Added detailed technical comparison with PHREEQC, NETPATH, and GWB, highlighting Hydrosheaf's advantages in sparse optimization and automated network inference.

### Enhancements

- **Endmember Database**: Added `nitrate_endmembers.json` with literature-validated isotopic signatures (Kendall 1998).
- **Hybrid Inference Logic**: Robust fallback mechanism ensures seamless operation whether isotope data is present or absent.

## 0.2.0 (2025-01-09)

### Major Extensions

- **Reactive Transport**: Added kinetic validation of inverse results using Arrhenius-corrected rates and Damköhler number analysis.
- **3D Flow Networks**: Implemented 3D graph inference with vertical anisotropy and a Bayesian topographic prior for flow direction without head data.
- **Temporal Dynamics**: Added time-series support, residence time estimation via cross-correlation center-of-mass, and seasonal decomposition.
- **Uncertainty Quantification**: Integrated Bayesian MCMC (NUTS) for reaction extent posteriors and bias-corrected bootstrap (BCa) for confidence intervals.

### Enhancements

- **Nitrate Discrimination**: Added robust "Low Nitrate" gating and Compositional Data Analysis (CoDA) for manure vs. fertilizer distinction.
- **Thermodynamic Constraints**: Expanded PHREEQC integration to enforce saturation index (SI) bounds on reaction fitting.
- **Documentation**: Substantially updated `Technical Reference` with mathematically rigorous, unit-tested examples.

### Verification

- Added comprehensive test suite `tests/` merging previous scattered tests.
- Verified all numerical examples in documentation against `tests/test_doc_examples.py`.

## 0.1.0

- Initial package scaffold, core transport/reaction fitting, CLI, and tests.
- PHREEQC integration scaffolding with SI-based constraints and bounds.
- Endmembers JSON loader and constraint-aware outputs.
- Probabilistic edge inference with head/DTW/topography fallback and edge confidence metadata.
- Optional isotope penalties with LMWL configuration and diagnostics.
