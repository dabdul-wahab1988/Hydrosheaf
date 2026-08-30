**Problem.** Groundwater inverse calibration can return a small residual without recovering the parameters that generated the observations, yet identifiability is rarely reported as a first-class calibration output.

**Objectives.** We quantify parameter recovery, uncertainty calibration and optimal measurement design in a PEST-style Gauss-Levenberg-Marquardt engine using production HydroSheaf transport and PHREEQC kinetic adapters.

**Method.** Three prospectively locked controlled-synthetic protocols used a noise-whitened log10-parameter Fisher information matrix, 250 paired replicates per strategy, percentile bootstrap intervals, an independently implemented grid-converged numerical truth, and 24 untouched MODFLOW 6/MODPATH 7 aquifer cases.

**Results.** All parameter-specific and joint D-, A-, E- and balanced criteria selected the same 50 d front observation; the anticipated target split was not confirmed. In the matched model, median dispersivity absolute log10 error fell from 0.2154 to 0.0173 and decay error from 0.0182 to 0.0146. Under independent numerical truth the same observation improved dispersivity (0.8262 to 0.1674) but degraded decay (0.1367 to 0.1541). In the PHREEQC kinetic adapter, one and four residence-time designs were both rank-one because the rate law identifies the product of rate constant and reactive surface area; an independent surface-area observation restored rank two. Linearised 95% coverage reached only 0.75-0.91 across the matched-model designs.

**Significance.** Measurement placement is parameter- and model-class-specific; structural product confounding requires independent parameter information. Identifiability diagnostics, not residuals alone, should accompany calibration reporting.
