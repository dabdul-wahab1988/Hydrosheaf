# E3 PHREEQC Forward Validation Report

Run timestamp: 2026-05-05T23:15:46Z

## Status
Completed live PHREEQC runs for the benchmark edges available in the 'complete' scenario.

## Summary Metrics

| Metric | Value |
| :--- | :--- |
| Edges Processed | 280 |
| Mean RMSE (mmol/L) | 1.413539 |
| Median RMSE (mmol/L) | 1.348509 |
| Mean NSE | 0.4713 |
| Median NSE | 0.5860 |

## Mineral Saturation Snapshot (Median SI)

| Phase | Median SI |
| :--- | :--- |
| Calcite | 2.32 |
| Dolomite | 4.25 |
| Gypsum | -6.44 |
| Halite | -6.86 |

## Interpretation
The live PHREEQC validation supports thermodynamic feasibility for Hydrosheaf's inverse reaction inferences.
The NSE and RMSE values indicate moderate forward predictive fit, not near-perfect agreement.
These diagnostics should be interpreted as a feasibility and residual check rather than a full reactive-transport calibration.
