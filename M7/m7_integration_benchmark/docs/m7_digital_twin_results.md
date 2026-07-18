# M7 Operational Synthetic-Twin Extension — Results

Generated deterministically from 24 independent process/observation
realisations with 80 ensemble members. Forecasts are prequential:
the forecast issued at month *t* cannot use observations after month *t*.

## Locked headline diagnostics

- Three-month head RMSE skill versus the open-loop twin:
  **11.9%** (replicate percentile interval
  -11.3% to 22.4%);
  the updated twin was better in **87.5%**
  of independent realisations.
- Three-month head RMSE skill versus oracle persistence:
  **39.4%**. Oracle persistence is intentionally
  conservative because it knows the exact current synthetic state.
- Empirical coverage of the nominal 90% three-month head interval:
  **84.7%**. Raw ensemble coverage was
  **34.6%**; the mean rolling
  monitoring-residual spread multiplier available at issue time was
  **3.78**. The correction is reported because
  raw EnKF uncertainty was substantially under-dispersed.
- Three-month head RMSE skill versus the updated wrong-topology control:
  **0.9%**.

## Interpretation

This extension demonstrates an operational mechanism: sparse observations update a
graph-constrained state ensemble, and that updated ensemble can be evaluated on
future states that were unavailable at issue time. It does **not** establish a field
digital twin, aquifer-specific fidelity, management safety, or transfer to Ghana.

The hidden truth and operational model intentionally differ in topology, coefficients,
nonlinearity, stochastic forcing, and an unmodelled chemistry pulse. This reduces the
inverse crime but does not reproduce all structural errors present in a real aquifer.
See `m7_digital_twin_protocol.md` for the claim boundary and required field extension.
