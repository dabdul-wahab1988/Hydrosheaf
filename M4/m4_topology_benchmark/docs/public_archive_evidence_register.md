# Public Archive Evidence Register

This register documents the validation tiers, evidence levels, and scientific boundaries for each public MODPATH archive ingested by the Hydrosheaf M4 pipeline.

---

## Tier 1: USGS Savage Municipal Well (Milford, NH)

- **Source Reference**: DOI: [10.5066/F7J102FK](https://doi.org/10.5066/F7J102FK)
- **MODPATH Version**: 5 (MODFLOW-2005)
- **Evidence Level**: **Validated (Topology and Source-Receptor Projection)**
- **Scale**: 3,000 particles, 410,769 pathline points, 174 reference edges.
- **Computed Metrics**:
  - Edge Precision, Recall, F1
  - Direction Agreement Rate
  - Capture-Envelope Convex Hull IoU and Source Cell Jaccard Index
  - Bootstrap Confidence Intervals (2,000 resamples)
  - Travel-Time Rank Correlation (Spearman $\rho$, Kendall $\tau$)
- **Allowed Claims & Guardrails**:
  - *Allowed*: Validates that Hydrosheaf reconstructs graph topology and receptor projections consistent with MODPATH pathlines.
  - *Guardrail*: Does not validate field geochemistry, capture-zone polygons, or travel times (since time definitions are not harmonised).

---

## Tier 2: Great Miami River Basin (Ohio)

- **Source Reference**: DOI: [10.5066/P9X4C9R6](https://doi.org/10.5066/P9X4C9R6)
- **MODPATH Version**: 7 (MODFLOW 6)
- **Evidence Level**: **Confirmation (MODFLOW 6 / MODPATH 7 compatibility)**
- **Scale**: 154 particles, 1,494 pathline points, 68 reference edges.
- **Computed Metrics**:
  - Edge Precision, Recall, F1
  - Direction Agreement Rate
  - Travel-Time Rank Correlation (Spearman $\rho = 0.999$, Kendall $\tau = 0.997$)
- **Allowed Claims & Guardrails**:
  - *Allowed*: Confirms compatible ingestion of MODFLOW 6 / MODPATH 7 endpoint and pathline datasets.
  - *Guardrail*: Capture-zone envelopes are omitted because release points do not form closed receptor point clouds.

---

## Tier 3: Long Island Nitrogen (New York)

- **Source Reference**: DOI: [10.5066/P97VFXZ4](https://doi.org/10.5066/P97VFXZ4)
- **MODPATH Version**: 7 (MODFLOW 6)
- **Evidence Level**: **Scalability Supplementary (Stub)**
- **Scale**: 0 particles, 0 pathline points (Stub mode).
- **Computed Metrics**: None (Wrote empty stubs).
- **Allowed Claims & Guardrails**:
  - *Allowed*: Used as a placeholder for pipeline scalability testing.
  - *Guardrail*: Contains no validation metrics; no scientific claims allowed.
