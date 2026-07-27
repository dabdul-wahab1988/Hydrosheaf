# M6 manuscript-ready asset plan

## Recommended six-figure main narrative

1. **Dataset and evidence ceiling:** existing Figure 1, explicitly describing
   Northern Ghana coordinates as masked and Tier 4 as the maximum M6
   chemistry/metadata tier.
2. **Field-transfer workflow:** existing Figure 2.
3. **Northern Ghana reaction-class stability:** existing Figure 3.
4. **Diagnostic ablation:** existing Figure 4, the principal M6 result.
5. **Within-campaign seasonal hold-forward:** Figure 5. Graph ridge is
   superior to persistence but not distinguishable from the expanding
   mean-delta baseline.
6. **Limitation map:** existing Figure 6, showing conservative equivalence-class
   reporting and external-data non-identifiability.

Extended Data Figure 1 reports external transfer. Extended Data Figures 2 and
3 answer validation and circularity objections rather than advancing the main
field narrative.

## Main tables

1. Dataset readiness and claim strength.
2. Variable availability.
3. Northern Ghana evidence-tier ablation.
4. External sparse-data transfer.
5. Northern Ghana truth-free seasonal hold-forward.

The `Tier 4` label must always be expanded as “maximum M6
chemistry/metadata tier”. It is not a complete field age–head–screen evidence
tier.

## Reproduction

```powershell
.venv\Scripts\python.exe M6\m6_field_transfer_benchmark\scripts\run_m6_q1.py
```

The R Figure 1 uses `sf` and `ggspatial` with the unsimplified 2021
geoBoundaries ADM1 Shapefile (16 Ghana regions) and an aligned national outline
derived from those polygons. Geometry and masked well coordinates are
transformed to WGS84 / UTM zone 30N (EPSG:32630); the map includes a CRS-aware
scale bar and true-north arrow and makes the coordinate mask explicit. The R
suite is authoritative for Figures 1–4 and 6, all three Extended Data figures
and Supplementary Figures S1–S11. Python Figure 5 uses the same two-column
dimensions, accessible palette and vector/raster export conventions.
