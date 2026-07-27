# Objective 6 — data-limited Ghanaian application

## Revised objective

> Apply the framework and its component diagnostics to Ghanaian aquifer
> datasets to determine which integrated interpretations are supportable under
> available data and which remain non-identifiable.

The corresponding research question is:

> Which components and integrated interpretations of the proposed framework
> are supportable from available Ghanaian aquifer data, and which remain
> non-identifiable?

## Scope clarification

The Northern Ghana workbook does not contain environmental age tracers,
screen-top/screen-bottom intervals, repeated hydraulic heads, exact public
coordinates, independently verified flow paths or reaction-mechanism truth.
Stable water isotopes provide recharge/source evidence, not residence-time
truth. Elevation minus the single static-water-level field is at most a
single-occasion head proxy.

The term **Tier 4** in M6 means the highest tier of the M6
chemistry–isotope–metadata reaction-diagnostic ladder. It does not mean that all
HydroSheaf pillars are observed in the field.

## What the Ghanaian data support

| Component diagnostic | Available evidence and result | Supportable interpretation |
|---|---|---|
| Data readiness | Northern Ghana has 320 wet/dry samples from 160 wells. Independent charge-balance screening gives 294 quantitative, 19 screening and 7 exploratory records. | The chemistry is suitable for quantitative class-level diagnostics after explicit QC filtering. |
| Seasonal chemistry prediction | In 140 complete quantitative well pairs revealed in 20 fixed, disclosed arbitrary batches (no real sampling-date sequence exists in the canonical data), graph-ridge MAE was lower than persistence by 0.173 log1p units, 95% CI [0.158, 0.188]. Its advantage over an expanding mean-delta model was 0.00512, CI [−0.00174, 0.01191]. | Within-campaign wet-to-dry chemistry is predictable under a disclosed arbitrary revelation order. Incremental value of the self-generated connectivity graph over a simple seasonal baseline is not established. |
| Reaction-family plausibility | At the maximum M6 tier, all 160 Northern Ghana wells were classified as partially identifiable; conservative reporting admitted 7.1 reaction alternatives on average rather than 5.4 single best-fit reactions. | Family/equivalence-class interpretations are supportable; unique mineral reactions are not. |
| Measurement-value audit | Removing Sr/SiO2 changed the non-identifiable fraction from 0.6% at Tier 3 to 51.9% at Tier 2, while mean MRS remained approximately 71. | Sr/SiO2 are load-bearing corroborative measurements; fit/reliability scores alone can conceal structural non-identifiability. |
| Edge-set sensitivity | Comparing three Hydrosheaf-generated edge sets against each other (chemistry-kNN as reference; no imported graph is used), total-variation distance was 0.12 for geographic-nearest and 0.05 for random-perturbed edges. | Network-level process attribution depends on the assumed edge set. Exact field connectivity is not identified. |
| External sparse-data transfer | With saturation indices computed independently for both external datasets (not withheld as unavailable), Talensi was 37.2% non-identifiable and charge-balance limited, materially lower than its matched Northern Ghana reference (54.2%); Lower Anayari was 95.3% non-identifiable, principally for silicate-dominant edges without Sr/SiO2, a family saturation indices do not corroborate. | The framework can diagnose when sparse datasets are screening-only and identify the next useful measurement; it does not validate mechanisms there. |
| Stable-isotope use | d18O and d2H are present for the Northern Ghana campaign. | Recharge/source and mixing plausibility can be discussed; groundwater age cannot be inferred from these isotope fields. |

## What remains non-identifiable

Under the available Ghanaian data, the following cannot be validated:

1. groundwater residence time or an age distribution;
2. exact directed field flow paths;
3. screen-resolved vertical connectivity;
4. unique reaction mechanisms;
5. a coupled age–topology–reaction field truth; and
6. an operational Ghanaian digital twin with validated state updating and
   forecasting.

An earlier revision of this workflow additionally read a different, antecedent
study's own inferred `Dominant_Process` and `Aquifer_Evolution_Label` fields, plus an
imported graph-edge sheet, for the same boreholes. Those fields, and the workbook that
supplied them, have been removed from this study entirely rather than retained as a
concordance reference (`DECISIONS.md`).

## How M6 and M7 now fit together

- **M6** is the Ghanaian data-limited component-diagnostic application:
  readiness, reaction equivalence, tier ablation, edge sensitivity, external
  transfer and limitation mapping.
- **M7.2** provides the strict within-campaign seasonal chemistry hold-forward
  test and shows that graph ridge is not distinguishable from the stronger
  expanding mean-delta baseline.
- **M7.3** uses fresh independent synthetic truth to test when combining
  hydraulics, age and chemistry helps, is redundant or creates false
  confidence.
- The optional digital-twin experiment remains a synthetic methods extension;
  it is not a Ghanaian field digital twin.

## Objective-level conclusion

Objective 6 is on course under the revised wording. The data support a
conservative field application that identifies reaction-family plausibility,
measurement value, seasonal predictability and sensitivity to assumed
connectivity. They do not support a complete integrated field reconstruction.
That distinction is the result of the objective, not a deficiency to hide.
