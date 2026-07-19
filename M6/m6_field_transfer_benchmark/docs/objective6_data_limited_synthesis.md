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
| Seasonal chemistry prediction | In 138 complete quantitative well pairs and 20 sequential August issue batches, graph-ridge MAE was lower than persistence by 0.168 log1p units, 95% CI [0.152, 0.185]. Its advantage over an expanding mean-delta model was 0.00584, CI [−0.00120, 0.01296]. | Within-campaign wet-to-dry chemistry is predictable. Incremental value of the processed graph over a simple seasonal baseline is not established. |
| Reaction-family plausibility | At the maximum M6 tier, all 160 Northern Ghana wells were classified as partially identifiable; conservative reporting admitted 7.3 reaction alternatives on average rather than 5.8 single best-fit reactions. | Family/equivalence-class interpretations are supportable; unique mineral reactions are not. |
| Measurement-value audit | Removing Sr/SiO2 changed the non-identifiable fraction from 0.6% at Tier 3 to 51.9% at Tier 2, while mean MRS remained approximately 71. | Sr/SiO2 are load-bearing corroborative measurements; fit/reliability scores alone can conceal structural non-identifiability. |
| Edge-set sensitivity | Network composition changed relative to the provided graph: total-variation distance was 0.093 for geographic-nearest, 0.171 for random-perturbed and 0.229 for chemistry-kNN edges. | Network-level process attribution depends on the assumed edge set. Exact field connectivity is not identified. |
| External sparse-data transfer | Talensi was 36% non-identifiable and charge-balance limited; Lower Anayari was 96.5% non-identifiable, principally for silicate-dominant edges without Sr/SiO2. | The framework can diagnose when sparse datasets are screening-only and identify the next useful measurement; it does not validate mechanisms there. |
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

The processed graph and prior `Dominant_Process` or
`Aquifer_Evolution_Label` fields are analytical references. Agreement with them
is concordance, not accuracy against independent truth.

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
