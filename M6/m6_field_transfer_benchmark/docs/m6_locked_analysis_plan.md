# M6 Locked Analysis Plan (Stage B1)

Locked 9 July 2026. This plan is the contract for the M6 field-transfer analysis. All
artifacts must trace to an experiment defined here.

## Reused, frozen upstream assets (no re-fitting)
- `m5_common.py` primitives: `fit_inverse` (FISTA sparse inverse), `reaction_matrix`,
  `matrix_diagnostics`, `equivalence_classes`, `support_from_extents`,
  `precision_recall_f1`, `jaccard`, `thermodynamic_bounds`, `MOLAR_MASS_G_MOL`,
  `ION_ORDER` (Ca,Mg,Na,K,HCO3,Cl,SO4,NO3,F,Fe,PO4), 16-reaction dictionary.
- **Frozen M5 MRS calibration** `M5/.../results/mrs_calibration_model.json` — the
  *transferred* identifiability classifier. M6 applies it unchanged; it is never re-fit.
  This is the scientific core of a *transfer* study.

## Inference unit and residual
- Concentrations converted mg/L → **mmol/L** (divide by molar mass).
- **Northern Ghana**: within-well seasonal residual `r = x_dry − x_wet` (mmol/L),
  one inference unit per well (160). Thermodynamic bounds from workbook SI
  (Calcite/Dolomite/Gypsum/Halite) — no live PHREEQC needed.
- **Edges** (Exp 4, and external sets): residual `r = x_target − x_source` along a
  candidate edge; unbounded non-negative fit where SI absent.

## Evidence tiers (reaction panel + evidence lifters)
Reaction-basis ions available per tier drive `matrix_diagnostics`; isotopes/Sr/SiO₂/SI/
metadata/edges are **evidence lifters/gates**, not linear-panel members.
- Tier 0: majors {Ca,Mg,Na,K,HCO3,Cl,SO4,NO3} (8-ion panel)
- Tier 1: +δ18O/δ2H (recharge/mixing evidence)
- Tier 2: +F (9-ion panel; enables fluorite/apatite)
- Tier 3: +Sr,SiO₂ (weathering endmember evidence)
- Tier 4: +aquifer/geology/lithology metadata, SI (thermo gate), QC class, graph edges

Native tier caps: Northern Ghana = Tier 4; Manu = Tier 2; Talensi = Tier 1.

## Identifiability classification (transferred)
Per inference unit: 9-feature vector = [rank_skill, coherence_skill, support_stability,
predictive_skill, thermodynamic_skill, reconstruction_rmse, heldout_rmse, n_measured_ions,
l1_norm]; standardise by frozen means/scales; prepend 1; `class = RESOLUTION_CLASSES[
argmax(design·class_coefficients)]`; if class==identifiable but any selected reaction sits
in an ambiguous equivalence class → downgrade to `equivalence_class`. MRS =
100·clip(design·reliability_coefficients,0,1). support_stability from bootstrap over
5% analytical noise (n=64); predictive_skill from leave-one-ion held-out RMSE.

## Experiments → outputs
- **E1 Readiness**: per-dataset tier attainment, missingness, CBE/Data_Class →
  `m6_dataset_readiness.csv`, `m6_variable_availability.csv`, `m6_cbe_distribution.csv`.
- **E2 NG seasonal×aquifer transfer**: 160 wells at Tier 4; reaction-class support,
  identifiability class, MRS by aquifer×season; wet→dry concordance →
  `m6_ng_field_pairs.csv`, `m6_ng_class_support.csv`, `m6_ng_aquifer_season_summary.csv`.
- **E3 Tier ablation**: re-run NG at Tiers 4→0 removing F/isotopes/Sr-SiO₂/SI/metadata/
  edges; identifiability-class migration, process-label flips →
  `m6_tier_ablation.csv`, `m6_tier_ablation_transitions.csv`.
- **E4 Edge uncertainty**: 4 edge sets (provided graph, chemistry-kNN, geographic-nearest,
  random/perturbed) + within-well; inferred process-network stability →
  `m6_edge_sensitivity.csv`, `m6_edge_network_summary.csv`.
- **E5 External transfer**: Talensi (Tier 1) + Manu (Tier 2) via chemistry-kNN edges;
  transfer success, uncertainty increase vs NG at matched tier, failure modes →
  `m6_external_transfer.csv`, `m6_external_summary.csv`.
- **E6 Limitation map**: synthesise per-process-class identifiable/weak/non-identifiable
  across datasets; conservative gated vs single-best-fit → `m6_limitation_map.csv`,
  `m6_conservative_vs_bestfit.csv`.

## Honesty boundary (enforced)
Mendeley `Dominant_Process`, `Aquifer_Evolution_Label`, `Graph_Edges`, and heuristic
residence times are **inferred references, never validation targets**. All "accuracy"
against them is reported as *concordance with prior labels*, not ground truth.

## Determinism
Fixed RNG seed (1234) for bootstrap/edge perturbation. No `Date.now`/wallclock in outputs.
