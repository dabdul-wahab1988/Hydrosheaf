# Supplementary Information

**Hydrosheaf: A Reproducible Computational Framework for Inferring Directed Groundwater Hydrochemical Evolution Networks from Sparse Multitracer Data**

Dickson Abdul-Wahab, Ebenezer Aquisman Asare, Dickson Adomako, Gibrilla Abass, Samuel Ganyaglo

This file accompanies the revised manuscript CAGEO-D-26-00847 (Computers and Geosciences). All values below were regenerated from the locked pipeline (git commit `463e1ce`, Python 3.14 environment, `M2/m2_benchmark` and `M3/m3_age_benchmark` runners); see the manuscript's metric audit appendix for provenance.

## Supplementary Table S1. Input data schema

| Data type | Required field | Unit | Required/optional | QC rule | Used in module |
|---|---|---|---|---|---|
| Location | latitude/longitude or x/y | decimal degrees or projected | required | valid coordinate range | graph |
| Elevation/head | elevation, head_m | m | recommended | non-missing preferred | topology |
| Major ions | Ca, Mg, Na, K, HCO3, Cl, SO4, NO3 | mmol/L or mg/L | required | charge-balance screening | chemistry |
| Trace tracers | F, Fe, PO4, Br, NH4 | mg/L or ug/L | optional | detection-limit handling | logic gates |
| Isotopes | 3H, 14C, d18O, d2H | TU, pmC, permil | optional/enhanced | tracer range check | age |
| Field parameters | pH, EC, temperature, DO/Eh | field units | recommended | field range check | geochemistry |

## Supplementary Table S2. Reaction dictionary (benchmark configuration)

| Process family | Reaction/process term | Chemical signature | PHREEQC/SI rule | Allowed direction |
|---|---|---|---|---|
| carbonate dissolution | calcite | Ca, HCO3 | SI/forward check when available | nonnegative |
| carbonate dissolution | dolomite | Ca, HCO3, Mg | SI/forward check when available | nonnegative |
| evaporite dissolution | gypsum | Ca, SO4 | SI/forward check when available | nonnegative |
| evaporite dissolution | halite | Cl, Na | SI/forward check when available | nonnegative |
| silicate weathering | albite | HCO3, Na | SI/forward check when available | nonnegative |
| redox process | pyrite oxidation (pyrite_net) | Fe, SO4 | SI/forward check when available | nonnegative |
| nitrate input | NO3src | NO3 | SI/forward check when available | nonnegative |
| redox process | denitrification (denit) | HCO3, NO3 | SI/forward check when available | nonnegative |
| redox process | sulfate_reduction, iron_reduction | SO4, Fe | SI/forward check when available | nonnegative |
| ion exchange | CaNa_exch | Ca, Na | SI/forward check when available | signed |
| ion exchange | NaCa_exch | Ca, Na | SI/forward check when available | nonnegative |
| ion exchange | MgNa_exch | Mg, Na | SI/forward check when available | signed |
| ion exchange | NaMg_exch | Mg, Na | SI/forward check when available | nonnegative |

Note: the benchmark dictionary is 14 reactions × 11 ions (Ca, Mg, Na, K, HCO3, Cl, SO4, NO3, F, Fe, Br per the locked configuration). The dictionary is rank-deficient (rank 8 of 11, condition number infinite) — see Supplementary Tables S6-S8 and Section 4.2 of the main text.

## Supplementary Table S3. Locked parameter defaults (main text Section 2 and Table S3)

| Parameter | Value/range | Used in module | Scientific rationale |
|---|---|---|---|
| n_realisations | 100 | benchmark | Monte Carlo recovery ensemble |
| lambda_l1 (default) | 0.002 | sparse inversion | default sparse-recovery penalty |
| AICc-selected lambda | 0.0483 (both field sites) | site model selection | AICc minimum on the regularization path (generic dictionary) |
| lambda_l2 (Tikhonov) | 0.0 (ridge floor ~1e-10 added) | sparse inversion | fixed stabilisation, not optimised with lambda_1 |
| major_ion_rel_sigma | 0.04 | synthetic data | configurable default for routine major-ion analytical precision |
| isotope_sigma_permil | 0.50 | isotope validation | configurable default for stable-isotope measurement noise |
| transport_models_enabled | evap, mix | edge fitting | minimum transport alternatives; both fitted, better-fitting candidate selected by the combined weighted objective (transport + residual chemistry + L1/EC-TDS/isotope/kinetic penalties) |
| edge_max_neighbors | 3 | graph construction | bounds candidates (~3n), avoids long-range false edges, retains sparse-network coverage |
| edge_radius_km (search radius) | 5.0 km | graph construction | a priori reproducibility-fixed benchmark default, not a site-calibrated hydraulic length scale; graph-stage sensitivity at 3, 5, 7.5 and 10 km is reported in Table S14, and field topology is conditional on this choice |
| edge_p_min (edge confidence threshold) | 0.75 | graph construction | minimum accepted directional confidence p_ij = Phi((h_i - h_j) / sigma_dh) |
| edge_gradient_min (minimum hydraulic gradient) | 1e-4 | graph construction | directional consistency floor; pairs below it but inside edge_radius_km are kept separately as lateral dispersive mixing candidates with p_ij shrunk toward 0.5 |
| edge_depth_mismatch | 20 m | graph construction | screen-depth difference beyond which an edge is flagged as depth-mismatched |
| sigma_meas / sigma_elev / sigma_topo (head uncertainty by evidence tier) | 0.5 / 1.0 / 10.0 m | graph construction | per-node head sigma for measured head, elevation proxy and topographic fallback; sets sigma_dh = sqrt(sigma_i^2 + sigma_j^2) in p_ij when no head covariance matrix is supplied |
| PSI trials | 30 default | field uncertainty | edge process stability under input perturbation |
| SI threshold (diagnostic band) | 0.2 | PHREEQC/forward checks | thermodynamic feasibility screening |
| age-ordering min_downstream_increase (epsilon) | 0.0 y | age constraint | ordering evaluated on posterior intervals; interval overlap grades severity |
| age severe log10 threshold | 0.3 | age constraint | non-overlapping reversals beyond this are severe (failure codes) |
| tau_agreement_tolerance | 0.4 | tracer inversion | relative spread threshold for tau-ambiguous flags |

## Supplementary Table S4. Pilot metadata

| Site | Sample ID | Latitude | Longitude | Elevation (m) | Aquifer setting | Land-use/mining influence | Available tracers | Completeness score |
|---|---|---|---|---|---|---|---|---|
| Talensi | TGW1 | 10.71 | 0.82 | 217 | Crystalline basement | Mining/Agri | stable isotopes; no nuclear tracers | 1.00 |
| Talensi | TGW2 | 10.71 | 0.81 | 222 | Crystalline basement | Mining/Agri | stable isotopes; no nuclear tracers | 1.00 |
| Lower Anayari | PUB1 | 10.91 | -1.06 | 217 | Alluvial/Basement | Agriculture | stable isotopes; no nuclear tracers | 1.00 |
| Lower Anayari | PUB2 | 10.91 | -1.06 | 217 | Alluvial/Basement | Agriculture | stable isotopes; no nuclear tracers | 1.00 |

## Supplementary Table S5. Field edge outputs (top edges by rank score)

| Edge ID | From node | To node | Distance (km) | Elevation/head relation | Process stability (PSI) | Age consistency | Chemical match R2 | Dominant reaction | Status |
|---|---|---|---|---|---|---|---|---|---|
| Talensi_TGW49->TGW36 | TGW49 | TGW36 | 3.19 | downgradient | 1.00 | field tracer absent | 1.00 | redox process | screening-level demonstration |
| Talensi_TGW34->TGW19 | TGW34 | TGW19 | 4.62 | downgradient | 1.00 | field tracer absent | 1.00 | nitrate input | screening-level demonstration |
| Talensi_TGW42->TGW45 | TGW42 | TGW45 | 3.36 | downgradient | 1.00 | field tracer absent | 1.00 | nitrate input | screening-level demonstration |
| Talensi_TGW3->TGW8 | TGW3 | TGW8 | 4.93 | downgradient | 1.00 | field tracer absent | 1.00 | redox process | screening-level demonstration |
| Talensi_TGW32->TGW12 | TGW32 | TGW12 | 4.68 | downgradient | 1.00 | field tracer absent | 1.00 | nitrate input | screening-level demonstration |
| Talensi_TGW44->TGW45 | TGW44 | TGW45 | 2.70 | downgradient | 1.00 | field tracer absent | 1.00 | nitrate input | screening-level demonstration |

These are the six highest-ranked field edges by rank score (chemistry R2 x PSI), regenerated from the canonical run; they are the same six edges reported in main-text Table 6. Ranking is strictly ordered with no ties, so the selection is reproducible.

Null edges: none. The seven zero-closure edges reported at the previous revision (Manu_9_11, Manu_10_11, Manu_12_11, Manu_19_11, Talensi_2_52, Talensi_3_52, Talensi_4_52), which converged on the sink wells NAB14 and THDW1 with chemistry R2 = 0.0 and all extents ~0, were artefacts of the unfiltered candidate graph used previously. They are rejected by the directional refinement and no unresolved null edge remains. Fourteen of the 258 retained edges return a negative chemistry R2 and are reported as poor fits rather than removed.

Direction convention and column definitions. Neither field site has hydraulic-head records, so the directional confidence p_ij of Section 2.3 is evaluated on the elevation proxy using the elevation-tier sigma of 1.0 m (Table S3). Every retained edge therefore carries a real directional confidence: the median is 0.93 at Talensi and 0.50 at Lower Anayari. The 0.50 at Lower Anayari is not a missing value but the exact value of Phi(0) - station elevations there are recorded as constants, so most Lower Anayari pairs fall below the gradient floor and are retained as lateral, dispersive mixing candidates rather than as gradient-supported directed edges. No retained edge at either site runs against the elevation proxy. The column headed "Process stability (PSI)" is the Monte Carlo process-stability index of Equation 7, identical to the "PSI probability" column of main-text Table 6; it is not the directional confidence and asserts nothing about hydraulic plausibility.

The column previously headed "Edge confidence" is renamed here to "Process stability (PSI)", because that is what the pipeline writes into it: the Monte Carlo process-stability index of Equation 7, identical to the "PSI probability" column of main-text Table 6. It is **not** a hydraulic confidence and is **not** derived from head data — no edge-confidence quantity is carried in the field results at all, because the directional confidence p_ij of Section 2.3 requires head or elevation uncertainties that these sites do not supply. A value of 1.00 therefore means that the dominant reaction was selected in every one of the 30 Monte Carlo perturbation trials for that edge; it carries no claim about the hydraulic plausibility of the connection. We thank the reviewer for pressing on this point: the original heading invited exactly the misreading that a confidence of 1.00 had been asserted for edges with no head data.

## Supplementary Table S6. Identifiability diagnostics: dictionary structure

| Metric | Value |
|---|---|
| n_reactions | 14 |
| n_ions | 11 |
| rank | 8 |
| rank deficiency | 3 |
| condition number | infinite |
| column pairs with |cosine| > 0.7 | 10 (calcite~anorthite 1.00; CaNa_exch~NaCa_exch 1.00; MgNa_exch~NaMg_exch 1.00; calcite~dolomite 0.95; pyrite_net~sulfate_reduction 1.00; etc.) |

## Supplementary Table S7. Extent- and family-level recovery (canonical, 100 realisations)

| Metric | Value |
|---|---|
| active reaction rows | 2,100 |
| active correlation R2 | 0.23 |
| active MAE | 0.37 mmol/L |
| active RMSE | 0.62 mmol/L |
| false-activation rate (>0.05 mmol/L) | 54.1% |
| active family rows | 1,700 |
| family correlation R2 | 0.16 |
| family MAE | 0.47 mmol/L |
| family sign-match rate | 36.5% |
| dominant-family hit rate | 48% |

## Supplementary Table S8. Dictionary sensitivity (leave-one-out, 5 locked realisations)

| Dictionary variant | n active rows | active corr-R2 | active MAE (mmol/L) |
|---|---|---|---|
| full dictionary | 441 | 0.09 | 0.40 |
| minus calcite | 291 | 0.02 | 0.32 |
| minus dolomite | 341 | 0.14 | 0.43 |
| minus gypsum | 400 | 0.08 | 0.42 |
| minus halite | 391 | 0.10 | 0.41 |
| minus albite | 441 | 0.17 | 0.34 |
| minus pyrite_oxidation_aerobic | 441 | 0.05 | 0.40 |

No single removal resolves the degeneracy; removing albite gives the largest improvement.

## Supplementary Table S8b. Reaction correlation and parameter-level uncertainty (recovered extents, 100 realisations)

Empirical correlations of recovered reaction extents across the synthetic benchmark (pairs with |r| > 0.5):

| Reaction A | Reaction B | Pearson r (recovered extents) |
|---|---|---|
| albite | pyrite_oxidation_aerobic | −0.715 |
| CaNa_exch | pyrite_oxidation_aerobic | 0.692 |
| calcite | dolomite | 0.661 |
| CaNa_exch | albite | −0.575 |
| gypsum | pyrite_oxidation_aerobic | −0.508 |

Parameter-level precision of recovered extents over active reaction rows (mean and standard deviation across 100 realisations; the SD is the empirical precision of the penalised solver, not a formal posterior interval — a fully Bayesian alternative is available, Supplementary Method S5):

| Reaction | n active rows | mean recovered extent (mmol/L) | SD (mmol/L) |
|---|---|---|---|
| calcite | 600 | +0.055 | 0.137 |
| dolomite | 400 | +0.083 | 0.109 |
| CaNa_exch | 400 | −0.057 | 0.168 |
| gypsum | 400 | −0.032 | 0.188 |
| halite | 300 | −0.019 | 0.186 |

The strong calcite–dolomite and albite–exchange correlations are the empirical signature of the dictionary degeneracy quantified in Tables S6–S8 and explain the limited extent-level recovery reported in Section 4.2.

## Supplementary Table S9. Topology baseline comparison (Savage MODPATH reference, 174 directed edges, DOI 10.5066/F7J102FK)

| Baseline construction | Candidates | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| head-gradient (elevation-as-head, downhill, k=2) — framework no-prior | 302 | 147 | 155 | 27 | 0.49 | 0.84 | 0.62 |
| all-pairs elevation-drop | 452 | 147 | 305 | 27 | 0.33 | 0.84 | 0.47 |
| proximity kNN (k=2) | 306 | 0 | 306 | 174 | 0.00 | 0.00 | 0.00 |
| conservative-tracer ordering | — | — | — | — | — | — | not evaluable on this archive (no hydrochemistry) |

## Supplementary Table S10. PSI separation for degenerate reaction pairs (site-aggregated, canonical)

| Site | Pair | PSI (member A) | PSI (member B) | Separation |
|---|---|---|---|---|
| Lower Anayari | CaNa_exch ~ NaCa_exch | 0.73 | 0.46 | CaNa > NaCa by 0.27 |
| Talensi | CaNa_exch ~ NaCa_exch | 0.83 | 0.42 | CaNa > NaCa by 0.41 |
| Lower Anayari | calcite ~ dolomite | 0.05 | 0.002 | both near zero |
| Talensi | calcite ~ dolomite | 0.06 | 0.04 | both near zero |

Values are the site-aggregated per-mineral inclusion probabilities (mean Monte Carlo inclusion probability across that site's edges), identical to the values plotted in Figure 7. PSI separates the directional exchange pair (CaNa_exch consistently higher than NaCa_exch at both sites), whereas both carbonate members remain near zero, so the metric does not falsely promote either carbonate when carbonates are not stably active. PSI measures stability to input perturbation, not the correctness of process attribution.

## Supplementary Table S11. Multi-endmember mixing + reaction stress test (canonical)

| Metric | Value |
|---|---|
| edges x realisations | 10 x 20 |
| true mixing fraction (two endmembers) | 0.35 (f1 = 0.20, f2 = 0.15) |
| true halite extent | 0.40 mmol/L |
| true calcite extent | 0.20 mmol/L |
| mix model selected | 73.1% of edges |
| median recovered mixing fraction | 0.24 |
| median chemistry R2 | 0.999 |
| median recovered halite extent | 0.04 mmol/L |
| median recovered calcite extent | 0.00 mmol/L |
| median family evaporite extent | 0.23 mmol/L |

Single-endmember transport absorbs most of the reaction signal when the true mixing is two-endmember: mixing–reaction coupling is quantified and discussed in Section 5.5.

## Supplementary Table S12. Age-ordering with wide posterior intervals (canonical, synthetic benchmark network)

The age-ordering audit is run on the 13 directed edges of the synthetic benchmark network in each of 100 realisations (1,300 edge-checks in total); the fractions below are means over the 100 realisations.

| Metric | Value |
|---|---|
| directed edges per realisation | 13 (× 100 realisations = 1,300 edge-checks) |
| violation fraction (point-estimate ordering) | 15.2% |
| severe violations (non-overlapping reversal beyond log10 0.3) | 2.4% |
| share of violations with overlapping intervals (retained, flagged "unresolved at stated uncertainty") | 84.3% |
| age-order consistency index | 0.85 |

## Supplementary Table S13. Runtime scaling of candidate-edge construction

| n_nodes | n_candidate_edges | wall_seconds | edges_per_node |
|---|---|---|---|
| 10 | 24 | 0.0001 | 2.40 |
| 20 | 54 | 0.0003 | 2.70 |
| 40 | 114 | 0.0010 | 2.85 |
| 80 | 234 | 0.0035 | 2.92 |
| 160 | 474 | 0.0145 | 2.96 |
| 320 | 954 | 0.0619 | 2.98 |

Candidate counts grow ~linearly (~3n) thanks to the maximum-neighbour pruning (default 3). Measured end-to-end wall times for the canonical runs on a 16-logical-core workstation: the full 100-realisation M2 synthetic benchmark (4 scenarios x 4 topology variants, including generation and evaluation) completed the realisation loop in ~1.7 minutes (start-to-finish including setup under 5 minutes), and the 1,272-row M3 public-age benchmark completed in ~7 minutes (16 processes).

## Supplementary Table S14. Graph-stage sensitivity to the edge search radius

The documented probabilistic builder and sheaf-inspired selection stage were rerun at four search radii. Chemistry fitting was not rerun, so these are topology-sensitivity results, not new field-performance estimates.

| Site | Radius (km) | Candidate edges | Retained edges | Primary | Lateral | Jaccard vs 5 km |
|---|---:|---:|---:|---:|---:|---:|
| Lower Anayari (Manu) | 3.0 | 265 | 117 | 12 | 105 | 0.442 |
| Lower Anayari (Manu) | 5.0 | 422 | 121 | 21 | 100 | 1.000 |
| Lower Anayari (Manu) | 7.5 | 480 | 123 | 45 | 78 | 0.671 |
| Lower Anayari (Manu) | 10.0 | 485 | 123 | 50 | 73 | 0.616 |
| Talensi | 3.0 | 89 | 85 | 78 | 7 | 0.291 |
| Talensi | 5.0 | 150 | 137 | 131 | 6 | 1.000 |
| Talensi | 7.5 | 177 | 151 | 143 | 8 | 0.258 |
| Talensi | 10.0 | 192 | 162 | 152 | 10 | 0.182 |

The 5-km setting is therefore a reproducibility-fixed benchmark default. The changes in retained topology show why field metrics should be interpreted conditionally and recalculated when a site-specific hydrogeological length scale is justified.

## Supplementary Method S1. Multi-criteria edge selection and sheaf-cohomology diagnostics (derivation)

Let the groundwater sampling network be the directed graph G = (V, E) with nodes V and candidate directed edges E. Each node v carries a stalk s_v in R^d, the local observation vector of d measured quantities (major-ion concentrations in mmol/L, stable-isotope values in permil, and hydraulic head or elevation in m). Each directed edge e = (u -> v) carries an affine restriction map

  f_e(s_u) = alpha_e * s_u + delta_e,

where alpha_e is an edge scale parameter and delta_e the expected offset along each dimension (e.g., the sign-consistent change in head or age). The maps are prescribed from the transport and age evidence for that edge: alpha_e absorbs bulk concentration scaling (evaporation factor or mixing dilution on the anchor ions) and delta_e encodes expected directional changes (head drop, age increase). The edge consistency residual is

  rho(e) = sqrt( sum_d w_d ( alpha_e s_u,d + delta_e,d - s_v,d )^2 ),

with weights w_d on observation dimensions (higher weight on conservative tracers). A candidate edge is consistent when rho(e) is small relative to the observation uncertainty.

Global consistency aggregates residuals over retained edges:

  rho_G = (1/|E|) sum_{e in E} rho(e)^2.

**Edge selection.** Each candidate edge receives a weighted multi-criteria score from the four evidence streams, and the retained set is chosen per node by a soft (inverse-temperature) selection, `sheaf_soft_beta` default 2.0, capped at `edge_max_neighbors`. This step is a weighted multi-criteria score and is implemented as such (`hydrosheaf/sheaf/topology_refine.py`); it involves no coboundary and no Laplacian, and we do not claim it outperforms an equivalent hand-designed score.

**Coboundary and cohomology.** The sheaf structure is a diagnostic layer computed on the retained graph (`hydrosheaf/sheaf/cohomology.py`). Stacking the edge relations gives the affine coboundary operator D with one row per (edge, observation dimension):

  D[(e,d), (u,d)] = sqrt(w_e) * alpha_e,   D[(e,d), (v,d)] = -sqrt(w_e),   b[(e,d)] = -sqrt(w_e) * delta_{e,d}

so that D has shape (|E|·d) × (|V|·d) and the affine edge relations are represented by Dx = b. This is a cellular-sheaf coboundary matrix with d-dimensional node stalks and an upstream restriction map alpha_e I. The reported quantities are the homogeneous nullity H0 = dim ker D, the first cohomology dimension H1 = |E|·d − rank D, the affine obstruction energy min_x ||Dx − b||², the per-edge leverage (the increase in obstruction energy attributable to each edge), and the maximum cycle obstruction. H0 is therefore the dimension of the homogeneous nullspace, not by itself the dimension of an affine global-section set. An exact affine global section exists only when the obstruction energy is zero within numerical tolerance; if the system is consistent, its affine solution set has dimension H0, whereas a positive obstruction energy means that the affine solution set is empty. In particular, H0 = 0 means that the homogeneous system has only the zero solution, not that no affine node assignment exists. We do not compute a spectral decomposition of any Laplacian: localisation is performed by the leave-one-edge-out leverage statistic, not by eigenvectors.

**What this adds.** Edge selection is a weighted multi-criteria score and adds nothing over an equivalent well-designed score. The coboundary layer adds a network-level closure diagnostic: it tests whether the affine right-hand side b is compatible with the retained relations, measures the residual obstruction when it is not, and localises that obstruction with leave-one-edge-out leverage. H0 describes the homogeneous nullspace, not a yes/no test for affine solvability. Section 4.6 reports positive obstruction energies for both field networks, so neither has an exact affine global section at the reported numerical tolerance.

## Supplementary Method S2. Transport decomposition details

For each retained edge, the upstream concentration vector x_u is transformed by one of two operators. Under the evaporation model, x_T = gamma x_u with gamma >= 1. Under the mixing model, x_T = (1 - f) x_u + f x_e with 0 <= f <= 1 and x_e a configured endmember vector. Both gamma and f are fitted by unconstrained weighted least squares, with the transport fit weighted towards the conservative anchor signals (EC/TDS, Cl, Br); the better-fitting candidate is selected per edge as the one minimising the combined weighted objective — the transport and residual-chemistry mismatch plus the L1, EC/TDS, isotope and kinetic penalties (identical to the rule stated in Section 2.6). The transport-corrected residual r = x_v - x_T is the chemistry attributed to reactive processes alone. Caveats: Cl and EC are treated as conservative anchors only where halite dissolution, agricultural chloride input, or anthropogenic salinity is not indicated; the multi-endmember stress test (Table S11) quantifies the coupling limit when the true mixing involves two endmembers.

## Supplementary Method S3. 3H/3He apparent ages

Where 3H/3He data are available, tritiogenic helium ingrowth is used to compute apparent ages without reliance on historical input reconstruction: the tritiogenic 3He concentration is estimated from measured 3He/4He and neon excess-air corrections, and the apparent age follows from the 3H decay equation with the tritium half-life of 12.32 yr (Lucas & Unterweger, 2000).

## Supplementary Method S4. Network Bayesian age update

Node posterior ages are obtained from single-node LPM inversion expressed as probability distributions. The directed graph imposes ordering constraints tau_j >= tau_i - epsilon for each retained edge (u_i -> v_j); the network Bayesian update conditions each node posterior on all incident constraints via a sequential factor-graph update, shrinking posterior width where constraints are informative. Constraint evaluation is interval-aware: violations whose intervals overlap are flagged "not resolved at stated uncertainty" and excluded from the severe set; non-overlapping reversals beyond the log10 severity threshold (0.3) receive a severe age-coherence flag. The current audit does not apply a continuous overlap-proportional weight, and any downstream exclusion or penalty based on these flags is not applied to the field results reported here.

## Supplementary Method S5. Bayesian alternative

A fully Bayesian edge-fitting alternative is supported (Markov-chain Monte Carlo over reaction extents with the same dictionary and gates), providing posterior extent distributions and R-hat/ESS diagnostics; the results reported in the main text use the penalised weighted least-squares solver with Monte Carlo input perturbation for the PSI analysis.

## Provenance manifest: schema and example

The pipeline writes a JSON provenance manifest for each run, recording the run timestamp, configuration path, output paths, and the effective parameters. Schema (fields written by the benchmark runners):

```json
{
  "run_utc": "ISO-8601 UTC timestamp",
  "config": "path to locked configuration file",
  "output": "path to primary output CSV",
  "summary": "path to summary CSV (when produced)",
  "max_rows": "row cap or null",
  "age_steps": "age discretisation steps (when applicable)",
  "scenario_ids": ["scenario identifiers executed"],
  "include_withdrawn": false,
  "n_output_rows": 1272,
  "pairwise_baseline": "scenario id used as pairwise baseline"
}
```

Example from the canonical public-age run (reproduced in this revision, 2026-08-19):

```json
{
  "run_utc": "2026-08-19T14:42:45.290054+00:00",
  "config": "M3/m3_age_benchmark/configs/design_matrix.yaml",
  "output": "M3/m3_age_benchmark/results/m3_tracerlpm_parity_agefractions_full.csv",
  "summary": "M3/m3_age_benchmark/results/m3_tracerlpm_parity_agefractions_full_summary.csv",
  "max_rows": null,
  "age_steps": 90,
  "scenario_ids": ["tracerlpm_parity_agefractions"],
  "include_withdrawn": false,
  "n_output_rows": 1272
}
```

## Supplementary Figure S1. Public USGS residence-time parity check (identifiability-gated)

Among 1,272 rows, 356 fits passed the identifiability gate (916 flagged non-identifiable); among identifiable fits, median absolute log10 error = 0.022, log10 RMSE = 0.235, 91.9% within a factor of two and 98.6% within a factor of ten (log10 R2 = 0.960). These are the same values quoted in main-text Section 4.4, Table 2, Table 4 and the Figure 5b annotation. The USGS ages and age fractions are reported model outputs, so this is screening-level parity rather than direct true-age validation.

## Supplementary Figure S2. PHREEQC-compatible forward-feasibility proxy

Panel (a) shows the RMSE distribution (median 0.02 mmol/L); panel (b) relates RMSE to Nash-Sutcliffe efficiency (median NSE = 1.00); panel (c) shows the thermodynamic-feasible fraction (0.88). The simulator uses locked PHREEQC-consistent saturation fields; live PHREEQC kinetic execution remains pending.

## Supplementary Figure S3. Field-scale chemical discovery fit

Edge-level chemistry R2 for the retained field graphs: median R2 = 0.53 for Lower Anayari (121 edges) and 0.82 for Talensi (137 edges), overall median 0.70 over 258 retained edges. These fits describe hydrochemical closure of a field demonstration and do not provide independent process-truth validation.
