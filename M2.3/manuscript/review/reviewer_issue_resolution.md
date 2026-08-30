# Reviewer issue resolution ledger

Issues raised in `manuscript_reviewer_report.md`, and their disposition. Four
issues were classified Critical; three were resolved with data in hand and one
remains open.

## Critical

| ID | Issue | Disposition |
|---|---|---|
| C1 | The 100-realisation benchmark is 100 noise replicates of one fixed scenario, not 100 independent scenarios; scored rows overstate independent information. | **Resolved.** `derive_identifiability.py` quantifies the replication structure: truth does not vary across realisations for any target, and the independent units are 6 transport parameters, 14 site residence times and 21 active edge-by-reaction combinations. Methods 3.4 now states the design explicitly. Results 4.2 reports statistics at both row and clustered level; clustered active-term R² is 0.011 against a row-level −0.023, and clustered transport error 0.064 against 0.058. The qualitative conclusions hold at the correct level. |
| C2 | Regularisation weights governing the reaction inversion are unreported, so a structural failure cannot be distinguished from a badly tuned penalty. | **Superseded, and partially open.** The structural analysis under C3 removes the dependence of the central claim on hyperparameters, since a null space cannot be closed by any penalty. The weights themselves are still not recorded in the shipped benchmark configuration and are not reported; this is now stated as a limitation rather than silently omitted. A sensitivity sweep remains a useful addition but is no longer load-bearing. |
| C3 | The claim that point reaction extents are "not recoverable" generalised from one solver on one scenario. | **Resolved, and the claim is now stronger.** The stoichiometry matrix was assembled from the declared dictionary: 14 reaction terms, 8 binding ions, rank 8, null-space dimension 6, condition number 13.65. Six explicit null directions were computed with maximum residual ion change at machine precision, and are chemically interpretable (cation-exchange pairs; calcite against albite and dolomite; sulfate reduction against calcite and gypsum). Non-identifiability is now established as a property of the forward map. Reported in new Table 5 and rebuilt into Results 4.2, Discussion 5.1 and the Conclusions. |
| C4 | The 54.1% leakage figure had no baseline and the 0.05 mmol/L threshold was unjustified. | **Resolved.** Leakage is now reported against detection of genuinely active terms at six thresholds from 0.01 to 0.50 mmol/L. At every threshold the inversion asserts extents on truth-zero terms more often than it detects active ones (54.1% against 42.1% at 0.05 mmol/L). The threshold is justified against a propagated analytical-noise floor of approximately 0.035 mmol/L derived from the benchmark's 3.5% major-ion relative sigma. |

## Major

| ID | Issue | Disposition |
|---|---|---|
| M1 | Abstract over length and number-dominated. | Resolved. Reduced to 272 words and restructured to lead with the identifiability finding. |
| M2 | No stated hypotheses; HARKing exposure on the negative result. | **Open.** Requires an author statement on what was predicted before the benchmark was run. Cannot be resolved from the record. |
| M3 | Literature anchored pre-2002; thin post-2022. | **Open.** Requires a literature search and 3 to 5 additions. |
| M4 | Competence-matched baselines never specified. | Resolved. Both baselines are now described in Methods 3.1 with their observation contracts and what they were not permitted to see. |
| M5 | Six locked cases with no power justification. | Resolved by reframing. Discussion 5.6 now states that six cases are too few to discriminate a passing from a failing method and that the programme is better read as a working test of the claim architecture. |
| M6 | Talensi charge imbalance diagnosed without checking the source publication. | **Open.** Requires comparison of the project file against the published tables. |
| M7 | Field edge denominators switch between 208 and 161 without explanation. | Resolved. Results 4.5 now states that the 47 unlabelled edges are those fitted with an evaporation model, which references no local end-member, and reports their closure distribution. |
| M8 | External age comparison discards 47% of rows without comment. | Resolved. `public_age_attrition.csv` characterises the attrition: it is entirely reference-side, none lost on the estimate side, but uneven across aquifer groups from 9.1% glacial to 75.8% western unconsolidated, with unpaired sites deeper. Reported in Results 4.4. |
| M9 | Figure 2(a) is a coordinate scatter, not a map. | **Resolved.** Replaced with a projected map carrying Ghana's national and first-level administrative boundaries, neighbouring countries for context, labelled study regions, a graticule, scale bar, north arrow and a locator inset showing the study frame within Ghana. Boundaries come from geoBoundaries (CC BY 4.0), fetched once by `fetch_boundaries.R` and cached to a local GeoPackage so the figure renders offline. Replacing the scatter also exposed a coordinate error the scatter had hidden; see below. |
| M10 | 28 references is thin. | **Open.** Follows from M3. |
| M11 | Section 4.1 read as a capability inventory. | Resolved. Cut to a single paragraph; the inventory lives in Table 2. |
| M12 | The measurement-value result was under-sold. | Resolved. Promoted to Discussion 5.5 with the structural explanation of why strontium and silica add independent constraints while isotopes and fluoride do not. |

## Found while resolving M9

Replacing the coordinate scatter with a real map exposed a data error that the
scatter had hidden. All 63 Talensi longitudes were stored as positive, placing
every sample outside Ghana; Talensi District lies near 0.8 degrees west. The sign
is negated in the derivation script, the source file is left unmodified, and the
correction is recorded in `field_provenance.json`, Table 4 and Section 2.2.
Because reflection preserves pairwise distance exactly, no distance-based result
changed. A point-in-polygon gate (`check_coordinates.R`) now runs in the pipeline
so that this class of error cannot recur silently, and it also flags three
Northern Ghana samples sitting 1.4 to 2.0 km beyond a generalised boundary
polygon, which is within cartographic tolerance.

An earlier version of the FIG-2 validation record claimed that coordinates had
been checked against the Ghana boundary. That claim rested on visual inspection
and was wrong. It has been corrected.

## Minor

Mixed-class residence-time degradation now discussed in Results 4.2. Terminology
fixed on "in-sample chemistry closure". Remaining open: Figure 4(b)
overplotting, Table 4 framing, and the dense locked-programme paragraph.

## Submission blockers, unchanged

Public repository DOI and URL are placeholders; Computers & Geosciences requires
a live documented public repository at submission. Clean-tree regeneration of the
locked run is outstanding. Named CRediT roles are outstanding. Redistribution
terms for the two published field datasets are unconfirmed.
