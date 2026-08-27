# Response to Reviewers — CAGEO-D-26-00847

**Manuscript:** Hydrosheaf: A Reproducible Computational Framework for Inferring Directed Groundwater Hydrochemical Evolution Networks from Sparse Multitracer Data

**Journal:** Computers and Geosciences

---

## Addendum — Response to the second-round review (appended)

We thank the reviewer for the careful re-assessment. The following concrete corrections were applied in response to the second-round report; each is now reflected in the revised manuscript, the Supplementary Information, and the regenerated figures.

1. **Formal sheaf claim (R2-M1, the central methodological concern).** We corrected the earlier wording rather than preserving a misleading simplification. Section 2.4 and Supplementary Method S1 now distinguish edge selection from the retained-graph diagnostic. Edge selection is explicitly a weighted multi-criteria score. The diagnostic uses d-dimensional node stalks, affine edge restriction maps and a coboundary matrix D with right-hand side b; it reports the homogeneous nullity, affine obstruction energy and leave-one-edge-out leverage. No graph or sheaf Laplacian is formed or spectrally decomposed. The text also states that the sheaf layer is a structural network-consistency diagnostic and is not claimed to outperform an equivalent weighted score for individual edge decisions.

2. **Figure/text metric discrepancies (R2-m1), and an explanation the reviewer is owed.** The reviewer is correct that the file reviewed showed Figure 3D with Point R² = 0.52 and Figure 5a with MAE = 183.3 y. On investigation the cause was not the analysis but the packaging: **the DOCX submitted for re-review was built from an earlier snapshot of the sources and therefore predated corrections that had already been applied.** The regenerated figures had carried R² = −0.03 and Median AE = 17.9 y before that DOCX was produced. We apologise for the wasted reviewer effort; that is precisely the failure mode the reviewer warned about, and we have removed the possibility of repeating it (item 13 below). Figure 3D reports pointwise isotope R² = −0.03 with edge-mean R² = 0.99, and Figure 5a reports median absolute error 17.9 y, both matching their captions.

   Two genuine figure/text discrepancies did survive into the current sources and are now fixed. (a) The Figure 5b annotation recomputed the parity R² from the plotted values, which are floored at 0.01 yr for the log axes, and so printed 0.97 where the canonical benchmark summary — and the text, Table 2 and Table 4 — report log10 R² = 0.960. The figure generator now reads `log10_r2` directly from the benchmark summary file, so the annotation cannot drift from the text again. (b) Section 4.2, Section 5.3 and the Figure 3 caption quoted an isotope "noise amplitude qA of 0.07‰" against the figure's "noise σ_Δ = 0.71 permil". The figure was right: the quantity is the noise on an isotope *difference* between the two nodes of an edge, √2 × the 0.5‰ per-node analytical sigma = 0.71‰. The 0.07 was a transposition introduced in our own audit table and inherited by the text; both are corrected, and the manuscript now states the derivation explicitly.

   Table 6 was regenerated. We also found and fixed a defect in the table generator that had made this table look worse than it is: the comparison between an edge's dominant process and its most stable PSI family was a substring test, so the label "Carbonates" never matched "carbonate dissolution" and edge Talensi_48_61 was reported as a mismatch when the two in fact agree. One of the six edges now has agreement between its dominant-process label and its most stable PSI family; the other five genuinely disagree. The note to Table 6 explains that the two columns answer different questions and that the disagreement is retained as a diagnostic of dictionary degeneracy, so no extent-level claim is made.

3. **PHREEQC Methods/Results contradiction (Part C).** Section 2.7 now describes only what was executed: a PHREEQC-compatible mass-balance proxy using locked saturation fields, with an explicit statement that live PHREEQC kinetic execution is not part of the current release (consistent with Sections 4.2 and 4.5).

4. **Corrupted reference entry (Part C).** The reference list is repaired: Lucas & Unterweger (2000), Harte (2021) and Jurgens et al. (2022) are distinct, alphabetised entries. Tracing the cause showed that the same fault had also eaten the author-and-year text out of five *in-text* citations — leaving fragments such as "(… release, DOI 10.5066/F7J102FK)" — and had broken a sentence in Code and data availability, which ran the Harte citation straight into the runtime discussion. All six sites are repaired. A bidirectional citation audit now confirms 68 reference entries, 68 cited, and no in-text citation without a matching entry.

5. **Age-overlap denominator (Part C; R1-M12).** Supplementary Table S12 now states that the audit runs on the 13 directed edges of the synthetic benchmark network in each of 100 realisations (1,300 edge-checks), and that the reported fractions are means over those realisations. The incorrect "field age network" label was corrected to "synthetic benchmark age network".

6. **Graph-construction parameters and the distance-decay length scale (R1-M7).** Pressed a second time on this point, we went back to the source rather than to our own previous description, and found that our Methods text was wrong in three respects. Section 2.3 is rewritten to describe what the code executes. (i) Candidates are **not** ranked by haversine distance: they are ranked by the directional confidence p_ij = Φ((h_i − h_j)/σ_Δh), with distance only as the tie-break. That expression is now given explicitly in Section 2.3. (ii) The search radius **does** have a fixed default — `edge_radius_km = 5.0 km` — contrary to our previous statement that it had none; it is now listed in Table S3 with its default and its basis. (iii) The Gaussian distance-decay term the reviewer asked about is real, but it lives in the optional three-dimensional probabilistic builder, where its length scale is not a separate σ_dist at all but the same `edge_radius_km`: P(d) = exp(−d²/2r²), giving P = 0.61 at one radius and 0.14 at two. The manuscript now says so, and states that all reported results use the two-dimensional builder in which distance enters only as the radius cut-off and the ranking tie-break, never as a multiplicative decay. Table S3 additionally now carries `edge_p_min` (0.75), `edge_gradient_min` (10⁻⁴), `edge_depth_mismatch` (20 m), and the per-tier head sigmas (0.5 / 1.0 / 10.0 m) that set σ_Δh. We are grateful for the persistence on this comment; it exposed a Methods description that would not have permitted independent repetition.

7. **Transport-model selection consistency (R1-M9b).** Section 2.6, Table S3, and Supplementary Method S2 now state the identical rule: both models are fitted and the better-fitting candidate minimises the combined weighted objective (transport + residual chemistry + L1/EC-TDS/isotope/kinetic penalties).

8. **Null-edge wording (R1-M13).** The seven sink-terminus edges are now described consistently as unresolved/failed fits (chemistry R² = 0.0, no reaction set explains the residual) in both the main text and the SI, rather than as evidence of successful null closure.

9. **PSI pair reconciliation and caveat (R1-M15).** Supplementary Table S10 and Sections 4.6/5.5 now quote the site-aggregated per-mineral PSI identical to Figure 7 (CaNa_exch 0.73 vs NaCa_exch 0.46 at Lower Anayari; 0.83 vs 0.42 at Talensi; both carbonates near zero, values regenerated on the corrected field topology), and we added the explicit caveat that PSI measures stability to input perturbation, not the correctness of process attribution.

10. **Reaction-correlation and parameter-uncertainty diagnostics (R2-M3).** A new Supplementary Table S8b reports the empirical correlation matrix of recovered extents (albite–pyrite −0.72; CaNa_exch–pyrite +0.69; calcite–dolomite +0.66; etc.) and per-reaction mean/SD of recovered extents over the 100 realisations, addressing the missing uncertainty diagnostic.

11. **Claim softening (Part D).** The abstract, Section 5.1, Section 5.3 and the Conclusions were revised so that the supported claim class is the *stability ranking of process-family hypotheses under the specified perturbation model* — explicitly not reliable identification of the true reaction family — and the field-management sentence now says high-PSI hypotheses can be *prioritised for targeted confirmation and monitoring design*, not management-level process attribution. Section 5.4 and Supplementary Table S5 now explain that the field "direction" is an elevation proxy (no hydraulic-head data) and that retained "upgradient/flat" edges survive on chemical ordering alone.

12. **Release/DOI (R2-m2).** The manuscript continues to state the pinned commit (463e1ce), environment, tests, CI and reproduction scripts; minting the archived-release DOI and tagging the public release remain actions that the authors will perform upon acceptance, and this is stated explicitly.

13. **Two further corrections the re-review prompted, in Supplementary Table S5.** The reviewer asked how retained "upgradient/flat" field edges could carry an edge confidence of 1.00 at sites with no hydraulic-head data. The answer turned out to be that the column was mislabelled: the pipeline writes the **process-stability index** into it — the same quantity as the "PSI probability" column of Table 6 — and the field results carry no edge-confidence quantity at all, because the directional confidence p_ij requires head or elevation uncertainties these sites do not supply. The column is renamed "Process stability (PSI)" and the note now states plainly that 1.00 means the dominant reaction was selected in all 30 Monte Carlo trials and asserts nothing about hydraulic plausibility. Separately, the six edges listed in Table S5 did not reproduce from the canonical run; they are replaced by the six highest-ranked edges by rank score (chemistry R² × PSI), which are now identical to those in Table 6 — the ranking is strictly ordered with no ties, so the selection is reproducible.

14. **A programmatic audit, so that item 2 cannot recur.** The reviewer's central procedural criticism — that we made metric consistency a feature of the revision and then shipped figures contradicting the text — is fair, and we have addressed it with machinery rather than with another manual pass. Two checks now ship with the reproduction materials and are run as the last step of building the submission package:
    - `audit/number_audit.py` re-reads every canonical metric from its result file and asserts that each rendered form appears in the manuscript, the Supplementary Information and the generated tables; it also verifies that Table 6 and Table S5 agree edge-for-edge, that every Table 6 row appears verbatim in the manuscript, that the age-ordering percentages are integer counts over the stated 1,300-edge-check denominator, and that a list of superseded values cannot reappear. It currently makes 131 assertions and passes.
    - `audit/verify_docx.py` reads the built DOCX files back and asserts that the corrected strings are present, that the stale strings are absent, that each DOCX is newer than its markdown source, and that the expected figures are embedded. This is the check that would have caught the stale build described in item 2.

    Three supporting reproducibility defects were fixed at the same time: the Figure 3A strip-plot jitter was unseeded (so the figure was not byte-reproducible between runs, despite the paper's reproducibility claim) and is now seeded; the figure-embedding step was not idempotent and silently doubled every figure if run twice, and now refuses to run on a document that already carries figures; and Supplementary Figure S1 quoted the parity statistics at lower precision (log10 RMSE 0.24) than the main text (0.235), which is now aligned.

15. **Field graph-direction logic (third-round comment), and the sheaf construction.** The reviewer asked us to state exactly which computational branch produced the Ghana field edges and how they passed the filters stated in Section 2.3. We traced this in the source rather than in our own prose, and the answer is that they passed none of those filters, because none was applied. The field edges were built by a plain two-dimensional Euclidean k-nearest-neighbour construction over (longitude, latitude) with k = 3 including self; elevation was loaded into the sample record and then never used, so edge direction was index order rather than gradient. We confirmed this by reconstructing the edge set independently — it reproduces the reported graph exactly, 82/82 at Lower Anayari and 126/126 at Talensi. Instrumenting the pipeline confirmed the rest: the refinement stage reported `sheaf_refinement: not_requested` and returned exactly as many edges as it was given, so **no candidate was ever rejected**; 61 edges ran uphill on the elevation proxy, and 148 of the 208 contained their own reverse. The reviewer's suspicion was correct and understated.

    We have corrected this at the source rather than in the wording. The field tier now runs through the construction Section 2.3 describes — the probabilistic builder with elevation as the head proxy, then the refinement — and every field number is regenerated. The graph stage now produces 572 candidates and retains 258 (121 Lower Anayari, 137 Talensi), rejecting 314. **No retained edge runs against the elevation proxy at either site**, reciprocal pairs fall from 148 to 56, and every retained edge carries a real directional confidence p_ij (median 0.93 at Talensi; 0.50 at Lower Anayari, which is exactly Phi(0) — station elevations there are recorded as constants, so most Lower Anayari edges fall below the gradient floor and are retained explicitly as lateral mixing candidates). Median chemistry R2 is now 0.70 overall (Lower Anayari 0.53, Talensi 0.82) and median PSI 0.97. The seven zero-closure sink edges that prompted R1-M13 do not survive the directional refinement and are gone; fourteen edges now return a negative chemistry R2 and are reported as poor fits rather than removed. An automated check now asserts that no retained field edge opposes the elevation proxy, so this cannot silently regress.

16. **What the sheaf construction actually is (R2-M1, revisited).** Checking the code showed that the previous answer was wrong in both directions. There is no eigendecomposition anywhere in the framework, so the earlier spectral-localisation claim is removed. The implemented diagnostic is an affine cellular-sheaf coboundary: each node has a d-dimensional stalk, each directed edge has the restriction map alpha_e I and offset delta_e, and the stacked matrix D and right-hand side b encode the affine relations. The reported H0 is dim ker D, the homogeneous nullity; it is not automatically the dimension of an affine global-section set. The affine system has an exact global section only when min ||Dx - b|| squared is zero within numerical tolerance. The field obstruction energies are positive: 0.295 at Lower Anayari and 0.128 at Talensi, so neither field network has an exact affine global section. H0 = 0 at Lower Anayari means a trivial homogeneous nullspace, while H0 = 10 at Talensi means a ten-dimensional homogeneous nullspace, not a ten-dimensional family of exact affine assignments. Localisation is by leave-one-edge-out obstruction leverage, not by eigenvectors. Edge selection remains a weighted multi-criteria score, as acknowledged in response to R1-M8.

17. **A coordinate error found while testing the field elevations, and a reproducibility fix.** Following the reviewer's question we checked the field coordinates against independent elevation data rather than against our own description. Two further defects surfaced. First, **the Talensi longitudes were recorded with the wrong sign**, a transcription error the authors had already identified when preparing the site map and which is now corrected at source in the archived data (+0.618 to +0.815 rather than -0.618 to -0.815), placing the site in Togo rather than in the Upper East Region; sampling SRTM at the recorded coordinates disagreed with the surveyed elevations by -80.7 m, and at the corrected coordinates by -5.5 m, which is ordinary SRTM performance. The coordinates are corrected. No reported result changes, and we verified this rather than assuming it: pairwise great-circle distance depends on longitude only through sin^2(dlon/2), which is even, so a uniform sign flip mirrors the layout without altering any inter-well distance; re-running the pipeline returns an identical edge set with a maximum difference of 5.7e-08 across all numeric columns. Any geolocated figure is regenerated. Second, the process-stability Monte Carlo drew from an unseeded global generator, so the reported PSI family distribution drifted between identical reruns; both scripts are now seeded and PSI output is byte-reproducible.

18. **We tested whether a DEM could resolve the Lower Anayari directions, and it cannot.** Because Lower Anayari elevations are village constants, we sampled per-well elevations from SRTM and ASTER at the surveyed coordinates, using Talensi -- which has genuine per-well survey -- as a control to measure the DEM's real vertical error in this terrain. That error is **RMSE 7.91 m** (SRTM; ASTER is worse at 10.25 m), implying 11.2 m of noise on a pairwise elevation difference. The within-village DEM differences have a median of 6.5 m. At the elevation-tier sigma currently configured (1.0 m) the DEM would appear to "recover" **99 of the 100** lateral edges as confidently directed; measured against the DEM's own noise, **none** of the 100 has a difference exceeding twice that noise. We therefore rejected the DEM substitution: it would have manufactured a directed network from differences that are demonstrably below the resolution of the data. The village-mean DEM does agree with the recorded values (r = 0.773), so the recorded elevations are approximately correct, merely coarse. Lower Anayari is accordingly reported as a directed sub-network on its between-village edges plus within-village edges retained as undirected lateral mixing candidates, and the analysis is archived in `M2/m2_benchmark/scripts/revision/rev_dem_lateral_recovery.py`.

---

We thank the Editor and both reviewers for their careful reading and constructive comments. The revision addresses every point raised. The most consequential change is a full provenance and consistency audit of every reported metric: all results were regenerated from a single locked pipeline (benchmark scripts, fixed configuration, pinned commit), and the manuscript text, tables, and figures now quote identical values throughout (see the metric audit in `manuscript/revision/metric_audit.md`). Where the audit changed a reported number, we state the corrected value explicitly below.

We also discovered and corrected a serious reporting error in the submitted version: the no-prior topology result (F1 = 0.86) reported in Section 4.3 had no surviving computational source in the repository or its history (a provenance audit found no code path producing 215/168/47/6). The revision replaces it with the reproducible result (F1 = 0.62) and reframes the topology claims accordingly (see our response to R1-M11 and R2-M2).

---

## Editor comments

### EDITOR-1 — Editor summary of reviewer concerns

> The two reviewers agree that the manuscript addresses an important problem in hydrogeochemistry and view positively the integration of graph-based pathway inference, inverse geochemical modelling, thermodynamic screening, uncertainty analysis, and reproducible software practices into a single workflow. However, both reviewers raise substantial concerns regarding the methodological foundations and validation of the framework. A central issue is the role of the proposed sheaf-based formulation, which is currently not defined with sufficient rigor. The reviewers request a clearer explanation of what the sheaf framework contributes beyond a weighted multi-criteria consistency score and whether it is essential to the methodology. The validation strategy also requires strengthening. In particular, the topology validation appears partly circular because MODPATH information is used both as prior information and as a benchmark for evaluation. Reviewers recommend giving greater emphasis to the no-prior results and clarifying precisely what aspects of connectivity are being recovered. Additional concerns relate to the treatment of non-uniqueness in inverse geochemical modelling, the simplicity of the synthetic benchmarks, and the limited independent validation of the Ghana case studies. Finally, the reviewers identify a number of inconsistencies between figures, tables, and reported metrics, and request stronger documentation of default parameters, reproducibility materials, and computational performance.

**Response:** We address all six themes. (1) The sheaf construction is now defined precisely (stalks, restriction maps, consistency residual, Laplacian localisation) in Section 2.4 and Supplementary Method S1, with an explicit statement of what it adds over a weighted multi-criteria score (R1-M8, R2-M1). (2) The topology validation is de-circularised: the prior-assisted F1 = 1.00 is relabelled a physics-prior ingestion-fidelity check, and the independent no-prior result (F1 = 0.62) is now the primary topology claim, with baseline graph-construction comparators (R1-M11, R2-M2). (3) Non-uniqueness in inverse modelling is quantified with new identifiability diagnostics (dictionary rank, condition number, collinearity, leave-one-out sensitivity) and discussed with respect to PSI behaviour on degenerate pairs (R1-M15, R2-M3). (4) The synthetic benchmark now includes a multi-endmember mixing-with-reactions stress test and honest reporting of extent-level recovery limits (R2-M4). (5) The Ghana demonstration is strengthened with site-specific context and an explicit evidence-availability statement (R2-M5). (6) A full audit reconciled every metric across text, tables, and figures (R1-M10, R2-m1); all parameter defaults are documented in Table S3; reproducibility materials (pinned commit, environment, tests, scripts, runtime scaling) are described in the Code and data availability section (R2-m2).

---

## Reviewer 1

### R1-M1 — Hyphenation of compound modifiers

> 1. Lines 40, 53, 80, and throughout. The terms "semi arid", "non uniqueness", "data limited", "up gradient", "down gradient", and "time averaged" appear without hyphens throughout the manuscript. Please hyphenate these compound modifiers consistently, as in "semi-arid", "non-uniqueness", and "data-limited".

**Response:** Agreed. All compound modifiers are hyphenated consistently throughout the revised manuscript and supplementary material (semi-arid, non-uniqueness, data-limited, up-gradient, down-gradient, time-averaged). A grep audit confirms no unhyphenated instances remain.

**Manuscript location:** throughout.

### R1-M2 — Distinction from combinatorial inverse modelling

> 2. Lines 48-51. The authors state that PHREEQC and NETPATH depend entirely on the analyst's prior selection of an initial and final water pair assumed to represent a connected upgradient-downgradient pathway. This is accurate and well-cited. However, there is no mention of existing approaches that have attempted to partially automate or systematise this selection, such as combinatorial inverse modelling as implemented by Manu et al. (2023), already cited in the reference list. How does Hydrosheaf's graph-based automation differ from or improve upon combinatorial approaches in practice? Please elaborate on this distinction, as in here it is central to the novelty claim.

**Response:** Agreed, and we thank the reviewer for the opportunity to sharpen the novelty claim. The revised Introduction and Section 5.2 now make the distinction explicit. Combinatorial inverse modelling (Manu et al., 2023) systematises the *pair selection* step by enumerating candidate initial–final water pairs and retaining those with chemically feasible solutions. Hydrosheaf differs in four practical respects: (i) candidate edges are scored against hydraulic, conservative-tracer, isotope, and age-order evidence jointly, so pair selection is constrained by multiple independent evidence streams rather than by mass-balance feasibility alone; (ii) non-reactive transport (evaporation or mixing) is separated from residual chemical change before reaction fitting, so the inverse problem is posed on the residual rather than on raw concentration differences; (iii) all inferred reactions are screened by PHREEQC-derived saturation-index gates and ranked by a Monte Carlo process-stability index (PSI); and (iv) the result is a directed network with per-edge confidence, uncertainty-ranked process hypotheses, and full provenance logging, rather than a set of independently feasible pair solutions.

**Manuscript location:** Section 1 (Introduction), Section 5.2.

### R1-M3 — Acknowledging wrapper-based combinations

> 3. Lines 62-73. The claim that Hydrosheaf fills a gap not addressed by any single existing tool is well-argued in general terms. However, existing Python-based scripted workflows combining PHREEQC wrappers (e.g., phreeqpy or iPhreeqc bindings) with graph analysis libraries (e.g., NetworkX) could, in principle, partially replicate components of this framework. The authors should briefly acknowledge whether such combinations have been attempted in the literature and why they do not constitute the integrated reproducible solution that Hydrosheaf provides.

**Response:** We acknowledge the point. The revised Introduction and Discussion now state that Python wrapper combinations (e.g., phreeqpy or iPhreeqc bindings with NetworkX) can in principle replicate individual components. To our knowledge, no published, versioned workflow combines transport decomposition, sparse inverse reaction fitting, thermodynamic gating, uncertainty diagnostics, and a multi-tier validation suite under a single configurable pipeline, and we have found no published report of such an integrated combination; the revised text states this carefully and restricts the claim to what we can support.

**Manuscript location:** Section 1, Section 5.2.

### R1-M4 — Reframing contribution six

> 4. Lines 87-100. The six computational contributions are clearly stated and form an effective summary of novelty. I would encourage the authors to revisit contribution six, which claims multi-tiered reproducible benchmarks as a computational contribution. Reproducible benchmarking is good scientific practice but it is not, in itself, a methodological contribution to the framework. Consider reframing this more as a validation strategy rather than listing it alongside the algorithmic contributions.

**Response:** Agreed. The contribution list now comprises five algorithmic contributions (directed-graph formulation; sheaf-style consistency scoring; sparse inverse fitting with thermodynamic gates; process-stability indexing; transport–reaction decomposition), followed by a separate statement describing the multi-tier validation strategy as an evaluation design rather than an algorithmic contribution.

**Manuscript location:** Section 1 (end of Introduction).

### R1-M5 — Consistent HCO3 notation

> 5. Lines 131-136. HCO3 is written without a charge superscript in several places in the main text and appears as HCO3- in others. Please ensure consistent chemical notation throughout. The same applies to ionic species in the stoichiometric dictionary description.

**Response:** Agreed. Charged species are written consistently with superscript charges throughout the main text, equations, and the stoichiometric dictionary description (HCO3⁻, SO4²⁻, NO3⁻, Ca²⁺, Mg²⁺, Na⁺, K⁺, Cl⁻).

**Manuscript location:** throughout, including Section 2.2 and Table S2.

### R1-M6 — meq/L vs mmol/L in the charge balance error

> 6. Lines 135-137. The charge balance error (Equation 1) uses concentration c in meq/L in the denominator and numerator, but line 135 states that concentrations are converted to mmol/L as the common unit basis. This creates a notational inconsistency that must be resolved. Are the charge balance calculations performed in meq/L on unconverted data, or after unit conversion? Please clarify.

**Response:** Agreed, and clarified. Concentrations are converted to mmol/L as the pipeline's common unit basis. Charge balance requires *equivalent* concentrations, so the CBE is evaluated on equivalent concentrations in meq/L (meq/L = z × mmol/L) computed from the converted data. Equation (1) and the surrounding text now state this convention explicitly.

**Manuscript location:** Section 2.2, Equation (1).

### R1-M7 — Default distance threshold

> 7. Lines 143-156 (Section 2.3). The configurable distance threshold used to generate candidate directed edges is described but no default value or recommended range is provided anywhere in the main text. This is a critical parameter, since too large a threshold generates physically implausible long-distance edges, while too small a threshold may miss valid connections in sparse datasets. What is the default value, and what is the scientific basis for its choice?

**Response:** Agreed, and we have corrected the description to match the implementation. Candidate-edge construction uses a configurable maximum-neighbour rule (default `edge_max_neighbors = 3`) with a Gaussian distance-decay term (configurable scale), a minimum accepted edge-confidence (default 0.75), and a minimum hydraulic gradient (default 10⁻⁴). The default of three downstream neighbours per node is a deliberate compromise: it bounds candidate counts to ~O(3n) (measured runtime scaling in the SI), avoids the physically implausible long-range edge sets that large thresholds generate, and retains enough candidates to avoid missing valid connections in sparse networks. The full default set and rationale are listed in Table S3.

**Manuscript location:** Section 2.3, Table S3, SI runtime table.

### R1-M8 — Rigorous sheaf definition and added value

> 8. Lines 157-176 (Section 2.4). The sheaf-inspired consistency framework is presented as the conceptual centrepiece of the paper... It is not clear from the text what the sheaf contributes beyond a weighted multi-criteria consistency check that could be described in more conventional terms. Could the authors clarify, in language accessible to a hydrogeology readership, what the sheaf framework adds in practice over a simpler weighted scoring of these same four criteria? If the benefit is mainly conceptual or organisational rather than producing materially different edge decisions, it would be helpful to state that plainly.

**Response:** We thank the reviewer for this challenge and have rewritten Section 2.4 accordingly, with the full construction in Supplementary Method S1. The text now separates edge selection from the retained-graph diagnostic. Edge selection is a weighted multi-criteria score. The diagnostic defines node stalks, affine restriction maps, the coboundary matrix D and its right-hand side b, and reports the homogeneous nullity, affine obstruction energy and leave-one-edge-out leverage. No graph or sheaf Laplacian is formed or spectrally decomposed. We state plainly that the added value is structural and network-level, not a claim that individual edge decisions outperform an equivalent weighted score.

**Manuscript location:** Section 2.4, Supplementary Method S1.

### R1-M9 — Defaults for epsilon, transport-model selection, lambda_2

> 9. Lines 179-224. Several parameters that affect the results are introduced without a stated default value or a basis for their choice... The most important is the tolerance epsilon in the age-ordering constraint (line 189-190)... In the transport decomposition (lines 199-207), it is also unclear how the choice between the evaporation model and the mixing model is made... To a lesser extent, the regularisation in Equation 6 would benefit from clarification of whether lambda_2 is fixed or optimised alongside lambda_1 and what its default value is. I would ask the authors to state the defaults for these parameters, ideally collected in Table S3, and to give a brief rationale for the age-ordering tolerance and the transport-model selection in particular.

**Response:** Agreed. All defaults are now stated in the text and collected in Table S3. (i) Age-ordering tolerance: the ordering criterion is evaluated as τ_j ≥ τ_i + ε with ε = 0 years by default (any downstream-age deficit counts as a violation candidate), while posterior intervals are used to grade severity: violations whose intervals overlap are flagged "not resolved at stated uncertainty" and are not counted as severe; only non-overlapping reversals beyond a log10 severity threshold of 0.3 receive a severe age-coherence flag; the audit does not apply a continuous overlap-proportional weight or an automatic refinement penalty. (ii) Transport-model selection: both the evaporation (γ) and mixing (f) models are fitted per edge, and the better-fitting candidate is selected algorithmically as the one minimising the combined weighted objective — the transport and residual-chemistry mismatch plus the L1, EC/TDS, isotope, and kinetic penalties; the transport fit itself is weighted towards the conservative anchor constraints (EC/TDS, Cl⁻, Br⁻). The criterion is stated in Section 2.6. (iii) λ₂: the Tikhonov coefficient is fixed (default 0.0) with a numerical ridge floor added for stability; it is not optimised alongside λ₁, which is selected by AICc over a grid.

**Manuscript location:** Sections 2.4-2.6, Table S3.

### R1-M10 — R2 inconsistency (0.86 vs 0.74)

> 10. Lines 316-332. The main text reports recovery R2 = 0.86 for active reaction extents, but Panel B of Figure 3 shows R2-active = 0.74 in the figure annotation. These values are inconsistent and must be reconciled. Please check all reported metrics against the figures and tables and correct any discrepancies before resubmission.

**Response:** Agreed, and this triggered the full metric audit. Investigation showed that the submitted text (R² = 0.86) had no surviving computational source, and the figure annotation (0.74) came from an analysis snapshot produced before a transport-selection correction (commit bd7c8aa, May 2026) that removed an artificial bias favouring the evaporation model. All benchmark results were regenerated from the current locked pipeline, and the canonical value for active-reaction recovery is now correlation R² = 0.23 with MAE = 0.37 mmol/L and RMSE = 0.62 mmol/L across 2,100 active reaction instances, with a false-activation rate of 54.1% (Section 4.2, Table 2, Figure 3B — identical values everywhere). The revision reports this honestly and explains the cause with the identifiability diagnostics requested in R2-M3: the reaction dictionary is rank-deficient (rank 8 of 11 ions, infinite condition number), so extent-level attribution is inherently limited. We recognise this is a substantial reduction from the submitted claim; it reflects the corrected pipeline, and we thank the reviewers for prompting the audit.

**Manuscript location:** Section 4.2, Table 2, Figure 3B, Section 5.5.

### R1-M11 — Topology validation reframing

> 11. Lines 335-347. A couple points need addressing here. First, the perfect topology result (F1 = 1.00) is obtained in prior-assisted mode, where MODPATH endpoint pairs are supplied as inputs and then recovered as outputs. This is closer to a pipeline integrity check than to an independent test of inference, so the no-prior result (F1 = 0.86) should carry the primary emphasis. Second, please be precise about what the benchmark demonstrates. The agreement with MODPATH establishes recovery of well-to-well endpoint connectivity, not the pathline geometry or travel time that MODPATH produces... Finally, the specific USGS MODPATH model or dataset used as the reference is never identified. What model is this, and is it a published, publicly accessible dataset? Please provide the full citation and confirm availability.

**Response:** Agreed on all three points. (a) The prior-assisted F1 = 1.00 is now labelled explicitly as a physics-prior ingestion-fidelity check, and the independent no-prior result carries the primary emphasis. We must also correct the submitted no-prior number: a provenance audit found no computational source for F1 = 0.86 (215 candidates / 168 TP / 47 FP / 6 FN) in the repository or its history. The reproducible result, identical in method to the dedicated M4 topology benchmark, is 302 candidates with TP = 147, FP = 155, FN = 27, precision = 0.49, recall = 0.84, F1 = 0.62 (head-gradient construction: elevation-as-head, downhill, k = 2). This is a high-recall result with substantial overconnection (FP > TP), reflecting the benchmark's information limit: in the outlet-convergent Savage system, true positives and false positives converge on the same outlet cells, and the distinction depends on internal permeability and pumping that coordinates and head proxies cannot carry. Inferred edges are therefore screening-level hypotheses requiring corroboration. (b) We now state precisely that agreement with MODPATH demonstrates recovery of well-to-well endpoint connectivity, not pathline geometry, travel time, or porosity-dependent transport. (c) The reference is identified and cited: the USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5 archive (Harte, 2021, U.S. Geological Survey data release, https://doi.org/10.5066/F7J102FK), publicly available. Baseline comparators are provided in the SI (all-pairs elevation-drop F1 = 0.47; proximity kNN F1 = 0.00; conservative-tracer ordering not evaluable on this archive because it contains no hydrochemistry).

**Manuscript location:** Sections 3.3, 4.3, Table 5, Figure 2, SI (baseline table), References.

### R1-M12 — Wide age ranges and the age-ordering constraint

> 12. Lines 349-362. The inferred age ranges for the older waters are very wide... What is missing is how these wide ranges feed back into the age-ordering constraint in Section 2.4. If two connected nodes have age estimates that overlap by this much, the ordering criterion becomes hard to apply. Please explain how the framework handles this case. Is the edge retained, rejected, or flagged with reduced confidence?

**Response:** Thank you for this question; it addresses a real design point and the framework's behaviour is now described and quantified. The ordering criterion is evaluated on posterior intervals rather than point estimates. When the intervals of a connected pair overlap, the edge is retained and the violation is flagged as "not resolved at stated uncertainty"; such cases are excluded from the severe-violation set. Non-overlapping reversals beyond the log10 severity threshold of 0.3 receive an age-coherence failure flag. The current audit reports these flags but applies no continuous overlap-proportional weight and no automatic exclusion or penalty to the field results. Quantified on the synthetic network, 15.2% of directed edges violate the point-estimate ordering, 84.3% of those violations are interval-overlap cases, and 2.4% of all edges are severe non-overlapping reversals. The age-order consistency index is 0.85.

**Manuscript location:** Section 2.5, Section 4.4, SI (age-overlap statistics).

### R1-M13 — Zero reaction extents in Table 6

> 13. Lines 382-400. Table 6 reports reaction extents of exactly 0.00 mmol/L for three Lower Anayari edges and one Talensi edge, yet with R2 = 1.00 and PSI = 1.00 in every case. A reaction extent of zero means the solver found no net chemical change along these edges, so a PSI of 1.00, which indicates maximum robustness, is difficult to interpret. Does PSI = 1.00 here mean the framework robustly detected that no reaction occurred, or is there an error in how the table values were computed? Please clarify this and revise Table 6 accordingly.

**Response:** The reviewer identified a real defect in the submitted Table 6, and the audit traced its cause. The seven zero-closure rows belonged to the superseded unfiltered as-run candidate graph, not to the canonical documented field output. After the directional refinement, the canonical run contains 572 candidates, 258 retained edges and zero unresolved null edges. Table 6 therefore lists six explained, top-ranked Talensi edges; one edge has agreement between the dominant extent label and the most stable PSI family, while five genuinely disagree. The table note explains the distinction and retains the disagreement as a dictionary-degeneracy diagnostic rather than treating PSI = 1 as evidence of a zero reaction.

**Manuscript location:** Section 4.6, Table 6, Table S5, Section 5.5.

### R1-M14 — MODPATH: benchmark vs prior tension

> 14. Lines 419-437. The characterisation of MODFLOW and MODPATH is broadly fair... However, the argument contains a tension that should be addressed. MODPATH is described here as too demanding for the target settings, yet it is used in Section 4.3 as the reference truth for topology validation, and the optional MODPATH priors are precisely what produce the perfect topology result. The same tool cannot be the standard the framework is meant to replace and the standard it is validated against without this being made explicit... Please reframe the comparison so that the strength of the topology claim matches what an inferred edge can actually support, and clarify the role of MODPATH as benchmark versus prior.

**Response:** Agreed, and the tension is now addressed explicitly in Sections 4.3 and 5.2. The position is: MODPATH's calibrated-input demands (conductivity field, boundary conditions, effective porosity) motivate Hydrosheaf's design for settings where those inputs do not exist. Where a MODPATH archive already exists, it plays two distinct, clearly separated roles: (i) as an optional physics prior, its ingestion is verified as a pipeline-integrity check (F1 = 1.00, labelled as such and explicitly not independent inference); and (ii) as an independent connectivity benchmark, the no-prior inference (F1 = 0.62) is the actual test of the framework's inference. We also now state that an inferred edge is conceptually weaker than an advective pathline: agreement with MODPATH endpoints demonstrates well-to-well connectivity, not pathline geometry, travel time, or porosity-dependent transport (cf. Meyer et al., 2018; Baker et al., 2025).

**Manuscript location:** Sections 4.3, 5.2.

### R1-M15 — Does PSI separate degenerate pairs?

> 15. Lines 471-486. You acknowledge that residual non-uniqueness persists for mineral pairs with similar stoichiometric signatures, such as calcite versus dolomite or Ca-Na versus Na-Ca exchange... What the discussion does not address is whether the PSI metric helps in practice. Does it separate these degenerate solutions, or does it assign similarly high stability to both members of the pair? This is a practically important point and should be discussed.

**Response:** We now answer this empirically. A dedicated diagnostic (SI Table) computed per-reaction Monte Carlo inclusion probabilities on the field edges for both degenerate pairs. For CaNa_exch versus NaCa_exch, PSI **does separate** the members: the mean absolute PSI difference per edge is 0.85 (Lower Anayari) and 0.76 (Talensi), with the dominant member at PSI ≈ 0.69-0.75 versus 0.39-0.45 for the other — the perturbation envelope consistently prefers one exchange direction per edge. For calcite versus dolomite, both PSI values are near zero on the sampled field edges (carbonates are not stably activated there), which is consistent with the field PSI family distribution (Carbonates 8/208 edges). The Discussion states the practical implication: PSI separates *directional* degeneracies (exchange pairs) but cannot rescue *compositional* degeneracies where neither member is active; extent-level attribution for near-collinear dissolution pairs remains limited and is quantified by the identifiability diagnostics.

**Manuscript location:** Section 4.6, Section 5.5, SI (PSI pair-separation table).

### R1-m1 — HCO3- in equations

> 1. Line 131. Write HCO3- consistently with charge superscript throughout the text and in all equations.

**Response:** Agreed; done throughout text and equations (see R1-M5).

### R1-m2 — 3H half-life citation

> 2. Line 180. The 3H half-life is stated as 12.32 yr. Please cite the primary source for this value rather than stating.

**Response:** Agreed. The value now cites the primary source: Lucas, L.L., & Unterweger, M.P. (2000), Comprehensive review and critical evaluation of the half-life of tritium, Journal of Research of the National Institute of Standards and Technology, 105(4), 541-549, https://doi.org/10.6028/jres.105.043 (4500 ± 8 days ≈ 12.32 yr).

**Manuscript location:** Section 2.1.2, References.

### R1-m3 — Analytical sigma defaults

> 3. Lines 247-248. The default analytical relative sigma of 4% for major ions and 0.5‰ for stable isotopes are stated without citation. Are these values taken from standard analytical practice in the literature, or are they arbitrary defaults? Please either cite a basis for these values or note explicitly that they are configurable defaults.

**Response:** Clarified. These are configurable defaults representing typical routine analytical precision for major-ion and stable-isotope measurements reported in the cited field studies; they are not presented as universal analytical constants. The text and Table S3 now state this explicitly.

**Manuscript location:** Section 2.8, Table S3.

### R1-m4 — Provenance manifest schema

> 4. Lines 265-272. The provenance manifest is described as a JSON file. It would be helpful to include an example or schema for this file in the Supplementary Information, as it is central to the reproducibility claim.

**Response:** Agreed. The SI now includes the JSON provenance-manifest schema and a real example generated by the canonical run.

**Manuscript location:** Supplementary Information (provenance manifest schema), Section 2.9.

### R1-m5 — Table 2 Ghana row

> 5. Table 2. The Ghana field demonstration row reports median R2 = 0.99 and PSI = 1.00 as the main metric. In light of the zero reaction extents discussed in the major comments, these values as currently reported are not interpretable at face value and should be revised or accompanied by an explanatory note.

**Response:** Agreed. The Ghana row now reports the corrected canonical values: median chemistry R² = 0.71 (Lower Anayari 0.60, Talensi 0.77) over 208 edges and median PSI = 1.00, with a note that seven null edges (unresolved, R² = 0.0) are excluded from process discovery (see R1-M13).

**Manuscript location:** Table 2, Section 4.6.

### R1-m6 — Compound noun strings

> 6. Throughout. The writing makes frequent use of compound noun strings such as "multi evidence consistency residual", "process interpretable mass transfers", and "uncertainty aware hydrochemical process networks". While not grammatically incorrect, some of these constructions are difficult to parse on first reading. Please review for clarity, particularly in the abstract and conclusions.

**Response:** Agreed. The abstract, conclusions, and discussion were rewritten with plainer constructions (e.g., "network-wide consistency residual", "process-attributable mass transfers", "uncertainty-ranked hydrochemical process networks"), and a readability pass was applied across the manuscript.

**Manuscript location:** Abstract, Conclusions, throughout.

---

## Reviewer 2

### R2-M1 — Precise sheaf definitions

> 1. The manuscript uses the terminology of sheaves, restriction maps, and sheaf consistency, but the actual implementation appears to be closer to a weighted multivariate residual or affine consistency score between candidate upstream and downstream nodes. The authors should define precisely what constitutes the sheaf, what the stalks and restriction maps are, how they are estimated or prescribed, and how the sheaf Laplacian is constructed and used in the algorithm. At present, the term "sheaf" risks sounding nominal rather than technically essential.

**Response:** Agreed; see our response to R1-M8. Section 2.4 and Supplementary Method S1 now define the implemented affine coboundary diagnostic precisely and state that no Laplacian is constructed. H0 is the homogeneous nullity, while affine solvability is tested separately by the obstruction energy. The field obstruction energies are positive, so neither field network is reported as having an exact affine global section.

**Manuscript location:** Section 2.4, Supplementary Method S1.

### R2-M2 — No-prior emphasis and baselines

> 2. The topology validation is not fully convincing. In the MODPATH prior-assisted case, MODPATH-derived information appears to be used both as a prior and as the reference for evaluation, which may lead to circular validation. The no-prior case should be emphasized as the primary test, and additional baseline comparisons should be provided, such as hydraulic-gradient-only, distance-based, elevation-based, or conservative-tracer-based graph construction methods.

**Response:** Agreed on both counts. (i) The no-prior case is now the primary topology test (F1 = 0.62, reproducible, identical method to the M4 benchmark), and the prior-assisted case is labelled an ingestion-fidelity check. (ii) Baseline graph-construction rules are compared against the same 174-edge reference in the SI: elevation-as-head downhill k = 2 (the framework construction, F1 = 0.62), all-pairs elevation-drop (F1 = 0.47, same recall 0.84 but double the false positives), proximity-only k-nearest k = 2 (F1 = 0.00, 306 edges, no true positives), and conservative-tracer ordering (not evaluable on the Savage archive, which contains no hydrochemistry; the framework's conservative-tracer evidence is evaluated in the synthetic benchmark instead). These comparisons show that the framework's neighbour-pruned, evidence-scored construction outperforms naive alternatives and that proximity alone carries no flow-connectivity information in this system.

**Manuscript location:** Sections 3.3, 4.3, Table 5, SI (baseline comparison table).

### R2-M3 — Identifiability diagnostics

> 3. The inverse geochemical modelling remains vulnerable to non-uniqueness and overfitting. The synthetic tests appear to use the same reaction dictionary for data generation and inversion, which may overestimate performance. The authors should provide more diagnostics on rank, condition number, parameter uncertainty, reaction correlation, and sensitivity to the reaction dictionary.

**Response:** Agreed, and these diagnostics are now reported (Section 4.2 and SI). (i) Dictionary diagnostics: the benchmark's 14-reaction × 11-ion dictionary has rank 8 (deficiency 3) and an infinite condition number; 10 column pairs have |cosine| > 0.7, including calcite~anorthite and CaNa_exch~NaCa_exch at |cosine| = 1.00. (ii) Reaction correlation: the empirical correlation matrix of recovered extents is provided in the SI. (iii) Dictionary sensitivity: leave-one-out refits over 5 locked realisations show that removing albite improves active-recovery R² from 0.09 to 0.17 and MAE from 0.40 to 0.34 mmol/L, but no single removal resolves the degeneracy, indicating a systemic identifiability limit rather than one problematic reaction. (iv) We also acknowledge the reviewer's point that generation and inversion share the reaction dictionary; the benchmark is therefore best interpreted as an identifiability stress test under idealised conditions, and the paper no longer claims near-perfect recovery. Honest extent-level performance (R² = 0.23, MAE 0.37 mmol/L, 54% false-activation rate) and family-level performance (R² = 0.16, dominant-family hit rate 48%) are reported.

**Manuscript location:** Section 4.2, Section 5.5, SI (identifiability tables), Figure 3B.

### R2-M4 — Transport correction realism

> 4. The transport correction uses simple evaporation or one-endmember mixing before fitting residual reactions. In real aquifers, mixing and reaction are coupled, and chloride or EC cannot always be treated as conservative anchors, especially in evaporite-bearing, mining-impacted, or anthropogenically affected systems. The authors should clarify how the framework distinguishes halite dissolution, evapoconcentration, agricultural chloride input, and mixing. The framework should also be tested against more realistic multi-endmember and reaction-mixing scenarios.

**Response:** Agreed. Three changes: (i) Section 2.6 now states the anchor-constraint caveat explicitly — Cl⁻ and EC are treated as conservative anchors only where halite dissolution, agricultural chloride input, or anthropogenic salinity is not indicated, and the transport-fit diagnostics surface cases where the anchor assumption fails. (ii) A new synthetic stress test (Section 4.2, SI) generates downstream chemistry from two-endmember mixing (f₁ = 0.20, f₂ = 0.15) with simultaneous halite (0.40 mmol/L) and calcite (0.20 mmol/L) extents; the single-endmember transport stage selects mixing on 73% of edges with median f = 0.24 (true total 0.35) and absorbs most of the reaction signal (median recovered halite extent 0.04 mmol/L, calcite 0.00) while still reproducing the chemistry (median R² = 0.999). This demonstrates quantitatively that when transport is misspecified, reaction extents are absorbed into the transport stage, and that halite dissolution, evapoconcentration, and mixing cannot be fully distinguished from major ions alone. (iii) The Discussion now states the practical consequences and recommends independent tracers (e.g., 87Sr/86Sr, boron or nitrate isotopes) to separate these processes.

**Manuscript location:** Section 2.6, Section 4.2, Section 5.5, SI (multi-endmember scenario).

### R2-M5 — Ghana site-specific evidence

> 5. The field demonstration is currently more a proof of workflow execution than an independent validation of hydrogeological interpretation. High chemical reconstruction performance does not by itself confirm the inferred flow paths or reactions. The Ghana cases should be supported by more site-specific evidence, including hydrogeological setting, hydraulic-head data, lithology/mineralogy, sampling density, isotope or trace-element constraints, and independent support for inferred processes such as ion exchange, gypsum dissolution, fluorite dissolution, denitrification, or nitrate input.

**Response:** Agreed, and the demonstration is repositioned and strengthened. Section 4.6 and the SI now provide the site context from the published record: the Lower Anayari alluvial/basement system with agricultural land use (Abdul-Wahab et al., 2021) and the Talensi mining-influenced crystalline basement system (Song et al., 2024), including lithology, land use, sampling densities (41 and 63 wells, with 121 and 137 retained graph edges in the canonical run), available tracers (major ions and stable isotopes; no nuclear tracers, no hydraulic-head records), and the absence of independent process-truth labels. The inferred process families (exchange-dominated at Lower Anayari; gypsum/fluorite and redox signals at Talensi) are presented as PSI-ranked hypotheses consistent with this setting, and the text states explicitly which independent evidence would be required to confirm them (repeat sampling, mineralogical analysis, head data, or dedicated tracer tests). We also note that the corrected median chemistry R² = 0.71 (not 0.94 as submitted) is reported, and that high reconstruction performance is explicitly not equated with process-truth validation.

**Manuscript location:** Section 4.6, Section 5.4, SI (Tables S4-S5), Discussion.

### R2-m1 — Internal consistency audit

> 1. Some results are internally inconsistent or insufficiently explained. Several metrics in the text, tables, and figure captions appear inconsistent. For example, different parts of the manuscript report different values for synthetic performance, Ghana field performance, and PSI. In Table 6, some reaction extents are reported as 0.00 while PSI is 1.00 and the interpretation refers to evaporite or carbonate signals; this needs explanation or correction. The authors should audit all tables, figures, captions, and textual summaries for consistency.

**Response:** Agreed; the full audit is summarised in the metric audit table (manuscript/revision/metric_audit.md) and is implemented throughout. Every reported metric was regenerated from a single locked pipeline, and text, tables, and figures now quote identical values; the generators were patched so tables and figures read the same canonical result files. The specific inconsistencies the reviewer lists were traced to mixed analysis snapshots (transport-selection correction; field mixing-endmember fix; M3 identifiability gate) and are corrected: synthetic reaction recovery (Section 4.2/Table 2/Fig. 3B all report R² = 0.23), Ghana performance (the canonical documented run reports median chemistry R² = 0.70, with 0.53 at Lower Anayari and 0.82 at Talensi), and Table 6 (rebuilt with zero unresolved null edges and one-agree/five-disagree process labels, see R1-M13). The Table 6 "interpretation" column now always matches the dominant-process column.

**Manuscript location:** throughout; metric audit in SI Appendix.

### R2-m2 — Software release, DOI, environment, tests, runtime

> 2. Because this is a computational framework paper, the authors should provide a versioned software release, DOI or archived repository, commit hash, environment file, test dataset, installation instructions, automated tests, and scripts reproducing every figure and table. Runtime and scalability should also be reported, since candidate edge construction may scale approximately with the square of the number of samples unless pruned.

**Response:** Agreed. The manuscript now distinguishes the public source repository from the exact analysis snapshot. The reported results are identified by commit 463e1ce on branch codex/m3-correctness, but that revision and the revised M2 field package were not present on the public origin when this package was frozen. The revision materials include a pre-release SHA-256 reproducibility manifest, the Python 3.14.6 environment export, locked configurations, test datasets, installation instructions, automated tests, CI workflow and scripts for reproducing the figures, tables and radius sensitivity. Runtime and scalability are reported in the SI. The authors must publish the exact snapshot as a versioned release and add a persistent DOI before resubmission; until then, complete public reproduction of the reported field package remains an author-held release action.

**Manuscript location:** Code and data availability, Section 2.9, SI (runtime table, environment, provenance schema).

---

We thank both reviewers again for the thorough and constructive reviews. All reviewer comments are reproduced above verbatim; the revised manuscript and response letter now describe the implemented diagnostics and canonical field output consistently; publication of the exact immutable release and DOI remains an author action.
