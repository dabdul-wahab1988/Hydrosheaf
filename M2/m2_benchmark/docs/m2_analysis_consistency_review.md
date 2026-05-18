# M2 Analysis Consistency Review

Updated: 2026-05-17

## Reviewer Gate

The M2 figure package is usable only with conservative captions. The synthetic benchmark, MODPATH topology check, M3 public-age screen, PHREEQC proxy, and Ghana field demonstration are different evidence tiers and must not be collapsed into one claim of full validation.

## Corrections Applied

- Fig. 5 now prefers `M3/m3_age_benchmark/results/m3_phase4_screened_full_results.csv` for the public USGS age panel, with the older M2 E1 CSV kept only as fallback.
- M3 was rerun at the canonical full setting: `--full --age-steps 35 --scenario screened_dgm_gases`.
- The Ghana field runner now configures a transparent local dilute mixing endmember, so field edges actually compare `evap` against `mix`.
- Fig. 3 wording was changed from proof language to synthetic benchmark evidence.
- Fig. 3B now reports all benchmark reactions and treats signed exchange extents as active when `abs(true extent) > 0`.
- Fig. 7 retains detailed process labels and removes guarantee language from the plot title.
- The figure/table evidence registers now point to the PNG Fig. 1 output, M3 public-age output, and the existing external validation plan rather than missing or removed artefacts.

## Observed Metrics And Reasons

### Fig. 3B Reaction Recovery

Current all-reaction recovery is moderate, not uniformly strong:

- active rows: `2100`
- inactive rows: `4900`
- active `R2`: `0.743`
- active MAE: `0.226 mmol/L`
- inactive false-positive rate above `0.05 mmol/L`: `35.3%`

The high false-positive rate is concentrated in inactive `albite` and `CaNa_exch` terms. This is expected from the current inverse setup because several reactions share ion signatures with carbonate/Na/Ca residual structure, while `lambda_l1 = 0.002` is permissive enough to keep small compensating extents. These rows should be interpreted as sparse-recovery leakage, not confident mineral identification.

### Fig. 3D Isotope Recovery

Pointwise isotope-shift recovery is weaker than edge-mean recovery because the synthetic isotope noise propagates into an edge-difference uncertainty of about `0.71 permil`. Many true shifts are zero, so pointwise noise dominates the scatter. The defensible claim is edge-mean/noise-aware consistency, not point-level isotope proof.

### Fig. 5 Public Age Validation

The refreshed M3 full screened public-age benchmark is weak-to-moderate:

- rows: `1272`
- median absolute log10 error: `0.383`
- log10 RMSE: `1.065`
- within factor 2: `0.410`
- within factor 10: `0.678`

The reason is not a missing run anymore. The result reflects heterogeneous public tracer support, ambiguous tritium/bomb-pulse solutions, young-gas contamination or correction sensitivity, and old-groundwater uncertainty in `14C`/`4He`. Fig. 5 should say screening-level public agreement.

### Fig. 4 Ghana Field Network

The earlier all-`evap` field result was a code-path issue: the field config enabled `mix`, but did not provide a mixing endmember, so the mix model had no candidate to evaluate. After adding a local dilute endmember:

- field edges: `208`
- transport models: `mix=121`, `evap=87`
- median chemistry `R2`: `0.991`

This improves internal consistency, but the Ghana network is still generated from coordinates/head proxies and has no independent process-truth graph.

## Required Claim Discipline

- Use `benchmark evidence`, `screening check`, `topology-only comparison`, `proxy`, and `field demonstration`.
- Avoid `proof`, `guarantee`, `validated process truth`, and `TracerLPM equivalence`.
- Keep PHREEQC as proxy evidence until a live PHREEQC backend is configured and rerun.
