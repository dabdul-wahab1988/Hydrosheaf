# M2 Analysis Consistency Review

Updated: 2026-05-19

## Reviewer Gate

The M2 figure package is usable only with conservative captions. The synthetic benchmark, MODPATH topology check, M3 public-age parity benchmark, PHREEQC proxy, and Ghana field demonstration are different evidence tiers and must not be collapsed into one claim of full validation.

## Corrections Applied

- Fig. 5 and Supplementary Fig. S1 now prefer `M3/m3_age_benchmark/results/m3_tracerlpm_parity_agefractions_full.csv` for the public USGS age panel, with older M3/M2 outputs kept only as fallbacks.
- M3 was rerun at the canonical full setting: `--full --age-steps 90 --scenario tracerlpm_parity_agefractions`.
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

The refreshed M3 age-fraction-constrained parity benchmark is the strongest independent public-age benchmark currently available:

- rows: `1272`; finite log-error rows: `1249`
- median absolute log10 error: `0.167`
- log10 RMSE: `0.740`
- log10 R2: `0.681`
- within factor 2: `0.614`
- within factor 10: `0.859`

This result uses USGS reported age fractions as TTD-shape constraints, so Fig. 5 should describe parity to USGS reported model outputs, not direct true-age validation or TracerLPM equivalence.

### Fig. 4 Ghana Field Network

The earlier all-`evap` field result was a code-path issue: the field config enabled `mix`, but did not provide a mixing endmember, so the mix model had no candidate to evaluate. After adding a local dilute endmember:

- field edges: `208`
- transport models: `mix=121`, `evap=87`
- median chemistry `R2`: `0.991`

This improves internal consistency, but the Ghana network is still generated from coordinates/head proxies and has no independent process-truth graph.

## Required Claim Discipline

- Use `benchmark evidence`, `reported-model parity`, `topology-only comparison`, `proxy`, and `field demonstration`.
- Avoid `proof`, `guarantee`, `validated process truth`, and `TracerLPM equivalence`.
- Keep PHREEQC as proxy evidence until a live PHREEQC backend is configured and rerun.
