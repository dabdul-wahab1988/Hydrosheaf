# M6 validation and robustness diagnostics

This document maps each Q1 reviewer concern to the concrete code and locked result
that addresses it. All numbers are reproducible via `python scripts/run_m6_q1.py`
(seed 1234). New evidence lives in `run_m6_synthetic_validation.py`,
`run_m6_robustness_diagnostics.py` and
`synthetic_reaction_truth_model.py`.

## Concern 1 (MAJOR) — Circularity: the tier-collapse was built into the gate
**What was done.** (a) Sr and SiO2 were promoted from hand-imposed evidence gates to
*genuine reaction-basis species* (`synthetic_reaction_truth_model.py`):
silicate weathering releases
SiO2, carbonate/gypsum release Sr. Their diagnostic value is now linear algebra.
(b) The field tier-ablation was re-run with the evidence gate **ON vs OFF**, alongside
gate-independent structural diagnostics.

**Result (`m6_field_gate_structural.csv`, Extended Data Figure 3a).**
| tier | gated non-id. | 95% CI | ungated non-id. | rank | silicate coherence |
|---|---:|---:|---:|---:|---:|
| T0 majors | 59.4% | 51.9–67.5 | **0.0%** | 8 | 0.707 |
| T1 +iso | 51.9% | 45.0–60.0 | **0.0%** | 8 | 0.707 |
| T2 +F | 51.9% | 44.4–60.0 | **0.0%** | 9 | 0.707 |
| T3 +Sr/SiO2 | 0.6% | 0.0–1.9 | 0.0% | 11 | 0.500 |
| T4 full | 0.0% | — | 0.0% | 11 | 0.500 |

**Honest conclusion.** The field "collapse" is driven *entirely by the conservative
evidence gate* — the classifier alone never returns non-identifiable. The gate is now
reported as an explicit, tunable **conservatism prior**, not a discovery. Its *direction*
is independently validated by (i) structural rank rising 8→11 and silicate-vs-non-silicate
coherence falling 0.707→0.500 when Sr/SiO2 are added, and (ii) the synthetic ground-truth
recovery below. Its *magnitude* (52%) is a prior choice and is labelled as such.

## Concern 2 (MAJOR) — No independent validation
**What was done.** `run_m6_synthetic_validation.py` forward-generates 360 synthetic wells
from a dilute recharge endmember with **known reactions and extents** (plus SiO2/Sr
signatures, evaporation, 3% noise), then runs the identical M6 inference at each tier and
scores recovery against truth.

**Result (`m6_synthetic_recovery_by_tier.csv`, Extended Data Figure 2).**
- Dominant-process accuracy vs truth improves with Sr/SiO2 where theory predicts:
  **cation-exchange 0.45→0.66**, carbonate 1.00, silicate ~0.98; evaporite ~0.80
  (Cl-confounded, expected).
- **Exact-mineral F1 stays limited at every tier (~0.49–0.74)** — mechanism is not
  uniquely recoverable regardless of tracers, which is exactly why M6 reports at the
  class/family level.
- Structural rank 8→9→11; silicate coherence 0.707→0.500 (gate-independent).

## Concern 3 (MAJOR) — Classifier may not discriminate ("everything partial")
**Result (`m6_mrs_discrimination.csv`, `m6_synthetic_class_spread.csv`,
Extended Data Figure 2d).**
- Northern Ghana field MRS is genuinely narrow: mean 70.8, sd 3.5, IQR 4.5.
- But on synthetic wells with varying information the true identifiability class spans the
  full range — **identifiable 35% / partially 30% / non 23% / equivalence 12%**.

**Honest conclusion.** The framework *does* discriminate; Northern Ghana's uniform
"partially identifiable" outcome reflects a **hydrochemically homogeneous, good-QC,
dilute-silicate-weathering aquifer population**, not a failure of the classifier. This is
now stated explicitly rather than presented as a property of the method.

## Concern 4 (MODERATE) — Transport-correction and inference-unit sensitivity
**Result (`m6_transport_sensitivity.csv`, `m6_transport_agreement.csv`,
Extended Data Figure 3b).**
- **No correction is degenerate**: 93% of wells collapse to evapoconcentration (the pure
  concentration artifact), agreement with the Cl treatment only 5.6%.
- The two principled treatments — Cl-conservative within-well and per-region
  recharge-endmember (no independent aquifer-type classification exists for these
  boreholes) — agree on the dominant family for **64.4%** of wells and both give a
  silicate-weathering-dominated network (71% vs 70%). Mean MRS is stable (71–72) across all.

**Conclusion.** A transport correction is *necessary*, and the qualitative conclusion
(silicate weathering dominant) is robust to the choice between principled corrections.

## Concern 5 (MODERATE) — Talensi charge-balance
**Result (`m6_talensi_cbe_diagnosis.csv`).** Median CBE −29%; cations 2.66 vs anions 4.49
meq/L (2.12 meq anion excess); HCO3 is 73.5% of anions. Treating alkalinity as CaCO3
*worsens* the balance (−36%), ruling that out. The imbalance points to an unmeasured cation
pool (no trace-metal cations reported for this artisanal-mining catchment) or cation
under-recovery. **Talensi is retained as a screening-only failure-mode vignette, not a
quantitative result.**

## Concern 6 (MODERATE) — Novelty / altitude
Framing (not code): the paper is positioned as decision-support — the synthetic validation
identifies *which* extra measurement resolves *which* ambiguity (Sr/SiO2 → cation-exchange
vs silicate), converting the equifinality point into an actionable monitoring recommendation.

## Minor
- **Confidence intervals** added to the tier non-identifiable fractions (bootstrap, B=400;
  Extended Data Figure 3a ribbon).
- **Limitation-map wording** clarified: silicate is non-identifiable in external sets
  *because Sr/SiO2 are unmeasured*, not because silicate is intrinsically unresolvable.
- **S4 colour collision** noted; metrics and process legends are separated by panel.

## Net effect on the claims
The central, defensible message is now stronger and honest:
1. **Goodness-of-fit and support stability do not reveal information loss** (MRS/stability
   flat while structure and gated class change) — unchanged, and central.
2. **Sr/SiO2 have real, validated, but targeted value** (cation-exchange vs silicate
   disambiguation; structural rank), demonstrated against known truth — no longer circular.
3. **Exact mechanism is not field-identifiable at any tier**, so conservative class/family
   reporting is the correct policy — now validated, not asserted.
4. **The evidence gate is an explicit conservatism prior**, directionally validated and
   magnitude-labelled — the circularity is disclosed and bounded.
