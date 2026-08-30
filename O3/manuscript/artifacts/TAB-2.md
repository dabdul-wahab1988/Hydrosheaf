| Benchmark | Capture-type metric | Value [95% CI] | Correctness-type metric | Value [95% CI] | CIs overlap? |
|---|---|---|---|---|:---:|
| Topology | Recall (true MODPATH edges recovered), n=174 | 0.845 [0.784, 0.891] | Precision (inferred edges that are true), n=302 | 0.487 [0.431, 0.543] | No |
| Age / residence time | Reportability rate, n=1272 (fitted rows admitted by the reference-free guard) | 0.259 [0.235, 0.283] | Within-factor-2 agreement, n=675, held-out cross-validated, uncalibrated | 0.564 [0.527, 0.601] | No |
| Reaction | Recall, pooled across 16 reactions, n=681 (true active phases recovered) | 0.599 [0.562, 0.635] | Precision, macro mean across 240 scenarios (1 minus false-discovery rate) | 0.639 [0.598, 0.681] | Yes |

Notes. Intervals are 95% Wilson score intervals computed from the underlying
counts (`O3/analysis/python/derive_headline_metrics.py`); the reaction-layer
precision interval is inverted from the reaction benchmark's own published
false-discovery-rate 95% CI rather than recomputed from scratch. Topology:
independent Head-gradient scenario against the USGS Savage MODPATH
reference; the informed structural floor (receptor-set only, zero hydraulic
information) reaches F1 0.552 and is not independent evidence. Age: capture
and correctness are two different kinds of rate, not a literal
recall/precision pair -- reportability is the fraction of attempted fits
that produce an actionable claim at all (an abstention design, not a
defect: see Discussion); within-factor-2 is the accuracy of the claims that
are produced, scored on an independent held-out cross-validation split.
Reaction: recall is pooled (micro-averaged) from
`tableS6_reaction_confusion_matrices.csv` because no macro-averaged recall is
published; the macro-averaged precision (0.639, from
`table1_comparative_inverse_performance.csv`) is the reaction benchmark's own
published headline correctness value. A pooled-precision cross-check (0.586,
95% CI [0.549, 0.622]) is reported in
`manuscript/artifacts/data/evidence_discrepancies.csv` and differs from the
macro value only because of the averaging convention, not a computational
disagreement. **The reaction layer's capture and correctness intervals
overlap; only the topology and age layers show a gap that is clearly
resolved at the 95% level (see Results Section 4.2 and Discussion Section
5.1).**
