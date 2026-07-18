# M7 — Coupled Integration Benchmark (thesis keystone)

Demonstrates the thesis's central claim: **integrating residence-time (age), graph
topology (connectivity) and inverse hydrogeochemistry (reactions), and testing their
mutual consistency, reduces interpretive non-uniqueness** in data-limited aquifers.

This is the result the field data cannot provide — no real aquifer supplies the true
flow network *and* true ages *and* true reactions at once — so it is demonstrated on a
controlled synthetic twin with known joint truth, driving the REAL Hydrosheaf modules.

See `docs/m7_integration_defensibility.md` for the full design and defensibility argument.

## M7.1 blind replicated stress test

The original locked M7 is retained for historical reproducibility. The stronger
M7.1 extension removes truth-informed candidate construction and tests the real
pipeline on independent aquifer realizations:

```bash
python scripts/run_m7_1_blind_benchmark.py
```

M7.1 separates truth generation from inference, freezes a regularized fusion
model and decision thresholds on 20 development aquifers, evaluates them once on
100 disjoint test aquifers, and runs an 8-aquifer PHREEQC/Bayesian-aging/topology-
posterior audit. It includes a truth-poisoning leakage test, strict held-out ions,
aquifer-level bootstrap intervals, calibration and graph-structure metrics, and a
machine-readable provenance manifest. See `docs/m7_1_blind_protocol.md` and
`docs/m7_1_blind_results.md`.

## What it drives (real framework, synthetic data)
- `hydrosheaf.nuclear.age_coherence.audit_graph_age_coherence` — age↔topology consistency audit
- `hydrosheaf.nuclear.lpm.convolve_input` + `input_history` — tritium forward model + single-node ages
- M5/M6 sparse inverse-reaction pillar (`m6_reactions.py`); chemistry falsifier via mass/Cl conservation
- (network-enhanced Bayesian dating is a dating-method question benchmarked in M3, not M7)

## Four integration tests
| Test | Coupling | Result file |
|---|---|---|
| T1 | age ↔ topology ( tritium age-coherence coupling) | `m7_age_gain.csv`, `m7_age_recovery.csv` |
| T2 | chemistry ↔ topology (falsifies impossible edges) | `m7_edge_classification.csv` |
| T3 | age ↔ chemistry (rejects age reversals) | `m7_edge_classification.csv` |
| T4 | integration gain (joint vs single streams) | `m7_integration_gain.csv`, `m7_trap_rejection.csv` |

## Reproduce
```bash
python scripts/run_m7_integration_benchmark.py          # deterministic, seed 1234
Rscript r_figures/plot_m7_integration_figures.R
```

## Optional operational synthetic-twin extension

The locked benchmark above remains the known-joint-truth integration test. A
separate extension demonstrates the update–forecast–verify loop required by an
operational digital twin:

```bash
python scripts/run_m7_digital_twin_extension.py
```

It writes only to `results/digital_twin/`, `figures/digital_twin/`, and
`docs/m7_digital_twin_results.md`; baseline M7 outputs are hashed before and after
the run. The extension uses transient graph states, sparse asynchronous monitoring,
sequential ensemble Kalman updating, 1/3/6-month prequential forecasts, uncertainty
calibration, and wrong-topology/shuffled-observation/sensor-dropout controls.

The scientifically accurate term is **operational synthetic-aquifer twin**, not a
validated field digital twin. See `docs/m7_digital_twin_protocol.md` for the locked
evaluation protocol, anti-inverse-crime design, permitted claims, and field
requirements.

## Honesty boundary
A controlled capability/mechanism demonstration, not field validation. Chemistry acts as
a falsifier (rejects impossible edges) not a confirmer (equifinality); the integration
value is the streams' complementary blind spots — only the joint test rejects every trap
type. The Ghana field application (M6) remains the realistic, chemistry-dominant transfer.

The replicated M7.1 result supersedes the locked single-realization M7 for
performance claims. In M7.1, naive equal-weight integration is inferior to the
hydraulic-spatial baseline; only development-trained logistic fusion improves edge
classification, and it does not improve every metric. The fitted age coefficient
is negative, neither topology nor Bayesian-age posteriors pass their convergence
gates, and the PHREEQC constraints are non-binding in the heavy audit. This is
synthetic capability evidence and a model-development diagnostic, not field
validation or a validated digital twin.
