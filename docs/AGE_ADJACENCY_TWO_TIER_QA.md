# Two-tier age/direct-adjacency QA protocol

**Protocol identifier:** `age-adjacency-two-tier-v1`  
**Status:** implementation contract for new runs; it does not rewrite locked results.

## Scientific decision

The age module has two different estimands and they must be run as separate
tiers:

| Tier | Permitted information | Direct-adjacency decision | Scientific claim ceiling |
| --- | --- | --- | --- |
| `T1_temporal_order` | Endpoint age estimates and their uncertainty | Always `ABSTAIN` for directness; age can score only downstream-older compatibility | Temporal ordering/direction diagnostic |
| `T2_segment_transport` | T1 plus independently supplied edge-specific direct and indirect travel-time/RTD hypotheses | Score only when both hypotheses, uncertainty, and provenance are present; otherwise `ABSTAIN` | Bounded direct-versus-indirect discrimination under the declared transport model |

Endpoint-age ordering is a reachability-compatible signal. It cannot identify
the cover relation because `u -> w -> v` and `u -> v` can have the same endpoint
age ordering. T1 therefore contains the two-hop skip as an explicit
identifiability control, but must not convert temporal compatibility into a
direct-edge probability.

T2 compares the observed increment

\[
  \Delta a = a_v-a_u
\]

with two *pre-specified* hypotheses. With endpoint covariance and process
uncertainty,

\[
  \sigma_\Delta^2 = \sigma_u^2+\sigma_v^2-2\operatorname{Cov}(a_u,a_v)
                      +\sigma_{process}^2.
\]

The directness evidence is

\[
  \log BF_{direct:indirect} =
  \log p(\Delta a\mid T_{direct})-
  \log p(\Delta a\mid T_{indirect}).
\]

The Bayes factor is conditional on the supplied segment length, velocity/RTD,
dispersion, mixing, and error model. It is not a universal probability that a
well pair is directly connected.

## Package contract

Every run is a new immutable package JSON with the following top-level fields:

```text
schema, protocol_id, run_id, tier, claim_tier, metadata, units,
provenance, ages, transport_hypotheses, predictions, truth_artifact
```

Required canonical units are fixed and machine-checked:

```text
age_years              = years
age_sigma_years        = years
age_covariance_years2  = years^2
travel_time_years      = years
travel_sigma_years     = years
probability            = 1
```

An age record contains `case_id`, `node_id`, `age_years`,
`age_sigma_years`, `age_status`, and `source_id`. `age_status` is one of
`observed`, `interval`, `censored`, `missing`, or `invalid`. Interval and
censored records must retain finite lower and upper limits; they must not be
silently converted to point ages.

Each prediction contains `case_id`, `edge_id`, `u`, `v`,
`direction_probability`, `prediction_status`, `adjacency_status`, and
`identifiability_stratum`. A prediction status is exactly `scored`, `ABSTAIN`,
or `invalid`. A directness score is represented by
`direct_probability` and `log_bayes_factor_direct_vs_indirect` only when
`adjacency_status = scored`.

The allowed identifiability strata are:

- `endpoint_age_order_only`: T1-compatible ordering information;
- `edge_transport_comparison`: both independent T2 hypotheses are available;
- `transport_censored`: a model horizon, sink, boundary, or missing segment
  prevents a complete comparison;
- `age_observation_censored`: age information is interval/censored;
- `unidentifiable`: the available data do not support either directness
  hypothesis.

T2 transport records contain `case_id`, `edge_id`, `u`, `v`, `hypothesis`,
`travel_time_years`, `travel_sigma_years`, `evidence_source_id`,
`independent_of_endpoint_age`, `path_basis`, and `intermediate_nodes`.
The direct record has no intermediate nodes; the indirect record must identify
at least one intermediate node. A scored T2 row requires both records,
strictly positive hypothesis-specific uncertainty, unequal direct/indirect
hypotheses, and `independent_of_endpoint_age = true`.

## Truth sealing and leakage gate

Truth is a separate JSON sidecar, referenced by relative path and SHA-256 in
`truth_artifact`. It must declare `sealed = true`, `role = evaluation_only`,
the same `run_id`, and complete-graph denominators for direct edges and
reachable ordered pairs. The prediction table must not contain a relation
label, `is_true_edge`, graph closure, reachable-pair label, process truth, or
any renamed truth-bearing alias.

The package metadata must declare:

```text
truth_sealed = true
truth_access_mode = evaluation_only
candidate_set_frozen = true
prediction_input_columns != evaluation_only_columns
```

The QA audit rejects truth-bearing keys in predictions and model-input
declarations. Labels are joined only by an evaluation-side scorer after the
prediction is frozen. Development calibration and threshold selection must
precede locked-test scoring; the audit does not select either.

## Claim tiers

`claim_tier` is an evidence label, not a model setting:

- `direction_diagnostic`: T1 temporal compatibility only;
- `controlled_synthetic_component`: an independently generated, truth-sealed
  T2 result under a declared generator;
- `calibrated_model_reference`: a calibrated MODFLOW/MODPATH or tracer-model
  reference, not independent field truth;
- `field_transfer_screening`: observed field transfer without independent
  adjacency truth.

An independent synthetic package must use
`controlled_synthetic_component`, a sealed truth sidecar, and
`integrated_scoring_allowed = true`. Observed-field and conceptual packages
cannot set integrated scoring true without an independent reference.

## Aiken emulation restriction

The Aiken County release is a calibrated model reference. Its CFC ages and
MODPATH paths can support age--transport concordance and parser/component
checks, but the release does not supply independent well-to-well direct-edge
truth. An Aiken package must therefore declare:

```text
source_kind = calibrated_model_reference
claim_tier = calibrated_model_reference
aiken_emulation = true
integrated_scoring_allowed = false
truth_artifact.available = false
```

It must contain no scored directness prediction and must not be pooled into a
topology F1, direct-adjacency PR-AUC, or field-validation claim. Active or
right-censored MODPATH endpoints remain transport diagnostics; they are not
completed well-to-well truth paths. A forward run must also remain separate
from backward-to-recharge runs.

## Deterministic QA command

Run the dependency-light audit on a new package directory:

```powershell
.venv\Scripts\python.exe `
  M7\m7_nonuniqueness_benchmark\scripts\audit_age_adjacency_tiers.py `
  --package M7\m7_nonuniqueness_benchmark\results\RUN-AGE-ADJACENCY-YYYYMMDD-01\package.json `
  --output M7\m7_nonuniqueness_benchmark\results\RUN-AGE-ADJACENCY-YYYYMMDD-01\qa.json
```

The audit is a fail-closed structural/scientific contract check. It verifies
required fields, units, finite/range constraints, source hashes, safe truth
sidecar hashing, truth sealing, candidate uniqueness, endpoint coverage,
identifiability strata, T1/T2 rules, claim boundaries, and Aiken restrictions.
It does not certify a generator, transport model, causal interpretation, or
field validity; those remain scientific review gates.
