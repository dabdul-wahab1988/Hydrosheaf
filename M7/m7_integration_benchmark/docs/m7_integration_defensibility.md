# M7 Coupled Integration Benchmark — Design & Defensibility

The thesis's central claim is that integrating three interdependent evidence
streams — residence-time (age), graph topology (connectivity) and inverse
hydrogeochemistry (reactions) — and testing their **mutual consistency** reduces
interpretive non-uniqueness in data-limited aquifers. Every other tier (M3–M6)
validates one pillar in isolation and the Ghana field application runs
chemistry-only (no tracers, no hydraulic heads). M7 supplies the missing
keystone: a demonstration that the three streams *jointly constrain each other*.

## Why a controlled synthetic twin (and why that is defensible)
No real aquifer provides the true flow network **and** true ages **and** true
reactions simultaneously, so the *reduction* in non-uniqueness the thesis claims
cannot be measured on field data — only demonstrated against known joint truth.
A controlled twin is therefore not a fallback; it is the only instrument that can
measure the causal claim "integration → fewer false interpretations." This is
consistent with the thesis's declared computational/methodological scope
(Chapter 1 §1.5) and with the synthetic and MODPATH reference-truth tiers already
used in M3–M5.

Four safeguards keep it defensible:
1. **Anti-inverse-crime.** The generator differs from the inversion: evaporative
   concentration, 3–5% analytical noise, and an *unmodelled mixing perturbation*
   on Na/Cl/Mg that no dictionary reaction can represent. The pipeline never
   receives back its own assumptions.
2. **Negative controls (trap edges).** Planted false connections, each designed
   to defeat exactly one single stream:
   - *Trap A* — geometrically down-gradient and chemically plausible, but an age
     reversal (edge into a young, mineralised local-recharge node). Only the age
     stream rejects it.
   - *Trap B* — geometrically down-gradient and age-coherent, but chemically
     impossible (edge into an old, dilute fossil water requires net precipitation).
     Only the chemistry stream rejects it.
   - *Spurious* — distant / up-gradient random pairs (geometry rejects).
3. **Drives framework components.** The locked benchmark uses the single-node
   LPM and age-coherence audit plus the M5/M6 reaction pillar; it does not execute
   Bayesian network aging or the full blind topology selector. The separate M7.1
   stress test executes candidate inference, sheaf refinement, inverse chemistry,
   constrained topology posterior, PHREEQC, and Bayesian network aging under a
   truth–inference firewall.
4. **Capability, not field truth.** M7 demonstrates the integration *mechanism*;
   it is explicitly not a claim about any real aquifer. The Ghana field transfer
   (M6) remains the realistic, chemistry-dominant application.

## The four integration tests
- **T1 age ↔ topology.** Tritium ages (accurate where the tracer is informative:
  within-factor-2 = 0.89 for < 70 yr water; modern/pre-modern accuracy = 1.00) power
  the age-coherence audit, which separates true edges from age-reversal traps. The
  network-enhanced Bayesian dating *method* is benchmarked separately in M3; M7 is an
  integration benchmark, not a dating-method benchmark.
- **T2 chemistry ↔ topology.** Chemistry (held-out-ion falsification + net-dissolution
  direction) rejects chemically-impossible edges (Trap B) that geometry and age accept.
- **T3 age ↔ chemistry.** Age coherence rejects age-reversal edges (Trap A) that
  geometry and chemistry accept.
- **T4 integration gain.** The joint classifier (geometry AND age AND chemistry)
  vs each single stream: precision, recall, and false-connection (trap) acceptance.

## Honest framing of the result
Chemistry is a **falsifier, not a confirmer**: it cleanly rejects physically
impossible edges but, because the reaction dictionary is underdetermined
(equifinality), it cannot uniquely confirm a plausible edge. Age coherence is a
strong, physical topology discriminator. The integration value is that the streams
have **complementary blind spots**: each trap type is rejected by exactly one
stream, so only the joint test rejects all of them. The gain is a large reduction
in false connections (lower non-uniqueness) at a conservative recall cost — the
same "reduce overinterpretation" trade-off that characterises the whole framework.

Locked headline numbers are in `docs/m7_results.md` / `results/`.

These safeguards describe the locked, constructed-trap demonstration. For
performance claims, the replicated blind M7.1 protocol and its adverse controls
supersede this single realization. M7.1 shows that naive equal-weight coupling is
inferior to its hydraulic-spatial baseline; only development-trained logistic fusion
improves selected held-out metrics, and important chemistry, age, and calibration
limitations remain.
