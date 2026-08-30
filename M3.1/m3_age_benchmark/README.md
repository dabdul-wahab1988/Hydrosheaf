# M3 Age Benchmark

This folder contains the Hydrosheaf M3 nuclear-tracer validation benchmark. It
is designed to test claims about graph-regularized age inference rather than
assume that network information always improves age estimates.

Run from the repository root:

```powershell
python M3/m3_age_benchmark/scripts/run_m3_age_benchmark.py
```

The benchmark compares:

- single-node age estimates;
- weak, medium, and strong graph-regularized estimates;
- randomized negative-control graphs;
- scenarios where graph regularization improves inference;
- scenarios where graph regularization degrades inference.

Outputs are written to `results/`, `tables/`, `figures/`, and `docs/`.

## Identified-TTD evidence path

M3 also contains a development-stage, set-valued age-evidence path. It asks
which age fractions are supported by the tracer-response matrix instead of
forcing every sample into one best scalar age or LPM family. It provides:

- sharp lower and upper bounds with witness TTDs for declared age functionals;
- explicit `IDENTIFIED`, `PARTIALLY_IDENTIFIED`, and `ABSTAIN` states;
- one-sided support for censored observations and explicit exclusion of
  contaminated-mixture rows that lack linear interval semantics;
- leakage-guarded prediction of held-out tracer concentrations;
- graph compatibility/falsification without graph-based tightening; and
- tracer/well/time information-gain ranking only for ensembles with declared
  probability semantics.

The frozen development configuration is
`configs/identified_ttd_protocol.yaml`. A small run is:

```powershell
python M3/m3_age_benchmark/scripts/run_m3_identified_ttd_benchmark.py `
  --max-rows 10 `
  --output M3/m3_age_benchmark/results/m3_identified_ttd_development.csv
```

The runner refuses to overwrite an existing result by default and writes a
hash-bearing manifest beside the CSV. Reported USGS ages and age fractions are
excluded from fitting and threshold selection. This implementation does not,
by itself, validate true TTD recovery, graph benefit, field flow paths, or an
optimal field sampling policy. Those claims require the prospective protocol
and confirmatory evidence described in
`docs/m3_identified_ttd_oed_implementation_plan_20260730.md`.
