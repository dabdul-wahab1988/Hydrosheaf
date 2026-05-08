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
