# M4 Topology Benchmark

This folder contains a focused topology-validation benchmark for Hydrosheaf.
The benchmark asks:

> Can Hydrosheaf's reduced-order graph construction reproduce reference
> advective connectivity from MODPATH particle tracking under controlled
> benchmark conditions?

The benchmark keeps two modes separate:

- independent graph inference compared against MODPATH reference connectivity;
- MODPATH-informed graph priors used inside Hydrosheaf.

Run from the repository root:

```powershell
python M4/m4_topology_benchmark/scripts/run_m4_topology_benchmark.py
```

Outputs are written to `results/`, `tables/`, and `docs/`.
