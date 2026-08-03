| Benchmark | Scenario | Row type | Independent | External reference | Headline metric | Value |
|---|---|---|:---:|---|---|---:|
| Topology | Spatial only | independent | yes | MODPATH (Savage archive) | F1 | 0.000 |
| Topology | Sink-aware baseline | informed control | no | MODPATH (Savage archive) | F1 | 0.552 |
| Topology | Head gradient | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Hodge pruned (diagnostic) | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Proj. gradient (diagnostic) | independent | yes | MODPATH (Savage archive) | F1 | 0.618 |
| Topology | Head depth | independent | yes | MODPATH (Savage archive) | F1 | 0.616 |
| Topology | Hydrostratigraphic | independent | yes | MODPATH (Savage archive) | F1 | 0.614 |
| Topology | Negative: random | negative control | no | MODPATH (Savage archive) | F1 | 0.006 |
| Topology | Negative: wrong direction | negative control | no | MODPATH (Savage archive) | F1 | 0.000 |
| Topology | Negative: misrouted sink | negative control | no | MODPATH (Savage archive) | F1 | 0.138 |
| Topology | Negative: shortcut (inapplicable) | negative control | no | MODPATH (Savage archive) | F1 | n/a |
| Topology | Sparse node | sensitivity diagnostic | no | MODPATH (Savage archive) | mean F1, 20 trials | 0.445 |
| Topology | Prior: override / merge / only | prior-informed | no | MODPATH (Savage archive) | edges in/out | 302→302/329/174 |
| Age / residence time | TracerLPM strict parity | independent | yes | USGS TracerLPM Table 4 | within-factor-2 | 0.875 (n=329/1272) |
| Age / residence time | TracerLPM age fractions | non-independent | no | USGS reported age fractions | within-factor-2 | 0.917 (n=289/1272) |
| Age / residence time | Hydrosheaf model selection | independent | yes | USGS national age release | within-factor-2 | 0.673 (n=309/1272) |
| Age / residence time | Old-water He-4 branch | independent | yes | USGS national age release | within-factor-2 | 0.697 (n=66/1272) |
| Age / residence time | Graph regularisation, randomised control | negative control | yes | Randomised candidate graph | Δ log10 RMSE vs single-node | +0.655 |
| Reaction | Bounded LS | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.435 |
| Reaction | Lasso | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.554 |
| Reaction | Elastic Net | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.551 |
| Reaction | Thermo Elastic Net (legacy primary) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.551 |
| Reaction | Hydrosheaf Guarded (primary) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.563 |
| Reaction | Hydrosheaf-Core (evidence-gated) | independent | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.561 |
| Reaction | Conventional PHREEQC inverse (strict+relaxed fallback) | conventional baseline | yes | 240-scenario PHREEQC benchmark | Phase F1 | 0.510 |
