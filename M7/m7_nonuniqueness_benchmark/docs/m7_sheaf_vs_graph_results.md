# M7.4 Results: Sheaf structure versus competence-matched weighted graphs

## Confirmatory status

`RUN-M7-SHEAF-VS-GRAPH-20260729-01` completed under the hash-locked protocol
with 32 development cases, 64 untouched locked-test cases (16 per scenario),
and 10,000 paired case-block bootstrap replicates. Development and test seeds
were disjoint; all primary metrics were finite; the independent generator
imported no HydroSheaf code; and the identity-limit maximum raw residual
difference between the graph Laplacian and affine sheaf was exactly 0.0. The
execution/equivalence gate therefore passed.

## Overall performance

| Model | Mean PR-AUC | Mean Brier | Mean log loss | Mean selected F1 |
|---|---:|---:|---:|---:|
| Edge-local weighted graph | 0.6665 | 0.1772 | 0.5189 | 0.6219 |
| Identity graph Laplacian | 0.5908 | 0.1970 | 0.5778 | 0.5647 |
| Affine sheaf | 0.6762 | 0.1777 | 0.5306 | 0.6235 |
| Permuted-map sheaf | 0.5853 | 0.1992 | 0.5852 | 0.5506 |

Relative to the identity graph Laplacian, the affine sheaf increased PR-AUC by
0.0854 (95% CI 0.0666 to 0.1050), reduced Brier score by 0.0193 (difference
-0.0193, CI -0.0235 to -0.0152), reduced log loss by 0.0472 (difference
-0.0472, CI -0.0573 to -0.0372), and increased selected-edge F1 by 0.0588
(CI 0.0418 to 0.0759). Native maps also outperformed the within-case
permuted-map control in PR-AUC by 0.0909 (CI 0.0705 to 0.1117), showing that
the gain depended on the edge-specific restrictions rather than the marginal
distribution of gains and offsets.

The stronger edge-local weighted graph was not beaten overall. The affine-
sheaf PR-AUC difference was +0.0097 (CI -0.0054 to +0.0248); Brier-score
difference was +0.00053 (CI -0.00330 to +0.00443); selected-F1 difference was
+0.00161 (CI -0.01565 to +0.01867); and log loss was 0.0117 higher (CI 0.00083
to 0.0230). The full predeclared incremental-superiority gate therefore failed.

## Scenario interpretation

In the identity-limit stratum, the identity graph Laplacian and affine sheaf
were identical in every reported metric, as required. In heterogeneous affine
cases, the sheaf materially outperformed the identity Laplacian but the
edge-local weighted graph remained strongest, showing that non-identity maps
matter while global solving is not automatically superior to informative local
compatibility. In noisy/missing cases the sheaf exceeded both graph comparators
in mean PR-AUC, but that stratum was secondary and is not promoted to a general
claim.

The clearest uniquely global contribution appeared in the incompatible-cycle
stratum. There, sheaf-minus-edge-local-graph PR-AUC was +0.0483 (95% CI 0.0258
to 0.0711), Brier difference was -0.00552 (CI -0.00895 to -0.00231), and
planted-conflict localisation PR-AUC improved by 0.0689 (CI 0.0467 to 0.0915).
Against the identity Laplacian, conflict-localisation PR-AUC improved by 0.1098
(CI 0.0917 to 0.1278).

## Claim decision

The predeclared full incremental-superiority claim is **not allowed**, because
the sheaf did not improve both PR-AUC and Brier score over the strong edge-local
weighted graph across all locked-test scenarios. The predeclared conditional
conflict-localisation claim **is allowed**:

> In this prospectively locked controlled-synthetic benchmark, edge-specific
> affine restrictions supplied value beyond an ordinary identity-coupled graph
> by representing non-identity relations and localising planted global
> inconsistencies. This value was conditional: the affine sheaf was not
> superior to a strong edge-local weighted graph overall and reduced exactly to
> the ordinary graph formulation in the identity limit.

The result is a model-conditioned scalar, static graph capability test. It is
not field validation, does not establish general HydroSheaf superiority, and
does not evaluate temporal, spatial three-dimensional, or vadose-zone
capability.

