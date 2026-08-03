# Supplementary Figures and Tables

**Conditional evidence integration and the incremental contribution of sheaf structure in controlled-synthetic groundwater benchmarks**

The complete authoritative supplementary tables are the version-controlled CSV files in `tables/publication/`. The compact views below preserve the claim-bearing comparisons while keeping the Word and PDF tables legible; omitted columns and per-case or per-edge rows remain available without loss in the cited CSV. Tables S1-S6 derive from locked process-based integration outputs, Table S7 from the prospectively locked representation benchmark, Table S8 from the strict public-pipeline acceptance run, Table S9 from the prospectively locked estimator diagnostic, and Tables S10-S13 from the post-review family-wise, precision and selection audit. Table S14 is the auxiliary M7.6 controlled-synthetic M3-mechanism diagnostic.

## Figure S1

![](figures/supporting_validation/figure_s1_model_domain_map.png)

**Figure S1.** Locked synthetic model domain (representative realization 4101). The MODFLOW/MODPATH truth network, observation nodes and hydraulic heads are shown in synthetic model-space coordinates. This is not a geographic map.

## Table S1

**Table S1.** Evidence-panel performance in every locked process-based integration condition. Compact view of `tables/publication/tableS1_all_evidence_conditions.csv`.

| Condition          | Panel   |   PR-AUC |   Brier |   Log loss |   Mean entropy |
|:-------------------|:--------|---------:|--------:|-----------:|---------------:|
| native             | H       |   0.1765 |  0.1076 |     0.3467 |         0.5029 |
| native             | A       |   0.1112 |  0.1087 |     0.3558 |         0.524  |
| native             | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| native             | HA      |   0.1114 |  0.1076 |     0.3467 |         0.5042 |
| native             | HC      |   0.4846 |  0.088  |     0.2676 |         0.4221 |
| native             | AC      |   0.471  |  0.089  |     0.2736 |         0.4361 |
| native             | HAC     |   0.4805 |  0.088  |     0.2677 |         0.4215 |
| age_permuted       | H       |   0.1765 |  0.1076 |     0.3467 |         0.5029 |
| age_permuted       | A       |   0.1397 |  0.1179 |     0.4539 |         0.524  |
| age_permuted       | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| age_permuted       | HA      |   0.1825 |  0.1096 |     0.3588 |         0.4726 |
| age_permuted       | HC      |   0.4846 |  0.088  |     0.2676 |         0.4221 |
| age_permuted       | AC      |   0.3521 |  0.1016 |     0.3588 |         0.4236 |
| age_permuted       | HAC     |   0.4165 |  0.0914 |     0.2781 |         0.4014 |
| hydraulic_permuted | H       |   0.1302 |  0.1191 |     0.467  |         0.5029 |
| hydraulic_permuted | A       |   0.1112 |  0.1087 |     0.3558 |         0.524  |
| hydraulic_permuted | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| hydraulic_permuted | HA      |   0.1173 |  0.1155 |     0.4358 |         0.4713 |
| hydraulic_permuted | HC      |   0.3803 |  0.1032 |     0.3717 |         0.4054 |
| hydraulic_permuted | AC      |   0.471  |  0.089  |     0.2736 |         0.4361 |
| hydraulic_permuted | HAC     |   0.3971 |  0.1001 |     0.3482 |         0.3879 |
| joint_misspecified | H       |   0.1302 |  0.1191 |     0.467  |         0.5029 |
| joint_misspecified | A       |   0.1397 |  0.1179 |     0.4539 |         0.524  |
| joint_misspecified | C       |   0.4587 |  0.0944 |     0.2951 |         0.4551 |
| joint_misspecified | HA      |   0.1371 |  0.1201 |     0.463  |         0.4729 |
| joint_misspecified | HC      |   0.3803 |  0.1032 |     0.3717 |         0.4054 |
| joint_misspecified | AC      |   0.3521 |  0.1016 |     0.3588 |         0.4236 |
| joint_misspecified | HAC     |   0.3305 |  0.105  |     0.3682 |         0.3845 |

## Table S2

**Table S2.** Predeclared evidence contrasts for discrimination and calibration (10,000 case-block bootstrap replicates). Complete five-metric record: `tables/publication/tableS2_case_block_contrasts.csv`.

| Contrast                      | Metric   |   Difference |   CI low |   CI high |   Cases |
|:------------------------------|:---------|-------------:|---------:|----------:|--------:|
| native incremental age        | PR-AUC   |      -0.006  |  -0.0122 |   -0.0011 |      12 |
| native incremental age        | Brier    |       0.0001 |  -0      |    0.0001 |      12 |
| native incremental age        | Log loss |       0      |  -0.0001 |    0.0002 |      12 |
| native incremental chemistry  | PR-AUC   |       0.4471 |   0.3575 |    0.5401 |      12 |
| native incremental chemistry  | Brier    |      -0.0196 |  -0.0213 |   -0.0176 |      12 |
| native incremental chemistry  | Log loss |      -0.0791 |  -0.085  |   -0.072  |      12 |
| native incremental hydraulics | PR-AUC   |       0.0091 |   0.001  |    0.0198 |      12 |
| native incremental hydraulics | Brier    |      -0.001  |  -0.0012 |   -0.0008 |      12 |
| native incremental hydraulics | Log loss |      -0.0059 |  -0.007  |   -0.005  |      12 |
| permuted age increment        | PR-AUC   |      -0.0754 |  -0.1353 |   -0.0148 |      12 |
| permuted age increment        | Brier    |       0.0034 |   0.0012 |    0.0055 |      12 |
| permuted age increment        | Log loss |       0.0105 |   0.003  |    0.0176 |      12 |
| permuted hydraulic increment  | PR-AUC   |      -0.0686 |  -0.112  |   -0.0271 |      12 |
| permuted hydraulic increment  | Brier    |       0.011  |   0.0055 |    0.0164 |      12 |
| permuted hydraulic increment  | Log loss |       0.0745 |   0.0404 |    0.1091 |      12 |
| joint misspecification        | PR-AUC   |      -0.139  |  -0.2041 |   -0.0742 |      12 |
| joint misspecification        | Brier    |       0.0106 |   0.0055 |    0.0163 |      12 |
| joint misspecification        | Log loss |       0.073  |   0.0401 |    0.1073 |      12 |

## Table S3

**Table S3.** Topology-to-age sensitivity aggregated over twelve locked cases. Values are means except median ESS; the complete CSV also reports bias, entropy, order violations and importance weights. Complete 96-row case audit: `tables/publication/tableS3_topology_age_sensitivity.csv`.

| Regime      | Topology   |   n |   MAE (yr) |   Coverage |   Width (yr) |       ESS |   Stable |
|:------------|:-----------|----:|-----------:|-----------:|-------------:|----------:|---------:|
| Informative | None       |  12 |     2.7644 |     0.9306 |      12.5035 | 50000     |   1      |
| Informative | Partial    |  12 |     2.7399 |     0.9236 |      12.3584 | 48866.4   |   1      |
| Informative | Complete   |  12 |     2.7025 |     0.9444 |      12.2511 | 48301.1   |   1      |
| Informative | Reversed   |  12 |     2.9846 |     0.9167 |      12.5469 |  2558     |   1      |
| 3H only     | None       |  12 |     4.7501 |     0.9583 |      25.7014 | 50000     |   1      |
| 3H only     | Partial    |  12 |     4.6927 |     0.9583 |      25.2452 | 47183.5   |   1      |
| 3H only     | Complete   |  12 |     4.586  |     0.9722 |      24.7896 | 45117.1   |   1      |
| 3H only     | Reversed   |  12 |     4.6101 |     0.9167 |      21.396  |   227.926 |   0.3333 |

## Table S4

**Table S4.** Reaction-family non-uniqueness aggregated by benchmark tier and true process (C, core; E, enhanced). Complete 216-row edgewise audit: `tables/publication/tableS4_reaction_edge_nonuniqueness.csv`.

| Tier   | Process            |   n |   Accuracy |   True prob. |   Entropy |   Eff. families |
|:-------|:-------------------|----:|-----------:|-------------:|----------:|----------------:|
| C      | Carbonate weather. |  24 |        0   |       0      |    0.1155 |          1.1378 |
| E      | Carbonate weather. |  24 |        0   |       0      |    0.0902 |          1.1103 |
| C      | Carbonate precip.  |  12 |        0   |       0      |    0.0338 |          1.0417 |
| E      | Carbonate precip.  |  12 |        0   |       0      |    0      |          1      |
| C      | Silicate weather.  |  24 |        1   |       0.8639 |    0.349  |          1.4663 |
| E      | Silicate weather.  |  24 |        1   |       0.9505 |    0.1613 |          1.1987 |
| C      | Denitrification    |  24 |        1   |       0.9941 |    0.0191 |          1.0226 |
| E      | Denitrification    |  24 |        1   |       0.9694 |    0.0623 |          1.0805 |
| C      | Sulfate reduction  |  12 |        1   |       0.9974 |    0.0134 |          1.014  |
| E      | Sulfate reduction  |  12 |        1   |       0.9857 |    0.0577 |          1.0636 |
| C      | Iron reduction     |  12 |        0   |       0      |    0.3413 |          1.4493 |
| E      | Iron reduction     |  12 |        0.5 |       0.4688 |    0.5256 |          1.7147 |

## Table S5

**Table S5.** Conflict diagnostics under native and adverse evidence conditions. The complete CSV additionally reports error rates within flagged and concordant subsets: `tables/publication/tableS5_conflict_diagnostics.csv`.

| Condition      |   N |   Flagged |   Span |   Error |   Overconf. |
|:---------------|----:|----------:|-------:|--------:|------------:|
| Native         | 827 |         0 | 0.092  |  0.1306 |      0      |
| Age perm.      | 827 |         0 | 0.1193 |  0.1306 |      0.0036 |
| Hyd. perm.     | 827 |         0 | 0.1205 |  0.1306 |      0.0351 |
| Joint misspec. | 827 |         0 | 0.123  |  0.1306 |      0.0387 |

## Table S6

**Table S6.** Exact paired-permutation tests with Benjamini-Hochberg correction for the predeclared process-based integration contrast family. Complete record: `tables/publication/tableS6_multiplicity_correction.csv`.

| Contrast       | Metric   |   n |   Difference |   Exact p |   BH p |
|:---------------|:---------|----:|-------------:|----------:|-------:|
| Native +age    | PR-AUC   |  12 |      -0.006  |    0.0625 | 0.0703 |
| Native +age    | Brier    |  12 |       0.0001 |    0.1875 | 0.1957 |
| Native +age    | Log loss |  12 |       0      |    0.8491 | 0.8491 |
| Native +age    | Entropy  |  12 |      -0.0006 |    0.002  | 0.0043 |
| Native +chem.  | PR-AUC   |  12 |       0.4471 |    0.0005 | 0.0012 |
| Native +chem.  | Brier    |  12 |      -0.0196 |    0.0005 | 0.0012 |
| Native +chem.  | Log loss |  12 |      -0.0791 |    0.0005 | 0.0012 |
| Native +chem.  | Entropy  |  12 |      -0.0827 |    0.0005 | 0.0012 |
| Native +hyd.   | PR-AUC   |  12 |       0.0091 |    0.0645 | 0.0703 |
| Native +hyd.   | Brier    |  12 |      -0.001  |    0.0005 | 0.0012 |
| Native +hyd.   | Log loss |  12 |      -0.0059 |    0.0005 | 0.0012 |
| Native +hyd.   | Entropy  |  12 |      -0.0146 |    0.0005 | 0.0012 |
| Permuted +age  | PR-AUC   |  12 |      -0.0754 |    0.041  | 0.0492 |
| Permuted +age  | Brier    |  12 |       0.0034 |    0.0181 | 0.0241 |
| Permuted +age  | Log loss |  12 |       0.0105 |    0.0264 | 0.0333 |
| Permuted +age  | Entropy  |  12 |      -0.0207 |    0.0005 | 0.0012 |
| Permuted +hyd. | PR-AUC   |  12 |      -0.0686 |    0.0151 | 0.0214 |
| Permuted +hyd. | Brier    |  12 |       0.011  |    0.0049 | 0.0073 |
| Permuted +hyd. | Log loss |  12 |       0.0745 |    0.0029 | 0.005  |
| Permuted +hyd. | Entropy  |  12 |      -0.0482 |    0.0005 | 0.0012 |
| Joint misspec. | PR-AUC   |  12 |      -0.139  |    0.0044 | 0.007  |
| Joint misspec. | Brier    |  12 |       0.0106 |    0.0024 | 0.0045 |
| Joint misspec. | Log loss |  12 |       0.073  |    0.0024 | 0.0045 |
| Joint misspec. | Entropy  |  12 |      -0.0706 |    0.0005 | 0.0012 |

## Table S7

**Table S7.** Claim-bearing representation-benchmark sheaf-versus-graph contrasts, including the identity limit and incompatible-cycle diagnostics (10,000 case-block bootstrap replicates). Complete 120-row contrast matrix: `tables/publication/tableS7_sheaf_vs_graph_contrasts.csv`.

| Scenario       | Comparator   | Metric          |   Cases |   Difference |   CI low |   CI high |
|:---------------|:-------------|:----------------|--------:|-------------:|---------:|----------:|
| All            | Edge-local   | PR-AUC          |      64 |       0.0097 |  -0.0054 |    0.0248 |
| All            | Edge-local   | Brier           |      64 |       0.0005 |  -0.0033 |    0.0044 |
| All            | Edge-local   | Log loss        |      64 |       0.0117 |   0.0008 |    0.0232 |
| All            | Edge-local   | Selected F1     |      64 |       0.0016 |  -0.0156 |    0.0183 |
| All            | Edge-local   | Conflict PR-AUC |      16 |       0.0689 |   0.0466 |    0.0914 |
| All            | Identity     | PR-AUC          |      64 |       0.0854 |   0.0666 |    0.105  |
| All            | Identity     | Brier           |      64 |      -0.0193 |  -0.0235 |   -0.0152 |
| All            | Identity     | Log loss        |      64 |      -0.0472 |  -0.0573 |   -0.0372 |
| All            | Identity     | Selected F1     |      64 |       0.0588 |   0.0418 |    0.0758 |
| All            | Identity     | Conflict PR-AUC |      16 |       0.1098 |   0.0912 |    0.1277 |
| All            | Permuted-map | PR-AUC          |      64 |       0.0909 |   0.0705 |    0.1117 |
| All            | Permuted-map | Brier           |      64 |      -0.0215 |  -0.0263 |   -0.0169 |
| All            | Permuted-map | Log loss        |      64 |      -0.0546 |  -0.066  |   -0.0433 |
| All            | Permuted-map | Selected F1     |      64 |       0.0729 |   0.0491 |    0.0964 |
| All            | Permuted-map | Conflict PR-AUC |      16 |       0.052  |   0.0323 |    0.07   |
| Identity limit | Identity     | PR-AUC          |      16 |       0      |   0      |    0      |
| Identity limit | Identity     | Brier           |      16 |       0      |   0      |    0      |
| Identity limit | Identity     | Log loss        |      16 |       0      |   0      |    0      |
| Identity limit | Identity     | Selected F1     |      16 |       0      |   0      |    0      |
| Incompatible   | Edge-local   | PR-AUC          |      16 |       0.0483 |   0.0258 |    0.0712 |
| Incompatible   | Edge-local   | Conflict PR-AUC |      16 |       0.0689 |   0.0467 |    0.0912 |

## Table S8

**Table S8.** Strict public-pipeline system acceptance. The execution criterion passed, whereas a general full-sheaf incremental-performance claim did not.

### Condition summary

| Condition      |   Cases |   Recall |   PR-AUC |   Brier |   Log loss |   Selected F1 |
|:---------------|--------:|---------:|---------:|--------:|-----------:|--------------:|
| age permuted   |       6 |   0.9815 |   0.3211 |  0.2378 |     0.6953 |        0.4222 |
| full sheaf     |       6 |   0.9815 |   0.3075 |  0.2171 |     0.6243 |        0.4222 |
| hydraulic only |       6 |   0.9815 |   0.3272 |  0.6068 |     7.764  |        0.4222 |
| local age      |       6 |   0.9815 |   0.2488 |  0.2256 |     0.6433 |        0.4222 |

### Paired case-block contrasts

| Left       | Comparator     | Metric      |   Cases |   Difference |   CI low |   CI high |
|:-----------|:---------------|:------------|--------:|-------------:|---------:|----------:|
| full sheaf | hydraulic only | PR-AUC      |       6 |      -0.0197 |  -0.0355 |   -0.0039 |
| full sheaf | hydraulic only | Brier       |       6 |      -0.3897 |  -0.4019 |   -0.3801 |
| full sheaf | hydraulic only | Selected F1 |       6 |       0      |   0      |    0      |
| full sheaf | local age      | PR-AUC      |       6 |       0.0586 |   0.0386 |    0.0777 |
| full sheaf | local age      | Brier       |       6 |      -0.0085 |  -0.0098 |   -0.0072 |
| full sheaf | local age      | Selected F1 |       6 |       0      |   0      |    0      |
| full sheaf | age permuted   | PR-AUC      |       6 |      -0.0136 |  -0.0622 |    0.0347 |
| full sheaf | age permuted   | Brier       |       6 |      -0.0207 |  -0.0429 |    0.0014 |
| full sheaf | age permuted   | Selected F1 |       6 |       0      |   0      |    0      |

Complete record: `tables/publication/tableS8_public_pipeline_acceptance.csv`.

## Table S9

**Table S9.** Claim-bearing estimator-diagnostic robust/hybrid contrasts under the primary separately cross-fitted calibration regime. Differences are left minus comparator; lower Brier score and log loss are favourable. Estimator abbreviations are EL, edge-local; OG, original affine-global; OH, original hybrid; RG, robust affine-global; RH, robust hybrid; and PRH, permuted robust hybrid. PRAUC denotes precision-recall area under the curve. The full 560-row matrix, including the shared-calibrator regime and secondary metrics, is `tables/publication/tableS9_robust_hybrid_contrasts.csv`.

| Scenario         | Comparison   | Metric   |   n | Difference [95% CI]        |
|:-----------------|:-------------|:---------|----:|:---------------------------|
| All              | RH vs EL     | PRAUC    | 128 | +0.0200 [+0.0073, +0.0324] |
| All              | RH vs EL     | Brier    | 128 | -0.0015 [-0.0042, +0.0011] |
| All              | RH vs EL     | LogLoss  | 128 | +0.0033 [-0.0034, +0.0101] |
| All              | RH vs PRH    | PRAUC    | 128 | +0.0441 [+0.0321, +0.0568] |
| All              | RH vs PRH    | Brier    | 128 | -0.0123 [-0.0149, -0.0098] |
| All              | RH vs PRH    | LogLoss  | 128 | -0.0311 [-0.0376, -0.0250] |
| All              | OH vs OG     | PRAUC    | 128 | +0.0029 [+0.0010, +0.0047] |
| All              | OH vs OG     | Brier    | 128 | -0.0006 [-0.0011, -0.0002] |
| All              | OH vs OG     | LogLoss  | 128 | -0.0015 [-0.0025, -0.0005] |
| All              | RG vs OG     | PRAUC    | 128 | -0.0048 [-0.0097, +0.0005] |
| All              | RG vs OG     | Brier    | 128 | +0.0029 [+0.0019, +0.0038] |
| All              | RG vs OG     | LogLoss  | 128 | +0.0073 [+0.0049, +0.0097] |
| All              | RH vs OH     | PRAUC    | 128 | -0.0008 [-0.0056, +0.0044] |
| All              | RH vs OH     | Brier    | 128 | +0.0017 [+0.0008, +0.0025] |
| All              | RH vs OH     | LogLoss  | 128 | +0.0043 [+0.0021, +0.0065] |
| Identity limit   | RH vs EL     | PRAUC    |  32 | +0.0126 [-0.0086, +0.0338] |
| Identity limit   | RH vs EL     | Brier    |  32 | +0.0004 [-0.0042, +0.0049] |
| Identity limit   | RH vs EL     | LogLoss  |  32 | +0.0014 [-0.0087, +0.0113] |
| Heterog. affine  | RH vs EL     | PRAUC    |  32 | -0.0100 [-0.0299, +0.0121] |
| Heterog. affine  | RH vs EL     | Brier    |  32 | +0.0040 [-0.0011, +0.0087] |
| Heterog. affine  | RH vs EL     | LogLoss  |  32 | +0.0251 [+0.0121, +0.0373] |
| Incompat. cycles | RH vs EL     | PRAUC    |  32 | +0.0437 [+0.0125, +0.0751] |
| Incompat. cycles | RH vs EL     | Brier    |  32 | -0.0033 [-0.0088, +0.0020] |
| Incompat. cycles | RH vs EL     | LogLoss  |  32 | -0.0056 [-0.0194, +0.0073] |
| Noisy/missing    | RH vs EL     | PRAUC    |  32 | +0.0335 [+0.0110, +0.0559] |
| Noisy/missing    | RH vs EL     | Brier    |  32 | -0.0072 [-0.0127, -0.0017] |
| Noisy/missing    | RH vs EL     | LogLoss  |  32 | -0.0076 [-0.0222, +0.0061] |

## Table S10

**Table S10.** Selected representation-benchmark contrasts with simultaneous 95% intervals controlling all 120 published contrasts as one family. Complete record: `tables/publication/tableS10_m7_4_multiplicity_adjusted.csv`.

| Scenario            | Comparison                            | Metric                       |   n | Difference [simultaneous 95% CI]   | FWER support   |
|:--------------------|:--------------------------------------|:-----------------------------|----:|:-----------------------------------|:---------------|
| all                 | affine sheaf vs weighted graph        | pr_auc                       |  64 | +0.0097 [-0.0149, +0.0343]         | No             |
| all                 | affine sheaf vs weighted graph        | brier                        |  64 | +0.0005 [-0.0058, +0.0069]         | No             |
| all                 | affine sheaf vs weighted graph        | log_loss                     |  64 | +0.0117 [-0.0066, +0.0301]         | No             |
| all                 | affine sheaf vs graph laplacian       | pr_auc                       |  64 | +0.0854 [+0.0539, +0.1169]         | Yes            |
| all                 | affine sheaf vs graph laplacian       | brier                        |  64 | -0.0193 [-0.0261, -0.0126]         | Yes            |
| all                 | affine sheaf vs graph laplacian       | log_loss                     |  64 | -0.0472 [-0.0638, -0.0306]         | Yes            |
| all                 | affine sheaf vs affine sheaf permuted | pr_auc                       |  64 | +0.0909 [+0.0571, +0.1246]         | Yes            |
| all                 | affine sheaf vs affine sheaf permuted | brier                        |  64 | -0.0215 [-0.0292, -0.0139]         | Yes            |
| all                 | affine sheaf vs affine sheaf permuted | log_loss                     |  64 | -0.0546 [-0.0734, -0.0357]         | Yes            |
| incompatible_cycles | affine sheaf vs weighted graph        | pr_auc                       |  16 | +0.0483 [+0.0095, +0.0871]         | Yes            |
| incompatible_cycles | affine sheaf vs weighted graph        | conflict_localisation_pr_auc |  16 | +0.0689 [+0.0318, +0.1061]         | Yes            |

## Table S11

**Table S11.** Selected estimator-diagnostic contrasts with simultaneous 95% intervals controlling all 560 published contrasts as one family. The selected robust-hybrid arm had local weight 1.0 and is the local-first/global-fallback estimator in the main text. Complete record: `tables/publication/tableS11_m7_5_multiplicity_adjusted.csv`.

| Scenario   | Comparison                                     | Metric                       |   n | Difference [simultaneous 95% CI]   | FWER support   |
|:-----------|:-----------------------------------------------|:-----------------------------|----:|:-----------------------------------|:---------------|
| all        | robust hybrid vs edge local                    | pr_auc                       | 128 | +0.0200 [-0.0055, +0.0454]         | No             |
| all        | robust hybrid vs edge local                    | brier                        | 128 | -0.0015 [-0.0069, +0.0039]         | No             |
| all        | robust hybrid vs edge local                    | log_loss                     | 128 | +0.0033 [-0.0102, +0.0168]         | No             |
| all        | robust hybrid vs edge local                    | conflict_localisation_pr_auc |  32 | +0.0770 [+0.0389, +0.1151]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | pr_auc                       | 128 | +0.0441 [+0.0192, +0.0690]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | brier                        | 128 | -0.0123 [-0.0175, -0.0071]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | log_loss                     | 128 | -0.0311 [-0.0440, -0.0183]         | Yes            |
| all        | robust hybrid vs robust hybrid permuted        | conflict_localisation_pr_auc |  32 | +0.0505 [+0.0154, +0.0855]         | Yes            |
| all        | robust affine global vs original affine global | brier                        | 128 | +0.0029 [+0.0010, +0.0049]         | No             |
| all        | robust affine global vs original affine global | log_loss                     | 128 | +0.0073 [+0.0025, +0.0122]         | No             |

## Table S12

**Table S12.** Post-review empirical precision and future-replication planning audit (20,000 simulations). Margins were not prespecified or field validated; POST_TEST rows are not evidence for completed tests. Complete record: `tables/publication/tableS12_precision_and_power.csv`.

| Design                   | Metric         |   Source n |   Planned n |   Margin |   P(CI favourable) |   P(CI clears margin) | Status               |
|:-------------------------|:---------------|-----------:|------------:|---------:|-------------------:|----------------------:|:---------------------|
| Evidence panel           | pr_auc         |          6 |          12 |     0.02 |             1      |                1      | Development planning |
| Evidence panel           | brier          |          6 |          12 |     0.01 |             1      |                1      | Development planning |
| Evidence panel           | log_loss       |          6 |          12 |     0.02 |             1      |                1      | Development planning |
| Representation           | pr_auc         |         32 |          64 |     0.02 |             0.1928 |                0.0001 | Development planning |
| Representation           | brier          |         32 |          64 |     0.01 |             0.0138 |                0      | Development planning |
| Representation           | log_loss       |         32 |          64 |     0.02 |             0.0001 |                0      | Development planning |
| Estimator                | pr_auc         |         64 |         128 |     0.02 |             0.5475 |                0.0006 | Development planning |
| Estimator                | brier          |         64 |         128 |     0.01 |             0.1228 |                0      | Development planning |
| Estimator                | log_loss       |         64 |         128 |     0.02 |             0.0022 |                0      | Development planning |
| Topology-age: two-tracer | age MAE        |         12 |          12 |     0.25 |             1      |                0      | Replication planning |
| Topology-age: two-tracer | interval width |         12 |          12 |     0.5  |             1      |                0      | Replication planning |
| Topology-age: tritium    | age MAE        |         12 |          12 |     0.25 |             1      |                0      | Replication planning |
| Topology-age: tritium    | interval width |         12 |          12 |     0.5  |             1      |                1      | Replication planning |
| Reaction recovery        | modal accuracy |         12 |          12 |     0.1  |             0.9268 |                0.0004 | Replication planning |

## Table S13

**Table S13.** Public-pipeline selection and confusion-count audit. All generated candidates were retained; no scalar probability threshold was applied. Complete record: `tables/publication/tableS13_public_pipeline_selection.csv`.

| Condition      |   Candidates |   Selected |   Min probability |   TP |   FP |   Conditional FN |   End-to-end FN |   End-to-end F1 |   Candidate recall |
|:---------------|-------------:|-----------:|------------------:|-----:|-----:|-----------------:|----------------:|----------------:|-------------------:|
| age_permuted   |          198 |        198 |            0.0222 |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |
| full_sheaf     |          198 |        198 |            0.1921 |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |
| hydraulic_only |          198 |        198 |            0.0257 |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |
| local_age      |          198 |        198 |            0.317  |   53 |  145 |                0 |               1 |          0.4206 |             0.9815 |

## Table S14

**Table S14.** Auxiliary M7.6 controlled-synthetic M3-mechanism diagnostic. The complete result record is `tables/publication/tableS14_m7_6_m3_mechanism.csv`; this table is not field validation and does not identify the USGS cause.

| Contrast                                            |   Difference |   CI low |   CI high |   Cases |   Bootstrap | Decision                                           |
|:----------------------------------------------------|-------------:|---------:|----------:|--------:|------------:|:---------------------------------------------------|
| severe minus none full infeasibility                |       0.2882 |   0.2118 |    0.3646 |      12 |       10000 | shared nuisance increased full-panel infeasibility |
| reducing minus nonreducing cfc11 pair rate severe   |       0.7188 |   0.6667 |    0.75   |      12 |       10000 | CFC-11 contrast positive                           |
| reducing minus nonreducing cfc12 pair rate severe   |       0.7396 |   0.724  |    0.75   |      12 |       10000 | CFC-12 specificity control failed                  |
| reducing minus nonreducing tritium pair rate severe |       0.3229 |   0.1354 |    0.4896 |      12 |       10000 | positive association not selective to CFC family   |
| E added to N T2 MAE none                            |       0      |   0      |    0      |      12 |       10000 | binding isotope-age control passed                 |
