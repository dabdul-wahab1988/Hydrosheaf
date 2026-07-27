# Supplementary Figures and Tables

**Evidence integration reduces groundwater interpretive non-uniqueness only conditionally: an independent benchmark with adverse controls**

This document reports the supplementary figure and supplementary tables
referenced in the main manuscript and its Supplementary Methods. Tables S1
through S5 are transcribed directly from the locked M7.3 benchmark outputs
(`results/m7_3_locked/`); no locked value has been recomputed for those
tables. Table S6 reports an additional multiplicity-correction robustness
check computed from the same locked per-case metrics in response to
pre-submission peer review; it re-analyses existing locked values with an
exact permutation test and a Benjamini-Hochberg correction rather than
introducing any new experimental result. Full metric and column
definitions are given in Supplementary Methods.

## Figure S1

![Figure S1. Locked synthetic model domain (representative realization 4101). The MODFLOW/MODPATH truth network, observation nodes and hydraulic heads are shown in synthetic model-space coordinates. This is not a geographic map and must not be presented on a Ghana or other geographic basemap.](figures/supporting_validation/figure_s1_model_domain_map.png)

**Figure S1.** Locked synthetic model domain (representative realization
4101). The MODFLOW/MODPATH truth network, observation nodes and hydraulic
heads are shown in synthetic model-space coordinates. This is not a
geographic map and must not be presented on a Ghana or other geographic
basemap.

## Table S1

**Table S1.** Every evidence panel (H, A, C, HA, HC, AC, HAC) reported
under every native and adverse locked-test condition (native,
age_permuted, hydraulic_permuted, joint_misspecified). Metrics are PR-AUC,
ROC-AUC, Brier score, log loss, mean edge entropy, the fraction of edges
with probability between 0.1 and 0.9, mean expected edge count per case,
expected calibration error, and overconfident error fraction.

| condition | panel | pr_auc | roc_auc | brier | log_loss | mean_edge_entropy | ambiguous_edge_fraction | mean_expected_edges_per_case | expected_calibration_error | overconfident_error_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| native | H | 0.1765 | 0.6495 | 0.1076 | 0.3467 | 0.5029 | 0.74 | 8.8313 | 0.0059 | 0 |
| native | A | 0.1112 | 0.4523 | 0.1087 | 0.3558 | 0.524 | 0.7666 | 8.9395 | 0.0124 | 0 |
| native | C | 0.4587 | 0.8916 | 0.0944 | 0.2951 | 0.4551 | 0.5768 | 8.3079 | 0.0777 | 0 |
| native | HA | 0.1114 | 0.4533 | 0.1076 | 0.3467 | 0.5042 | 0.74 | 8.9004 | 0.0046 | 0 |
| native | HC | 0.4846 | 0.9089 | 0.088 | 0.2676 | 0.4221 | 0.5018 | 8.4146 | 0.0602 | 0 |
| native | AC | 0.471 | 0.9036 | 0.089 | 0.2736 | 0.4361 | 0.5127 | 8.4794 | 0.0692 | 0 |
| native | HAC | 0.4805 | 0.9077 | 0.088 | 0.2677 | 0.4215 | 0.5018 | 8.4263 | 0.0599 | 0 |
| age_permuted | H | 0.1765 | 0.6495 | 0.1076 | 0.3467 | 0.5029 | 0.74 | 8.8313 | 0.0059 | 0 |
| age_permuted | A | 0.1397 | 0.4846 | 0.1179 | 0.4539 | 0.524 | 0.7666 | 8.9395 | 0.0553 | 0.0339 |
| age_permuted | C | 0.4587 | 0.8916 | 0.0944 | 0.2951 | 0.4551 | 0.5768 | 8.3079 | 0.0777 | 0 |
| age_permuted | HA | 0.1825 | 0.634 | 0.1096 | 0.3588 | 0.4726 | 0.6481 | 8.0402 | 0.0139 | 0.0193 |
| age_permuted | HC | 0.4846 | 0.9089 | 0.088 | 0.2676 | 0.4221 | 0.5018 | 8.4146 | 0.0602 | 0 |
| age_permuted | AC | 0.3521 | 0.7841 | 0.1016 | 0.3588 | 0.4236 | 0.503 | 8.0302 | 0.0621 | 0.0278 |
| age_permuted | HAC | 0.4165 | 0.875 | 0.0914 | 0.2781 | 0.4014 | 0.4861 | 7.7208 | 0.0397 | 0.0036 |
| hydraulic_permuted | H | 0.1302 | 0.5008 | 0.1191 | 0.467 | 0.5029 | 0.74 | 8.8313 | 0.0642 | 0.0351 |
| hydraulic_permuted | A | 0.1112 | 0.4523 | 0.1087 | 0.3558 | 0.524 | 0.7666 | 8.9395 | 0.0124 | 0 |
| hydraulic_permuted | C | 0.4587 | 0.8916 | 0.0944 | 0.2951 | 0.4551 | 0.5768 | 8.3079 | 0.0777 | 0 |
| hydraulic_permuted | HA | 0.1173 | 0.4814 | 0.1155 | 0.4358 | 0.4713 | 0.6445 | 8.0083 | 0.0375 | 0.0351 |
| hydraulic_permuted | HC | 0.3803 | 0.7878 | 0.1032 | 0.3717 | 0.4054 | 0.4788 | 7.9093 | 0.064 | 0.0351 |
| hydraulic_permuted | AC | 0.471 | 0.9036 | 0.089 | 0.2736 | 0.4361 | 0.5127 | 8.4794 | 0.0692 | 0 |
| hydraulic_permuted | HAC | 0.3971 | 0.8037 | 0.1001 | 0.3482 | 0.3879 | 0.4426 | 7.3766 | 0.0741 | 0.0351 |
| joint_misspecified | H | 0.1302 | 0.5008 | 0.1191 | 0.467 | 0.5029 | 0.74 | 8.8313 | 0.0642 | 0.0351 |
| joint_misspecified | A | 0.1397 | 0.4846 | 0.1179 | 0.4539 | 0.524 | 0.7666 | 8.9395 | 0.0553 | 0.0339 |
| joint_misspecified | C | 0.4587 | 0.8916 | 0.0944 | 0.2951 | 0.4551 | 0.5768 | 8.3079 | 0.0777 | 0 |
| joint_misspecified | HA | 0.1371 | 0.4897 | 0.1201 | 0.463 | 0.4729 | 0.6481 | 8.0502 | 0.0722 | 0.052 |
| joint_misspecified | HC | 0.3803 | 0.7878 | 0.1032 | 0.3717 | 0.4054 | 0.4788 | 7.9093 | 0.064 | 0.0351 |
| joint_misspecified | AC | 0.3521 | 0.7841 | 0.1016 | 0.3588 | 0.4236 | 0.503 | 8.0302 | 0.0621 | 0.0278 |
| joint_misspecified | HAC | 0.3305 | 0.7613 | 0.105 | 0.3682 | 0.3845 | 0.4534 | 7.2063 | 0.0522 | 0.0387 |
## Table S2

**Table S2.** All case-block evidence contrasts computed for the
benchmark (10,000-replicate case-block bootstrap; `resampling_unit` is an
independent MODFLOW case), extending Table 3 of the main text to every
predeclared contrast and metric.

| contrast | condition | full_panel | baseline_panel | metric | mean_difference | ci95_low | ci95_high | n_cases | n_bootstrap | resampling_unit |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| native_incremental_age | native | HAC | HC | pr_auc | -0.006 | -0.0122 | -0.0011 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_age | native | HAC | HC | brier | 0.0001 | -0 | 0.0001 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_age | native | HAC | HC | log_loss | 0 | -0.0001 | 0.0002 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_age | native | HAC | HC | mean_edge_entropy | -0.0006 | -0.0009 | -0.0003 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_age | native | HAC | HC | overconfident_error_fraction | 0 | 0 | 0 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_chemistry | native | HAC | HA | pr_auc | 0.4471 | 0.3575 | 0.5401 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_chemistry | native | HAC | HA | brier | -0.0196 | -0.0213 | -0.0176 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_chemistry | native | HAC | HA | log_loss | -0.0791 | -0.085 | -0.072 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_chemistry | native | HAC | HA | mean_edge_entropy | -0.0827 | -0.0985 | -0.0653 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_chemistry | native | HAC | HA | overconfident_error_fraction | 0 | 0 | 0 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_hydraulics | native | HAC | AC | pr_auc | 0.0091 | 0.001 | 0.0198 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_hydraulics | native | HAC | AC | brier | -0.001 | -0.0012 | -0.0008 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_hydraulics | native | HAC | AC | log_loss | -0.0059 | -0.007 | -0.005 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_hydraulics | native | HAC | AC | mean_edge_entropy | -0.0146 | -0.0185 | -0.0118 | 12 | 10000 | independent_MODFLOW_case |
| native_incremental_hydraulics | native | HAC | AC | overconfident_error_fraction | 0 | 0 | 0 | 12 | 10000 | independent_MODFLOW_case |
| permuted_age_increment | age_permuted | HAC | HC | pr_auc | -0.0754 | -0.1353 | -0.0148 | 12 | 10000 | independent_MODFLOW_case |
| permuted_age_increment | age_permuted | HAC | HC | brier | 0.0034 | 0.0012 | 0.0055 | 12 | 10000 | independent_MODFLOW_case |
| permuted_age_increment | age_permuted | HAC | HC | log_loss | 0.0105 | 0.003 | 0.0176 | 12 | 10000 | independent_MODFLOW_case |
| permuted_age_increment | age_permuted | HAC | HC | mean_edge_entropy | -0.0207 | -0.0236 | -0.0175 | 12 | 10000 | independent_MODFLOW_case |
| permuted_age_increment | age_permuted | HAC | HC | overconfident_error_fraction | 0.0036 | 0 | 0.0073 | 12 | 10000 | independent_MODFLOW_case |
| permuted_hydraulic_increment | hydraulic_permuted | HAC | AC | pr_auc | -0.0686 | -0.112 | -0.0271 | 12 | 10000 | independent_MODFLOW_case |
| permuted_hydraulic_increment | hydraulic_permuted | HAC | AC | brier | 0.011 | 0.0055 | 0.0164 | 12 | 10000 | independent_MODFLOW_case |
| permuted_hydraulic_increment | hydraulic_permuted | HAC | AC | log_loss | 0.0745 | 0.0404 | 0.1091 | 12 | 10000 | independent_MODFLOW_case |
| permuted_hydraulic_increment | hydraulic_permuted | HAC | AC | mean_edge_entropy | -0.0482 | -0.0593 | -0.0372 | 12 | 10000 | independent_MODFLOW_case |
| permuted_hydraulic_increment | hydraulic_permuted | HAC | AC | overconfident_error_fraction | 0.035 | 0.0229 | 0.0483 | 12 | 10000 | independent_MODFLOW_case |
| joint_misspecification | joint_misspecified | HAC | C | pr_auc | -0.139 | -0.2041 | -0.0742 | 12 | 10000 | independent_MODFLOW_case |
| joint_misspecification | joint_misspecified | HAC | C | brier | 0.0106 | 0.0055 | 0.0163 | 12 | 10000 | independent_MODFLOW_case |
| joint_misspecification | joint_misspecified | HAC | C | log_loss | 0.073 | 0.0401 | 0.1073 | 12 | 10000 | independent_MODFLOW_case |
| joint_misspecification | joint_misspecified | HAC | C | mean_edge_entropy | -0.0706 | -0.0827 | -0.0586 | 12 | 10000 | independent_MODFLOW_case |
| joint_misspecification | joint_misspecified | HAC | C | overconfident_error_fraction | 0.0387 | 0.0278 | 0.0495 | 12 | 10000 | independent_MODFLOW_case |
## Table S3

**Table S3.** Per-case, per-tracer-regime topology-to-age sensitivity for
each of the four topology conditions (none, partial_true, complete_true,
reversed), including age MAE, bias, 95% coverage, mean interval width,
normalised marginal age entropy, true- and assumed-edge order-violation
probabilities, importance-sampling effective sample size (`importance_ess`)
and its fraction of 50,000 particles, log mean topology weight, and the
predeclared stability flag (`importance_stable_ess_ge_400`).

| seed | tracer_regime | graph_condition | n_nodes | n_assumed_edges | n_particles | order_scale_years | age_mae_years | age_bias_years | age_95_coverage | mean_interval_width_years | mean_normalized_age_entropy | true_edge_order_violation_probability | assumed_edge_order_violation_probability | importance_ess | importance_ess_fraction | log_mean_topology_weight | importance_stable_ess_ge_400 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5301 | informative | none | 12 | 0 | 50000 | 5 | 2.9214 | -0.1678 | 0.9167 | 12.6875 | 0.5693 | 0.0044 |  | 50000 | 1 | 0 | True |
| 5301 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.8986 | -0.1519 | 0.9167 | 12.5833 | 0.5684 | 0.0043 | 0.0005 | 49259.3832 | 0.9852 | -0.2668 | True |
| 5301 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.8404 | -0.1635 | 0.9167 | 12.3958 | 0.5659 | 0.0022 | 0.0022 | 48121.2754 | 0.9624 | -0.7499 | True |
| 5301 | informative | reversed | 12 | 9 | 50000 | 5 | 2.8358 | -1.1544 | 0.9167 | 12.8125 | 0.5672 | 0.0037 | 0.9955 | 564.3485 | 0.0113 | -29.0454 | True |
| 5301 | tritium_only | none | 12 | 0 | 50000 | 5 | 5.1513 | 2.3084 | 1 | 27.6042 | 0.6135 | 0.0138 |  | 50000 | 1 | 0 | True |
| 5301 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 5.0975 | 2.2587 | 1 | 27.2292 | 0.6099 | 0.0102 | 0.0036 | 47580.4461 | 0.9516 | -0.3779 | True |
| 5301 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 5.0061 | 2.196 | 1 | 26.7981 | 0.6057 | 0.0051 | 0.0051 | 45510.9776 | 0.9102 | -0.8728 | True |
| 5301 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.6102 | -1.102 | 1 | 19.7935 | 0.5514 | 0.0449 | 0.9545 | 32.8218 | 0.0007 | -28.462 | False |
| 5302 | informative | none | 12 | 0 | 50000 | 5 | 2.6515 | -0.8582 | 0.9167 | 13.375 | 0.5724 | 0.0014 |  | 50000 | 1 | 0 | True |
| 5302 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.6345 | -0.8491 | 0.9167 | 13.2292 | 0.5706 | 0.0011 | 0.0005 | 49229.9929 | 0.9846 | -0.3142 | True |
| 5302 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.6251 | -0.8457 | 0.9167 | 13.2083 | 0.5696 | 0.0007 | 0.0007 | 48916.4241 | 0.9783 | -0.6364 | True |
| 5302 | informative | reversed | 12 | 9 | 50000 | 5 | 3.6274 | -2.4487 | 0.9167 | 13.404 | 0.5725 | 0.0014 | 0.9984 | 954.4612 | 0.0191 | -31.7807 | True |
| 5302 | tritium_only | none | 12 | 0 | 50000 | 5 | 4.6564 | 1.2414 | 0.9167 | 29.4583 | 0.6173 | 0.0108 |  | 50000 | 1 | 0 | True |
| 5302 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 4.6233 | 1.1854 | 0.9167 | 29.0714 | 0.6132 | 0.0077 | 0.0023 | 47813.8089 | 0.9563 | -0.3592 | True |
| 5302 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 4.4996 | 1.1222 | 0.9167 | 28.6563 | 0.6096 | 0.0036 | 0.0036 | 46145.3991 | 0.9229 | -0.8023 | True |
| 5302 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 5.5508 | -2.2538 | 0.75 | 20.7349 | 0.5903 | 0.0366 | 0.9629 | 248.2882 | 0.005 | -29.7733 | False |
| 5303 | informative | none | 12 | 0 | 50000 | 5 | 2.6957 | -0.6043 | 1 | 12.5833 | 0.5653 | 0.0036 |  | 50000 | 1 | 0 | True |
| 5303 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.6236 | -0.5846 | 1 | 12.4167 | 0.5632 | 0.0023 | 0.002 | 48262.5067 | 0.9653 | -0.4927 | True |
| 5303 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.6289 | -0.5816 | 1 | 12.3167 | 0.5621 | 0.0018 | 0.0018 | 48119.3571 | 0.9624 | -0.8559 | True |
| 5303 | informative | reversed | 12 | 9 | 50000 | 5 | 2.6498 | -1.3407 | 0.9167 | 11.9792 | 0.5643 | 0.0041 | 0.9951 | 4105.8648 | 0.0821 | -26.3517 | True |
| 5303 | tritium_only | none | 12 | 0 | 50000 | 5 | 4.4584 | 0.7264 | 1 | 21.9792 | 0.5965 | 0.0163 |  | 50000 | 1 | 0 | True |
| 5303 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 4.395 | 0.6472 | 1 | 21.5 | 0.5914 | 0.0113 | 0.005 | 46533.0977 | 0.9307 | -0.5625 | True |
| 5303 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 4.266 | 0.5849 | 1 | 21.0223 | 0.5875 | 0.0063 | 0.0063 | 44967.1406 | 0.8993 | -1.068 | True |
| 5303 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.0731 | -0.6618 | 0.9167 | 18.3961 | 0.5905 | 0.0359 | 0.9627 | 754.014 | 0.0151 | -25.559 | True |
| 5304 | informative | none | 12 | 0 | 50000 | 5 | 2.2764 | -1.1384 | 0.9167 | 11.8542 | 0.5628 | 0.0027 |  | 50000 | 1 | 0 | True |
| 5304 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.2111 | -1.1201 | 0.9167 | 11.7083 | 0.5612 | 0.0021 | 0.001 | 48837.9073 | 0.9768 | -0.3843 | True |
| 5304 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.2209 | -1.1242 | 0.9167 | 11.5417 | 0.5594 | 0.0013 | 0.0013 | 48510.0849 | 0.9702 | -0.778 | True |
| 5304 | informative | reversed | 12 | 9 | 50000 | 5 | 2.9835 | -1.6999 | 0.9167 | 11.9947 | 0.5654 | 0.0036 | 0.9959 | 4668.1116 | 0.0934 | -26.4198 | True |
| 5304 | tritium_only | none | 12 | 0 | 50000 | 5 | 2.6069 | -0.0593 | 1 | 18.8333 | 0.5943 | 0.0144 |  | 50000 | 1 | 0 | True |
| 5304 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 2.5633 | -0.1217 | 1 | 18.3333 | 0.5898 | 0.0098 | 0.005 | 46937.871 | 0.9388 | -0.4937 | True |
| 5304 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 2.4886 | -0.1822 | 1 | 17.875 | 0.5858 | 0.0053 | 0.0053 | 45436.5489 | 0.9087 | -0.9505 | True |
| 5304 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.819 | -0.4341 | 1 | 18.8496 | 0.5716 | 0.0651 | 0.9342 | 75.9644 | 0.0015 | -25.2794 | False |
| 5305 | informative | none | 12 | 0 | 50000 | 5 | 3.253 | 0.0606 | 0.8333 | 13.8958 | 0.577 | 0.0044 |  | 50000 | 1 | 0 | True |
| 5305 | informative | partial_true | 12 | 5 | 50000 | 5 | 3.2041 | 0.0679 | 0.8333 | 13.7917 | 0.5757 | 0.0027 | 0.0034 | 48463.7215 | 0.9693 | -0.3978 | True |
| 5305 | informative | complete_true | 12 | 9 | 50000 | 5 | 3.1784 | 0.0774 | 0.8333 | 13.7083 | 0.5748 | 0.0023 | 0.0023 | 48127.2502 | 0.9625 | -0.6936 | True |
| 5305 | informative | reversed | 12 | 9 | 50000 | 5 | 3.0032 | -1.2698 | 0.8333 | 13.5103 | 0.5752 | 0.0043 | 0.9953 | 1491.3018 | 0.0298 | -31.2682 | True |
| 5305 | tritium_only | none | 12 | 0 | 50000 | 5 | 7.7008 | 4.3802 | 0.9167 | 34.75 | 0.6284 | 0.0192 |  | 50000 | 1 | 0 | True |
| 5305 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 7.614 | 4.33 | 1 | 34.2708 | 0.6248 | 0.0124 | 0.0107 | 46553.8231 | 0.9311 | -0.4929 | True |
| 5305 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 7.5053 | 4.271 | 1 | 33.875 | 0.6214 | 0.0083 | 0.0083 | 44830.2696 | 0.8966 | -0.904 | True |
| 5305 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 5.1836 | -0.6503 | 0.9167 | 25.6878 | 0.5392 | 0.0403 | 0.9523 | 38.2477 | 0.0008 | -31.5395 | False |
| 5306 | informative | none | 12 | 0 | 50000 | 5 | 2.8177 | 0.3147 | 1 | 13.1667 | 0.5713 | 0.0022 |  | 50000 | 1 | 0 | True |
| 5306 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.8107 | 0.3175 | 1 | 12.9583 | 0.5695 | 0.0017 | 0.0008 | 49175.9379 | 0.9835 | -0.2853 | True |
| 5306 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.7603 | 0.3269 | 1 | 12.8333 | 0.5684 | 0.001 | 0.001 | 48605.1488 | 0.9721 | -0.6402 | True |
| 5306 | informative | reversed | 12 | 9 | 50000 | 5 | 2.8517 | -0.8052 | 1 | 13.5138 | 0.5751 | 0.0028 | 0.9969 | 1730.1098 | 0.0346 | -29.9909 | True |
| 5306 | tritium_only | none | 12 | 0 | 50000 | 5 | 6.2893 | 3.7109 | 1 | 33.6458 | 0.6236 | 0.0119 |  | 50000 | 1 | 0 | True |
| 5306 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 6.2277 | 3.6372 | 1 | 33.1875 | 0.6197 | 0.0087 | 0.0024 | 47892.801 | 0.9579 | -0.3131 | True |
| 5306 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 6.1208 | 3.5815 | 1 | 32.7292 | 0.6158 | 0.0039 | 0.0039 | 45796.5365 | 0.9159 | -0.7692 | True |
| 5306 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.1661 | 0.8172 | 1 | 30.891 | 0.5943 | 0.0427 | 0.9572 | 207.5633 | 0.0042 | -30.5755 | False |
| 5307 | informative | none | 12 | 0 | 50000 | 5 | 3.4124 | 0.0535 | 0.8333 | 13.7083 | 0.5758 | 0.0096 |  | 50000 | 1 | 0 | True |
| 5307 | informative | partial_true | 12 | 5 | 50000 | 5 | 3.3683 | 0.0691 | 0.9167 | 13.6339 | 0.5747 | 0.0064 | 0.0074 | 47904.1384 | 0.9581 | -0.4599 | True |
| 5307 | informative | complete_true | 12 | 9 | 50000 | 5 | 3.3183 | 0.073 | 0.9167 | 13.5417 | 0.5737 | 0.0054 | 0.0054 | 47307.2804 | 0.9461 | -0.8258 | True |
| 5307 | informative | reversed | 12 | 9 | 50000 | 5 | 2.7849 | -1.2986 | 0.9167 | 13.8273 | 0.5763 | 0.0078 | 0.9912 | 1397.0779 | 0.0279 | -30.5794 | True |
| 5307 | tritium_only | none | 12 | 0 | 50000 | 5 | 7.786 | 4.0946 | 0.75 | 34.75 | 0.6267 | 0.0312 |  | 50000 | 1 | 0 | True |
| 5307 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 7.675 | 4.0454 | 0.75 | 34.2917 | 0.6229 | 0.0213 | 0.0215 | 45652.2945 | 0.913 | -0.5932 | True |
| 5307 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 7.5485 | 3.9741 | 0.75 | 33.8542 | 0.6193 | 0.016 | 0.016 | 43637.7024 | 0.8728 | -1.0718 | True |
| 5307 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.7223 | -0.2895 | 0.75 | 23.4318 | 0.5411 | 0.047 | 0.9515 | 34.1509 | 0.0007 | -31.4872 | False |
| 5308 | informative | none | 12 | 0 | 50000 | 5 | 2.8175 | -0.0671 | 0.9167 | 13.0417 | 0.5704 | 0.0047 |  | 50000 | 1 | 0 | True |
| 5308 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.8328 | -0.0651 | 0.9167 | 12.9375 | 0.5691 | 0.0046 | 0.0003 | 49516.4618 | 0.9903 | -0.2197 | True |
| 5308 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.7525 | -0.0549 | 0.9167 | 12.7917 | 0.5677 | 0.0024 | 0.0024 | 48265.6895 | 0.9653 | -0.6969 | True |
| 5308 | informative | reversed | 12 | 9 | 50000 | 5 | 2.5319 | -1.06 | 0.9167 | 12.6458 | 0.5705 | 0.0046 | 0.9948 | 2584.7145 | 0.0517 | -30.2671 | True |
| 5308 | tritium_only | none | 12 | 0 | 50000 | 5 | 4.5063 | 1.445 | 1 | 28.9583 | 0.6131 | 0.0145 |  | 50000 | 1 | 0 | True |
| 5308 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 4.4696 | 1.3837 | 0.9167 | 28.6042 | 0.6094 | 0.0119 | 0.0018 | 48178.8895 | 0.9636 | -0.2827 | True |
| 5308 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 4.3409 | 1.322 | 1 | 28.125 | 0.6052 | 0.0055 | 0.0055 | 45500.2716 | 0.91 | -0.8489 | True |
| 5308 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.3007 | -0.6485 | 1 | 22.1781 | 0.5979 | 0.0319 | 0.967 | 579.6478 | 0.0116 | -28.847 | True |
| 5309 | informative | none | 12 | 0 | 50000 | 5 | 2.5882 | -1.3799 | 0.9167 | 10.7083 | 0.551 | 0.0032 |  | 50000 | 1 | 0 | True |
| 5309 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.5597 | -1.3673 | 0.8333 | 10.5208 | 0.5489 | 0.0028 | 0.0008 | 49001.1146 | 0.98 | -0.3813 | True |
| 5309 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.5406 | -1.367 | 1 | 10.4375 | 0.5469 | 0.0016 | 0.0016 | 48368.3283 | 0.9674 | -0.9501 | True |
| 5309 | informative | reversed | 12 | 9 | 50000 | 5 | 2.9282 | -1.7204 | 0.9167 | 11 | 0.5546 | 0.0059 | 0.9936 | 8633.246 | 0.1727 | -23.5211 | True |
| 5309 | tritium_only | none | 12 | 0 | 50000 | 5 | 2.9465 | -0.2117 | 0.9167 | 16.7708 | 0.5826 | 0.0153 |  | 50000 | 1 | 0 | True |
| 5309 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 2.9585 | -0.2771 | 0.9167 | 16.2697 | 0.578 | 0.0122 | 0.0026 | 47514.3662 | 0.9503 | -0.4158 | True |
| 5309 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 2.8187 | -0.3572 | 1 | 15.6875 | 0.5727 | 0.0055 | 0.0055 | 45093.8969 | 0.9019 | -1.109 | True |
| 5309 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.0763 | -0.4298 | 0.9167 | 17.3719 | 0.5801 | 0.0439 | 0.9517 | 297.0407 | 0.0059 | -22.7947 | False |
| 5310 | informative | none | 12 | 0 | 50000 | 5 | 2.5983 | -1.1307 | 1 | 11.375 | 0.5562 | 0.0041 |  | 50000 | 1 | 0 | True |
| 5310 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.5797 | -1.1232 | 1 | 11.2083 | 0.5543 | 0.0032 | 0.0014 | 48561.2731 | 0.9712 | -0.4526 | True |
| 5310 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.5332 | -1.1196 | 1 | 11.0833 | 0.5526 | 0.002 | 0.002 | 48160.6321 | 0.9632 | -0.8857 | True |
| 5310 | informative | reversed | 12 | 9 | 50000 | 5 | 3.1975 | -1.5995 | 0.9167 | 11.7087 | 0.5594 | 0.0059 | 0.9936 | 4561.9815 | 0.0912 | -24.7873 | True |
| 5310 | tritium_only | none | 12 | 0 | 50000 | 5 | 3.224 | -0.2395 | 1 | 18.0208 | 0.5875 | 0.016 |  | 50000 | 1 | 0 | True |
| 5310 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 3.1706 | -0.3153 | 1 | 17.4583 | 0.5825 | 0.0114 | 0.0042 | 46693.3432 | 0.9339 | -0.5453 | True |
| 5310 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 3.1147 | -0.3815 | 1 | 16.9597 | 0.5781 | 0.0059 | 0.0059 | 45033.3495 | 0.9007 | -1.0681 | True |
| 5310 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.5414 | -0.6415 | 0.9167 | 19.1639 | 0.591 | 0.0375 | 0.9612 | 1346.4481 | 0.0269 | -23.9975 | True |
| 5311 | informative | none | 12 | 0 | 50000 | 5 | 2.4467 | -1.3108 | 1 | 11.2292 | 0.5565 | 0.0029 |  | 50000 | 1 | 0 | True |
| 5311 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.4661 | -1.3054 | 0.9167 | 11.0833 | 0.5544 | 0.0022 | 0.0011 | 48717.5762 | 0.9744 | -0.4278 | True |
| 5311 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.3998 | -1.3023 | 1 | 10.9792 | 0.553 | 0.0014 | 0.0014 | 48417.857 | 0.9684 | -0.8411 | True |
| 5311 | informative | reversed | 12 | 9 | 50000 | 5 | 2.9987 | -1.7544 | 0.9167 | 11.4583 | 0.56 | 0.003 | 0.9966 | 5912.4328 | 0.1182 | -25.2297 | True |
| 5311 | tritium_only | none | 12 | 0 | 50000 | 5 | 3.3624 | -0.1703 | 1 | 18.4792 | 0.5902 | 0.0158 |  | 50000 | 1 | 0 | True |
| 5311 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 3.291 | -0.2413 | 1 | 18.0597 | 0.5854 | 0.0116 | 0.004 | 46951.6455 | 0.939 | -0.4931 | True |
| 5311 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 3.1856 | -0.3058 | 1 | 17.5417 | 0.5811 | 0.006 | 0.006 | 45076.2973 | 0.9015 | -1.035 | True |
| 5311 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.7639 | -0.8973 | 0.8333 | 18.4051 | 0.5771 | 0.0421 | 0.9564 | 110.1378 | 0.0022 | -24.3984 | False |
| 5312 | informative | none | 12 | 0 | 50000 | 5 | 2.694 | -1.1152 | 0.9167 | 12.4167 | 0.567 | 0.0034 |  | 50000 | 1 | 0 | True |
| 5312 | informative | partial_true | 12 | 5 | 50000 | 5 | 2.6891 | -1.1055 | 0.9167 | 12.2292 | 0.5648 | 0.0028 | 0.001 | 48894.9819 | 0.9779 | -0.3749 | True |
| 5312 | informative | complete_true | 12 | 9 | 50000 | 5 | 2.6314 | -1.097 | 0.9167 | 12.1756 | 0.5637 | 0.0017 | 0.0017 | 48336.5781 | 0.9667 | -0.7775 | True |
| 5312 | informative | reversed | 12 | 9 | 50000 | 5 | 3.4223 | -1.9804 | 0.9167 | 12.7083 | 0.5707 | 0.0042 | 0.9955 | 2531.2806 | 0.0506 | -27.05 | True |
| 5312 | tritium_only | none | 12 | 0 | 50000 | 5 | 4.3137 | 1.2101 | 1 | 25.1667 | 0.6071 | 0.0161 |  | 50000 | 1 | 0 | True |
| 5312 | tritium_only | partial_true | 12 | 5 | 50000 | 5 | 4.2264 | 1.1424 | 1 | 24.6667 | 0.6029 | 0.0125 | 0.0034 | 47415.3811 | 0.9483 | -0.4118 | True |
| 5312 | tritium_only | complete_true | 12 | 9 | 50000 | 5 | 4.137 | 1.0849 | 1 | 24.3513 | 0.5989 | 0.0063 | 0.0063 | 45140.2629 | 0.9028 | -0.9512 | True |
| 5312 | tritium_only | reversed | 12 | 9 | 50000 | 5 | 4.514 | -0.5243 | 1 | 21.8487 | 0.5957 | 0.0405 | 0.9583 | 414.9644 | 0.0083 | -27.634 | True |
## Table S4

**Table S4.** Edgewise reaction support and entropy for every scored true
flow edge (12 locked test cases, 9 true edges per case, core and enhanced
panels), reporting the true and modal predicted reaction family, whether
the modal family was correct, the probability assigned to the true family,
family-support entropy, and the effective number of supported families.

| seed | edge_id | tier | true_family | true_process | n_bootstrap | modal_family | modal_family_correct | true_family_probability | family_support_entropy | effective_supported_families |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 5301 | MF5301_P0_M0->MF5301_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1391 | 1.1492 |
| 5301 | MF5301_P0_M0->MF5301_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1391 | 1.1492 |
| 5301 | MF5301_P0_M1->MF5301_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5301 | MF5301_P0_M1->MF5301_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5301 | MF5301_P0_M2->MF5301_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.6406 | 0.8739 | 2.3962 |
| 5301 | MF5301_P0_M2->MF5301_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9375 | 0.2338 | 1.2634 |
| 5301 | MF5301_P1_M0->MF5301_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5301 | MF5301_P1_M0->MF5301_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5301 | MF5301_P1_M1->MF5301_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5301 | MF5301_P1_M1->MF5301_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5301 | MF5301_P1_M2->MF5301_P1_M3 | core | iron_reduction | iron_reduction | 64 | sulfate_reduction | False | 0 | 0.6809 | 1.9756 |
| 5301 | MF5301_P1_M2->MF5301_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | iron_reduction | True | 0.7344 | 0.6383 | 1.8932 |
| 5301 | MF5301_P2_M0->MF5301_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9375 | 0.2338 | 1.2634 |
| 5301 | MF5301_P2_M0->MF5301_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5301 | MF5301_P2_M1->MF5301_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.2742 | 1.3154 |
| 5301 | MF5301_P2_M1->MF5301_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1391 | 1.1492 |
| 5301 | MF5301_P2_M2->MF5301_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5301 | MF5301_P2_M2->MF5301_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5302 | MF5302_P0_M0->MF5302_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1391 | 1.1492 |
| 5302 | MF5302_P0_M0->MF5302_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.0805 | 1.0838 |
| 5302 | MF5302_P0_M1->MF5302_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5302 | MF5302_P0_M1->MF5302_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5302 | MF5302_P0_M2->MF5302_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5302 | MF5302_P0_M2->MF5302_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5302 | MF5302_P1_M0->MF5302_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5302 | MF5302_P1_M0->MF5302_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5302 | MF5302_P1_M1->MF5302_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5302 | MF5302_P1_M1->MF5302_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5302 | MF5302_P1_M2->MF5302_P1_M3 | core | iron_reduction | iron_reduction | 64 | sulfate_reduction | False | 0 | 0.6809 | 1.9756 |
| 5302 | MF5302_P1_M2->MF5302_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | iron_reduction | True | 0.8281 | 0.4588 | 1.5822 |
| 5302 | MF5302_P2_M0->MF5302_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9375 | 0.2338 | 1.2634 |
| 5302 | MF5302_P2_M0->MF5302_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5302 | MF5302_P2_M1->MF5302_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5302 | MF5302_P2_M1->MF5302_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5302 | MF5302_P2_M2->MF5302_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 0.9844 | 0.0805 | 1.0838 |
| 5302 | MF5302_P2_M2->MF5302_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 0.9844 | 0.0805 | 1.0838 |
| 5303 | MF5303_P0_M0->MF5303_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.0805 | 1.0838 |
| 5303 | MF5303_P0_M0->MF5303_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5303 | MF5303_P0_M1->MF5303_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5303 | MF5303_P0_M1->MF5303_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5303 | MF5303_P0_M2->MF5303_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8438 | 0.5288 | 1.697 |
| 5303 | MF5303_P0_M2->MF5303_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8594 | 0.4552 | 1.5764 |
| 5303 | MF5303_P1_M0->MF5303_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5303 | MF5303_P1_M0->MF5303_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5303 | MF5303_P1_M1->MF5303_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5303 | MF5303_P1_M1->MF5303_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 0.9844 | 0.0805 | 1.0838 |
| 5303 | MF5303_P1_M2->MF5303_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.5445 | 1.7238 |
| 5303 | MF5303_P1_M2->MF5303_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0.4219 | 0.6809 | 1.9756 |
| 5303 | MF5303_P2_M0->MF5303_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.6875 | 0.6211 | 1.8609 |
| 5303 | MF5303_P2_M0->MF5303_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9531 | 0.1892 | 1.2083 |
| 5303 | MF5303_P2_M1->MF5303_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5303 | MF5303_P2_M1->MF5303_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5303 | MF5303_P2_M2->MF5303_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5303 | MF5303_P2_M2->MF5303_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5304 | MF5304_P0_M0->MF5304_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.3901 | 1.4771 |
| 5304 | MF5304_P0_M0->MF5304_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.3901 | 1.4771 |
| 5304 | MF5304_P0_M1->MF5304_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5304 | MF5304_P0_M1->MF5304_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5304 | MF5304_P0_M2->MF5304_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9062 | 0.3111 | 1.365 |
| 5304 | MF5304_P0_M2->MF5304_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5304 | MF5304_P1_M0->MF5304_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5304 | MF5304_P1_M0->MF5304_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5304 | MF5304_P1_M1->MF5304_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 0.9844 | 0.0805 | 1.0838 |
| 5304 | MF5304_P1_M1->MF5304_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5304 | MF5304_P1_M2->MF5304_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5304 | MF5304_P1_M2->MF5304_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0.1406 | 0.4061 | 1.5009 |
| 5304 | MF5304_P2_M0->MF5304_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.5781 | 0.8579 | 2.3581 |
| 5304 | MF5304_P2_M0->MF5304_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8438 | 0.4842 | 1.6229 |
| 5304 | MF5304_P2_M1->MF5304_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5304 | MF5304_P2_M1->MF5304_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5304 | MF5304_P2_M2->MF5304_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5304 | MF5304_P2_M2->MF5304_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5305 | MF5305_P0_M0->MF5305_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.2338 | 1.2634 |
| 5305 | MF5305_P0_M0->MF5305_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5305 | MF5305_P0_M1->MF5305_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5305 | MF5305_P0_M1->MF5305_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5305 | MF5305_P0_M2->MF5305_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8906 | 0.3452 | 1.4123 |
| 5305 | MF5305_P0_M2->MF5305_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9375 | 0.2338 | 1.2634 |
| 5305 | MF5305_P1_M0->MF5305_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5305 | MF5305_P1_M0->MF5305_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5305 | MF5305_P1_M1->MF5305_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 0.9844 | 0.0805 | 1.0838 |
| 5305 | MF5305_P1_M1->MF5305_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 0.9531 | 0.1892 | 1.2083 |
| 5305 | MF5305_P1_M2->MF5305_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.4826 | 1.6202 |
| 5305 | MF5305_P1_M2->MF5305_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | iron_reduction | True | 0.6875 | 0.6211 | 1.8609 |
| 5305 | MF5305_P2_M0->MF5305_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5305 | MF5305_P2_M0->MF5305_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5305 | MF5305_P2_M1->MF5305_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5305 | MF5305_P2_M1->MF5305_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5305 | MF5305_P2_M2->MF5305_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5305 | MF5305_P2_M2->MF5305_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5306 | MF5306_P0_M0->MF5306_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5306 | MF5306_P0_M0->MF5306_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5306 | MF5306_P0_M1->MF5306_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5306 | MF5306_P0_M1->MF5306_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5306 | MF5306_P0_M2->MF5306_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5306 | MF5306_P0_M2->MF5306_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5306 | MF5306_P1_M0->MF5306_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5306 | MF5306_P1_M0->MF5306_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5306 | MF5306_P1_M1->MF5306_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5306 | MF5306_P1_M1->MF5306_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5306 | MF5306_P1_M2->MF5306_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5306 | MF5306_P1_M2->MF5306_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0.2969 | 0.6082 | 1.8371 |
| 5306 | MF5306_P2_M0->MF5306_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.875 | 0.3768 | 1.4576 |
| 5306 | MF5306_P2_M0->MF5306_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8438 | 0.4334 | 1.5425 |
| 5306 | MF5306_P2_M1->MF5306_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.0805 | 1.0838 |
| 5306 | MF5306_P2_M1->MF5306_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.0805 | 1.0838 |
| 5306 | MF5306_P2_M2->MF5306_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5306 | MF5306_P2_M2->MF5306_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 0.8594 | 0.4061 | 1.5009 |
| 5307 | MF5307_P0_M0->MF5307_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1391 | 1.1492 |
| 5307 | MF5307_P0_M0->MF5307_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.0805 | 1.0838 |
| 5307 | MF5307_P0_M1->MF5307_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5307 | MF5307_P0_M1->MF5307_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5307 | MF5307_P0_M2->MF5307_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.7344 | 0.5789 | 1.784 |
| 5307 | MF5307_P0_M2->MF5307_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5307 | MF5307_P1_M0->MF5307_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5307 | MF5307_P1_M0->MF5307_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5307 | MF5307_P1_M1->MF5307_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5307 | MF5307_P1_M1->MF5307_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 0.9531 | 0.1892 | 1.2083 |
| 5307 | MF5307_P1_M2->MF5307_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.3452 | 1.4123 |
| 5307 | MF5307_P1_M2->MF5307_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0.3594 | 0.6531 | 1.9214 |
| 5307 | MF5307_P2_M0->MF5307_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5307 | MF5307_P2_M0->MF5307_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5307 | MF5307_P2_M1->MF5307_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5307 | MF5307_P2_M1->MF5307_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5307 | MF5307_P2_M2->MF5307_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5307 | MF5307_P2_M2->MF5307_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5308 | MF5308_P0_M0->MF5308_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5308 | MF5308_P0_M0->MF5308_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5308 | MF5308_P0_M1->MF5308_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5308 | MF5308_P0_M1->MF5308_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5308 | MF5308_P0_M2->MF5308_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9531 | 0.219 | 1.2449 |
| 5308 | MF5308_P0_M2->MF5308_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5308 | MF5308_P1_M0->MF5308_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5308 | MF5308_P1_M0->MF5308_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5308 | MF5308_P1_M1->MF5308_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5308 | MF5308_P1_M1->MF5308_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5308 | MF5308_P1_M2->MF5308_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.5623 | 1.7548 |
| 5308 | MF5308_P1_M2->MF5308_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | iron_reduction | True | 0.7969 | 0.5047 | 1.6565 |
| 5308 | MF5308_P2_M0->MF5308_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5308 | MF5308_P2_M0->MF5308_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5308 | MF5308_P2_M1->MF5308_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.2742 | 1.3154 |
| 5308 | MF5308_P2_M1->MF5308_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1391 | 1.1492 |
| 5308 | MF5308_P2_M2->MF5308_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5308 | MF5308_P2_M2->MF5308_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5309 | MF5309_P0_M0->MF5309_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5309 | MF5309_P0_M0->MF5309_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5309 | MF5309_P0_M1->MF5309_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5309 | MF5309_P0_M1->MF5309_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5309 | MF5309_P0_M2->MF5309_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9688 | 0.1391 | 1.1492 |
| 5309 | MF5309_P0_M2->MF5309_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5309 | MF5309_P1_M0->MF5309_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5309 | MF5309_P1_M0->MF5309_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 0.9844 | 0.0805 | 1.0838 |
| 5309 | MF5309_P1_M1->MF5309_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5309 | MF5309_P1_M1->MF5309_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 0.9375 | 0.2338 | 1.2634 |
| 5309 | MF5309_P1_M2->MF5309_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5309 | MF5309_P1_M2->MF5309_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0.0156 | 0.0805 | 1.0838 |
| 5309 | MF5309_P2_M0->MF5309_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8438 | 0.4334 | 1.5425 |
| 5309 | MF5309_P2_M0->MF5309_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9062 | 0.3111 | 1.365 |
| 5309 | MF5309_P2_M1->MF5309_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1892 | 1.2083 |
| 5309 | MF5309_P2_M1->MF5309_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1892 | 1.2083 |
| 5309 | MF5309_P2_M2->MF5309_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 0.875 | 0.3768 | 1.4576 |
| 5309 | MF5309_P2_M2->MF5309_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 0.5 | 0.6931 | 2 |
| 5310 | MF5310_P0_M0->MF5310_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5310 | MF5310_P0_M0->MF5310_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5310 | MF5310_P0_M1->MF5310_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0.4061 | 1.5009 |
| 5310 | MF5310_P0_M1->MF5310_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5310 | MF5310_P0_M2->MF5310_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.8594 | 0.4552 | 1.5764 |
| 5310 | MF5310_P0_M2->MF5310_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5310 | MF5310_P1_M0->MF5310_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5310 | MF5310_P1_M0->MF5310_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5310 | MF5310_P1_M1->MF5310_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5310 | MF5310_P1_M1->MF5310_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5310 | MF5310_P1_M2->MF5310_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.2338 | 1.2634 |
| 5310 | MF5310_P1_M2->MF5310_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0.1094 | 0.3452 | 1.4123 |
| 5310 | MF5310_P2_M0->MF5310_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9688 | 0.1391 | 1.1492 |
| 5310 | MF5310_P2_M0->MF5310_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5310 | MF5310_P2_M1->MF5310_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5310 | MF5310_P2_M1->MF5310_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5310 | MF5310_P2_M2->MF5310_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5310 | MF5310_P2_M2->MF5310_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 0.9375 | 0.2338 | 1.2634 |
| 5311 | MF5311_P0_M0->MF5311_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.1892 | 1.2083 |
| 5311 | MF5311_P0_M0->MF5311_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0.2338 | 1.2634 |
| 5311 | MF5311_P0_M1->MF5311_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5311 | MF5311_P0_M1->MF5311_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5311 | MF5311_P0_M2->MF5311_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9844 | 0.0805 | 1.0838 |
| 5311 | MF5311_P0_M2->MF5311_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 1 | 0 | 1 |
| 5311 | MF5311_P1_M0->MF5311_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5311 | MF5311_P1_M0->MF5311_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5311 | MF5311_P1_M1->MF5311_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5311 | MF5311_P1_M1->MF5311_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5311 | MF5311_P1_M2->MF5311_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.1892 | 1.2083 |
| 5311 | MF5311_P1_M2->MF5311_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | iron_reduction | True | 0.6875 | 0.6211 | 1.8609 |
| 5311 | MF5311_P2_M0->MF5311_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9375 | 0.2689 | 1.3086 |
| 5311 | MF5311_P2_M0->MF5311_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9219 | 0.3267 | 1.3864 |
| 5311 | MF5311_P2_M1->MF5311_P2_M2 | core | carbonate | carbonate_weathering | 64 | denitrification | False | 0 | 0.6435 | 1.9031 |
| 5311 | MF5311_P2_M1->MF5311_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | denitrification | False | 0 | 0.6931 | 2 |
| 5311 | MF5311_P2_M2->MF5311_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5311 | MF5311_P2_M2->MF5311_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5312 | MF5312_P0_M0->MF5312_P0_M1 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5312 | MF5312_P0_M0->MF5312_P0_M1 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5312 | MF5312_P0_M1->MF5312_P0_M2 | core | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5312 | MF5312_P0_M1->MF5312_P0_M2 | enhanced | carbonate | carbonate_precipitation | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5312 | MF5312_P0_M2->MF5312_P0_M3 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.5469 | 0.6887 | 1.9912 |
| 5312 | MF5312_P0_M2->MF5312_P0_M3 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.9531 | 0.1892 | 1.2083 |
| 5312 | MF5312_P1_M0->MF5312_P1_M1 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5312 | MF5312_P1_M0->MF5312_P1_M1 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5312 | MF5312_P1_M1->MF5312_P1_M2 | core | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5312 | MF5312_P1_M1->MF5312_P1_M2 | enhanced | sulfate_reduction | sulfate_reduction | 64 | sulfate_reduction | True | 1 | 0 | 1 |
| 5312 | MF5312_P1_M2->MF5312_P1_M3 | core | iron_reduction | iron_reduction | 64 | silicate_exchange | False | 0 | 0.3768 | 1.4576 |
| 5312 | MF5312_P1_M2->MF5312_P1_M3 | enhanced | iron_reduction | iron_reduction | 64 | iron_reduction | True | 0.5469 | 0.6887 | 1.9912 |
| 5312 | MF5312_P2_M0->MF5312_P2_M1 | core | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.7031 | 0.6694 | 1.9531 |
| 5312 | MF5312_P2_M0->MF5312_P2_M1 | enhanced | silicate_exchange | silicate_weathering | 64 | silicate_exchange | True | 0.7188 | 0.6922 | 1.9982 |
| 5312 | MF5312_P2_M1->MF5312_P2_M2 | core | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5312 | MF5312_P2_M1->MF5312_P2_M2 | enhanced | carbonate | carbonate_weathering | 64 | silicate_exchange | False | 0 | 0 | 1 |
| 5312 | MF5312_P2_M2->MF5312_P2_M3 | core | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
| 5312 | MF5312_P2_M2->MF5312_P2_M3 | enhanced | denitrification | denitrification | 64 | denitrification | True | 1 | 0 | 1 |
## Table S5

**Table S5.** The predeclared univariate probability-span conflict
diagnostic (threshold 0.50) under every locked-test condition, reporting
the number of scored edges, the conflict fraction, the mean univariate
probability span, the integrated (HAC) error rate on conflicting versus
concordant edges, and the integrated overconfident-error fraction. No
condition produced a non-zero conflict fraction, including joint
misspecification, illustrating that this single-threshold heuristic did
not detect the misspecification that paired discrimination and calibration
metrics detected clearly (Table 3).

| condition | n_edges | conflict_fraction | mean_univariate_probability_span | integrated_error_rate_conflict | integrated_error_rate_concordant | integrated_overconfident_error_fraction |
| --- | --- | --- | --- | --- | --- | --- |
| native | 827 | 0 | 0.092 |  | 0.1306 | 0 |
| age_permuted | 827 | 0 | 0.1193 |  | 0.1306 | 0.0036 |
| hydraulic_permuted | 827 | 0 | 0.1205 |  | 0.1306 | 0.0351 |
| joint_misspecified | 827 | 0 | 0.123 |  | 0.1306 | 0.0387 |

## Table S6

**Table S6.** Multiplicity-correction robustness check for the six
predeclared case-block contrasts across four metrics (24 tests), added in
response to pre-submission peer review. Each row reports an exact
two-sided sign-flip permutation p-value, computed directly from the
per-case paired differences (n=12, all 4096 sign assignments enumerated
exactly), and a Benjamini-Hochberg adjusted p-value across the full family
of 24 tests at q=0.05. This check does not alter any locked result; it
re-analyses the same per-case values reported in Table 3 with a method
appropriate for a twelve-case sample and corrected for testing 24
hypotheses at once.

| Contrast | Metric | Cases | Mean difference | Exact permutation p (two-sided) | Benjamini-Hochberg adjusted p | Survives q=0.05 |
| --- | --- | --- | --- | --- | --- | --- |
| Native age added | PR-AUC | 12 | -0.0060 | 0.0625 | 0.0703 | No |
| Native age added | Brier | 12 | 0.0001 | 0.1875 | 0.1957 | No |
| Native age added | Log loss | 12 | 0.0000 | 0.8491 | 0.8491 | No |
| Native age added | Entropy | 12 | -0.0006 | 0.0020 | 0.0043 | Yes |
| Native chemistry added | PR-AUC | 12 | 0.4471 | 0.0005 | 0.0012 | Yes |
| Native chemistry added | Brier | 12 | -0.0196 | 0.0005 | 0.0012 | Yes |
| Native chemistry added | Log loss | 12 | -0.0791 | 0.0005 | 0.0012 | Yes |
| Native chemistry added | Entropy | 12 | -0.0827 | 0.0005 | 0.0012 | Yes |
| Native hydraulics added | PR-AUC | 12 | 0.0091 | 0.0645 | 0.0703 | No |
| Native hydraulics added | Brier | 12 | -0.0010 | 0.0005 | 0.0012 | Yes |
| Native hydraulics added | Log loss | 12 | -0.0059 | 0.0005 | 0.0012 | Yes |
| Native hydraulics added | Entropy | 12 | -0.0146 | 0.0005 | 0.0012 | Yes |
| Permuted age added | PR-AUC | 12 | -0.0754 | 0.0410 | 0.0492 | Yes |
| Permuted age added | Brier | 12 | 0.0034 | 0.0181 | 0.0241 | Yes |
| Permuted age added | Log loss | 12 | 0.0105 | 0.0264 | 0.0333 | Yes |
| Permuted age added | Entropy | 12 | -0.0207 | 0.0005 | 0.0012 | Yes |
| Permuted hydraulics added | PR-AUC | 12 | -0.0686 | 0.0151 | 0.0214 | Yes |
| Permuted hydraulics added | Brier | 12 | 0.0110 | 0.0049 | 0.0073 | Yes |
| Permuted hydraulics added | Log loss | 12 | 0.0745 | 0.0029 | 0.0050 | Yes |
| Permuted hydraulics added | Entropy | 12 | -0.0482 | 0.0005 | 0.0012 | Yes |
| Age + hydraulics misspecified | PR-AUC | 12 | -0.1390 | 0.0044 | 0.0070 | Yes |
| Age + hydraulics misspecified | Brier | 12 | 0.0106 | 0.0024 | 0.0045 | Yes |
| Age + hydraulics misspecified | Log loss | 12 | 0.0730 | 0.0024 | 0.0045 | Yes |
| Age + hydraulics misspecified | Entropy | 12 | -0.0706 | 0.0005 | 0.0012 | Yes |
