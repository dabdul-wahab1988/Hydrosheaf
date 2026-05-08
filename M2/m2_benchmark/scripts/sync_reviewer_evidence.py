from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = BENCHMARK_ROOT / "tables"
DOC_DIR = BENCHMARK_ROOT / "docs"


def _read_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(BENCHMARK_ROOT / path)


def _scalar(df: pd.DataFrame, column: str, row: int = 0) -> float:
    return float(df.iloc[row][column])


def _first_status(df: pd.DataFrame, group: str = "all") -> pd.Series:
    if "group" in df.columns:
        matched = df[df["group"].astype(str) == group]
        if not matched.empty:
            return matched.iloc[0]
    return df.iloc[0]


def build_table4() -> pd.DataFrame:
    transport = _read_csv("results/transport_recovery.csv")
    reactions = _read_csv("results/reaction_recovery.csv")
    e1 = _first_status(_read_csv("external/usgs_age/results/usgs_age_validation_summary.csv"))
    e1b = _read_csv("external/usgs_age/results/e1b_joint_lpm_validation.csv")
    e2 = _read_csv("external/modpath/results/modpath_topology_summary.csv").iloc[0]
    e2b = _read_csv("external/modpath/results/modpath_geochemical_replay_summary.csv").iloc[0]
    e3 = _read_csv("results/phreeqc_forward_validation_summary.csv")
    e4 = _read_csv("external/northern_ghana/results/northern_ghana_validation_summary.csv").iloc[0]
    e5 = _first_status(_read_csv("external/dgmeta/results/dgmeta_validation_summary.csv"))
    e6 = _read_csv("external/e6_nonlinear/results/e6_nonlinear_validation.csv")

    e3_median = e3[e3.iloc[:, 0] == "median"].iloc[0]
    e3_mean = e3[e3.iloc[:, 0] == "mean"].iloc[0]
    e3_sum = e3[e3.iloc[:, 0] == "sum"].iloc[0]

    e1b_match = float(e1b["model_match"].mean())
    e1b_rmse = float(((e1b["log10_error"]) ** 2).mean() ** 0.5)

    rows = [
        [
            "Synthetic aquifer benchmark",
            "locked M2 ground truth",
            "transport and reaction extents",
            "median transport absolute error={:.3f}; median reaction error={:.3f}".format(
                float(transport["absolute_error"].median()),
                float(reactions["absolute_error_mmolL"].median()),
            ),
            "direct recovery of known processes",
            "SyntheticDataGuide.docx",
            "completed",
        ],
        [
            "Public tracer-age validation",
            "USGS public-supply aquifer groundwater-age data release, DOI 10.5066/P9W7T0DN",
            "published mean age and age fractions",
            "E1 all samples n={}; log10 age RMSE={:.3f}; median bias={:.3f}; CI coverage={:.3f}; E1b joint-LPM subset n={}; RMSE={:.3f}; family match={:.2%}".format(
                int(e1["n_samples"]),
                float(e1["log10_age_rmse"]),
                float(e1["median_log10_bias"]),
                float(e1["ci_coverage_fraction"]),
                len(e1b),
                e1b_rmse,
                e1b_match,
            ),
            "screening-level age-magnitude agreement; weak model-family reproduction is reported explicitly",
            "Jurgens et al. USGS data release 10.5066/P9W7T0DN",
            "completed with interpretation guardrail",
        ],
        [
            "MODFLOW/MODPATH topology validation",
            "USGS Savage Municipal Water-Supply Well MODFLOW-2005/MODPATH5 archive, DOI 10.5066/F7J102FK",
            "directed endpoint/pathline graph priors",
            "TP={}; FP={}; FN={}; F1={:.2f}; direction agreement={:.2f}".format(
                int(e2["true_positive_edges"]),
                int(e2["false_positive_edges"]),
                int(e2["false_negative_edges"]),
                float(e2["edge_f1"]),
                float(e2["direction_agreement_rate"]),
            ),
            "validates MODPATH-to-graph conversion, not geochemical inference",
            "Harte USGS data release 10.5066/F7J102FK",
            "completed",
        ],
        [
            "MODPATH-conditioned geochemical replay",
            "USGS Savage MODPATH topology with controlled injected chemistry",
            "transport gamma and reaction identity/extent on real MODPATH edges",
            "n={}; gamma MAE={:.3f}; reaction F1={:.2f}; extent MAE={:.3f}".format(
                int(e2b["n_edges"]),
                float(e2b["gamma_mae"]),
                float(e2b["reaction_f1"]),
                float(e2b["mean_extent_mae"]),
            ),
            "semi-synthetic integration test because the archive lacks paired chemistry at graph nodes",
            "E2b replay script and USGS Savage archive",
            "completed with scope guardrail",
        ],
        [
            "Live PHREEQC forward validation",
            "USGS PHREEQC version 3 examples and databases, DOI 10.3133/tm6A43",
            "major-ion residuals and saturation-index diagnostics",
            "n={}; median RMSE={:.2f} mmol/L; mean RMSE={:.2f} mmol/L; median NSE={:.2f}; mean NSE={:.2f}".format(
                int(e3_sum["phreeqc_ok"]),
                float(e3_median["rmse_mmolL"]),
                float(e3_mean["rmse_mmolL"]),
                float(e3_median["nse"]),
                float(e3_mean["nse"]),
            ),
            "supports thermodynamic feasibility with moderate forward predictive fit",
            "Parkhurst and Appelo 2013 USGS TM 6-A43",
            "completed with moderate-fit caveat",
        ],
        [
            "Nonlinear PHREEQC synthetic validation",
            "PHREEQC-generated nonlinear synthetic ground truth",
            "forward-validation NSE under thermodynamically consistent synthetic conditions",
            "n={}; median NSE={:.2f}; mean NSE={:.2f}".format(
                len(e6),
                float(e6["nse"].median()),
                float(e6["nse"].mean()),
            ),
            "shows improved fit under internally consistent nonlinear PHREEQC truth; does not prove perfect recovery",
            "E6 nonlinear synthetic validation script",
            "completed with no near-perfect claim",
        ],
        [
            "Data-limited pilot scenario",
            "Corrected Northern Ghana aquifer workbook",
            "field hydrochemistry, generated wet/dry graph edges, and process-fit diagnostics under sparse inputs",
            "n_wells={}; quantitative samples={}; candidate edge fits={}; top-20 minimum chemistry R2={:.3f}; median chemistry R2={:.3f}".format(
                int(e4["n_wells"]),
                int(e4["fit_samples_quantitative"]),
                int(e4["candidate_edge_fits"]),
                float(e4["top20_min_chemistry_r2"]),
                float(e4["median_chemistry_r2"]),
            ),
            "replaces the local Manu pilot with corrected field-hydrochemistry workflow evidence; generated graph is not independent process truth",
            "data/NorthenGhana/NorthernGhana.xlsx; public DOI/source URL not embedded in workbook",
            "completed field-hydrochemistry demonstration with generated-graph guardrail",
        ],
        [
            "DGMETA dissolved-gas validation",
            "USGS DGMETA example workbooks, DOI 10.5066/P9NQ1RFY",
            "recharge temperature and excess-air estimates",
            "n={}; temperature MAE={:.2f} C; excess-air MAE={:.2f} cm3/kg; large residual fraction={:.3f}".format(
                int(e5["n_rows"]),
                float(e5["temperature_mae_c"]),
                float(e5["excess_air_mae_cm3kg"]),
                float(e5["large_residual_fraction"]),
            ),
            "validates dissolved-gas fitting mechanics; exact DGMETA parity still needs macro-level verification",
            "USGS DGMETA examples 10.5066/P9NQ1RFY",
            "completed with DGMETA-equivalence guardrail",
        ],
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "benchmark",
            "data_source",
            "target_variable",
            "performance_metric",
            "expected_evidence",
            "key_reference",
            "m2_status",
        ],
    )


def build_table_s4() -> pd.DataFrame:
    rows = [
        [
            "M2 locked synthetic benchmark",
            "all required chemistry, tracer, topology, and ground-truth reaction fields",
            "local generated package",
            "synthetic recovery",
            "completed",
            "SyntheticDataGuide.docx",
        ],
        [
            "USGS public-supply groundwater-age data release",
            "site metadata, environmental tracers, dissolved gases, LPM outputs, age fractions",
            "local cached zip plus https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004",
            "external residence-time validation",
            "completed with screening-level interpretation",
            "10.5066/P9W7T0DN",
        ],
        [
            "USGS Savage MODFLOW/MODPATH model archive",
            "MODFLOW-2005 inputs/outputs, MODPATH5 endpoints/pathlines, MOC3D transport files",
            "local cached archive plus https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute",
            "external topology validation and MODPATH-conditioned geochemical replay",
            "completed",
            "10.5066/F7J102FK",
        ],
        [
            "USGS PHREEQC version 3 examples and databases",
            "speciation, batch reaction, inverse modelling, transport examples, thermodynamic databases",
            "https://www.usgs.gov/software/phreeqc-version-3",
            "live forward reaction validation",
            "completed with moderate-fit caveat",
            "10.3133/tm6A43",
        ],
        [
            "USGS DGMETA example workbooks",
            "noble-gas model outputs, recharge temperature, excess-air parameters",
            "USGS DGMETA example workbooks and cached parsed outputs",
            "dissolved-gas correction validation",
            "completed with DGMETA-equivalence guardrail",
            "10.5066/P9NQ1RFY",
        ],
        [
            "Corrected Northern Ghana aquifer workbook",
            "160 boreholes, 320 dry/wet hydrochemical samples, major ions, isotopes, Sr, SiO2, coordinates, depth, static water level, and distance covariates; no supplied edge or saturation-index sheet",
            "data/NorthenGhana/NorthernGhana.xlsx plus optional data/NorthenGhana/SI.pdf",
            "field hydrochemistry and generated sparse graph-edge demonstration",
            "completed with field-data and generated-graph guardrails; public DOI/source URL not embedded",
            "public source citation to be supplied in manuscript if available",
        ],
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "resource",
            "available_variables",
            "access_pathway",
            "hydrosheaf_validation_workflow",
            "status",
            "identifier",
        ],
    )


def build_workplan() -> pd.DataFrame:
    rows = [
        [
            "E1_public_tracer_age",
            "3.3 and Figure 4A",
            "USGS groundwater-age public-supply aquifer data",
            "10.5066/P9W7T0DN",
            "https://www.usgs.gov/data/data-distribution-groundwater-age-aquifers-used-public-supply-continental-united-states-2004",
            "Compare Hydrosheaf residence-time estimates with published mean ages and age fractions",
            "external/usgs_age/results/usgs_age_validation.csv; external/usgs_age/results/e1b_joint_lpm_validation.csv; figures/figure4a_public_tracer_age_agreement.png",
            "completed_with_guardrail",
            "E1 n=689, RMSE=1.230; E1b n=69, RMSE=0.956, model-family match=10.14%; interpret as screening-level age agreement.",
        ],
        [
            "E2_modpath_topology",
            "3.4 and Figure S2",
            "USGS Savage Municipal Water-Supply Well MODFLOW/MODPATH archive",
            "10.5066/F7J102FK",
            "https://www.usgs.gov/data/modflow-2005-modpath-and-moc3d-used-groundwater-flow-simulation-pathlines-analysis-and-solute",
            "Convert MODPATH endpoints/pathlines into graph priors and compare inferred edges with particle-tracking edges",
            "external/modpath/results/modpath_graph_priors.csv; external/modpath/results/modpath_topology_agreement.csv; figures/figure_s2_modpath_to_graph_prior_real_archive.png",
            "completed",
            "TP=174; FP=0; FN=0; F1=1.00; direction agreement=1.00; topology only, not geochemical inference.",
        ],
        [
            "E2b_modpath_geochemical_replay",
            "3.4 and Figure S2B",
            "USGS Savage MODPATH graph with injected known-truth chemistry",
            "10.5066/F7J102FK",
            "local MODPATH archive plus controlled replay script",
            "Inject known transport/reaction signals along real MODPATH edges and test Hydrosheaf recovery",
            "external/modpath/results/modpath_geochemical_replay.csv; figures/figure_s2b_modpath_conditioned_geochemical_replay.png",
            "completed_with_guardrail",
            "n=80; gamma MAE=0.005; reaction F1=1.00; extent MAE=0.002; semi-synthetic integration validation.",
        ],
        [
            "E3_live_phreeqc_forward",
            "3.5 and Figure S3",
            "USGS PHREEQC version 3 examples and databases",
            "10.3133/tm6A43",
            "https://pubs.usgs.gov/tm/06/a43/",
            "Run live PHREEQC forward validation for inferred reaction pathways",
            "results/phreeqc_forward_validation.csv; figures/figure_s3_phreeqc_forward_diagnostics_live.png",
            "completed_with_moderate_fit_caveat",
            "n=280; median RMSE=1.35 mmol/L; median NSE=0.59; mean NSE=0.47.",
        ],
        [
            "E4c_northern_ghana_field_hydrochemistry",
            "3.6",
            "Corrected Northern Ghana aquifer workbook",
            "public source citation to be supplied in manuscript if available",
            "data/NorthenGhana/NorthernGhana.xlsx",
            "Run Hydrosheaf on charge-balance-screened wet/dry field hydrochemistry using a generated coordinate/head-proxy graph",
            "external/northern_ghana/results/northern_ghana_samples.csv; external/northern_ghana/results/northern_ghana_edge_results.csv; figures/figure_s8_e4_sparse_pilot_network.png",
            "completed_field_demonstration_with_generated_graph_guardrail",
            "Replaces local Manu pilot. Field-hydrochemistry workflow evidence only; generated graph is not independent process truth, and there is no independent tracer-age, MODPATH, or PHREEQC truth.",
        ],
        [
            "E5_dgmeta_dissolved_gas",
            "Appendix and Table 4",
            "USGS DGMETA example workbooks",
            "10.5066/P9NQ1RFY",
            "https://doi.org/10.5066/P9NQ1RFY",
            "Validate dissolved-gas correction mechanics against DGMETA reference outputs",
            "external/dgmeta/results/dgmeta_hydrosheaf_comparison.csv; external/dgmeta/results/dgmeta_validation_summary.csv",
            "completed_with_guardrail",
            "n=409; temperature MAE=0.51 C; excess-air MAE=2.44 cm3/kg; exact DGMETA parity still requires macro-level verification.",
        ],
        [
            "E6_nonlinear_phreeqc_synthetic",
            "Supplementary nonlinear validation",
            "PHREEQC-generated nonlinear synthetic truth",
            "local generated package",
            "external/e6_nonlinear/results/e6_nonlinear_validation.csv",
            "Test forward fit when truth is generated by the same nonlinear thermodynamic class",
            "external/e6_nonlinear/results/e6_nonlinear_validation.csv",
            "completed_with_no_perfect_recovery_claim",
            "n=50; median NSE=0.88; mean NSE=0.85.",
        ],
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "validation_tier",
            "required_by_m2_section",
            "primary_dataset",
            "source_or_doi",
            "source_url",
            "hydrosheaf_task",
            "required_outputs",
            "status",
            "note",
        ],
    )


def write_docs(table4: pd.DataFrame) -> None:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    summary = [
        "# M2 Hydrosheaf Benchmark Results Summary",
        "",
        f"Reviewer evidence synchronized: {timestamp}",
        "",
        "This summary is intentionally conservative. The current artefacts support a defensible M2 validation package only when the stated guardrails are retained in captions, results text, and discussion.",
        "",
        "## Table 4 Snapshot",
        "",
        table4.to_markdown(index=False),
        "",
        "## Interpretation Guardrails",
        "",
        "- E1/E1b supports screening-level residence-time agreement, not full TracerLPM model-family equivalence.",
        "- E2 validates MODPATH-to-graph conversion; E2b is a semi-synthetic replay on real MODPATH topology because paired chemistry is absent.",
        "- E3 supports thermodynamic feasibility with moderate forward predictive fit; do not describe NSE around 0.59 as high.",
        "- E4c replaces the local Manu pilot with the corrected Northern Ghana workbook as field-hydrochemistry/data-limited workflow evidence; the graph is generated from coordinates and head proxy, not independent process truth.",
        "- E5 validates dissolved-gas fitting mechanics against DGMETA examples; exact DGMETA parity still needs workbook macro-level verification.",
        "- E6 shows strong but not near-perfect nonlinear synthetic performance.",
        "",
        "## Main Outputs",
        "",
        "- `docs/02_figures.md`",
        "- `docs/03_tables.md`",
        "- `tables/table4_validation_design_and_results.csv`",
        "- `tables/table_s4_benchmark_dataset_inventory.csv`",
        "- `figures/figure1_architecture_workflow.png`",
        "- `figures/figure4_external_validation_composite.png`",
        "- `figures/figure_s2_modpath_to_graph_prior_real_archive.png`",
        "- `figures/figure_s2b_modpath_conditioned_geochemical_replay.png`",
        "- `figures/figure_s3_phreeqc_forward_diagnostics_live.png`",
        "- `figures/figure_s5_e6_nonlinear_validation.png`",
        "- `figures/figure_s6_dgmeta_diagnostics.png`",
        "- `figures/figure_s7_e1_tracer_age_residuals.png`",
        "- `figures/figure_s8_e4_sparse_pilot_network.png`",
        "",
    ]
    (DOC_DIR / "m2_results_summary.md").write_text("\n".join(summary), encoding="utf-8")

    figures = [
        "# 02 Figures Evidence Register",
        "",
        f"Updated: {timestamp}",
        "",
        "| figure | artefact | evidence source | allowed manuscript claim | required guardrail | status |",
        "|:--|:--|:--|:--|:--|:--|",
        "| Figure 1 | `figures/figure1_architecture_workflow.png` | `tables/table1_module_architecture.csv` | Hydrosheaf integrates preprocessing, graph priors, LPM age estimation, sparse inverse reaction fitting, PHREEQC diagnostics, and uncertainty outputs. | Caption must not imply every optional external module is required for minimum operation. | revised |",
        "| Figure 2 | `figures/figure2_process_network.png` | `config/ground_truth.yaml`; synthetic node/edge tables | Shows the locked synthetic process network used for M2 recovery tests. | Treat as synthetic benchmark design, not a field aquifer map. | revised |",
        "| Figure 3 | `figures/figure3_synthetic_benchmark_recovery.png` | `results/transport_recovery.csv`; `results/reaction_recovery.csv`; `results/missing_data_sensitivity.csv`; `results/topology_robustness.csv` | Synthetic recovery is strong for locked transport/reaction cases and degrades under missing-data/topology perturbation scenarios. | Caption must match actual panels: transport, reaction, scenario bias, and candidate-edge robustness. | usable after caption correction |",
        "| Figure 4 | `figures/figure4_external_validation_composite.png` | E1, E2, E3, and E5 result files | Integrates the main external validation evidence: USGS tracer-age comparison, MODPATH topology agreement, live PHREEQC diagnostics, and DGMETA temperature agreement. | Must retain layer-specific limitations; this is not a claim of full reactive-transport equivalence. | added |",
        "| Figure 4A | `figures/figure4a_public_tracer_age_agreement.png` | `external/usgs_age/results/usgs_age_validation.csv`; `external/usgs_age/results/e1b_joint_lpm_validation.csv` | Public tracer-age comparison supports screening-level age-magnitude validation. | Use only if a standalone age-validation figure is needed; do not claim TracerLPM equivalence. | supplementary candidate |",
        "| Figure 5 | `figures/figure5_sensitivity_uncertainty.png` | `results/sensitivity_uncertainty_summary.csv` | Monte Carlo and degraded-input diagnostics quantify uncertainty sensitivity. | Claims must stay tied to the synthetic benchmark and scaled contributions. | revised |",
        "| Figure S1 | `figures/figure_s1_sparse_fitting_algorithm.png` | method implementation and reaction dictionary | Summarizes the sparse fitting workflow. | Keep as schematic, not validation evidence. | revised |",
        "| Figure S2 | `figures/figure_s2_modpath_to_graph_prior_real_archive.png` | `external/modpath/results/modpath_topology_summary.csv` | MODPATH endpoint/pathline conversion preserves source-receptor topology. | State explicitly that this validates topology only. | usable with visual tightening if journal requires |",
        "| Figure S2B | `figures/figure_s2b_modpath_conditioned_geochemical_replay.png` | `external/modpath/results/modpath_geochemical_replay_summary.csv` | Hydrosheaf recovers injected chemistry on real MODPATH topology. | Label as semi-synthetic replay, not a paired field validation. | usable with guarded interpretation |",
        "| Figure S3 | `figures/figure_s3_phreeqc_forward_diagnostics_live.png` | `results/phreeqc_forward_validation.csv`; `results/phreeqc_forward_validation_summary.csv` | Live PHREEQC diagnostics support thermodynamic feasibility with moderate forward fit. | Caption must describe RMSE/NSE/SI diagnostics, not time-series concentration evolution. | usable after caption correction |",
        "| Figure S4 | `figures/figure_s4_residence_time_network_update.png` | `results/age_inference_validation.csv`; `results/age_network_consistency.csv` | Network priors modify residence-time estimates relative to single-node age estimates in the synthetic benchmark. | Synthetic/internal residence-time evidence only. | added |",
        "| Figure S5 | `figures/figure_s5_e6_nonlinear_validation.png` | `external/e6_nonlinear/results/e6_nonlinear_validation.csv` | E6 supports nonlinear PHREEQC-consistent synthetic performance with median NSE around 0.88. | Supplementary consistency test; does not replace E3 and does not prove perfect recovery. | added |",
        "| Figure S6 | `figures/figure_s6_dgmeta_diagnostics.png` | `external/dgmeta/results/dgmeta_hydrosheaf_comparison.csv`; `external/dgmeta/results/dgmeta_validation_summary.csv` | DGMETA diagnostics show recharge-temperature and excess-air agreement patterns. | Exact DGMETA parity still requires macro/workbook option verification. | added |",
        "| Figure S7 | `figures/figure_s7_e1_tracer_age_residuals.png` | `external/usgs_age/results/usgs_age_validation.csv` | Tracer-age residuals are stratified by tracer support and age class. | Use to disclose performance heterogeneity, not to overstate age-model equivalence. | added |",
        "| Figure S8 | `figures/figure_s8_e4_sparse_pilot_network.png` | `external/northern_ghana/results/northern_ghana_edge_results.csv` | E4c shows corrected Northern Ghana field-hydrochemistry fits under wet/dry sparse-data conditions. | Generated coordinate/head-proxy graph only; no independent tracer-age, MODPATH, process-truth graph, or PHREEQC truth. | revised |",
        "",
    ]
    (DOC_DIR / "02_figures.md").write_text("\n".join(figures), encoding="utf-8")

    tables = [
        "# 03 Tables Evidence Register",
        "",
        f"Updated: {timestamp}",
        "",
        "| table | artefact | evidence source | allowed manuscript claim | required guardrail | status |",
        "|:--|:--|:--|:--|:--|:--|",
        "| Table 1 | `tables/table1_module_architecture.csv` | Hydrosheaf module design | Defines the package modules and validation roles. | Keep descriptive, not performance evidence. | usable |",
        "| Table 2 | `tables/table2_input_fields.csv` | M2 input specification | Separates minimum viable from enhanced inputs. | Maintain unit conventions and optional tracer wording. | usable |",
        "| Table 3 | `tables/table3_residence_time_options.csv` | Hydrosheaf tracer/LPM implementation | Documents tracer/LPM support and limits. | Must say full LPM-family fitting is a validation layer with identifiability limits. | usable with cautious wording |",
        "| Table 4 | `tables/table4_validation_design_and_results.csv` | E1-E6 reports and result CSVs | Central reviewer-facing validation table. | Do not remove caveats for E1b, E2b, E3, E4, E5, or E6. | revised |",
        "| Table 5 | `tables/table5_method_comparison.csv` | method comparison | Positions Hydrosheaf relative to PHREEQC/NETPATH, MODPATH, and TracerLPM. | Avoid claiming replacement of specialist tools. | usable |",
        "| Table S2 | `tables/table_s2_tracer_properties.csv` | tracer property literature and implementation | Summarizes tracer properties and uncertainty sources. | Keep as background, not validation result. | usable |",
        "| Table S3 | `tables/table_s3_reaction_dictionary.csv` | Hydrosheaf reaction dictionary | Lists available process labels and stoichiometry. | Ensure units remain mmol/L extents in text. | usable |",
        "| Table S4 | `tables/table_s4_benchmark_dataset_inventory.csv` | local and public validation datasets | Provides dataset provenance and current status. | Must stay synchronized with Table 4 and workplan. | revised |",
        "| External validation workplan | `tables/external_validation_workplan.csv` | E1-E6 reports | Records status and output paths for each validation tier. | Treat E4c as field-hydrochemistry demonstration, E2b as semi-synthetic, and E3 as moderate fit. | revised |",
        "",
    ]
    (DOC_DIR / "03_tables.md").write_text("\n".join(tables), encoding="utf-8")


def write_manifest() -> None:
    manifest_paths = []
    for folder in ["config", "data", "results", "tables", "figures", "docs", "scripts", "external"]:
        root = BENCHMARK_ROOT / folder
        if not root.exists():
            continue
        for path in sorted(root.rglob("*")):
            if not path.is_file():
                continue
            if "__pycache__" in path.parts or path.suffix == ".pyc":
                continue
            manifest_paths.append(path.relative_to(BENCHMARK_ROOT).as_posix())
    (BENCHMARK_ROOT / "MANIFEST.md").write_text(
        "# M2 Benchmark Manifest\n\n" + "\n".join(f"- `{path}`" for path in manifest_paths) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    table4 = build_table4()
    table4.to_csv(TABLE_DIR / "table4_validation_design_and_results.csv", index=False)
    build_table_s4().to_csv(TABLE_DIR / "table_s4_benchmark_dataset_inventory.csv", index=False)
    build_workplan().to_csv(TABLE_DIR / "external_validation_workplan.csv", index=False)
    write_docs(table4)
    write_manifest()
    print("Synchronized reviewer evidence tables, registers, summary, and manifest.")


if __name__ == "__main__":
    main()
