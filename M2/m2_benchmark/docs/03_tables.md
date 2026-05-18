# 03 Tables Evidence Register

Updated: 2026-05-17

| table | artefact | evidence source | allowed manuscript claim | required guardrail | status |
|:--|:--|:--|:--|:--|:--|
| Table 1 | `tables/table1_module_architecture.csv` | Hydrosheaf module design | Defines the package modules and validation roles. | Keep descriptive, not performance evidence. | usable |
| Table 2 | `tables/table2_input_fields.csv` | M2 input specification | Separates minimum viable from enhanced inputs. | Maintain unit conventions and optional tracer wording. | usable |
| Table 3 | `tables/table3_residence_time_options.csv` | Hydrosheaf tracer/LPM implementation | Documents tracer/LPM support and limits. | Must say full LPM-family fitting is a validation layer with identifiability limits. | usable with cautious wording |
| Table 4 | `tables/table4_validation_design_and_results.csv` | M2 synthetic outputs, M3 public-age output, MODPATH summary, PHREEQC proxy, Ghana field outputs | Central reviewer-facing validation-status table. | Keep statuses tier-specific; do not collapse proxy, synthetic, topology-only, and field-demonstration evidence into one validation claim. | revised |
| Table 5 | `tables/table5_method_comparison.csv` | method comparison | Positions Hydrosheaf relative to PHREEQC/NETPATH, MODPATH, and TracerLPM. | Avoid claiming replacement of specialist tools. | usable |
| Table S2 | `tables/table_s2_tracer_properties.csv` | tracer property literature and implementation | Summarizes tracer properties and uncertainty sources. | Keep as background, not validation result. | usable |
| Table S3 | `tables/table_s3_reaction_dictionary.csv` | Hydrosheaf reaction dictionary | Lists available process labels and stoichiometry. | Ensure units remain mmol/L extents in text. | usable |
| Table S4 | `tables/table_s4_benchmark_dataset_inventory.csv` | local and public validation datasets | Provides dataset provenance and current status. | Must stay synchronized with Table 4 and workplan. | revised |
| External validation plan | `docs/external_validation_plan.md` | E1-E4c plan and result-path expectations | Records status and output paths for each validation tier. | Treat public age as M3 screening evidence, MODPATH as topology-only, PHREEQC as proxy unless live rerun exists, and Ghana as field demonstration. | revised |
