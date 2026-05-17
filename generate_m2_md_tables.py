import pandas as pd
from pathlib import Path

def convert():
    table_dir = Path("M2/m2_benchmark/tables")
    out_dir = Path("M2/m2_benchmark/tables/Manuscript_Ready")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    mapping = {
        "table1_module_architecture.csv": "Table1_Module_Architecture.md",
        "table2_input_fields.csv": "Table2_Validation_Suite.md",
        "table3_residence_time_options.csv": "Table3_Global_Validation_Performance.md",
        "table4_validation_design_and_results.csv": "Table4_MRT_Accuracy.md",
        "table5_method_comparison.csv": "Table5_MODPATH_Agreement.md",
        "table_s3_reaction_dictionary.csv": "Table6_Discovery_and_PSI.md",
    }
    
    for csv_file, md_file in mapping.items():
        path = table_dir / csv_file
        if path.exists():
            df = pd.read_csv(path)
            md_content = df.to_markdown(index=False)
            (out_dir / md_file).write_text(md_content)
            print(f"Generated {out_dir / md_file}")
        else:
            print(f"Missing {path}")

if __name__ == "__main__":
    convert()
