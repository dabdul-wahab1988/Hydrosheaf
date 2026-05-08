import pandas as pd
from pathlib import Path

BENCHMARK_ROOT = Path(__file__).resolve().parent.parent

def clean_workplan():
    path = BENCHMARK_ROOT / "tables" / "external_validation_workplan.csv"
    if not path.exists(): return
    df = pd.read_csv(path)
    
    # Standard columns
    cols = [
        'validation_tier', 'required_by_m2_section', 'primary_dataset', 
        'source_or_doi', 'source_url', 'hydrosheaf_task', 
        'required_outputs', 'status', 'note'
    ]
    
    # Deduplicate by validation_tier, keeping the last one (most recent results)
    df = df.drop_duplicates(subset=['validation_tier'], keep='last')
    
    # Ensure all columns exist
    for c in cols:
        if c not in df.columns:
            df[c] = ""
            
    # Reorder
    df = df[cols]
    
    df.to_csv(path, index=False)
    print(f"Cleaned {path}")

def clean_table4():
    path = BENCHMARK_ROOT / "tables" / "table4_validation_design_and_results.csv"
    if not path.exists(): return
    df = pd.read_csv(path)
    
    # Deduplicate by benchmark name
    df = df.drop_duplicates(subset=['benchmark'], keep='last')
    
    df.to_csv(path, index=False)
    print(f"Cleaned {path}")

if __name__ == "__main__":
    clean_workplan()
    clean_table4()
