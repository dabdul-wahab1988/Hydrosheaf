import sys, traceback
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SCRIPT_DIR = REPO_ROOT / "M3" / "m3_age_benchmark" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import run_m3_usgs_benchmark as usgs
import pandas as pd

print("Loading USGS national dataset...", flush=True)
df = usgs.load_usgs_national_dataset()
print(f"Loaded {len(df)} rows.", flush=True)

results = []
for i, (_, row) in enumerate(df.iterrows(), start=1):
    if i == 1 or i % 100 == 0:
        print(f"  Processing row {i}/{len(df)}...", flush=True)
    try:
        res = usgs.fit_benchmark_row(row, age_steps=usgs.M3_DEFAULT_AGE_STEPS)
        results.append(res)
    except Exception as exc:
        print(f"ERROR at row {i}: {exc}", flush=True)
        traceback.print_exc()
        break

out_df = pd.DataFrame(results)
out_path = SCRIPT_DIR.parent / "results" / "m3_usgs_benchmark_results.csv"
out_df.to_csv(out_path, index=False)
print(f"Wrote {len(out_df)} rows to {out_path}", flush=True)
