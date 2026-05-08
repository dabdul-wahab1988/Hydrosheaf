import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os
import time
import warnings

# Suppress warnings
warnings.filterwarnings('ignore')

# Add Hydrosheaf to path
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.nuclear.joint_lpm import fit_lpm_models

BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = BENCHMARK_ROOT / "external" / "usgs_age" / "input" / "DataForNationalGroundwaterAge_1_1"
RESULT_DIR = BENCHMARK_ROOT / "external" / "usgs_age" / "results"

def run_e1b_validation():
    print("Starting E1b Joint-LPM Validation...")
    
    t4_path = INPUT_DIR / "Table_4_LPM_ModOut.txt"
    if not t4_path.exists():
        t4_path = BENCHMARK_ROOT / "external" / "usgs_age" / "input" / "Table_4_LPM_ModOut.txt"
        if not t4_path.exists():
            print(f"Error: Table_4_LPM_ModOut.txt not found.")
            return
            
    df4 = pd.read_csv(t4_path, sep='\t', low_memory=False, na_values=['na', 'NaN', ''])
    
    # Filter for valid LPM fits
    df4 = df4[df4['LPM_Name'].notna()]
    
    # Subsample for performance (grid search is expensive)
    if len(df4) > 100:
        df4 = df4.sample(100, random_state=42)
        print(f"Subsampled to 100 rows for performance.")
    
    print(f"Total reference fits to evaluate: {len(df4)}")
    
    results = []
    
    # Optimized fit parameters for validation speed
    fit_kwargs = {
        'age_steps': 40,
        'models': ["PFM", "EM", "DM", "EPM", "BMM-DM-DM"]
    }
    
    start_time = time.time()
    for i, (_, row) in enumerate(df4.iterrows()):
        if (i+1) % 10 == 0 or i == 0:
            elapsed = time.time() - start_time
            rate = (i+1) / elapsed if elapsed > 0 else 0
            print(f"  Processed {i+1}/{len(df4)} rows... ({rate:.1f} rows/s)")
            
        sample_id = row['SampleID']
        sample_date = str(row['SampleDate'])
        
        try:
            if '/' in sample_date:
                sample_year = float(sample_date.split('/')[-1])
            else:
                sample_year = float(pd.to_datetime(sample_date).year)
        except:
            sample_year = 2010.0
            
        obs_dict = {}
        tracers_found = []
        for j in range(1, 11):
            name = row.get(f'LPM_Tracer_Name_{j:02d}')
            if pd.isna(name) or name == 'na': continue
            val = row.get(f'LPM_Meas_Tracer_{j:02d}')
            err = row.get(f'LPM_Meas_Tracer_{j:02d}_Err')
            if pd.isna(val): continue
            
            hs_key = None
            if name == '3H': hs_key = 'tritium_TU'
            elif name == '3He(trit)': hs_key = 'he3_trit_TU'
            elif name == 'SF6': hs_key = 'sf6_pptv'
            elif name == '14C': hs_key = 'c14_pmc'
            elif name == 'CFC-11': hs_key = 'cfc11_pptv'
            elif name == 'CFC-12': hs_key = 'cfc12_pptv'
            elif name == 'CFC-113': hs_key = 'cfc113_pptv'
            
            if hs_key:
                obs_dict[hs_key] = float(val)
                tracers_found.append(name)
                if not pd.isna(err):
                    obs_dict[hs_key.replace('_', '_sigma_')] = float(err)

        if len(tracers_found) < 1:
            continue
            
        try:
            fits = fit_lpm_models(obs_dict, sample_year=sample_year, **fit_kwargs)
        except Exception:
            continue
            
        if not fits or not fits[0].converged:
            continue
            
        best = fits[0]
        
        if best.model.startswith('BMM'):
             age1 = best.parameters.get('mean_age_1_years', 0.0)
             age2 = best.parameters.get('mean_age_2_years', 0.0)
             f = best.parameters.get('binary_fraction', 0.5)
             hs_age = f * age1 + (1.0 - f) * age2
        else:
             hs_age = best.parameters.get('mean_age_years', 0.0)
             
        ref_model = str(row['LPM_Name'])
        ref_age = row['Rpt_TotAge_yrs']
        if pd.isna(ref_age):
             ref_age = row['LPM_Age_C1_yrs']
             
        try:
            ref_age = float(ref_age)
        except:
            continue
            
        results.append({
            "sample_id": sample_id,
            "ref_model": ref_model,
            "ref_age": ref_age,
            "hs_model": best.model,
            "hs_age": hs_age,
            "hs_rmse": best.rmse_standardized,
            "n_tracers": best.n_tracers,
            "tracers": ",".join(tracers_found)
        })
        
    if not results:
        print("No validation results generated.")
        return
        
    df_res = pd.DataFrame(results)
    df_res['log10_error'] = np.log10(df_res['hs_age'] + 1) - np.log10(df_res['ref_age'] + 1)
    rmse = np.sqrt(np.mean(df_res['log10_error']**2))
    bias = np.median(df_res['log10_error'])
    df_res['model_match'] = df_res['ref_model'].str.upper() == df_res['hs_model'].str.upper()
    match_rate = df_res['model_match'].mean()
    
    print(f"\nE1b Results Summary (n={len(df_res)}):")
    print(f"  Log10 Age RMSE: {rmse:.4f}")
    print(f"  Median Log10 Bias: {bias:.4f}")
    print(f"  Model Family Match Rate: {match_rate:.2%}")
    
    output_path = RESULT_DIR / "e1b_joint_lpm_validation.csv"
    df_res.to_csv(output_path, index=False)
    
    report_path = RESULT_DIR / "e1b_joint_lpm_report.md"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("# E1b Joint-LPM Validation Report\n\n")
        f.write(f"Run timestamp: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n\n")
        f.write("## Metrics\n\n")
        f.write(f"| Metric | Value |\n")
        f.write(f"| :--- | :--- |\n")
        f.write(f"| Samples Evaluated | {len(df_res)} |\n")
        f.write(f"| Log10 Age RMSE | {rmse:.4f} |\n")
        f.write(f"| Median Log10 Bias | {bias:.4f} |\n")
        f.write(f"| Model Family Match Rate | {match_rate:.2%} |\n\n")
        f.write("## Accuracy by Model Family\n\n")
        grouped = df_res.groupby('ref_model')['log10_error'].agg(['count', 'mean', 'std'])
        f.write(grouped.to_markdown() + "\n\n")
        f.write("## Interpretation\n")
        f.write("The Joint-LPM validation tests Hydrosheaf's ability to jointly fit multiple tracers ")
        f.write("using standard lumped parameter model families, compared to official USGS TracerLPM fits.\n")

    print(f"Validation complete. Results saved to {output_path}")

if __name__ == "__main__":
    run_e1b_validation()
