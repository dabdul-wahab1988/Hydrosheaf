import pandas as pd
import numpy as np
import sys
from pathlib import Path
import math

# Add hydrosheaf to path
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from hydrosheaf.nuclear.multi_tracer import build_atmospheric_tracer_input
import run_m3_usgs_benchmark as usgs

def _apparent_3h3he_age(tritium: float, he3_trit: float) -> float:
    try:
        tritium = float(tritium)
        he3_trit = float(he3_trit)
    except (TypeError, ValueError):
        return float("nan")
    if not math.isfinite(tritium) or not math.isfinite(he3_trit) or tritium <= 0 or he3_trit < 0:
        return float("nan")
    lambda_tritium = math.log(2.0) / 12.32
    return math.log1p(he3_trit / tritium) / lambda_tritium

def main():
    print("Loading datasets...")
    df = usgs.load_usgs_national_dataset()
    
    sf6_hist = build_atmospheric_tracer_input("SF6")
    sf6_years = sf6_hist.years
    sf6_values = sf6_hist.values
    
    # Pre-compute an inverted interpolation for SF6 PFM age
    def sf6_pfm_age(sample_year, sf6_conc):
        if not math.isfinite(sf6_conc) or sf6_conc <= 0:
            return float('nan')
        if sf6_conc >= np.max(sf6_values):
            return 0.0
        if sf6_conc <= np.min(sf6_values):
            return sample_year - np.min(sf6_years)
        
        try:
            recharge_year = np.interp(sf6_conc, sf6_values, sf6_years)
            return max(0.0, sample_year - recharge_year)
        except Exception:
            return float('nan')

    records = []
    
    for _, row in df.iterrows():
        site_id = row['site_id']
        ref_age = usgs._parse_age(row.get('reference_age_years'))
        if pd.isna(ref_age):
            continue
            
        sample_date = row.get('sample_date')
        if pd.isna(sample_date):
            sample_year = 2010.0 # fallback
        else:
            try:
                sample_year = pd.to_datetime(sample_date).year + pd.to_datetime(sample_date).dayofyear / 365.25
            except Exception:
                sample_year = 2010.0
                
        # 3H/3He
        raw_3h = row.get('raw_tritium_TU', float('nan'))
        raw_3he = row.get('raw_he3_trit_TU', float('nan'))
        corr_3h = row.get('tritium_TU', float('nan'))
        corr_3he = row.get('he3_trit_TU', float('nan'))
        
        raw_3h3he_age = _apparent_3h3he_age(raw_3h, raw_3he)
        corr_3h3he_age = _apparent_3h3he_age(corr_3h, corr_3he)
        
        # SF6
        raw_sf6 = row.get('raw_sf6_pptv', float('nan'))
        corr_sf6 = row.get('sf6_pptv', float('nan'))
        
        raw_sf6_age = sf6_pfm_age(sample_year, raw_sf6)
        corr_sf6_age = sf6_pfm_age(sample_year, corr_sf6)
        
        records.append({
            'site_id': site_id,
            'ref_age': ref_age,
            'age_class': usgs._age_class(ref_age),
            'sample_year': sample_year,
            'raw_3h': raw_3h,
            'corr_3h': corr_3h,
            'raw_3he': raw_3he,
            'corr_3he': corr_3he,
            'raw_3h3he_age': raw_3h3he_age,
            'corr_3h3he_age': corr_3h3he_age,
            'raw_sf6': raw_sf6,
            'corr_sf6': corr_sf6,
            'raw_sf6_age': raw_sf6_age,
            'corr_sf6_age': corr_sf6_age,
            'dgm_correction': row.get('dissolved_gas_correction', '')
        })
        
    out_df = pd.DataFrame(records)
    
    # 3H/3He errors
    out_df['raw_3h3he_err'] = np.abs(np.log10(out_df['raw_3h3he_age'] + 1) - np.log10(out_df['ref_age'] + 1))
    out_df['corr_3h3he_err'] = np.abs(np.log10(out_df['corr_3h3he_age'] + 1) - np.log10(out_df['ref_age'] + 1))
    
    # Use a small tolerance to define "improved"
    # To handle NaN comparisons safely, we compare explicitly
    improved_3h3he = []
    for r, c in zip(out_df['raw_3h3he_err'], out_df['corr_3h3he_err']):
        if pd.isna(r) or pd.isna(c):
            improved_3h3he.append(False)
        else:
            improved_3h3he.append(r < c - 1e-6)
    out_df['3h3he_improved_by_raw'] = improved_3h3he
    
    # SF6 errors
    out_df['raw_sf6_err'] = np.abs(np.log10(out_df['raw_sf6_age'] + 1) - np.log10(out_df['ref_age'] + 1))
    out_df['corr_sf6_err'] = np.abs(np.log10(out_df['corr_sf6_age'] + 1) - np.log10(out_df['ref_age'] + 1))
    
    improved_sf6 = []
    for r, c in zip(out_df['raw_sf6_err'], out_df['corr_sf6_err']):
        if pd.isna(r) or pd.isna(c):
            improved_sf6.append(False)
        else:
            improved_sf6.append(r < c - 1e-6)
    out_df['sf6_improved_by_raw'] = improved_sf6
    
    out_file = Path(__file__).parent.parent / "results" / "diagnostic_young_tracers.csv"
    out_df.to_csv(out_file, index=False)
    print(f"Wrote diagnostics to {out_file}")
    
    print("\n--- 3H/3He Summary ---")
    valid_3h3he = out_df.dropna(subset=['raw_3h3he_age', 'corr_3h3he_age', 'ref_age'])
    # Only keep where they differ
    diff_3h3he = valid_3h3he[abs(valid_3h3he['raw_3h3he_age'] - valid_3h3he['corr_3h3he_age']) > 1e-3]
    print(f"Rows with both raw and corrected 3H/3He ages: {len(valid_3h3he)}")
    print(f"Rows where raw and corrected 3H/3He ages differ: {len(diff_3h3he)}")
    if not diff_3h3he.empty:
        print(f"Raw 3H/3He error better than corrected: {diff_3h3he['3h3he_improved_by_raw'].sum()} / {len(diff_3h3he)}")
        print(f"Median log10 error Raw: {diff_3h3he['raw_3h3he_err'].median():.3f}")
        print(f"Median log10 error Corr: {diff_3h3he['corr_3h3he_err'].median():.3f}")

    print("\n--- SF6 Summary ---")
    valid_sf6 = out_df.dropna(subset=['raw_sf6_age', 'corr_sf6_age', 'ref_age'])
    diff_sf6 = valid_sf6[abs(valid_sf6['raw_sf6_age'] - valid_sf6['corr_sf6_age']) > 1e-3]
    print(f"Rows with both raw and corrected SF6 ages: {len(valid_sf6)}")
    print(f"Rows where raw and corrected SF6 ages differ: {len(diff_sf6)}")
    if not diff_sf6.empty:
        print(f"Raw SF6 error better than corrected: {diff_sf6['sf6_improved_by_raw'].sum()} / {len(diff_sf6)}")
        print(f"Median log10 error Raw: {diff_sf6['raw_sf6_err'].median():.3f}")
        print(f"Median log10 error Corr: {diff_sf6['corr_sf6_err'].median():.3f}")
        
    # Also evaluate just the concentrations vs standard
    print("\n--- Interpretation ---")
    print("This checks simple Piston Flow Model (apparent age) differences.")

if __name__ == '__main__':
    main()