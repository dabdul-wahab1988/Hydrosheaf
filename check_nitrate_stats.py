"""
Check Nitrate Statistics in Manu.xlsx
"""
import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

def main():
    try:
        df = pd.read_excel("manu.xlsx")
        
        # Check column name - user mentioned NO3
        if 'NO3' in df.columns:
            no3 = df['NO3']
            print(f"Nitrate (NO3) Stats [mg/L]:")
            print(f"  Count: {no3.count()}")
            print(f"  Min:   {no3.min():.2f}")
            print(f"  Max:   {no3.max():.2f}")
            print(f"  Mean:  {no3.mean():.2f}")
            print(f"  Std:   {no3.std():.2f}")
            
            # Check for outliers
            print("\nTop 5 Highest Nitrate Samples:")
            # Use Station or Sample ID
            id_col = 'Station' if 'Station' in df.columns else 'Sample ID'
            print(df.nlargest(5, 'NO3')[[id_col, 'NO3']])
            
        else:
            print("Column 'NO3' not found. Available columns:", df.columns.tolist())
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
