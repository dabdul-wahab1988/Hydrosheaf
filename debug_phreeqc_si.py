"""
Debug PHREEQC Saturation Indices directly.
"""
import sys
import os
from pathlib import Path
import traceback

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from phreeqpython import PhreeqPython
except ImportError:
    print("phreeqpython not installed.")
    sys.exit(1)

def main():
    print("Debug: Testing PhreeqPython...")
    
    # Path to database
    db_path = str(Path.cwd() / "hydrosheaf" / "databases" / "phreeqc.dat")
    if not os.path.exists(db_path):
        print(f"Error: Database not found at {db_path}")
        return

    print(f"Using database: {db_path}")

    pp = None
    try:
        # Check constructor: Does it take 'database'?
        # In phreeqpython, it often looks for 'phreeqc.dat' in current or internal dir.
        # But we can try to pass full path.
        pp = PhreeqPython(database=db_path)
    except Exception as e:
        print(f"Error initializing PhreeqPython with path: {e}")
        # Try default
        print("Trying default database...")
        try:
            pp = PhreeqPython() 
        except Exception as e2:
            print(f"Error initializing default PhreeqPython: {e2}")
            return

    if not pp:
        print("Failed to initialize PhreeqPython.")
        return

    # Define Solution L1
    # pH=6.4, Temp=25
    # Ca=40.64 mg/L
    # Mg=23.67 mg/L
    # Na=54.23 mg/L
    # K=4.5 mg/L
    # HCO3=296 mg/L
    # Cl=26.74 mg/L
    # SO4=12.74 mg/L
    # NO3=42.35 mg/L
    
    print("Adding solution L1...")
    try:
        sol1 = pp.add_solution({
            'pH': 6.4,
            'temp': 25.0,
            'units': 'mg/L',
            'Ca': 40.64,
            'Mg': 23.67,
            'Na': 54.23,
            'K': 4.5,
            'Alkalinity': 296.0, # Will try 'Alkalinity' (as CaCO3 usually)
                                 # Or 'C(4)'?
                                 # Let's try explicit species if phreeqpython supports key 'HCO3' -> C(4)
                                 # Or try standard 'Alkalinity'
            'Cl': 26.74,
            'S(6)': 12.74, # SO4
            'N(5)': 42.35  # NO3
        })
    except Exception as e:
        print(f"Error adding solution: {e}")
        traceback.print_exc()
        return
    
    # Let's try to get SI
    print("\nSaturation Indices (L1):")
    try:
        # Note: phreeqpython .si() method
        si_calcite = sol1.si("Calcite")
        si_gypsum = sol1.si("Gypsum")
        si_dolomite = sol1.si("Dolomite")
        
        print(f"  Calcite: {si_calcite}")
        print(f"  Gypsum: {si_gypsum}")
        print(f"  Dolomite: {si_dolomite}")
        
    except Exception as e:
        print(f"Error calculating SI: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    main()
