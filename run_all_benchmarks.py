#!/usr/bin/env python3
"""
Master script to run all PhD benchmarks for Hydrosheaf (M2-M5).
This ensures reproducibility and allows examiners to regenerate claims.
"""

import subprocess
import sys
import os

def run_script(path, args=None):
    print(f"\n{'='*60}")
    print(f"  Running: {path}")
    print(f"{'='*60}\n")
    
    cmd = [sys.executable, path]
    if args:
        cmd.extend(args)
    
    try:
        subprocess.check_call(cmd)
        print(f"\n[✓] SUCCESS: {path}\n")
    except subprocess.CalledProcessError as e:
        print(f"\n[✗] FAILED: {path} (Exit code: {e.returncode})\n")
        return False
    return True

def main():
    root = os.getcwd()
    
    benchmarks = [
        ("M2/m2_benchmark/scripts/run_m2_benchmark.py", ["--realisations", "10"]),
        ("M3/m3_age_benchmark/scripts/run_m3_age_benchmark.py", []),
        ("M4/m4_topology_benchmark/scripts/run_m4_topology_benchmark.py", []),
        ("M5/m5_inverse_reaction_benchmark/scripts/run_m5_inverse_reaction_benchmark.py", []),
    ]
    
    success_count = 0
    for script, args in benchmarks:
        script_path = os.path.join(root, script)
        if not os.path.exists(script_path):
            print(f"[!] WARNING: Script not found: {script_path}")
            continue
            
        if run_script(script_path, args):
            success_count += 1
            
    print(f"\n{'-'*60}")
    print(f"Benchmark Summary: {success_count}/{len(benchmarks)} completed successfully.")
    print(f"{'-'*60}")

if __name__ == "__main__":
    main()
