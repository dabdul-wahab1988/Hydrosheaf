
import os
import glob
from pathlib import Path

def clean_tex_artifacts(directory):
    extensions = [
        "*.aux", "*.toc", "*.out", "*.bcf", "*.run.xml", 
        "*.bbl", "*.blg", "*.log", "*.fls", "*.fdb_latexmk",
        "*.lof", "*.lot", "*.synctex.gz"
    ]
    
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"Directory {directory} does not exist.")
        return

    print(f"Cleaning TeX artifacts in {directory}...")
    for ext in extensions:
        for file_path in dir_path.glob(ext):
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path.name}")
            except Exception as e:
                print(f"Error deleting {file_path.name}: {e}")

if __name__ == "__main__":
    clean_tex_artifacts("hydrosheaf_manual")
