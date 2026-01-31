import os
from pathlib import Path

def generate_tree(root_path: Path, exclude_dirs: set):
    """Generates a string representation of the directory tree."""
    tree_str = "Project Structure:\n"
    tree_str += f"{root_path.name}/\n"
    
    for path in sorted(root_path.rglob('*')):
        # Skip excluded directories and their contents
        is_excluded = False
        for part in path.parts:
            if part in exclude_dirs:
                is_excluded = True
                break
        if is_excluded:
            continue
            
        depth = len(path.relative_to(root_path).parts)
        indent = "    " * depth
        if path.is_dir():
             tree_str += f"{indent}{path.name}/\n"
        elif path.suffix == '.py':
             tree_str += f"{indent}{path.name}\n"
             
    return tree_str

def bundle_files(root_dir=".", output_file="hydrosheaf_code_bundle.txt"):
    """
    Bundles all .py files in the directory into a single text file,
    preceded by a directory tree structure.
    """
    root_path = Path(root_dir).resolve()
    
    # Directories to exclude
    exclude_dirs = {
        '.git', '__pycache__', 'venv', 'env', '.pytest_cache', 
        'build', 'dist', 'hydrosheaf.egg-info', '.idea', '.vscode'
    }

    print(f"Bundling python files from: {root_path}")
    
    with open(output_file, "w", encoding="utf-8") as outfile:
        # 1. Write Directory Tree
        outfile.write("=" * 80 + "\n")
        outfile.write("PROJECT STRUCTURE\n")
        outfile.write("=" * 80 + "\n\n")
        outfile.write(generate_tree(root_path, exclude_dirs))
        outfile.write("\n" + "=" * 80 + "\n\n")

        # 2. Write File Contents
        file_count = 0
        for path in sorted(root_path.rglob('*.py')):
            # Skip if in excluded directory
            is_excluded = False
            for part in path.parts:
                if part in exclude_dirs:
                    is_excluded = True
                    break
            if is_excluded:
                continue

            rel_path = path.relative_to(root_path)
            
            outfile.write(f"\n{'=' * 80}\n")
            outfile.write(f"FILE: {rel_path}\n")
            outfile.write(f"{'=' * 80}\n\n")
            
            try:
                with open(path, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    outfile.write(content)
                    outfile.write("\n") # Ensure newline between files
                file_count += 1
                print(f"Added: {rel_path}")
            except Exception as e:
                outfile.write(f"Error reading file: {e}\n")
                print(f"Error reading {rel_path}: {e}")

    print(f"\nBundle created successfully: {output_file}")
    print(f"Total files bundled: {file_count}")

if __name__ == "__main__":
    bundle_files()
