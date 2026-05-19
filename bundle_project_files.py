import argparse
from pathlib import Path


EXCLUDE_DIRS = {
    ".git",
    "__pycache__",
    "venv",
    "env",
    ".pytest_cache",
    "build",
    "dist",
    "hydrosheaf.egg-info",
    ".idea",
    ".vscode",
}

M2_BUNDLE_PATHS = (
    "bundle_project_files.py",
    "generate_m2_md_tables.py",
    "hydrosheaf",
    "M2/m2_benchmark/scripts",
    "scripts/analysis",
    "M3/m3_age_benchmark/scripts",
    "M6/m6_robustness_benchmark/scripts",
)


def is_excluded(path: Path, exclude_dirs: set[str]) -> bool:
    return any(part in exclude_dirs for part in path.parts)


def iter_python_files(root_path: Path, include_paths: tuple[str, ...] | None, exclude_dirs: set[str]) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    roots = [root_path / rel for rel in include_paths] if include_paths else [root_path]

    for current_root in roots:
        if not current_root.exists():
            print(f"Warning: include path does not exist: {current_root}")
            continue

        candidates = current_root.rglob("*.py") if current_root.is_dir() else [current_root]
        for path in candidates:
            path = path.resolve()
            if path.suffix != ".py" or is_excluded(path, exclude_dirs):
                continue
            if path not in seen:
                seen.add(path)
                paths.append(path)

    return sorted(paths, key=lambda item: item.relative_to(root_path).as_posix().lower())


def generate_tree(root_path: Path, exclude_dirs: set[str], files: list[Path] | None = None):
    """Generates a string representation of the directory tree."""
    tree_str = "Project Structure:\n"
    tree_str += f"{root_path.name}/\n"

    if files is None:
        tree_paths = [path for path in sorted(root_path.rglob("*")) if not is_excluded(path, exclude_dirs)]
    else:
        rel_entries: set[Path] = set()
        for file_path in files:
            rel_path = file_path.relative_to(root_path)
            rel_entries.add(rel_path)
            rel_entries.update(rel_path.parents)
        rel_entries.discard(Path("."))
        tree_paths = [root_path / rel for rel in sorted(rel_entries, key=lambda item: item.as_posix().lower())]

    for path in tree_paths:
        depth = len(path.relative_to(root_path).parts)
        indent = "    " * depth
        if path.is_dir():
            tree_str += f"{indent}{path.name}/\n"
        elif path.suffix == ".py":
            tree_str += f"{indent}{path.name}\n"

    return tree_str


def bundle_files(root_dir=".", output_file="hydrosheaf_code_bundle.txt", include_paths: tuple[str, ...] | None = None):
    """
    Bundles .py files into a single text file, preceded by a directory tree.
    """
    root_path = Path(root_dir).resolve()
    files = iter_python_files(root_path, include_paths, EXCLUDE_DIRS)

    print(f"Bundling python files from: {root_path}")
    if include_paths:
        print("Profile paths:")
        for rel in include_paths:
            print(f"  - {rel}")

    with open(output_file, "w", encoding="utf-8") as outfile:
        # 1. Write Directory Tree
        outfile.write("=" * 80 + "\n")
        outfile.write("PROJECT STRUCTURE\n")
        outfile.write("=" * 80 + "\n\n")
        outfile.write(generate_tree(root_path, EXCLUDE_DIRS, files if include_paths else None))
        outfile.write("\n" + "=" * 80 + "\n\n")

        # 2. Write File Contents
        file_count = 0
        for path in files:
            rel_path = path.relative_to(root_path)

            outfile.write(f"\n{'=' * 80}\n")
            outfile.write(f"FILE: {rel_path}\n")
            outfile.write(f"{'=' * 80}\n\n")

            try:
                with open(path, "r", encoding="utf-8") as infile:
                    content = infile.read()
                    outfile.write(content)
                    outfile.write("\n")
                file_count += 1
                print(f"Added: {rel_path}")
            except Exception as e:
                outfile.write(f"Error reading file: {e}\n")
                print(f"Error reading {rel_path}: {e}")

    print(f"\nBundle created successfully: {output_file}")
    print(f"Total files bundled: {file_count}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Bundle Hydrosheaf Python source files into one text file.")
    parser.add_argument(
        "--profile",
        choices=["all", "m2"],
        default="all",
        help="Use 'm2' for M2 analysis and manuscript figure code.",
    )
    parser.add_argument("--root-dir", default=".", help="Repository root to scan.")
    parser.add_argument("--output", default=None, help="Output bundle file.")
    args = parser.parse_args()

    include_paths = M2_BUNDLE_PATHS if args.profile == "m2" else None
    output_file = args.output or ("m2_analysis_figure_code_bundle.txt" if args.profile == "m2" else "hydrosheaf_code_bundle.txt")
    bundle_files(args.root_dir, output_file, include_paths)


if __name__ == "__main__":
    main()
