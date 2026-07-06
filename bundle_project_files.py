import argparse
from pathlib import Path


EXCLUDE_DIRS = {
    ".git",
    "__pycache__",
    "venv",
    ".venv",
    "env",
    ".pytest_cache",
    "pymc_cache",
    "build",
    "dist",
    "hydrosheaf.egg-info",
    ".idea",
    ".vscode",
    "tests",
    "scratch",
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

M4_BUNDLE_PATHS = (
    "pyproject.toml",
    "README.md",
    "REPRODUCIBILITY.md",
    "bundle_project_files.py",
    "hydrosheaf",
    "M4/m4_topology_benchmark",
)

M4_MINIMAL_BUNDLE_PATHS = (
    "M4/m4_topology_benchmark/scripts",
    "hydrosheaf/graph/build.py",
    "hydrosheaf/physics/modpath.py",
    "hydrosheaf/validation/__init__.py",
    "hydrosheaf/validation/modpath_archive.py",
    "hydrosheaf/validation/topology.py",
)

DEFAULT_BUNDLE_EXTENSIONS = (".py",)
M4_BUNDLE_EXTENSIONS = (".py", ".md", ".yaml", ".yml", ".toml", ".json", ".csv", ".dat")


def is_excluded(path: Path, exclude_dirs: set[str]) -> bool:
    return any(part in exclude_dirs for part in path.parts)


def read_text_file(path: Path) -> str:
    try:
        with open(path, "r", encoding="utf-8") as infile:
            return infile.read()
    except UnicodeDecodeError:
        with open(path, "r", encoding="latin-1") as infile:
            return infile.read()


def normalise_bundle_content(content: str) -> str:
    """Strip trailing whitespace so generated bundle files pass git checks."""
    return "\n".join(line.rstrip() for line in content.splitlines())


def iter_bundle_files(
    root_path: Path,
    include_paths: tuple[str, ...] | None,
    exclude_dirs: set[str],
    extensions: tuple[str, ...] = DEFAULT_BUNDLE_EXTENSIONS,
) -> list[Path]:
    paths: list[Path] = []
    seen: set[Path] = set()
    allowed_extensions = {ext.lower() for ext in extensions}
    roots = [root_path / rel for rel in include_paths] if include_paths else [root_path]

    for current_root in roots:
        if not current_root.exists():
            print(f"Warning: include path does not exist: {current_root}")
            continue

        candidates = (
            (path for path in current_root.rglob("*") if path.is_file())
            if current_root.is_dir()
            else [current_root]
        )
        for path in candidates:
            path = path.resolve()
            if path.suffix.lower() not in allowed_extensions or is_excluded(path, exclude_dirs):
                continue
            if path not in seen:
                seen.add(path)
                paths.append(path)

    return sorted(paths, key=lambda item: item.relative_to(root_path).as_posix().lower())


def generate_tree(
    root_path: Path,
    exclude_dirs: set[str],
    files: list[Path] | None = None,
    extensions: tuple[str, ...] = DEFAULT_BUNDLE_EXTENSIONS,
):
    """Generates a string representation of the directory tree."""
    tree_str = "Project Structure:\n"
    tree_str += f"{root_path.name}/\n"
    allowed_extensions = {ext.lower() for ext in extensions}

    def tree_sort_key(path: Path) -> tuple[tuple[str, int], ...]:
        rel_path = path.relative_to(root_path)
        parts = rel_path.parts
        key_parts: list[tuple[str, int]] = []
        for idx, part in enumerate(parts):
            partial = root_path / Path(*parts[: idx + 1])
            key_parts.append((part.lower(), 0 if partial.is_dir() else 1))
        return tuple(key_parts)

    if files is None:
        tree_paths = [
            path
            for path in sorted(root_path.rglob("*"), key=tree_sort_key)
            if not is_excluded(path, exclude_dirs)
            and (path.is_dir() or path.suffix.lower() in allowed_extensions)
        ]
    else:
        rel_entries: set[Path] = set()
        for file_path in files:
            rel_path = file_path.relative_to(root_path)
            rel_entries.add(rel_path)
            rel_entries.update(rel_path.parents)
        rel_entries.discard(Path("."))
        tree_paths = sorted((root_path / rel for rel in rel_entries), key=tree_sort_key)

    for path in tree_paths:
        depth = len(path.relative_to(root_path).parts)
        indent = "    " * depth
        if path.is_dir():
            tree_str += f"{indent}{path.name}/\n"
        elif path.suffix.lower() in allowed_extensions:
            tree_str += f"{indent}{path.name}\n"

    return tree_str


def bundle_files(
    root_dir=".",
    output_file="hydrosheaf_code_bundle.txt",
    include_paths: tuple[str, ...] | None = None,
    extensions: tuple[str, ...] = DEFAULT_BUNDLE_EXTENSIONS,
):
    """
    Bundles source and text data files into a single text file, preceded by a directory tree.
    """
    root_path = Path(root_dir).resolve()
    files = iter_bundle_files(root_path, include_paths, EXCLUDE_DIRS, extensions)

    print(f"Bundling files from: {root_path}")
    print("Extensions: " + ", ".join(extensions))
    if include_paths:
        print("Profile paths:")
        for rel in include_paths:
            print(f"  - {rel}")

    with open(output_file, "w", encoding="utf-8") as outfile:
        # 1. Write Directory Tree
        outfile.write("=" * 80 + "\n")
        outfile.write("PROJECT STRUCTURE\n")
        outfile.write("=" * 80 + "\n\n")
        outfile.write(generate_tree(root_path, EXCLUDE_DIRS, files if include_paths else None, extensions))
        outfile.write("\n" + "=" * 80 + "\n\n")

        # 2. Write File Contents
        file_count = 0
        for path in files:
            rel_path = path.relative_to(root_path)

            outfile.write(f"\n{'=' * 80}\n")
            outfile.write(f"FILE: {rel_path}\n")
            outfile.write(f"{'=' * 80}\n\n")

            try:
                content = normalise_bundle_content(read_text_file(path))
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
    parser = argparse.ArgumentParser(description="Bundle Hydrosheaf source files into one text file.")
    parser.add_argument(
        "--profile",
        choices=["all", "m2", "m4", "m4-minimal"],
        default="all",
        help=(
            "Use 'm2' for M2 analysis/figure code, 'm4' for Hydrosheaf plus the "
            "M4 topology benchmark, or 'm4-minimal' for reviewer-facing M4 code only."
        ),
    )
    parser.add_argument("--root-dir", default=".", help="Repository root to scan.")
    parser.add_argument("--output", default=None, help="Output bundle file.")
    args = parser.parse_args()

    include_paths = None
    extensions = DEFAULT_BUNDLE_EXTENSIONS
    default_output = "hydrosheaf_code_bundle.txt"
    if args.profile == "m2":
        include_paths = M2_BUNDLE_PATHS
        default_output = "m2_analysis_figure_code_bundle.txt"
    elif args.profile == "m4":
        include_paths = M4_BUNDLE_PATHS
        extensions = M4_BUNDLE_EXTENSIONS
        default_output = "hydrosheaf_m4_benchmark_code_bundle.txt"
    elif args.profile == "m4-minimal":
        include_paths = M4_MINIMAL_BUNDLE_PATHS
        default_output = "hydrosheaf_m4_minimal_code_bundle.txt"

    output_file = args.output or default_output
    bundle_files(args.root_dir, output_file, include_paths, extensions)


if __name__ == "__main__":
    main()
