# Hydrosheaf Development Guide

This guide is for **developers** and **contributors** who want to build PEST++ from source or modify the codebase.

## Table of Contents

1. [Development Setup](#development-setup)
2. [Building PEST++ from Source](#building-pest-from-source)
3. [Running Tests](#running-tests)
4. [Contributing](#contributing)

---

## Development Setup

### Prerequisites

- Python >= 3.8
- Git
- Visual Studio 2022 Build Tools (for PEST++ compilation on Windows)
- CMake and Ninja (for PEST++ compilation)

### Clone and Install

```bash
git clone https://github.com/dabdul-wahab1988/Hydrosheaf.git
cd Hydrosheaf

# Install in development mode
pip install -e .

# Install with all optional dependencies
pip install -e .[all]

# Install development tools
pip install pytest pytest-cov black flake8 mypy
```

### Activate Development Environment

```bash
# Create virtual environment (recommended)
python -m venv venv

# Activate
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

# Install dependencies
pip install -e .[all]
```

---

## Building PEST++ from Source

**Note:** This is only necessary if you need to modify PEST++ itself. Most users should:
- Download pre-compiled binaries from [USGS/pestpp releases](https://github.com/usgs/pestpp/releases), OR
- Use Hydrosheaf's core inverse modeling (which doesn't require PEST++)

### Windows (Visual Studio 2022)

#### System Requirements

- Windows 10 or later
- Visual Studio 2022 with Build Tools
- CMake (usually included with VS 2022)
- Ninja (usually included with VS 2022)
- ~500 MB disk space

#### Build Steps

```batch
# Navigate to repo
cd C:\path\to\Hydrosheaf

# Run the build script
compile_pestpp.bat
```

**What the script does:**
1. Sets up Visual Studio 2022 compiler environment
2. Configures build with CMake and Ninja
3. Compiles PEST++ source code
4. Copies compiled `.exe` files to `bin/` directory

#### Manual Build (if script fails)

```batch
# Navigate to PEST++ source
cd pestpp

# Create build directory
mkdir build
cd build

# Configure with CMake
cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release ..

# Build with Ninja
ninja

# Copy binaries to root bin/ directory
copy *.exe ..\..\bin\
```

#### Troubleshooting

**Issue:** "vcvarsall.bat not found"
- **Solution:** Ensure Visual Studio 2022 Build Tools are installed. Check the path in `compile_pestpp.bat` matches your installation.

**Issue:** "CMake not found"
- **Solution:** Install CMake via Visual Studio Installer or download from https://cmake.org

**Issue:** "Ninja not found"
- **Solution:** Install Ninja via Visual Studio Installer or download from https://github.com/ninja-build/ninja/releases

### Linux/macOS

```bash
cd pestpp
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cp bin/* ../../bin/
```

---

## Running Tests

### Core Package Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_edge_fit.py -v

# Run with coverage report
python -m pytest tests/ --cov=hydrosheaf --cov-report=html

# Run PHREEQC-specific tests (requires PHREEQC)
python -m pytest tests/test_phreeqc*.py -v

# Run accuracy/regression tests
python -m pytest tests/test_accuracy*.py -v
```

---

## Development Workflows

### Adding a New Module

1. Create new module in `hydrosheaf/my_module/`
2. Add `__init__.py` with public exports
3. Write implementation files
4. Add corresponding tests in `tests/test_my_module.py`
5. Update documentation in `HYDROSHEAF_USER_MANUAL.md`
6. Run tests: `pytest tests/test_my_module.py`

### Modifying Core Functions

1. Make changes to function in `hydrosheaf/inference/`, `hydrosheaf/models/`, etc.
2. Run relevant tests to ensure no regressions
3. Update docstrings and type hints
4. Add/update tests if behavior changed
5. Update user manual if API changed
6. Commit with descriptive message

### Code Style

```bash
# Format code with Black
black hydrosheaf/ tests/

# Check linting with Flake8
flake8 hydrosheaf/ tests/

# Type checking with MyPy
mypy hydrosheaf/
```

---

## Contributing

### Before Submitting a Pull Request

1. **Fork the repository** and create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following code style guidelines

3. **Add/update tests** for new functionality:
   ```bash
   pytest tests/
   ```

4. **Update documentation:**
   - Update docstrings
   - Update `HYDROSHEAF_USER_MANUAL.md` if API changed
   - Update `README.md` if user-facing

5. **Run quality checks:**
   ```bash
   black hydrosheaf/ tests/
   flake8 hydrosheaf/ tests/
   mypy hydrosheaf/
   pytest tests/ --cov=hydrosheaf
   ```

6. **Commit with clear messages:**
   ```bash
   git commit -m "feat: Add new feature description"
   git commit -m "fix: Fix issue with X"
   git commit -m "docs: Update documentation"
   ```

7. **Push and create Pull Request:**
   ```bash
   git push origin feature/your-feature-name
   ```

### Commit Message Format

- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation changes
- `test:` Test additions/modifications
- `refactor:` Code refactoring (no behavior change)
- `perf:` Performance improvements
- `chore:` Build, CI, dependencies

---

## Documentation

### Updating User Manual

All code changes that affect user-facing APIs should update `HYDROSHEAF_USER_MANUAL.md`:

```markdown
### New Function
**Function:** `my_function(arg1, arg2, config)`  
**Path:** `hydrosheaf.my_module`  
**Returns:** `MyResult`  
**Description:** Detailed description of what it does.
```

### Writing Docstrings

```python
def my_function(x: List[float], y: List[float]) -> float:
    """
    Short description.
    
    Longer description if needed, explaining the algorithm or use case.
    
    Parameters
    ----------
    x : List[float]
        Description of x
    y : List[float]
        Description of y
    
    Returns
    -------
    float
        Description of return value
        
    Examples
    --------
    >>> my_function([1, 2, 3], [4, 5, 6])
    42.0
    """
```

---

## CI/CD Pipeline

The repository uses GitHub Actions for:
- Running pytest on all commits
- Type checking with mypy
- Code coverage reporting
- Building documentation

See `.github/workflows/` for configuration details.

---

## Getting Help

- **Questions about code?** Open a GitHub Discussion
- **Found a bug?** Open a GitHub Issue with reproduction steps
- **Want to contribute?** See [Contributing](#contributing) section above
- **Documentation questions?** Check `HYDROSHEAF_USER_MANUAL.md`

---

## Project Structure

```
Hydrosheaf/
├── hydrosheaf/           # Main Python package
│   ├── api.py            # High-level entry points
│   ├── cli.py            # Command-line interface
│   ├── config.py         # Configuration defaults
│   ├── inference/        # Core inverse solver
│   ├── graph/            # Network topology
│   ├── sheaf/            # Sheaf refinement
│   ├── models/           # Chemical/physical models
│   ├── phreeqc/          # PHREEQC integration
│   ├── temporal/         # Time-series analysis
│   ├── uncertainty/      # Uncertainty quantification
│   ├── calibration/      # PEST++ integration
│   └── ...
├── tests/                # Test suite
├── pestpp/               # PEST++ source code (optional)
├── docs/                 # Additional documentation
├── HYDROSHEAF_USER_MANUAL.md  # User reference
├── README.md             # User overview
└── DEVELOPMENT.md        # This file
```

---

## Performance Profiling

```bash
# Profile with cProfile
python -m cProfile -o results.prof your_script.py

# View results
python -m pstats results.prof
```

---

## Release Process

1. Update version in `hydrosheaf/__init__.py` and `pyproject.toml`
2. Update `CHANGELOG.md` with new features/fixes
3. Commit: `git commit -m "chore: Bump version to x.y.z"`
4. Tag: `git tag vx.y.z`
5. Push: `git push origin main --tags`
6. GitHub Actions will build and upload to PyPI

---

**Last Updated:** January 2026  
**Maintainer:** Dickson Abdul-Wahab
