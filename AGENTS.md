# Hydrosheaf AGENTS.md

This file guides agentic coding assistants working in this repository.
It summarizes build/lint/test commands and local coding conventions.

## Repository Overview

- Core Python package lives in `hydrosheaf/`.
- Tests live in `tests/` (pytest + unittest-style tests).

## Cursor/Copilot Rules

- No `.cursor/rules/`, `.cursorrules`, or `.github/copilot-instructions.md` files were found.
- If these are added later, follow them and update this document.

## Environment Prerequisites

- Python >= 3.8
- Optional: PHREEQC and related Python bindings for thermodynamic constraints

## Core Package: Install & Run

- Install core package: `pip install .`
- Install with plotting support: `pip install .[plot]`
- Install with PHREEQC support: `pip install .[phreeqc]`
- Install everything: `pip install .[all]`
- Run CLI: `hydrosheaf --help`

## Core Package: Tests

- Run all tests: `python -m pytest tests/`
- Run a single test file: `python -m pytest tests/test_edge_fit.py`
- Run a single test case: `python -m pytest tests/test_edge_fit.py::EdgeFitTests::test_edge_fit_evap_and_halite`
- PHREEQC-specific tests (requires PHREEQC): `python -m pytest tests/test_phreeqc*.py`
- Accuracy/regression tests: `python -m pytest tests/test_accuracy*.py`
- Documentation examples: `python -m pytest tests/test_doc_examples.py`


## General Coding Guidelines

- Prefer clear, minimal changes that match existing patterns.
- Avoid introducing new dependencies unless required.
- Keep public APIs stable and update tests when behavior changes.
- Preserve domain terminology (ion names, mineral labels, hydrologic terms).
- Use descriptive names; avoid one-letter variables outside small loops.
- Avoid inline comments unless necessary; favor readable code.
- Do not add license headers or large boilerplate blocks.

## Python Style

- Use standard library imports first, third-party next, local imports last.
- Group imports with a blank line between groups.
- Keep imports sorted alphabetically within each group.
- Use `from __future__ import annotations` when a module already uses it.
- Use type hints for public functions and complex data structures.
- Prefer `List`, `Dict`, `Tuple`, `Optional` from `typing` for clarity.
- Favor dataclasses for configuration/state containers.
- Keep functions short and focused; extract helpers for repeated logic.
- Use `snake_case` for functions/variables and `PascalCase` for classes.
- Use ALL_CAPS for constants (see `DEFAULT_ION_ORDER`).
- Prefer explicit error handling with clear messages (`ValueError`, `TypeError`).
- Avoid raising generic `Exception` unless there is no better option.
- Validate inputs near boundaries (API adapters, public entry points).
- Keep docstrings brief and informative for public helpers.

## Data & Units Conventions

- Internal concentration units are mmol/L.
- Ion order defaults to `Config.DEFAULT_ION_ORDER` (10-ion order).
- Respect configured `Config.ion_order` and `Config.weights` when fitting.

## Testing Conventions

- Tests use pytest runner with unittest-style test classes.
- Follow existing naming patterns: `test_*.py`, `*Tests` classes.
- For new tests, keep fixtures local to the test file unless reused.
- Avoid slow tests unless necessary; separate integration tests clearly.
- When adding new features, update or add targeted tests.

## Error Handling & Validation

- Use explicit validation for required inputs (see `api.validate_required_inputs`).
- Prefer returning structured results and diagnostics rather than printing.
- Maintain consistent exception messages (full sentence, actionable context).
- Do not silently swallow errors in public APIs.

## Logging & Diagnostics

- The core package primarily uses returned objects for diagnostics.
- If adding logging, keep it lightweight and optional.



## When Updating This File

- Update commands if tooling changes.
- Add Cursor/Copilot rules if they appear.
- Keep the document concise and around ~150 lines.
