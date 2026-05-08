"""Run E1 public groundwater-age validation for the M2 benchmark.

The E1 tier compares Hydrosheaf single-node residence-time estimates against
published TracerLPM-style outputs from the USGS public-supply aquifer
groundwater-age data release (DOI 10.5066/P9W7T0DN).

The script is deliberately conservative: it will not generate Figure 4A or
validation metrics unless it can read actual source tables from the USGS release.
If ScienceBase is unreachable, it writes a failure report with the exact blocker
and leaves empty schema files for downstream reproducibility checks.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import traceback
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover - handled in report if plotting unavailable.
    plt = None

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.nuclear.multi_tracer import infer_multi_tracer_age


BENCHMARK_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = BENCHMARK_ROOT / "external" / "usgs_age" / "input"
RESULT_DIR = BENCHMARK_ROOT / "external" / "usgs_age" / "results"
FIGURE_DIR = BENCHMARK_ROOT / "figures"

CATALOG_URL = (
    "https://data.usgs.gov/datacatalog/data/"
    "USGS:62857962d34e3bef0c9a6db8"
)
METADATA_XML_URL = (
    "https://data.usgs.gov/datacatalog/metadata/"
    "USGS.62857962d34e3bef0c9a6db8.xml"
)
SCIENCEBASE_ITEM_JSON = (
    "https://www.sciencebase.gov/catalog/item/"
    "62857962d34e3bef0c9a6db8?fields=files&format=json"
)
SCIENCEBASE_FILE_ZIP = (
    "https://www.sciencebase.gov/catalog/file/get/"
    "62857962d34e3bef0c9a6db8"
)

VALIDATION_COLUMNS = [
    "site_id",
    "sample_id",
    "sample_year",
    "reference_mean_age_years",
    "hydrosheaf_age_years",
    "hydrosheaf_ci_low_years",
    "hydrosheaf_ci_high_years",
    "hydrosheaf_method",
    "tritium_TU",
    "he3_trit_TU",
    "sf6_pptv",
    "cfc11_pptv",
    "cfc12_pptv",
    "cfc113_pptv",
    "c14_pmc",
    "supported_tracer_count",
    "supported_tracers",
    "reference_anthropocene_fraction",
    "reference_holocene_fraction",
    "reference_pleistocene_fraction",
    "log10_error",
    "inside_hydrosheaf_ci",
]


@dataclass
class SourceAttempt:
    url: str
    target: Optional[str]
    status: str
    detail: str


@dataclass
class TableCandidate:
    name: str
    path: str
    frame: pd.DataFrame


def clean_name(value: Any) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^0-9a-zA-Z]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text


def ensure_dirs() -> None:
    for path in (INPUT_DIR, RESULT_DIR, FIGURE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def fetch_url(url: str, target: Path, timeout: int = 90) -> SourceAttempt:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "Hydrosheaf-M2-E1-validation/1.0",
            "Accept": "*/*",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            content = response.read()
            status = getattr(response, "status", 200)
            ctype = response.headers.get("content-type", "")
        if status != 200:
            return SourceAttempt(url, str(target), "failed", f"HTTP {status}")
        if not content:
            return SourceAttempt(url, str(target), "failed", "empty response body")
        if "text/html" in ctype.lower() and target.suffix.lower() not in {".html", ".htm"}:
            snippet = content[:160].decode("utf-8", errors="replace")
            return SourceAttempt(
                url,
                str(target),
                "failed",
                f"unexpected HTML response instead of data file: {snippet}",
            )
        target.write_bytes(content)
        return SourceAttempt(url, str(target), "ok", f"{len(content)} bytes")
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
        return SourceAttempt(url, str(target), "failed", repr(exc))


def download_sources() -> List[SourceAttempt]:
    attempts: List[SourceAttempt] = []
    attempts.append(fetch_url(METADATA_XML_URL, INPUT_DIR / "usgs_age_metadata.xml"))

    item_json = INPUT_DIR / "sciencebase_item_files.json"
    attempts.append(fetch_url(SCIENCEBASE_ITEM_JSON, item_json))
    if item_json.exists():
        try:
            payload = json.loads(item_json.read_text(encoding="utf-8"))
            for item_file in payload.get("files", []) or []:
                name = item_file.get("name") or item_file.get("title")
                uri = item_file.get("downloadUri") or item_file.get("url")
                if not name or not uri:
                    continue
                if not re.search(r"\.(csv|txt|xlsx?|zip)$", name, flags=re.I):
                    continue
                attempts.append(fetch_url(uri, INPUT_DIR / name))
        except Exception as exc:
            attempts.append(
                SourceAttempt(
                    str(item_json),
                    None,
                    "failed",
                    f"could not parse ScienceBase file manifest: {exc!r}",
                )
            )

    archive = INPUT_DIR / "sciencebase_usgs_age_release.zip"
    attempts.append(fetch_url(SCIENCEBASE_FILE_ZIP, archive))
    if archive.exists() and zipfile.is_zipfile(archive):
        extract_dir = INPUT_DIR / "sciencebase_usgs_age_release"
        extract_dir.mkdir(exist_ok=True)
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(extract_dir)
        attempts.append(
            SourceAttempt(
                str(archive),
                str(extract_dir),
                "ok",
                f"extracted {len(zf.namelist())} files",
            )
        )
    elif archive.exists():
        attempts.append(
            SourceAttempt(
                str(archive),
                None,
                "failed",
                "downloaded file is not a valid zip archive",
            )
        )

    return attempts


def candidate_files(input_dir: Path) -> List[Path]:
    files: List[Path] = []
    for suffix in ("*.csv", "*.txt", "*.tsv", "*.xls", "*.xlsx"):
        files.extend(input_dir.rglob(suffix))
    return [
        path
        for path in sorted(files)
        if path.name.lower()
        not in {"readme.md", "sciencebase_item_files.json", "usgs_age_metadata.xml"}
        and "readme" not in path.name.lower()
    ]


def read_table_file(path: Path) -> List[TableCandidate]:
    candidates: List[TableCandidate] = []
    suffix = path.suffix.lower()
    try:
        if suffix in {".xlsx", ".xls"}:
            sheets = pd.read_excel(path, sheet_name=None)
            for sheet_name, frame in sheets.items():
                if not frame.empty:
                    candidates.append(
                        TableCandidate(f"{path.name}:{sheet_name}", str(path), normalize_columns(frame))
                    )
        else:
            sep = "\t" if suffix in {".txt", ".tsv"} else ","
            try:
                frame = pd.read_csv(path, sep=sep, low_memory=False)
            except UnicodeDecodeError:
                frame = pd.read_csv(path, sep=sep, low_memory=False, encoding="cp1252")
            if not frame.empty:
                candidates.append(TableCandidate(path.name, str(path), normalize_columns(frame)))
    except Exception as exc:
        print(f"[warn] Could not read {path}: {exc}")
    return candidates


def normalize_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out.columns = [clean_name(c) for c in out.columns]
    return out


def discover_tables(input_dir: Path) -> List[TableCandidate]:
    candidates: List[TableCandidate] = []
    for path in candidate_files(input_dir):
        candidates.extend(read_table_file(path))
    return candidates


def score_table(frame: pd.DataFrame, patterns: Iterable[str]) -> int:
    cols = list(frame.columns)
    score = 0
    for pattern in patterns:
        if any(re.search(pattern, col, flags=re.I) for col in cols):
            score += 1
    return score


def pick_table(candidates: List[TableCandidate], patterns: Iterable[str]) -> Optional[TableCandidate]:
    scored = [
        (score_table(candidate.frame, patterns), len(candidate.frame), candidate)
        for candidate in candidates
    ]
    scored = [item for item in scored if item[0] > 0]
    if not scored:
        return None
    scored.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return scored[0][2]


def first_matching_column(frame: pd.DataFrame, patterns: Iterable[str]) -> Optional[str]:
    for pattern in patterns:
        for col in frame.columns:
            if re.search(pattern, col, flags=re.I):
                return col
    return None


def choose_join_column(left: pd.DataFrame, right: pd.DataFrame) -> Optional[str]:
    priority = [
        r"^sample_id$",
        r"sample.*id",
        r"^site_id$",
        r"site.*id",
        r"station.*id",
        r"well.*id",
        r"usgs.*id",
        r"local.*id",
    ]
    common = set(left.columns).intersection(right.columns)
    for pattern in priority:
        for col in common:
            if re.search(pattern, col, flags=re.I):
                return col
    return None


def to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(
        series.astype(str).str.replace("<", "", regex=False).str.replace(">", "", regex=False),
        errors="coerce",
    )


def decimal_year(value: Any, default: float = 2011.0) -> float:
    if pd.isna(value):
        return default
    try:
        num = float(value)
        if 1900 <= num <= 2100:
            return num
    except Exception:
        pass
    parsed = pd.to_datetime(value, errors="coerce")
    if pd.isna(parsed):
        return default
    start = pd.Timestamp(year=parsed.year, month=1, day=1)
    end = pd.Timestamp(year=parsed.year + 1, month=1, day=1)
    return float(parsed.year + (parsed - start).total_seconds() / (end - start).total_seconds())


def infer_sample_age(row: pd.Series) -> Dict[str, Any]:
    sample_year = decimal_year(row.get("sample_year"), default=2011.0)
    estimate = infer_multi_tracer_age(
        {
            "tritium_TU": row.get("tritium_TU"),
            "tritium_sigma_TU": row.get("tritium_sigma_TU"),
            "he3_trit_TU": row.get("he3_trit_TU"),
            "he3_trit_sigma_TU": row.get("he3_trit_sigma_TU"),
            "sf6_pptv": row.get("sf6_pptv"),
            "sf6_sigma_pptv": row.get("sf6_sigma_pptv"),
            "cfc11_pptv": row.get("cfc11_pptv"),
            "cfc11_sigma_pptv": row.get("cfc11_sigma_pptv"),
            "cfc12_pptv": row.get("cfc12_pptv"),
            "cfc12_sigma_pptv": row.get("cfc12_sigma_pptv"),
            "cfc113_pptv": row.get("cfc113_pptv"),
            "cfc113_sigma_pptv": row.get("cfc113_sigma_pptv"),
            "c14_pmc": row.get("c14_pmc"),
            "c14_sigma_pmc": row.get("c14_sigma_pmc"),
            "he4_ccpg": row.get("he4_ccpg"),
        },
        sample_year=sample_year,
    )
    return {
        "hydrosheaf_age_years": estimate["age_years"],
        "hydrosheaf_ci_low_years": estimate["ci_low_years"],
        "hydrosheaf_ci_high_years": estimate["ci_high_years"],
        "hydrosheaf_method": estimate["method"],
        "supported_tracer_count": estimate["n_tracers"],
        "supported_tracers": ";".join(estimate["tracers"]),
    }


def build_validation_dataset(candidates: List[TableCandidate]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    age_table = pick_table(
        candidates,
        [
            r"mean.*age",
            r"age.*mean",
            r"anthropocene",
            r"holocene",
            r"pleistocene",
            r"lpm",
        ],
    )
    tracer_table = pick_table(
        candidates,
        [
            r"tritium|h3|3h",
            r"carbon.*14|c14|14c|pmc",
            r"sf6",
            r"helium|he3|he4",
        ],
    )
    site_table = pick_table(
        candidates,
        [
            r"principal.*aquifer",
            r"aquifer",
            r"latitude|longitude|lat|lon",
            r"well.*depth",
        ],
    )

    diagnostics: Dict[str, Any] = {
        "candidate_tables": [
            {"name": c.name, "path": c.path, "rows": len(c.frame), "columns": list(c.frame.columns)}
            for c in candidates
        ],
        "age_table": age_table.name if age_table else None,
        "tracer_table": tracer_table.name if tracer_table else None,
        "site_table": site_table.name if site_table else None,
    }

    if age_table is None or tracer_table is None:
        return pd.DataFrame(columns=VALIDATION_COLUMNS), diagnostics

    join_col = choose_join_column(age_table.frame, tracer_table.frame)
    if join_col is None:
        diagnostics["join_error"] = "No common sample/site identifier column found."
        return pd.DataFrame(columns=VALIDATION_COLUMNS), diagnostics

    merged = age_table.frame.merge(
        tracer_table.frame,
        on=join_col,
        how="inner",
        suffixes=("_age", "_tracer"),
    )
    if site_table is not None:
        site_join = choose_join_column(merged, site_table.frame)
        if site_join is not None:
            merged = merged.merge(site_table.frame, on=site_join, how="left", suffixes=("", "_site"))

    ref_age_col = first_matching_column(
        merged,
        [
            r"^rpt_totage_yrs$",
            r"^totage_yrs$",
            r"tot.*age.*yrs",
            r"total.*age",
            r"mean.*age",
            r"age.*mean",
            r"meanage",
            r"mean.*residence",
            r"mrt",
        ],
    )
    tritium_col = first_matching_column(merged, [r"^3h_tu$", r"tritium.*tu", r"tritium", r"^h3$", r"3h"])
    tritium_sigma_col = first_matching_column(merged, [r"^3h_err_tu$", r"tritium.*err", r"h3.*err"])
    he3_col = first_matching_column(merged, [r"^3he_trit_tu$", r"he3.*trit", r"3he.*trit"])
    he3_sigma_col = first_matching_column(merged, [r"^3he_trit_err_tu$", r"he3.*err", r"3he.*err"])
    sf6_col = first_matching_column(merged, [r"^sf6_pptv$", r"sf6"])
    sf6_sigma_col = first_matching_column(merged, [r"^sf6_err_pptv$", r"sf6.*err"])
    cfc11_col = first_matching_column(merged, [r"^cfc_?11_pptv$", r"cfc.*11"])
    cfc11_sigma_col = first_matching_column(merged, [r"^cfc_?11_err_pptv$", r"cfc.*11.*err"])
    cfc12_col = first_matching_column(merged, [r"^cfc_?12_pptv$", r"cfc.*12"])
    cfc12_sigma_col = first_matching_column(merged, [r"^cfc_?12_err_pptv$", r"cfc.*12.*err"])
    cfc113_col = first_matching_column(merged, [r"^cfc_?113_pptv$", r"cfc.*113"])
    cfc113_sigma_col = first_matching_column(merged, [r"^cfc_?113_err_pptv$", r"cfc.*113.*err"])
    c14_col = first_matching_column(merged, [r"^14c_pmc$", r"carbon.*14.*pmc", r"c14.*pmc", r"14c.*pmc", r"pmc"])
    c14_sigma_col = first_matching_column(merged, [r"^14c_err_pmc$", r"c14.*err", r"14c.*err"])
    he4_col = first_matching_column(merged, [r"^4he_ccpg$", r"he4.*ccpg", r"4he"])
    sample_year_col = first_matching_column(
        merged,
        [r"sample.*year", r"sample.*date", r"date.*sample", r"collection.*date"],
    )
    sample_id_col = first_matching_column(merged, [r"sample.*id", r"^sample_id$"]) or join_col
    site_id_col = first_matching_column(merged, [r"site.*id", r"station.*id", r"well.*id", r"usgs.*id"])

    diagnostics.update(
        {
            "join_column": join_col,
            "reference_age_column": ref_age_col,
            "tritium_column": tritium_col,
            "he3_tritium_column": he3_col,
            "sf6_column": sf6_col,
            "cfc11_column": cfc11_col,
            "cfc12_column": cfc12_col,
            "cfc113_column": cfc113_col,
            "c14_column": c14_col,
            "sample_year_column": sample_year_col,
        }
    )
    if ref_age_col is None:
        diagnostics["reference_error"] = "No published mean-age column was detected."
        return pd.DataFrame(columns=VALIDATION_COLUMNS), diagnostics

    working = pd.DataFrame()
    working["site_id"] = merged[site_id_col] if site_id_col else merged[join_col]
    working["sample_id"] = merged[sample_id_col]
    working["sample_year"] = (
        merged[sample_year_col].map(decimal_year) if sample_year_col else 2011.0
    )
    working["reference_mean_age_years"] = to_numeric(merged[ref_age_col])
    working["tritium_TU"] = to_numeric(merged[tritium_col]) if tritium_col else np.nan
    working["tritium_sigma_TU"] = to_numeric(merged[tritium_sigma_col]) if tritium_sigma_col else np.nan
    working["he3_trit_TU"] = to_numeric(merged[he3_col]) if he3_col else np.nan
    working["he3_trit_sigma_TU"] = to_numeric(merged[he3_sigma_col]) if he3_sigma_col else np.nan
    working["sf6_pptv"] = to_numeric(merged[sf6_col]) if sf6_col else np.nan
    working["sf6_sigma_pptv"] = to_numeric(merged[sf6_sigma_col]) if sf6_sigma_col else np.nan
    working["cfc11_pptv"] = to_numeric(merged[cfc11_col]) if cfc11_col else np.nan
    working["cfc11_sigma_pptv"] = to_numeric(merged[cfc11_sigma_col]) if cfc11_sigma_col else np.nan
    working["cfc12_pptv"] = to_numeric(merged[cfc12_col]) if cfc12_col else np.nan
    working["cfc12_sigma_pptv"] = to_numeric(merged[cfc12_sigma_col]) if cfc12_sigma_col else np.nan
    working["cfc113_pptv"] = to_numeric(merged[cfc113_col]) if cfc113_col else np.nan
    working["cfc113_sigma_pptv"] = to_numeric(merged[cfc113_sigma_col]) if cfc113_sigma_col else np.nan
    working["c14_pmc"] = to_numeric(merged[c14_col]) if c14_col else np.nan
    working["c14_sigma_pmc"] = to_numeric(merged[c14_sigma_col]) if c14_sigma_col else np.nan
    working["he4_ccpg"] = to_numeric(merged[he4_col]) if he4_col else np.nan

    fraction_specs = {
        "reference_anthropocene_fraction": [r"anthropocene", r"young", r"modern"],
        "reference_holocene_fraction": [r"holocene"],
        "reference_pleistocene_fraction": [r"pleistocene"],
    }
    for output_col, patterns in fraction_specs.items():
        source_col = first_matching_column(merged, patterns)
        working[output_col] = to_numeric(merged[source_col]) if source_col else np.nan

    inferred = working.apply(infer_sample_age, axis=1, result_type="expand")
    validation = pd.concat([working, inferred], axis=1)
    validation = validation[
        pd.notna(validation["reference_mean_age_years"])
        & pd.notna(validation["hydrosheaf_age_years"])
        & (validation["reference_mean_age_years"] > 0)
        & (validation["hydrosheaf_age_years"] > 0)
    ].copy()

    validation["log10_error"] = np.log10(validation["hydrosheaf_age_years"]) - np.log10(
        validation["reference_mean_age_years"]
    )
    validation["inside_hydrosheaf_ci"] = (
        (validation["reference_mean_age_years"] >= validation["hydrosheaf_ci_low_years"])
        & (validation["reference_mean_age_years"] <= validation["hydrosheaf_ci_high_years"])
    )

    return validation.reindex(columns=VALIDATION_COLUMNS), diagnostics


def summarize_validation(validation: pd.DataFrame, diagnostics: Dict[str, Any]) -> pd.DataFrame:
    if validation.empty:
        return pd.DataFrame(
            [
                {
                    "status": "blocked_no_validation_rows",
                    "n_samples": 0,
                    "log10_age_rmse": np.nan,
                    "median_log10_bias": np.nan,
                    "ci_coverage_fraction": np.nan,
                    "age_table": diagnostics.get("age_table"),
                    "tracer_table": diagnostics.get("tracer_table"),
                    "note": diagnostics.get("join_error")
                    or diagnostics.get("reference_error")
                    or "No readable source tables with both reference ages and supported tracers were available.",
                }
            ]
        )

    grouped = []
    for label, frame in [("all", validation)] + list(validation.groupby("hydrosheaf_method")):
        errors = frame["log10_error"].dropna()
        grouped.append(
            {
                "status": "completed",
                "group": label,
                "n_samples": int(len(frame)),
                "log10_age_rmse": float(np.sqrt(np.mean(errors**2))) if len(errors) else np.nan,
                "median_log10_bias": float(np.median(errors)) if len(errors) else np.nan,
                "median_abs_log10_error": float(np.median(np.abs(errors))) if len(errors) else np.nan,
                "ci_coverage_fraction": float(frame["inside_hydrosheaf_ci"].mean()),
                "age_table": diagnostics.get("age_table"),
                "tracer_table": diagnostics.get("tracer_table"),
                "note": "Hydrosheaf estimate uses upgraded multi-tracer support: 3H, 3H/3He, 14C, SF6, CFC-11, CFC-12, and CFC-113 where available.",
            }
        )
    return pd.DataFrame(grouped)


def plot_figure4a(validation: pd.DataFrame) -> Optional[Path]:
    if validation.empty or plt is None:
        return None
    fig_path = FIGURE_DIR / "figure4a_public_tracer_age_agreement.png"
    fig, ax = plt.subplots(figsize=(7, 6))
    def plot_group(method: Any) -> str:
        text = str(method)
        if "3H/3He" in text:
            return "includes 3H/3He"
        if "SF6" in text or "CFC" in text:
            return "includes SF6/CFC"
        if "14C" in text:
            return "14C only"
        if "3H" in text:
            return "3H only"
        return "other"

    plot_data = validation.copy()
    plot_data["plot_group"] = plot_data["hydrosheaf_method"].map(plot_group)
    groups = ["14C only", "3H only", "includes 3H/3He", "includes SF6/CFC", "other"]
    colors = {
        "14C only": "#dc2626",
        "3H only": "#2563eb",
        "includes 3H/3He": "#16a34a",
        "includes SF6/CFC": "#7c3aed",
        "other": "#475569",
    }
    for group in groups:
        subset = plot_data[plot_data["plot_group"] == group]
        if subset.empty:
            continue
        ax.scatter(
            subset["reference_mean_age_years"],
            subset["hydrosheaf_age_years"],
            s=24,
            alpha=0.72,
            label=group,
            color=colors[group],
            edgecolor="white",
            linewidth=0.35,
        )
    low = max(0.1, min(validation["reference_mean_age_years"].min(), validation["hydrosheaf_age_years"].min()) * 0.5)
    high = max(validation["reference_mean_age_years"].max(), validation["hydrosheaf_age_years"].max()) * 2.0
    grid = np.logspace(np.log10(low), np.log10(high), 200)
    ax.plot(grid, grid, color="#0f172a", linewidth=1.2, label="1:1")
    ax.fill_between(grid, grid / 10.0, grid * 10.0, color="#94a3b8", alpha=0.18, label="+/- 1 log10 year")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_xlabel("Published TracerLPM mean age (years)")
    ax.set_ylabel("Hydrosheaf single-node age estimate (years)")
    ax.set_title("Figure 4A. Public tracer-age validation")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)
    return fig_path


def write_failure_report(
    source_attempts: List[SourceAttempt],
    diagnostics: Dict[str, Any],
    summary: pd.DataFrame,
) -> Path:
    report_path = RESULT_DIR / "e1_failure_report.md"
    attempts_md = "\n".join(
        f"- `{a.status}` {a.url} -> {a.target or 'n/a'}: {a.detail}"
        for a in source_attempts
    )
    candidate_count = len(diagnostics.get("candidate_tables", []))
    note = str(summary.iloc[0].get("note", "")) if not summary.empty else ""
    report = [
        "# E1 USGS Public Tracer-Age Validation Report",
        "",
        f"Run timestamp UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Status",
        "",
        "Blocked. No validation metrics or Figure 4A were generated because the actual USGS source tables could not be read from this machine.",
        "",
        "## Source Confirmed",
        "",
        "- Data release DOI: `10.5066/P9W7T0DN`.",
        f"- Catalog URL: `{CATALOG_URL}`.",
        "- The accessible USGS metadata confirms seven tables containing site information, published mean age and age fractions, environmental tracer concentrations, and TracerLPM calibration details for 1,279 public-supply aquifer sites.",
        "",
        "## Download Attempts",
        "",
        attempts_md or "- No download attempts were made.",
        "",
        "## Local Data Discovery",
        "",
        f"- Candidate source tables found in `external/usgs_age/input`: {candidate_count}.",
        f"- Diagnostic note: {note}",
        "",
        "## What Is Needed To Complete E1",
        "",
        "1. Retry the script when ScienceBase is reachable, or manually download the data release from the DOI page.",
        "2. Place the seven USGS release tables in `external/usgs_age/input` as CSV/XLSX/TXT files.",
        "3. Re-run `python M2\\m2_benchmark\\scripts\\run_e1_usgs_age_validation.py --no-download`.",
        "",
        "The script will then generate `usgs_age_validation.csv`, `usgs_age_validation_summary.csv`, and `figures/figure4a_public_tracer_age_agreement.png` from the actual public tables.",
        "",
    ]
    report_path.write_text("\n".join(report), encoding="utf-8")
    return report_path


def write_success_report(
    validation: pd.DataFrame,
    summary: pd.DataFrame,
    diagnostics: Dict[str, Any],
    figure_path: Optional[Path],
) -> Path:
    report_path = RESULT_DIR / "e1_validation_report.md"
    summary_md = summary.to_markdown(index=False)
    methods = validation["hydrosheaf_method"].value_counts().to_dict()
    method_md = "\n".join(f"- `{method}`: {count} samples" for method, count in methods.items())
    report = [
        "# E1 USGS Public Tracer-Age Validation Report",
        "",
        f"Run timestamp UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Status",
        "",
        "Completed from the local USGS release files supplied for DOI `10.5066/P9W7T0DN`.",
        "",
        "## Source Tables Used",
        "",
        f"- Age/reference table: `{diagnostics.get('age_table')}`",
        f"- Tracer table: `{diagnostics.get('tracer_table')}`",
        f"- Site table: `{diagnostics.get('site_table')}`",
        "",
        "## Outputs",
        "",
        f"- Validation rows: `{RESULT_DIR / 'usgs_age_validation.csv'}`",
        f"- Summary metrics: `{RESULT_DIR / 'usgs_age_validation_summary.csv'}`",
        f"- Figure 4A: `{figure_path}`",
        "",
        "## Sample Counts By Hydrosheaf Tracer Method",
        "",
        method_md,
        "",
        "## Metrics",
        "",
        summary_md,
        "",
        "## Interpretation Guardrail",
        "",
        "This E1 run is a public-data agreement test for Hydrosheaf's upgraded single-node multi-tracer age inference. The implementation now uses 3H/3He closed-system ingrowth and SF6/CFC atmospheric-equivalent screening ages in addition to 3H and 14C. The USGS reference ages were fitted with the full TracerLPM workflow and site-specific corrections, so remaining disagreement should be interpreted as evidence for future work on local recharge histories, dissolved-gas corrections, and mixture-model fitting.",
        "",
    ]
    report_path.write_text("\n".join(report), encoding="utf-8")
    return report_path


def update_external_status(
    status: str,
    note: str,
    performance_metric: Optional[str] = None,
) -> None:
    workplan = BENCHMARK_ROOT / "tables" / "external_validation_workplan.csv"
    if not workplan.exists():
        pass
    else:
        df = pd.read_csv(workplan)
        if "status" in df.columns:
            mask = df["validation_tier"] == "E1_public_tracer_age"
            if mask.any():
                df.loc[mask, "status"] = status
                if "note" not in df.columns:
                    df["note"] = ""
                df.loc[mask, "note"] = note
                df.to_csv(workplan, index=False)

    table4 = BENCHMARK_ROOT / "tables" / "table4_validation_design_and_results.csv"
    if table4.exists():
        df = pd.read_csv(table4)
        if "m2_status" in df.columns:
            mask = df["benchmark"] == "Public tracer-age validation"
            if mask.any():
                df.loc[mask, "m2_status"] = status
                if performance_metric:
                    df.loc[mask, "performance_metric"] = performance_metric
                if "run_note" not in df.columns:
                    df["run_note"] = ""
                df.loc[mask, "run_note"] = note
                df.to_csv(table4, index=False)

    table_s4 = BENCHMARK_ROOT / "tables" / "table_s4_benchmark_dataset_inventory.csv"
    if table_s4.exists():
        df = pd.read_csv(table_s4)
        if "status" in df.columns:
            mask = df["resource"] == "USGS public-supply groundwater-age data release"
            if mask.any():
                df.loc[mask, "status"] = status
                if "note" not in df.columns:
                    df["note"] = ""
                df.loc[mask, "note"] = note
                df.to_csv(table_s4, index=False)


def write_manifest(
    source_attempts: List[SourceAttempt],
    diagnostics: Dict[str, Any],
    validation: pd.DataFrame,
    summary: pd.DataFrame,
    figure_path: Optional[Path],
) -> Path:
    def json_safe(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(k): json_safe(v) for k, v in value.items()}
        if isinstance(value, list):
            return [json_safe(v) for v in value]
        if isinstance(value, tuple):
            return [json_safe(v) for v in value]
        if isinstance(value, (np.integer,)):
            return int(value)
        if isinstance(value, (np.floating, float)):
            number = float(value)
            return number if math.isfinite(number) else None
        if pd.isna(value):
            return None
        return value

    manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "validation_tier": "E1_public_tracer_age",
        "source_doi": "10.5066/P9W7T0DN",
        "source_catalog_url": CATALOG_URL,
        "source_attempts": [attempt.__dict__ for attempt in source_attempts],
        "diagnostics": diagnostics,
        "outputs": {
            "validation_csv": str(RESULT_DIR / "usgs_age_validation.csv"),
            "summary_csv": str(RESULT_DIR / "usgs_age_validation_summary.csv"),
            "figure4a": str(figure_path) if figure_path else None,
            "success_report": str(RESULT_DIR / "e1_validation_report.md"),
            "failure_report": str(RESULT_DIR / "e1_failure_report.md"),
        },
        "n_validation_rows": int(len(validation)),
        "summary": summary.to_dict(orient="records"),
    }
    out = RESULT_DIR / "e1_usgs_age_source_manifest.json"
    out.write_text(json.dumps(json_safe(manifest), indent=2, allow_nan=False), encoding="utf-8")
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Use only files already present in external/usgs_age/input.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_dirs()
    source_attempts: List[SourceAttempt] = []
    if not args.no_download:
        source_attempts = download_sources()

    try:
        candidates = discover_tables(INPUT_DIR)
        validation, diagnostics = build_validation_dataset(candidates)
        validation.to_csv(RESULT_DIR / "usgs_age_validation.csv", index=False)
        summary = summarize_validation(validation, diagnostics)
        summary.to_csv(RESULT_DIR / "usgs_age_validation_summary.csv", index=False)
        figure_path = plot_figure4a(validation)
        write_manifest(source_attempts, diagnostics, validation, summary, figure_path)

        if validation.empty:
            report_path = write_failure_report(source_attempts, diagnostics, summary)
            update_external_status(
                "blocked_source_tables_unavailable",
                "Source metadata reached, but ScienceBase/source data tables were not readable locally.",
            )
            print(f"E1 blocked; report written to {report_path}")
            return 0

        failure_report = RESULT_DIR / "e1_failure_report.md"
        if failure_report.exists():
            failure_report.unlink()
        success_report = write_success_report(validation, summary, diagnostics, figure_path)
        all_row = summary[summary["group"] == "all"].iloc[0]
        metric_text = (
            f"n={int(all_row['n_samples'])}; "
            f"log10 age RMSE={float(all_row['log10_age_rmse']):.3f}; "
            f"median log10 bias={float(all_row['median_log10_bias']):.3f}; "
            f"CI coverage={float(all_row['ci_coverage_fraction']):.3f}"
        )
        update_external_status(
            "completed",
            f"Generated {len(validation)} public tracer-age validation rows. Report: {success_report.name}.",
            metric_text,
        )
        print(f"E1 completed with {len(validation)} validation rows.")
        if figure_path:
            print(f"Figure 4A written to {figure_path}")
        return 0
    except Exception as exc:
        validation = pd.DataFrame(columns=VALIDATION_COLUMNS)
        validation.to_csv(RESULT_DIR / "usgs_age_validation.csv", index=False)
        summary = pd.DataFrame(
            [
                {
                    "status": "failed_exception",
                    "n_samples": 0,
                    "log10_age_rmse": np.nan,
                    "median_log10_bias": np.nan,
                    "ci_coverage_fraction": np.nan,
                    "note": repr(exc),
                }
            ]
        )
        summary.to_csv(RESULT_DIR / "usgs_age_validation_summary.csv", index=False)
        diagnostics = {"exception": repr(exc), "traceback": traceback.format_exc()}
        write_manifest(source_attempts, diagnostics, validation, summary, None)
        report_path = write_failure_report(source_attempts, diagnostics, summary)
        update_external_status("failed_exception", repr(exc))
        print(f"E1 failed; report written to {report_path}")
        print(traceback.format_exc())
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
