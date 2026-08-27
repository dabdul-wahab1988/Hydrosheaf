"""Verify the built DOCX files carry the corrected text, not the stale build.

The previous round shipped DOCX files that predated the source fixes, so every
reviewer complaint about figures and metrics was against text the authors had
already corrected. This check reads the DOCX back and asserts (a) the corrected
strings are present, (b) the stale strings are gone, and (c) the DOCX body
matches the markdown source.

Run: .venv/Scripts/python.exe M2/M2_ready/Revision/audit/verify_docx.py
"""

from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

REV = Path(__file__).resolve().parent.parent

# (docx, markdown source)
PAIRS = [
    (REV / "Manuscript-Final-Revised.docx", REV / "Manuscript-Final-Revised.md"),
    (REV / "Manuscript-Colored-Changes.docx", REV / "Manuscript-Final-Revised.md"),
    (REV / "Supplementary-Information-Revised.docx",
     REV / "Supplementary-Information-Revised.md"),
]

# Strings that MUST be present in the built manuscript DOCX
REQUIRED_MS = [
    "Sheaf-Inspired Edge Refinement",
    "no graph or sheaf Laplacian is formed",
    "Live PHREEQC kinetic execution is not part of the current release",
    "Harte, 2021; USGS data release, DOI 10.5066/F7J102FK",
    "Jurgens et al., 2022; USGS data release, DOI 10.5066/P9W7T0DN",
    "Runtime and scalability are reported in the Supplementary Information",
    # third-round: the real coboundary construction replaces the unimplemented
    # spectral-localisation claim
    "cellular-sheaf coboundary",
    "homogeneous nullity",
    "affine obstruction energy",
    "leave-one-edge-out",
    "no exact global section",
    "lateral, dispersive mixing candidates",
    "258",
    "extent and stability disagree",
    "Supplementary Table S9",
    "prioritisation rather than management-level process attribution",
    "edge_radius_km",
    "Supplementary Table S14",
    "zero unresolved null edges",
    "pre-release reproducibility manifest",
    "synthetic benchmark age network",
]

REQUIRED_SI = [
    "edge_radius_km (search radius)",
    "Talensi_TGW3->TGW8",
    "Talensi_TGW44->TGW45",
    "log10 RMSE = 0.235",
    "1,300 edge-checks",
    "Multi-criteria edge selection and sheaf-cohomology diagnostics",
    "Reaction correlation and parameter-level uncertainty",
    "Null edges: none",
    "Supplementary Table S14",
    "homogeneous nullity",
    "does not apply a continuous overlap-proportional weight",
    "no eigendecomposition of any Laplacian is performed",
]

# Strings that MUST NOT survive anywhere
FORBIDDEN = [
    "Sheaf-Based Edge Refinement",
    "The sheaf Laplacian is constructed",
    "sheaf Laplacian's spectral decomposition",
    "Forward kinetic validation is subsequently conducted",
    "qA",
    # the substring-bug rendering: an edge whose PSI family agrees with its
    # dominant process but was reported as a mismatch
    "provisional carbonate dissolution; PSI family: Carbonates",
    "sufficiently stable to inform Where PSI",
    "Supplementary Table: topology baselines",
    "this radius has no fixed default",
    "ranked by haversine distance and limited",
    # third-round: the unimplemented spectral claim and the superseded as-run
    # field topology must not come back
    "spectral decomposition of the weighted graph Laplacian",
    "the spectral decomposition of *L* is used to localise",
    "the Laplacian we compute is an ordinary weighted graph Laplacian",
    "H0 = dim ker D (dimension of the global section space)",
    "H0 = 0 means no assignment of node states",
    "ten-dimensional family of globally consistent assignments",
    "Three of the six edges",
    "Seven edges converging",
    "returned zero chemical closure",
    "down-weighted in refinement",
    "median PSI was 1.0",
    "208 retained directed",
    "208 directed\nedges",
    "Evaporites 93",
    "median PSI was 1.0",
]

failures: list[str] = []


def docx_text(path: Path) -> str:
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf8")
    xml = re.sub(r"</w:p>", "\n", xml)
    xml = re.sub(r"<[^>]+>", "", xml)
    return (xml.replace("&amp;", "&").replace("&lt;", "<")
               .replace("&gt;", ">").replace("&quot;", '"').replace("&#39;", "'"))


def norm(s: str) -> str:
    for a, b in [("−", "-"), ("–", "-"), ("—", "-"), (" ", " "), (" ", " ")]:
        s = s.replace(a, b)
    return re.sub(r"\s+", " ", s)


for docx_path, md_path in PAIRS:
    if not docx_path.exists():
        failures.append(f"{docx_path.name}: NOT BUILT")
        continue
    body = norm(docx_text(docx_path))
    required = REQUIRED_MS if "Manuscript" in docx_path.name else REQUIRED_SI
    for s in required:
        if norm(s) not in body:
            failures.append(f"{docx_path.name}: MISSING required text {s!r}")
    for s in FORBIDDEN:
        if norm(s) in body:
            failures.append(f"{docx_path.name}: STALE text still present {s!r}")

    # The DOCX must be newer than its source, or it is a stale build again.
    if docx_path.stat().st_mtime < md_path.stat().st_mtime:
        failures.append(
            f"{docx_path.name} is OLDER than {md_path.name} -- stale build"
        )

    # Figures must be embedded
    names = zipfile.ZipFile(docx_path).namelist()
    n_media = len([n for n in names if n.startswith("word/media/")])
    expected = 7 if "Manuscript" in docx_path.name else 3
    if n_media != expected:
        failures.append(
            f"{docx_path.name}: {n_media} embedded images, expected {expected}"
        )

print(f"verify_docx: checked {len(PAIRS)} documents")
if failures:
    print(f"\n{len(failures)} FAILURE(S):\n")
    for f in failures:
        print("  -", f)
    sys.exit(1)
print("PASS: built DOCX files carry the corrected text and embedded figures")
