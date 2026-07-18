"""Shared engine for the M6 Hydrosheaf field-transfer & robustness benchmark.

Reuses the frozen M5 inverse-reaction primitives and the frozen M5 Mechanism
Resolution Score (MRS) calibration. Nothing here re-fits M5: M6 is a *transfer*
study, so the upstream calibration is applied unchanged to real Ghanaian data.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

# --- locate repo + reuse M5 primitives ---------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[3]
M5_SCRIPTS = REPO_ROOT / "M5" / "m5_inverse_reaction_benchmark" / "scripts"
M5_RESULTS = REPO_ROOT / "M5" / "m5_inverse_reaction_benchmark" / "results"
import sys
for p in (str(REPO_ROOT), str(M5_SCRIPTS)):
    if p not in sys.path:
        sys.path.insert(0, p)

import m5_common as m5  # noqa: E402  (fit_inverse, matrix_diagnostics, ...)

SEED = 1234
BENCH_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCH_DIR / "results"

DATA = {
    "northern_ghana": REPO_ROOT / "data" / "NorthenGhana" / "Aquifers_Dataset_Mendeley.xlsx",
    "talensi": REPO_ROOT / "data" / "Talensi_MiningArea" / "talensi.csv",
    "manu": REPO_ROOT / "data" / "LowerAnayari" / "manu.csv",
}

# Reaction-basis ions (subset of M5 ION_ORDER used by the linear panel)
MAJOR_IONS = ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3"]
CHARGE = {  # signed charge for CBE (meq/L)
    "Ca": +2, "Mg": +2, "Na": +1, "K": +1,
    "HCO3": -1, "Cl": -1, "SO4": -2, "NO3": -1, "F": -1,
}

# Evidence-tier ladder ---------------------------------------------------------
# reaction-panel ions per tier; evidence lifters listed separately
TIERS = {
    "tier0_majors":        {"ions": MAJOR_IONS,          "lifters": []},
    "tier1_isotopes":      {"ions": MAJOR_IONS,          "lifters": ["iso"]},
    "tier2_fluoride":      {"ions": MAJOR_IONS + ["F"],  "lifters": ["iso"]},
    "tier3_sr_sio2":       {"ions": MAJOR_IONS + ["F"],  "lifters": ["iso", "sr_sio2"]},
    "tier4_full_metadata": {"ions": MAJOR_IONS + ["F"],  "lifters": ["iso", "sr_sio2", "si", "meta", "edges"]},
}
TIER_ORDER = list(TIERS)

# Reaction family -> interpretable process label (for concordance, NOT truth) --
FAMILY_OF = {r.label: r.family for r in m5.REACTIONS}
PROCESS_LABEL = {
    "silicate": "Silicate weathering / recharge mixing",
    "carbonate": "Carbonate weathering",
    "ion_exchange": "Cation exchange / Na enrichment",
    "evaporite": "Evapoconcentration / salinity",
    "anthropogenic": "Agro-domestic nitrate loading",
    "redox": "Redox transformation",
    "trace_mineral": "Fluoride mobilisation",
}

SI_COL_TO_PHASE = {"Calcite_SI": "Calcite", "Dolomite_SI": "Dolomite",
                   "Gypsum_SI": "Gypsum", "Halite_SI": "Halite"}


# --- data harmonisation -------------------------------------------------------
def _to_mmol(series_mgL: pd.Series, ion: str) -> pd.Series:
    return pd.to_numeric(series_mgL, errors="coerce") / m5.MOLAR_MASS_G_MOL[ion]


def load_northern_ghana() -> pd.DataFrame:
    """320 seasonal samples with full metadata + SI (Tier 4 capable)."""
    f = DATA["northern_ghana"]
    hs = pd.read_excel(f, sheet_name="Hydrochemistry_Seasonal")
    wn = pd.read_excel(f, sheet_name="Wells_Nodes")[
        ["Well_ID", "Aquifer_Type", "Geology_Group", "Lithology", "Land_Use",
         "Latitude", "Longitude"]
    ]
    df = hs.merge(wn, on="Well_ID", how="left", suffixes=("", "_node"))
    out = pd.DataFrame({"dataset": "northern_ghana", "sample_id": df["Sample_ID"],
                        "site_id": df["Well_ID"], "season": df["Season"]})
    for ion in MAJOR_IONS + ["F"]:
        col = f"{ion}_mg_L"
        out[ion] = _to_mmol(df[col], ion) if col in df else np.nan
    out["Fe"] = np.nan
    out["PO4"] = np.nan
    out["Sr_mgL"] = pd.to_numeric(df.get("Sr_mg_L"), errors="coerce")
    out["SiO2_mgL"] = pd.to_numeric(df.get("SiO2_mg_L"), errors="coerce")
    out["d18O"] = pd.to_numeric(df.get("d18O_permil"), errors="coerce")
    out["d2H"] = pd.to_numeric(df.get("d2H_permil"), errors="coerce")
    for c in SI_COL_TO_PHASE:
        out[c] = pd.to_numeric(df.get(c), errors="coerce")
    out["Aquifer_Type"] = df["Aquifer_Type"]
    out["Geology_Group"] = df.get("Geology_Group")
    out["Lithology"] = df.get("Lithology")
    out["Latitude"] = pd.to_numeric(df.get("Latitude"), errors="coerce")
    out["Longitude"] = pd.to_numeric(df.get("Longitude"), errors="coerce")
    out["Data_Class"] = df.get("Data_Class")
    out["CBE_provided"] = pd.to_numeric(df.get("Charge_Balance_Error_pct"), errors="coerce")
    out["Dominant_Process_prior"] = df.get("Dominant_Process")
    out["Facies"] = df.get("Hydrochemical_Facies")
    return out


def load_talensi() -> pd.DataFrame:
    """63 samples, majors + Fe + isotopes; no F/Sr/SiO2/season (Tier 1)."""
    df = pd.read_csv(DATA["talensi"])
    out = pd.DataFrame({"dataset": "talensi", "sample_id": df["Code"],
                        "site_id": df["Code"], "season": "unknown"})
    m = {"HCO3": "HCO3", "Na": "Na", "K": "K", "Ca": "Ca", "Mg": "Mg",
         "Cl": "Cl", "SO4": "SO4", "NO3": "NO3"}
    for ion, col in m.items():
        out[ion] = _to_mmol(df[col], ion)
    out["F"] = np.nan
    out["Fe"] = _to_mmol(df["Fe"], "Fe")
    out["PO4"] = np.nan
    out["Sr_mgL"] = np.nan
    out["SiO2_mgL"] = np.nan
    out["d18O"] = pd.to_numeric(df["d18O"], errors="coerce")
    out["d2H"] = pd.to_numeric(df["d2H"], errors="coerce")
    for c in SI_COL_TO_PHASE:
        out[c] = np.nan
    out["Aquifer_Type"] = "Birimian basement (mining area)"
    out["Latitude"] = pd.to_numeric(df["Latitude"], errors="coerce")
    out["Longitude"] = pd.to_numeric(df["Longitude"], errors="coerce")
    out["Data_Class"] = np.nan
    out["Dominant_Process_prior"] = np.nan
    return out


def load_manu() -> pd.DataFrame:
    """41 samples, majors + F + Fe + isotopes; no Sr/SiO2/season (Tier 2)."""
    df = pd.read_csv(DATA["manu"])
    out = pd.DataFrame({"dataset": "manu", "sample_id": df["Sample ID"],
                        "site_id": df["Sample ID"], "season": "unknown"})
    m = {"Na": "Na", "K": "K", "Ca": "Ca", "Mg": "Mg", "F": "F",
         "Cl": "Cl", "HCO3": "HCO3", "NO3": "NO3", "SO4": "SO4"}
    for ion, col in m.items():
        out[ion] = _to_mmol(df[col], ion)
    out["Fe"] = _to_mmol(df["Fe"], "Fe")
    out["PO4"] = np.nan
    out["Sr_mgL"] = np.nan
    out["SiO2_mgL"] = np.nan
    out["d18O"] = pd.to_numeric(df["d18O"], errors="coerce")
    out["d2H"] = pd.to_numeric(df["d2H"], errors="coerce")
    for c in SI_COL_TO_PHASE:
        out[c] = np.nan
    out["Aquifer_Type"] = "Regolith/alluvial shallow (Lower Anayari)"
    out["Latitude"] = pd.to_numeric(df["Y coordinate"], errors="coerce")
    out["Longitude"] = pd.to_numeric(df["X coordinate"], errors="coerce")
    out["Data_Class"] = np.nan
    out["Dominant_Process_prior"] = np.nan
    return out


def load_all() -> dict[str, pd.DataFrame]:
    return {"northern_ghana": load_northern_ghana(),
            "talensi": load_talensi(), "manu": load_manu()}


# --- quality control ----------------------------------------------------------
def charge_balance_error(row: Mapping[str, float]) -> float:
    """Independent CBE (%) in meq/L from harmonised mmol/L ions."""
    cat = sum(max(row.get(i, 0) or 0, 0) * CHARGE[i] for i in ["Ca", "Mg", "Na", "K"])
    an = sum(abs(row.get(i, 0) or 0) * abs(CHARGE[i])
             for i in ["HCO3", "Cl", "SO4", "NO3", "F"] if pd.notna(row.get(i)))
    if cat + an == 0:
        return np.nan
    return 100.0 * (cat - an) / (cat + an)


def cbe_class(cbe: float) -> str:
    a = abs(cbe)
    if not np.isfinite(cbe):
        return "unusable"
    if a <= 5:
        return "quantitative"
    if a <= 10:
        return "screening"
    return "exploratory"


# --- inference core -----------------------------------------------------------
def conservative_factor(x_source: Mapping[str, float],
                        x_target: Mapping[str, float]) -> float:
    """Evapoconcentration/dilution factor from conservative Cl (clipped)."""
    cs, ct = x_source.get("Cl", np.nan), x_target.get("Cl", np.nan)
    if not (np.isfinite(cs) and np.isfinite(ct)) or cs <= 0:
        return 1.0
    return float(np.clip(ct / cs, 0.2, 5.0))


def residual_vector(x_source: Mapping[str, float], x_target: Mapping[str, float],
                    ions: Sequence[str], transport_correct: bool = True) -> np.ndarray:
    """Reactive residual over an ion panel (mmol/L).

    With transport_correct, the conservative component is removed using Cl as a
    conservative tracer: x_transport = x_source * (Cl_target/Cl_source), so the
    residual isolates non-conservative water-rock reactions. Cl itself is the
    reference, so halite is intentionally Cl-confounded (documented limitation).
    """
    f = conservative_factor(x_source, x_target) if transport_correct else 1.0
    return np.array([(x_target.get(i, np.nan) or np.nan)
                     - f * (x_source.get(i, np.nan) or np.nan) for i in ions], dtype=float)


def _penalty_for_panel(ions: Sequence[str], high: float = 20.0) -> np.ndarray:
    """A reaction is eligible only if every ion it uses is measured in the panel.
    Ineligible reactions get a near-prohibitive L1 penalty so they cannot be
    spuriously selected (e.g. fluorite when F is not measured)."""
    ionset = set(ions)
    scales = np.ones(len(m5.REACTIONS), dtype=float)
    for idx, rxn in enumerate(m5.REACTIONS):
        used = {ion for ion, coef in rxn.stoich.items() if coef != 0}
        if not used.issubset(ionset):
            scales[idx] = high
    return scales


def _upstream_si(target_row: Mapping[str, float]) -> dict[str, float]:
    si = {}
    for col, phase in SI_COL_TO_PHASE.items():
        v = target_row.get(col)
        if v is not None and np.isfinite(v):
            si[phase] = float(v)
    return si


RESOLUTION_ORDER = ["non_identifiable", "partially_identifiable",
                    "equivalence_class", "identifiable"]
# families whose mechanism cannot be *reported* as resolved without an
# independent corroborating tracer being measured and consistent
TRACER_GATED = {"silicate", "carbonate", "evaporite", "trace_mineral"}


def evidence_corroborated(family: str, row: Mapping[str, float],
                          lifters: Sequence[str]) -> bool:
    """Is the dominant process family independently corroborated by an
    available evidence lifter? Families inferable from major ions alone
    (ion_exchange, anthropogenic, redox) return True by construction."""
    L = set(lifters)
    if family == "silicate":
        return "sr_sio2" in L and np.isfinite(row.get("SiO2_mgL", np.nan)) and row["SiO2_mgL"] >= 10.0
    if family == "carbonate":
        si_ok = "si" in L and np.isfinite(row.get("Calcite_SI", np.nan))
        sr_ok = "sr_sio2" in L and np.isfinite(row.get("Sr_mgL", np.nan)) and row["Sr_mgL"] >= 0.10
        return si_ok or sr_ok
    if family == "evaporite":
        if "iso" not in L:
            return False
        d18, d2 = row.get("d18O", np.nan), row.get("d2H", np.nan)
        return np.isfinite(d18) and np.isfinite(d2) and (d2 - 8.0 * d18) < 8.0
    if family == "trace_mineral":
        f = row.get("F", np.nan)
        return np.isfinite(f) and f > 0.05          # F must be measured (Tier 2+)
    return True


def evidence_lift(base_class: str, family: str, row: Mapping[str, float],
                  lifters: Sequence[str]) -> tuple[str, bool]:
    """Downgrade a tracer-gated attribution one resolution level when its
    corroborating tracer is unavailable or inconsistent."""
    corrob = evidence_corroborated(family, row, lifters)
    if family in TRACER_GATED and not corrob:
        idx = max(0, RESOLUTION_ORDER.index(base_class) - 1)
        return RESOLUTION_ORDER[idx], corrob
    return base_class, corrob


@dataclass
class UnitResult:
    support: set
    extents: np.ndarray
    rmse: float
    l1_norm: float
    n_ions: int
    dominant_family: str
    dominant_process: str
    rank_skill: float
    coherence_skill: float
    support_stability: float
    predictive_skill: float
    thermodynamic_skill: float
    heldout_rmse: float
    base_resolution_class: str
    resolution_class: str
    evidence_corroborated: bool
    mrs: float
    thermo_violations: int


class TransferClassifier:
    """Applies the FROZEN M5 MRS calibration to M6 field units (unchanged)."""

    def __init__(self, path: Path = M5_RESULTS / "mrs_calibration_model.json"):
        d = json.loads(Path(path).read_text())
        self.features = d["features"]
        self.means = np.array(d["feature_means"], float)
        self.scales = np.array(d["feature_scales"], float)
        self.classes = d["resolution_classes"]
        self.class_coef = np.array(d["class_coefficients"], float)      # (10,4)
        self.reliab = np.array(d["reliability_coefficients"], float)    # (10,)

    def classify(self, feats: Sequence[float], ambiguous: bool) -> tuple[str, float]:
        design = np.concatenate([[1.0], (np.asarray(feats, float) - self.means) / self.scales])
        cls = self.classes[int(np.argmax(design @ self.class_coef))]
        mrs = 100.0 * float(np.clip(design @ self.reliab, 0.0, 1.0))
        if cls == "identifiable" and ambiguous:
            cls = "equivalence_class"
        return cls, mrs


def fit_unit(x_source: Mapping[str, float], x_target: Mapping[str, float],
             ions: Sequence[str], si: Mapping[str, float] | None,
             clf: TransferClassifier, class_map: Mapping[str, str],
             rng: np.random.Generator, n_boot: int = 48,
             lifters: Sequence[str] = (), evidence_row: Mapping[str, float] | None = None,
             heldout: bool = True, apply_gate: bool = True,
             transport_correct: bool = True) -> UnitResult | None:
    r = residual_vector(x_source, x_target, ions, transport_correct=transport_correct)
    if not np.all(np.isfinite(r)):
        return None
    method = "thermo" if si else "nonneg"
    pen = _penalty_for_panel(ions)
    fit = m5.fit_inverse(r, ions, method=method, lambda_l1=0.01, lambda_l2=0.001,
                         upstream_si=si or {}, penalty_scales=pen)
    extents = np.asarray(fit["extents"], float)
    support = m5.support_from_extents(extents)
    diag = m5.matrix_diagnostics(ions)

    # bootstrap support stability (5% analytical noise)
    base = support
    jac = []
    for _ in range(n_boot):
        rn = r + rng.normal(0, 0.05 * np.maximum(np.abs(r), 1e-6))
        fb = m5.fit_inverse(rn, ions, method=method, lambda_l1=0.01, lambda_l2=0.001,
                            upstream_si=si or {}, penalty_scales=pen)
        jac.append(m5.jaccard(base, m5.support_from_extents(np.asarray(fb["extents"], float))))
    stability = float(np.mean(jac)) if jac else 1.0

    # leave-one-ion held-out prediction
    ho = []
    if heldout and len(ions) > 3:
        for k in range(len(ions)):
            keep = [j for j in range(len(ions)) if j != k]
            sub = [ions[j] for j in keep]
            fk = m5.fit_inverse(r[keep], sub, method=method, lambda_l1=0.01,
                                lambda_l2=0.001, upstream_si=si or {},
                                penalty_scales=_penalty_for_panel(sub))
            pred_full = m5.reaction_matrix(ions).T @ np.asarray(fk["extents"], float)
            ho.append(abs(pred_full[k] - r[k]))
    heldout_rmse = float(np.sqrt(np.mean(np.square(ho)))) if ho else float(fit["rmse"])

    rank_skill = diag["rank"] / max(1.0, min(diag["n_ions"], diag["n_equivalence_classes"]))
    coherence_skill = 1.0 - min(1.0, diag["class_collapsed_coherence"])
    predictive_skill = math.exp(-np.mean(ho) / (0.05 + float(np.mean(np.abs(r))))) if ho else 0.5
    thermo_skill = 1.0 - min(1.0, float(fit["thermodynamic_violations"]) / max(len(support), 1))

    feats = [rank_skill, coherence_skill, stability, predictive_skill, thermo_skill,
             float(fit["rmse"]), heldout_rmse, len(ions), float(fit["l1_norm"])]
    ambiguous = any(sum(1 for v in class_map.values() if v == class_map[s]) > 1 for s in support)
    base_cls, mrs = clf.classify(feats, ambiguous)

    fam_extent: dict[str, float] = {}
    for lab in support:
        fam_extent[FAMILY_OF[lab]] = fam_extent.get(FAMILY_OF[lab], 0.0) + abs(
            extents[m5.REACTION_LABELS.index(lab)])
    dom_fam = max(fam_extent, key=fam_extent.get) if fam_extent else "none"

    erow = evidence_row if evidence_row is not None else x_target
    if apply_gate:
        cls, corrob = evidence_lift(base_cls, dom_fam, erow, lifters)
    else:
        cls, corrob = base_cls, evidence_corroborated(dom_fam, erow, lifters)
    return UnitResult(
        support=support, extents=extents, rmse=float(fit["rmse"]),
        l1_norm=float(fit["l1_norm"]), n_ions=len(ions), dominant_family=dom_fam,
        dominant_process=PROCESS_LABEL.get(dom_fam, "unresolved"),
        rank_skill=rank_skill, coherence_skill=coherence_skill,
        support_stability=stability, predictive_skill=predictive_skill,
        thermodynamic_skill=thermo_skill, heldout_rmse=heldout_rmse,
        base_resolution_class=base_cls, resolution_class=cls,
        evidence_corroborated=bool(corrob), mrs=mrs,
        thermo_violations=int(fit["thermodynamic_violations"]))


# --- edge-set generators ------------------------------------------------------
def _ionic_load(row: Mapping[str, float]) -> float:
    return float(np.nansum([abs(row.get(i, 0) or 0) * abs(CHARGE.get(i, 1))
                            for i in MAJOR_IONS]))


def chemistry_knn_edges(samples: pd.DataFrame, k: int = 3) -> list[tuple]:
    """kNN by standardised major-ion chemistry; direction by increasing ionic load."""
    feats = samples[MAJOR_IONS].apply(pd.to_numeric, errors="coerce")
    feats = feats.fillna(feats.mean())
    z = (feats - feats.mean()) / feats.std(ddof=0).replace(0, 1)
    ids = samples["sample_id"].tolist()
    load = {ids[i]: _ionic_load(samples.iloc[i]) for i in range(len(samples))}
    edges = []
    Z = z.to_numpy()
    for i in range(len(ids)):
        d = np.sqrt(((Z - Z[i]) ** 2).sum(1))
        order = np.argsort(d)
        for j in order[1:k + 1]:
            a, b = ids[i], ids[int(j)]
            src, tgt = (a, b) if load[a] <= load[b] else (b, a)
            edges.append((src, tgt))
    return sorted(set(edges))


def geographic_edges(samples: pd.DataFrame, k: int = 3) -> list[tuple]:
    lat = pd.to_numeric(samples["Latitude"], errors="coerce").to_numpy()
    lon = pd.to_numeric(samples["Longitude"], errors="coerce").to_numpy()
    ids = samples["sample_id"].tolist()
    load = {ids[i]: _ionic_load(samples.iloc[i]) for i in range(len(samples))}
    edges = []
    for i in range(len(ids)):
        d = np.sqrt((lat - lat[i]) ** 2 + (lon - lon[i]) ** 2)
        order = np.argsort(np.where(np.isfinite(d), d, np.inf))
        for j in order[1:k + 1]:
            a, b = ids[i], ids[int(j)]
            src, tgt = (a, b) if load[a] <= load[b] else (b, a)
            edges.append((src, tgt))
    return sorted(set(edges))


def random_edges(samples: pd.DataFrame, n: int, rng: np.random.Generator) -> list[tuple]:
    ids = samples["sample_id"].tolist()
    edges = set()
    while len(edges) < min(n, len(ids) * (len(ids) - 1)):
        a, b = rng.choice(ids, 2, replace=False)
        edges.add((a, b))
    return sorted(edges)


def get_class_map() -> dict[str, str]:
    class_map, _ = m5.equivalence_classes()
    return class_map


if __name__ == "__main__":
    data = load_all()
    for name, df in data.items():
        print(f"{name}: {df.shape[0]} samples, seasons={df['season'].unique().tolist()}")
