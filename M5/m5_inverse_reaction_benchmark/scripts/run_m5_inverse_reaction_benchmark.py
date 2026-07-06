"""Run the complete M5 identifiability-aware inverse-reaction benchmark."""
from __future__ import annotations

import argparse
import json
import math
import platform
import re
import sys
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hydrosheaf.isotopes import (  # noqa: E402
    compute_d_excess,
    evaporation_index,
)
from hydrosheaf.models.ec_tds import fit_linear_model  # noqa: E402
from hydrosheaf.models.evidence_lifted import (  # noqa: E402
    evidence_lifted_resolution,
)
from hydrosheaf.models.exchange import (  # noqa: E402
    classify_exchange_direction,
    compute_cai_indices,
)
from hydrosheaf.models.gibbs import compute_gibbs_metrics  # noqa: E402
from hydrosheaf.models.redox import classify_redox  # noqa: E402
from m5_common import (  # noqa: E402
    ION_ORDER,
    MOLAR_MASS_G_MOL,
    PHREEQC_TOTALS,
    REACTIONS,
    REACTION_BY_LABEL,
    REACTION_LABELS,
    class_set,
    equivalence_classes,
    finite_or_blank,
    find_phreeqc,
    fit_inverse,
    jaccard,
    matrix_diagnostics,
    precision_recall_f1,
    reaction_matrix,
    run_phreeqc,
    support_from_extents,
)


BENCHMARK_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BENCHMARK_DIR / "results"
TABLES_DIR = BENCHMARK_DIR / "tables"
DOCS_DIR = BENCHMARK_DIR / "docs"
PHREEQC_DIR = BENCHMARK_DIR / "phreeqc_inputs"
PHREEQC_INVERSE_DIR = PHREEQC_DIR / "inverse_baseline"
GHANA_WORKBOOK = REPO_ROOT / "data" / "NorthenGhana" / "Aquifers_Dataset_Mendeley.xlsx"
LEGACY_NORTHERN_GHANA_WORKBOOK = (
    REPO_ROOT / "data" / "NorthenGhana" / "NorthernGhana.xlsx"
)
TALENSI_CSV = REPO_ROOT / "data" / "Talensi_MiningArea" / "talensi.csv"
LOWER_ANAYARI_CSV = REPO_ROOT / "data" / "LowerAnayari" / "manu.csv"
RANDOM_SEED = 20250615
N_SCENARIOS_PER_ARCHETYPE = 60
NOISE_LEVELS = [0.0, 0.03, 0.08]
SUPPORT_THRESHOLD = 0.015
PHREEQC_INVERSE_UNCERTAINTY = 0.05
PHREEQC_INVERSE_FALLBACK_UNCERTAINTY = 0.20
PHREEQC_INVERSE_MAX_M5_PHASES = 10
DETECTION_THRESHOLD_GRID = [0.015, 0.020, 0.025, 0.030, 0.035, 0.040]
HYDROSHEAF_CORE_METHOD = "hydrosheaf_core"
GUARDED_METHOD = "hydrosheaf_guarded"
PRIMARY_METHOD = GUARDED_METHOD
LEGACY_PRIMARY_METHOD = "thermo_elastic_net"
SI_THRESHOLD = 0.1
CORE_EVIDENCE_MIN_SCALE = 0.35
CORE_EVIDENCE_MAX_SCALE = 3.00
DATA_TIER_PANEL = "full_11"
DATA_TIER_NOISE_LEVEL = 0.03
ION_CHARGES = {
    "Ca": 2.0,
    "Mg": 2.0,
    "Na": 1.0,
    "K": 1.0,
    "HCO3": 1.0,
    "Cl": 1.0,
    "SO4": 2.0,
    "NO3": 1.0,
    "F": 1.0,
    "Fe": 2.0,
    "PO4": 3.0,
}
DATA_TIER_DIAGNOSTICS = {
    "core": [],
    "plus_lite": ["SiO2", "Sr", "water_isotope_evaporation"],
    "enhanced": [
        "SiO2",
        "Sr",
        "water_isotope_evaporation",
        "Br",
        "DO",
        "DOC",
        "d34S_SO4",
        "d18O_SO4",
        "d15N_NO3",
        "d18O_NO3",
    ],
}
DIAGNOSTIC_SIGNATURES = {
    "calcite": {"Sr": 0.70},
    "dolomite": {"Sr": 0.55},
    "albite": {"SiO2": 0.85},
    "anorthite": {"SiO2": 1.00, "Sr": 0.20},
    "k_feldspar": {"SiO2": 0.75},
    "gypsum": {
        "Sr": 0.20,
        "water_isotope_evaporation": 0.15,
        "d34S_SO4": -0.65,
        "d18O_SO4": 0.75,
    },
    "halite": {"water_isotope_evaporation": 0.25, "Br": 1.00},
    "fluorite": {"Sr": 0.15},
    "apatite": {"Sr": 0.15},
    "pyrite_oxidation": {"DO": 0.80, "d34S_SO4": 0.90, "d18O_SO4": -0.45},
    "NO3src": {"d15N_NO3": -0.45, "d18O_NO3": 0.55},
    "denit": {"DOC": 0.80, "DO": -0.70, "d15N_NO3": 0.95, "d18O_NO3": 0.75},
    "CaNa_exch": {},
    "NaCa_exch": {},
    "MgNa_exch": {},
    "NaMg_exch": {},
}
MRS_FEATURES = [
    "rank_skill",
    "coherence_skill",
    "support_stability",
    "predictive_skill",
    "thermodynamic_skill",
    "reconstruction_rmse_mmolL",
    "heldout_rmse_mmolL",
    "n_measured_ions",
    "l1_norm",
]
RESOLUTION_CLASSES = [
    "non_identifiable",
    "partially_identifiable",
    "equivalence_class",
    "identifiable",
]

PANELS = {
    "full_11": list(ION_ORDER),
    "major_8": ["Ca", "Mg", "Na", "K", "HCO3", "Cl", "SO4", "NO3"],
    "no_alkalinity": [ion for ion in ION_ORDER if ion != "HCO3"],
    "no_redox_trace": [
        ion for ion in ION_ORDER if ion not in {"NO3", "Fe", "PO4"}
    ],
    "core_5": ["Ca", "Mg", "Na", "HCO3", "Cl"],
}

PHREEQC_INVERSE_PHASES = {
    "calcite": (
        "M5_calcite",
        ["M5_calcite", "    Ca(HCO3)2 = Ca+2 + 2 HCO3-", "    log_k 0.0"],
    ),
    "dolomite": (
        "M5_dolomite",
        ["M5_dolomite", "    CaMg(HCO3)4 = Ca+2 + Mg+2 + 4 HCO3-", "    log_k 0.0"],
    ),
    "albite": (
        "M5_albite",
        ["M5_albite", "    NaHCO3 = Na+ + HCO3-", "    log_k 0.0"],
    ),
    "anorthite": (
        "M5_anorthite",
        ["M5_anorthite", "    Ca(HCO3)2 = Ca+2 + 2 HCO3-", "    log_k 0.0"],
    ),
    "k_feldspar": (
        "M5_kfeldspar",
        ["M5_kfeldspar", "    KHCO3 = K+ + HCO3-", "    log_k 0.0"],
    ),
    "gypsum": (
        "M5_gypsum",
        ["M5_gypsum", "    CaSO4:2H2O = Ca+2 + SO4-2 + 2 H2O", "    log_k 0.0"],
    ),
    "halite": (
        "M5_halite",
        ["M5_halite", "    NaCl = Na+ + Cl-", "    log_k 0.0"],
    ),
    "fluorite": (
        "M5_fluorite",
        ["M5_fluorite", "    CaF2 = Ca+2 + 2 F-", "    log_k 0.0"],
    ),
    "apatite": (
        "M5_apatite",
        ["M5_apatite", "    Ca5(PO4)3F = 5 Ca+2 + 3 PO4-3 + F-", "    log_k 0.0"],
    ),
    "pyrite_oxidation": (
        "M5_pyriteox",
        [
            "M5_pyriteox",
            "    FeS2 + 3.5 O2 + H2O = Fe+2 + 2 SO4-2 + 2 H+",
            "    log_k 0.0",
        ],
    ),
    "NO3src": (
        "M5_NO3src",
        ["M5_NO3src", "    HNO3 = H+ + NO3-", "    log_k 0.0"],
    ),
    "denit": (
        "M5_denit",
        [
            "M5_denit",
            (
                "    CH2O + 0.8 NO3- = 0.4 N2 + 0.8 HCO3- "
                "+ 0.2 CO2 + 0.6 H2O"
            ),
            "    log_k 0.0",
        ],
    ),
}

PHREEQC_INVERSE_PHASE_TO_REACTION = {
    phase: label for label, (phase, _) in PHREEQC_INVERSE_PHASES.items()
}

ARCHETYPES = {
    "carbonate": {
        "pH": 7.15,
        "temp": 24.0,
        "Ca": 1.50,
        "Mg": 0.75,
        "Na": 0.80,
        "K": 0.08,
        "HCO3": 3.50,
        "Cl": 0.55,
        "SO4": 0.22,
        "NO3": 0.18,
        "F": 0.015,
        "Fe": 0.002,
        "PO4": 0.012,
    },
    "crystalline": {
        "pH": 6.75,
        "temp": 27.0,
        "Ca": 0.48,
        "Mg": 0.30,
        "Na": 1.25,
        "K": 0.19,
        "HCO3": 1.90,
        "Cl": 0.70,
        "SO4": 0.32,
        "NO3": 0.25,
        "F": 0.045,
        "Fe": 0.012,
        "PO4": 0.018,
    },
    "evaporitic": {
        "pH": 7.55,
        "temp": 30.0,
        "Ca": 2.20,
        "Mg": 1.30,
        "Na": 6.50,
        "K": 0.28,
        "HCO3": 2.40,
        "Cl": 7.40,
        "SO4": 2.20,
        "NO3": 0.38,
        "F": 0.055,
        "Fe": 0.002,
        "PO4": 0.010,
    },
    "mixed": {
        "pH": 7.05,
        "temp": 28.0,
        "Ca": 1.05,
        "Mg": 0.62,
        "Na": 2.70,
        "K": 0.16,
        "HCO3": 2.65,
        "Cl": 2.50,
        "SO4": 0.85,
        "NO3": 0.42,
        "F": 0.035,
        "Fe": 0.006,
        "PO4": 0.015,
    },
}

ARCHETYPE_POOLS = {
    "carbonate": [
        "calcite",
        "dolomite",
        "gypsum",
        "NO3src",
        "denit",
        "CaNa_exch",
    ],
    "crystalline": [
        "albite",
        "anorthite",
        "k_feldspar",
        "fluorite",
        "apatite",
        "pyrite_oxidation",
        "CaNa_exch",
        "NaCa_exch",
    ],
    "evaporitic": [
        "halite",
        "gypsum",
        "calcite",
        "NO3src",
        "CaNa_exch",
        "MgNa_exch",
    ],
    "mixed": list(REACTION_LABELS),
}

EXTENT_RANGES = {
    "apatite": (0.003, 0.015),
    "pyrite_oxidation": (0.008, 0.080),
    "NO3src": (0.030, 0.350),
    "denit": (0.025, 0.150),
    "CaNa_exch": (0.025, 0.220),
    "NaCa_exch": (0.025, 0.220),
    "MgNa_exch": (0.025, 0.180),
    "NaMg_exch": (0.025, 0.180),
}


def ensure_directories() -> None:
    for directory in (
        RESULTS_DIR,
        TABLES_DIR,
        DOCS_DIR,
        PHREEQC_DIR,
        PHREEQC_INVERSE_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def _finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _sample_from_vector(
    vector: Sequence[float],
    measured_ions: Sequence[str],
    *,
    equivalents: bool,
) -> dict[str, float]:
    measured = set(measured_ions)
    sample: dict[str, float] = {}
    tds = 0.0
    for ion in ION_ORDER:
        if ion not in measured:
            continue
        value = _finite_float(vector[ION_ORDER.index(ion)])
        if value is None:
            continue
        value = max(value, 0.0)
        tds += value * MOLAR_MASS_G_MOL[ion]
        sample[ion] = value * ION_CHARGES[ion] if equivalents else value
    sample["TDS"] = tds
    return sample


def _safe_gibbs_metrics(sample_meq: Mapping[str, float]) -> dict[str, object]:
    required = {"Na", "Ca", "Cl", "HCO3"}
    if not required.issubset(sample_meq):
        return {
            "tds": sample_meq.get("TDS", np.nan),
            "ratio_na_naca": np.nan,
            "ratio_cl_clhco3": np.nan,
            "gibbs_dominance": "unknown",
            "gibbs_evap_weight": np.nan,
            "gibbs_mix_weight": np.nan,
            "gibbs_evap_penalty": np.nan,
        }
    return compute_gibbs_metrics(sample_meq)


def _safe_exchange_direction(sample_meq: Mapping[str, float]) -> tuple[str, float, float]:
    if not {"Na", "Cl"}.issubset(sample_meq):
        return "unknown", np.nan, np.nan
    if not {"SO4", "HCO3", "NO3"} & set(sample_meq):
        return "unknown", np.nan, np.nan
    cai1, cai2 = compute_cai_indices(sample_meq)
    return (
        classify_exchange_direction(sample_meq),
        float(cai1) if cai1 is not None else np.nan,
        float(cai2) if cai2 is not None else np.nan,
    )


def _safe_redox_state(sample_mmol: Mapping[str, float]) -> str:
    if "NO3" not in sample_mmol and "Fe" not in sample_mmol:
        return "unknown"
    return classify_redox(sample_mmol)


def _ion_delta(
    upstream: Mapping[str, float],
    downstream: Mapping[str, float],
    ion: str,
) -> float:
    if ion not in upstream or ion not in downstream:
        return np.nan
    return float(downstream[ion]) - float(upstream[ion])


def _reaction_alignment(
    label: str,
    residual: np.ndarray,
    panel_ions: Sequence[str],
) -> float:
    vector = reaction_matrix(panel_ions)[REACTION_LABELS.index(label)]
    denominator = float(np.linalg.norm(vector) * np.linalg.norm(residual))
    if denominator <= 1e-12:
        return 0.0
    return float(np.dot(vector, residual) / denominator)


def _core_scale_from_score(score: float) -> float:
    return float(
        np.clip(
            math.exp(-2.0 * (score - 0.5)),
            CORE_EVIDENCE_MIN_SCALE,
            CORE_EVIDENCE_MAX_SCALE,
        )
    )


def hydrosheaf_core_evidence(
    upstream_vector: Sequence[float],
    observed_residual: Sequence[float],
    panel_ions: Sequence[str],
    upstream_si: Mapping[str, float] | None,
) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Build sparse-data Hydrosheaf evidence gates for the M5 reaction dictionary."""
    upstream_full = np.asarray(upstream_vector, dtype=float)
    residual = np.asarray(observed_residual, dtype=float)
    measured = list(panel_ions)
    downstream_full = upstream_full.copy()
    for offset, ion in enumerate(panel_ions):
        index = ION_ORDER.index(ion)
        downstream_full[index] = max(upstream_full[index] + residual[offset], 1e-10)

    upstream_mmol = _sample_from_vector(upstream_full, measured, equivalents=False)
    downstream_mmol = _sample_from_vector(downstream_full, measured, equivalents=False)
    upstream_meq = _sample_from_vector(upstream_full, measured, equivalents=True)
    downstream_meq = _sample_from_vector(downstream_full, measured, equivalents=True)

    gibbs_u = _safe_gibbs_metrics(upstream_meq)
    gibbs_d = _safe_gibbs_metrics(downstream_meq)
    exchange_u, cai1_u, cai2_u = _safe_exchange_direction(upstream_meq)
    exchange_d, cai1_d, cai2_d = _safe_exchange_direction(downstream_meq)
    redox_u = _safe_redox_state(upstream_mmol)
    redox_d = _safe_redox_state(downstream_mmol)
    si_values = upstream_si or {}

    tds_delta = float(downstream_meq.get("TDS", np.nan)) - float(
        upstream_meq.get("TDS", np.nan)
    )
    ca_delta = _ion_delta(upstream_mmol, downstream_mmol, "Ca")
    mg_delta = _ion_delta(upstream_mmol, downstream_mmol, "Mg")
    na_delta = _ion_delta(upstream_mmol, downstream_mmol, "Na")
    k_delta = _ion_delta(upstream_mmol, downstream_mmol, "K")
    hco3_delta = _ion_delta(upstream_mmol, downstream_mmol, "HCO3")
    cl_delta = _ion_delta(upstream_mmol, downstream_mmol, "Cl")
    so4_delta = _ion_delta(upstream_mmol, downstream_mmol, "SO4")
    no3_delta = _ion_delta(upstream_mmol, downstream_mmol, "NO3")
    f_delta = _ion_delta(upstream_mmol, downstream_mmol, "F")
    fe_delta = _ion_delta(upstream_mmol, downstream_mmol, "Fe")
    po4_delta = _ion_delta(upstream_mmol, downstream_mmol, "PO4")

    penalty_scales: list[float] = []
    rows: list[dict[str, object]] = []
    direct_exchange = {"CaNa_exch", "MgNa_exch"}
    reverse_exchange = {"NaCa_exch", "NaMg_exch"}

    for reaction in REACTIONS:
        label = reaction.label
        family = reaction.family
        score = 0.50
        flags: list[str] = []
        alignment = _reaction_alignment(label, residual, panel_ions)
        if alignment >= 0.55:
            score += 0.20
            flags.append("residual_strongly_aligned")
        elif alignment >= 0.25:
            score += 0.10
            flags.append("residual_weakly_aligned")
        elif alignment <= -0.25:
            score -= 0.18
            flags.append("residual_opposed")

        gibbs_downstream = str(gibbs_d["gibbs_dominance"])
        if gibbs_downstream == "rock" and family in {
            "carbonate",
            "silicate",
            "trace_mineral",
        }:
            score += 0.08
            flags.append("gibbs_rock_support")
        elif gibbs_downstream == "evaporation" and family == "evaporite":
            score += 0.14
            flags.append("gibbs_evaporation_support")
        elif gibbs_downstream == "precipitation" and family == "evaporite":
            score -= 0.10
            flags.append("gibbs_precipitation_evaporite_conflict")

        if label in direct_exchange | reverse_exchange:
            if exchange_d == "freshening":
                score += 0.22 if label in direct_exchange else -0.22
                flags.append(
                    "cai_freshening_support"
                    if label in direct_exchange
                    else "cai_freshening_conflict"
                )
            elif exchange_d == "salinization":
                score += 0.22 if label in reverse_exchange else -0.22
                flags.append(
                    "cai_salinization_support"
                    if label in reverse_exchange
                    else "cai_salinization_conflict"
                )

        if label == "NO3src":
            if math.isfinite(no3_delta) and no3_delta > 0:
                score += 0.18
                flags.append("nitrate_increase_support")
            if redox_d == "reducing":
                score -= 0.10
                flags.append("reducing_no3_source_conflict")
        elif label == "denit":
            if (
                math.isfinite(no3_delta)
                and math.isfinite(hco3_delta)
                and no3_delta < 0
                and hco3_delta > 0
            ):
                score += 0.20
                flags.append("denitrification_stoichiometry_support")
            elif math.isfinite(no3_delta) and no3_delta > 0:
                score -= 0.16
                flags.append("denitrification_no3_conflict")
        elif label == "pyrite_oxidation":
            if (
                math.isfinite(so4_delta)
                and math.isfinite(fe_delta)
                and so4_delta > 0
                and fe_delta >= -1e-8
            ):
                score += 0.14
                flags.append("pyrite_products_support")
            if redox_d == "reducing":
                score -= 0.16
                flags.append("reducing_pyrite_oxidation_conflict")

        if label == "halite" and math.isfinite(na_delta) and math.isfinite(cl_delta):
            if na_delta > 0 and cl_delta > 0:
                score += 0.16
                flags.append("halite_na_cl_support")
        elif label == "gypsum" and math.isfinite(ca_delta) and math.isfinite(so4_delta):
            if ca_delta > 0 and so4_delta > 0:
                score += 0.14
                flags.append("gypsum_ca_so4_support")
        elif label in {"calcite", "dolomite", "anorthite"}:
            carbonate_cation = ca_delta
            if label == "dolomite" and math.isfinite(mg_delta):
                carbonate_cation = min(ca_delta, mg_delta)
            if (
                math.isfinite(carbonate_cation)
                and math.isfinite(hco3_delta)
                and carbonate_cation > 0
                and hco3_delta > 0
            ):
                score += 0.10
                flags.append("carbonate_alkalinity_support")
        elif label == "albite" and math.isfinite(na_delta) and math.isfinite(hco3_delta):
            if na_delta > 0 and hco3_delta > 0:
                score += 0.10
                flags.append("albite_na_alkalinity_support")
        elif label == "k_feldspar" and math.isfinite(k_delta) and math.isfinite(hco3_delta):
            if k_delta > 0 and hco3_delta > 0:
                score += 0.10
                flags.append("k_feldspar_k_alkalinity_support")
        elif label == "fluorite" and math.isfinite(ca_delta) and math.isfinite(f_delta):
            if ca_delta > 0 and f_delta > 0:
                score += 0.10
                flags.append("fluorite_ca_f_support")
        elif label == "apatite" and math.isfinite(po4_delta) and math.isfinite(f_delta):
            if po4_delta > 0 and f_delta > 0:
                score += 0.10
                flags.append("apatite_po4_f_support")

        if reaction.si_phase:
            si_value = _finite_float(si_values.get(reaction.si_phase))
            if si_value is not None:
                if si_value < -SI_THRESHOLD and alignment > 0:
                    score += 0.10
                    flags.append("undersaturated_dissolution_support")
                elif si_value > SI_THRESHOLD and alignment > 0:
                    score -= 0.12
                    flags.append("saturated_dissolution_conflict")
                elif si_value > SI_THRESHOLD and alignment < 0:
                    score += 0.08
                    flags.append("saturated_precipitation_support")

        if family == "evaporite" and math.isfinite(tds_delta):
            if tds_delta > 30.0:
                score += 0.04
                flags.append("salinity_increase_support")
            elif tds_delta < -30.0:
                score -= 0.04
                flags.append("salinity_decrease_conflict")

        score = float(np.clip(score, 0.05, 0.95))
        penalty_scale = _core_scale_from_score(score)
        penalty_scales.append(penalty_scale)
        rows.append(
            {
                "reaction": label,
                "family": family,
                "panel_ions": ";".join(panel_ions),
                "hydrosheaf_core_evidence_score": score,
                "penalty_scale": penalty_scale,
                "residual_alignment": alignment,
                "gibbs_upstream": gibbs_u["gibbs_dominance"],
                "gibbs_downstream": gibbs_downstream,
                "exchange_upstream": exchange_u,
                "exchange_downstream": exchange_d,
                "redox_upstream": redox_u,
                "redox_downstream": redox_d,
                "cai1_upstream": cai1_u,
                "cai1_downstream": cai1_d,
                "cai2_upstream": cai2_u,
                "cai2_downstream": cai2_d,
                "tds_delta_mg_L": tds_delta,
                "evidence_flags": ";".join(flags),
            }
        )
    return np.asarray(penalty_scales, dtype=float), rows


def _sum_major_ions_mg(row: pd.Series) -> float:
    total = 0.0
    for ion in ION_ORDER:
        column = f"{ion}_mg_L"
        value = _finite_float(row.get(column))
        if value is not None:
            total += value
    return total


def _calibrate_field_ec_tds(hydro: pd.DataFrame) -> dict[str, tuple[float, float]]:
    ec_totals: list[float] = []
    tds_totals: list[float] = []
    ec_values: list[float] = []
    tds_values: list[float] = []
    for _, row in hydro.iterrows():
        total = _sum_major_ions_mg(row)
        if total <= 0:
            continue
        ec = _finite_float(row.get("EC_uS_cm"))
        tds = _finite_float(row.get("TDS_mg_L"))
        if ec is not None:
            ec_totals.append(total)
            ec_values.append(ec)
        if tds is not None:
            tds_totals.append(total)
            tds_values.append(tds)
    ec_model = fit_linear_model(ec_totals, ec_values) if ec_values else (0.0, 0.0)
    tds_model = fit_linear_model(tds_totals, tds_values) if tds_values else (0.0, 0.0)
    return {"ec": ec_model, "tds": tds_model}


def _field_optional_tracer_metrics(wet: pd.Series, dry: pd.Series) -> dict[str, object]:
    metrics: dict[str, object] = {}
    tracer_count = 0
    for column in ["SiO2_mg_L", "Sr_mg_L", "d18O_permil", "d2H_permil"]:
        wet_value = _finite_float(wet.get(column))
        dry_value = _finite_float(dry.get(column))
        short = column.replace("_mg_L", "").replace("_permil", "")
        metrics[f"wet_{short}"] = wet_value if wet_value is not None else np.nan
        metrics[f"dry_{short}"] = dry_value if dry_value is not None else np.nan
        if wet_value is not None and dry_value is not None:
            metrics[f"delta_{short}"] = dry_value - wet_value
            tracer_count += 1
        else:
            metrics[f"delta_{short}"] = np.nan

    wet_d18 = _finite_float(wet.get("d18O_permil"))
    wet_d2h = _finite_float(wet.get("d2H_permil"))
    dry_d18 = _finite_float(dry.get("d18O_permil"))
    dry_d2h = _finite_float(dry.get("d2H_permil"))
    if None not in (wet_d18, wet_d2h, dry_d18, dry_d2h):
        wet_excess = compute_d_excess(float(wet_d18), float(wet_d2h))
        dry_excess = compute_d_excess(float(dry_d18), float(dry_d2h))
        wet_evap = evaporation_index(float(wet_d18), float(wet_d2h), 10.0, 8.0)
        dry_evap = evaporation_index(float(dry_d18), float(dry_d2h), 10.0, 8.0)
        metrics.update(
            {
                "wet_d_excess": wet_excess,
                "dry_d_excess": dry_excess,
                "delta_d_excess": dry_excess - wet_excess,
                "wet_gmwl_deviation": wet_evap,
                "dry_gmwl_deviation": dry_evap,
                "delta_gmwl_deviation": dry_evap - wet_evap,
            }
        )
    else:
        metrics.update(
            {
                "wet_d_excess": np.nan,
                "dry_d_excess": np.nan,
                "delta_d_excess": np.nan,
                "wet_gmwl_deviation": np.nan,
                "dry_gmwl_deviation": np.nan,
                "delta_gmwl_deviation": np.nan,
            }
        )
    metrics["plus_lite_tracer_count"] = tracer_count
    return metrics


def _numeric_field_value(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        censored = text.startswith("<")
        text = text.lstrip("<>").strip()
        try:
            number = float(text)
        except ValueError:
            return None
        return 0.5 * number if censored else number
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _field_vector_from_row(
    row: pd.Series,
    column_map: Mapping[str, str],
) -> tuple[np.ndarray, list[str]]:
    values: list[float] = []
    measured: list[str] = []
    for ion in ION_ORDER:
        column = column_map.get(ion)
        value = _numeric_field_value(row.get(column)) if column else None
        if value is not None:
            values.append(value / MOLAR_MASS_G_MOL[ion])
            measured.append(ion)
        else:
            values.append(np.nan)
    return np.asarray(values, dtype=float), measured


def _field_major_ion_sum_mg(
    row: pd.Series,
    column_map: Mapping[str, str],
) -> float:
    total = 0.0
    for ion in ION_ORDER:
        column = column_map.get(ion)
        value = _numeric_field_value(row.get(column)) if column else None
        if value is not None:
            total += max(value, 0.0)
    return total


def _field_optional_diagnostics_from_rows(
    upstream: pd.Series,
    downstream: pd.Series,
    optional_columns: Mapping[str, str],
) -> tuple[list[str], np.ndarray]:
    diagnostics: list[str] = []
    observed: list[float] = []
    si_up = _numeric_field_value(upstream.get(optional_columns.get("SiO2", "")))
    si_down = _numeric_field_value(downstream.get(optional_columns.get("SiO2", "")))
    if si_up is not None and si_down is not None:
        diagnostics.append("SiO2")
        observed.append((si_down - si_up) / 10.0)

    sr_up = _numeric_field_value(upstream.get(optional_columns.get("Sr", "")))
    sr_down = _numeric_field_value(downstream.get(optional_columns.get("Sr", "")))
    if sr_up is not None and sr_down is not None:
        diagnostics.append("Sr")
        observed.append((sr_down - sr_up) / 0.1)

    d18_col = optional_columns.get("d18O")
    d2h_col = optional_columns.get("d2H")
    if d18_col and d2h_col:
        up_d18 = _numeric_field_value(upstream.get(d18_col))
        up_d2h = _numeric_field_value(upstream.get(d2h_col))
        down_d18 = _numeric_field_value(downstream.get(d18_col))
        down_d2h = _numeric_field_value(downstream.get(d2h_col))
        if None not in (up_d18, up_d2h, down_d18, down_d2h):
            upstream_evap = evaporation_index(
                float(up_d18), float(up_d2h), 10.0, 8.0
            )
            downstream_evap = evaporation_index(
                float(down_d18), float(down_d2h), 10.0, 8.0
            )
            diagnostics.append("water_isotope_evaporation")
            observed.append(downstream_evap - upstream_evap)
    return diagnostics, np.asarray(observed, dtype=float)


def _external_sample_records(
    frame: pd.DataFrame,
    *,
    dataset: str,
    id_column: str,
    x_column: str,
    y_column: str,
    z_column: str,
    column_map: Mapping[str, str],
    optional_columns: Mapping[str, str],
    metadata_columns: Sequence[str] = (),
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        sample_id = str(row.get(id_column, "")).strip()
        x = _numeric_field_value(row.get(x_column))
        y = _numeric_field_value(row.get(y_column))
        if not sample_id or x is None or y is None:
            continue
        vector, measured = _field_vector_from_row(row, column_map)
        if len(measured) < 6:
            continue
        records.append(
            {
                "dataset": dataset,
                "sample_id": sample_id,
                "x": float(x),
                "y": float(y),
                "z": _numeric_field_value(row.get(z_column)),
                "row": row,
                "vector": vector,
                "measured": measured,
                "major_sum_mg_L": _field_major_ion_sum_mg(row, column_map),
                "column_map": dict(column_map),
                "optional_columns": dict(optional_columns),
                "metadata": {
                    column: row.get(column, "")
                    for column in metadata_columns
                    if column in row.index
                },
            }
        )
    return records


def _nearest_external_edges(
    records: Sequence[Mapping[str, object]],
    *,
    neighbours: int = 2,
) -> list[tuple[Mapping[str, object], Mapping[str, object]]]:
    if len(records) < 2:
        return []
    from scipy.spatial import cKDTree

    coords = np.asarray([[float(item["x"]), float(item["y"])] for item in records])
    tree = cKDTree(coords)
    k = min(len(records), neighbours + 1)
    pairs: set[tuple[int, int]] = set()
    for index, item in enumerate(records):
        _, neighbours_idx = tree.query([float(item["x"]), float(item["y"])], k=k)
        for other in np.atleast_1d(neighbours_idx):
            other_index = int(other)
            if index == other_index:
                continue
            pairs.add(tuple(sorted((index, other_index))))

    oriented: list[tuple[Mapping[str, object], Mapping[str, object]]] = []
    for left_index, right_index in sorted(pairs):
        left = records[left_index]
        right = records[right_index]
        left_z = _finite_float(left.get("z"))
        right_z = _finite_float(right.get("z"))
        if left_z is not None and right_z is not None and left_z != right_z:
            upstream, downstream = (
                (left, right) if left_z > right_z else (right, left)
            )
        else:
            left_sum = float(left.get("major_sum_mg_L", 0.0))
            right_sum = float(right.get("major_sum_mg_L", 0.0))
            upstream, downstream = (
                (left, right) if left_sum <= right_sum else (right, left)
            )
        oriented.append((upstream, downstream))
    return oriented


def _delta_consistency_score(
    observed_delta: float | None,
    predicted_delta: float,
    scale_floor: float,
) -> float:
    if observed_delta is None or not math.isfinite(predicted_delta):
        return np.nan
    scale = scale_floor + abs(observed_delta)
    return float(math.exp(-abs(observed_delta - predicted_delta) / scale))


def _diagnostic_signature_matrix(diagnostics: Sequence[str]) -> np.ndarray:
    return np.asarray(
        [
            [
                float(DIAGNOSTIC_SIGNATURES.get(reaction.label, {}).get(name, 0.0))
                for name in diagnostics
            ]
            for reaction in REACTIONS
        ],
        dtype=float,
    )


def _simulate_optional_diagnostics(
    scenario: Mapping[str, object],
    tier: str,
    noise_level: float = DATA_TIER_NOISE_LEVEL,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    diagnostics = list(DATA_TIER_DIAGNOSTICS[tier])
    if not diagnostics:
        return diagnostics, np.asarray([], dtype=float), np.asarray([], dtype=float)
    signature = _diagnostic_signature_matrix(diagnostics)
    clean = signature.T @ np.asarray(scenario["true_extents"], dtype=float)
    seed = RANDOM_SEED + 500000 + int(scenario["scenario_index"]) * 113
    seed += int(noise_level * 1000) + 7919 * len(diagnostics)
    rng = np.random.default_rng(seed)
    sigma = noise_level * np.maximum(np.abs(clean), 0.02)
    observed = clean + rng.normal(0.0, sigma)
    return diagnostics, clean, observed


def _optional_diagnostic_evidence(
    diagnostics: Sequence[str],
    observed: Sequence[float],
) -> tuple[np.ndarray, list[dict[str, object]]]:
    if not diagnostics:
        rows = [
            {
                "reaction": reaction.label,
                "optional_evidence_score": 0.5,
                "optional_penalty_scale": 1.0,
                "optional_diagnostic_alignment": np.nan,
                "optional_support_alignment": np.nan,
                "diagnostics_with_nonzero_signature": 0,
            }
            for reaction in REACTIONS
        ]
        return np.ones(len(REACTIONS), dtype=float), rows

    observed_vector = np.asarray(observed, dtype=float)
    observed_norm = float(np.linalg.norm(observed_vector))
    signature = _diagnostic_signature_matrix(diagnostics)
    scales: list[float] = []
    rows: list[dict[str, object]] = []
    for reaction, vector in zip(REACTIONS, signature):
        signature_norm = float(np.linalg.norm(vector))
        if signature_norm <= 1e-12 or observed_norm <= 1e-12:
            raw_alignment = np.nan
            support_alignment = 0.0
        else:
            raw_alignment = float(
                np.dot(vector, observed_vector) / (signature_norm * observed_norm)
            )
            support_alignment = abs(raw_alignment) if reaction.signed else raw_alignment
        score = float(np.clip(0.5 + 0.35 * support_alignment, 0.05, 0.95))
        scale = float(np.clip(math.exp(-2.2 * (score - 0.5)), 0.25, 4.0))
        scales.append(scale)
        rows.append(
            {
                "reaction": reaction.label,
                "optional_evidence_score": score,
                "optional_penalty_scale": scale,
                "optional_diagnostic_alignment": raw_alignment,
                "optional_support_alignment": support_alignment,
                "diagnostics_with_nonzero_signature": int(np.count_nonzero(vector)),
            }
        )
    return np.asarray(scales, dtype=float), rows


def _combine_core_and_optional_penalties(
    core_scales: Sequence[float],
    optional_scales: Sequence[float],
) -> np.ndarray:
    combined = np.asarray(core_scales, dtype=float) * np.asarray(
        optional_scales, dtype=float
    )
    return np.clip(combined, 0.20, 5.00)


def _combine_evidence_scores(core_score: object, optional_score: object) -> float:
    core = _finite_float(core_score)
    optional = _finite_float(optional_score)
    core = 0.5 if core is None else core
    optional = 0.5 if optional is None else optional
    return float(np.clip(core + optional - 0.5, 0.05, 0.95))


def _ambiguous_class_members(class_map: Mapping[str, str]) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = {}
    for reaction, class_id in class_map.items():
        grouped.setdefault(class_id, []).append(reaction)
    return {
        class_id: sorted(members)
        for class_id, members in grouped.items()
        if len(members) > 1
    }


def _evidence_lifted_resolution_frame(
    evidence: pd.DataFrame,
    class_map: Mapping[str, str],
    *,
    score_column: str,
    group_columns: Sequence[str],
    evidence_source: str,
) -> pd.DataFrame:
    if evidence.empty or score_column not in evidence.columns:
        return pd.DataFrame()

    ambiguous = _ambiguous_class_members(class_map)
    rows: list[dict[str, object]] = []
    group_keys: str | list[str]
    group_keys = list(group_columns) if len(group_columns) > 1 else group_columns[0]
    for key_values, group in evidence.groupby(group_keys, dropna=False):
        if len(group_columns) == 1:
            key_tuple = (key_values,)
        elif isinstance(key_values, tuple):
            key_tuple = key_values
        else:
            key_tuple = tuple(key_values)
        group_context = dict(zip(group_columns, key_tuple))
        by_reaction = group.drop_duplicates("reaction").set_index("reaction")
        for class_id, members in ambiguous.items():
            available = [member for member in members if member in by_reaction.index]
            if len(available) != len(members):
                continue
            scores = {
                member: float(by_reaction.loc[member, score_column])
                for member in members
            }
            resolution = evidence_lifted_resolution(
                members,
                scores,
                class_id=class_id,
            )
            row = {
                **group_context,
                "evidence_source": evidence_source,
                **resolution.as_row(),
            }
            if "true_active" in by_reaction.columns:
                true_members = [
                    member
                    for member in members
                    if bool(by_reaction.loc[member, "true_active"])
                ]
                row["true_active_members"] = ";".join(true_members)
                row["n_true_active_members"] = len(true_members)
                row["top_member_true_active"] = resolution.top_member in true_members
            if "recovered_active" in by_reaction.columns:
                recovered_members = [
                    member
                    for member in members
                    if bool(by_reaction.loc[member, "recovered_active"])
                ]
                row["recovered_active_members"] = ";".join(recovered_members)
                row["n_recovered_active_members"] = len(recovered_members)
                row["top_member_recovered_active"] = (
                    resolution.top_member in recovered_members
                )
            rows.append(row)
    return pd.DataFrame(rows)


def write_core_evidence_lifted_resolution(
    class_map: Mapping[str, str],
    evidence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if evidence is None:
        evidence_path = RESULTS_DIR / "hydrosheaf_core_evidence.csv"
        if not evidence_path.exists():
            return pd.DataFrame()
        evidence = pd.read_csv(evidence_path)
    resolution = _evidence_lifted_resolution_frame(
        evidence,
        class_map,
        score_column="hydrosheaf_core_evidence_score",
        group_columns=[
            "scenario_id",
            "archetype",
            "replicate",
            "noise_level",
            "transport_error_level",
            "panel",
            "method",
        ],
        evidence_source="hydrosheaf_core_major_ions",
    )
    resolution.to_csv(
        RESULTS_DIR / "hydrosheaf_core_evidence_lifted_resolution.csv",
        index=False,
    )
    return resolution


def write_field_evidence_lifted_resolution(
    class_map: Mapping[str, str],
    evidence: pd.DataFrame,
) -> pd.DataFrame:
    resolution = _evidence_lifted_resolution_frame(
        evidence,
        class_map,
        score_column="hydrosheaf_core_evidence_score",
        group_columns=[
            "well_id",
            "aquifer",
            "geology_group",
            "lithology",
            "region",
            "district",
        ],
        evidence_source="ghana_field_hydrosheaf_core",
    )
    resolution.to_csv(
        RESULTS_DIR / "ghana_evidence_lifted_resolution.csv",
        index=False,
    )
    return resolution


def _draw_support(
    archetype: str,
    count: int,
    rng: np.random.Generator,
) -> list[str]:
    pool = list(ARCHETYPE_POOLS[archetype])
    rng.shuffle(pool)
    selected: list[str] = []
    opposite_pairs = [
        {"CaNa_exch", "NaCa_exch"},
        {"MgNa_exch", "NaMg_exch"},
    ]
    for label in pool:
        if any(label in pair and any(item in pair for item in selected) for pair in opposite_pairs):
            continue
        selected.append(label)
        if len(selected) == count:
            return selected
    return selected


def _draw_extent(label: str, rng: np.random.Generator) -> float:
    low, high = EXTENT_RANGES.get(label, (0.025, 0.280))
    extent = float(rng.uniform(low, high))
    if REACTION_BY_LABEL[label].signed and rng.random() < 0.22:
        extent *= -1.0
    return extent


def generate_scenarios() -> list[dict[str, object]]:
    rng = np.random.default_rng(RANDOM_SEED)
    matrix = reaction_matrix()
    scenarios: list[dict[str, object]] = []
    scenario_index = 0
    for archetype, base in ARCHETYPES.items():
        base_vector = np.asarray([base[ion] for ion in ION_ORDER], dtype=float)
        for replicate in range(N_SCENARIOS_PER_ARCHETYPE):
            n_active = 1 + (replicate % 5)
            support = _draw_support(archetype, n_active, rng)
            extents = np.zeros(len(REACTIONS), dtype=float)
            for label in support:
                extents[REACTION_LABELS.index(label)] = _draw_extent(label, rng)
            delta = matrix.T @ extents
            scale = 1.0
            while np.any(base_vector + scale * delta < 1e-8):
                scale *= 0.7
                if scale < 0.05:
                    raise RuntimeError(f"Could not construct positive scenario {scenario_index}.")
            extents *= scale
            transport_error_level = [0.0, 0.02, 0.05][replicate % 3]
            scenarios.append(
                {
                    "scenario_id": f"M5S{scenario_index + 1:04d}",
                    "scenario_index": scenario_index,
                    "archetype": archetype,
                    "replicate": replicate + 1,
                    "n_true_reactions": len(support),
                    "true_support": ";".join(sorted(support)),
                    "transport_error_level": transport_error_level,
                    "topology_confidence": float(rng.uniform(0.55, 0.98)),
                    "residence_time_confidence": float(rng.uniform(0.50, 0.96)),
                    "base": dict(base),
                    "true_extents": extents,
                }
            )
            scenario_index += 1
    return scenarios


def _solution_lines(solution_id: int, base: Mapping[str, float]) -> list[str]:
    return [
        f"SOLUTION {solution_id} upstream",
        "    units mmol/kgw",
        f"    temp {base['temp']:.6f}",
        f"    pH {base['pH']:.6f}",
        "    pe 4",
        f"    Ca {base['Ca']:.12g}",
        f"    Mg {base['Mg']:.12g}",
        f"    Na {base['Na']:.12g}",
        f"    K {base['K']:.12g}",
        f"    Cl {base['Cl']:.12g}",
        f"    S(6) {base['SO4']:.12g} as SO4",
        f"    N(5) {base['NO3']:.12g} as NO3",
        f"    F {base['F']:.12g}",
        f"    Fe {base['Fe']:.12g}",
        f"    P {base['PO4']:.12g} as PO4",
        f"    Alkalinity {base['HCO3']:.12g} as HCO3",
        "    -water 1.0",
    ]


def _reaction_components(extents: Sequence[float]) -> dict[str, float]:
    components: dict[str, float] = {}
    for reaction, extent in zip(REACTIONS, extents):
        for formula, coefficient in reaction.phreeqc_components:
            components[formula] = components.get(formula, 0.0) + (
                float(extent) * coefficient / 1000.0
            )
    return {
        formula: value
        for formula, value in components.items()
        if abs(value) > 1e-14
    }


def build_phreeqc_input(
    scenarios: Sequence[Mapping[str, object]],
    selected_name: str,
) -> str:
    phases = sorted(
        {
            reaction.si_phase
            for reaction in REACTIONS
            if reaction.si_phase is not None
        }
    )
    totals = " ".join(PHREEQC_TOTALS.values())
    lines: list[str] = []
    for index, scenario in enumerate(scenarios, start=1):
        lines.extend(_solution_lines(index, scenario["base"]))
        if index == 1:
            lines.extend(
                [
                    "SELECTED_OUTPUT 1",
                    f"    -file {selected_name}",
                    "    -reset false",
                    "    -high_precision true",
                    "    -simulation true",
                    "    -solution true",
                    "    -pH true",
                    "    -charge_balance true",
                    f"    -totals {totals}",
                    f"    -saturation_indices {' '.join(phases)}",
                ]
            )
        lines.append("END")
        lines.append("")
        lines.append(f"USE solution {index}")
        lines.append(f"REACTION {index}")
        for formula, coefficient in sorted(
            _reaction_components(scenario["true_extents"]).items()
        ):
            lines.append(f"    {formula} {coefficient:.14g}")
        lines.append("    1 moles in 1 steps")
        lines.append(f"SAVE solution {10000 + index}")
        lines.append("END")
        lines.append("")
    return "\n".join(lines)


def run_live_phreeqc(
    scenarios: list[dict[str, object]],
) -> tuple[pd.DataFrame, float, Path, Path]:
    executable, database = find_phreeqc()
    input_path = PHREEQC_DIR / "m5_factorial_benchmark.pqi"
    output_path = PHREEQC_DIR / "m5_factorial_benchmark.out"
    selected_path = PHREEQC_DIR / "m5_factorial_selected.tsv"
    input_path.write_text(
        build_phreeqc_input(scenarios, selected_path.name),
        encoding="ascii",
    )
    if selected_path.exists():
        selected_path.unlink()
    elapsed = run_phreeqc(input_path, output_path, database, executable)
    raw = pd.read_csv(selected_path, sep="\t")
    raw.columns = [str(column).strip() for column in raw.columns]
    if len(raw) != 2 * len(scenarios):
        raise RuntimeError(
            f"Expected {2 * len(scenarios)} PHREEQC rows, found {len(raw)}."
        )

    matrix = reaction_matrix()
    rows: list[dict[str, object]] = []
    for index, scenario in enumerate(scenarios):
        upstream_row = raw.iloc[2 * index]
        downstream_row = raw.iloc[2 * index + 1]
        upstream = np.asarray(
            [float(upstream_row[PHREEQC_TOTALS[ion]]) * 1000.0 for ion in ION_ORDER]
        )
        downstream = np.asarray(
            [float(downstream_row[PHREEQC_TOTALS[ion]]) * 1000.0 for ion in ION_ORDER]
        )
        true_delta = downstream - upstream
        expected_delta = matrix.T @ np.asarray(scenario["true_extents"], dtype=float)
        generation_rmse = float(np.sqrt(np.mean((true_delta - expected_delta) ** 2)))
        upstream_si: dict[str, float] = {}
        downstream_si: dict[str, float] = {}
        for reaction in REACTIONS:
            if not reaction.si_phase:
                continue
            heading = f"si_{reaction.si_phase}"
            upstream_si[reaction.si_phase] = float(upstream_row.get(heading, np.nan))
            downstream_si[reaction.si_phase] = float(downstream_row.get(heading, np.nan))
        scenario["upstream"] = upstream
        scenario["downstream"] = downstream
        scenario["true_delta"] = true_delta
        scenario["upstream_si"] = upstream_si
        scenario["downstream_si"] = downstream_si
        scenario["generation_rmse"] = generation_rmse
        row: dict[str, object] = {
            "scenario_id": scenario["scenario_id"],
            "archetype": scenario["archetype"],
            "replicate": scenario["replicate"],
            "n_true_reactions": scenario["n_true_reactions"],
            "true_support": scenario["true_support"],
            "transport_error_level": scenario["transport_error_level"],
            "topology_confidence": scenario["topology_confidence"],
            "residence_time_confidence": scenario["residence_time_confidence"],
            "generation_rmse_mmolL": generation_rmse,
            "upstream_pH": float(upstream_row["pH"]),
            "downstream_pH": float(downstream_row["pH"]),
            "upstream_charge_balance_raw": float(upstream_row["charge"]),
            "downstream_charge_balance_raw": float(downstream_row["charge"]),
        }
        for ion, up_value, down_value, delta_value in zip(
            ION_ORDER, upstream, downstream, true_delta
        ):
            row[f"upstream_{ion}"] = up_value
            row[f"downstream_{ion}"] = down_value
            row[f"delta_{ion}"] = delta_value
        for reaction, extent in zip(REACTIONS, scenario["true_extents"]):
            row[f"true_extent_{reaction.label}"] = float(extent)
        for phase, value in upstream_si.items():
            row[f"upstream_si_{phase}"] = value
        for phase, value in downstream_si.items():
            row[f"downstream_si_{phase}"] = value
        rows.append(row)
    return pd.DataFrame(rows), elapsed, executable, database


def hydrate_scenarios_from_phreeqc(
    scenarios: Sequence[dict[str, object]],
    phreeqc_frame: pd.DataFrame,
) -> None:
    """Attach PHREEQC chemistry vectors to generated scenario definitions."""
    by_id = phreeqc_frame.set_index("scenario_id")
    for scenario in scenarios:
        row = by_id.loc[scenario["scenario_id"]]
        scenario["upstream"] = np.asarray(
            [row[f"upstream_{ion}"] for ion in ION_ORDER], dtype=float
        )
        scenario["downstream"] = np.asarray(
            [row[f"downstream_{ion}"] for ion in ION_ORDER], dtype=float
        )
        scenario["true_delta"] = np.asarray(
            [row[f"delta_{ion}"] for ion in ION_ORDER], dtype=float
        )
        scenario["upstream_si"] = {
            reaction.si_phase: float(
                row.get(f"upstream_si_{reaction.si_phase}", np.nan)
            )
            for reaction in REACTIONS
            if reaction.si_phase
        }


def _phreeqc_inverse_phase_definition(labels: Sequence[str]) -> list[str]:
    lines = ["PHASES"]
    for label in labels:
        _, phase_lines = PHREEQC_INVERSE_PHASES[label]
        lines.extend(phase_lines)
    lines.append("END")
    return lines


def _inverse_solution_lines(
    solution_id: int,
    name: str,
    scenario: Mapping[str, object],
    row: pd.Series,
    prefix: str,
) -> list[str]:
    base = scenario["base"]
    return [
        f"SOLUTION {solution_id} {name}",
        "    units mmol/kgw",
        f"    temp {base['temp']:.6f}",
        f"    pH {row[f'{prefix}_pH']:.8g}",
        "    pe 4",
        f"    Ca {row[f'{prefix}_Ca']:.12g}",
        f"    Mg {row[f'{prefix}_Mg']:.12g}",
        f"    Na {row[f'{prefix}_Na']:.12g}",
        f"    K {row[f'{prefix}_K']:.12g}",
        f"    Cl {row[f'{prefix}_Cl']:.12g}",
        f"    S(6) {row[f'{prefix}_SO4']:.12g} as SO4",
        f"    N(5) {row[f'{prefix}_NO3']:.12g} as NO3",
        f"    F {row[f'{prefix}_F']:.12g}",
        f"    Fe {row[f'{prefix}_Fe']:.12g}",
        f"    P {row[f'{prefix}_PO4']:.12g} as PO4",
        f"    Alkalinity {row[f'{prefix}_HCO3']:.12g} as HCO3",
        "    -water 1.0",
    ]


def _projection_score(label: str, delta: np.ndarray) -> float:
    vector = reaction_matrix()[REACTION_LABELS.index(label)]
    denom = float(vector @ vector)
    if denom <= 0.0:
        return 0.0
    return abs(float(delta @ vector) / denom)


def _phreeqc_inverse_candidates(
    scenario: Mapping[str, object],
) -> tuple[list[str], list[str]]:
    pool = list(ARCHETYPE_POOLS[str(scenario["archetype"])])
    m5_labels = [
        label for label in pool if label in PHREEQC_INVERSE_PHASES
    ]
    for label in ("halite", "NO3src", "denit"):
        if label not in m5_labels:
            m5_labels.append(label)
    m5_labels = list(dict.fromkeys(m5_labels))
    if len(m5_labels) > PHREEQC_INVERSE_MAX_M5_PHASES:
        delta = np.asarray(scenario["true_delta"], dtype=float)
        priority = [
            label for label in ("halite", "NO3src", "denit") if label in m5_labels
        ]
        scored = sorted(
            (label for label in m5_labels if label not in priority),
            key=lambda label: _projection_score(label, delta),
            reverse=True,
        )
        m5_labels = (priority + scored)[:PHREEQC_INVERSE_MAX_M5_PHASES]

    exchange_phases: list[str] = []
    if any(label in pool for label in ("CaNa_exch", "NaCa_exch")):
        exchange_phases.extend(["CaX2", "NaX"])
    if any(label in pool for label in ("MgNa_exch", "NaMg_exch")):
        exchange_phases.extend(["MgX2", "NaX"])
    return m5_labels, list(dict.fromkeys(exchange_phases))


def _phreeqc_inverse_fallback_candidates(
    m5_labels: Sequence[str],
    exchange_phases: Sequence[str],
) -> tuple[list[str], list[str]]:
    balance_labels = ["calcite", "dolomite", "gypsum", "halite", "NO3src", "denit"]
    labels = list(dict.fromkeys([*m5_labels, *balance_labels]))
    exchanges = list(dict.fromkeys([*exchange_phases, "CaX2", "NaX", "MgX2"]))
    return labels, exchanges


def _build_phreeqc_inverse_input(
    scenario: Mapping[str, object],
    row: pd.Series,
    m5_labels: Sequence[str],
    exchange_phases: Sequence[str],
    uncertainty: float = PHREEQC_INVERSE_UNCERTAINTY,
) -> str:
    phase_names = [
        PHREEQC_INVERSE_PHASES[label][0] for label in m5_labels
    ] + list(exchange_phases)
    lines = _phreeqc_inverse_phase_definition(m5_labels)
    lines.append("")
    lines.extend(_inverse_solution_lines(1, "upstream", scenario, row, "upstream"))
    lines.extend(["END", ""])
    lines.extend(_inverse_solution_lines(2, "downstream", scenario, row, "downstream"))
    lines.extend(["END", ""])
    lines.extend(
        [
            "INVERSE_MODELING 1",
            "    -solutions 1 2",
            f"    -uncertainty {uncertainty}",
            "    -range",
            "    -phases",
        ]
    )
    for phase in phase_names:
        lines.append(f"        {phase}")
    lines.extend(["END", ""])
    return "\n".join(lines)


def _parse_phreeqc_inverse_output(
    output_text: str,
    phase_names: Sequence[str],
) -> tuple[list[dict[str, object]], dict[str, int]]:
    phase_name_set = set(phase_names)
    models: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    in_phase_table = False
    transfer_pattern = re.compile(
        r"^\s*(\S+)\s+([-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)"
    )
    for line in output_text.splitlines():
        if line.startswith("Phase mole transfers:"):
            if current is not None:
                models.append(current)
            current = {"phase_transfers_mol": {}, "minimal": False}
            in_phase_table = True
            continue
        if current is not None and "Model contains minimum number of phases" in line:
            current["minimal"] = True
        if not in_phase_table or current is None:
            continue
        if (
            not line.strip()
            or line.lstrip().startswith("Redox")
            or line.lstrip().startswith("Sum of")
            or line.startswith("=")
        ):
            in_phase_table = False
            continue
        match = transfer_pattern.match(line)
        if not match:
            continue
        phase, transfer = match.group(1), float(match.group(2))
        if phase in phase_name_set:
            transfers = current["phase_transfers_mol"]
            assert isinstance(transfers, dict)
            transfers[phase] = transfer
    if current is not None:
        models.append(current)

    counts = {"models_found": 0, "minimal_models_found": 0}
    model_match = re.search(r"Number of models found:\s+(\d+)", output_text)
    minimal_match = re.search(
        r"Number of minimal models found:\s+(\d+)", output_text
    )
    if model_match:
        counts["models_found"] = int(model_match.group(1))
    if minimal_match:
        counts["minimal_models_found"] = int(minimal_match.group(1))
    return models, counts


def _exchange_support_from_transfers(
    transfers: Mapping[str, float],
    scenario: Mapping[str, object],
) -> set[str]:
    threshold_mol = SUPPORT_THRESHOLD / 1000.0
    delta = np.asarray(scenario["true_delta"], dtype=float)
    support: set[str] = set()
    if abs(transfers.get("CaX2", 0.0)) >= threshold_mol or abs(
        transfers.get("NaX", 0.0)
    ) >= threshold_mol:
        vector = reaction_matrix()[REACTION_LABELS.index("CaNa_exch")]
        coefficient = float(delta @ vector) / float(vector @ vector)
        if abs(coefficient) >= SUPPORT_THRESHOLD:
            support.add("CaNa_exch" if coefficient > 0 else "NaCa_exch")
    if abs(transfers.get("MgX2", 0.0)) >= threshold_mol or abs(
        transfers.get("NaX", 0.0)
    ) >= threshold_mol:
        vector = reaction_matrix()[REACTION_LABELS.index("MgNa_exch")]
        coefficient = float(delta @ vector) / float(vector @ vector)
        if abs(coefficient) >= SUPPORT_THRESHOLD:
            support.add("MgNa_exch" if coefficient > 0 else "NaMg_exch")
    return support


def _phreeqc_inverse_support(
    model: Mapping[str, object],
    scenario: Mapping[str, object],
) -> set[str]:
    threshold_mol = SUPPORT_THRESHOLD / 1000.0
    transfers = model.get("phase_transfers_mol", {})
    if not isinstance(transfers, Mapping):
        return set()
    support = {
        PHREEQC_INVERSE_PHASE_TO_REACTION[phase]
        for phase, transfer in transfers.items()
        if phase in PHREEQC_INVERSE_PHASE_TO_REACTION
        and abs(float(transfer)) >= threshold_mol
    }
    support.update(_exchange_support_from_transfers(transfers, scenario))
    return support


def _support_metric_row(
    scenario: Mapping[str, object],
    support: set[str],
    class_map: Mapping[str, str],
) -> dict[str, float | str]:
    truth_support = support_from_extents(
        np.asarray(scenario["true_extents"], dtype=float), SUPPORT_THRESHOLD
    )
    precision, recall, f1, false_discovery = precision_recall_f1(
        truth_support, support
    )
    class_precision, class_recall, class_f1, class_false_discovery = (
        precision_recall_f1(
            class_set(truth_support, class_map),
            class_set(support, class_map),
        )
    )
    return {
        "selected_support": ";".join(sorted(support)),
        "selected_classes": ";".join(sorted(class_set(support, class_map))),
        "phase_precision": precision,
        "phase_recall": recall,
        "phase_f1": f1,
        "phase_false_discovery_rate": false_discovery,
        "class_precision": class_precision,
        "class_recall": class_recall,
        "class_f1": class_f1,
        "class_false_discovery_rate": class_false_discovery,
    }


def _prefixed_metrics(
    prefix: str,
    metrics: Mapping[str, float | str],
) -> dict[str, float | str]:
    return {f"{prefix}_{key}": value for key, value in metrics.items()}


def run_phreeqc_inverse_baseline(
    scenarios: Sequence[Mapping[str, object]],
    phreeqc_frame: pd.DataFrame,
    class_map: Mapping[str, str],
    executable: Path,
    database: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    """Run a conventional PHREEQC INVERSE_MODELING baseline on clean pairs."""
    by_id = phreeqc_frame.set_index("scenario_id")
    scenario_rows: list[dict[str, object]] = []
    model_rows: list[dict[str, object]] = []
    total_elapsed = 0.0
    for stale_path in PHREEQC_INVERSE_DIR.glob("M5S*.*"):
        if stale_path.suffix.lower() in {".pqi", ".out", ".log"}:
            stale_path.unlink()
    for index, scenario in enumerate(scenarios, start=1):
        scenario_id = str(scenario["scenario_id"])
        row = by_id.loc[scenario_id]
        m5_labels, exchange_phases = _phreeqc_inverse_candidates(scenario)
        active_labels = list(m5_labels)
        active_exchange_phases = list(exchange_phases)
        active_uncertainty = PHREEQC_INVERSE_UNCERTAINTY
        active_setup = "strict_5pct"
        error_message = ""
        strict_error_message = ""
        models: list[dict[str, object]] = []
        counts = {"models_found": 0, "minimal_models_found": 0}
        attempts = [
            (
                "strict_5pct",
                m5_labels,
                exchange_phases,
                PHREEQC_INVERSE_UNCERTAINTY,
            )
        ]
        fallback_labels, fallback_exchange_phases = (
            _phreeqc_inverse_fallback_candidates(m5_labels, exchange_phases)
        )
        attempts.append(
            (
                "relaxed_20pct",
                fallback_labels,
                fallback_exchange_phases,
                PHREEQC_INVERSE_FALLBACK_UNCERTAINTY,
            )
        )
        for setup_name, labels, exchanges, uncertainty in attempts:
            active_labels = list(labels)
            active_exchange_phases = list(exchanges)
            active_uncertainty = uncertainty
            active_setup = setup_name
            phase_names = [
                PHREEQC_INVERSE_PHASES[label][0] for label in labels
            ] + list(exchanges)
            input_path = PHREEQC_INVERSE_DIR / f"{scenario_id}_{setup_name}.pqi"
            output_path = PHREEQC_INVERSE_DIR / f"{scenario_id}_{setup_name}.out"
            input_path.write_text(
                _build_phreeqc_inverse_input(
                    scenario, row, labels, exchanges, uncertainty
                ),
                encoding="ascii",
            )
            try:
                elapsed = run_phreeqc(input_path, output_path, database, executable)
                total_elapsed += elapsed
                output_text = output_path.read_text(
                    encoding="utf-8", errors="replace"
                )
                models, counts = _parse_phreeqc_inverse_output(
                    output_text, phase_names
                )
                error_message = ""
                break
            except RuntimeError as exc:
                message = str(exc).splitlines()[0]
                if setup_name == "strict_5pct":
                    strict_error_message = message
                error_message = message

        minimal_models = [
            model for model in models if bool(model.get("minimal"))
        ] or models
        minimal_model_ids = {id(model) for model in minimal_models}
        first_support: set[str] = set()
        first_minimal_seen = False
        minimal_union_support: set[str] = set()
        oracle_support: set[str] = set()
        oracle_phase_f1 = -1.0
        oracle_class_f1 = -1.0
        for model_index, model in enumerate(models, start=1):
            support = _phreeqc_inverse_support(model, scenario)
            metrics = _support_metric_row(scenario, support, class_map)
            if id(model) in minimal_model_ids:
                minimal_union_support.update(support)
                if not first_minimal_seen:
                    first_support = set(support)
                    first_minimal_seen = True
                if (
                    float(metrics["class_f1"]) > oracle_class_f1
                    or (
                        float(metrics["class_f1"]) == oracle_class_f1
                        and float(metrics["phase_f1"]) > oracle_phase_f1
                    )
                ):
                    oracle_class_f1 = float(metrics["class_f1"])
                    oracle_phase_f1 = float(metrics["phase_f1"])
                    oracle_support = set(support)
            transfers = model.get("phase_transfers_mol", {})
            transfer_summary = ""
            if isinstance(transfers, Mapping):
                transfer_summary = ";".join(
                    f"{phase}:{1000.0 * float(value):.6g}"
                    for phase, value in sorted(transfers.items())
                    if abs(float(value)) >= SUPPORT_THRESHOLD / 1000.0
                )
            model_rows.append(
                {
                    "scenario_id": scenario_id,
                    "archetype": scenario["archetype"],
                    "model_index": model_index,
                    "minimal": bool(model.get("minimal")),
                    "phase_transfers_mmolL": transfer_summary,
                    **metrics,
                }
            )

        base_row: dict[str, object] = {
            "scenario_id": scenario_id,
            "archetype": scenario["archetype"],
            "replicate": scenario["replicate"],
            "n_true_reactions": scenario["n_true_reactions"],
            "true_support": scenario["true_support"],
            "phreeqc_inverse_setup": active_setup,
            "phreeqc_inverse_uncertainty": active_uncertainty,
            "strict_error_message": strict_error_message,
            "candidate_reactions": ";".join(active_labels),
            "candidate_exchange_phases": ";".join(active_exchange_phases),
            "candidate_phase_count": len(active_labels) + len(active_exchange_phases),
            "phreeqc_inverse_success": counts["models_found"] > 0,
            "models_found": counts["models_found"],
            "minimal_models_found": counts["minimal_models_found"],
            "parsed_models": len(models),
            "parsed_minimal_models": sum(
                1 for model in models if bool(model.get("minimal"))
            ),
            "error_message": error_message,
        }
        base_row.update(
            _prefixed_metrics(
                "first_minimal",
                _support_metric_row(scenario, first_support, class_map),
            )
        )
        base_row.update(
            _prefixed_metrics(
                "minimal_union",
                _support_metric_row(scenario, minimal_union_support, class_map),
            )
        )
        base_row.update(
            _prefixed_metrics(
                "oracle_minimal",
                _support_metric_row(scenario, oracle_support, class_map),
            )
        )
        scenario_rows.append(base_row)
        if index % 40 == 0:
            print(f"M5: PHREEQC inverse baseline parsed {index} scenarios...")
    return pd.DataFrame(scenario_rows), pd.DataFrame(model_rows), total_elapsed


def _load_or_run_phreeqc_inverse_baseline(
    scenarios: Sequence[Mapping[str, object]],
    phreeqc_frame: pd.DataFrame,
    class_map: Mapping[str, str],
    executable: Path,
    database: Path,
    *,
    reuse_existing: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    baseline_path = RESULTS_DIR / "phreeqc_inverse_baseline.csv"
    models_path = RESULTS_DIR / "phreeqc_inverse_baseline_models.csv"
    if reuse_existing and baseline_path.exists() and models_path.exists():
        print("M5: reusing conventional PHREEQC inverse-model baseline...")
        return pd.read_csv(baseline_path), pd.read_csv(models_path), float("nan")
    print("M5: running conventional PHREEQC inverse-model baseline...")
    phreeqc_inverse, phreeqc_inverse_models, elapsed = run_phreeqc_inverse_baseline(
        scenarios, phreeqc_frame, class_map, executable, database
    )
    phreeqc_inverse.to_csv(baseline_path, index=False)
    phreeqc_inverse_models.to_csv(models_path, index=False)
    return phreeqc_inverse, phreeqc_inverse_models, elapsed


def observed_residual(
    scenario: Mapping[str, object],
    noise_level: float,
) -> np.ndarray:
    seed = RANDOM_SEED + int(scenario["scenario_index"]) * 101 + int(noise_level * 1000)
    rng = np.random.default_rng(seed)
    upstream = np.asarray(scenario["upstream"], dtype=float)
    downstream = np.asarray(scenario["downstream"], dtype=float)
    up_sigma = noise_level * np.maximum(upstream, 0.02)
    down_sigma = noise_level * np.maximum(downstream, 0.02)
    observed_up = upstream + rng.normal(0.0, up_sigma)
    observed_down = downstream + rng.normal(0.0, down_sigma)
    confounder_level = float(scenario["transport_error_level"])
    confounder = rng.normal(
        0.0,
        confounder_level * np.maximum(upstream, 0.05),
    )
    return observed_down - observed_up + confounder


def evaluate_fit(
    scenario: Mapping[str, object],
    extents: np.ndarray,
    observed: np.ndarray,
    panel_ions: Sequence[str],
    class_map: Mapping[str, str],
    support_threshold: float = SUPPORT_THRESHOLD,
) -> dict[str, float | str]:
    truth_extents = np.asarray(scenario["true_extents"], dtype=float)
    truth_support = support_from_extents(truth_extents, SUPPORT_THRESHOLD)
    predicted_support = support_from_extents(extents, support_threshold)
    precision, recall, f1, false_discovery = precision_recall_f1(
        truth_support, predicted_support
    )
    truth_classes = class_set(truth_support, class_map)
    predicted_classes = class_set(predicted_support, class_map)
    class_precision, class_recall, class_f1, class_false_discovery = precision_recall_f1(
        truth_classes, predicted_classes
    )
    active_intersection = truth_support & predicted_support
    direction_accuracy = (
        np.mean(
            [
                np.sign(extents[REACTION_LABELS.index(label)])
                == np.sign(truth_extents[REACTION_LABELS.index(label)])
                for label in active_intersection
            ]
        )
        if active_intersection
        else np.nan
    )
    full_prediction = reaction_matrix().T @ extents
    panel_indices = [ION_ORDER.index(ion) for ion in panel_ions]
    hidden_indices = [index for index in range(len(ION_ORDER)) if index not in panel_indices]
    reconstruction_rmse = float(
        np.sqrt(np.mean((observed[panel_indices] - full_prediction[panel_indices]) ** 2))
    )
    heldout_rmse = (
        float(
            np.sqrt(
                np.mean(
                    (
                        np.asarray(scenario["true_delta"])[hidden_indices]
                        - full_prediction[hidden_indices]
                    )
                    ** 2
                )
            )
        )
        if hidden_indices
        else np.nan
    )
    active_truth = np.abs(truth_extents) >= SUPPORT_THRESHOLD
    relative_bias = (
        float(
            np.mean(
                (extents[active_truth] - truth_extents[active_truth])
                / np.maximum(np.abs(truth_extents[active_truth]), 1e-12)
            )
        )
        if np.any(active_truth)
        else np.nan
    )
    return {
        "phase_precision": precision,
        "phase_recall": recall,
        "phase_f1": f1,
        "false_discovery_rate": false_discovery,
        "class_precision": class_precision,
        "class_recall": class_recall,
        "class_f1": class_f1,
        "class_false_discovery_rate": class_false_discovery,
        "direction_accuracy": direction_accuracy,
        "extent_rmse_mmolL": float(np.sqrt(np.mean((extents - truth_extents) ** 2))),
        "extent_mae_mmolL": float(np.mean(np.abs(extents - truth_extents))),
        "relative_extent_bias": relative_bias,
        "reconstruction_rmse_mmolL": reconstruction_rmse,
        "heldout_rmse_mmolL": heldout_rmse,
        "truth_support": ";".join(sorted(truth_support)),
        "selected_support": ";".join(sorted(predicted_support)),
        "truth_classes": ";".join(sorted(truth_classes)),
        "selected_classes": ";".join(sorted(predicted_classes)),
    }


def calibrate_hyperparameters(
    scenarios: Sequence[Mapping[str, object]],
    class_map: Mapping[str, str],
) -> tuple[dict[str, float], pd.DataFrame]:
    training = [
        scenario
        for scenario in scenarios
        if scenario["archetype"] != "mixed" and int(scenario["replicate"]) <= 12
    ]
    lambda_grid = np.geomspace(0.0003, 0.08, 9)
    l2_grid = [0.005, 0.02, 0.08]
    rows: list[dict[str, object]] = []
    for lambda_l2 in l2_grid:
        for lambda_l1 in lambda_grid:
            metrics_by_threshold: dict[float, list[dict[str, float | str]]] = {
                threshold: [] for threshold in DETECTION_THRESHOLD_GRID
            }
            for scenario in training:
                observed = observed_residual(scenario, 0.03)
                fit = fit_inverse(
                    observed,
                    ION_ORDER,
                    "thermo_elastic_net",
                    float(lambda_l1),
                    float(lambda_l2),
                    scenario["upstream_si"],
                )
                extents = np.asarray(fit["extents"])
                for threshold in DETECTION_THRESHOLD_GRID:
                    metrics_by_threshold[threshold].append(
                        evaluate_fit(
                            scenario,
                            extents,
                            observed,
                            ION_ORDER,
                            class_map,
                            support_threshold=threshold,
                        )
                    )
            for threshold, metrics in metrics_by_threshold.items():
                frame = pd.DataFrame(metrics)
                score = (
                    frame["class_f1"].mean()
                    + 0.5 * frame["phase_f1"].mean()
                    - 0.30 * frame["false_discovery_rate"].mean()
                    - 0.10 * frame["extent_rmse_mmolL"].mean()
                    - 0.05 * frame["reconstruction_rmse_mmolL"].mean()
                )
                rows.append(
                    {
                        "lambda_l1": float(lambda_l1),
                        "lambda_l2": float(lambda_l2),
                        "support_threshold_mmolL": float(threshold),
                        "n_training_scenarios": len(training),
                        "mean_phase_f1": frame["phase_f1"].mean(),
                        "mean_class_f1": frame["class_f1"].mean(),
                        "mean_phase_precision": frame["phase_precision"].mean(),
                        "mean_phase_recall": frame["phase_recall"].mean(),
                        "mean_false_discovery_rate": frame[
                            "false_discovery_rate"
                        ].mean(),
                        "mean_extent_rmse_mmolL": frame[
                            "extent_rmse_mmolL"
                        ].mean(),
                        "mean_reconstruction_rmse_mmolL": frame[
                            "reconstruction_rmse_mmolL"
                        ].mean(),
                        "selection_score": score,
                    }
                )
    selection = pd.DataFrame(rows).sort_values(
        ["selection_score", "lambda_l1"],
        ascending=[False, False],
    )
    best = selection.iloc[0]
    return {
        "lambda_l1": float(best["lambda_l1"]),
        "lambda_l2": float(best["lambda_l2"]),
        "support_threshold_mmolL": float(best["support_threshold_mmolL"]),
    }, selection


def _hyperparameters_from_selection(selection: pd.DataFrame) -> dict[str, float]:
    best = selection.sort_values(
        ["selection_score", "lambda_l1"],
        ascending=[False, False],
    ).iloc[0]
    return {
        "lambda_l1": float(best["lambda_l1"]),
        "lambda_l2": float(best["lambda_l2"]),
        "support_threshold_mmolL": float(
            best.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
        ),
    }


def _replace_method_rows(
    existing_path: Path,
    new_rows: pd.DataFrame,
    method: str,
) -> pd.DataFrame:
    if existing_path.exists():
        existing = pd.read_csv(existing_path)
        if "method" in existing.columns:
            existing = existing[existing["method"] != method]
        combined = pd.concat([existing, new_rows], ignore_index=True)
    else:
        combined = new_rows.copy()
    combined.to_csv(existing_path, index=False)
    return combined


def run_factorial_fits(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
    method_names: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    method_settings = {
        "bounded_ls": {
            "fit_method": "bounded_ls",
            "lambda_l1": 0.0,
            "lambda_l2": 0.0,
            "support_threshold": SUPPORT_THRESHOLD,
        },
        "lasso": {
            "fit_method": "lasso",
            "lambda_l1": hyperparameters["lambda_l1"],
            "lambda_l2": 0.0,
            "support_threshold": SUPPORT_THRESHOLD,
        },
        "elastic_net": {
            "fit_method": "elastic_net",
            "lambda_l1": hyperparameters["lambda_l1"],
            "lambda_l2": hyperparameters["lambda_l2"],
            "support_threshold": SUPPORT_THRESHOLD,
        },
        LEGACY_PRIMARY_METHOD: {
            "fit_method": LEGACY_PRIMARY_METHOD,
            "lambda_l1": hyperparameters["lambda_l1"],
            "lambda_l2": hyperparameters["lambda_l2"],
            "support_threshold": SUPPORT_THRESHOLD,
        },
        GUARDED_METHOD: {
            "fit_method": LEGACY_PRIMARY_METHOD,
            "lambda_l1": hyperparameters["lambda_l1"],
            "lambda_l2": hyperparameters["lambda_l2"],
            "support_threshold": hyperparameters.get(
                "support_threshold_mmolL", SUPPORT_THRESHOLD
            ),
        },
        HYDROSHEAF_CORE_METHOD: {
            "fit_method": LEGACY_PRIMARY_METHOD,
            "lambda_l1": hyperparameters["lambda_l1"],
            "lambda_l2": hyperparameters["lambda_l2"],
            "support_threshold": hyperparameters.get(
                "support_threshold_mmolL", SUPPORT_THRESHOLD
            ),
            "use_core_evidence": True,
        },
    }
    if method_names is not None:
        selected_methods = set(method_names)
        method_settings = {
            method: settings
            for method, settings in method_settings.items()
            if method in selected_methods
        }
        missing = selected_methods - set(method_settings)
        if missing:
            raise ValueError(f"Unknown M5 method(s): {sorted(missing)}")
    fit_rows: list[dict[str, object]] = []
    recovery_rows: list[dict[str, object]] = []
    heldout_rows: list[dict[str, object]] = []
    core_evidence_rows: list[dict[str, object]] = []
    full_matrix = reaction_matrix()
    for scenario in scenarios:
        for noise_level in NOISE_LEVELS:
            observed = observed_residual(scenario, noise_level)
            for panel_name, panel_ions in PANELS.items():
                indices = [ION_ORDER.index(ion) for ion in panel_ions]
                hidden = [ion for ion in ION_ORDER if ion not in panel_ions]
                fit_cache: dict[tuple[object, ...], dict[str, object]] = {}
                for method, settings in method_settings.items():
                    fit_method = str(settings["fit_method"])
                    lambda_l1 = float(settings["lambda_l1"])
                    lambda_l2 = float(settings["lambda_l2"])
                    support_threshold = float(settings["support_threshold"])
                    penalty_scales = None
                    evidence_rows: list[dict[str, object]] = []
                    if bool(settings.get("use_core_evidence", False)):
                        penalty_scales, evidence_rows = hydrosheaf_core_evidence(
                            np.asarray(scenario["upstream"], dtype=float),
                            observed[indices],
                            panel_ions,
                            scenario["upstream_si"],
                        )
                    penalty_signature = (
                        tuple(np.round(penalty_scales, 6))
                        if penalty_scales is not None
                        else ()
                    )
                    fit_key = (fit_method, lambda_l1, lambda_l2, penalty_signature)
                    if fit_key not in fit_cache:
                        fit_cache[fit_key] = fit_inverse(
                            observed[indices],
                            panel_ions,
                            fit_method,
                            lambda_l1,
                            lambda_l2,
                            scenario["upstream_si"],
                            SI_THRESHOLD,
                            penalty_scales=penalty_scales,
                        )
                    fit = fit_cache[fit_key]
                    extents = np.asarray(fit["extents"], dtype=float)
                    metrics = evaluate_fit(
                        scenario,
                        extents,
                        observed,
                        panel_ions,
                        class_map,
                        support_threshold=support_threshold,
                    )
                    row = {
                        "scenario_id": scenario["scenario_id"],
                        "archetype": scenario["archetype"],
                        "replicate": scenario["replicate"],
                        "noise_level": noise_level,
                        "transport_error_level": scenario["transport_error_level"],
                        "panel": panel_name,
                        "n_measured_ions": len(panel_ions),
                        "method": method,
                        "fit_method": fit_method,
                        "lambda_l1": lambda_l1,
                        "lambda_l2": lambda_l2,
                        "support_threshold_mmolL": support_threshold,
                        "solver_rmse_mmolL": fit["rmse"],
                        "l1_norm": fit["l1_norm"],
                        "iterations": fit["iterations"],
                        "converged": fit["converged"],
                        "runtime_ms": fit["runtime_ms"],
                        "thermodynamic_violations": fit[
                            "thermodynamic_violations"
                        ],
                        "topology_confidence": scenario["topology_confidence"],
                        "residence_time_confidence": scenario[
                            "residence_time_confidence"
                        ],
                        **metrics,
                    }
                    fit_rows.append(row)
                    truth_extents = np.asarray(scenario["true_extents"], dtype=float)
                    if evidence_rows:
                        for evidence, reaction, truth, recovered in zip(
                            evidence_rows,
                            REACTIONS,
                            truth_extents,
                            extents,
                        ):
                            core_evidence_rows.append(
                                {
                                    "scenario_id": scenario["scenario_id"],
                                    "archetype": scenario["archetype"],
                                    "replicate": scenario["replicate"],
                                    "noise_level": noise_level,
                                    "transport_error_level": scenario[
                                        "transport_error_level"
                                    ],
                                    "panel": panel_name,
                                    "method": method,
                                    "true_active": abs(truth) >= SUPPORT_THRESHOLD,
                                    "recovered_active": abs(recovered)
                                    >= support_threshold,
                                    "true_extent_mmolL": truth,
                                    "recovered_extent_mmolL": recovered,
                                    "equivalence_class": class_map[reaction.label],
                                    **evidence,
                                }
                            )
                    for reaction, truth, recovered in zip(
                        REACTIONS, truth_extents, extents
                    ):
                        recovery_rows.append(
                            {
                                "scenario_id": scenario["scenario_id"],
                                "archetype": scenario["archetype"],
                                "noise_level": noise_level,
                                "panel": panel_name,
                                "method": method,
                                "reaction": reaction.label,
                                "family": reaction.family,
                                "true_extent_mmolL": truth,
                                "recovered_extent_mmolL": recovered,
                                "absolute_error_mmolL": abs(recovered - truth),
                                "true_active": abs(truth) >= SUPPORT_THRESHOLD,
                                "recovered_active": abs(recovered)
                                >= support_threshold,
                                "direction_correct": (
                                    np.sign(recovered) == np.sign(truth)
                                    if abs(truth) >= SUPPORT_THRESHOLD
                                    and abs(recovered) >= support_threshold
                                    else ""
                                ),
                                "equivalence_class": class_map[reaction.label],
                                "support_threshold_mmolL": support_threshold,
                            }
                        )
                    full_prediction = full_matrix.T @ extents
                    for ion in hidden:
                        ion_index = ION_ORDER.index(ion)
                        heldout_rows.append(
                            {
                                "scenario_id": scenario["scenario_id"],
                                "archetype": scenario["archetype"],
                                "noise_level": noise_level,
                                "panel": panel_name,
                                "method": method,
                                "fit_method": fit_method,
                                "support_threshold_mmolL": support_threshold,
                                "heldout_ion": ion,
                                "true_delta_mmolL": np.asarray(
                                    scenario["true_delta"]
                                )[ion_index],
                                "predicted_delta_mmolL": full_prediction[ion_index],
                                "absolute_error_mmolL": abs(
                                    np.asarray(scenario["true_delta"])[ion_index]
                                    - full_prediction[ion_index]
                                ),
                            }
                        )
    return (
        pd.DataFrame(fit_rows),
        pd.DataFrame(recovery_rows),
        pd.DataFrame(heldout_rows),
        pd.DataFrame(core_evidence_rows),
    )


def run_data_tier_experiment(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compare sparse-data, Plus-lite, and Enhanced evidence tiers.

    Plus-lite and Enhanced diagnostics are generated as a controlled synthetic
    measurement-design experiment from reaction truth plus measurement noise.
    They are not field-measured reaction truth.
    """
    panel_ions = PANELS[DATA_TIER_PANEL]
    indices = [ION_ORDER.index(ion) for ion in panel_ions]
    support_threshold = float(
        hyperparameters.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
    )
    metric_rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    diagnostic_rows: list[dict[str, object]] = []
    for scenario in scenarios:
        observed = observed_residual(scenario, DATA_TIER_NOISE_LEVEL)
        core_scales, core_evidence_rows = hydrosheaf_core_evidence(
            np.asarray(scenario["upstream"], dtype=float),
            observed[indices],
            panel_ions,
            scenario["upstream_si"],
        )
        core_by_reaction = {
            str(row["reaction"]): row for row in core_evidence_rows
        }
        for tier, diagnostics in DATA_TIER_DIAGNOSTICS.items():
            diagnostic_names, clean_values, observed_values = _simulate_optional_diagnostics(
                scenario,
                tier,
                DATA_TIER_NOISE_LEVEL,
            )
            optional_scales, optional_rows = _optional_diagnostic_evidence(
                diagnostic_names,
                observed_values,
            )
            optional_by_reaction = {
                str(row["reaction"]): row for row in optional_rows
            }
            if tier == "core":
                penalty_scales = core_scales
            else:
                penalty_scales = _combine_core_and_optional_penalties(
                    core_scales,
                    optional_scales,
                )
            fit = fit_inverse(
                observed[indices],
                panel_ions,
                LEGACY_PRIMARY_METHOD,
                hyperparameters["lambda_l1"],
                hyperparameters["lambda_l2"],
                scenario["upstream_si"],
                SI_THRESHOLD,
                penalty_scales=penalty_scales,
            )
            extents = np.asarray(fit["extents"], dtype=float)
            metrics = evaluate_fit(
                scenario,
                extents,
                observed,
                panel_ions,
                class_map,
                support_threshold=support_threshold,
            )
            metric_rows.append(
                {
                    "scenario_id": scenario["scenario_id"],
                    "archetype": scenario["archetype"],
                    "replicate": scenario["replicate"],
                    "data_tier": tier,
                    "panel": DATA_TIER_PANEL,
                    "noise_level": DATA_TIER_NOISE_LEVEL,
                    "n_major_ions": len(panel_ions),
                    "n_optional_diagnostics": len(diagnostic_names),
                    "optional_diagnostics": ";".join(diagnostic_names),
                    "diagnostic_source": (
                        "major-ion evidence only"
                        if tier == "core"
                        else "controlled synthetic optional diagnostics"
                    ),
                    "solver_rmse_mmolL": fit["rmse"],
                    "l1_norm": fit["l1_norm"],
                    "iterations": fit["iterations"],
                    "converged": fit["converged"],
                    "thermodynamic_violations": fit[
                        "thermodynamic_violations"
                    ],
                    **metrics,
                }
            )
            truth_extents = np.asarray(scenario["true_extents"], dtype=float)
            for reaction, truth, recovered, combined_scale in zip(
                REACTIONS,
                truth_extents,
                extents,
                penalty_scales,
            ):
                core_row = core_by_reaction[reaction.label]
                optional_row = optional_by_reaction[reaction.label]
                combined_evidence_score = _combine_evidence_scores(
                    core_row["hydrosheaf_core_evidence_score"],
                    optional_row["optional_evidence_score"],
                )
                evidence_rows.append(
                    {
                        "scenario_id": scenario["scenario_id"],
                        "archetype": scenario["archetype"],
                        "replicate": scenario["replicate"],
                        "data_tier": tier,
                        "reaction": reaction.label,
                        "family": reaction.family,
                        "equivalence_class": class_map[reaction.label],
                        "true_extent_mmolL": truth,
                        "recovered_extent_mmolL": recovered,
                        "true_active": abs(truth) >= SUPPORT_THRESHOLD,
                        "recovered_active": abs(recovered) >= support_threshold,
                        "core_evidence_score": core_row[
                            "hydrosheaf_core_evidence_score"
                        ],
                        "core_penalty_scale": core_row["penalty_scale"],
                        "optional_evidence_score": optional_row[
                            "optional_evidence_score"
                        ],
                        "combined_evidence_score": combined_evidence_score,
                        "optional_penalty_scale": optional_row[
                            "optional_penalty_scale"
                        ],
                        "combined_penalty_scale": float(combined_scale),
                        "optional_diagnostic_alignment": optional_row[
                            "optional_diagnostic_alignment"
                        ],
                        "optional_support_alignment": optional_row[
                            "optional_support_alignment"
                        ],
                        "diagnostics_with_nonzero_signature": optional_row[
                            "diagnostics_with_nonzero_signature"
                        ],
                    }
                )
            for name, clean, measured in zip(
                diagnostic_names,
                clean_values,
                observed_values,
            ):
                diagnostic_rows.append(
                    {
                        "scenario_id": scenario["scenario_id"],
                        "archetype": scenario["archetype"],
                        "replicate": scenario["replicate"],
                        "data_tier": tier,
                        "diagnostic": name,
                        "clean_value": float(clean),
                        "observed_value": float(measured),
                        "absolute_noise": float(measured - clean),
                    }
                )
    return (
        pd.DataFrame(metric_rows),
        pd.DataFrame(evidence_rows),
        pd.DataFrame(diagnostic_rows),
    )


def write_data_tier_outputs(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> pd.DataFrame:
    print("M5: running Core/Plus-lite/Enhanced data-tier experiment...")
    data_tiers, data_tier_evidence, data_tier_diagnostics = run_data_tier_experiment(
        scenarios,
        hyperparameters,
        class_map,
    )
    data_tiers.to_csv(RESULTS_DIR / "data_tier_experiment.csv", index=False)
    data_tier_evidence.to_csv(
        RESULTS_DIR / "data_tier_reaction_evidence.csv", index=False
    )
    data_tier_resolution = _evidence_lifted_resolution_frame(
        data_tier_evidence,
        class_map,
        score_column="combined_evidence_score",
        group_columns=[
            "scenario_id",
            "archetype",
            "replicate",
            "data_tier",
        ],
        evidence_source="data_tier_combined_evidence",
    )
    data_tier_resolution.to_csv(
        RESULTS_DIR / "data_tier_evidence_lifted_resolution.csv",
        index=False,
    )
    data_tier_diagnostics.to_csv(
        RESULTS_DIR / "data_tier_optional_diagnostics.csv", index=False
    )
    return data_tiers


def add_mechanism_resolution_scores(
    fits: pd.DataFrame,
    class_map: Mapping[str, str],
) -> pd.DataFrame:
    diagnostics = {name: matrix_diagnostics(ions) for name, ions in PANELS.items()}
    def support_set(value: object) -> set[str]:
        if pd.isna(value) or not str(value).strip():
            return set()
        return set(str(value).split(";"))

    support_lookup: dict[tuple[object, ...], list[set[str]]] = {}
    for key, group in fits.groupby(["scenario_id", "noise_level", "panel"]):
        support_lookup[key] = [support_set(value) for value in group["selected_support"]]
    rows: list[dict[str, object]] = []
    primary = fits[fits["method"] == PRIMARY_METHOD].copy()
    if primary.empty:
        primary = fits[fits["method"] == LEGACY_PRIMARY_METHOD].copy()
    for _, row in primary.iterrows():
        key = (row["scenario_id"], row["noise_level"], row["panel"])
        supports = support_lookup[key]
        selected = support_set(row["selected_support"])
        stability = float(np.mean([jaccard(selected, support) for support in supports]))
        diagnostic = diagnostics[str(row["panel"])]
        rank_skill = diagnostic["rank"] / max(
            1.0,
            min(
                diagnostic["n_ions"],
                diagnostic["n_equivalence_classes"],
            ),
        )
        coherence_skill = 1.0 - min(
            1.0, diagnostic["class_collapsed_coherence"]
        )
        scale = 0.05 + abs(float(row["extent_mae_mmolL"]))
        predictive_error = (
            float(row["heldout_rmse_mmolL"])
            if np.isfinite(row["heldout_rmse_mmolL"])
            else float(row["reconstruction_rmse_mmolL"])
        )
        predictive_skill = math.exp(-predictive_error / scale)
        thermo_skill = 1.0 - min(
            1.0,
            float(row["thermodynamic_violations"]) / max(len(selected), 1),
        )
        uncalibrated_mrs = 100.0 * (
            0.25 * rank_skill
            + 0.15 * coherence_skill
            + 0.25 * stability
            + 0.20 * predictive_skill
            + 0.15 * thermo_skill
        )
        true_class = (
            "identifiable"
            if row["phase_f1"] >= 0.80
            else "equivalence_class"
            if row["class_f1"] >= 0.80
            else "partially_identifiable"
            if row["class_f1"] >= 0.50
            else "non_identifiable"
        )
        rows.append(
            {
                **row.to_dict(),
                "rank_skill": rank_skill,
                "coherence_skill": coherence_skill,
                "support_stability": stability,
                "predictive_skill": predictive_skill,
                "thermodynamic_skill": thermo_skill,
                "uncalibrated_mrs": uncalibrated_mrs,
                "true_resolution_class": true_class,
            }
        )
    frame = pd.DataFrame(rows)
    features = frame[MRS_FEATURES].copy()
    features["heldout_rmse_mmolL"] = features["heldout_rmse_mmolL"].fillna(
        features["reconstruction_rmse_mmolL"]
    )
    values = np.nan_to_num(features.to_numpy(dtype=float))
    training = frame["archetype"].ne("mixed").to_numpy()
    means = values[training].mean(axis=0)
    scales = values[training].std(axis=0)
    scales[scales == 0.0] = 1.0
    design = np.column_stack([np.ones(len(values)), (values - means) / scales])
    penalty = np.eye(design.shape[1])
    penalty[0, 0] = 0.0
    outcomes = np.column_stack(
        [
            frame["true_resolution_class"].eq(label).astype(float)
            for label in RESOLUTION_CLASSES
        ]
    )
    class_coefficients = np.linalg.solve(
        design[training].T @ design[training] + penalty,
        design[training].T @ outcomes[training],
    )
    class_scores = design @ class_coefficients
    predicted = np.asarray(RESOLUTION_CLASSES)[np.argmax(class_scores, axis=1)]
    reliability_coefficients = np.linalg.solve(
        design[training].T @ design[training] + penalty,
        design[training].T @ frame.loc[training, "class_f1"].to_numpy(),
    )
    calibrated = np.clip(design @ reliability_coefficients, 0.0, 1.0)
    frame["mechanism_resolution_score"] = 100.0 * calibrated
    frame["predicted_resolution_class"] = predicted
    frame["resolution_class_correct"] = (
        frame["true_resolution_class"] == frame["predicted_resolution_class"]
    )
    frame["edge_conditioned_mrs"] = frame["mechanism_resolution_score"] * np.sqrt(
        frame["topology_confidence"] * frame["residence_time_confidence"]
    )
    calibration = {
        "training_archetypes": ["carbonate", "crystalline", "evaporitic"],
        "heldout_archetype": "mixed",
        "features": MRS_FEATURES,
        "feature_means": means.tolist(),
        "feature_scales": scales.tolist(),
        "resolution_classes": RESOLUTION_CLASSES,
        "class_coefficients": class_coefficients.tolist(),
        "reliability_coefficients": reliability_coefficients.tolist(),
        "ridge_penalty": 1.0,
    }
    (RESULTS_DIR / "mrs_calibration_model.json").write_text(
        json.dumps(calibration, indent=2), encoding="utf-8"
    )
    return frame


def run_regularization_paths(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
) -> pd.DataFrame:
    selected = [
        next(
            scenario
            for scenario in scenarios
            if scenario["archetype"] == archetype and scenario["replicate"] == 1
        )
        for archetype in ARCHETYPES
    ]
    rows: list[dict[str, object]] = []
    for scenario in selected:
        observed = observed_residual(scenario, 0.03)
        for lambda_l1 in np.geomspace(0.0001, 0.20, 24):
            fit = fit_inverse(
                observed,
                ION_ORDER,
                "elastic_net",
                float(lambda_l1),
                hyperparameters["lambda_l2"],
                scenario["upstream_si"],
            )
            extents = np.asarray(fit["extents"])
            row = {
                "scenario_id": scenario["scenario_id"],
                "archetype": scenario["archetype"],
                "lambda_l1": float(lambda_l1),
                "lambda_l2": hyperparameters["lambda_l2"],
                "rmse_mmolL": fit["rmse"],
                "n_selected": len(
                    support_from_extents(
                        extents,
                        hyperparameters.get(
                            "support_threshold_mmolL", SUPPORT_THRESHOLD
                        ),
                    )
                ),
            }
            for label, extent in zip(REACTION_LABELS, extents):
                row[label] = extent
            rows.append(row)
    return pd.DataFrame(rows)


def run_next_best_measurement(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> pd.DataFrame:
    selected = [
        scenario for scenario in scenarios if int(scenario["replicate"]) <= 10
    ]
    support_threshold = float(
        hyperparameters.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
    )
    rows: list[dict[str, object]] = []
    for scenario in selected:
        observed = observed_residual(scenario, 0.03)
        full_fit = fit_inverse(
            observed,
            ION_ORDER,
            "thermo_elastic_net",
            hyperparameters["lambda_l1"],
            hyperparameters["lambda_l2"],
            scenario["upstream_si"],
        )
        full_metrics = evaluate_fit(
            scenario,
            np.asarray(full_fit["extents"]),
            observed,
            ION_ORDER,
            class_map,
            support_threshold=support_threshold,
        )
        full_classes = set(str(full_metrics["selected_classes"]).split(";"))
        for ion in ION_ORDER:
            panel = [value for value in ION_ORDER if value != ion]
            indices = [ION_ORDER.index(value) for value in panel]
            fit = fit_inverse(
                observed[indices],
                panel,
                "thermo_elastic_net",
                hyperparameters["lambda_l1"],
                hyperparameters["lambda_l2"],
                scenario["upstream_si"],
            )
            metrics = evaluate_fit(
                scenario,
                np.asarray(fit["extents"]),
                observed,
                panel,
                class_map,
                support_threshold=support_threshold,
            )
            omitted_classes = set(str(metrics["selected_classes"]).split(";"))
            hidden_error = float(metrics["heldout_rmse_mmolL"])
            class_gain = float(full_metrics["class_f1"]) - float(metrics["class_f1"])
            support_change = 1.0 - jaccard(full_classes, omitted_classes)
            scale = 0.02 + abs(
                float(np.asarray(scenario["true_delta"])[ION_ORDER.index(ion)])
            )
            value_score = hidden_error / scale + max(0.0, class_gain) + support_change
            rows.append(
                {
                    "scenario_id": scenario["scenario_id"],
                    "archetype": scenario["archetype"],
                    "candidate_measurement": ion,
                    "heldout_absolute_error_mmolL": hidden_error,
                    "class_f1_without_measurement": metrics["class_f1"],
                    "class_f1_with_measurement": full_metrics["class_f1"],
                    "realised_class_f1_gain": class_gain,
                    "support_change": support_change,
                    "measurement_value_score": value_score,
                }
            )
    return pd.DataFrame(rows)


def run_bootstrap_stability(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
) -> pd.DataFrame:
    selected = [
        scenario for scenario in scenarios if int(scenario["replicate"]) <= 3
    ]
    counts: dict[tuple[str, str], int] = {}
    extent_values: dict[tuple[str, str], list[float]] = {}
    n_bootstrap = 40
    support_threshold = float(
        hyperparameters.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
    )
    for scenario in selected:
        for bootstrap in range(n_bootstrap):
            rng = np.random.default_rng(
                RANDOM_SEED + 900000 + int(scenario["scenario_index"]) * 100 + bootstrap
            )
            truth = np.asarray(scenario["true_delta"])
            scale = 0.03 * np.maximum(
                np.asarray(scenario["downstream"]), 0.02
            )
            residual = truth + rng.normal(0.0, scale)
            fit = fit_inverse(
                residual,
                ION_ORDER,
                "thermo_elastic_net",
                hyperparameters["lambda_l1"],
                hyperparameters["lambda_l2"],
                scenario["upstream_si"],
            )
            for label, extent in zip(REACTION_LABELS, fit["extents"]):
                key = (str(scenario["scenario_id"]), label)
                counts[key] = counts.get(key, 0) + (
                    abs(float(extent)) >= support_threshold
                )
                extent_values.setdefault(key, []).append(float(extent))
    rows = []
    scenario_lookup = {str(item["scenario_id"]): item for item in selected}
    for (scenario_id, label), values in extent_values.items():
        rows.append(
            {
                "scenario_id": scenario_id,
                "archetype": scenario_lookup[scenario_id]["archetype"],
                "reaction": label,
                "n_bootstrap": n_bootstrap,
                "support_threshold_mmolL": support_threshold,
                "support_frequency": counts[(scenario_id, label)] / n_bootstrap,
                "mean_extent_mmolL": np.mean(values),
                "sd_extent_mmolL": np.std(values, ddof=1),
            }
        )
    return pd.DataFrame(rows)


def run_thermodynamic_threshold_sensitivity(
    scenarios: Sequence[Mapping[str, object]],
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> pd.DataFrame:
    selected = [
        scenario for scenario in scenarios if int(scenario["replicate"]) <= 8
    ]
    support_threshold = float(
        hyperparameters.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
    )
    rows: list[dict[str, object]] = []
    for threshold in [0.0, 0.1, 0.3]:
        for scenario in selected:
            observed = observed_residual(scenario, 0.03)
            fit = fit_inverse(
                observed,
                ION_ORDER,
                "thermo_elastic_net",
                hyperparameters["lambda_l1"],
                hyperparameters["lambda_l2"],
                scenario["upstream_si"],
                threshold,
            )
            metrics = evaluate_fit(
                scenario,
                np.asarray(fit["extents"]),
                observed,
                ION_ORDER,
                class_map,
                support_threshold=support_threshold,
            )
            rows.append(
                {
                    "scenario_id": scenario["scenario_id"],
                    "archetype": scenario["archetype"],
                    "si_threshold": threshold,
                    "support_threshold_mmolL": support_threshold,
                    "thermodynamic_violations": fit[
                        "thermodynamic_violations"
                    ],
                    **metrics,
                }
            )
    return pd.DataFrame(rows)


def load_ghana_pairs() -> tuple[pd.DataFrame, pd.DataFrame]:
    hydro = pd.read_excel(GHANA_WORKBOOK, sheet_name="Hydrochemistry_Seasonal")
    wells = pd.read_excel(GHANA_WORKBOOK, sheet_name="Wells_Nodes")
    hydro["Sampling_Date"] = pd.to_datetime(hydro["Sampling_Date"])
    years = sorted(hydro["Sampling_Date"].dropna().dt.year.unique().tolist())
    if years != [2025]:
        raise ValueError(f"Northern Ghana sampling years must be [2025], found {years}.")
    metadata = wells[
        [
            "Well_ID",
            "Aquifer_Type",
            "Geology_Group",
            "Lithology",
            "Region",
            "District",
        ]
    ].drop_duplicates("Well_ID")
    hydro = hydro.drop(columns=["Aquifer_Type"], errors="ignore").merge(
        metadata,
        on=["Well_ID", "Region"],
        how="left",
    )
    return hydro, wells


def _ghana_vector(row: pd.Series) -> tuple[np.ndarray, list[str]]:
    values: list[float] = []
    measured: list[str] = []
    for ion in ION_ORDER:
        column = f"{ion}_mg_L"
        if column in row.index and pd.notna(row[column]):
            values.append(float(row[column]) / MOLAR_MASS_G_MOL[ion])
            measured.append(ion)
        else:
            values.append(np.nan)
    return np.asarray(values, dtype=float), measured


def run_ghana_field_demonstration(
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hydro, _ = load_ghana_pairs()
    support_threshold = float(
        hyperparameters.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
    )
    calibration_path = RESULTS_DIR / "mrs_calibration_model.json"
    calibration = json.loads(calibration_path.read_text(encoding="utf-8"))
    feature_means = np.asarray(calibration["feature_means"], dtype=float)
    feature_scales = np.asarray(calibration["feature_scales"], dtype=float)
    class_coefficients = np.asarray(
        calibration["class_coefficients"], dtype=float
    )
    reliability_coefficients = np.asarray(
        calibration["reliability_coefficients"], dtype=float
    )
    quantitative = hydro[
        hydro["Data_Class"].astype(str).str.contains(
            "Quantitative inverse modelling", case=False, na=False
        )
    ]
    ec_tds_models = _calibrate_field_ec_tds(quantitative)
    pair_rows: list[dict[str, object]] = []
    extent_rows: list[dict[str, object]] = []
    heldout_rows: list[dict[str, object]] = []
    class_rows: list[dict[str, object]] = []
    field_evidence_rows: list[dict[str, object]] = []
    for well_id, group in quantitative.groupby("Well_ID"):
        seasons = {
            str(row["Season"]).strip().lower(): row
            for _, row in group.iterrows()
        }
        if "wet" not in seasons or "dry" not in seasons:
            continue
        wet = seasons["wet"]
        dry = seasons["dry"]
        wet_vector, wet_measured = _ghana_vector(wet)
        dry_vector, dry_measured = _ghana_vector(dry)
        measured = [
            ion
            for ion in ION_ORDER
            if ion in wet_measured and ion in dry_measured
        ]
        if len(measured) < 6:
            continue
        indices = [ION_ORDER.index(ion) for ion in measured]
        residual = dry_vector[indices] - wet_vector[indices]
        upstream_si = {
            "Calcite": float(wet.get("Calcite_SI", np.nan)),
            "Dolomite": float(wet.get("Dolomite_SI", np.nan)),
            "Gypsum": float(wet.get("Gypsum_SI", np.nan)),
            "Halite": float(wet.get("Halite_SI", np.nan)),
        }
        penalty_scales, evidence_rows = hydrosheaf_core_evidence(
            wet_vector,
            residual,
            measured,
            upstream_si,
        )
        method_fits = {}
        for method, l2 in [
            ("lasso", 0.0),
            ("elastic_net", hyperparameters["lambda_l2"]),
            ("thermo_elastic_net", hyperparameters["lambda_l2"]),
        ]:
            method_fits[method] = fit_inverse(
                residual,
                measured,
                method,
                hyperparameters["lambda_l1"],
                l2,
                upstream_si,
            )
        method_fits[GUARDED_METHOD] = method_fits[LEGACY_PRIMARY_METHOD]
        method_fits[HYDROSHEAF_CORE_METHOD] = fit_inverse(
            residual,
            measured,
            LEGACY_PRIMARY_METHOD,
            hyperparameters["lambda_l1"],
            hyperparameters["lambda_l2"],
            upstream_si,
            SI_THRESHOLD,
            penalty_scales=penalty_scales,
        )
        fit = method_fits[PRIMARY_METHOD]
        extents = np.asarray(fit["extents"])
        selected = support_from_extents(extents, support_threshold)
        selected_classes = class_set(selected, class_map)
        supports = [
            support_from_extents(value["extents"], support_threshold)
            for value in method_fits.values()
        ]
        stability = float(np.mean([jaccard(selected, support) for support in supports]))
        diagnostic = matrix_diagnostics(measured)
        leave_one_errors = []
        leave_one_changes = []
        for heldout in measured:
            panel = [ion for ion in measured if ion != heldout]
            panel_indices = [measured.index(ion) for ion in panel]
            leave_fit = fit_inverse(
                residual[panel_indices],
                panel,
                "thermo_elastic_net",
                hyperparameters["lambda_l1"],
                hyperparameters["lambda_l2"],
                upstream_si,
            )
            prediction = reaction_matrix(measured).T @ np.asarray(
                leave_fit["extents"]
            )
            heldout_index = measured.index(heldout)
            error = abs(residual[heldout_index] - prediction[heldout_index])
            leave_support = support_from_extents(
                leave_fit["extents"], support_threshold
            )
            support_change = 1.0 - jaccard(selected, leave_support)
            leave_one_errors.append(error)
            leave_one_changes.append(support_change)
            heldout_rows.append(
                {
                    "well_id": well_id,
                    "aquifer": wet["Aquifer_Type"],
                    "lithology": wet["Lithology"],
                    "heldout_ion": heldout,
                    "observed_delta_mmolL": residual[heldout_index],
                    "predicted_delta_mmolL": prediction[heldout_index],
                    "absolute_error_mmolL": error,
                    "support_change": support_change,
                    "measurement_value_score": error
                    / (0.02 + abs(residual[heldout_index]))
                    + support_change,
                }
            )
        mean_heldout = float(np.mean(leave_one_errors))
        rank_skill = diagnostic["rank"] / max(
            1.0,
            min(
                diagnostic["n_ions"],
                diagnostic["n_equivalence_classes"],
            ),
        )
        coherence_skill = 1.0 - min(
            1.0, diagnostic["class_collapsed_coherence"]
        )
        predictive_skill = math.exp(
            -mean_heldout / (0.05 + float(np.mean(np.abs(residual))))
        )
        thermo_skill = 1.0 - min(
            1.0,
            float(fit["thermodynamic_violations"]) / max(len(selected), 1),
        )
        uncalibrated_mrs = 100.0 * (
            0.25 * rank_skill
            + 0.15 * coherence_skill
            + 0.25 * stability
            + 0.20 * predictive_skill
            + 0.15 * thermo_skill
        )
        feature_values = np.asarray(
            [
                rank_skill,
                coherence_skill,
                stability,
                predictive_skill,
                thermo_skill,
                fit["rmse"],
                mean_heldout,
                len(measured),
                fit["l1_norm"],
            ],
            dtype=float,
        )
        calibrated_design = np.concatenate(
            [[1.0], (feature_values - feature_means) / feature_scales]
        )
        mrs = 100.0 * float(
            np.clip(calibrated_design @ reliability_coefficients, 0.0, 1.0)
        )
        ambiguous = any(
            sum(1 for value in class_map.values() if value == class_id) > 1
            for class_id in selected_classes
        )
        resolution = RESOLUTION_CLASSES[
            int(np.argmax(calibrated_design @ class_coefficients))
        ]
        if resolution == "identifiable" and ambiguous:
            resolution = "equivalence_class"
        wet_sum_mg = _sum_major_ions_mg(wet)
        dry_sum_mg = _sum_major_ions_mg(dry)
        ec_slope, ec_intercept = ec_tds_models["ec"]
        tds_slope, tds_intercept = ec_tds_models["tds"]
        predicted_wet_ec = ec_slope * wet_sum_mg + ec_intercept
        predicted_dry_ec = ec_slope * dry_sum_mg + ec_intercept
        predicted_wet_tds = tds_slope * wet_sum_mg + tds_intercept
        predicted_dry_tds = tds_slope * dry_sum_mg + tds_intercept
        wet_ec = _finite_float(wet.get("EC_uS_cm"))
        dry_ec = _finite_float(dry.get("EC_uS_cm"))
        wet_tds = _finite_float(wet.get("TDS_mg_L"))
        dry_tds = _finite_float(dry.get("TDS_mg_L"))
        observed_ec_delta = (
            dry_ec - wet_ec if wet_ec is not None and dry_ec is not None else None
        )
        observed_tds_delta = (
            dry_tds - wet_tds if wet_tds is not None and dry_tds is not None else None
        )
        optional_tracers = _field_optional_tracer_metrics(wet, dry)
        evidence_scores = [float(row["hydrosheaf_core_evidence_score"]) for row in evidence_rows]
        conflict_count = sum(
            "conflict" in str(row["evidence_flags"]) for row in evidence_rows
        )
        pair_rows.append(
            {
                "well_id": well_id,
                "aquifer": wet["Aquifer_Type"],
                "geology_group": wet["Geology_Group"],
                "lithology": wet["Lithology"],
                "region": wet["Region"],
                "district": wet["District"],
                "wet_date": pd.Timestamp(wet["Sampling_Date"]).date().isoformat(),
                "dry_date": pd.Timestamp(dry["Sampling_Date"]).date().isoformat(),
                "n_measured_ions": len(measured),
                "support_threshold_mmolL": support_threshold,
                "selected_support": ";".join(sorted(selected)),
                "selected_classes": ";".join(sorted(selected_classes)),
                "support_stability": stability,
                "mean_heldout_rmse_mmolL": mean_heldout,
                "mean_leave_one_support_change": float(
                    np.mean(leave_one_changes)
                ),
                "mean_hydrosheaf_core_evidence_score": float(
                    np.mean(evidence_scores)
                ),
                "max_hydrosheaf_core_evidence_score": float(
                    np.max(evidence_scores)
                ),
                "hydrosheaf_core_conflict_count": int(conflict_count),
                "observed_ec_delta_uS_cm": (
                    observed_ec_delta if observed_ec_delta is not None else np.nan
                ),
                "predicted_ec_delta_uS_cm": predicted_dry_ec - predicted_wet_ec,
                "ec_delta_consistency_score": _delta_consistency_score(
                    observed_ec_delta,
                    predicted_dry_ec - predicted_wet_ec,
                    50.0,
                ),
                "observed_tds_delta_mg_L": (
                    observed_tds_delta if observed_tds_delta is not None else np.nan
                ),
                "predicted_tds_delta_mg_L": predicted_dry_tds - predicted_wet_tds,
                "tds_delta_consistency_score": _delta_consistency_score(
                    observed_tds_delta,
                    predicted_dry_tds - predicted_wet_tds,
                    25.0,
                ),
                "uncalibrated_mrs": uncalibrated_mrs,
                "mechanism_resolution_score": mrs,
                "resolution_class": resolution,
                "solver_rmse_mmolL": fit["rmse"],
                "thermodynamic_violations": fit[
                    "thermodynamic_violations"
                ],
                "converged": fit["converged"],
                **optional_tracers,
            }
        )
        for evidence, reaction, extent in zip(evidence_rows, REACTIONS, extents):
            field_evidence_rows.append(
                {
                    "well_id": well_id,
                    "aquifer": wet["Aquifer_Type"],
                    "geology_group": wet["Geology_Group"],
                    "lithology": wet["Lithology"],
                    "region": wet["Region"],
                    "district": wet["District"],
                    "selected": reaction.label in selected,
                    "extent_mmolL": extent,
                    "equivalence_class": class_map[reaction.label],
                    **evidence,
                }
            )
        for reaction, extent in zip(REACTIONS, extents):
            extent_rows.append(
                {
                    "well_id": well_id,
                    "aquifer": wet["Aquifer_Type"],
                    "lithology": wet["Lithology"],
                    "reaction": reaction.label,
                    "family": reaction.family,
                    "extent_mmolL": extent,
                    "active": abs(extent) >= support_threshold,
                    "equivalence_class": class_map[reaction.label],
                    "support_threshold_mmolL": support_threshold,
                }
            )
        for class_id in sorted(set(class_map.values())):
            class_rows.append(
                {
                    "well_id": well_id,
                    "aquifer": wet["Aquifer_Type"],
                    "lithology": wet["Lithology"],
                    "equivalence_class": class_id,
                    "selected": class_id in selected_classes,
                }
            )
    return (
        pd.DataFrame(pair_rows),
        pd.DataFrame(extent_rows),
        pd.DataFrame(heldout_rows),
        pd.DataFrame(class_rows),
        pd.DataFrame(field_evidence_rows),
    )


def _external_field_edge_outputs(
    *,
    dataset: str,
    edge_id: str,
    upstream: Mapping[str, object],
    downstream: Mapping[str, object],
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    support_threshold = float(
        hyperparameters.get("support_threshold_mmolL", SUPPORT_THRESHOLD)
    )
    upstream_vector = np.asarray(upstream["vector"], dtype=float)
    downstream_vector = np.asarray(downstream["vector"], dtype=float)
    measured = [
        ion
        for ion in ION_ORDER
        if ion in upstream["measured"] and ion in downstream["measured"]
    ]
    if len(measured) < 6:
        return [], []

    indices = [ION_ORDER.index(ion) for ion in measured]
    residual = downstream_vector[indices] - upstream_vector[indices]
    core_scales, core_rows = hydrosheaf_core_evidence(
        upstream_vector,
        residual,
        measured,
        {},
    )
    optional_diagnostics, optional_observed = _field_optional_diagnostics_from_rows(
        upstream["row"],
        downstream["row"],
        upstream["optional_columns"],
    )
    optional_scales, optional_rows = _optional_diagnostic_evidence(
        optional_diagnostics,
        optional_observed,
    )
    combined_scales = _combine_core_and_optional_penalties(
        core_scales,
        optional_scales,
    )
    core_by_reaction = {str(row["reaction"]): row for row in core_rows}
    optional_by_reaction = {str(row["reaction"]): row for row in optional_rows}

    pair_rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    for data_tier, penalty_scales in [
        ("core", core_scales),
        ("available_plus_lite", combined_scales),
    ]:
        fit = fit_inverse(
            residual,
            measured,
            LEGACY_PRIMARY_METHOD,
            hyperparameters["lambda_l1"],
            hyperparameters["lambda_l2"],
            {},
            SI_THRESHOLD,
            penalty_scales=penalty_scales,
        )
        extents = np.asarray(fit["extents"], dtype=float)
        selected = support_from_extents(extents, support_threshold)
        selected_classes = class_set(selected, class_map)
        evidence_scores = []
        for reaction, extent in zip(REACTIONS, extents):
            core_row = core_by_reaction[reaction.label]
            optional_row = optional_by_reaction[reaction.label]
            optional_score = float(optional_row["optional_evidence_score"])
            combined_score = (
                float(core_row["hydrosheaf_core_evidence_score"])
                if data_tier == "core"
                else _combine_evidence_scores(
                    core_row["hydrosheaf_core_evidence_score"],
                    optional_score,
                )
            )
            evidence_scores.append(combined_score)
            evidence_rows.append(
                {
                    "dataset": dataset,
                    "edge_id": edge_id,
                    "data_tier": data_tier,
                    "upstream_id": upstream["sample_id"],
                    "downstream_id": downstream["sample_id"],
                    "reaction": reaction.label,
                    "family": reaction.family,
                    "equivalence_class": class_map[reaction.label],
                    "selected": reaction.label in selected,
                    "extent_mmolL": float(extent),
                    "core_evidence_score": core_row[
                        "hydrosheaf_core_evidence_score"
                    ],
                    "optional_evidence_score": optional_score,
                    "combined_evidence_score": combined_score,
                    "core_penalty_scale": core_row["penalty_scale"],
                    "optional_penalty_scale": optional_row[
                        "optional_penalty_scale"
                    ],
                    "combined_penalty_scale": float(
                        penalty_scales[REACTION_LABELS.index(reaction.label)]
                    ),
                    "optional_diagnostic_alignment": optional_row[
                        "optional_diagnostic_alignment"
                    ],
                    "optional_support_alignment": optional_row[
                        "optional_support_alignment"
                    ],
                    "diagnostics_with_nonzero_signature": optional_row[
                        "diagnostics_with_nonzero_signature"
                    ],
                    "panel_ions": ";".join(measured),
                    "optional_diagnostics": ";".join(optional_diagnostics),
                    "evidence_flags": core_row["evidence_flags"],
                }
            )
        metadata = {
            f"meta_{key}": value for key, value in dict(upstream["metadata"]).items()
        }
        pair_rows.append(
            {
                "dataset": dataset,
                "edge_id": edge_id,
                "data_tier": data_tier,
                "upstream_id": upstream["sample_id"],
                "downstream_id": downstream["sample_id"],
                "n_measured_ions": len(measured),
                "measured_ions": ";".join(measured),
                "n_optional_diagnostics": len(optional_diagnostics),
                "optional_diagnostics": ";".join(optional_diagnostics),
                "support_threshold_mmolL": support_threshold,
                "selected_support": ";".join(sorted(selected)),
                "selected_classes": ";".join(sorted(selected_classes)),
                "mean_combined_evidence_score": float(np.mean(evidence_scores)),
                "max_combined_evidence_score": float(np.max(evidence_scores)),
                "solver_rmse_mmolL": fit["rmse"],
                "l1_norm": fit["l1_norm"],
                "thermodynamic_violations": fit["thermodynamic_violations"],
                "converged": fit["converged"],
                "upstream_major_sum_mg_L": upstream["major_sum_mg_L"],
                "downstream_major_sum_mg_L": downstream["major_sum_mg_L"],
                **metadata,
            }
        )
    return pair_rows, evidence_rows


def _legacy_northern_ghana_external_edges() -> list[
    tuple[str, Mapping[str, object], Mapping[str, object]]
]:
    dry = pd.read_excel(LEGACY_NORTHERN_GHANA_WORKBOOK, sheet_name="Dry")
    wet = pd.read_excel(LEGACY_NORTHERN_GHANA_WORKBOOK, sheet_name="Wet")
    dry = dry.copy()
    wet = wet.copy()
    dry["_sample_id"] = dry["Well_ID"].astype(str) + "_dry"
    wet["_sample_id"] = wet["Well_ID"].astype(str) + "_wet"
    column_map = {ion: f"{ion}_mg_L" for ion in ION_ORDER}
    optional_columns = {
        "SiO2": "SiO2_mg_L",
        "Sr": "Sr_mg_L",
        "d18O": "d18O_permil",
        "d2H": "d2H_permil",
    }
    metadata = ["Well_ID", "Region", "District", "Community_Code"]
    wet_records = _external_sample_records(
        wet,
        dataset="NorthernGhana.xlsx",
        id_column="_sample_id",
        x_column="Longitude",
        y_column="Latitude",
        z_column="Elevation_m",
        column_map=column_map,
        optional_columns=optional_columns,
        metadata_columns=metadata,
    )
    dry_records = _external_sample_records(
        dry,
        dataset="NorthernGhana.xlsx",
        id_column="_sample_id",
        x_column="Longitude",
        y_column="Latitude",
        z_column="Elevation_m",
        column_map=column_map,
        optional_columns=optional_columns,
        metadata_columns=metadata,
    )
    wet_by_well = {
        str(record["metadata"]["Well_ID"]): record for record in wet_records
    }
    dry_by_well = {
        str(record["metadata"]["Well_ID"]): record for record in dry_records
    }
    edges = []
    for well_id in sorted(set(wet_by_well) & set(dry_by_well)):
        edges.append((f"{well_id}_wet_to_dry", wet_by_well[well_id], dry_by_well[well_id]))
    return edges


def _spatial_external_edges_from_csv(
    path: Path,
    *,
    dataset: str,
    id_column: str,
    x_column: str,
    y_column: str,
    z_column: str,
    column_map: Mapping[str, str],
    optional_columns: Mapping[str, str],
    metadata_columns: Sequence[str],
) -> list[tuple[str, Mapping[str, object], Mapping[str, object]]]:
    frame = pd.read_csv(path)
    records = _external_sample_records(
        frame,
        dataset=dataset,
        id_column=id_column,
        x_column=x_column,
        y_column=y_column,
        z_column=z_column,
        column_map=column_map,
        optional_columns=optional_columns,
        metadata_columns=metadata_columns,
    )
    edges = []
    for index, (upstream, downstream) in enumerate(
        _nearest_external_edges(records),
        start=1,
    ):
        edge_id = (
            f"{dataset}_{index:03d}_"
            f"{upstream['sample_id']}_to_{downstream['sample_id']}"
        )
        edges.append((edge_id, upstream, downstream))
    return edges


def run_external_field_transfer(
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    edge_specs: list[tuple[str, Mapping[str, object], Mapping[str, object]]] = []
    edge_specs.extend(_legacy_northern_ghana_external_edges())
    edge_specs.extend(
        _spatial_external_edges_from_csv(
            TALENSI_CSV,
            dataset="Talensi",
            id_column="Code",
            x_column="Longitude",
            y_column="Latitude",
            z_column="Elevation",
            column_map={
                "Ca": "Ca",
                "Mg": "Mg",
                "Na": "Na",
                "K": "K",
                "HCO3": "HCO3",
                "Cl": "Cl",
                "SO4": "SO4",
                "NO3": "NO3",
                "Fe": "Fe",
            },
            optional_columns={"d18O": "d18O", "d2H": "d2H"},
            metadata_columns=["Town", "Code"],
        )
    )
    edge_specs.extend(
        _spatial_external_edges_from_csv(
            LOWER_ANAYARI_CSV,
            dataset="LowerAnayari",
            id_column="Sample ID",
            x_column="X coordinate",
            y_column="Y coordinate",
            z_column="Elevation",
            column_map={
                "Ca": "Ca",
                "Mg": "Mg",
                "Na": "Na",
                "K": "K",
                "HCO3": "HCO3",
                "Cl": "Cl",
                "SO4": "SO4",
                "NO3": "NO3",
                "F": "F",
                "Fe": "Fe",
            },
            optional_columns={"d18O": "d18O", "d2H": "d2H"},
            metadata_columns=["Sample ID", "Station"],
        )
    )

    pair_rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    for edge_id, upstream, downstream in edge_specs:
        rows, evidence = _external_field_edge_outputs(
            dataset=str(upstream["dataset"]),
            edge_id=edge_id,
            upstream=upstream,
            downstream=downstream,
            hyperparameters=hyperparameters,
            class_map=class_map,
        )
        pair_rows.extend(rows)
        evidence_rows.extend(evidence)
    pairs = pd.DataFrame(pair_rows)
    evidence = pd.DataFrame(evidence_rows)
    resolution = _evidence_lifted_resolution_frame(
        evidence,
        class_map,
        score_column="combined_evidence_score",
        group_columns=["dataset", "edge_id", "data_tier"],
        evidence_source="external_field_transfer",
    )
    return pairs, evidence, resolution


def write_external_field_transfer_outputs(
    hyperparameters: Mapping[str, float],
    class_map: Mapping[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print("M5: running external field ELRI transfer on NorthernGhana, Talensi, and Lower Anayari...")
    pairs, evidence, resolution = run_external_field_transfer(
        hyperparameters,
        class_map,
    )
    pairs.to_csv(RESULTS_DIR / "external_field_transfer_pairs.csv", index=False)
    evidence.to_csv(RESULTS_DIR / "external_field_reaction_evidence.csv", index=False)
    resolution.to_csv(
        RESULTS_DIR / "external_field_evidence_lifted_resolution.csv",
        index=False,
    )
    return pairs, evidence, resolution


def write_diagnostics_and_classes() -> tuple[dict[str, str], pd.DataFrame]:
    class_map, class_rows = equivalence_classes()
    pd.DataFrame(class_rows).to_csv(
        RESULTS_DIR / "equivalence_classes.csv", index=False
    )
    diagnostic_rows = []
    for panel, ions in PANELS.items():
        diagnostic_rows.append(
            {
                "panel": panel,
                "ions": ";".join(ions),
                **{
                    key: finite_or_blank(value)
                    for key, value in matrix_diagnostics(ions).items()
                },
            }
        )
    diagnostics = pd.DataFrame(diagnostic_rows)
    diagnostics.to_csv(RESULTS_DIR / "identifiability_diagnostics.csv", index=False)
    return class_map, diagnostics


def write_summary(
    scenarios: Sequence[Mapping[str, object]],
    phreeqc_frame: pd.DataFrame,
    fits: pd.DataFrame,
    phreeqc_inverse: pd.DataFrame,
    mrs: pd.DataFrame,
    data_tiers: pd.DataFrame,
    field_pairs: pd.DataFrame,
    hyperparameters: Mapping[str, float],
    phreeqc_elapsed: float,
    phreeqc_inverse_elapsed: float,
    executable: Path,
    database: Path,
) -> None:
    primary = fits[
        (fits["method"] == PRIMARY_METHOD)
        & (fits["panel"] == "full_11")
        & (fits["noise_level"] == 0.03)
    ]
    legacy_primary = fits[
        (fits["method"] == LEGACY_PRIMARY_METHOD)
        & (fits["panel"] == "full_11")
        & (fits["noise_level"] == 0.03)
    ]
    guarded = fits[
        (fits["method"] == GUARDED_METHOD)
        & (fits["panel"] == "full_11")
        & (fits["noise_level"] == 0.03)
    ]
    core_method = fits[
        (fits["method"] == HYDROSHEAF_CORE_METHOD)
        & (fits["panel"] == "full_11")
        & (fits["noise_level"] == 0.03)
    ]
    if primary.empty:
        primary = guarded.copy() if not guarded.empty else legacy_primary.copy()
    low_residual = primary[
        primary["reconstruction_rmse_mmolL"]
        <= primary["reconstruction_rmse_mmolL"].quantile(0.25)
    ]
    low_residual_wrong = float((low_residual["phase_f1"] < 0.80).mean())
    holdout = mrs[mrs["archetype"] == "mixed"]
    inverse_success = (
        float(phreeqc_inverse["phreeqc_inverse_success"].mean())
        if "phreeqc_inverse_success" in phreeqc_inverse
        else float("nan")
    )
    evidence_path = RESULTS_DIR / "hydrosheaf_core_evidence.csv"
    if evidence_path.exists():
        core_evidence = pd.read_csv(evidence_path)
        core_evidence_primary = core_evidence[
            (core_evidence["panel"] == "full_11")
            & (core_evidence["noise_level"] == 0.03)
        ]
    else:
        core_evidence_primary = pd.DataFrame()
    if data_tiers.empty:
        data_tier_summary: dict[str, dict[str, float]] = {}
    else:
        data_tier_summary = {
            str(tier): {
                "mean_phase_f1": float(group["phase_f1"].mean()),
                "mean_class_f1": float(group["class_f1"].mean()),
                "mean_false_discovery_rate": float(
                    group["false_discovery_rate"].mean()
                ),
                "mean_extent_rmse_mmolL": float(group["extent_rmse_mmolL"].mean()),
                "n_optional_diagnostics": float(
                    group["n_optional_diagnostics"].max()
                ),
            }
            for tier, group in data_tiers.groupby("data_tier")
        }
    data_tier_resolution_path = RESULTS_DIR / "data_tier_evidence_lifted_resolution.csv"
    if data_tier_resolution_path.exists():
        data_tier_resolution = pd.read_csv(data_tier_resolution_path)
        data_tier_resolution_summary = {
            str(tier): {
                "mean_evidence_lifted_resolution_index": float(
                    group["evidence_lifted_resolution_index"].mean()
                ),
                "median_evidence_lifted_resolution_index": float(
                    group["evidence_lifted_resolution_index"].median()
                ),
                "conditionally_preferred_or_resolved_fraction": float(
                    group["resolution_status"]
                    .isin(
                        [
                            "conditionally_preferred",
                            "evidence_lifted_resolved",
                        ]
                    )
                    .mean()
                ),
                "n_class_evaluations": int(len(group)),
            }
            for tier, group in data_tier_resolution.groupby("data_tier")
        }
    else:
        data_tier_resolution_summary = {}
    core_resolution_path = RESULTS_DIR / "hydrosheaf_core_evidence_lifted_resolution.csv"
    if core_resolution_path.exists():
        core_resolution = pd.read_csv(core_resolution_path)
        core_resolution_primary = core_resolution[
            (core_resolution["panel"] == "full_11")
            & (core_resolution["noise_level"] == 0.03)
        ]
        core_elri_mean = (
            float(core_resolution_primary["evidence_lifted_resolution_index"].mean())
            if not core_resolution_primary.empty
            else None
        )
    else:
        core_elri_mean = None
    external_resolution_path = (
        RESULTS_DIR / "external_field_evidence_lifted_resolution.csv"
    )
    external_pairs_path = RESULTS_DIR / "external_field_transfer_pairs.csv"
    if external_resolution_path.exists():
        external_resolution = pd.read_csv(external_resolution_path)
        external_field_elri_summary: dict[str, dict[str, dict[str, float]]] = {}
        for keys, group in external_resolution.groupby(["dataset", "data_tier"]):
            dataset, tier = keys
            external_field_elri_summary.setdefault(str(dataset), {})[str(tier)] = {
                "mean_evidence_lifted_resolution_index": float(
                    group["evidence_lifted_resolution_index"].mean()
                ),
                "median_evidence_lifted_resolution_index": float(
                    group["evidence_lifted_resolution_index"].median()
                ),
                "conditionally_preferred_or_resolved_fraction": float(
                    group["resolution_status"]
                    .isin(
                        [
                            "conditionally_preferred",
                            "evidence_lifted_resolved",
                        ]
                    )
                    .mean()
                ),
                "n_class_evaluations": int(len(group)),
            }
    else:
        external_field_elri_summary = {}
    if external_pairs_path.exists():
        external_pairs = pd.read_csv(external_pairs_path)
        external_field_edge_counts = {
            str(dataset): int(group["edge_id"].nunique())
            for dataset, group in external_pairs.groupby("dataset")
        }
    else:
        external_pairs = pd.DataFrame()
        external_field_edge_counts = {}
    inverse_strict_success = (
        float(
            (
                phreeqc_inverse["phreeqc_inverse_success"]
                & phreeqc_inverse["phreeqc_inverse_setup"].eq("strict_5pct")
            ).mean()
        )
        if {"phreeqc_inverse_success", "phreeqc_inverse_setup"}.issubset(
            phreeqc_inverse.columns
        )
        else float("nan")
    )
    inverse_relaxed_fraction = (
        float(phreeqc_inverse["phreeqc_inverse_setup"].eq("relaxed_20pct").mean())
        if "phreeqc_inverse_setup" in phreeqc_inverse
        else float("nan")
    )
    summary = {
        "generated_at": pd.Timestamp.now(tz="UTC").isoformat(),
        "random_seed": RANDOM_SEED,
        "phreeqc_version": "3.7.3-15968",
        "phreeqc_executable": str(executable),
        "phreeqc_database": str(database),
        "phreeqc_runtime_seconds": (
            float(phreeqc_elapsed) if np.isfinite(phreeqc_elapsed) else None
        ),
        "phreeqc_inverse_runtime_seconds": (
            float(phreeqc_inverse_elapsed)
            if np.isfinite(phreeqc_inverse_elapsed)
            else None
        ),
        "n_phreeqc_scenarios": len(scenarios),
        "n_benchmark_fits": len(fits),
        "n_phreeqc_inverse_scenarios": len(phreeqc_inverse),
        "n_field_pairs": len(field_pairs),
        "archetypes": sorted(ARCHETYPES),
        "noise_levels": NOISE_LEVELS,
        "panels": PANELS,
        "selected_hyperparameters": dict(hyperparameters),
        "primary_method": PRIMARY_METHOD,
        "guarded_method": GUARDED_METHOD,
        "hydrosheaf_core_method": HYDROSHEAF_CORE_METHOD,
        "legacy_primary_method": LEGACY_PRIMARY_METHOD,
        "maximum_phreeqc_generation_rmse_mmolL": float(
            phreeqc_frame["generation_rmse_mmolL"].max()
        ),
        "primary_mean_phase_f1": float(primary["phase_f1"].mean()),
        "primary_mean_class_f1": float(primary["class_f1"].mean()),
        "primary_mean_false_discovery_rate": float(
            primary["false_discovery_rate"].mean()
        ),
        "primary_mean_extent_rmse_mmolL": float(
            primary["extent_rmse_mmolL"].mean()
        ),
        "legacy_thermo_mean_phase_f1": float(legacy_primary["phase_f1"].mean()),
        "legacy_thermo_mean_class_f1": float(legacy_primary["class_f1"].mean()),
        "legacy_thermo_mean_false_discovery_rate": float(
            legacy_primary["false_discovery_rate"].mean()
        ),
        "guarded_mean_phase_f1": (
            float(guarded["phase_f1"].mean()) if not guarded.empty else None
        ),
        "guarded_mean_class_f1": (
            float(guarded["class_f1"].mean()) if not guarded.empty else None
        ),
        "guarded_mean_false_discovery_rate": (
            float(guarded["false_discovery_rate"].mean()) if not guarded.empty else None
        ),
        "hydrosheaf_core_mean_phase_f1": (
            float(core_method["phase_f1"].mean()) if not core_method.empty else None
        ),
        "hydrosheaf_core_mean_class_f1": (
            float(core_method["class_f1"].mean()) if not core_method.empty else None
        ),
        "hydrosheaf_core_mean_false_discovery_rate": (
            float(core_method["false_discovery_rate"].mean())
            if not core_method.empty
            else None
        ),
        "hydrosheaf_core_mean_evidence_score": (
            float(core_evidence_primary["hydrosheaf_core_evidence_score"].mean())
            if not core_evidence_primary.empty
            else None
        ),
        "hydrosheaf_core_supported_true_active_fraction": (
            float(
                core_evidence_primary[
                    core_evidence_primary["hydrosheaf_core_evidence_score"] >= 0.60
                ]["true_active"].mean()
            )
            if not core_evidence_primary.empty
            else None
        ),
        "data_tier_summary": data_tier_summary,
        "data_tier_evidence_lifted_resolution_summary": (
            data_tier_resolution_summary
        ),
        "hydrosheaf_core_mean_evidence_lifted_resolution_index": core_elri_mean,
        "external_field_edge_counts": external_field_edge_counts,
        "external_field_evidence_lifted_resolution_summary": (
            external_field_elri_summary
        ),
        "phreeqc_inverse_success_fraction": inverse_success,
        "phreeqc_inverse_strict_success_fraction": inverse_strict_success,
        "phreeqc_inverse_relaxed_fallback_fraction": inverse_relaxed_fraction,
        "phreeqc_inverse_mean_models_found": float(
            phreeqc_inverse["models_found"].mean()
        ) if "models_found" in phreeqc_inverse else None,
        "phreeqc_inverse_mean_minimal_models_found": float(
            phreeqc_inverse["minimal_models_found"].mean()
        ) if "minimal_models_found" in phreeqc_inverse else None,
        "phreeqc_inverse_first_minimal_mean_phase_f1": float(
            phreeqc_inverse["first_minimal_phase_f1"].mean()
        ) if "first_minimal_phase_f1" in phreeqc_inverse else None,
        "phreeqc_inverse_first_minimal_mean_class_f1": float(
            phreeqc_inverse["first_minimal_class_f1"].mean()
        ) if "first_minimal_class_f1" in phreeqc_inverse else None,
        "phreeqc_inverse_minimal_union_mean_class_f1": float(
            phreeqc_inverse["minimal_union_class_f1"].mean()
        ) if "minimal_union_class_f1" in phreeqc_inverse else None,
        "phreeqc_inverse_oracle_minimal_mean_class_f1": float(
            phreeqc_inverse["oracle_minimal_class_f1"].mean()
        ) if "oracle_minimal_class_f1" in phreeqc_inverse else None,
        "low_residual_wrong_phase_fraction": low_residual_wrong,
        "mixed_holdout_mrs_classification_accuracy": float(
            holdout["resolution_class_correct"].mean()
        ),
        "field_resolution_class_counts": (
            field_pairs["resolution_class"].value_counts().to_dict()
            if "resolution_class" in field_pairs
            else {}
        ),
        "field_median_core_evidence_score": (
            float(field_pairs["mean_hydrosheaf_core_evidence_score"].median())
            if "mean_hydrosheaf_core_evidence_score" in field_pairs
            else None
        ),
        "field_median_tds_consistency_score": (
            float(field_pairs["tds_delta_consistency_score"].median())
            if "tds_delta_consistency_score" in field_pairs
            else None
        ),
        "field_pairs_with_plus_lite_tracers": (
            int((field_pairs["plus_lite_tracer_count"] > 0).sum())
            if "plus_lite_tracer_count" in field_pairs
            else 0
        ),
        "claim_guardrail": (
            "Hydrosheaf is evaluated as an identifiability-aware sparse linear "
            "inverse reaction model with PHREEQC-generated truth and "
            "thermodynamic screening, not as a fully coupled nonlinear "
            "PHREEQC inverse solver."
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    (RESULTS_DIR / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    def fmt3(value: object) -> str:
        number = _finite_float(value)
        return f"{number:.3f}" if number is not None else "NA"

    def tier_text(tier: str) -> str:
        item = dict(summary["data_tier_summary"]).get(tier, {})
        return (
            f"{fmt3(item.get('mean_class_f1'))} class F1, "
            f"{fmt3(item.get('mean_false_discovery_rate'))} FDR"
        )

    def tier_elri_text(tier: str) -> str:
        item = dict(
            summary["data_tier_evidence_lifted_resolution_summary"]
        ).get(tier, {})
        return (
            f"{fmt3(item.get('mean_evidence_lifted_resolution_index'))} mean ELRI, "
            f"{fmt3(item.get('conditionally_preferred_or_resolved_fraction'))} "
            "conditionally preferred/resolved"
        )

    def external_elri_text(dataset: str) -> str:
        item = (
            dict(summary["external_field_evidence_lifted_resolution_summary"])
            .get(dataset, {})
            .get("available_plus_lite", {})
        )
        edges = dict(summary["external_field_edge_counts"]).get(dataset, 0)
        return (
            f"{dataset}: {edges} edges, "
            f"{fmt3(item.get('median_evidence_lifted_resolution_index'))} "
            "median ELRI"
        )

    lines = [
        "# M5 Complete Analysis Summary",
        "",
        f"- Live PHREEQC scenarios: {summary['n_phreeqc_scenarios']}.",
        f"- Factorial inverse fits: {summary['n_benchmark_fits']}.",
        f"- Northern Ghana quantitative wet-dry pairs: {summary['n_field_pairs']}.",
        (
            "- Maximum PHREEQC-to-stoichiometric generation RMSE: "
            f"{summary['maximum_phreeqc_generation_rmse_mmolL']:.3e} mmol/L."
        ),
        (
            "- Primary Hydrosheaf Guarded 3% noise full-panel phase F1 / "
            "equivalence-class F1: "
            f"{summary['primary_mean_phase_f1']:.3f} / "
            f"{summary['primary_mean_class_f1']:.3f}."
        ),
        (
            "- Primary Hydrosheaf Guarded false-discovery rate: "
            f"{summary['primary_mean_false_discovery_rate']:.3f}; "
            "Hydrosheaf-Core evidence-gated comparator class F1 / "
            "false-discovery rate: "
            f"{fmt3(summary['hydrosheaf_core_mean_class_f1'])} / "
            f"{fmt3(summary['hydrosheaf_core_mean_false_discovery_rate'])}; "
            "legacy thermodynamic elastic-net false-discovery rate: "
            f"{summary['legacy_thermo_mean_false_discovery_rate']:.3f}."
        ),
        (
            "- Hydrosheaf-Core mean reaction-evidence score: "
            f"{fmt3(summary['hydrosheaf_core_mean_evidence_score'])}; "
            "fraction of high-evidence synthetic reaction rows that were truly "
            "active: "
            f"{fmt3(summary['hydrosheaf_core_supported_true_active_fraction'])}."
        ),
        (
            "- Data-tier experiment at 3% noise: Core "
            f"{tier_text('core')}; Plus-lite {tier_text('plus_lite')}; "
            f"Enhanced {tier_text('enhanced')}."
        ),
        (
            "- Evidence-lifted resolution of ambiguous reaction classes: Core "
            f"{tier_elri_text('core')}; Plus-lite {tier_elri_text('plus_lite')}; "
            f"Enhanced {tier_elri_text('enhanced')}. ELRI quantifies "
            "within-class evidence separation, not new mass-balance uniqueness."
        ),
        (
            "- Conventional PHREEQC inverse-model baseline success fraction "
            "(strict 5% plus relaxed 20% fallback): "
            f"{summary['phreeqc_inverse_success_fraction']:.3f}; strict-only "
            f"success fraction: {summary['phreeqc_inverse_strict_success_fraction']:.3f}; "
            "first-minimal "
            "phase F1 / equivalence-class F1: "
            f"{summary['phreeqc_inverse_first_minimal_mean_phase_f1']:.3f} / "
            f"{summary['phreeqc_inverse_first_minimal_mean_class_f1']:.3f}."
        ),
        (
            "- PHREEQC inverse-model multiplicity: "
            f"{summary['phreeqc_inverse_mean_models_found']:.2f} feasible models "
            "and "
            f"{summary['phreeqc_inverse_mean_minimal_models_found']:.2f} minimal "
            "models per scenario on average."
        ),
        (
            "- Fraction of lowest-residual fits with phase F1 below 0.80: "
            f"{summary['low_residual_wrong_phase_fraction']:.3f}."
        ),
        (
            "- Mixed-archetype held-out MRS classification accuracy: "
            f"{summary['mixed_holdout_mrs_classification_accuracy']:.3f}."
        ),
        (
            "- Ghana median Hydrosheaf-Core evidence score / TDS consistency "
            "score: "
            f"{fmt3(summary['field_median_core_evidence_score'])} / "
            f"{fmt3(summary['field_median_tds_consistency_score'])}; "
            "pairs with optional SiO2/Sr/isotope support: "
            f"{summary['field_pairs_with_plus_lite_tracers']}."
        ),
        (
            "- External field ELRI transfer: "
            f"{external_elri_text('NorthernGhana.xlsx')}; "
            f"{external_elri_text('Talensi')}; "
            f"{external_elri_text('LowerAnayari')}. These are field "
            "plausibility audits, not reaction-truth validation."
        ),
        "",
        "## Claim Guardrail",
        "",
        str(summary["claim_guardrail"]),
        "",
        "The Northern Ghana component is a chemistry-only seasonal transfer "
        "demonstration without independent flow-path, groundwater-age, or "
        "reaction-truth validation.",
        "",
    ]
    with (DOCS_DIR / "m5_results_summary.md").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        handle.write("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--skip-field",
        action="store_true",
        help="Skip the Northern Ghana workbook analysis.",
    )
    parser.add_argument(
        "--reuse-synthetic",
        action="store_true",
        help="Reuse existing synthetic outputs and run only field/summary stages.",
    )
    parser.add_argument(
        "--reuse-phreeqc",
        action="store_true",
        help="Reuse phreeqc_ground_truth.csv and regenerate synthetic model outputs.",
    )
    args = parser.parse_args()
    if args.reuse_synthetic and args.reuse_phreeqc:
        parser.error("--reuse-synthetic and --reuse-phreeqc are mutually exclusive.")
    ensure_directories()
    print("M5: building structural diagnostics and equivalence classes...")
    class_map, _ = write_diagnostics_and_classes()
    if args.reuse_synthetic:
        print("M5: reusing existing synthetic benchmark outputs...")
        scenarios = generate_scenarios()
        phreeqc_frame = pd.read_csv(RESULTS_DIR / "phreeqc_ground_truth.csv")
        hydrate_scenarios_from_phreeqc(scenarios, phreeqc_frame)
        selection = pd.read_csv(RESULTS_DIR / "hyperparameter_selection.csv")
        hyperparameters = _hyperparameters_from_selection(selection)
        fits = pd.read_csv(RESULTS_DIR / "benchmark_fits.csv")
        mrs = add_mechanism_resolution_scores(fits, class_map)
        mrs.to_csv(
            RESULTS_DIR / "mechanism_resolution_scores.csv", index=False
        )
        data_tiers = write_data_tier_outputs(scenarios, hyperparameters, class_map)
        phreeqc_elapsed = float("nan")
        executable, database = find_phreeqc()
        phreeqc_inverse, _, phreeqc_inverse_elapsed = (
            _load_or_run_phreeqc_inverse_baseline(
                scenarios,
                phreeqc_frame,
                class_map,
                executable,
                database,
                reuse_existing=True,
            )
        )
    elif args.reuse_phreeqc:
        print("M5: reusing PHREEQC ground truth and regenerating model outputs...")
        scenarios = generate_scenarios()
        phreeqc_frame = pd.read_csv(RESULTS_DIR / "phreeqc_ground_truth.csv")
        hydrate_scenarios_from_phreeqc(scenarios, phreeqc_frame)
        phreeqc_elapsed = float("nan")
        executable, database = find_phreeqc()
        selection_path = RESULTS_DIR / "hyperparameter_selection.csv"
        if selection_path.exists():
            print("M5: reusing locked sparse-model hyperparameters...")
            selection = pd.read_csv(selection_path)
            hyperparameters = _hyperparameters_from_selection(selection)
        else:
            print("M5: calibrating guarded sparse-model hyperparameters...")
            hyperparameters, selection = calibrate_hyperparameters(scenarios, class_map)
            selection.to_csv(selection_path, index=False)
        print(
            "M5: selected lambda_l1="
            f"{hyperparameters['lambda_l1']:.6g}, "
            f"lambda_l2={hyperparameters['lambda_l2']:.6g}, "
            "support_threshold="
            f"{hyperparameters['support_threshold_mmolL']:.3f} mmol/L."
        )
        fits_path = RESULTS_DIR / "benchmark_fits.csv"
        recovery_path = RESULTS_DIR / "reaction_recovery.csv"
        heldout_path = RESULTS_DIR / "heldout_ion_results.csv"
        if fits_path.exists() and recovery_path.exists() and heldout_path.exists():
            print("M5: refreshing Hydrosheaf-Core method rows only...")
            core_fits, core_recovery, core_heldout, core_evidence = run_factorial_fits(
                scenarios,
                hyperparameters,
                class_map,
                method_names=[HYDROSHEAF_CORE_METHOD],
            )
            fits = _replace_method_rows(fits_path, core_fits, HYDROSHEAF_CORE_METHOD)
            recovery = _replace_method_rows(
                recovery_path, core_recovery, HYDROSHEAF_CORE_METHOD
            )
            heldout = _replace_method_rows(
                heldout_path, core_heldout, HYDROSHEAF_CORE_METHOD
            )
        else:
            print("M5: running guarded factorial model comparisons...")
            fits, recovery, heldout, core_evidence = run_factorial_fits(
                scenarios, hyperparameters, class_map
            )
            fits.to_csv(fits_path, index=False)
            recovery.to_csv(recovery_path, index=False)
            heldout.to_csv(heldout_path, index=False)
        core_evidence.to_csv(
            RESULTS_DIR / "hydrosheaf_core_evidence.csv", index=False
        )
        print("M5: calculating Mechanism Resolution Scores...")
        mrs = add_mechanism_resolution_scores(fits, class_map)
        mrs.to_csv(RESULTS_DIR / "mechanism_resolution_scores.csv", index=False)
        data_tiers = write_data_tier_outputs(scenarios, hyperparameters, class_map)
        phreeqc_inverse, _, phreeqc_inverse_elapsed = (
            _load_or_run_phreeqc_inverse_baseline(
                scenarios,
                phreeqc_frame,
                class_map,
                executable,
                database,
                reuse_existing=True,
            )
        )
        print("M5: reusing or refreshing auxiliary synthetic analyses...")
        regularization_path = RESULTS_DIR / "regularization_paths.csv"
        if not regularization_path.exists():
            regularization = run_regularization_paths(scenarios, hyperparameters)
            regularization.to_csv(regularization_path, index=False)
        next_best_path = RESULTS_DIR / "next_best_measurement.csv"
        if not next_best_path.exists():
            next_best = run_next_best_measurement(
                scenarios, hyperparameters, class_map
            )
            next_best.to_csv(next_best_path, index=False)
        bootstrap_path = RESULTS_DIR / "bootstrap_support_stability.csv"
        if not bootstrap_path.exists():
            bootstrap = run_bootstrap_stability(scenarios, hyperparameters)
            bootstrap.to_csv(bootstrap_path, index=False)
        thermo_sensitivity_path = RESULTS_DIR / "thermodynamic_threshold_sensitivity.csv"
        if not thermo_sensitivity_path.exists():
            thermo_sensitivity = run_thermodynamic_threshold_sensitivity(
                scenarios, hyperparameters, class_map
            )
            thermo_sensitivity.to_csv(thermo_sensitivity_path, index=False)
    else:
        print("M5: generating 240 live-PHREEQC reaction scenarios...")
        scenarios = generate_scenarios()
        phreeqc_frame, phreeqc_elapsed, executable, database = run_live_phreeqc(
            scenarios
        )
        phreeqc_frame.to_csv(RESULTS_DIR / "phreeqc_ground_truth.csv", index=False)
        print("M5: calibrating sparse-model hyperparameters on non-mixed archetypes...")
        hyperparameters, selection = calibrate_hyperparameters(scenarios, class_map)
        selection.to_csv(RESULTS_DIR / "hyperparameter_selection.csv", index=False)
        print(
            "M5: selected lambda_l1="
            f"{hyperparameters['lambda_l1']:.6g}, "
            f"lambda_l2={hyperparameters['lambda_l2']:.6g}, "
            "support_threshold="
            f"{hyperparameters['support_threshold_mmolL']:.3f} mmol/L."
        )
        print("M5: running factorial model, noise, and missing-ion comparisons...")
        fits, recovery, heldout, core_evidence = run_factorial_fits(
            scenarios, hyperparameters, class_map
        )
        fits.to_csv(RESULTS_DIR / "benchmark_fits.csv", index=False)
        recovery.to_csv(RESULTS_DIR / "reaction_recovery.csv", index=False)
        heldout.to_csv(RESULTS_DIR / "heldout_ion_results.csv", index=False)
        core_evidence.to_csv(
            RESULTS_DIR / "hydrosheaf_core_evidence.csv", index=False
        )
        print("M5: calculating Mechanism Resolution Scores...")
        mrs = add_mechanism_resolution_scores(fits, class_map)
        mrs.to_csv(RESULTS_DIR / "mechanism_resolution_scores.csv", index=False)
        data_tiers = write_data_tier_outputs(scenarios, hyperparameters, class_map)
        phreeqc_inverse, _, phreeqc_inverse_elapsed = (
            _load_or_run_phreeqc_inverse_baseline(
                scenarios,
                phreeqc_frame,
                class_map,
                executable,
                database,
                reuse_existing=False,
            )
        )
        print("M5: running regularization, measurement-value, and stability analyses...")
        regularization = run_regularization_paths(scenarios, hyperparameters)
        regularization.to_csv(RESULTS_DIR / "regularization_paths.csv", index=False)
        next_best = run_next_best_measurement(
            scenarios, hyperparameters, class_map
        )
        next_best.to_csv(RESULTS_DIR / "next_best_measurement.csv", index=False)
        bootstrap = run_bootstrap_stability(scenarios, hyperparameters)
        bootstrap.to_csv(
            RESULTS_DIR / "bootstrap_support_stability.csv", index=False
        )
        thermo_sensitivity = run_thermodynamic_threshold_sensitivity(
            scenarios, hyperparameters, class_map
        )
        thermo_sensitivity.to_csv(
            RESULTS_DIR / "thermodynamic_threshold_sensitivity.csv", index=False
        )
    write_core_evidence_lifted_resolution(class_map)
    if args.skip_field:
        field_pairs = pd.DataFrame()
    else:
        print("M5: running the 2025 Northern Ghana chemistry-only demonstration...")
        field_pairs, field_extents, field_heldout, field_classes, field_evidence = (
            run_ghana_field_demonstration(hyperparameters, class_map)
        )
        field_pairs.to_csv(RESULTS_DIR / "ghana_field_pairs.csv", index=False)
        field_extents.to_csv(RESULTS_DIR / "ghana_field_reaction_extents.csv", index=False)
        field_heldout.to_csv(RESULTS_DIR / "ghana_field_heldout_ions.csv", index=False)
        field_classes.to_csv(RESULTS_DIR / "ghana_field_class_support.csv", index=False)
        field_evidence.to_csv(
            RESULTS_DIR / "ghana_field_hydrosheaf_core_evidence.csv", index=False
        )
        write_field_evidence_lifted_resolution(class_map, field_evidence)
        write_external_field_transfer_outputs(hyperparameters, class_map)
    write_summary(
        scenarios,
        phreeqc_frame,
        fits,
        phreeqc_inverse,
        mrs,
        data_tiers,
        field_pairs,
        hyperparameters,
        phreeqc_elapsed,
        phreeqc_inverse_elapsed,
        executable,
        database,
    )
    print("M5 complete analysis finished. Generate tables and figures next.")


if __name__ == "__main__":
    main()
