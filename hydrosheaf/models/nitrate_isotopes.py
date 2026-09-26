"""
Dual Isotope Mixing Model for Nitrate Source Apportionment.
"""

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Default path to endmembers DB
DEFAULT_DB_PATH = (
    Path(__file__).resolve().parent.parent / "databases" / "nitrate_endmembers.json"
)


@dataclass
class IsotopeSample:
    d15N: float
    d18O: float
    B: Optional[float] = None  # Boron concentration (ug/L or mg/L)
    d11B: Optional[float] = None  # delta-11B (permil)
    ln_B: Optional[float] = None  # ln(B in ug/L)
    sample_id: Optional[str] = None

    def __post_init__(self) -> None:
        if self.ln_B is None and self.B is not None:
            val = float(self.B)
            if val > 0:
                self.ln_B = math.log(val)


@dataclass
class SourceIsotopes:
    name: str
    d15N_mean: float
    d15N_std: float
    d18O_mean: float
    d18O_std: float
    d15N_d18O_corr: Optional[float] = None
    d15N_d18O_cov: Optional[float] = None
    ln_B_mean: Optional[float] = None
    ln_B_std: Optional[float] = None
    d11B_mean: Optional[float] = None
    d11B_std: Optional[float] = None
    B_mean: Optional[float] = None
    B_std: Optional[float] = None

    def covariance_d15N_d18O(self) -> float:
        if self.d15N_d18O_cov is not None:
            return float(self.d15N_d18O_cov)
        if self.d15N_d18O_corr is None:
            return 0.0
        corr = float(self.d15N_d18O_corr)
        corr = max(-0.999999, min(0.999999, corr))
        return corr * float(self.d15N_std) * float(self.d18O_std)


DEFAULT_PROCESS_CONSTRAINTS: Dict[str, Any] = {
    "denitrification_enabled": True,
    "denitrification_min_extent": 0.05,
    "denitrification_strength": 1.5,
    "denitrification_slope": 0.65,
    "denitrification_slope_sigma": 0.25,
    "denitrification_perp_sigma": 4.0,
    "nitrification_enabled": True,
    "nitrification_strength": 1.0,
    "nitrification_sigma": 4.0,
    "nitrification_o2_o18": 23.5,
}

DEFAULT_ISOTOPE_QC_CONFIG: Dict[str, float] = {
    "min_top_probability": 0.55,
    "min_top_margin": 0.12,
    "max_normalized_entropy": 0.90,
    "min_tail_probability": 0.02,
    "max_sensitivity_tv": 0.35,
    "sensitivity_delta": 0.5,
    "min_top_ci_gap": 0.0,
}


def _parse_optional_float(value: object) -> Optional[float]:
    if value is None:
        return None
    return float(value)


def _gaussian_score(value: float, mean: float, sigma: float) -> float:
    sigma = max(float(sigma), 1e-12)
    z = (float(value) - float(mean)) / sigma
    return math.exp(-0.5 * z * z)


def _parse_process_config(config: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    merged = dict(DEFAULT_PROCESS_CONSTRAINTS)
    if config:
        merged.update(config)
    return merged


def _parse_qc_config(config: Optional[Dict[str, Any]]) -> Dict[str, float]:
    merged: Dict[str, float] = dict(DEFAULT_ISOTOPE_QC_CONFIG)
    if config:
        for key, value in config.items():
            if key in merged and isinstance(value, (int, float)):
                merged[key] = float(value)
    return merged


def _normalize_probabilities(
    source_probs: Dict[str, float], sources: List[SourceIsotopes]
) -> Dict[str, float]:
    if not sources:
        return {}
    values: Dict[str, float] = {}
    total = 0.0
    for src in sources:
        value = max(float(source_probs.get(src.name, 0.0)), 0.0)
        values[src.name] = value
        total += value
    if total <= 0.0:
        uniform = 1.0 / float(len(sources))
        return {src.name: uniform for src in sources}
    return {name: val / total for name, val in values.items()}


def _summary_moments(
    source_probs: Dict[str, float], sources: List[SourceIsotopes]
) -> Dict[str, float]:
    mean_n = 0.0
    mean_o = 0.0
    var_n_mix = 0.0
    var_o_mix = 0.0
    cov_no_mix = 0.0
    for src in sources:
        p = float(source_probs.get(src.name, 0.0))
        mu_n = float(src.d15N_mean)
        mu_o = float(src.d18O_mean)
        var_n = max(float(src.d15N_std) ** 2, 1e-24)
        var_o = max(float(src.d18O_std) ** 2, 1e-24)
        cov_no = float(src.covariance_d15N_d18O())
        mean_n += p * mu_n
        mean_o += p * mu_o
        var_n_mix += p * p * var_n
        var_o_mix += p * p * var_o
        cov_no_mix += p * p * cov_no

    var_n = max(var_n_mix, 1e-24)
    var_o = max(var_o_mix, 1e-24)
    cov_no = cov_no_mix
    max_abs_cov = math.sqrt(var_n * var_o) * 0.999999
    cov_no = max(-max_abs_cov, min(max_abs_cov, cov_no))
    return {
        "mean_d15N": mean_n,
        "mean_d18O": mean_o,
        "var_d15N": var_n,
        "var_d18O": var_o,
        "cov_d15N_d18O": cov_no,
    }


def compute_isotope_qc_diagnostics(
    sample: IsotopeSample,
    sources: List[SourceIsotopes],
    source_probs: Dict[str, float],
    qc_config: Optional[Dict[str, Any]] = None,
    prior_probs: Optional[Dict[str, float]] = None,
    posterior_samples: Optional[Any] = None,
    ci_lower: Optional[Dict[str, float]] = None,
    ci_upper: Optional[Dict[str, float]] = None,
) -> Tuple[Dict[str, float], List[str]]:
    if not sources:
        return {}, []

    conf = _parse_qc_config(qc_config)
    probs = _normalize_probabilities(source_probs, sources)
    ranked = sorted(probs.items(), key=lambda kv: kv[1], reverse=True)
    top_name, top_prob = ranked[0]
    second_prob = ranked[1][1] if len(ranked) > 1 else 0.0
    top_margin = max(top_prob - second_prob, 0.0)
    n_sources = len(ranked)

    entropy = 0.0
    for _, p in ranked:
        p_safe = max(p, 1e-12)
        entropy -= p_safe * math.log(p_safe)
    normalized_entropy = entropy / math.log(n_sources) if n_sources > 1 else 0.0
    effective_sources = math.exp(entropy)

    moments = _summary_moments(probs, sources)
    dx_n = float(sample.d15N) - moments["mean_d15N"]
    dx_o = float(sample.d18O) - moments["mean_d18O"]
    det = (
        moments["var_d15N"] * moments["var_d18O"]
        - moments["cov_d15N_d18O"] * moments["cov_d15N_d18O"]
    )
    det = max(det, 1e-24)
    mahal = (
        moments["var_d18O"] * dx_n * dx_n
        - 2.0 * moments["cov_d15N_d18O"] * dx_n * dx_o
        + moments["var_d15N"] * dx_o * dx_o
    ) / det
    mahal = max(mahal, 0.0)
    tail_probability = math.exp(-0.5 * mahal)
    rmse = math.sqrt(0.5 * (dx_n * dx_n + dx_o * dx_o))

    delta = max(float(conf["sensitivity_delta"]), 1e-6)
    base_probs = probs
    perturbed = [
        IsotopeSample(float(sample.d15N) + delta, float(sample.d18O)),
        IsotopeSample(float(sample.d15N) - delta, float(sample.d18O)),
        IsotopeSample(float(sample.d15N), float(sample.d18O) + delta),
        IsotopeSample(float(sample.d15N), float(sample.d18O) - delta),
    ]
    tv_values: List[float] = []
    for perturbed_sample in perturbed:
        perturbed_probs = compute_isotope_prob(
            perturbed_sample, sources, prior_probs=prior_probs
        )
        normalized_perturbed = _normalize_probabilities(perturbed_probs, sources)
        tv = 0.5 * sum(
            abs(base_probs[src.name] - normalized_perturbed[src.name]) for src in sources
        )
        tv_values.append(tv)
    sensitivity_tv_mean = float(sum(tv_values) / len(tv_values))
    sensitivity_tv_max = float(max(tv_values))

    top_ci_gap = float("nan")
    if len(ranked) > 1:
        second_name = ranked[1][0]
        if posterior_samples is not None:
            arr = posterior_samples
            try:
                arr_shape = arr.shape
                if len(arr_shape) >= 2 and arr_shape[-1] == len(sources):
                    arr2d = np.asarray(arr).reshape(-1, len(sources))
                    idx_top = [src.name for src in sources].index(top_name)
                    idx_second = [src.name for src in sources].index(second_name)
                    top_low = float(np.percentile(arr2d[:, idx_top], 2.5))
                    second_high = float(np.percentile(arr2d[:, idx_second], 97.5))
                    top_ci_gap = top_low - second_high
            except Exception:
                pass
        if math.isnan(top_ci_gap) and ci_lower and ci_upper:
            if top_name in ci_lower and second_name in ci_upper:
                top_ci_gap = float(ci_lower[top_name]) - float(ci_upper[second_name])

    flags: List[str] = []
    if (
        top_prob < conf["min_top_probability"]
        or top_margin < conf["min_top_margin"]
        or normalized_entropy > conf["max_normalized_entropy"]
    ):
        flags.append("qc_low_identifiability")
    if tail_probability < conf["min_tail_probability"]:
        flags.append("qc_posterior_predictive_mismatch")
    if sensitivity_tv_max > conf["max_sensitivity_tv"]:
        flags.append("qc_high_sensitivity")
    if not math.isnan(top_ci_gap) and top_ci_gap < conf["min_top_ci_gap"]:
        flags.append("qc_posterior_overlap")

    diagnostics = {
        "top_probability": float(top_prob),
        "top_margin": float(top_margin),
        "normalized_entropy": float(normalized_entropy),
        "effective_sources": float(effective_sources),
        "predicted_d15N": float(moments["mean_d15N"]),
        "predicted_d18O": float(moments["mean_d18O"]),
        "posterior_predictive_rmse": float(rmse),
        "posterior_predictive_tail_probability": float(tail_probability),
        "sensitivity_tv_mean": float(sensitivity_tv_mean),
        "sensitivity_tv_max": float(sensitivity_tv_max),
    }
    if not math.isnan(top_ci_gap):
        diagnostics["top_ci_gap"] = float(top_ci_gap)
    if flags:
        diagnostics["n_qc_flags"] = float(len(flags))
    return diagnostics, flags


def compute_process_prior_probs(
    sample: IsotopeSample,
    sources: List[SourceIsotopes],
    base_prior_probs: Optional[Dict[str, float]] = None,
    water_d18O: Optional[float] = None,
    denitrification_extent: Optional[float] = None,
    process_config: Optional[Dict[str, Any]] = None,
) -> Tuple[Dict[str, float], List[str], Dict[str, float]]:
    if not sources:
        return {}, [], {}

    conf = _parse_process_config(process_config)
    n_sources = len(sources)
    priors: Dict[str, float] = {}
    weights: Dict[str, float] = {src.name: 1.0 for src in sources}
    flags: List[str] = []
    diagnostics: Dict[str, float] = {}

    for src in sources:
        if base_prior_probs and src.name in base_prior_probs:
            priors[src.name] = max(float(base_prior_probs[src.name]), 1e-12)
        else:
            priors[src.name] = 1.0 / n_sources

    denit_on = bool(conf.get("denitrification_enabled", True))
    denit_extent = (
        float(denitrification_extent)
        if denitrification_extent is not None
        else None
    )
    if (
        denit_on
        and denit_extent is not None
        and denit_extent >= float(conf.get("denitrification_min_extent", 0.05))
    ):
        slope = float(conf.get("denitrification_slope", 0.65))
        slope_sigma = float(conf.get("denitrification_slope_sigma", 0.25))
        perp_sigma = float(conf.get("denitrification_perp_sigma", 4.0))
        strength = float(conf.get("denitrification_strength", 1.5))
        denit_scale = max(0.0, denit_extent)
        denit_scores: List[float] = []

        for src in sources:
            dx_n = float(sample.d15N) - float(src.d15N_mean)
            dx_o = float(sample.d18O) - float(src.d18O_mean)
            if dx_n <= 0.0:
                continue
            obs_slope = dx_o / max(dx_n, 1e-12)
            slope_score = _gaussian_score(obs_slope, slope, slope_sigma)
            perp = abs(dx_o - slope * dx_n) / math.sqrt(1.0 + slope * slope)
            perp_score = _gaussian_score(perp, 0.0, perp_sigma)
            score = slope_score * perp_score
            denit_scores.append(score)
            weights[src.name] *= 1.0 + strength * denit_scale * score

        if denit_scores and max(denit_scores) > 0.3:
            flags.append("denitrification_trajectory_match")
        diagnostics["denitrification_extent"] = denit_extent
        diagnostics["denitrification_best_score"] = (
            max(denit_scores) if denit_scores else 0.0
        )

    nitrif_on = bool(conf.get("nitrification_enabled", True))
    water_o18 = float(water_d18O) if water_d18O is not None else None
    if nitrif_on and water_o18 is not None:
        o2_o18 = float(conf.get("nitrification_o2_o18", 23.5))
        nitrif_sigma = float(conf.get("nitrification_sigma", 4.0))
        nitrif_strength = float(conf.get("nitrification_strength", 1.0))
        expected_nitrif = (2.0 / 3.0) * water_o18 + (1.0 / 3.0) * o2_o18
        sample_score = _gaussian_score(float(sample.d18O), expected_nitrif, nitrif_sigma)

        for src in sources:
            source_sigma = max(float(src.d18O_std), 1e-12)
            source_alignment = _gaussian_score(
                float(src.d18O_mean), expected_nitrif, source_sigma
            )
            weights[src.name] *= 1.0 + nitrif_strength * sample_score * source_alignment

        if sample_score > 0.3:
            flags.append("nitrification_pathway_match")
        diagnostics["nitrification_expected_d18O"] = expected_nitrif
        diagnostics["nitrification_sample_score"] = sample_score

    unnormalized: Dict[str, float] = {}
    total = 0.0
    for src in sources:
        value = max(priors[src.name] * weights[src.name], 1e-18)
        unnormalized[src.name] = value
        total += value

    if total <= 0.0:
        uniform = 1.0 / n_sources
        return ({src.name: uniform for src in sources}, flags, diagnostics)

    adjusted = {name: val / total for name, val in unnormalized.items()}
    if flags:
        flags.append("process_constraints_applied")
    return adjusted, flags, diagnostics


@lru_cache(maxsize=8)
def _load_isotope_endmembers_cached(
    path_str: str,
    include_subsources: bool = False,
) -> List[SourceIsotopes]:
    """Load endmember definitions from JSON."""
    path = Path(path_str)
    if not path.exists():
        return []

    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    raw_sources = data.get("sources", {})
    has_subsources = "Septic_Sewage" in raw_sources or "Animal_Manure" in raw_sources
    has_legacy_manure = "Manure" in raw_sources

    sources = []
    for name, params in raw_sources.items():
        if has_subsources and has_legacy_manure:
            if include_subsources:
                if name == "Manure":
                    continue
            else:
                if name in ("Septic_Sewage", "Animal_Manure"):
                    continue
        ln_b_params = params.get("ln_B")
        d11b_params = params.get("d11B")
        b_params = params.get("B")

        ln_b_mean = _parse_optional_float(ln_b_params.get("mean")) if isinstance(ln_b_params, dict) else None
        ln_b_std = _parse_optional_float(ln_b_params.get("std")) if isinstance(ln_b_params, dict) else None

        b_mean = _parse_optional_float(b_params.get("mean")) if isinstance(b_params, dict) else None
        b_std = _parse_optional_float(b_params.get("std")) if isinstance(b_params, dict) else None
        if ln_b_mean is None and b_mean is not None and b_mean > 0:
            ln_b_mean = math.log(b_mean)
            ln_b_std = 0.5 if b_std is None else max(0.1, b_std / b_mean)

        d11b_mean = _parse_optional_float(d11b_params.get("mean")) if isinstance(d11b_params, dict) else None
        d11b_std = _parse_optional_float(d11b_params.get("std")) if isinstance(d11b_params, dict) else None

        sources.append(
            SourceIsotopes(
                name=name,
                d15N_mean=params["d15N"]["mean"],
                d15N_std=params["d15N"]["std"],
                d18O_mean=params["d18O"]["mean"],
                d18O_std=params["d18O"]["std"],
                d15N_d18O_corr=_parse_optional_float(
                    params.get("correlation", {}).get("d15N_d18O")
                    if isinstance(params.get("correlation"), dict)
                    else params.get("d15N_d18O_corr")
                ),
                d15N_d18O_cov=_parse_optional_float(
                    params.get("covariance", {}).get("d15N_d18O")
                    if isinstance(params.get("covariance"), dict)
                    else params.get("d15N_d18O_cov")
                ),
                ln_B_mean=ln_b_mean,
                ln_B_std=ln_b_std,
                d11B_mean=d11b_mean,
                d11B_std=d11b_std,
                B_mean=b_mean,
                B_std=b_std,
            )
        )
    return sources


def load_isotope_endmembers(
    path: Path = DEFAULT_DB_PATH,
    include_subsources: bool = False,
) -> List[SourceIsotopes]:
    # Return a shallow copy to keep call-site behavior mutable-safe.
    return list(_load_isotope_endmembers_cached(str(path), include_subsources=include_subsources))


def load_endmember_database(path: Path = DEFAULT_DB_PATH) -> Dict[str, SourceIsotopes]:
    """Return dictionary of all endmembers including specialized trace subsources."""
    all_sources = list(_load_isotope_endmembers_cached(str(path), include_subsources=True))
    db = {s.name: s for s in all_sources}
    legacy = list(_load_isotope_endmembers_cached(str(path), include_subsources=False))
    for s in legacy:
        if s.name not in db:
            db[s.name] = s
    if "Manure" not in db and "Animal_Manure" in db:
        db["Manure"] = db["Animal_Manure"]
    return db


def compute_isotope_prob(
    sample: IsotopeSample,
    sources: List[SourceIsotopes],
    prior_probs: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    """
    Compute posterior probability of each source using a multivariate normal likelihood.

    P(Source|Sample) is proportional to P(Sample|Source) * P(Source).

    When Boron (ln_B) and/or delta-11B are provided in the sample and present in sources,
    evaluates a 3D or 4D multivariate likelihood:
        [d15N, d18O, ln_B, d11B]^T ~ N(mu_mix, Sigma_mix)
    which disentangles septic sewage from animal manure. Otherwise evaluates exact 2D likelihood.
    """
    if not sources:
        return {}

    posteriors: Dict[str, float] = {}
    total_likelihood = 0.0
    likelihoods: Dict[str, float] = {}
    x_n = float(sample.d15N)
    x_o = float(sample.d18O)

    # Check for active Boron tracers
    has_b = sample.ln_B is not None and all(s.ln_B_mean is not None for s in sources)
    has_d11b = sample.d11B is not None and all(s.d11B_mean is not None for s in sources)
    x_b = float(sample.ln_B) if has_b and sample.ln_B is not None else 0.0
    x_d11b = float(sample.d11B) if has_d11b and sample.d11B is not None else 0.0

    k_dim = 2 + (1 if has_b else 0) + (1 if has_d11b else 0)
    norm_const = (2.0 * math.pi) ** (k_dim / 2.0)

    for src in sources:
        var_n = max(float(src.d15N_std) ** 2, 1e-24)
        var_o = max(float(src.d18O_std) ** 2, 1e-24)
        cov_no = float(src.covariance_d15N_d18O())
        max_abs_cov = math.sqrt(var_n * var_o) * 0.999999
        cov_no = max(-max_abs_cov, min(max_abs_cov, cov_no))

        det_2d = var_n * var_o - cov_no * cov_no
        if det_2d <= 1e-30:
            det_2d = 1e-30

        dx_n = x_n - float(src.d15N_mean)
        dx_o = x_o - float(src.d18O_mean)
        inv00 = var_o / det_2d
        inv11 = var_n / det_2d
        inv01 = -cov_no / det_2d
        quad = inv00 * dx_n * dx_n + 2.0 * inv01 * dx_n * dx_o + inv11 * dx_o * dx_o

        det_total = det_2d

        if has_b and src.ln_B_mean is not None:
            std_b = max(float(src.ln_B_std if src.ln_B_std is not None else 0.5), 1e-6)
            var_b = std_b ** 2
            dx_b = x_b - float(src.ln_B_mean)
            quad += (dx_b * dx_b) / var_b
            det_total *= var_b

        if has_d11b and src.d11B_mean is not None:
            std_d11b = max(float(src.d11B_std if src.d11B_std is not None else 3.0), 1e-6)
            var_d11b = std_d11b ** 2
            dx_d11b = x_d11b - float(src.d11B_mean)
            quad += (dx_d11b * dx_d11b) / var_d11b
            det_total *= var_d11b

        exponent = -0.5 * quad
        norm = 1.0 / (norm_const * math.sqrt(det_total))
        lik = norm * math.exp(exponent)

        # Apply prior
        prior = 1.0 / len(sources)
        if prior_probs and src.name in prior_probs:
            prior = prior_probs[src.name]

        unnormalized = lik * prior
        likelihoods[src.name] = unnormalized
        total_likelihood += unnormalized

    # 2. Normalize
    if total_likelihood <= 0:
        # Fallback to uniform if sample is extremely far from all sources
        uniform_p = 1.0 / len(sources)
        return {src.name: uniform_p for src in sources}

    for name, val in likelihoods.items():
        posteriors[name] = val / total_likelihood

    return posteriors
