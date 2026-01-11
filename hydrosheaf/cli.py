"""Command-line interface for hydrosheaf."""

import argparse
import csv
from typing import Dict, Iterable, List, Mapping, Tuple

from .config import Config, DEFAULT_ION_ORDER
from .data.endmembers import load_endmembers_json
from .data.minerals import MINERAL_LIBRARY
from .data.schema import build_endmember_vectors, normalize_samples
from .data.units import mgL_to_mmolL, meqL_to_mmolL
from .inference.network_fit import fit_network, infer_edges
from .models.ec_tds import calibrate_ec_tds
from .isotopes import fit_lmwl
from .outputs.export import export_edge_results_csv, export_edge_results_json
from .outputs.interpret import print_interpretation_report
from .outputs.temporal import temporal_edge_results_to_rows
from .graph.types import Edge
import json


def _parse_weights(value: str) -> List[float]:
    weights = [float(item) for item in value.split(",")]
    if len(weights) != 10:
        raise ValueError("weights must have 10 comma-separated values.")
    return weights


def _parse_ion_order(value: str) -> List[str]:
    order = [item.strip() for item in value.split(",") if item.strip()]
    if len(order) != 10:
        raise ValueError("ion-order must have 10 comma-separated values.")
    return order


def _parse_endmembers(values: Iterable[str]) -> Dict[str, List[float]]:
    endmembers: Dict[str, List[float]] = {}
    for entry in values:
        if ":" not in entry:
            raise ValueError("endmember format must be name:v1,v2,...")
        name, raw = entry.split(":", 1)
        vector = [float(item) for item in raw.split(",")]
        if len(vector) != 10:
            raise ValueError(f"endmember '{name}' must have 10 values.")
        endmembers[name] = vector
    return endmembers


def _parse_layer_minerals(values: Iterable[str]) -> Dict[int, List[str]]:
    mapping: Dict[int, List[str]] = {}
    for entry in values:
        if ":" not in entry:
            raise ValueError("--layer-minerals entries must be like '1:calcite,gypsum'")
        layer_raw, minerals_raw = entry.split(":", 1)
        try:
            layer_idx = int(layer_raw.strip())
        except ValueError as exc:
            raise ValueError(f"Invalid layer index in --layer-minerals: {layer_raw}") from exc
        minerals = [m.strip() for m in minerals_raw.split(",") if m.strip()]
        if not minerals:
            continue
        mapping[layer_idx] = minerals
    return mapping


def _read_samples(path: str) -> List[Dict[str, object]]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows: List[Dict[str, object]] = []
        for row in reader:
            parsed: Dict[str, object] = {}
            for key, value in row.items():
                if value is None:
                    continue
                parsed[key] = value if value != "" else None
            rows.append(parsed)
    return rows


def _read_edges(path: str) -> List[Edge]:
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        edges: List[Edge] = []
        for row in reader:
            if "u" not in row or "v" not in row:
                raise ValueError("edges CSV must include 'u' and 'v' columns.")
            edge_id = row.get("edge_id") or f"{row['u']}->{row['v']}"
            attrs: Dict[str, object] = {}
            for key, value in row.items():
                if key in {"edge_id", "u", "v"}:
                    continue
                if value in (None, ""):
                    continue
                try:
                    attrs[key] = float(value)
                except (TypeError, ValueError):
                    attrs[key] = value
            edges.append(Edge(edge_id=edge_id, u=row["u"], v=row["v"], attrs=attrs))
    return edges


def _required_columns(ion_order: List[str]) -> List[str]:
    return ["sample_id", "site_id", *ion_order, "EC", "TDS", "pH"]


def _validate_required(samples: List[Mapping[str, object]], ion_order: List[str]) -> None:
    required = _required_columns(ion_order)
    for entry in samples:
        missing = [col for col in required if col not in entry]
        if missing:
            raise ValueError(f"Sample {entry.get('sample_id', '?')} missing columns: {missing}")


def _convert_samples(
    samples: List[Dict[str, object]],
    ion_order: List[str],
    unit: str,
) -> List[Dict[str, object]]:
    if unit == "mmol/L":
        return samples
    converted: List[Dict[str, object]] = []
    for entry in samples:
        updated = dict(entry)
        for ion in ion_order:
            if ion not in updated:
                continue
            value = updated[ion]
            if value is None:
                continue
            value = float(value)
            if unit == "mg/L":
                updated[ion] = mgL_to_mmolL(value, ion)
            elif unit == "meq/L":
                updated[ion] = meqL_to_mmolL(value, ion)
        converted.append(updated)
    return converted


def _convert_endmembers(
    endmembers: Dict[str, List[float]],
    ion_order: List[str],
    unit: str,
) -> Dict[str, List[float]]:
    if unit == "mmol/L":
        return endmembers
    converted: Dict[str, List[float]] = {}
    for name, vector in endmembers.items():
        new_vector: List[float] = []
        for value, ion in zip(vector, ion_order):
            if unit == "mg/L":
                new_vector.append(mgL_to_mmolL(value, ion))
            elif unit == "meq/L":
                new_vector.append(meqL_to_mmolL(value, ion))
        converted[name] = new_vector
    return converted


def main() -> None:
    parser = argparse.ArgumentParser(description="Run hydrosheaf edge inference.")
    parser.add_argument("--samples", help="Path to samples CSV.")
    parser.add_argument("--edges", help="Path to edges CSV.")
    parser.add_argument("--physics-priors", type=str, help="Path to physics priors CSV/JSON (u,v,p_uv,tt_mean_days,...).")
    parser.add_argument(
        "--physics-priors-mode",
        choices=["override", "merge", "only"],
        default="override",
        help="How to apply physics priors to edges (override attrs, merge edges, or use only physics edges).",
    )
    parser.add_argument("--modpath-endpoints", type=str, help="Path to MODPATH endpoints file (requires flopy).")
    parser.add_argument("--modpath-max-match-distance", type=float, default=500.0, help="Max distance for matching endpoints to node coords (default: 500).")
    parser.add_argument("--modpath-time-unit-days", type=float, default=1.0, help="Convert MODPATH time units to days (default: 1).")
    parser.add_argument("--output", help="Output file path.")
    parser.add_argument("--format", choices=["json", "csv"], default="json")
    parser.add_argument("--interpret", action="store_true", help="Print a geochemical interpretation report after inference.")
    parser.add_argument("--report-only", type=str, help="Print interpretation report from an existing JSON results file and exit.")
    parser.add_argument("--lambda-sparse", type=float, default=0.0)
    parser.add_argument("--lambda-l1", type=float, default=0.0)
    parser.add_argument("--allow-signed", action="store_true")
    parser.add_argument("--reaction-max-iter", type=int, default=300)
    parser.add_argument("--reaction-tol", type=float, default=1e-6)
    parser.add_argument("--charge-balance-limit", type=float, default=0.1)
    parser.add_argument("--ec-tds-penalty-limit", type=float, default=0.0)
    parser.add_argument("--ec-tds-penalty-enabled", action="store_true")
    parser.add_argument(
        "--missing-policy",
        choices=["skip", "impute_zero"],
        default="skip",
        help="How to handle missing ions.",
    )
    parser.add_argument(
        "--detection-policy",
        choices=["half", "zero", "value", "drop"],
        default="half",
        help="How to handle detection limits like '<0.01'.",
    )
    parser.add_argument("--eta-ec", type=float, default=0.0)
    parser.add_argument("--eta-tds", type=float, default=0.0)
    parser.add_argument("--weights", type=str)
    parser.add_argument(
        "--signed-reaction",
        action="append",
        default=[],
        help="Allow signed extent for a reaction label (repeatable).",
    )
    parser.add_argument(
        "--ion-order",
        type=str,
        help="Comma-separated ion order (10 entries).",
    )
    parser.add_argument(
        "--unit",
        choices=["mmol/L", "mg/L", "meq/L"],
        default="mmol/L",
        help="Input sample units for ions.",
    )
    parser.add_argument(
        "--infer-edges",
        action="store_true",
        help="Infer edges from coordinates when edges CSV is not supplied.",
    )
    parser.add_argument(
        "--infer-edges-method",
        choices=["probabilistic", "probabilistic_map", "simple"],
        default="probabilistic",
        help="Edge inference strategy.",
    )
    parser.add_argument("--max-neighbors", type=int, default=1)
    parser.add_argument("--allow-uphill", action="store_true")
    parser.add_argument("--head-key", type=str, default="hydraulic_head")
    parser.add_argument("--elevation-key", type=str, default="elevation")
    parser.add_argument("--edge-p-min", type=float, default=0.75)
    parser.add_argument("--edge-radius-km", type=float, default=5.0)
    parser.add_argument("--edge-max-neighbors", type=int, default=3)
    parser.add_argument("--edge-sigma-meas", type=float, default=0.5)
    parser.add_argument("--edge-sigma-dtw", type=float, default=1.0)
    parser.add_argument("--edge-sigma-elev", type=float, default=1.0)
    parser.add_argument("--edge-sigma-topo", type=float, default=10.0)
    parser.add_argument(
        "--edge-head-inference",
        choices=["heuristic", "bayesian", "bayesian_mcmc"],
        default="heuristic",
        help="Head estimation strategy for probabilistic edge inference.",
    )
    parser.add_argument("--edge-dtw-prior-mu", type=float, default=5.0, help="Prior mean for depth-to-water (m).")
    parser.add_argument("--edge-dtw-prior-sigma", type=float, default=5.0, help="Prior std for depth-to-water (m).")
    parser.add_argument("--edge-head-prior-mu", type=float, default=0.0, help="Fallback prior mean head (m) if elevation missing.")
    parser.add_argument("--edge-head-prior-sigma", type=float, default=1000.0, help="Fallback prior std head (m) if elevation missing.")
    parser.add_argument("--edge-topo-sigma-depth", type=float, default=5.0, help="3D topographic prior depth uncertainty (m).")
    parser.add_argument("--edge-topo-reject-p", type=float, default=0.1, help="Reject if P(flow downhill) below this threshold.")
    parser.add_argument("--edge-map-weight", type=float, default=1.0, help="MAP prior weight on edge probability (requires probabilistic_map).")
    parser.add_argument("--edge-map-candidate-multiplier", type=int, default=5, help="Candidate expansion factor before MAP selection.")
    parser.add_argument("--edge-map-p-min", type=float, default=0.1, help="Minimum candidate edge probability for MAP selection.")
    parser.add_argument("--edge-gradient-min", type=float, default=1e-4)
    parser.add_argument("--edge-head-key", type=str, default="head_meas")
    parser.add_argument("--edge-dtw-key", type=str, default="dtw")
    parser.add_argument("--edge-elevation-key", type=str, default="elevation")
    parser.add_argument("--edge-aquifer-key", type=str, default="aquifer_unit")
    parser.add_argument("--edge-screen-depth-key", type=str, default="screen_depth")
    parser.add_argument("--edge-well-depth-key", type=str, default="well_depth")
    parser.add_argument("--edge-depth-mismatch", type=float, default=20.0)
    parser.add_argument("--isotope-enabled", action="store_true")
    parser.add_argument("--isotope-weight", type=float, default=1.0)
    parser.add_argument("--isotope-d-excess-weight", type=float, default=0.0)
    parser.add_argument("--isotope-d18o-key", type=str, default="18O")
    parser.add_argument("--isotope-d2h-key", default="2H", help="Column name for delta-2H.")
    parser.add_argument("--auto-lmwl", action="store_true", help="Automatically calibrate LMWL from samples.")
    parser.add_argument("--iso-consistency-weight", type=float, default=0.0, help="Weight for isotope-chemistry consistency penalty.")
    
    # Mineral Library Options
    parser.add_argument("--lmwl-a", type=float)
    parser.add_argument("--lmwl-b", type=float)
    parser.add_argument("--fit-lmwl", action="store_true")
    parser.add_argument(
        "--endmember",
        action="append",
        default=[],
        help="Mixing endmember name:v1,v2,... (10 values).",
    )
    parser.add_argument(
        "--endmember-id",
        action="append",
        default=[],
        help="Mixing endmember site_id from samples (repeatable).",
    )
    parser.add_argument(
        "--endmembers-json",
        type=str,
        help="Path to endmembers.json (see proposal).",
    )
    parser.add_argument(
        "--calibrate-ec-tds",
        action="store_true",
        help="Fit EC/TDS linear models from samples.",
    )
    parser.add_argument("--phreeqc-enabled", action="store_true")
    parser.add_argument("--phreeqc-mode", choices=["phreeqpython", "subprocess"], default="phreeqpython")
    parser.add_argument("--phreeqc-database", type=str, default=Config().phreeqc_database)
    parser.add_argument("--phreeqc-executable", type=str, default="")
    parser.add_argument("--temp-default-c", type=float, default=25.0)
    parser.add_argument("--si-threshold", type=float, default=0.2)
    parser.add_argument("--constraints-hard", action="store_true")
    
    # Nitrate Source Discrim
    parser.add_argument("--nitrate-source-enabled", action="store_true", help="Enable nitrate source discrimination (manure vs fertilizer).")
    parser.add_argument("--nitrate-source-min-conc", type=float, help="Minimum nitrate concentration (mg/L) to attempt discrimination (default: 10).")

    # Uncertainty quantification arguments
    parser.add_argument("--uncertainty", choices=["none", "bootstrap", "bayesian", "monte_carlo"], default="none", help="Uncertainty quantification method")
    parser.add_argument("--bootstrap-n", type=int, default=1000, help="Bootstrap resamples (default: 1000)")
    parser.add_argument("--bootstrap-ci", choices=["percentile", "bca"], default="percentile", help="Bootstrap CI method")
    parser.add_argument("--bayesian-samples", type=int, default=5000, help="MCMC samples per chain (default: 5000)")
    parser.add_argument("--bayesian-chains", type=int, default=4, help="MCMC chains (default: 4)")
    parser.add_argument("--bayesian-accept", type=float, default=0.95, help="Target acceptance rate (default: 0.95)")
    parser.add_argument("--monte-carlo-n", type=int, default=1000, help="Monte Carlo samples (default: 1000)")
    parser.add_argument("--input-uncertainty", type=float, default=5.0, help="Input uncertainty percentage (default: 5.0)")
    parser.add_argument("--uncertainty-seed", type=int, help="Random seed for uncertainty quantification")

    # Temporal dynamics arguments
    parser.add_argument("--temporal-data", type=str, help="Path to time-series CSV")
    parser.add_argument("--temporal-enabled", action="store_true", help="Enable temporal dynamics analysis")
    parser.add_argument("--temporal-output", type=str, help="Optional output path for temporal results JSON.")
    parser.add_argument("--temporal-window", type=int, default=365, help="Analysis window in days (default: 365)")
    parser.add_argument("--temporal-min-samples", type=int, default=3, help="Minimum samples per node (default: 3)")
    parser.add_argument("--temporal-interp", choices=["linear", "spline", "nearest"], default="linear", help="Interpolation method")
    parser.add_argument("--temporal-frequency", type=int, default=30, help="Interpolation grid spacing in days (default: 30)")
    parser.add_argument(
        "--residence-method",
        choices=["gradient", "cross_correlation", "bayesian_lag", "ttd", "tracer_decay"],
        help="Residence time estimation method",
    )
    parser.add_argument(
        "--residence-tracer",
        type=str,
        default="Cl",
        help="Tracer for temporal residence time via cross-correlation (e.g., Cl, 18O, 2H, or 'auto' / comma-list).",
    )
    parser.add_argument("--tau-agreement-tolerance", type=float, default=0.4, help="Relative tolerance for multi-tracer tau disagreement (default: 0.4).")
    parser.add_argument("--tau-min-peak-corr", type=float, default=0.2, help="Minimum correlation peak to accept a tracer (default: 0.2).")
    parser.add_argument("--tau-max-relative-uncertainty", type=float, default=1.5, help="Reject tracer if uncertainty/tau exceeds this (default: 1.5).")
    parser.add_argument("--tau-max-uncertainty-days", type=float, default=180.0, help="Reject tracer if abs uncertainty exceeds this (default: 180).")
    parser.add_argument("--tau-physics-blend-threshold", type=float, default=0.5, help="Blend tau with physics prior if relative diff exceeds this (default: 0.5).")
    parser.add_argument("--ttd-grid-dt-days", type=float, default=30.0, help="TTD convolution grid step in days (default: 30).")
    parser.add_argument("--ttd-max-lag-days", type=float, default=365.0, help="TTD maximum lag support in days (default: 365).")
    parser.add_argument("--ttd-smoothness-lambda", type=float, default=0.0, help="TTD smoothness penalty strength (default: 0).")
    parser.add_argument("--ttd-min-r2", type=float, default=0.2, help="Minimum R^2 to accept a tracer TTD fit (default: 0.2).")
    parser.add_argument("--ttd-attenuation-k-max", type=float, default=0.02, help="Max k (1/day) for attenuation grid search (default: 0.02).")
    parser.add_argument("--ttd-attenuation-k-steps", type=int, default=6, help="Number of attenuation k grid points (default: 6).")
    parser.add_argument("--bayes-lag-grid-dt-days", type=float, default=5.0, help="Bayesian lag tau grid step in days (default: 5).")
    parser.add_argument("--bayes-lag-max-lag-days", type=float, default=365.0, help="Bayesian lag max lag support in days (default: 365).")
    parser.add_argument("--bayes-lag-prior-sigma-multiplier", type=float, default=1.0, help="Multiplies the physics sigma for tau prior width (default: 1).")
    parser.add_argument("--bayes-lag-min-pairs", type=int, default=5, help="Minimum overlap pairs for Bayesian lag (default: 5).")
    parser.add_argument("--residence-k", type=float, default=1.0, help="Hydraulic conductivity m/day for gradient method (default: 1.0)")
    parser.add_argument("--residence-porosity", type=float, default=0.2, help="Effective porosity for gradient method (default: 0.2)")
    parser.add_argument("--residence-coupling", action="store_true", help="Couple static edge chemistry sparsity to residence time (short tau => stronger sparsity).")
    parser.add_argument("--residence-ref-days", type=float, default=30.0, help="Reference residence time (days) for residence-coupling scaling.")
    parser.add_argument(
        "--layer-minerals",
        action="append",
        default=[],
        help="Layer-specific minerals as 'layer:calcite,gypsum,...' (repeatable). Uses downstream node layer_key.",
    )

    # Vadose-zone physics (optional, open 1D Richards column) -> physics priors
    parser.add_argument("--vadose-enabled", action="store_true", help="Enable vadose 1D Richards column to generate travel-time priors.")
    parser.add_argument("--vadose-forcing", type=str, help="Vadose forcing CSV (timestamp,P_mm_day,ET0_mm_day[,I_mm_day]).")
    parser.add_argument("--vadose-profile", type=str, help="Vadose profile JSON/YAML (layers, depth_m, optional root_depth_m).")
    parser.add_argument("--vadose-links", type=str, help="Vadose links CSV mapping u->v edges with u_depth_m (lysimeter depth).")
    parser.add_argument("--vadose-water-table", type=str, help="Optional water table CSV (timestamp,water_table_depth_m) for lower boundary.")
    parser.add_argument("--vadose-lower-boundary", choices=["free_drainage", "constant_head_from_wt"], default="free_drainage")
    parser.add_argument("--vadose-dz", type=float, default=0.1, help="Vadose grid spacing dz (m).")
    parser.add_argument("--vadose-kc", type=float, default=1.0, help="Crop coefficient Kc for ET0->ETp.")
    parser.add_argument("--vadose-evap-frac", type=float, default=0.3, help="Fraction of ETp treated as surface evaporation (rest is transpiration sink).")
    parser.add_argument("--vadose-priors-output", type=str, help="Optional path to write vadose-derived physics priors CSV.")
    parser.add_argument("--vadose-details-output", type=str, help="Optional path to write vadose priors details JSON (TTD kernels, tau series).")
    parser.add_argument("--vadose-calibrate", action="store_true", help="Calibrate vadose Ks multiplier and Kc against observed theta.")
    parser.add_argument("--vadose-observations", type=str, help="Vadose observations CSV for calibration (timestamp,depth_m,theta).")
    parser.add_argument("--vadose-calibrate-max-nfev", type=int, default=25, help="Max function evaluations for vadose calibration (default: 25).")
    parser.add_argument("--vadose-calibrate-output", type=str, help="Optional path to write vadose calibration result JSON.")
    parser.add_argument("--vadose-no3-loading", type=str, help="Optional nitrate loading CSV (timestamp,NO3_mg_L) to predict breakthrough to water table.")
    parser.add_argument("--vadose-no3-k", type=float, default=0.0, help="First-order nitrate attenuation rate k (1/day) for vadose breakthrough (default: 0).")
    parser.add_argument("--vadose-no3-output", type=str, help="Optional output CSV for vadose nitrate breakthrough curves (long format).")
    parser.add_argument("--vadose-no3-summary-output", type=str, help="Optional output JSON for vadose nitrate breakthrough summaries per edge.")

    # Reactive transport validation arguments
    parser.add_argument("--validate-forward", action="store_true", help="Run forward RT validation")
    parser.add_argument("--rt-simulator", choices=["phreeqc_kinetic", "mt3dms"], default="phreeqc_kinetic", help="Reactive transport simulator")
    parser.add_argument("--rt-time-steps", type=int, default=100, help="Number of kinetic time steps (default: 100)")
    parser.add_argument("--rt-rmse-threshold", type=float, default=1.0, help="RMSE threshold for consistency (default: 1.0)")
    parser.add_argument("--rt-nse-threshold", type=float, default=0.5, help="NSE threshold for consistency (default: 0.5)")
    parser.add_argument("--rt-residence-time", type=float, help="Default residence time in days (if not computed)")
    parser.add_argument("--rt-custom-rates", type=str, help="Path to custom rate laws YAML")

    # 3D network arguments
    parser.add_argument("--3d", action="store_true", dest="network_3d", help="Enable 3D network")
    parser.add_argument("--z-key", type=str, default="screen_depth", help="Z coordinate column name")
    parser.add_argument("--anisotropy", type=float, default=0.1, help="Vertical anisotropy factor α_v (default: 0.1)")
    parser.add_argument("--layer-file", type=str, help="Layer definition YAML file")
    parser.add_argument("--aquitard-p", type=float, default=0.3, help="Aquitard leakage probability (default: 0.3)")
    parser.add_argument("--screen-overlap-threshold", type=float, default=5.0, help="Screen overlap threshold in meters (default: 5.0)")

    # Mineral Library Options
    parser.add_argument(
        "--minerals",
        type=str,
        help="Comma-separated list of active minerals (e.g., 'calcite,albite,pyrite_oxidation_aerobic'). by default it uses a standard set.",
    )
    parser.add_argument(
        "--list-minerals",
        action="store_true",
        help="Print available minerals in the library and exit.",
    )

    args = parser.parse_args()

    if args.list_minerals:
        print("Available Minerals in Library:")
        for name, stoich in sorted(MINERAL_LIBRARY.items()):
            formula = ", ".join(f"{ion}:{coeff}" for ion, coeff in stoich.items())
            print(f"  - {name}: {formula}")
        return

    if args.report_only:
        with open(args.report_only, 'r') as f:
            data = json.load(f)
        print_interpretation_report(data)
        return

    if not args.samples or not args.output:
        parser.error("--samples and --output are required unless --report-only or --list-minerals is used.")

    ion_order = _parse_ion_order(args.ion_order) if args.ion_order else DEFAULT_ION_ORDER.copy()
    samples = _read_samples(args.samples)
    _validate_required(samples, ion_order)

    # Optional physics priors (file-based and/or MODPATH endpoints).
    physics_priors = []
    if args.physics_priors:
        from .physics.priors import load_physics_priors

        physics_priors.extend(load_physics_priors(args.physics_priors))
    if args.modpath_endpoints:
        from .physics.modpath import node_coords_from_samples, priors_from_modpath_endpoints

        node_coords = node_coords_from_samples(samples, node_id_key="site_id", x_key="x", y_key="y", z_key=None)
        if not node_coords:
            raise ValueError("--modpath-endpoints requires samples to include x and y coordinates in the same system as MODPATH.")
        physics_priors.extend(
            priors_from_modpath_endpoints(
                args.modpath_endpoints,
                node_coords,
                max_match_distance=args.modpath_max_match_distance,
                time_unit_days=args.modpath_time_unit_days,
            )
        )

    # Optional vadose-zone physics priors (Richards column -> travel time summaries).
    vadose_details = None
    if args.vadose_enabled:
        if not args.vadose_forcing or not args.vadose_profile or not args.vadose_links:
            raise ValueError("--vadose-enabled requires --vadose-forcing, --vadose-profile, and --vadose-links")
        from .vadose.io import load_forcing_csv, load_links_csv, load_profile, load_water_table_csv
        from .vadose.contracts import VadoseRunConfig
        from .vadose.calibrate import calibrate_ks_and_kc, load_theta_observations_csv, scale_profile_ks
        from .vadose.run import build_vadose_edge_priors

        forcing = load_forcing_csv(args.vadose_forcing)
        profile = load_profile(args.vadose_profile)
        links = load_links_csv(args.vadose_links)
        wt_series = load_water_table_csv(args.vadose_water_table) if args.vadose_water_table else None
        vadose_cfg = VadoseRunConfig(
            dz_m=float(args.vadose_dz),
            kc=float(args.vadose_kc),
            evaporation_fraction=float(args.vadose_evap_frac),
            lower_boundary=str(args.vadose_lower_boundary),
        )
        if args.vadose_calibrate:
            if not args.vadose_observations:
                raise ValueError("--vadose-calibrate requires --vadose-observations (timestamp,depth_m,theta)")
            from dataclasses import replace

            obs = load_theta_observations_csv(args.vadose_observations)
            cal = calibrate_ks_and_kc(
                profile,
                forcing,
                obs,
                config=vadose_cfg,
                water_table_depth_m=wt_series,
                max_nfev=int(args.vadose_calibrate_max_nfev),
            )
            profile = scale_profile_ks(profile, cal.ks_multiplier)
            vadose_cfg = replace(vadose_cfg, kc=float(cal.kc))
            if args.vadose_calibrate_output:
                with open(args.vadose_calibrate_output, "w", encoding="utf-8") as handle:
                    json.dump(
                        {
                            "ks_multiplier": cal.ks_multiplier,
                            "kc": cal.kc,
                            "cost": cal.cost,
                            "nfev": cal.nfev,
                            "status": cal.status,
                            "message": cal.message,
                        },
                        handle,
                        indent=2,
                    )
        vadose_edge_priors, vadose_details = build_vadose_edge_priors(
            profile,
            forcing,
            links,
            config=vadose_cfg,
            water_table_depth_m=wt_series,
        )
        physics_priors.extend([p.to_physics_prior() for p in vadose_edge_priors])
        if args.vadose_priors_output:
            # Write as a physics-priors CSV compatible with hydrosheaf/physics/priors.py
            import csv as _csv

            with open(args.vadose_priors_output, "w", newline="", encoding="utf-8") as handle:
                writer = _csv.DictWriter(
                    handle,
                    fieldnames=["u", "v", "p_uv", "tt_mean_days", "tt_std_days", "tt_p10_days", "tt_p90_days", "source"],
                )
                writer.writeheader()
                for p in vadose_edge_priors:
                    writer.writerow(
                        {
                            "u": p.u,
                            "v": p.v,
                            "p_uv": p.p_uv,
                            "tt_mean_days": p.tt_mean_days,
                            "tt_std_days": p.tt_std_days,
                            "tt_p10_days": p.tt_p10_days,
                            "tt_p90_days": p.tt_p90_days,
                            "source": p.source,
                        }
                    )
        if args.vadose_details_output and vadose_details is not None:
            with open(args.vadose_details_output, "w", encoding="utf-8") as handle:
                json.dump(vadose_details, handle, indent=2)

        # Optional nitrate breakthrough curves to the water table (per vadose link).
        if args.vadose_no3_loading and (args.vadose_no3_output or args.vadose_no3_summary_output):
            from datetime import datetime
            from .vadose.nitrate import load_no3_loading_csv, predict_no3_breakthrough
            import csv as _csv

            loading = load_no3_loading_csv(args.vadose_no3_loading)
            run_meta = (vadose_details or {}).get("_vadose_run", {})
            run_ts = [datetime.strptime(str(s), "%Y-%m-%d") for s in (run_meta.get("timestamps") or [])]
            run_recharge = [float(v) for v in (run_meta.get("recharge_m_day") or [])]
            if not run_ts or len(run_ts) != len(run_recharge):
                raise ValueError("Vadose run metadata missing timestamps/recharge for NO3 breakthrough; re-run with --vadose-enabled.")

            summaries = []
            rows = []
            for p in vadose_edge_priors:
                edge_id = f"{p.u}->{p.v}"
                edge_detail = (vadose_details or {}).get(edge_id) or {}
                tau_grid = edge_detail.get("ttd_tau_days") or []
                g = edge_detail.get("ttd_g") or []
                if not tau_grid or not g:
                    continue
                points, summary = predict_no3_breakthrough(
                    edge_id=edge_id,
                    ttd_tau_days=tau_grid,
                    ttd_g=g,
                    timestamps=run_ts,
                    recharge_m_day=run_recharge,
                    loading=loading,
                    k_1_day=float(args.vadose_no3_k),
                    c_in_units="mg/L",
                )
                summaries.append(
                    {
                        "edge_id": summary.edge_id,
                        "k_1_day": summary.k_1_day,
                        "c_in_units": summary.c_in_units,
                        "total_mass_delivered": summary.total_mass_delivered,
                        "total_mass_no_attenuation": summary.total_mass_no_attenuation,
                        "attenuated_fraction": summary.attenuated_fraction,
                        "peak_c_wt": summary.peak_c_wt,
                        "peak_time": summary.peak_time,
                    }
                )
                for pt in points:
                    rows.append(
                        {
                            "edge_id": edge_id,
                            "timestamp": pt.timestamp.strftime("%Y-%m-%d"),
                            "no3_c_wt": pt.c_wt,
                            "recharge_m_day": pt.recharge_m_day,
                            "no3_flux_mass_per_m2_day": pt.flux_mass_per_m2_day,
                        }
                    )

            if args.vadose_no3_output and rows:
                with open(args.vadose_no3_output, "w", newline="", encoding="utf-8") as handle:
                    writer = _csv.DictWriter(
                        handle,
                        fieldnames=["edge_id", "timestamp", "no3_c_wt", "recharge_m_day", "no3_flux_mass_per_m2_day"],
                    )
                    writer.writeheader()
                    writer.writerows(rows)
            if args.vadose_no3_summary_output:
                with open(args.vadose_no3_summary_output, "w", encoding="utf-8") as handle:
                    json.dump(summaries, handle, indent=2)

    config = Config(
        ion_order=ion_order,
        lambda_sparse=args.lambda_sparse,
        lambda_l1=args.lambda_l1,
        allow_signed_reactions=args.allow_signed,
        reaction_max_iter=args.reaction_max_iter,
        reaction_tol=args.reaction_tol,
        charge_balance_limit=args.charge_balance_limit,
        ec_tds_penalty_limit=args.ec_tds_penalty_limit,
        ec_tds_penalty_enabled=args.ec_tds_penalty_enabled,
        missing_policy=args.missing_policy,
        detection_limit_policy=args.detection_policy,
        eta_ec=args.eta_ec,
        eta_tds=args.eta_tds,
        unit="mmol/L",
        phreeqc_enabled=args.phreeqc_enabled,
        phreeqc_mode=args.phreeqc_mode,
        phreeqc_database=args.phreeqc_database,
        phreeqc_executable=args.phreeqc_executable,
        temp_default_c=args.temp_default_c,
        si_threshold_tau=args.si_threshold,
        constraints_hard=args.constraints_hard,
        edge_p_min=args.edge_p_min,
        edge_radius_km=args.edge_radius_km,
        edge_max_neighbors=args.edge_max_neighbors,
        edge_sigma_meas=args.edge_sigma_meas,
        edge_sigma_dtw=args.edge_sigma_dtw,
        edge_sigma_elev=args.edge_sigma_elev,
        edge_sigma_topo=args.edge_sigma_topo,
        edge_head_inference=args.edge_head_inference,
        edge_dtw_prior_mu=args.edge_dtw_prior_mu,
        edge_dtw_prior_sigma=args.edge_dtw_prior_sigma,
        edge_head_prior_mu=args.edge_head_prior_mu,
        edge_head_prior_sigma=args.edge_head_prior_sigma,
        edge_topo_sigma_depth=args.edge_topo_sigma_depth,
        edge_topo_reject_p=args.edge_topo_reject_p,
        edge_map_prior_weight=args.edge_map_weight,
        edge_map_candidate_multiplier=args.edge_map_candidate_multiplier,
        edge_map_p_min=args.edge_map_p_min,
        edge_gradient_min=args.edge_gradient_min,
        edge_head_key=args.edge_head_key,
        edge_dtw_key=args.edge_dtw_key,
        edge_elevation_key=args.edge_elevation_key,
        edge_aquifer_key=args.edge_aquifer_key,
        edge_screen_depth_key=args.edge_screen_depth_key,
        edge_well_depth_key=args.edge_well_depth_key,
        edge_depth_mismatch=args.edge_depth_mismatch,
        isotope_enabled=args.isotope_enabled,
        isotope_weight=args.isotope_weight,
        isotope_d_excess_weight=args.isotope_d_excess_weight,
        isotope_d18o_key=args.isotope_d18o_key,
        nitrate_source_enabled=args.nitrate_source_enabled,
        uncertainty_method=args.uncertainty,
        bootstrap_n_resamples=args.bootstrap_n,
        bootstrap_ci_method=args.bootstrap_ci,
        bayesian_n_samples=args.bayesian_samples,
        bayesian_n_chains=args.bayesian_chains,
        bayesian_target_accept=args.bayesian_accept,
        monte_carlo_n_samples=args.monte_carlo_n,
        input_uncertainty_pct=args.input_uncertainty,
        temporal_enabled=args.temporal_enabled,
        temporal_window_days=args.temporal_window,
        temporal_min_samples=args.temporal_min_samples,
        temporal_interpolation_method=args.temporal_interp,
        temporal_frequency_days=args.temporal_frequency,
        residence_time_method=args.residence_method if args.residence_method else "cross_correlation",
        residence_time_tracer=args.residence_tracer,
        residence_time_hydraulic_k=args.residence_k,
        residence_time_porosity=args.residence_porosity,
        tau_agreement_tolerance=args.tau_agreement_tolerance,
        tau_min_peak_corr=args.tau_min_peak_corr,
        tau_max_relative_uncertainty=args.tau_max_relative_uncertainty,
        tau_max_uncertainty_days=args.tau_max_uncertainty_days,
        tau_physics_blend_threshold=args.tau_physics_blend_threshold,
        ttd_grid_dt_days=args.ttd_grid_dt_days,
        ttd_max_lag_days=args.ttd_max_lag_days,
        ttd_smoothness_lambda=args.ttd_smoothness_lambda,
        ttd_min_r2=args.ttd_min_r2,
        ttd_attenuation_k_max=args.ttd_attenuation_k_max,
        ttd_attenuation_k_steps=args.ttd_attenuation_k_steps,
        bayes_lag_grid_dt_days=args.bayes_lag_grid_dt_days,
        bayes_lag_max_lag_days=args.bayes_lag_max_lag_days,
        bayes_lag_prior_sigma_multiplier=args.bayes_lag_prior_sigma_multiplier,
        bayes_lag_min_pairs=args.bayes_lag_min_pairs,
        residence_time_coupling_enabled=args.residence_coupling,
        residence_time_reference_days=args.residence_ref_days,
        reactive_transport_validation=args.validate_forward,
        rt_simulator=args.rt_simulator,
        rt_n_time_steps=args.rt_time_steps,
        rt_consistency_rmse_threshold=args.rt_rmse_threshold,
        rt_consistency_nse_threshold=args.rt_nse_threshold,
        rt_default_residence_time=args.rt_residence_time if args.rt_residence_time else 30.0,
        rt_custom_rates_file=args.rt_custom_rates if args.rt_custom_rates else "",
        network_3d_enabled=args.network_3d,
        z_coordinate_key=args.z_key,
        vertical_anisotropy=args.anisotropy,
        aquitard_leakage_p=args.aquitard_p,
        screen_overlap_threshold=args.screen_overlap_threshold,
    )
    if args.layer_minerals:
        config.layer_mineral_map = _parse_layer_minerals(args.layer_minerals)
    if args.nitrate_source_min_conc is not None:
        config.nitrate_source_min_mg_L = args.nitrate_source_min_conc
    config.isotope_d2h_key = args.isotope_d2h_key
    config.auto_lmwl = args.auto_lmwl
    config.isotope_consistency_weight = args.iso_consistency_weight
    
    # Edges / topology
    edges = None
    if args.edges:
        edges = _read_edges(args.edges)
    elif args.infer_edges:
        edge_attr_overrides = None
        if physics_priors:
            edge_attr_overrides = {f"{p.u}->{p.v}": p.attrs() for p in physics_priors}
        edges = infer_edges(
            samples,
            max_neighbors=args.max_neighbors,
            allow_uphill=args.allow_uphill,
            head_key=args.head_key,
            elevation_key=args.elevation_key,
            method=args.infer_edges_method,
            config=config,
            edge_attr_overrides=edge_attr_overrides,
        )
    elif physics_priors and args.physics_priors_mode == "only":
        edges = []
    else:
        raise ValueError("Provide --edges, use --infer-edges, or use --physics-priors-mode only with priors.")

    if physics_priors:
        from .physics.priors import apply_physics_priors

        edges = apply_physics_priors(list(edges or []), physics_priors, mode=args.physics_priors_mode)
        if not edges:
            raise ValueError("Physics priors produced no usable edges (check node ids and inputs).")
    if args.weights:
        config.weights = _parse_weights(args.weights)
    if args.eta_ec or args.eta_tds:
        config.ec_tds_penalty_enabled = True
    if args.isotope_enabled:
        config.isotope_enabled = True
    if args.endmember:
        config.mixing_endmembers = _parse_endmembers(args.endmember)
    samples = normalize_samples(samples, ion_order, config.detection_limit_policy)
    if args.endmember_id:
        endmember_map = build_endmember_vectors(samples, args.endmember_id, ion_order)
        missing_ids = [end_id for end_id in args.endmember_id if end_id not in endmember_map]
        if missing_ids:
            raise ValueError(f"Endmember IDs not found in samples: {missing_ids}")
        for end_id, vector in endmember_map.items():
            config.mixing_endmembers.setdefault(end_id, vector)
    if args.endmembers_json:
        endmember_map, meta = load_endmembers_json(args.endmembers_json)
        if meta.get("ion_order") and meta.get("ion_order") != ion_order:
            raise ValueError("endmembers.json ion_order does not match CLI ion_order")
        for end_id, vector in endmember_map.items():
            config.mixing_endmembers.setdefault(end_id, vector)
    if args.signed_reaction:
        config.signed_reaction_labels = args.signed_reaction
    
    if args.minerals:
        config.active_minerals = [m.strip() for m in args.minerals.split(",") if m.strip()]


    if args.unit != "mmol/L":
        samples = _convert_samples(samples, ion_order, args.unit)
        if config.mixing_endmembers:
            config.mixing_endmembers = _convert_endmembers(
                config.mixing_endmembers,
                ion_order,
                args.unit,
            )

    if args.calibrate_ec_tds:
        calibrate_ec_tds(samples, config)
        config.ec_tds_penalty_enabled = True

    if args.lmwl_a is not None and args.lmwl_b is not None:
        config.lmwl_a = args.lmwl_a
        config.lmwl_b = args.lmwl_b
        config.lmwl_defined = True
        config.isotope_enabled = True
    elif args.fit_lmwl:
        a, b = fit_lmwl(samples, d18o_key=config.isotope_d18o_key, d2h_key=config.isotope_d2h_key)
        config.lmwl_a = a
        config.lmwl_b = b
        config.lmwl_defined = True
        config.isotope_enabled = True

    if config.isotope_enabled and config.auto_lmwl:
        try:
            intercept, slope = fit_lmwl(
                samples,
                d18o_key=config.isotope_d18o_key,
                d2h_key=config.isotope_d2h_key
            )
            config.lmwl_a = intercept
            config.lmwl_b = slope
            config.lmwl_defined = True
            print(f"Auto-calibrated LMWL: d2H = {slope:.2f} * d18O + {intercept:.2f}")
        except Exception as e:
            print(f"Warning: Failed to auto-calibrate LMWL: {e}")

    temporal_results = []
    temporal_by_edge = {}
    tau_overrides = None
    if args.temporal_enabled:
        if not args.temporal_data:
            raise ValueError("--temporal-enabled requires --temporal-data")
        from .temporal.time_series import load_time_series_csv
        from .temporal.temporal_edge_fit import fit_temporal_edge
        from .graph.build import build_edges

        # Prefer site_id to match inferred edge node ids; fall back to other common keys.
        temporal_nodes = load_time_series_csv(args.temporal_data, ion_order, node_id_column="site_id")
        if not temporal_nodes:
            temporal_nodes = load_time_series_csv(args.temporal_data, ion_order, node_id_column="node_id")
        if not temporal_nodes:
            temporal_nodes = load_time_series_csv(args.temporal_data, ion_order, node_id_column="sample_id")

        tau_overrides = {}
        for edge in build_edges(edges):
            if edge.u not in temporal_nodes or edge.v not in temporal_nodes:
                continue
            hydraulic_params = {
                "lmwl_intercept": float(getattr(config, "lmwl_a", 10.0)),
                "lmwl_slope": float(getattr(config, "lmwl_b", 8.0)),
                "agreement_tolerance": float(getattr(config, "tau_agreement_tolerance", 0.4)),
                "min_peak_corr": float(getattr(config, "tau_min_peak_corr", 0.2)),
                "max_relative_uncertainty": float(getattr(config, "tau_max_relative_uncertainty", 1.5)),
                "max_uncertainty_days": float(getattr(config, "tau_max_uncertainty_days", 180.0)),
                "grid_dt_days": float(getattr(config, "ttd_grid_dt_days", 30.0)),
                "max_lag_days": float(getattr(config, "ttd_max_lag_days", 365.0)),
                "smoothness_lambda": float(getattr(config, "ttd_smoothness_lambda", 0.0)),
                "ttd_min_r2": float(getattr(config, "ttd_min_r2", 0.2)),
                "attenuation_k_max": float(getattr(config, "ttd_attenuation_k_max", 0.02)),
                "attenuation_k_steps": int(getattr(config, "ttd_attenuation_k_steps", 6)),
                "bayes_lag_grid_dt_days": float(getattr(config, "bayes_lag_grid_dt_days", 5.0)),
                "bayes_lag_max_lag_days": float(getattr(config, "bayes_lag_max_lag_days", 365.0)),
                "bayes_lag_prior_sigma_multiplier": float(getattr(config, "bayes_lag_prior_sigma_multiplier", 1.0)),
                "bayes_lag_min_pairs": int(getattr(config, "bayes_lag_min_pairs", 5)),
            }
            # Provide physics fallback/prior when edge geometry supports it (typical when edges are inferred).
            attrs = edge.attrs or {}
            # If physics travel-time priors are available, use them directly as a tau prior for bayesian_lag.
            tau_mu = attrs.get("physics_tau_mean_days", attrs.get("edge_residence_time_days"))
            tau_sigma = attrs.get("physics_tau_std_days")
            tau_p10 = attrs.get("physics_tau_p10_days")
            tau_p90 = attrs.get("physics_tau_p90_days")
            if tau_mu is not None:
                try:
                    hydraulic_params["tau_prior_mu_days"] = float(tau_mu)  # type: ignore[arg-type]
                except (TypeError, ValueError):
                    pass
            if tau_sigma is not None:
                try:
                    hydraulic_params["tau_prior_sigma_days"] = float(tau_sigma)  # type: ignore[arg-type]
                except (TypeError, ValueError):
                    pass
            elif tau_p10 is not None and tau_p90 is not None:
                try:
                    p10 = float(tau_p10)  # type: ignore[arg-type]
                    p90 = float(tau_p90)  # type: ignore[arg-type]
                    # Approximate sigma from central 80% interval for a normal: width ≈ 2.563*sigma.
                    hydraulic_params["tau_prior_sigma_days"] = max(1e-9, (p90 - p10) / 2.563)
                except (TypeError, ValueError):
                    pass

            distance_m = attrs.get("distance_3d_m")
            if distance_m is None and attrs.get("distance_km") is not None:
                try:
                    distance_m = float(attrs.get("distance_km")) * 1000.0
                except (TypeError, ValueError):
                    distance_m = None
            if distance_m is not None:
                try:
                    hydraulic_params["distance_m"] = float(distance_m)
                except (TypeError, ValueError):
                    pass
            delta_h = attrs.get("delta_h")
            if delta_h is not None and distance_m:
                try:
                    hydraulic_params["gradient"] = abs(float(delta_h)) / max(1e-9, float(distance_m))
                except (TypeError, ValueError):
                    pass
            h_grad = attrs.get("horizontal_gradient")
            if "gradient" not in hydraulic_params and h_grad is not None:
                try:
                    hydraulic_params["gradient"] = abs(float(h_grad))
                except (TypeError, ValueError):
                    pass
            hydraulic_params["K_m_day"] = float(getattr(config, "residence_time_hydraulic_k", 1.0))
            hydraulic_params["porosity"] = float(getattr(config, "residence_time_porosity", 0.2))
            # Allow per-edge overrides if provided as edge attributes (e.g., from an edges CSV).
            for key, out_key in [("K_m_day", "K_m_day"), ("porosity", "porosity"), ("distance_m", "distance_m"), ("gradient", "gradient")]:
                if key in attrs and attrs.get(key) not in (None, ""):
                    try:
                        hydraulic_params[out_key] = float(attrs.get(key))  # type: ignore[arg-type]
                    except (TypeError, ValueError):
                        pass

            hydraulic_params["physics_blend_threshold"] = float(getattr(config, "tau_physics_blend_threshold", 0.5))

            temporal_res = fit_temporal_edge(
                temporal_nodes[edge.u],
                temporal_nodes[edge.v],
                config,
                edge_id=edge.edge_id,
                hydraulic_params=hydraulic_params,
            )
            temporal_results.append(temporal_res)
            temporal_by_edge[edge.edge_id] = temporal_res
            tau_overrides[edge.edge_id] = temporal_res.residence_time_days

    results = fit_network(samples, edges, config, residence_time_overrides=tau_overrides)

    # Attach temporal summaries to the main per-edge output rows (when available).
    if temporal_by_edge:
        for res in results:
            t = temporal_by_edge.get(res.edge_id)
            if t is None:
                continue
            res.temporal_residence_time_days = t.residence_time_days
            res.temporal_residence_time_method = t.residence_time_method
            res.temporal_residence_time_uncertainty = t.residence_time_uncertainty
            res.temporal_transport_model = t.transport_model
            res.temporal_gamma_mean = t.gamma_mean
            res.temporal_gamma_std = t.gamma_std
            res.temporal_f_mean = t.f_mean
            res.temporal_f_std = t.f_std
            res.temporal_reaction_extents_mean = list(t.reaction_extents_mean or [])
            res.temporal_reaction_extents_std = list(t.reaction_extents_std or [])
            res.temporal_total_residual_norm = t.total_residual_norm
            res.temporal_n_time_points = len(t.timestamps or [])
            res.temporal_residence_time_flags = list(t.residence_time_flags or [])
            res.temporal_residence_time_details = dict(t.residence_time_details or {})

    if args.format == "csv":
        export_edge_results_csv(results, args.output)
    else:
        export_edge_results_json(results, args.output)

    if args.interpret:
        # Convert results to list of dicts for the interpreter if it's not already
        from .outputs.tables import edge_results_table
        table_data = edge_results_table(results)
        print_interpretation_report(table_data)

    # Optional temporal sidecar output (time-varying edge fits)
    if args.temporal_enabled:
        temporal_out = args.temporal_output
        if not temporal_out:
            base = args.output
            if base.lower().endswith(".json"):
                temporal_out = base[:-5] + ".temporal.json"
            else:
                temporal_out = base + ".temporal.json"
        with open(temporal_out, "w", encoding="utf-8") as handle:
            json.dump(temporal_edge_results_to_rows(temporal_results), handle, indent=2)


if __name__ == "__main__":
    main()
