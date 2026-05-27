"""Read MODFLOW ASCII formatted head (.fhd) files and compute grid-aware diagnostics.

Parses layer-by-layer head blocks and returns a mapping from 1-based
sequential MODFLOW cell index to head value, suitable for joining onto
a node-mapping DataFrame via the "cell_{int}" node_id convention.

Also provides continuous head-gradient computation via finite differences
on the structured MODFLOW grid, plus MODFLOW 6 binary .hds support.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from statistics import mode
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Grid geometry
# ---------------------------------------------------------------------------


@dataclass
class GridGeometry:
    """MODFLOW structured-grid geometry needed for finite-difference gradients.

    Parameters
    ----------
    ncol, nrow, nlay:
        Number of columns, rows, and layers.
    dx, dy:
        Cell width and height in model length units (e.g. feet for Savage).
        When the grid has variable spacing these may be set to the *core*
        uniform spacing; all nodes in the node-mapping DataFrame are assumed
        to lie within this uniform region.
    rotation_deg:
        Grid rotation angle in degrees counter-clockwise from the projected
        coordinate system axes.  MODFLOW ANGROT convention (positive = CCW).
    origin_x, origin_y:
        Model origin (lower-left corner of the grid) in projected coordinates.
    cell_x, cell_y:
        Optional (ncol, nrow) arrays of cell-centre coordinates in model space.
        When *None*, they are built from *dx* / *dy* assuming uniform spacing.
    """

    ncol: int
    nrow: int
    nlay: int = 1
    dx: float = 0.0
    dy: float = 0.0
    rotation_deg: float = 0.0
    origin_x: float = 0.0
    origin_y: float = 0.0
    cell_x: Optional[np.ndarray] = None  # (ncol, nrow)
    cell_y: Optional[np.ndarray] = None  # (ncol, nrow)

    def __post_init__(self) -> None:
        if self.cell_x is None or self.cell_y is None:
            xc = np.arange(self.ncol, dtype=np.float64) * self.dx + self.dx / 2.0
            yc = np.arange(self.nrow, dtype=np.float64) * self.dy + self.dy / 2.0
            self.cell_x, self.cell_y = np.meshgrid(xc, yc, indexing="ij")


# ---------------------------------------------------------------------------
# FHD header detection
# ---------------------------------------------------------------------------


def _try_parse_header(tokens: list[str]) -> Optional[Tuple[int, int, int]]:
    """Return (ncol, nrow, ilay) from a MODFLOW .fhd header line.

    Handles two common formats:

    Classic MODFLOW-2000:
        KSTP  KPER  PERTIM  TOTIM  NCOL  NROW  ILAY

    MODFLOW-NWT / MODFLOW-2005 (extended):
        KSTP  KPER  PERTIM  TOTIM  TEXT  NCOL  NROW  ILAY  [FORMAT]

    The key insight: NCOL, NROW, ILAY are always the last run of positive
    integers on the line, immediately before any optional format string.
    We collect integers by scanning left-to-right after skipping the first
    two (KSTP, KPER) and any floats / text tokens; the first three
    consecutive positive integers encountered after any non-numeric text
    are NCOL, NROW, ILAY.
    """
    if len(tokens) < 7:
        return None
    try:
        kstp = int(tokens[0])
        kper = int(tokens[1])
        if kstp <= 0 or kper <= 0:
            return None
    except (ValueError, IndexError):
        return None

    # Collect the first run of positive integers that appear after the
    # float/text tokens (PERTIM, TOTIM, and optional TEXT label).
    int_run: list[int] = []
    for tok in tokens[2:]:
        try:
            v = int(tok)
            if v > 0:
                int_run.append(v)
            else:
                # Zero or negative integer breaks the run
                if int_run:
                    break
        except ValueError:
            # Float (PERTIM / TOTIM) or text (HEAD / format string)
            if int_run:
                # We already have integers — a non-int signals the end
                break
            # Otherwise keep scanning past floats / text tokens

    if len(int_run) >= 3:
        ncol, nrow, ilay = int_run[0], int_run[1], int_run[2]
        return ncol, nrow, ilay
    return None


# ---------------------------------------------------------------------------
# FHD parser
# ---------------------------------------------------------------------------


def parse_fhd(
    path: str | Path,
    hdry_upper: float = 1e28,
    hdry_lower: Optional[float] = None,
) -> Dict[int, float]:
    """Parse a MODFLOW ASCII formatted head file (.fhd).

    Returns a dict mapping 1-based sequential MODFLOW cell index to head
    value (in whatever units the model uses — feet for the Savage NH model).

    Dry / inactive cells are omitted from the result.  Two exclusion rules
    apply (either triggers exclusion):
    - ``val >= hdry_upper`` — standard MODFLOW HDRY marker (default 1e28).
    - ``val <= hdry_lower`` — negative HDRY marker, e.g. -300 ft as used by
      MODFLOW-NWT models.  Defaults to None (no lower threshold).

    The sequential cell index matches MODFLOW's own numbering convention:
    cells are counted layer-by-layer, row-by-row, column-by-column (1-based).
    Cell n in layer L, row R, column C satisfies:
        n = (L-1)*NROW*NCOL + (R-1)*NCOL + C
    which is exactly the order they appear in the .fhd output.
    """
    path = Path(path)
    heads: Dict[int, float] = {}

    # Auto-detect negative HDRY from the first few data lines if not provided.
    # For MODFLOW-NWT the typical HDRY is -300.0; for MODFLOW-2000 it is 1e30.
    # We scan the first block of the file and if all non-header values in
    # that block are exactly one repeated negative number, we treat it as HDRY.
    _hdry_lower = hdry_lower

    with path.open("r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]

    if _hdry_lower is None:
        _hdry_lower = _detect_hdry_lower(lines)

    i = 0
    while i < len(lines):
        tokens = lines[i].split()
        parsed = _try_parse_header(tokens)

        if parsed is None:
            i += 1
            continue

        ncol, nrow, ilay = parsed

        # Skip optional "HEAD IN LAYER..." description line (classic format)
        i += 1
        if i < len(lines) and lines[i].lstrip().upper().startswith("HEAD"):
            i += 1

        # Read ncol * nrow head values (may span multiple lines)
        n_cells = ncol * nrow
        vals: list[float] = []
        while len(vals) < n_cells and i < len(lines):
            # Stop if next line looks like a header (don't consume it)
            next_tokens = lines[i].split()
            if _try_parse_header(next_tokens) is not None:
                break
            for tok in next_tokens:
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
            i += 1

        # Map values to 1-based sequential cell index for this layer
        layer_offset = (ilay - 1) * ncol * nrow
        for j, val in enumerate(vals[:n_cells]):
            dry = val >= hdry_upper
            if _hdry_lower is not None:
                dry = dry or val <= _hdry_lower
            if not dry:
                heads[layer_offset + j + 1] = val

    return heads


def _detect_hdry_lower(lines: list[str]) -> Optional[float]:
    """Heuristically detect a negative HDRY value from the first data block."""
    # Read up to 500 numeric values after the first header line
    in_data = False
    vals: list[float] = []
    for line in lines[:500]:
        toks = line.split()
        if not in_data:
            if _try_parse_header(toks) is not None:
                in_data = True
            continue
        for tok in toks:
            try:
                vals.append(float(tok))
            except ValueError:
                pass
        if len(vals) >= 200:
            break

    if not vals:
        return None
    negative = [v for v in vals if v < 0]
    if not negative:
        return None
    # If ≥80% of sampled values are the same negative number → treat as HDRY
    try:
        candidate = mode(round(v, 2) for v in negative)
    except Exception:
        return None
    count = sum(1 for v in vals if round(v, 2) == candidate)
    if count / len(vals) >= 0.50:
        # Unit-agnostic buffer: 1% of |HDRY|, floored at 0.5 model-units
        buffer = max(0.5, abs(float(candidate)) * 0.01)
        return float(candidate) + buffer
    return None


def cell_heads_to_sample_field(
    node_df: pd.DataFrame,
    head_map: Dict[int, float],
    field_name: str = "hydraulic_head",
) -> pd.DataFrame:
    """Attach MODFLOW head values to a node-mapping DataFrame.

    Parameters
    ----------
    node_df:
        DataFrame with a ``node_id`` column containing values like ``"cell_130749"``.
    head_map:
        Output of :func:`parse_fhd`.
    field_name:
        Column name to add (default ``"hydraulic_head"``).

    Returns
    -------
    DataFrame with *field_name* added (NaN where no head is available).
    """
    out = node_df.copy()
    cell_nums = (
        out["node_id"]
        .str.extract(r"cell_(\d+)", expand=False)
        .astype("Int64")
    )
    out[field_name] = cell_nums.map(head_map)
    return out


# ---------------------------------------------------------------------------
# Continuous head gradient from structured grid
# ---------------------------------------------------------------------------

def _cell_to_ijk(cell_idx: int, ncol: int, nrow: int) -> Tuple[int, int, int]:
    """Convert 1-based MODFLOW cell index to 0-based (col, row, layer)."""
    idx = cell_idx - 1
    cells_per_layer = ncol * nrow
    lay = idx // cells_per_layer
    remainder = idx % cells_per_layer
    row = remainder // ncol
    col = remainder % ncol
    return col, row, lay


def _head_array_from_map(
    head_map: Dict[int, float],
    ncol: int,
    nrow: int,
    nlay: int,
) -> np.ndarray:
    """Build (ncol, nrow, nlay) float32 array from cell-index → head dict.

    Cells not present in *head_map* are set to NaN.
    """
    arr = np.full((ncol, nrow, nlay), np.nan, dtype=np.float32)
    for cell_idx, val in head_map.items():
        col, row, lay = _cell_to_ijk(cell_idx, ncol, nrow)
        if 0 <= col < ncol and 0 <= row < nrow and 0 <= lay < nlay:
            arr[col, row, lay] = float(val)
    return arr


def _gaussian_smooth_layer(
    layer: np.ndarray, sigma: float
) -> np.ndarray:
    """Gaussian-smooth a 2D layer, preserving NaN mask."""
    if sigma <= 0:
        return layer
    try:
        from scipy.ndimage import gaussian_filter  # type: ignore[import-untyped]
    except ImportError:
        return layer

    mask = ~np.isnan(layer)
    if mask.sum() < 3:
        return layer
    filled = np.where(mask, layer, 0.0)
    smoothed = gaussian_filter(filled, sigma=sigma, mode="nearest")
    weights = gaussian_filter(mask.astype(np.float32), sigma=sigma, mode="nearest")
    weights = np.where(weights < 1e-10, 1.0, weights)
    return np.where(mask, smoothed / weights, np.nan)


def _safe_gradient_1d(
    arr: np.ndarray, i: int, step: float, axis_size: int
) -> float:
    """Compute df/dx at index *i* with NaN-safe fallback.

    Prefers centred difference; falls back to one-sided if a neighbour is
    NaN, and returns ``math.nan`` when no valid neighbour exists.
    """
    left_ok = i > 0 and not np.isnan(arr[i - 1])
    right_ok = i < axis_size - 1 and not np.isnan(arr[i + 1])

    if left_ok and right_ok:
        return float((arr[i + 1] - arr[i - 1]) / (2.0 * step))
    if right_ok:
        return float((arr[i + 1] - arr[i]) / step)
    if left_ok:
        return float((arr[i] - arr[i - 1]) / step)
    return float("nan")


def compute_head_gradient(
    head_map: Dict[int, float],
    grid: GridGeometry,
    sigma: float = 1.0,
) -> Dict[int, Tuple[float, float, float]]:
    """Compute continuous hydraulic-head gradient at each active cell.

    Uses centred finite differences on the MODFLOW structured grid, with
    optional Gaussian pre-smoothing to suppress difference-amplified noise.

    Parameters
    ----------
    head_map:
        Output of :func:`parse_fhd` (or :func:`parse_hds`).
        Maps 1-based MODFLOW cell index → head value.
    grid:
        Grid geometry descriptor.
    sigma:
        Gaussian smoothing radius in cell-width units before differencing.
        Set to 0 to disable smoothing.

    Returns
    -------
    dict mapping 1-based cell index → ``(gx, gy, gz)`` where each component
    is in head-units per model-length-unit (e.g. ft/ft for Savage).
    The gradient is rotated into world (projected) coordinates — *gx* points
    along the projected x-axis, *gy* along projected y.
    *gz* is always 0.0 for single-layer grids.
    """
    head_array = _head_array_from_map(head_map, grid.ncol, grid.nrow, grid.nlay)

    # Gaussian smoothing per layer (horizontal only)
    if sigma > 0:
        for k in range(grid.nlay):
            head_array[:, :, k] = _gaussian_smooth_layer(head_array[:, :, k], sigma)

    # Rotation matrix for grid → world
    theta = math.radians(grid.rotation_deg)
    cos_t = math.cos(theta)
    sin_t = math.sin(theta)

    gradient_map: Dict[int, Tuple[float, float, float]] = {}

    for cell_idx in head_map:
        col, row, lay = _cell_to_ijk(cell_idx, grid.ncol, grid.nrow)
        if not (0 <= col < grid.ncol and 0 <= row < grid.nrow and 0 <= lay < grid.nlay):
            continue

        h_center = head_array[col, row, lay]
        if np.isnan(h_center):
            continue

        # Gradient in grid coordinates
        gx_grid = _safe_gradient_1d(head_array[:, row, lay], col, grid.dx, grid.ncol)
        gy_grid = _safe_gradient_1d(head_array[col, :, lay], row, grid.dy, grid.nrow)

        if np.isnan(gx_grid) or np.isnan(gy_grid):
            continue

        # Rotate to world coordinates (ANGROT positive = CCW from world axes)
        gx_world = gx_grid * cos_t - gy_grid * sin_t
        gy_world = gx_grid * sin_t + gy_grid * cos_t
        gz = 0.0

        # Vertical gradient if multi-layer
        if grid.nlay > 1:
            gz = _safe_gradient_1d(head_array[col, row, :], lay, 1.0, grid.nlay)

        gradient_map[cell_idx] = (float(gx_world), float(gy_world), float(gz))

    return gradient_map


def map_gradient_to_nodes(
    gradient_map: Dict[int, Tuple[float, float, float]],
    node_df: pd.DataFrame,
) -> Dict[str, Tuple[float, float, float]]:
    """Join cell-level gradient vectors onto node IDs.

    Parameters
    ----------
    gradient_map:
        Output of :func:`compute_head_gradient`.
    node_df:
        DataFrame with ``node_id`` column containing ``"cell_{int}"`` values.

    Returns
    -------
    ``{node_id: (gx, gy, gz)}`` for every node whose cell index exists in
    *gradient_map*.  Nodes without matching gradient data are omitted.
    """
    result: Dict[str, Tuple[float, float, float]] = {}
    cell_nums = (
        node_df["node_id"]
        .str.extract(r"cell_(\d+)", expand=False)
        .astype("Int64")
    )
    valid_mask = cell_nums.notna()
    valid_cells = cell_nums[valid_mask].astype(int)
    valid_ids = node_df.loc[valid_mask, "node_id"].astype(str)
    for cell_num, node_id in zip(valid_cells, valid_ids):
        grad = gradient_map.get(cell_num)
        if grad is not None:
            result[node_id] = grad
    return result


# ---------------------------------------------------------------------------
# MODFLOW 6 binary head support (.hds)
# ---------------------------------------------------------------------------


def parse_hds(
    path: str | Path,
    kstpkper: Optional[Tuple[int, int]] = None,
) -> Dict[int, float]:
    """Parse a MODFLOW 6 binary head-save file (.hds) via ``flopy``.

    Returns the same ``{cell_index: head_value}`` dict format as
    :func:`parse_fhd`, using MODFLOWʼs 1-based sequential cell numbering.

    Parameters
    ----------
    path:
        Path to the binary .hds file.
    kstpkper:
        Optional ``(kstp, kper)`` tuple to select a specific time step and
        stress period.  When *None*, the last time step is used.

    Returns
    -------
    dict mapping 1-based cell index → head value.
    Dry / inactive cells are excluded (flopy reports them as NaN).

    Raises
    ------
    ImportError
        If ``flopy`` is not installed.
    FileNotFoundError
        If *path* does not exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Head-save file not found: {path}")

    try:
        from flopy.utils import HeadFile  # type: ignore[import-untyped]
    except ImportError:
        raise ImportError(
            "flopy is required to read MODFLOW 6 binary .hds files. "
            "Install with: pip install flopy"
        )

    hf = HeadFile(str(path))
    if kstpkper is not None:
        kstp, kper = kstpkper
        head_array = hf.get_data(kstpkper=(kstp, kper))
    else:
        # Use the last record
        # NOTE: get_recordarray() loads metadata for all time steps.
        # For very large transient models this may be memory-intensive.
        # Consider using flopy's HeadFile.get_kstpkper() if record count is known.
        records = hf.get_recordarray()
        if records.shape[0] == 0:
            return {}
        head_array = hf.get_data(kstpkper=(records[-1]["kstp"], records[-1]["kper"]))

    # flopy returns a 2D (nrow, ncol) or 3D (nlay, nrow, ncol) array
    # Convert to the 1-based sequential cell index convention
    heads: Dict[int, float] = {}
    arr = np.asarray(head_array, dtype=np.float64)

    if arr.ndim == 2:
        # Single layer: (nrow, ncol)
        nrow, ncol = arr.shape
        for r in range(nrow):
            for c in range(ncol):
                val = arr[r, c]
                if not np.isnan(val) and np.isfinite(val):
                    cell = r * ncol + c + 1
                    heads[cell] = float(val)
    elif arr.ndim == 3:
        # Multi-layer: (nlay, nrow, ncol)
        nlay, nrow, ncol = arr.shape
        for l in range(nlay):
            for r in range(nrow):
                for c in range(ncol):
                    val = arr[l, r, c]
                    if not np.isnan(val) and np.isfinite(val):
                        cell = l * nrow * ncol + r * ncol + c + 1
                        heads[cell] = float(val)

    return heads


# ---------------------------------------------------------------------------
# Cell-by-cell flow (.cbc) — optional, best-effort
# ---------------------------------------------------------------------------


def try_parse_cbc(path: str | Path) -> Optional[np.ndarray]:
    """Attempt to parse MODFLOW cell-by-cell flow output (.cbc).

    Returns a structured NumPy array with columns ``(cell_from, cell_to, face,
    flow_rate)`` or *None* if parsing fails or the file is unavailable.

    Face-flow data provides the definitive physical flow-direction signal:
    a face with ``Q >> 0`` is one where water actually moves, regardless of
    the head-gradient sign.

    .. deprecated::
        Tested as an edge-direction prior on the Savage NH benchmark and found
        to be ineffective for topology inference: the dominant-outflow direction
        on a face is anti-correlated with particle trajectories (2.3% alignment,
        worse than useless).  MODPATH tracks magnitude-weighted flow on each
        face, not just the sign of the net flux, so the CBC net-flux direction
        does not reliably indicate the preferred downstream cell.

        This function is retained for diagnostic inspection of model budgets
        and may be useful for other model geometries or non-trajectory
        applications, but it should **not** be used as an edge-direction prior.

    Parameters
    ----------
    path:
        Path to the binary .cbc (or .ccf) file.

    Returns
    -------
    Structured array or *None*.
    """
    path = Path(path)
    if not path.exists():
        return None

    try:
        from flopy.utils import CellBudgetFile  # type: ignore[import-untyped]
    except ImportError:
        return None

    try:
        cbc = CellBudgetFile(str(path))
        records = cbc.get_data(text="FLOW JA FACE", full3D=False)
        if records is None or len(records) == 0:
            # Try alternative budget text labels
            for text in ("FLOW JA FACE", "FLOW-JA-FACE", "FLOW RIGHT FACE"):
                records = cbc.get_data(text=text, full3D=False)
                if records is not None and len(records) > 0:
                    break
            if records is None or len(records) == 0:
                return None

        # Collect face-flow records
        rows: list[Tuple[int, int, int, float]] = []
        for rec in records:
            if rec.ndim == 1:
                for i, q in enumerate(rec):
                    if abs(q) > 1e-15:
                        rows.append((0, i + 1, 0, float(q)))
            elif rec.ndim == 2:
                nconn = rec.shape[1]
                for i in range(rec.shape[0]):
                    for j in range(nconn):
                        q = float(rec[i, j])
                        if abs(q) > 1e-15:
                            rows.append((i + 1, j + 1, 0, q))

        if not rows:
            return None

        dtype = np.dtype([
            ("cell_from", np.int32),
            ("cell_to", np.int32),
            ("face", np.int32),
            ("flow_rate", np.float64),
        ])
        return np.array(rows, dtype=dtype)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Grid-geometry construction from node coordinates
# ---------------------------------------------------------------------------


def build_grid_geometry_from_params(
    ncol: int,
    nrow: int,
    nlay: int = 1,
    dx: float = 0.0,
    dy: float = 0.0,
    rotation_deg: float = 0.0,
    origin_x: float = 0.0,
    origin_y: float = 0.0,
) -> GridGeometry:
    """Build a :class:`GridGeometry` from explicit grid parameters.

    This is the simplest constructor for cases where the MODFLOW grid
    structure is known (e.g. from the .lst listing file).

    .. deprecated::
        This is a trivial wrapper around :class:`GridGeometry`.  Prefer
        constructing ``GridGeometry(...)`` directly.

    Parameters
    ----------
    ncol, nrow, nlay:
        Grid dimensions.
    dx, dy:
        Uniform cell width / height in model length units.
    rotation_deg:
        ANGROT grid rotation (degrees CCW from projected axes).
    origin_x, origin_y:
        Grid origin in projected coordinates.

    Returns
    -------
    Populated :class:`GridGeometry`.
    """
    return GridGeometry(
        ncol=ncol,
        nrow=nrow,
        nlay=nlay,
        dx=dx,
        dy=dy,
        rotation_deg=rotation_deg,
        origin_x=origin_x,
        origin_y=origin_y,
    )


def build_grid_geometry_from_nodes(
    node_df: pd.DataFrame,
    ncol: int,
    nrow: int,
    nlay: int = 1,
    dx: float = 0.0,
    dy: float = 0.0,
    rotation_deg: float = 0.0,
    origin_x: float = 0.0,
    origin_y: float = 0.0,
) -> GridGeometry:
    """Build :class:`GridGeometry`, optionally inferring spacing from node positions.

    When *dx* or *dy* are 0, the function attempts to estimate uniform cell
    size from the median spacing of node x / y coordinates in *node_df*.

    Parameters
    ----------
    node_df:
        DataFrame with ``x``, ``y`` columns in model coordinates.
    ncol, nrow, nlay:
        Grid dimensions.
    dx, dy:
        Cell size; 0 triggers auto-detection from *node_df*.
    rotation_deg, origin_x, origin_y:
        As in :func:`build_grid_geometry_from_params`.

    Returns
    -------
    Populated :class:`GridGeometry`.
    """
    _dx = dx
    _dy = dy

    if _dx <= 0 and "x" in node_df.columns:
        x_vals = node_df["x"].dropna().sort_values().values
        if len(x_vals) >= 2:
            diffs = np.diff(np.unique(x_vals))
            diffs = diffs[diffs > 1e-6]
            if len(diffs) > 0:
                _dx = float(np.median(diffs))

    if _dy <= 0 and "y" in node_df.columns:
        y_vals = node_df["y"].dropna().sort_values().values
        if len(y_vals) >= 2:
            diffs = np.diff(np.unique(y_vals))
            diffs = diffs[diffs > 1e-6]
            if len(diffs) > 0:
                _dy = float(np.median(diffs))

    return GridGeometry(
        ncol=ncol,
        nrow=nrow,
        nlay=nlay,
        dx=_dx,
        dy=_dy,
        rotation_deg=rotation_deg,
        origin_x=origin_x,
        origin_y=origin_y,
    )
