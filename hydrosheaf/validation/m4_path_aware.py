"""Truth-blind path-aware physical evidence for the M4-C benchmark.

M4-C is deliberately a physical routing layer, not a MODPATH imitation.  It
builds a local transition model on the public structured MODFLOW grid from:

* the archived hydraulic-head field;
* the head-gradient direction implied by that field;
* absolute CBC face-flow activity as a path-capacity proxy; and
* public CBC well extraction context as a soft endpoint preference.

The sign of a CBC face record is *not* used as a standalone flow direction.
The previous benchmark audit showed that a naive dominant-face direction can
be anti-aligned with particle trajectories.  Head-gradient direction is
therefore the directional evidence; CBC face magnitudes only describe how
active a local connection is.

The module returns path features for every endpoint pair supplied by the
caller.  It never accepts, imports, or inspects MODPATH reference edges.  A
path not reached by the finite beam is represented by zero support rather
than being removed from the candidate universe.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any

import numpy as np

from hydrosheaf.physics.modflow_head import GridGeometry, compute_head_gradient

M4_C = "M4-C_path_aware"


def _cell_to_ijk(cell_id: int, ncol: int, nrow: int) -> tuple[int, int, int]:
    index = int(cell_id) - 1
    cells_per_layer = int(ncol) * int(nrow)
    layer, remainder = divmod(index, cells_per_layer)
    row, col = divmod(remainder, int(ncol))
    return col, row, layer


def _ijk_to_cell(col: int, row: int, layer: int, ncol: int, nrow: int) -> int:
    return int(layer) * int(nrow) * int(ncol) + int(row) * int(ncol) + int(col) + 1


def _cell_from_node_id(node_id: object) -> int:
    text = str(node_id)
    if "_" not in text:
        raise ValueError(f"M4-C node IDs must end in a MODFLOW cell number: {text!r}")
    try:
        cell = int(text.rsplit("_", 1)[1])
    except ValueError as exc:
        raise ValueError(f"M4-C node ID has no integer cell suffix: {text!r}") from exc
    if cell <= 0:
        raise ValueError(f"M4-C cell IDs must be positive: {text!r}")
    return cell


def _face_array(budget_file: Any, record_name: str, shape: tuple[int, int, int]) -> np.ndarray:
    try:
        blocks = budget_file.get_data(text=record_name, full3D=False)
    except Exception:
        return np.zeros(shape, dtype=np.float64)
    for block in reversed(blocks or ()):
        array = np.asarray(block, dtype=np.float64)
        if array.ndim == 3:
            if array.shape != shape:
                raise ValueError(
                    f"CBC {record_name!r} has shape {array.shape}; expected {shape}."
                )
            return np.nan_to_num(array, nan=0.0, posinf=0.0, neginf=0.0)
    return np.zeros(shape, dtype=np.float64)


def aggregate_savage_cbc_face_activity(
    cbc_path: str | Path,
    *,
    nlay: int,
    nrow: int,
    ncol: int,
) -> tuple[dict[int, float], dict[str, object]]:
    """Aggregate absolute public CBC face-flow activity by MODFLOW cell.

    Each face magnitude is assigned to both cells incident to that face.  The
    result is a non-directional activity measure.  Direction is supplied by
    :func:`compute_head_gradient` in the transition model below.
    """

    from flopy.utils import CellBudgetFile

    path = Path(cbc_path)
    if not path.exists():
        raise FileNotFoundError(path)
    shape = (int(nlay), int(nrow), int(ncol))
    budget_file = CellBudgetFile(str(path))
    right = _face_array(budget_file, "FLOW RIGHT FACE", shape)
    front = _face_array(budget_file, "FLOW FRONT FACE", shape)
    lower = _face_array(budget_file, "FLOW LOWER FACE", shape)

    activity = np.abs(right) + np.abs(front) + np.abs(lower)
    if ncol > 1:
        activity[:, :, 1:] += np.abs(right[:, :, :-1])
    if nrow > 1:
        activity[:, 1:, :] += np.abs(front[:, :-1, :])
    if nlay > 1:
        activity[1:, :, :] += np.abs(lower[:-1, :, :])

    values = activity[np.isfinite(activity) & (activity > 0.0)]
    scale = float(np.median(values)) if values.size else 1.0
    scale = max(scale, 1e-12)
    by_cell: dict[int, float] = {}
    for layer in range(int(nlay)):
        for row in range(int(nrow)):
            for col in range(int(ncol)):
                cell = _ijk_to_cell(col, row, layer, ncol, nrow)
                by_cell[cell] = float(activity[layer, row, col])
    metadata = {
        "cbc_path": str(path),
        "record_names_used": [
            "FLOW RIGHT FACE",
            "FLOW FRONT FACE",
            "FLOW LOWER FACE",
        ],
        "face_activity_definition": (
            "sum of absolute face-flow magnitudes, assigning each face to both incident cells"
        ),
        "activity_scale_median_nonzero": scale,
        "direction_from_cbc_sign": False,
        "truth_blind": True,
    }
    return by_cell, metadata


@dataclass(frozen=True)
class PathAwareConfig:
    """Predeclared M4-C routing parameters.

    These values are physical-routing hyperparameters, not Savage-fitted
    thresholds.  The endpoint decision threshold is deliberately kept in the
    benchmark runner so ranking metrics can be reported independently.
    """

    head_smoothing_sigma_cells: float = 1.0
    transition_temperature: float = 0.75
    gradient_alignment_weight: float = 2.0
    head_drop_weight: float = 1.5
    face_activity_weight: float = 0.25
    pumping_sink_weight: float = 0.35
    vertical_transition_penalty: float = 0.50
    path_decay_per_step: float = 0.998
    beam_width: int = 32
    branch_width: int = 3
    max_steps: int = 600
    min_log_support: float = -700.0
    vertical_layer_distance_proxy: float = 1.0

    def __post_init__(self) -> None:
        if self.transition_temperature <= 0.0:
            raise ValueError("M4-C transition_temperature must be positive.")
        if not 0.0 < self.path_decay_per_step <= 1.0:
            raise ValueError("M4-C path_decay_per_step must be in (0, 1].")
        if self.beam_width < 1 or self.branch_width < 1 or self.max_steps < 1:
            raise ValueError("M4-C beam_width, branch_width, and max_steps must be positive.")
        if self.vertical_layer_distance_proxy <= 0.0:
            raise ValueError("M4-C vertical_layer_distance_proxy must be positive.")

    def to_dict(self) -> dict[str, object]:
        return {
            "head_smoothing_sigma_cells": self.head_smoothing_sigma_cells,
            "transition_temperature": self.transition_temperature,
            "gradient_alignment_weight": self.gradient_alignment_weight,
            "head_drop_weight": self.head_drop_weight,
            "face_activity_weight": self.face_activity_weight,
            "pumping_sink_weight": self.pumping_sink_weight,
            "vertical_transition_penalty": self.vertical_transition_penalty,
            "path_decay_per_step": self.path_decay_per_step,
            "beam_width": self.beam_width,
            "branch_width": self.branch_width,
            "max_steps": self.max_steps,
            "min_log_support": self.min_log_support,
            "vertical_layer_distance_proxy": self.vertical_layer_distance_proxy,
        }


@dataclass(frozen=True)
class PathAwareResult:
    """Truth-blind path features and routing audit information."""

    edge_features: Mapping[str, Mapping[str, object]]
    source_diagnostics: tuple[Mapping[str, object], ...]
    head_gradient_count: int
    head_drop_scale: float
    face_activity_metadata: Mapping[str, object]
    config: PathAwareConfig

    def to_dict(self) -> dict[str, object]:
        return {
            "n_edge_features": len(self.edge_features),
            "n_source_diagnostics": len(self.source_diagnostics),
            "head_gradient_count": self.head_gradient_count,
            "head_drop_scale": self.head_drop_scale,
            "face_activity_metadata": dict(self.face_activity_metadata),
            "config": self.config.to_dict(),
            "truth_blind": True,
        }


@dataclass(frozen=True)
class _Transition:
    cell: int
    log_probability: float
    activity_support: float
    vertical: bool


@dataclass(frozen=True)
class _PathState:
    cell: int
    log_probability: float
    path: tuple[int, ...]
    activity_sum: float


class _GridTransitionModel:
    def __init__(
        self,
        *,
        head_map: Mapping[int, float],
        gradient_map: Mapping[int, tuple[float, float, float]],
        grid: GridGeometry,
        activity_by_cell: Mapping[int, float],
        activity_scale: float,
        sink_by_cell: Mapping[int, float],
        sink_scale: float,
        config: PathAwareConfig,
    ) -> None:
        self.heads = {int(key): float(value) for key, value in head_map.items()}
        self.gradients = dict(gradient_map)
        self.grid = grid
        self.activity = {int(key): float(value) for key, value in activity_by_cell.items()}
        self.activity_scale = max(float(activity_scale), 1e-12)
        self.sinks = {int(key): max(0.0, float(value)) for key, value in sink_by_cell.items()}
        self.sink_scale = max(float(sink_scale), 1e-12)
        self.config = config
        self.head_drop_scale = self._estimate_head_drop_scale()
        self._cache: dict[int, tuple[_Transition, ...]] = {}

    def _estimate_head_drop_scale(self) -> float:
        drops: list[float] = []
        for cell in self.heads:
            for neighbour, _vertical, _sign in self._candidate_cells(cell):
                if neighbour in self.heads:
                    drops.append(abs(self.heads[cell] - self.heads[neighbour]))
        if not drops:
            return 1.0
        return max(float(np.median(np.asarray(drops, dtype=float))), 1e-6)

    def _candidate_cells(self, cell: int) -> tuple[tuple[int, bool, int], ...]:
        col, row, layer = _cell_to_ijk(cell, self.grid.ncol, self.grid.nrow)
        candidates: list[tuple[int, bool, int]] = []
        if col > 0:
            candidates.append((_ijk_to_cell(col - 1, row, layer, self.grid.ncol, self.grid.nrow), False, -1))
        if col + 1 < self.grid.ncol:
            candidates.append((_ijk_to_cell(col + 1, row, layer, self.grid.ncol, self.grid.nrow), False, 1))
        if row > 0:
            candidates.append((_ijk_to_cell(col, row - 1, layer, self.grid.ncol, self.grid.nrow), False, -1))
        if row + 1 < self.grid.nrow:
            candidates.append((_ijk_to_cell(col, row + 1, layer, self.grid.ncol, self.grid.nrow), False, 1))
        if layer > 0:
            candidates.append((_ijk_to_cell(col, row, layer - 1, self.grid.ncol, self.grid.nrow), True, -1))
        if layer + 1 < self.grid.nlay:
            candidates.append((_ijk_to_cell(col, row, layer + 1, self.grid.ncol, self.grid.nrow), True, 1))
        return tuple(candidates)

    def _horizontal_unit_vector(self, axis: str, sign: int) -> tuple[float, float]:
        theta = math.radians(float(self.grid.rotation_deg))
        cosine = math.cos(theta)
        sine = math.sin(theta)
        if axis == "x":
            return sign * cosine, sign * sine
        return -sign * sine, sign * cosine

    def _score_transition(self, cell: int, neighbour: int, vertical: bool, sign: int) -> tuple[float, float]:
        current_head = self.heads[cell]
        neighbour_head = self.heads[neighbour]
        drop_score = math.tanh((current_head - neighbour_head) / self.head_drop_scale)
        gx, gy, gz = self.gradients.get(cell, (0.0, 0.0, 0.0))
        if vertical:
            gradient_norm = abs(float(gz))
            alignment = (-(float(gz)) * float(sign) / gradient_norm) if gradient_norm > 1e-12 else 0.0
        else:
            col_a, row_a, _ = _cell_to_ijk(cell, self.grid.ncol, self.grid.nrow)
            col_b, row_b, _ = _cell_to_ijk(neighbour, self.grid.ncol, self.grid.nrow)
            axis = "x" if col_a != col_b else "y"
            unit_x, unit_y = self._horizontal_unit_vector(axis, sign)
            norm = math.hypot(float(gx), float(gy))
            alignment = (-(float(gx) * unit_x + float(gy) * unit_y) / norm) if norm > 1e-12 else 0.0
        activity = math.log1p(max(0.0, self.activity.get(neighbour, 0.0)) / self.activity_scale)
        activity_support = min(1.0, activity / math.log1p(10.0))
        sink = math.log1p(max(0.0, self.sinks.get(neighbour, 0.0)) / self.sink_scale)
        sink_support = min(1.0, sink / math.log1p(10.0))
        score = (
            self.config.gradient_alignment_weight * alignment
            + self.config.head_drop_weight * drop_score
            + self.config.face_activity_weight * activity_support
            + self.config.pumping_sink_weight * sink_support
            - (self.config.vertical_transition_penalty if vertical else 0.0)
        )
        return float(score), float(activity_support)

    def transitions(self, cell: int) -> tuple[_Transition, ...]:
        cached = self._cache.get(int(cell))
        if cached is not None:
            return cached
        scored: list[tuple[int, float, float, bool]] = []
        for neighbour, vertical, sign in self._candidate_cells(int(cell)):
            if neighbour not in self.heads or neighbour not in self.gradients:
                continue
            score, activity_support = self._score_transition(int(cell), neighbour, vertical, sign)
            scored.append((neighbour, score, activity_support, vertical))
        if not scored:
            self._cache[int(cell)] = ()
            return ()
        scored.sort(key=lambda item: (item[1], -item[0]), reverse=True)
        scored = scored[: int(self.config.branch_width)]
        logits = np.asarray([item[1] for item in scored], dtype=float)
        logits = (logits - float(np.max(logits))) / float(self.config.transition_temperature)
        probabilities = np.exp(logits)
        probabilities /= float(np.sum(probabilities))
        transitions = tuple(
            _Transition(
                cell=item[0],
                log_probability=float(math.log(max(probability, 1e-300))),
                activity_support=item[2],
                vertical=item[3],
            )
            for item, probability in zip(scored, probabilities)
        )
        self._cache[int(cell)] = transitions
        return transitions


def _sink_context(
    observations: Iterable[Mapping[str, object]],
) -> tuple[dict[int, float], float]:
    sink_by_cell: dict[int, float] = {}
    values: list[float] = []
    for row in observations:
        try:
            cell = _cell_from_node_id(row.get("site_id", ""))
        except ValueError:
            continue
        sink = 0.0
        for key in ("well_rate", "river_leakage", "head_boundary_flux"):
            try:
                rate = float(row.get(key, 0.0) or 0.0)
            except (TypeError, ValueError):
                rate = 0.0
            sink += max(0.0, -rate)
        sink_by_cell[cell] = sink
        if sink > 0.0:
            values.append(sink)
    return sink_by_cell, max(float(np.median(np.asarray(values))), 1e-12) if values else 1.0


def _trace_source(
    source_cell: int,
    model: _GridTransitionModel,
) -> tuple[dict[int, _PathState], dict[str, object]]:
    config = model.config
    initial = _PathState(source_cell, 0.0, (source_cell,), 0.0)
    beam: list[_PathState] = [initial]
    best: dict[int, _PathState] = {source_cell: initial}
    for step in range(1, int(config.max_steps) + 1):
        by_cell: dict[int, _PathState] = {}
        for state in beam:
            visited = set(state.path)
            for transition in model.transitions(state.cell):
                if transition.cell in visited:
                    continue
                log_probability = (
                    state.log_probability
                    + transition.log_probability
                    + math.log(float(config.path_decay_per_step))
                )
                if log_probability < float(config.min_log_support):
                    continue
                candidate = _PathState(
                    transition.cell,
                    log_probability,
                    state.path + (transition.cell,),
                    state.activity_sum + transition.activity_support,
                )
                previous = by_cell.get(transition.cell)
                if previous is None or candidate.log_probability > previous.log_probability:
                    by_cell[transition.cell] = candidate
        if not by_cell:
            break
        beam = sorted(
            by_cell.values(),
            key=lambda state: (state.log_probability, -state.cell),
            reverse=True,
        )[: int(config.beam_width)]
        for state in beam:
            previous = best.get(state.cell)
            if previous is None or state.log_probability > previous.log_probability:
                best[state.cell] = state
    reachable = [state for cell, state in best.items() if cell != source_cell]
    max_log = max((state.log_probability for state in reachable), default=float("-inf"))
    diagnostics = {
        "source_cell": source_cell,
        "visited_cells_best_beam": len(best),
        "max_reached_steps": max((len(state.path) - 1 for state in best.values()), default=0),
        "max_log_support": max_log if math.isfinite(max_log) else None,
        "beam_terminated_early": len(beam) == 0,
    }
    return best, diagnostics


def build_path_aware_features(
    observations: Iterable[Mapping[str, object]],
    *,
    head_map: Mapping[int, float],
    cbc_path: str | Path,
    grid: GridGeometry,
    config: PathAwareConfig | None = None,
) -> PathAwareResult:
    """Build M4-C path features for all directed observation pairs.

    The returned mapping contains one record for every directed pair.  The
    input observations, head map, CBC file, and grid are the only inputs used;
    no reference topology is accepted by this API.
    """

    resolved = config or PathAwareConfig()
    rows = [dict(row) for row in observations]
    rows.sort(key=lambda row: str(row.get("site_id", "")))
    node_cells = {
        str(row["site_id"]): _cell_from_node_id(row["site_id"])
        for row in rows
    }
    if len(node_cells) != len(rows):
        raise ValueError("M4-C observations require unique site_id values.")
    missing_cells = [
        (node_id, cell)
        for node_id, cell in node_cells.items()
        if cell not in head_map
    ]
    if missing_cells:
        raise ValueError(f"M4-C head map is missing benchmark cells: {missing_cells[:5]}")

    activity_by_cell, activity_metadata = aggregate_savage_cbc_face_activity(
        cbc_path,
        nlay=grid.nlay,
        nrow=grid.nrow,
        ncol=grid.ncol,
    )
    activity_values = np.asarray(
        [activity_by_cell.get(int(cell), 0.0) for cell in head_map], dtype=float
    )
    nonzero_activity = activity_values[activity_values > 0.0]
    activity_scale = float(np.median(nonzero_activity)) if nonzero_activity.size else 1.0
    sink_by_cell, sink_scale = _sink_context(rows)
    gradients = compute_head_gradient(
        {int(key): float(value) for key, value in head_map.items()},
        grid,
        sigma=float(resolved.head_smoothing_sigma_cells),
    )
    model = _GridTransitionModel(
        head_map=head_map,
        gradient_map=gradients,
        grid=grid,
        activity_by_cell=activity_by_cell,
        activity_scale=activity_scale,
        sink_by_cell=sink_by_cell,
        sink_scale=sink_scale,
        config=resolved,
    )
    source_diagnostics: list[Mapping[str, object]] = []
    edge_features: dict[str, Mapping[str, object]] = {}
    for source_id in sorted(node_cells):
        source_cell = node_cells[source_id]
        best, diagnostics = _trace_source(source_cell, model)
        source_diagnostics.append({"source_id": source_id, **diagnostics})
        reachable = [state for cell, state in best.items() if cell != source_cell]
        max_log = max((state.log_probability for state in reachable), default=float("-inf"))
        for target_id in sorted(node_cells):
            if target_id == source_id:
                continue
            target_cell = node_cells[target_id]
            state = best.get(target_cell)
            if state is None:
                record = {
                    "path_reachability": 0.0,
                    "path_relative_support": 0.0,
                    "path_log_support": float(resolved.min_log_support),
                    "path_steps": 0.0,
                    "path_transition_efficiency": 0.0,
                    "path_activity_support": 0.0,
                    "path_endpoint_sink_support": float(
                        1.0 if model.sinks.get(target_cell, 0.0) > 0.0 else 0.0
                    ),
                    "path_source_cell": source_cell,
                    "path_target_cell": target_cell,
                    "path_reached": False,
                }
            else:
                steps = max(1, len(state.path) - 1)
                relative = (
                    math.exp(max(-700.0, min(0.0, state.log_probability - max_log)))
                    if math.isfinite(max_log)
                    else 0.0
                )
                record = {
                    "path_reachability": 1.0,
                    "path_relative_support": float(relative),
                    "path_log_support": float(state.log_probability),
                    "path_steps": float(steps),
                    "path_transition_efficiency": float(
                        math.exp(max(-700.0, state.log_probability / steps))
                    ),
                    "path_activity_support": float(state.activity_sum / steps),
                    "path_endpoint_sink_support": float(
                        1.0 if model.sinks.get(target_cell, 0.0) > 0.0 else 0.0
                    ),
                    "path_source_cell": source_cell,
                    "path_target_cell": target_cell,
                    "path_reached": True,
                }
            edge_features[f"{source_id}->{target_id}"] = record
    activity_metadata = {
        **dict(activity_metadata),
        "activity_scale_used": activity_scale,
        "pumping_sink_scale_used": sink_scale,
        "n_active_head_cells": len(head_map),
        "n_gradient_cells": len(gradients),
        "vertical_geometry": (
            "unit layer-index proxy; no public physical layer thickness was supplied"
        ),
    }
    return PathAwareResult(
        edge_features=edge_features,
        source_diagnostics=tuple(source_diagnostics),
        head_gradient_count=len(gradients),
        head_drop_scale=model.head_drop_scale,
        face_activity_metadata=activity_metadata,
        config=resolved,
    )


__all__ = [
    "M4_C",
    "PathAwareConfig",
    "PathAwareResult",
    "aggregate_savage_cbc_face_activity",
    "build_path_aware_features",
]
