"""
Project storage utilities (filesystem-backed JSON).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import threading
import uuid


PROJECTS_DIR = Path(__file__).resolve().parent.parent / "project_storage"
_LOCK = threading.Lock()
_PROJECTS: Dict[str, Dict[str, Any]] = {}


def _now_iso() -> str:
    return datetime.utcnow().isoformat()


def _project_path(project_id: str) -> Path:
    return PROJECTS_DIR / f"{project_id}.json"


def _normalize_project(project: Dict[str, Any]) -> Dict[str, Any]:
    project.setdefault("analysis_jobs", [])
    project.setdefault("datasets", [])
    project.setdefault("results", [])
    project.setdefault("inputs", [])
    return project


def load_projects() -> None:
    PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
    for path in PROJECTS_DIR.glob("*.json"):
        try:
            with path.open("r", encoding="utf-8") as handle:
                data = json.load(handle)
            if "id" not in data:
                data["id"] = path.stem
            _PROJECTS[data["id"]] = _normalize_project(data)
        except Exception:
            # Skip invalid/corrupt files rather than crashing startup.
            continue


def list_projects() -> List[Dict[str, Any]]:
    return list(_PROJECTS.values())


def get_project(project_id: str) -> Optional[Dict[str, Any]]:
    return _PROJECTS.get(project_id)


def create_project(name: str, description: Optional[str]) -> Dict[str, Any]:
    with _LOCK:
        project_id = str(uuid.uuid4())[:12]
        now = _now_iso()
        project = _normalize_project(
            {
                "id": project_id,
                "name": name,
                "description": description,
                "created_at": now,
                "updated_at": now,
            }
        )
        _PROJECTS[project_id] = project
        _write_project(project)
        return project


def update_project(project: Dict[str, Any]) -> None:
    with _LOCK:
        project["updated_at"] = _now_iso()
        _PROJECTS[project["id"]] = _normalize_project(project)
        _write_project(project)


def delete_project(project_id: str) -> None:
    with _LOCK:
        if project_id in _PROJECTS:
            del _PROJECTS[project_id]
        path = _project_path(project_id)
        if path.exists():
            path.unlink()


def add_analysis_job(project_id: str, job_id: str) -> Optional[Dict[str, Any]]:
    project = get_project(project_id)
    if not project:
        return None
    if job_id not in project["analysis_jobs"]:
        project["analysis_jobs"].append(job_id)
    update_project(project)
    return project


def add_result(project_id: str, result: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    project = get_project(project_id)
    if not project:
        return None
    result = dict(result)
    result["saved_at"] = _now_iso()
    project["results"].append(result)
    update_project(project)
    return project


def add_input(project_id: str, input_payload: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    project = get_project(project_id)
    if not project:
        return None
    payload = dict(input_payload)
    payload["saved_at"] = _now_iso()
    project["inputs"].append(payload)
    dataset_id = payload.get("dataset_id")
    if dataset_id and dataset_id not in project["datasets"]:
        project["datasets"].append(dataset_id)
    update_project(project)
    return project


def _write_project(project: Dict[str, Any]) -> None:
    PROJECTS_DIR.mkdir(parents=True, exist_ok=True)
    path = _project_path(project["id"])
    with path.open("w", encoding="utf-8") as handle:
        json.dump(project, handle, indent=2, sort_keys=True)


load_projects()
