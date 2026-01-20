"""
Projects Router - Handles project management and result export
"""

from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
from datetime import datetime
import json
import io
import csv
import re
import zipfile

from .. import project_store

router = APIRouter()


class ProjectCreate(BaseModel):
    """Create a new project"""

    name: str = Field(..., min_length=1, max_length=100)
    description: Optional[str] = None


class Project(BaseModel):
    """Project model"""

    id: str
    name: str
    description: Optional[str]
    created_at: str
    updated_at: str
    analysis_jobs: List[str] = []
    datasets: List[str] = []
    results: List[Dict[str, Any]] = []
    inputs: List[Dict[str, Any]] = []


def _safe_name(value: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned or fallback


def _csv_from_records(records: List[Dict[str, Any]]) -> str:
    if not records:
        return ""
    fieldnames = sorted({key for record in records for key in record.keys()})
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    for record in records:
        writer.writerow({key: record.get(key) for key in fieldnames})
    return output.getvalue()


def _build_text_report(project: Dict[str, Any]) -> str:
    report_lines = [
        "=" * 60,
        "HYDROSHEAF PROJECT REPORT",
        "=" * 60,
        "",
        f"Project Name: {project['name']}",
        f"Description: {project.get('description', 'N/A')}",
        f"Created: {project['created_at']}",
        f"Exported: {datetime.utcnow().isoformat()}",
        "",
        "-" * 60,
        "FRAMEWORK INFORMATION",
        "-" * 60,
        "Hydrosheaf: Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry",
        "",
        "Authors:",
        "  - Dickson Abdul-Wahab",
        "  - Ebenezer Aquisman Asare",
        "  - Abdul Rashid Dickson",
        "",
        "-" * 60,
        "ANALYSIS SUMMARY",
        "-" * 60,
        f"Total Analyses: {len(project.get('analysis_jobs', []))}",
        f"Total Results Saved: {len(project.get('results', []))}",
        "",
    ]

    results = project.get("results", [])
    if results:
        report_lines.append("-" * 60)
        report_lines.append("DETAILED RESULTS")
        report_lines.append("-" * 60)

        for i, result in enumerate(results, 1):
            report_lines.append(f"\n--- Result {i} ---")
            report_lines.append(f"Name: {result.get('name', 'Unnamed')}")
            report_lines.append(f"Type: {result.get('analysis_type', 'Unknown')}")
            report_lines.append(f"Saved: {result.get('saved_at', 'Unknown')}")

            if result.get("transport_model"):
                report_lines.append("\nTransport Model:")
                tm = result["transport_model"]
                report_lines.append(
                    f"  Dominant Process: {tm.get('dominant_process', 'N/A')}"
                )
                report_lines.append(
                    f"  Evaporation Fraction: {tm.get('evaporation_fraction', 'N/A')}"
                )

            if result.get("reactions"):
                report_lines.append("\nReactions:")
                for rxn in result["reactions"]:
                    report_lines.append(
                        f"  - {rxn.get('mineral', 'Unknown')}: "
                        f"{rxn.get('rate_mmol_L', 'N/A')} mmol/L ({rxn.get('direction', 'N/A')})"
                    )
    else:
        report_lines.append("No results saved in this project yet.")

    report_lines.append("")
    report_lines.append("=" * 60)
    report_lines.append("END OF REPORT")
    report_lines.append("=" * 60)

    return "\n".join(report_lines)


@router.post("/create")
async def create_project(project: ProjectCreate):
    """Create a new project"""
    created = project_store.create_project(project.name, project.description)
    return {
        "project_id": created["id"],
        "name": created["name"],
        "message": "Project created successfully. All your analysis results will be saved to this project.",
    }


@router.get("/list")
async def list_projects():
    """List all projects"""
    return [
        {
            "id": p["id"],
            "name": p["name"],
            "description": p["description"],
            "created_at": p["created_at"],
            "analysis_count": len(p["analysis_jobs"]),
            "has_results": len(p.get("results", [])) > 0,
        }
        for p in project_store.list_projects()
    ]


@router.get("/{project_id}")
async def get_project(project_id: str):
    """Get project details"""
    project = project_store.get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    return project


@router.post("/{project_id}/add-analysis")
async def add_analysis_to_project(project_id: str, job_id: str):
    """Add an analysis job to a project"""
    project = project_store.add_analysis_job(project_id, job_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    return {"message": "Analysis added to project"}


@router.post("/{project_id}/add-result")
async def add_result_to_project(project_id: str, result: Dict[str, Any]):
    """Add a result to a project"""
    project = project_store.add_result(project_id, result)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    return {"message": "Result saved to project"}


@router.get("/{project_id}/export")
async def export_project(project_id: str):
    """Export project results as JSON (can be converted to PDF on frontend)"""
    project = project_store.get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    export_data = {
        "project_name": project["name"],
        "description": project["description"],
        "created_at": project["created_at"],
        "exported_at": datetime.utcnow().isoformat(),
        "framework": "Hydrosheaf - Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry",
        "authors": [
            "Dickson Abdul-Wahab",
            "Ebenezer Aquisman Asare",
            "Abdul Rashid Dickson",
        ],
        "analysis_summary": {
            "total_analyses": len(project["analysis_jobs"]),
            "total_results": len(project.get("results", [])),
        },
        "results": project.get("results", []),
        "inputs": project.get("inputs", []),
    }

    return export_data


@router.get("/{project_id}/download")
async def download_project_report(project_id: str):
    """Download project report as a text document"""
    project = project_store.get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    report_content = _build_text_report(project)

    # Return as downloadable file
    return StreamingResponse(
        io.BytesIO(report_content.encode("utf-8")),
        media_type="text/plain",
        headers={
            "Content-Disposition": f"attachment; filename={project['name'].replace(' ', '_')}_report.txt"
        },
    )


@router.get("/{project_id}/download-complete")
async def download_complete_project(project_id: str):
    """Download complete project export as a ZIP archive"""
    project = project_store.get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    safe_name = _safe_name(project.get("name", ""), f"project_{project_id}")
    exported_at = datetime.utcnow().isoformat()

    metadata = {
        "project_id": project["id"],
        "project_name": project["name"],
        "description": project.get("description"),
        "created_at": project.get("created_at"),
        "updated_at": project.get("updated_at"),
        "exported_at": exported_at,
        "framework": "Hydrosheaf - Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry",
        "authors": [
            "Dickson Abdul-Wahab",
            "Ebenezer Aquisman Asare",
            "Abdul Rashid Dickson",
        ],
        "analysis_summary": {
            "total_analyses": len(project.get("analysis_jobs", [])),
            "total_results": len(project.get("results", [])),
            "total_inputs": len(project.get("inputs", [])),
        },
    }

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        root = safe_name
        zf.writestr(f"{root}/README.txt", "Hydrosheaf project export bundle.\n")
        zf.writestr(f"{root}/project_metadata.json", json.dumps(metadata, indent=2))

        inputs = project.get("inputs", [])
        zf.writestr(f"{root}/data/inputs.json", json.dumps(inputs, indent=2))
        for idx, payload in enumerate(inputs, 1):
            run_dir = f"{root}/data/run_{idx}"
            samples = payload.get("samples", [])
            edges = payload.get("edges", [])
            config = payload.get("config", {})
            run_meta = {
                "analysis_name": payload.get("analysis_name"),
                "analysis_type": payload.get("analysis_type"),
                "analysis_job_id": payload.get("analysis_job_id"),
                "dataset_id": payload.get("dataset_id"),
                "dataset_name": payload.get("dataset_name"),
                "dataset_description": payload.get("dataset_description"),
                "edge_source": payload.get("edge_source"),
                "saved_at": payload.get("saved_at"),
            }
            zf.writestr(f"{run_dir}/metadata.json", json.dumps(run_meta, indent=2))
            zf.writestr(f"{run_dir}/samples.json", json.dumps(samples, indent=2))
            csv_samples = _csv_from_records(samples)
            if csv_samples:
                zf.writestr(f"{run_dir}/samples.csv", csv_samples)
            zf.writestr(f"{run_dir}/edges.json", json.dumps(edges, indent=2))
            csv_edges = _csv_from_records(edges)
            if csv_edges:
                zf.writestr(f"{run_dir}/edges.csv", csv_edges)
            zf.writestr(f"{run_dir}/config.json", json.dumps(config, indent=2))

        results = project.get("results", [])
        zf.writestr(f"{root}/results/results.json", json.dumps(results, indent=2))
        result_summaries = []
        for result in results:
            transport = (
                result.get("transport_model", {})
                if isinstance(result.get("transport_model"), dict)
                else {}
            )
            metadata_block = (
                result.get("metadata", {})
                if isinstance(result.get("metadata"), dict)
                else {}
            )
            result_summaries.append(
                {
                    "name": result.get("name"),
                    "analysis_type": result.get("analysis_type"),
                    "job_id": result.get("job_id"),
                    "saved_at": result.get("saved_at"),
                    "transport_dominant_process": transport.get("dominant_process"),
                    "transport_evaporation_fraction": transport.get(
                        "evaporation_fraction"
                    ),
                    "reaction_count": (
                        len(result.get("reactions", []))
                        if isinstance(result.get("reactions"), list)
                        else None
                    ),
                    "edges_analyzed": metadata_block.get("edges_analyzed"),
                    "samples_analyzed": metadata_block.get("samples_analyzed"),
                }
            )

        csv_results = _csv_from_records(result_summaries)
        if csv_results:
            zf.writestr(f"{root}/results/results_summary.csv", csv_results)

        report_content = _build_text_report(project)
        zf.writestr(f"{root}/reports/summary_report.txt", report_content)

    buffer.seek(0)
    return StreamingResponse(
        buffer,
        media_type="application/zip",
        headers={"Content-Disposition": f"attachment; filename={safe_name}.zip"},
    )


@router.delete("/{project_id}")
async def delete_project(project_id: str):
    """Delete a project"""
    project = project_store.get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    project_store.delete_project(project_id)
    return {"message": "Project deleted successfully"}
