"""
Projects Router - Handles project management and result export
"""

from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
from datetime import datetime
import uuid
import json
import io

router = APIRouter()

# In-memory storage for projects
projects: Dict[str, Dict] = {}


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


@router.post("/create")
async def create_project(project: ProjectCreate):
    """Create a new project"""
    project_id = str(uuid.uuid4())[:12]
    
    now = datetime.utcnow().isoformat()
    projects[project_id] = {
        "id": project_id,
        "name": project.name,
        "description": project.description,
        "created_at": now,
        "updated_at": now,
        "analysis_jobs": [],
        "datasets": [],
        "results": [],
    }
    
    return {
        "project_id": project_id,
        "name": project.name,
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
        for p in projects.values()
    ]


@router.get("/{project_id}")
async def get_project(project_id: str):
    """Get project details"""
    if project_id not in projects:
        raise HTTPException(status_code=404, detail="Project not found")
    
    return projects[project_id]


@router.post("/{project_id}/add-analysis")
async def add_analysis_to_project(project_id: str, job_id: str):
    """Add an analysis job to a project"""
    if project_id not in projects:
        raise HTTPException(status_code=404, detail="Project not found")
    
    projects[project_id]["analysis_jobs"].append(job_id)
    projects[project_id]["updated_at"] = datetime.utcnow().isoformat()
    
    return {"message": "Analysis added to project"}


@router.post("/{project_id}/add-result")
async def add_result_to_project(project_id: str, result: Dict[str, Any]):
    """Add a result to a project"""
    if project_id not in projects:
        raise HTTPException(status_code=404, detail="Project not found")
    
    result["saved_at"] = datetime.utcnow().isoformat()
    projects[project_id]["results"].append(result)
    projects[project_id]["updated_at"] = datetime.utcnow().isoformat()
    
    return {"message": "Result saved to project"}


@router.get("/{project_id}/export")
async def export_project(project_id: str):
    """Export project results as JSON (can be converted to PDF on frontend)"""
    if project_id not in projects:
        raise HTTPException(status_code=404, detail="Project not found")
    
    project = projects[project_id]
    
    export_data = {
        "project_name": project["name"],
        "description": project["description"],
        "created_at": project["created_at"],
        "exported_at": datetime.utcnow().isoformat(),
        "framework": "Hydrosheaf - Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry",
        "authors": [
            "Dickson Abdul-Wahab",
            "Ebenezer Aquisman Asare", 
            "Abdul Rashid Dickson"
        ],
        "analysis_summary": {
            "total_analyses": len(project["analysis_jobs"]),
            "total_results": len(project.get("results", [])),
        },
        "results": project.get("results", []),
    }
    
    return export_data


@router.get("/{project_id}/download")
async def download_project_report(project_id: str):
    """Download project report as a text document"""
    if project_id not in projects:
        raise HTTPException(status_code=404, detail="Project not found")
    
    project = projects[project_id]
    
    # Generate a formatted text report
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
        f"Total Analyses: {len(project['analysis_jobs'])}",
        f"Total Results Saved: {len(project.get('results', []))}",
        "",
    ]
    
    # Add results
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
            
            if result.get('transport_model'):
                report_lines.append("\nTransport Model:")
                tm = result['transport_model']
                report_lines.append(f"  Dominant Process: {tm.get('dominant_process', 'N/A')}")
                report_lines.append(f"  Evaporation Fraction: {tm.get('evaporation_fraction', 'N/A')}")
            
            if result.get('reactions'):
                report_lines.append("\nReactions:")
                for rxn in result['reactions']:
                    report_lines.append(f"  - {rxn.get('mineral', 'Unknown')}: {rxn.get('rate_mmol_L', 'N/A')} mmol/L ({rxn.get('direction', 'N/A')})")
    else:
        report_lines.append("No results saved in this project yet.")
    
    report_lines.append("")
    report_lines.append("=" * 60)
    report_lines.append("END OF REPORT")
    report_lines.append("=" * 60)
    
    report_content = "\n".join(report_lines)
    
    # Return as downloadable file
    return StreamingResponse(
        io.BytesIO(report_content.encode('utf-8')),
        media_type="text/plain",
        headers={
            "Content-Disposition": f"attachment; filename={project['name'].replace(' ', '_')}_report.txt"
        }
    )


@router.delete("/{project_id}")
async def delete_project(project_id: str):
    """Delete a project"""
    if project_id not in projects:
        raise HTTPException(status_code=404, detail="Project not found")
    
    del projects[project_id]
    return {"message": "Project deleted successfully"}
