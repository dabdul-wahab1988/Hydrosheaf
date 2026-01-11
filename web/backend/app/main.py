"""
Hydrosheaf Web Application - FastAPI Backend
Provides REST API endpoints for the Hydrosheaf hydrogeochemistry framework.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pathlib import Path

from .routers import analysis, network, samples, projects

app = FastAPI(
    title="Hydrosheaf API",
    description="REST API for Hydrosheaf Groundwater Hydrogeochemistry Framework",
    version="1.0.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc",
)

# CORS middleware for frontend development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:5174", "http://127.0.0.1:5174"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(analysis.router, prefix="/api/analysis", tags=["Analysis"])
app.include_router(network.router, prefix="/api/network", tags=["Network"])
app.include_router(samples.router, prefix="/api/samples", tags=["Samples"])
app.include_router(projects.router, prefix="/api/projects", tags=["Projects"])


@app.get("/api/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "service": "Hydrosheaf API",
        "version": "1.0.0",
    }


@app.get("/api/")
async def root():
    """Root API endpoint with welcome message"""
    return {
        "message": "Welcome to Hydrosheaf API",
        "description": "Sheaf-Theoretic Methods in Groundwater Hydrogeochemistry",
        "documentation": "/api/docs",
    }
