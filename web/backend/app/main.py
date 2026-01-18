"""
Hydrosheaf Web Application - FastAPI Backend
Provides REST API endpoints for the Hydrosheaf hydrogeochemistry framework.
"""

import os
import asyncio
from datetime import datetime, timedelta
from contextlib import asynccontextmanager

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse
from pathlib import Path

try:
    from slowapi import Limiter, _rate_limit_exceeded_handler
    from slowapi.util import get_remote_address
    from slowapi.errors import RateLimitExceeded
    RATE_LIMITING_AVAILABLE = True
except ImportError:
    RATE_LIMITING_AVAILABLE = False

from .routers import analysis, network, samples, projects, demo
from .routers.samples import init_demo_data
from .database import init_db
from .websocket_manager import manager as ws_manager

# Configuration constants
MAX_JOB_AGE_HOURS = int(os.getenv("MAX_JOB_AGE_HOURS", "24"))
CLEANUP_INTERVAL_SECONDS = int(os.getenv("CLEANUP_INTERVAL_SECONDS", "3600"))

# CORS origins from environment variable or default to dev origins
CORS_ORIGINS = os.getenv(
    "CORS_ORIGINS",
    "http://localhost:5173,http://127.0.0.1:5173,http://localhost:5174,http://127.0.0.1:5174"
).split(",")


async def cleanup_old_jobs():
    """Background task to remove old jobs and networks from memory."""
    while True:
        await asyncio.sleep(CLEANUP_INTERVAL_SECONDS)
        try:
            cutoff = datetime.utcnow() - timedelta(hours=MAX_JOB_AGE_HOURS)

            # Cleanup analysis jobs
            jobs_to_remove = []
            for job_id, job in analysis.analysis_jobs.items():
                created_at = job.get('created_at')
                if created_at:
                    try:
                        job_time = datetime.fromisoformat(created_at.replace('Z', '+00:00').replace('+00:00', ''))
                        if job_time < cutoff:
                            jobs_to_remove.append(job_id)
                    except (ValueError, TypeError):
                        pass

            for job_id in jobs_to_remove:
                del analysis.analysis_jobs[job_id]

            # Cleanup old networks (keep max 100)
            if len(network.networks) > 100:
                # Remove oldest networks (simple FIFO)
                network_ids = list(network.networks.keys())
                for net_id in network_ids[:-100]:
                    del network.networks[net_id]

            if jobs_to_remove:
                print(f"[Cleanup] Removed {len(jobs_to_remove)} old analysis jobs")

        except Exception as e:
            print(f"[Cleanup] Error during cleanup: {e}")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Manage application lifespan - start background tasks on startup."""
    # Initialize database
    init_db()
    print("[Startup] SQLite database initialized")

    # Initialize demo data
    init_demo_data()
    print("[Startup] Demo data initialized")

    # Start cleanup task
    cleanup_task = asyncio.create_task(cleanup_old_jobs())
    print(f"[Startup] Job cleanup task started (interval: {CLEANUP_INTERVAL_SECONDS}s, max age: {MAX_JOB_AGE_HOURS}h)")
    yield
    # Cancel cleanup task on shutdown
    cleanup_task.cancel()
    try:
        await cleanup_task
    except asyncio.CancelledError:
        pass


app = FastAPI(
    title="Hydrosheaf API",
    description="REST API for Hydrosheaf Groundwater Hydrogeochemistry Framework",
    version="1.0.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc",
    lifespan=lifespan,
)

# Rate limiting setup
if RATE_LIMITING_AVAILABLE:
    limiter = Limiter(key_func=get_remote_address)
    app.state.limiter = limiter
    app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)
    print("[Startup] Rate limiting enabled")
else:
    limiter = None
    print("[Startup] Rate limiting disabled (slowapi not installed)")

# CORS middleware for frontend development
app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(analysis.router, prefix="/api/analysis", tags=["Analysis"])
app.include_router(network.router, prefix="/api/network", tags=["Network"])
app.include_router(samples.router, prefix="/api/samples", tags=["Samples"])
app.include_router(projects.router, prefix="/api/projects", tags=["Projects"])
app.include_router(demo.router, prefix="/api/demo", tags=["Demo Objectives"])


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


@app.websocket("/ws/analysis/{job_id}")
async def websocket_analysis_progress(websocket: WebSocket, job_id: str):
    """
    WebSocket endpoint for real-time analysis progress updates.

    Connect to ws://localhost:8000/ws/analysis/{job_id} to receive
    progress updates for a specific analysis job.

    Messages are JSON with format:
    - Progress: {"type": "progress", "job_id": "...", "progress": 50, "step": "...", "status": "running"}
    - Completion: {"type": "completion", "job_id": "...", "status": "completed", "progress": 100}
    """
    await ws_manager.connect(websocket, job_id)
    try:
        # Keep connection alive, wait for client messages or disconnection
        while True:
            # Wait for any message from client (ping/pong or close)
            try:
                data = await websocket.receive_text()
                # Client can send "ping" to keep connection alive
                if data == "ping":
                    await websocket.send_text('{"type": "pong"}')
            except WebSocketDisconnect:
                break
    finally:
        await ws_manager.disconnect(websocket, job_id)
