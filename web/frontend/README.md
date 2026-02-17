# Hydrosheaf Web Frontend

React + Vite frontend for the Hydrosheaf web application.

## Setup

```bash
cd web/frontend
npm install
npm run dev
```

The dev server runs on `http://localhost:5173` and expects the backend at `http://localhost:8000`.
Override the backend URL with `VITE_API_BASE` or `VITE_API_HOST`.

## Build & Lint

```bash
npm run build
npm run lint
npm run preview
```

## Demo Objectives Page

The demo objectives page queries `/api/demo` endpoints for research objectives 1–4.
These endpoints require `hydrosheaf_synthetic_csv/` at the repo root and the backend running.

## Environment Variables

- `VITE_API_BASE`: Full backend API base URL (defaults to `http://localhost:8000/api`).
- `VITE_API_HOST`: Host override for WebSocket connections (defaults to window host).
- `VITE_APP_VERSION`: Optional app version string shown in UI.

## Notes

- The frontend assumes mg/L inputs; conversions happen in the backend adapter.
- The demo objectives data uses synthetic CSVs; see `examples/demo_research_objectives.py`.

