from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from app.api.auth import router as auth_router
from app.api.jobs import router as jobs_router
from app.core.config import settings
from app.core.db import init_db
from app.services.job_runner import job_runner


@asynccontextmanager
async def lifespan(app: FastAPI):
    init_db()
    job_runner.start()
    yield
    job_runner.stop()


app = FastAPI(title="AlphaGEM Web Backend", version="0.1.0", lifespan=lifespan)
app.include_router(auth_router, prefix="/api/auth")
app.include_router(jobs_router, prefix="/api")
frontend_dir = settings.repo_root / "web-frontend"
assets_dir = frontend_dir / "assets"
if assets_dir.exists():
    app.mount("/assets", StaticFiles(directory=assets_dir), name="assets")


@app.get("/healthz")
def healthz() -> dict[str, str]:
    return {"status": "ok", "repoRoot": str(settings.repo_root)}


@app.get("/")
def landing() -> FileResponse:
    return FileResponse(frontend_dir / "index.html")


@app.get("/en")
def landing_en() -> FileResponse:
    return FileResponse(frontend_dir / "index-en.html")


@app.get("/en/")
def landing_en_slash() -> FileResponse:
    return FileResponse(frontend_dir / "index-en.html")


@app.get("/AlphaGEM")
def alphagem_app() -> FileResponse:
    return FileResponse(frontend_dir / "alphagem.html")


@app.get("/AlphaGEM/")
def alphagem_app_slash() -> FileResponse:
    return FileResponse(frontend_dir / "alphagem.html")


@app.get("/en/AlphaGEM")
def alphagem_app_en() -> FileResponse:
    return FileResponse(frontend_dir / "alphagem-en.html")


@app.get("/en/AlphaGEM/")
def alphagem_app_en_slash() -> FileResponse:
    return FileResponse(frontend_dir / "alphagem-en.html")
