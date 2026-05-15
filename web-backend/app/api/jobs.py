from __future__ import annotations

import shutil
from pathlib import Path

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from fastapi.responses import FileResponse

from app.core.auth import get_request_identity
from app.core.config import settings
from app.core.db import ensure_user
from app.schemas import ArtifactResponse, JobListResponse, JobResponse, LogResponse, UserResponse
from app.services.job_runner import job_runner
from app.services.repository import create_job, get_artifact_for_job, get_job_for_user, list_artifacts, list_jobs_for_user

router = APIRouter()

VALID_MODES = {"structure_alignment", "plmsearch"}
VALID_REFNAMES = {"yeast", "human", "ecoli", "strco", "synechocystis"}
VALID_MEDIA = {"min", "full"}


def _job_to_response(job: dict) -> JobResponse:
    artifacts = [ArtifactResponse.model_validate(item) for item in list_artifacts(job["id"])]
    if job["status"] == "queued":
        user_message = "Your job has been received and is waiting to start."
    elif job["status"] == "running":
        user_message = "Your job is currently running."
    elif job["status"] == "succeeded":
        user_message = "Your job finished successfully. Files are ready to download."
    elif job["status"] == "failed":
        user_message = "This job did not finish successfully."
    else:
        user_message = None
    payload = {
        "id": job["id"],
        "public_name": job["public_name"],
        "mode": job["mode"],
        "refname": job["refname"],
        "growthMedium": job["growth_medium"],
        "status": job["status"],
        "created_at": job["created_at"],
        "started_at": job["started_at"],
        "finished_at": job["finished_at"],
        "tmscore": job["tmscore"],
        "coverage": job["coverage"],
        "plddt": job["plddt"],
        "queue_position": job["queue_position"],
        "user_message": user_message,
        "artifacts": [artifact.model_dump() for artifact in artifacts],
    }
    return JobResponse.model_validate(payload)


async def _store_upload(destination: Path, upload: UploadFile) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as handle:
        shutil.copyfileobj(upload.file, handle)
    await upload.close()
    return destination


def _require_job(job_id: str, user_id: str) -> dict:
    job = get_job_for_user(job_id, user_id)
    if not job:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")
    return job


@router.get("/me", response_model=UserResponse)
def me(identity: dict = Depends(get_request_identity)) -> UserResponse:
    user = ensure_user(identity["email"], role=identity["role"], auth_provider=identity["provider"])
    return UserResponse(
        email=user["email"],
        display_name=user["display_name"],
    )


@router.get("/jobs", response_model=JobListResponse)
def list_jobs(identity: dict = Depends(get_request_identity)) -> JobListResponse:
    user = ensure_user(identity["email"], role=identity["role"], auth_provider=identity["provider"])
    jobs = [_job_to_response(job) for job in list_jobs_for_user(user["id"])]
    return JobListResponse(jobs=jobs)


@router.post("/jobs", response_model=JobResponse, status_code=status.HTTP_201_CREATED)
async def create_job_endpoint(
    public_name: str = Form(...),
    mode: str = Form(...),
    refname: str = Form(...),
    grothmedium: str = Form("min"),
    tmscore: float | None = Form(default=None),
    coverage: float | None = Form(default=None),
    plddt: float | None = Form(default=None),
    fasta: UploadFile = File(...),
    maplist: UploadFile | None = File(default=None),
    structure_archive: UploadFile | None = File(default=None),
    identity: dict = Depends(get_request_identity),
) -> JobResponse:
    if mode not in VALID_MODES:
        raise HTTPException(status_code=400, detail="Invalid mode")
    if refname not in VALID_REFNAMES:
        raise HTTPException(status_code=400, detail="Invalid refname")
    if grothmedium not in VALID_MEDIA:
        raise HTTPException(status_code=400, detail="Invalid grothmedium")
    if mode == "structure_alignment" and (maplist is None or structure_archive is None):
        raise HTTPException(status_code=400, detail="structure_alignment requires maplist and structure_archive")

    user = ensure_user(identity["email"], role=identity["role"], auth_provider=identity["provider"])
    queue_position = len(list_jobs_for_user(user["id"]))
    email = user["email"]
    sanitized_name = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in public_name.strip()) or "job"
    internal_name = f"job_{sanitized_name}_{abs(hash(public_name + email)) % 10_000_000}"
    working_dir = settings.jobs_root / internal_name
    inputs_dir = working_dir / "input"
    logs_dir = working_dir / "logs"
    fasta_path = await _store_upload(inputs_dir / (fasta.filename or "genome.fasta"), fasta)

    maplist_path: Path | None = None
    structure_path: Path | None = None

    if maplist is not None:
        maplist_path = await _store_upload(inputs_dir / (maplist.filename or "maplist.xlsx"), maplist)
    if structure_archive is not None:
        archive_path = await _store_upload(inputs_dir / (structure_archive.filename or "structure.zip"), structure_archive)
        structure_path = inputs_dir / "structure"
        structure_path.mkdir(parents=True, exist_ok=True)
        try:
            shutil.unpack_archive(str(archive_path), str(structure_path))
        except shutil.ReadError as exc:
            raise HTTPException(status_code=400, detail=f"Unsupported structure archive: {exc}") from exc

    log_path = logs_dir / "job.log"
    job = create_job(
        user_id=user["id"],
        public_name=public_name,
        internal_name=internal_name,
        mode=mode,
        refname=refname,
        growth_medium=grothmedium,
        working_dir=working_dir,
        log_path=log_path,
        fasta_path=fasta_path,
        maplist_path=maplist_path,
        structure_path=structure_path,
        tmscore=tmscore,
        coverage=coverage,
        plddt=plddt,
        queue_position=queue_position,
    )
    job_runner.enqueue(job["id"])
    return _job_to_response(job)


@router.get("/jobs/{job_id}", response_model=JobResponse)
def get_job_endpoint(job_id: str, identity: dict = Depends(get_request_identity)) -> JobResponse:
    user = ensure_user(identity["email"], role=identity["role"], auth_provider=identity["provider"])
    job = _require_job(job_id, user["id"])
    return _job_to_response(job)


@router.get("/jobs/{job_id}/logs", response_model=LogResponse)
def get_job_logs(job_id: str, identity: dict = Depends(get_request_identity)) -> LogResponse:
    user = ensure_user(identity["email"], role=identity["role"], auth_provider=identity["provider"])
    job = _require_job(job_id, user["id"])
    log_path = Path(job["log_path"])
    text = log_path.read_text(encoding="utf-8") if log_path.exists() else ""
    return LogResponse(job_id=job_id, text=text[-6000:])


@router.get("/jobs/{job_id}/artifacts/{artifact_id}")
def download_artifact(job_id: str, artifact_id: str, identity: dict = Depends(get_request_identity)) -> FileResponse:
    user = ensure_user(identity["email"], role=identity["role"], auth_provider=identity["provider"])
    _require_job(job_id, user["id"])
    artifact = get_artifact_for_job(artifact_id, job_id)
    if not artifact:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Artifact not found")
    path = Path(artifact["stored_path"])
    if not path.exists():
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Artifact file missing")
    return FileResponse(path, filename=artifact["display_name"])
