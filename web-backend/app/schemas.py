from __future__ import annotations

from datetime import datetime

from pydantic import BaseModel, Field


class UserResponse(BaseModel):
    email: str
    display_name: str | None = None


class AuthRequest(BaseModel):
    email: str
    password: str
    display_name: str | None = None


class AuthResponse(BaseModel):
    access_token: str
    token_type: str = "bearer"
    user: UserResponse


class ArtifactResponse(BaseModel):
    id: str
    kind: str
    display_name: str
    size_bytes: int
    created_at: datetime


class JobResponse(BaseModel):
    id: str
    public_name: str
    mode: str
    refname: str
    growth_medium: str = Field(alias="growthMedium")
    status: str
    created_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    tmscore: float | None = None
    coverage: float | None = None
    plddt: float | None = None
    queue_position: int | None = None
    user_message: str | None = None
    artifacts: list[ArtifactResponse] = []

    model_config = {"populate_by_name": True}


class JobListResponse(BaseModel):
    jobs: list[JobResponse]


class LogResponse(BaseModel):
    job_id: str
    text: str
