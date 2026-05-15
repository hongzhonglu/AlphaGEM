from __future__ import annotations

import sqlite3
import uuid
from pathlib import Path

from app.core.db import connect, utc_now


def row_to_dict(row: sqlite3.Row | None) -> dict | None:
    return dict(row) if row is not None else None


def create_job(
    *,
    user_id: str,
    public_name: str,
    internal_name: str,
    mode: str,
    refname: str,
    growth_medium: str,
    working_dir: Path,
    log_path: Path,
    fasta_path: Path,
    maplist_path: Path | None,
    structure_path: Path | None,
    tmscore: float | None,
    coverage: float | None,
    plddt: float | None,
    queue_position: int,
) -> dict:
    job_id = str(uuid.uuid4())
    now = utc_now()
    with connect() as conn:
        conn.execute(
            """
            INSERT INTO jobs (
                id, user_id, public_name, internal_name, mode, refname, growth_medium,
                status, created_at, working_dir, log_path, fasta_path, maplist_path,
                structure_path, tmscore, coverage, plddt, queue_position
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, 'queued', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                job_id,
                user_id,
                public_name,
                internal_name,
                mode,
                refname,
                growth_medium,
                now,
                str(working_dir),
                str(log_path),
                str(fasta_path),
                str(maplist_path) if maplist_path else None,
                str(structure_path) if structure_path else None,
                tmscore,
                coverage,
                plddt,
                queue_position,
            ),
        )
        return get_job(job_id)


def list_jobs_for_user(user_id: str) -> list[dict]:
    with connect() as conn:
        rows = conn.execute(
            "SELECT * FROM jobs WHERE user_id = ? ORDER BY created_at DESC",
            (user_id,),
        ).fetchall()
        return [dict(row) for row in rows]


def get_job(job_id: str) -> dict | None:
    with connect() as conn:
        row = conn.execute("SELECT * FROM jobs WHERE id = ?", (job_id,)).fetchone()
        return row_to_dict(row)


def get_job_for_user(job_id: str, user_id: str) -> dict | None:
    with connect() as conn:
        row = conn.execute(
            "SELECT * FROM jobs WHERE id = ? AND user_id = ?",
            (job_id, user_id),
        ).fetchone()
        return row_to_dict(row)


def update_job(job_id: str, **fields: object) -> dict | None:
    if not fields:
        return get_job(job_id)

    assignments = ", ".join(f"{key} = ?" for key in fields)
    values = list(fields.values()) + [job_id]
    with connect() as conn:
        conn.execute(f"UPDATE jobs SET {assignments} WHERE id = ?", values)
        return get_job(job_id)


def add_artifact(job_id: str, kind: str, stored_path: Path, display_name: str, size_bytes: int) -> dict:
    artifact_id = str(uuid.uuid4())
    now = utc_now()
    with connect() as conn:
        conn.execute(
            """
            INSERT INTO job_artifacts (id, job_id, kind, stored_path, display_name, size_bytes, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (artifact_id, job_id, kind, str(stored_path), display_name, size_bytes, now),
        )
        row = conn.execute("SELECT * FROM job_artifacts WHERE id = ?", (artifact_id,)).fetchone()
        return dict(row)


def list_artifacts(job_id: str) -> list[dict]:
    with connect() as conn:
        rows = conn.execute(
            "SELECT * FROM job_artifacts WHERE job_id = ? ORDER BY created_at ASC",
            (job_id,),
        ).fetchall()
        return [dict(row) for row in rows]


def get_artifact_for_job(artifact_id: str, job_id: str) -> dict | None:
    with connect() as conn:
        row = conn.execute(
            "SELECT * FROM job_artifacts WHERE id = ? AND job_id = ?",
            (artifact_id, job_id),
        ).fetchone()
        return row_to_dict(row)
