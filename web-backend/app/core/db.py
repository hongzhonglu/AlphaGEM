from __future__ import annotations

import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Iterator

from app.core.auth import hash_password
from app.core.config import settings


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def ensure_directories() -> None:
    settings.jobs_root.mkdir(parents=True, exist_ok=True)
    settings.uploads_root.mkdir(parents=True, exist_ok=True)
    settings.var_root.mkdir(parents=True, exist_ok=True)


def get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(settings.database_path, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


@contextmanager
def connect() -> Iterator[sqlite3.Connection]:
    conn = get_connection()
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def _table_columns(conn: sqlite3.Connection, table_name: str) -> set[str]:
    rows = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
    return {row["name"] for row in rows}


def _ensure_column(conn: sqlite3.Connection, table_name: str, column_sql: str) -> None:
    column_name = column_sql.split(" ", 1)[0]
    if column_name not in _table_columns(conn, table_name):
        conn.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_sql}")


def init_db() -> None:
    ensure_directories()
    with connect() as conn:
        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS users (
                id TEXT PRIMARY KEY,
                email TEXT NOT NULL UNIQUE,
                display_name TEXT,
                role TEXT NOT NULL DEFAULT 'user',
                auth_provider TEXT NOT NULL DEFAULT 'local',
                password_hash TEXT,
                created_at TEXT NOT NULL,
                last_seen_at TEXT NOT NULL
            );

            CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY,
                user_id TEXT NOT NULL,
                public_name TEXT NOT NULL,
                internal_name TEXT NOT NULL UNIQUE,
                mode TEXT NOT NULL,
                refname TEXT NOT NULL,
                growth_medium TEXT NOT NULL,
                status TEXT NOT NULL,
                created_at TEXT NOT NULL,
                started_at TEXT,
                finished_at TEXT,
                working_dir TEXT NOT NULL,
                output_dir TEXT,
                log_path TEXT NOT NULL,
                fasta_path TEXT NOT NULL,
                maplist_path TEXT,
                structure_path TEXT,
                tmscore REAL,
                coverage REAL,
                plddt REAL,
                exit_code INTEGER,
                error_message TEXT,
                queue_position INTEGER
            );

            CREATE TABLE IF NOT EXISTS job_artifacts (
                id TEXT PRIMARY KEY,
                job_id TEXT NOT NULL,
                kind TEXT NOT NULL,
                stored_path TEXT NOT NULL,
                display_name TEXT NOT NULL,
                size_bytes INTEGER NOT NULL,
                created_at TEXT NOT NULL,
                FOREIGN KEY(job_id) REFERENCES jobs(id)
            );
            """
        )

        _ensure_column(conn, "users", "display_name TEXT")
        _ensure_column(conn, "users", "role TEXT NOT NULL DEFAULT 'user'")
        _ensure_column(conn, "users", "auth_provider TEXT NOT NULL DEFAULT 'local'")
        _ensure_column(conn, "users", "password_hash TEXT")
        _ensure_column(conn, "users", "created_at TEXT")
        _ensure_column(conn, "users", "last_seen_at TEXT")

        seed_admin_user(conn)


def seed_admin_user(conn: sqlite3.Connection) -> None:
    email = settings.seed_admin_email.lower()
    row = conn.execute("SELECT * FROM users WHERE email = ?", (email,)).fetchone()
    now = utc_now()
    if row:
        updates: list[tuple[str, object]] = []
        if row["role"] != "admin":
            updates.append(("role", "admin"))
        if not row["password_hash"]:
            updates.append(("password_hash", hash_password(settings.seed_admin_password)))
        if row["auth_provider"] != "local":
            updates.append(("auth_provider", "local"))
        if updates:
            assignments = ", ".join(f"{name} = ?" for name, _ in updates)
            values = [value for _, value in updates] + [email]
            conn.execute(f"UPDATE users SET {assignments} WHERE email = ?", values)
        return

    conn.execute(
        """
        INSERT INTO users (id, email, display_name, role, auth_provider, password_hash, created_at, last_seen_at)
        VALUES (?, ?, ?, 'admin', 'local', ?, ?, ?)
        """,
        (
            f"user_{uuid.uuid4().hex}",
            email,
            "AlphaGEM Admin",
            hash_password(settings.seed_admin_password),
            now,
            now,
        ),
    )


def ensure_user(email: str, *, role: str = "user", auth_provider: str = "local") -> sqlite3.Row:
    with connect() as conn:
        row = conn.execute("SELECT * FROM users WHERE email = ?", (email,)).fetchone()
        now = utc_now()
        if row:
            conn.execute("UPDATE users SET last_seen_at = ? WHERE email = ?", (now, email))
            row = conn.execute("SELECT * FROM users WHERE email = ?", (email,)).fetchone()
            return row

        user_id = f"user_{uuid.uuid4().hex}"
        conn.execute(
            """
            INSERT INTO users (id, email, display_name, role, auth_provider, password_hash, created_at, last_seen_at)
            VALUES (?, ?, ?, ?, ?, NULL, ?, ?)
            """,
            (user_id, email, email.split("@")[0], role, auth_provider, now, now),
        )
        row = conn.execute("SELECT * FROM users WHERE email = ?", (email,)).fetchone()
        return row


def create_local_user(email: str, password: str, display_name: str | None = None, role: str = "user") -> sqlite3.Row:
    with connect() as conn:
        existing = conn.execute("SELECT * FROM users WHERE email = ?", (email,)).fetchone()
        if existing:
            raise ValueError("User already exists")

        now = utc_now()
        user_id = f"user_{uuid.uuid4().hex}"
        conn.execute(
            """
            INSERT INTO users (id, email, display_name, role, auth_provider, password_hash, created_at, last_seen_at)
            VALUES (?, ?, ?, ?, 'local', ?, ?, ?)
            """,
            (
                user_id,
                email,
                display_name or email.split("@")[0],
                role,
                hash_password(password),
                now,
                now,
            ),
        )
        row = conn.execute("SELECT * FROM users WHERE email = ?", (email,)).fetchone()
        return row


def get_user_by_email(email: str) -> sqlite3.Row | None:
    with connect() as conn:
        return conn.execute("SELECT * FROM users WHERE email = ?", (email.lower(),)).fetchone()


def touch_user(email: str) -> sqlite3.Row | None:
    with connect() as conn:
        now = utc_now()
        conn.execute("UPDATE users SET last_seen_at = ? WHERE email = ?", (now, email.lower()))
        return conn.execute("SELECT * FROM users WHERE email = ?", (email.lower(),)).fetchone()
