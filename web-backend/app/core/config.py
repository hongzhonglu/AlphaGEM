from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    repo_root: Path
    storage_root: Path
    jobs_root: Path
    uploads_root: Path
    var_root: Path
    database_path: Path
    max_upload_size_mb: int
    default_user_email: str
    auth_secret: str
    token_ttl_seconds: int
    seed_admin_email: str
    seed_admin_password: str


def _path_from_env(name: str, default: Path) -> Path:
    value = os.getenv(name)
    return Path(value).expanduser().resolve() if value else default


REPO_ROOT = Path(__file__).resolve().parents[3]
STORAGE_ROOT = _path_from_env("ALPHAGEM_STORAGE_ROOT", REPO_ROOT / "storage")
VAR_ROOT = _path_from_env("ALPHAGEM_VAR_ROOT", REPO_ROOT / "var")

settings = Settings(
    repo_root=REPO_ROOT,
    storage_root=STORAGE_ROOT,
    jobs_root=STORAGE_ROOT / "jobs",
    uploads_root=STORAGE_ROOT / "uploads",
    var_root=VAR_ROOT,
    database_path=_path_from_env("ALPHAGEM_DB_PATH", VAR_ROOT / "alphagem.db"),
    max_upload_size_mb=int(os.getenv("ALPHAGEM_MAX_UPLOAD_MB", "512")),
    default_user_email=os.getenv("ALPHAGEM_DEFAULT_USER_EMAIL", "local-dev@alphagem"),
    auth_secret=os.getenv("ALPHAGEM_AUTH_SECRET", "change-me-before-public-deploy"),
    token_ttl_seconds=int(os.getenv("ALPHAGEM_TOKEN_TTL_SECONDS", "86400")),
    seed_admin_email=os.getenv("ALPHAGEM_ADMIN_EMAIL", "admin@alphagem.local"),
    seed_admin_password=os.getenv("ALPHAGEM_ADMIN_PASSWORD", "AlphaGEM_admin_123"),
)
