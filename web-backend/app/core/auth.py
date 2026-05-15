from __future__ import annotations

import base64
import hashlib
import hmac
import json
import secrets
import time
from typing import Any

from fastapi import Header, HTTPException, status

from app.core.config import settings


def hash_password(password: str, salt: str | None = None) -> str:
    salt_value = salt or secrets.token_hex(16)
    digest = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt_value.encode("utf-8"), 100_000)
    return f"{salt_value}${digest.hex()}"


def verify_password(password: str, stored_hash: str | None) -> bool:
    if not stored_hash or "$" not in stored_hash:
        return False
    salt, _ = stored_hash.split("$", 1)
    expected = hash_password(password, salt)
    return hmac.compare_digest(expected, stored_hash)


def _b64url(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).decode("utf-8").rstrip("=")


def _unb64url(data: str) -> bytes:
    padding = "=" * (-len(data) % 4)
    return base64.urlsafe_b64decode(data + padding)


def create_access_token(email: str, role: str) -> str:
    payload = {
        "email": email.lower(),
        "role": role,
        "exp": int(time.time()) + settings.token_ttl_seconds,
        "iss": "alphagem-local-auth",
    }
    payload_raw = json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    payload_token = _b64url(payload_raw)
    signature = hmac.new(settings.auth_secret.encode("utf-8"), payload_token.encode("utf-8"), hashlib.sha256).hexdigest()
    return f"{payload_token}.{signature}"


def decode_access_token(token: str) -> dict[str, Any]:
    try:
        payload_token, signature = token.split(".", 1)
    except ValueError as exc:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token format") from exc

    expected = hmac.new(settings.auth_secret.encode("utf-8"), payload_token.encode("utf-8"), hashlib.sha256).hexdigest()
    if not hmac.compare_digest(expected, signature):
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid token signature")

    payload = json.loads(_unb64url(payload_token))
    if int(payload.get("exp", 0)) < int(time.time()):
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Token expired")
    return payload


def get_request_identity(
    authorization: str | None = Header(default=None),
    cf_access_email: str | None = Header(default=None, alias="Cf-Access-Authenticated-User-Email"),
    x_dev_user_email: str | None = Header(default=None, alias="X-Dev-User-Email"),
) -> dict[str, str]:
    if authorization and authorization.lower().startswith("bearer "):
        token = authorization.split(" ", 1)[1].strip()
        payload = decode_access_token(token)
        return {"email": payload["email"], "role": payload["role"], "provider": "local"}

    if cf_access_email:
        return {"email": cf_access_email.lower(), "role": "user", "provider": "cloudflare"}

    if x_dev_user_email:
        return {"email": x_dev_user_email.lower(), "role": "user", "provider": "dev-header"}

    if settings.default_user_email:
        return {"email": settings.default_user_email.lower(), "role": "user", "provider": "default"}

    raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Missing authenticated user")
