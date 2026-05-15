from __future__ import annotations

from fastapi import APIRouter, HTTPException, status

from app.core.auth import create_access_token, verify_password
from app.core.db import create_local_user, get_user_by_email, touch_user
from app.schemas import AuthRequest, AuthResponse, UserResponse

router = APIRouter()


def _user_response(row) -> UserResponse:
    return UserResponse(
        email=row["email"],
        display_name=row["display_name"],
    )


@router.post("/register", response_model=AuthResponse, status_code=status.HTTP_201_CREATED)
def register(payload: AuthRequest) -> AuthResponse:
    email = payload.email.lower().strip()
    if "@" not in email:
        raise HTTPException(status_code=400, detail="Please enter a valid email address")
    if len(payload.password) < 8:
        raise HTTPException(status_code=400, detail="Password must be at least 8 characters long")
    try:
        user = create_local_user(email=email, password=payload.password, display_name=payload.display_name)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    token = create_access_token(user["email"], user["role"])
    return AuthResponse(access_token=token, user=_user_response(user))


@router.post("/login", response_model=AuthResponse)
def login(payload: AuthRequest) -> AuthResponse:
    email = payload.email.lower().strip()
    user = get_user_by_email(email)
    if not user or not verify_password(payload.password, user["password_hash"]):
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Incorrect email or password")
    user = touch_user(email)
    token = create_access_token(user["email"], user["role"])
    return AuthResponse(access_token=token, user=_user_response(user))
