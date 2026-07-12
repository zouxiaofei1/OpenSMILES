"""Settings routes: API key configuration."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter
from pydantic import BaseModel, Field

from server import deps

router = APIRouter(prefix="/api/v1/settings", tags=["settings"])


class ApiKeyBody(BaseModel):
    """Body for setting Anthropic API key."""

    api_key: str = Field(default="")


@router.get("/api-key")
def get_api_key() -> dict[str, Any]:
    """Return public status of configured API key (no full key)."""
    return deps.get_secrets().public_status()


@router.post("/api-key")
def post_api_key(body: ApiKeyBody) -> dict[str, Any]:
    """Store API key and return public status."""
    deps.get_secrets().set_key(body.api_key)
    return deps.get_secrets().public_status()