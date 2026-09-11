"""Settings API: read/update persistent app settings (LLM config + theme).

GET  /api/v1/settings -> redacted view (api_key as configured + last-4 hint)
POST /api/v1/settings -> persist updates; api_key empty means keep current.
"""
from __future__ import annotations

from typing import Any

from fastapi import APIRouter
from pydantic import BaseModel, Field

from server.backend.settings_store import get_store

router = APIRouter(prefix="/api/v1", tags=["settings"])


class SettingsBody(BaseModel):
    model: str | None = Field(default=None, min_length=1, max_length=200)
    base_url: str | None = Field(default=None, min_length=1, max_length=500)
    api_key: str | None = Field(default=None, max_length=500)  # 空/省略 = 保留原值
    concurrency: int | None = Field(default=None, ge=1, le=64)
    theme: str | None = Field(default=None, pattern="^(dark|light)$")


@router.get("/settings")
def get_settings() -> dict[str, Any]:
    return {"ok": True, "settings": get_store().public()}


@router.post("/settings")
def save_settings(body: SettingsBody) -> dict[str, Any]:
    store = get_store()
    updates: dict[str, Any] = {}
    if body.model is not None:
        updates["model"] = body.model.strip()
    if body.base_url is not None:
        updates["base_url"] = body.base_url.strip()
    if body.api_key is not None and body.api_key.strip():
        updates["api_key"] = body.api_key.strip()  # 空字符串 → 不更新
    if body.concurrency is not None:
        updates["concurrency"] = body.concurrency
    if body.theme is not None:
        updates["theme"] = body.theme
    saved = store.save(updates)
    return {"ok": True, "settings": store.public(saved)}
