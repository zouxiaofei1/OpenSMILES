"""History-version support: git commit listing + worktree/cache housekeeping.

GET  /api/v1/git/commits — HEAD-lineage commit list (used by the topbar picker)
POST /api/v1/history/gc  — sweep leftover history worktrees (optional)
"""
from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Query

from server import history_store

router = APIRouter(prefix="/api/v1", tags=["history"])


@router.get("/git/commits")
def git_commits(refresh: bool = Query(False)) -> dict[str, Any]:
    """Return all HEAD-lineage commits: {ok, head, commits:[{hash,short,date,subject}]}."""
    return history_store.list_commits(force=refresh)


@router.post("/history/gc")
def history_gc() -> dict[str, Any]:
    """Remove leftover history worktrees (startup sweep)."""
    history_store.cleanup_orphans()
    return {"ok": True}
