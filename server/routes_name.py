"""Namer API: SMILES → bilingual NameResult JSON."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter
from pydantic import BaseModel, Field

from server.deps import get_namer, name_result_dict

router = APIRouter(prefix="/api/v1", tags=["name"])


class NameBody(BaseModel):
    """Request body for POST /name."""

    smiles: str = Field(..., min_length=1)


@router.post("/name")
def name_smiles(body: NameBody) -> dict[str, Any]:
    """Run SMILESNNamer.name and return serialized NameResult."""
    result = get_namer().name(body.smiles)
    return name_result_dict(result)
