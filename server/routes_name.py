"""Namer API: SMILES → bilingual NameResult JSON."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter
from pydantic import BaseModel, Field

from server.atom_ids_svg import build_atom_ids_svg
from server.deps import get_namer, name_result_dict
from server.locants_svg import build_locants_svg

router = APIRouter(prefix="/api/v1", tags=["name"])


class NameBody(BaseModel):
    """Request body for POST /name."""

    smiles: str = Field(..., min_length=1)


class LocantsBody(BaseModel):
    """Request body for POST /name/locants-svg."""

    smiles: str = Field(..., min_length=1)


@router.post("/name")
def name_smiles(body: NameBody) -> dict[str, Any]:
    """Run SMILESNNamer.name and return serialized NameResult."""
    result = get_namer().name(body.smiles)
    return name_result_dict(result)


@router.post("/name/locants-svg")
def name_locants_svg(body: LocantsBody) -> dict[str, Any]:
    """SMILES → 带 L4 编号标注的结构图 SVG;编号不可用返回 ok:false。"""
    res = build_locants_svg(body.smiles)
    return {"ok": True, **res} if res else {"ok": False}


@router.post("/name/atom-ids-svg")
def name_atom_ids_svg(body: LocantsBody) -> dict[str, Any]:
    """SMILES → 标注原子/SSSR 环索引的结构图 SVG;不可用返回 ok:false。"""
    res = build_atom_ids_svg(body.smiles)
    return {"ok": True, **res} if res else {"ok": False}
