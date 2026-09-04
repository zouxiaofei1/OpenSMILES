"""Namer API: SMILES → bilingual NameResult JSON."""

from __future__ import annotations

import json
import threading
from collections import OrderedDict
from pathlib import Path
from typing import Any

import httpx
from fastapi import APIRouter
from pydantic import BaseModel, Field

from server.atom_ids_svg import build_atom_ids_svg
from server.deps import get_namer, name_result_dict
from server.locants_svg import build_locants_svg

router = APIRouter(prefix="/api/v1", tags=["name"])

ROOT = Path(__file__).resolve().parents[1]
GOLD_SOURCE = ROOT / "data" / "merged_benchmark.json"

# 懒加载的基准索引: canonical SMILES → gold 记录。合并了合并基准中 4010 行的
# {english_name, chinese_name, tier, source}, 供 Namer 页命名时附带标准答案。
_gold_index: dict[str, dict[str, Any]] | None = None
_gold_lock = threading.Lock()


def _build_gold_index() -> dict[str, dict[str, Any]]:
    """Read merged_benchmark.json → {canonical SMILES: gold record}. First wins on dup."""
    from rdkit import Chem

    idx: dict[str, dict[str, Any]] = {}
    try:
        rows = json.loads(GOLD_SOURCE.read_text(encoding="utf-8"))
    except Exception:
        return idx
    if not isinstance(rows, list):
        return idx
    for r in rows:
        smi = (r.get("smiles") or "").strip()
        try:
            mol = Chem.MolFromSmiles(smi)
            cs = Chem.MolToSmiles(mol) if mol is not None else ""
        except Exception:
            continue
        if not cs or cs in idx:
            continue
        idx[cs] = {
            "id": r.get("id"),
            "en": r.get("english_name") or "",
            "zh": r.get("chinese_name") or "",
            "tier": r.get("tier"),
            "source": r.get("source"),
        }
    return idx


def lookup_gold(smiles: str) -> dict[str, Any] | None:
    """Find the merged_benchmark record for a SMILES (canonical match), else None.

    Index is built once on first call. Never raises: any failure → None so a
    gold lookup never disturbs the naming main flow.
    """
    global _gold_index
    with _gold_lock:
        if _gold_index is None:
            _gold_index = _build_gold_index()
        idx = _gold_index
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles((smiles or "").strip())
        if mol is None:
            return None
        rec = idx.get(Chem.MolToSmiles(mol))
    except Exception:
        return None
    return rec if rec is not None else None


# ── PubChem 查询: PUG-REST property/IUPACName（PubChem 页 2.1.1 IUPAC Name）──────
PUG_IUPAC_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{ik}"
    "/property/IUPACName/JSON"
)
PUG_USER_AGENT = "ChemAgentNamer/0.1 (local IUPAC-name lookup)"

# 以 InChIKey 缓存稳定结果(命中/未收录); 网络异常不缓存。上限防无限增长。
_pc_cache: OrderedDict[str, dict[str, Any]] = OrderedDict()
_pc_cache_lock = threading.Lock()
_PC_CACHE_MAX = 512


def _pc_cache_get(key: str) -> dict[str, Any] | None:
    with _pc_cache_lock:
        return _pc_cache.get(key)


def _pc_cache_put(key: str, val: dict[str, Any]) -> None:
    with _pc_cache_lock:
        _pc_cache[key] = val
        _pc_cache.move_to_end(key)
        while len(_pc_cache) > _PC_CACHE_MAX:
            _pc_cache.popitem(last=False)


def _pc_query(ik: str) -> tuple[bool, dict[str, Any]]:
    """GET PUG-REST IUPACName for an InChIKey. Returns (ok, normalized dict).

    ok=False means a non-recoverable failure (network/timeout/bad response);
    the caller must NOT cache that. Found/not-found both come back ok=True.
    """
    try:
        r = httpx.get(
            PUG_IUPAC_URL.format(ik=ik),
            timeout=20.0,
            headers={"User-Agent": PUG_USER_AGENT},
        )
    except httpx.HTTPError as exc:
        return False, {"ok": False, "error": f"PubChem 请求失败: {exc.__class__.__name__}"}
    if r.status_code == 404:
        return True, {"ok": True, "found": False, "cid": None, "iupac": ""}
    if r.status_code != 200:
        return False, {"ok": False, "error": f"PubChem HTTP {r.status_code}"}
    try:
        data = r.json()
        if "Fault" in data:  # PUG 偶发以 200 返回 Fault
            return True, {"ok": True, "found": False, "cid": None, "iupac": ""}
        props = data.get("PropertyTable", {}).get("Properties") or []
    except Exception:
        return False, {"ok": False, "error": "PubChem 响应解析失败"}
    if not props:
        return True, {"ok": True, "found": False, "cid": None, "iupac": ""}
    rec = props[0]
    cid = rec.get("CID")
    iupac = (rec.get("IUPACName") or "").strip()
    return True, {
        "ok": True,
        "found": bool(iupac),
        "cid": cid,
        "iupac": iupac,
        "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else None,
    }


def _pubchem_ik_lookup(ik: str) -> dict[str, Any]:
    """Cached PubChem lookup by InChIKey; network failures bypass the cache."""
    cached = _pc_cache_get(ik)
    if cached is not None:
        return dict(cached)
    ok, result = _pc_query(ik)
    if ok:
        _pc_cache_put(ik, result)
    return result


def pubchem_iupac(smiles: str) -> dict[str, Any]:
    """SMILES → PubChem 2.1.1 computed IUPACName via PUG-REST.

    InChIKey (not the raw SMILES) is the query key to dodge slash/ring-closure
    URL-encoding issues. When the exact stereoisomer has no PubChem record, retry
    once without stereochemistry (PubChem usually indexes the racemic/parent
    form), flagging the fallback via `racemic`. Never raises.
    """
    from rdkit import Chem
    from rdkit.Chem import inchi

    try:
        mol = Chem.MolFromSmiles((smiles or "").strip())
    except Exception:
        mol = None
    if mol is None:
        return {"ok": False, "error": "无法解析 SMILES"}
    try:
        ik = inchi.MolToInchiKey(mol)
    except Exception:
        return {"ok": False, "error": "无法生成 InChIKey"}

    res = _pubchem_ik_lookup(ik)
    if res["ok"] and not res["found"] and res.get("cid") is None and ik:
        # 精确立体未收录 → 去立体重查父/外消旋记录
        try:
            flat = Chem.RemoveStereochemistry(Chem.Mol(mol))
            ik2 = inchi.MolToInchiKey(flat)
        except Exception:
            ik2 = ""
        if ik2 and ik2 != ik:
            res2 = _pubchem_ik_lookup(ik2)
            if res2["ok"] and res2.get("cid") is not None:
                res2 = dict(res2, racemic=True)
                return res2
    return res


class NameBody(BaseModel):
    """Request body for POST /name."""

    smiles: str = Field(..., min_length=1)


class LocantsBody(BaseModel):
    """Request body for POST /name/locants-svg."""

    smiles: str = Field(..., min_length=1)
    orient: bool = True  # 开关: 稠环按 preferred_orientation 水平摆放


@router.post("/name")
def name_smiles(body: NameBody) -> dict[str, Any]:
    """Run SMILESNNamer.name and return serialized NameResult + merged_benchmark gold."""
    result = get_namer().name(body.smiles)
    payload = name_result_dict(result)
    payload["gold"] = lookup_gold(body.smiles)
    return payload


class PubchemBody(BaseModel):
    """Request body for POST /name/pubchem-iupac."""

    smiles: str = Field(..., min_length=1)


@router.post("/name/pubchem-iupac")
def name_pubchem_iupac(body: PubchemBody) -> dict[str, Any]:
    """SMILES → PubChem 2.1.1 (computed) IUPAC name via PUG-REST proxy."""
    return pubchem_iupac(body.smiles)


@router.post("/name/locants-svg")
def name_locants_svg(body: LocantsBody) -> dict[str, Any]:
    """SMILES → 带 L4 编号标注的结构图 SVG;编号不可用返回 ok:false。"""
    res = build_locants_svg(body.smiles, orient=body.orient)
    return {"ok": True, **res} if res else {"ok": False}


@router.post("/name/atom-ids-svg")
def name_atom_ids_svg(body: LocantsBody) -> dict[str, Any]:
    """SMILES → 标注原子/SSSR 环索引的结构图 SVG;不可用返回 ok:false。"""
    res = build_atom_ids_svg(body.smiles)
    return {"ok": True, **res} if res else {"ok": False}
