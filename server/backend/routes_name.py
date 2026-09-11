"""Namer API: SMILES → bilingual NameResult JSON."""

from __future__ import annotations

import difflib
import json
import threading
from collections import OrderedDict
from pathlib import Path
from typing import Any, Literal

from server.backend.legacy_engines import name_result as legacy_name_result
from server.backend.ml_engine import name_result as ml_name_result

import httpx
from fastapi import APIRouter
from pydantic import BaseModel, Field

from benchmarks.preview_metrics import similarity
from server.backend.atom_ids_svg import build_atom_ids_svg
from server.backend.deps import get_namer, name_result_dict
from server.backend.locants_svg import build_locants_svg

router = APIRouter(prefix="/api/v1", tags=["name"])

ROOT = Path(__file__).resolve().parents[2]
GOLD_SOURCE = ROOT / "benchmarks" / "merged_benchmark.json"

# 懒加载的基准索引(同一份 merged_benchmark.json 一次解析产出两张表):
#   canonical SMILES → gold 记录 — Namer 页命名时附带标准答案
#   {id 字符串 → 行记录}       — 输入框可输 chebi-1 / tiers-1 这类 id 直接定位到某行
_gold_index: dict[str, dict[str, Any]] | None = None
_id_index: dict[str, dict[str, Any]] | None = None
_index_lock = threading.Lock()


def _record(r: dict[str, Any]) -> dict[str, Any]:
    """merged_benchmark 行 → 对外记录(含 smiles, 供 id 解析出真实 SMILES)。"""
    return {
        "id": r.get("id"),
        "smiles": (r.get("smiles") or "").strip(),
        "en": r.get("english_name") or "",
        "zh": r.get("chinese_name") or "",
        "tier": r.get("tier"),
        "source": r.get("source"),
    }


def _load_indexes() -> None:
    """Read merged_benchmark.json once → fill canonical-SMILES & exact-id maps.

    两表均 first-wins 去重; 文件读/解析失败时置空表, 不抛异常, 保证查询路径安全。
    """
    global _gold_index, _id_index
    try:
        rows = json.loads(GOLD_SOURCE.read_text(encoding="utf-8"))
    except Exception:
        rows = []
    if not isinstance(rows, list):
        rows = []
    gold: dict[str, dict[str, Any]] = {}
    ids: dict[str, dict[str, Any]] = {}
    from rdkit import Chem

    for r in rows:
        rec = _record(r)
        rid = rec["id"]
        if rid and rid not in ids:
            ids[rid] = rec
        smi = rec["smiles"]
        try:
            mol = Chem.MolFromSmiles(smi)
            cs = Chem.MolToSmiles(mol) if mol is not None else ""
        except Exception:
            cs = ""
        if cs and cs not in gold:
            gold[cs] = rec
    _gold_index = gold
    _id_index = ids


def _ensure_indexes() -> None:
    """Build both maps once (double-checked under lock); subsequent calls are cheap."""
    global _gold_index, _id_index
    if _gold_index is not None and _id_index is not None:
        return
    with _index_lock:
        if _gold_index is None or _id_index is None:
            _load_indexes()


def lookup_gold(smiles: str) -> dict[str, Any] | None:
    """Find the merged_benchmark record for a SMILES (canonical match), else None.

    Index is built once on first call. Never raises: any failure → None so a
    gold lookup never disturbs the naming main flow.
    """
    _ensure_indexes()
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles((smiles or "").strip())
        if mol is None:
            return None
        rec = _gold_index.get(Chem.MolToSmiles(mol))
    except Exception:
        return None
    return rec if rec is not None else None


def lookup_id_smiles(text: str | None) -> dict[str, Any] | None:
    """Exact merged_benchmark id match (e.g. "chebi-1" → record), else None."""
    _ensure_indexes()
    return _id_index.get((text or "").strip()) if text is not None else None


def _id_or_text(text: str) -> str:
    """Input that is a benchmark id → that row's stored SMILES, else text unchanged.

    Lets the whole name-family of endpoints accept either a SMILES or an id like
    chebi-1 / tiers-1. Non-matching input passes through so the normal SMILES
    parse-error path still handles garbage.
    """
    rec = lookup_id_smiles(text)
    if rec and rec["smiles"]:
        return rec["smiles"]
    return text


# ── 引擎名 ↔ 基准名 差异(Namer 页 gold 卡的并排高亮)──────────────────────────
# 相似度阈值: 归一化后引擎名与基准名达到该相似度才高亮差异。命得基本对时差异
# 才是可看的细节(如 (furan-3-yl) 少一层括号); 错得离谱时整名标红没有信息量。
GOLD_DIFF_MIN_SIM = 0.70


def _name_similarity(pred: str, gold: str, lang: str) -> float:
    """归一化后的字符级相似度 0..1, 用于判断差异是否值得高亮。

    与 benchmark 预览页共用同一实现(preview_metrics), 两边百分比一致;
    gold 非空(调用处已保证)故不会取到 None。"""
    return similarity(pred, gold, lang) or 0.0


def _diff_sides(pred: str, gold: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """逐字符 diff → (引擎侧片段, gold 侧片段), 片段形如 {"t": 文本, "same": bool}。

    在原始串(未归一)上比对, 前端高亮的就是名称原貌; 相邻同类片段合并以减少 DOM 节点。
    """
    left: list[dict[str, Any]] = []
    right: list[dict[str, Any]] = []

    def push(acc: list[dict[str, Any]], text: str, same: bool) -> None:
        if not text:
            return
        if acc and acc[-1]["same"] == same:
            acc[-1]["t"] += text
        else:
            acc.append({"t": text, "same": same})

    matcher = difflib.SequenceMatcher(None, pred, gold, autojunk=False)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        same = tag == "equal"
        push(left, pred[i1:i2], same)
        push(right, gold[j1:j2], same)
    return left, right


def build_gold_diff(
    pred_en: str, pred_zh: str, gold: dict[str, Any] | None
) -> dict[str, Any] | None:
    """引擎名 ↔ 基准名 的相似度与逐字符差异; 无 gold 返回 None。

    每种语言各给 {similarity, show, pred, gold}: show 由相似度是否 ≥ 阈值决定,
    前端据此只在"命得基本对"时高亮细节。gold 缺该语言名称则跳过该语言。
    """
    if not gold:
        return None
    out: dict[str, Any] = {}
    for lang, pred in (("en", pred_en), ("zh", pred_zh)):
        ref = (gold.get(lang) or "").strip()
        if not ref:
            continue
        pred = (pred or "").strip()
        sim = _name_similarity(pred, ref, lang)
        left, right = _diff_sides(pred, ref)
        out[lang] = {
            "similarity": round(sim, 4),
            "show": sim >= GOLD_DIFF_MIN_SIM,
            "pred": left,
            "gold": right,
        }
    return out or None


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
    """Request body for POST /name. `smiles` 可以是 SMILES 或 merged_benchmark id。

    `engine` 选命名引擎: "src"=现役规则引擎; "v2"/"v3"=tools 历史规则引擎;
    "ml"=本地神经网络 SMILES2IUPAC(仅英文)。
    """

    smiles: str = Field(..., min_length=1)
    engine: Literal["src", "v2", "v3", "ml"] = "src"


class ResolveBody(BaseModel):
    """Request body for POST /name/resolve."""

    text: str = Field(..., min_length=1)


class LocantsBody(BaseModel):
    """Request body for POST /name/locants-svg. `smiles` 亦容忍基准 id。"""

    smiles: str = Field(..., min_length=1)
    orient: bool = True  # 开关: 稠环按 preferred_orientation 水平摆放


@router.post("/name")
def name_smiles(body: NameBody) -> dict[str, Any]:
    """Run src/v3 engine .name and return serialized result + merged_benchmark gold.

    入参若是 merged_benchmark id(chebi-1 / tiers-1)先解析为该行 SMILES 再命名。
    gold 与引擎无关(基准答案), 两个引擎都附上便于对照。
    """
    smi = _id_or_text(body.smiles)
    if body.engine == "ml":
        payload = ml_name_result(smi)
    elif body.engine == "src":
        result = get_namer().name(smi)
        payload = name_result_dict(result)
    else:  # v2 / v3 历史规则引擎
        payload = legacy_name_result(smi, body.engine)
    payload["engine"] = body.engine
    gold = lookup_gold(smi)
    payload["gold"] = gold
    # 命名成功才谈得上差异; 失败时前端本就隐藏 gold 卡, 无需白算一遍 diff。
    payload["gold_diff"] = (
        build_gold_diff(payload.get("en") or "", payload.get("zh") or "", gold)
        if gold and payload.get("success")
        else None
    )
    return payload


@router.post("/name/resolve")
def resolve_name(body: ResolveBody) -> dict[str, Any]:
    """Resolve an input that may be a merged_benchmark id into the SMILES to use.

    命中 → {"ok": true, "kind": "id", "id", "smiles"};否则原样放行
    (kind "smiles"), 由下游照常按 SMILES 解析(坏输入报错路径不变)。
    """
    text = (body.text or "").strip()
    rec = lookup_id_smiles(text)
    if rec and rec["smiles"]:
        return {"ok": True, "kind": "id", "id": rec["id"], "smiles": rec["smiles"]}
    return {"ok": True, "kind": "smiles", "smiles": text}


class PubchemBody(BaseModel):
    """Request body for POST /name/pubchem-iupac."""

    smiles: str = Field(..., min_length=1)


@router.post("/name/pubchem-iupac")
def name_pubchem_iupac(body: PubchemBody) -> dict[str, Any]:
    """SMILES → PubChem 2.1.1 (computed) IUPAC name via PUG-REST proxy."""
    return pubchem_iupac(_id_or_text(body.smiles))


@router.post("/name/locants-svg")
def name_locants_svg(body: LocantsBody) -> dict[str, Any]:
    """SMILES → 带 L4 编号标注的结构图 SVG;编号不可用返回 ok:false。"""
    res = build_locants_svg(_id_or_text(body.smiles), orient=body.orient)
    return {"ok": True, **res} if res else {"ok": False}


@router.post("/name/atom-ids-svg")
def name_atom_ids_svg(body: LocantsBody) -> dict[str, Any]:
    """SMILES → 标注原子/SSSR 环索引的结构图 SVG;不可用返回 ok:false。"""
    res = build_atom_ids_svg(_id_or_text(body.smiles))
    return {"ok": True, **res} if res else {"ok": False}
