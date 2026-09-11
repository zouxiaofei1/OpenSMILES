"""预览行的派生指标：稠环判定、RDKit 复杂度、重原子数、预测↔金标相似度。

由预览生成器 (benchmark_preview_parallel) 与 server 的旧缓存补算路径
(routes_benchmark) 共用，保证"新生成的行"与"补算出来的行"字段完全一致。
相似度取归一化后的字符级 ratio，与 Namer 页 gold 卡 (routes_name) 同口径，
两个页面对同一对名字给出的百分比不会打架。

字段: fused(稠环) / cx(BertzCT 复杂度) / hac(重原子数) / sim_en / sim_zh(相似度, 无金标为 None)。
"""

from __future__ import annotations

import difflib
from typing import Any

# 派生字段全集; 旧缓存缺任一键即触发一次性补算。
DERIVED_KEYS = ("fused", "cx", "hac", "sim_en", "sim_zh")


def is_fused(smiles: str) -> bool:
    """是否含稠合环系: 至少两个环共享 ≥2 原子(P-25.3.3 稠环范畴)。"""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    rings = mol.GetRingInfo().AtomRings()
    for i in range(len(rings)):
        ri = set(rings[i])
        for j in range(i + 1, len(rings)):
            if len(ri & set(rings[j])) >= 2:
                return True
    return False


def complexity(smiles: str) -> float:
    """BertzCT 分子复杂度指数; SMILES 无法解析时记 0。"""
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0.0
    try:
        return float(Descriptors.BertzCT(mol))
    except Exception:
        return 0.0


def heavy_atoms(smiles: str) -> int:
    """重原子数(不含氢); SMILES 无法解析时记 0。"""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return int(mol.GetNumHeavyAtoms())


def _normalized(name: str, lang: str) -> str:
    """按 benchmark 判分口径归一: EN 小写/统一连字符, ZH 统一括号种类。"""
    from namepredict.constants import normalize_en, normalize_zh

    s = name or ""
    return normalize_zh(s) if lang == "zh" else normalize_en(s)


def similarity(pred: str, gold: str, lang: str) -> float | None:
    """预测名 ↔ 金标名 的字符级相似度 0..1; 无金标名返回 None(不适用)。"""
    if not (gold or "").strip():
        return None
    a = _normalized(pred, lang)
    b = _normalized(gold, lang)
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a, b, autojunk=False).ratio()


def row_derived(smiles: str, pred_en: str, pred_zh: str, gold_en: str, gold_zh: str) -> dict[str, Any]:
    """一行的派生字段字典，可直接 **展开进预览行。"""
    return {
        "fused": is_fused(smiles),
        "cx": complexity(smiles),
        "hac": heavy_atoms(smiles),
        "sim_en": similarity(pred_en, gold_en, "en"),
        "sim_zh": similarity(pred_zh, gold_zh, "zh"),
    }


def row_derived_of(row: dict[str, Any]) -> dict[str, Any]:
    """同上，但入参是已生成好的预览行(补算旧缓存时用)。"""
    return row_derived(
        str(row.get("s") or ""),
        str(row.get("en") or ""),
        str(row.get("zh") or ""),
        str(row.get("ge") or ""),
        str(row.get("gz") or ""),
    )


def has_derived(rows: list[dict[str, Any]]) -> bool:
    """rows 是否已带全部派生字段(空 rows 视为已带, 无需补算)。"""
    return not rows or all(k in rows[0] for k in DERIVED_KEYS)
