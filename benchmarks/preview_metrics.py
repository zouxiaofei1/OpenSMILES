"""预览行的派生指标：稠环判定、RDKit 复杂度、重原子数、预测↔金标相似度与差异。

由预览生成器 (benchmark_preview_parallel) 与 server 的旧缓存补算路径
(routes_benchmark) 共用，保证"新生成的行"与"补算出来的行"字段完全一致。
相似度取归一化后的字符级 ratio，与 Namer 页 gold 卡同口径，
两个页面对同一对名字给出的百分比不会打架。
逐字符差异 (build_gold_diff) 同样归此一处，Namer 页与 benchmark 预览页的
行内高亮由同一个函数产出，不会出现分割位置不一致。

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
    from namepredict.tools.re import normalize_en, normalize_zh

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


def _diff_sides(pred: str, gold: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """逐字符 diff → (预测侧片段, 基准侧片段), 片段形如 {"t": 文本, "same": bool}。

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


# 相似度阈值: 归一化后预测名与基准名达到该相似度, Namer 页才高亮差异。命得基本对时
# 差异才是可看的细节(如 (furan-3-yl) 少一层括号); 错得离谱时整名标红没有信息量。
# benchmark 预览页忽略它(双击即全量高亮), 只作 Namer 页的开关。
GOLD_DIFF_MIN_SIM = 0.70


def build_gold_diff(
    pred_en: str, pred_zh: str, gold: dict[str, Any] | None
) -> dict[str, Any] | None:
    """预测名 ↔ 基准名 的相似度与逐字符差异; 无基准返回 None。

    每种语言各给 {similarity, show, pred, gold}: show 由相似度是否 ≥ 阈值决定,
    前端据此决定要不要高亮; gold 缺该语言名称则跳过该语言。
    """
    if not gold:
        return None
    out: dict[str, Any] = {}
    for lang, pred in (("en", pred_en), ("zh", pred_zh)):
        ref = (gold.get(lang) or "").strip()
        if not ref:
            continue
        pred = (pred or "").strip()
        sim = similarity(pred, ref, lang) or 0.0
        left, right = _diff_sides(pred, ref)
        out[lang] = {
            "similarity": round(sim, 4),
            "show": sim >= GOLD_DIFF_MIN_SIM,
            "pred": left,
            "gold": right,
        }
    return out or None


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
