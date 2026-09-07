"""L5 整分子磷酸（kind=phosphate）命名 worker：由 numbered 母体计数 + o_side 臂组装双语名。"""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.layer5.stems import _metal_en_prefix, _metal_zh_suffix


def _tail_en(h: int) -> str:
    """英文磷酸词尾（含酸式 H 词）：dihydrogen/hydrogen phosphate 或 phosphate。"""
    if h == 2:
        return "dihydrogen phosphate"
    if h == 1:
        return "hydrogen phosphate"
    return "phosphate"


def _hyd_zh(h: int) -> str:
    """中文酸式 H 词：二氢 / 氢（h=0 空串）。"""
    if h == 2:
        return "二氢"
    if h == 1:
        return "氢"
    return ""


def _arm_zh_root(zh: str) -> str:
    """臂中文词干：简单单字根（甲基→甲/苯基→苯）去基；复合/带位次/立体的臂词干保留。"""
    if (zh.endswith("基") and "-" not in zh and not zh.startswith("(")
            and not any(ch.isdigit() for ch in zh)):
        return zh[:-1]
    return zh


def _group_arms(arms: list[dict]) -> list[tuple[str, str, int]]:
    """按臂 radical 名分组（相同 en 合并计数），按 en 字母序返回 [(en, zh, m), …]。"""
    table: dict[str, list] = {}
    for s in arms:
        en = (s.get("en") or "").strip()
        zh = (s.get("zh") or "").strip()
        if not en:
            continue
        row = table.get(en)
        if row is None:
            table[en] = [en, zh, 1]
        else:
            row[2] += 1
    return [tuple(table[k]) for k in sorted(table)]


def _salt_name(n_oh: int, salt: dict) -> tuple[str, str] | None:
    """纯磷酸碱金属盐（k=0, d>0）：metal [dihydrogen|hydrogen] phosphate / 磷酸[二氢|氢]{金属}。"""
    metal_en = _metal_en_prefix(salt)
    metal_zh = _metal_zh_suffix(salt)
    if not metal_en or not metal_zh:
        return None
    return (f"{metal_en} {_tail_en(n_oh)}", f"磷酸{_hyd_zh(n_oh)}{metal_zh}")


def phosphate_names(numbered: dict) -> tuple[str, str] | None:
    """组装磷酸整名；不支持形态返回 None。"""
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "phosphate":
        return None
    n_oh = int(parent.get("n_oh") or 0)
    n_om = int(parent.get("n_om") or 0)
    salt = parent.get("salt_meta") or {}
    arms = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]

    if n_om > 0:
        if arms or not salt.get("metal"):
            return None
        return _salt_name(n_oh, salt)
    if not arms:
        if n_oh == 3 and not salt:
            return ("phosphoric acid", "磷酸")
        return None
    if salt:
        return None  # 中性酯不应带金属（盐+酯混合不在本轮）
    # 中性磷酸酯：k = len(arms) > 0，h = 3 − k
    groups = _group_arms(arms)
    if not groups:
        return None
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for en, zh, m in groups:
        if m > 1:
            m_en, m_zh = MULT_EN.get(m), MULT_ZH.get(m)
            en_w = f"{m_en}{en}" if m_en else None
            zh_w = f"{m_zh}{_arm_zh_root(zh)}" if m_zh else None
        else:
            en_w, zh_w = en, _arm_zh_root(zh)
        if not en_w or not zh_w:
            return None
        en_parts.append(en_w)
        zh_parts.append(zh_w)
    tail = _tail_en(n_oh)
    en = f"{' '.join(en_parts)} {tail}"
    zh = f"磷酸{_hyd_zh(n_oh)}{''.join(zh_parts)}酯"
    return (en, zh)
