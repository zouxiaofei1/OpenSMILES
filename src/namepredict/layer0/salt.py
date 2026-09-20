"""L0 盐解离：金属阳离子/卤化物/氢卤酸反离子；返回有机 mol 与盐元数据供 L5。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.constants import (
    HALIDE_EN, HALIDE_HX_EN, HALIDE_HX_ZH, HALIDE_ZH, HALO_Z,
    METAL_ION_EN, METAL_ION_ZH, MULT_EN, MULT_ZH,
)


def _single_atom(mol: Mol):
    """单原子片段返回该原子，否则 None。"""
    if mol.GetNumAtoms() != 1:
        return None
    return mol.GetAtomWithIdx(0)

def _metal_ion(mol: Mol) -> str | None:
    """单原子金属阳离子 → 英文金属名（P-71.2）；非金属或多原子返回 None。"""
    a = _single_atom(mol)
    if a is None or a.GetFormalCharge() < 1:
        return None
    return METAL_ION_EN.get(a.GetAtomicNum())

def _halide_ion(mol: Mol) -> int | None:
    """单原子卤素阴离子 [X-] → 原子序数；否则 None。"""
    a = _single_atom(mol)
    if a is None or a.GetFormalCharge() != -1 or a.GetAtomicNum() not in HALO_Z:
        return None
    return a.GetAtomicNum() if a.GetTotalNumHs() == 0 else None

def _hydrogen_halide(mol: Mol) -> int | None:
    """中性卤化氢 HX（卤素带 1 个 H）→ 原子序数；否则 None。"""
    a = _single_atom(mol)
    if a is None or a.GetFormalCharge() != 0 or a.GetAtomicNum() not in HALO_Z:
        return None
    return a.GetAtomicNum() if a.GetTotalNumHs() == 1 else None


def _mult_word(words: dict[int, str], n: int, mult: dict) -> str | None:
    """按份数加数量前缀（n=1 不加）；缺词表返回 None。"""
    if n == 1:
        return words.get(1) or None
    m = mult.get(n)
    base = words.get(1)
    return f"{m}{base}" if m and base else None


def _from_frags(frags: tuple[Mol, ...]) -> tuple[Mol, dict] | None:
    """从片段中提取单一有机分子与盐元数据；不满足则返回 None。"""
    metals: list[str] = []
    anions: list[int] = []   # 卤素阴离子 X⁻ 的原子序数
    hx: list[int] = []       # 中性卤化氢 HX 的原子序数
    organics: list[Mol] = []
    for f in frags:
        m = _metal_ion(f)
        if m is not None:
            metals.append(m)
            continue
        z = _halide_ion(f)
        if z is not None:
            anions.append(z)
            continue
        z = _hydrogen_halide(f)
        if z is not None:
            hx.append(z)
            continue
        organics.append(f)
    if not organics or len({Chem.MolToSmiles(f, isomericSmiles=False) for f in organics}) != 1:
        return None  # 无有机片段，或存在多种有机片段：非简单盐
    organic = organics[0]
    n_org = len(organics)
    charge = sum(a.GetFormalCharge() for a in organic.GetAtoms())
    if metals:  # 金属阳离子：金属种类须唯一，且不与卤素反离子共存
        if anions or hx or len(set(metals)) != 1 or charge >= 0:
            return None
        meta = {"metal": metals[0], "metal_zh": METAL_ION_ZH[metals[0]],
                "n_metal": len(metals), "n_org": n_org}
    elif anions:  # 卤素阴离子 X⁻：有机物须为阳离子（P-71.2 有机阳离子 + 卤离子）
        if hx or len(set(anions)) != 1 or charge <= 0:
            return None
        en = _mult_word({1: HALIDE_EN[anions[0]]}, len(anions), MULT_EN)
        zh = _mult_word({1: HALIDE_ZH[anions[0]]}, len(anions), MULT_ZH)
        if en is None or zh is None:
            return None
        meta = {"halide": en, "halide_zh": zh, "n_org": n_org}
    else:  # 中性卤化氢 HX：有机物须为中性碱（氢卤酸盐）
        if len(set(hx)) != 1 or charge != 0:
            return None
        en = _mult_word({1: HALIDE_HX_EN[hx[0]]}, len(hx), MULT_EN)
        zh = _mult_word({1: HALIDE_HX_ZH[hx[0]]}, len(hx), MULT_ZH)
        if en is None or zh is None:
            return None
        meta = {"acid_salt": en, "acid_salt_zh": zh, "n_org": n_org}
    return organic, meta


def dissociate_salt(mol: Mol) -> tuple[Mol, dict]:
    """返回有机 mol 与盐元数据；非简单盐时 meta 为空。"""
    if len(Chem.GetMolFrags(mol)) < 2:
        return mol, {}
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    hit = _from_frags(frags)
    return hit if hit is not None else (mol, {})
