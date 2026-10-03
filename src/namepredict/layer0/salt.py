"""L0 盐解离：金属阳离子/卤化物/氢卤酸反离子；返回有机 mol 与盐元数据供 L5。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol

from namepredict.constants import (
    Cl, HALIDE_EN, HALIDE_HX_EN, HALIDE_HX_ZH, HALIDE_ZH, HALO_Z,
    METAL_ION_EN, METAL_ION_ZH, MULT_EN, MULT_ZH,
)

# 多原子无机反离子（P-71.2）：片段规范 SMILES → (en, zh)；阳离子并入金属通道、阴离子并入卤离子通道。
POLY_CATION = {"[NH4+]": ("azanium", "铵")}
POLY_ANION = {
    "[OH-]": ("hydroxide", "氢氧化物"),
    "O=[N+]([O-])[O-]": ("nitrate", "硝酸盐"),
    "[O-][Cl+3]([O-])([O-])[O-]": ("perchlorate", "高氯酸盐"),
    "[N-]=[N+]=[N-]": ("azide", "叠氮化物"),
    "F[P-](F)(F)(F)(F)F": ("hexafluorophosphate", "六氟磷酸盐"),
}

# 卤离子名 → 氢卤酸加合物名（P-71.3）：中性有机碱 + X⁻ 成盐时按 HX 命名（吡啶盐酸盐）
_HX_BY_HALIDE_EN = {HALIDE_EN[z]: (HALIDE_HX_EN[z], HALIDE_HX_ZH[z]) for z in HALIDE_EN}


def _frag_role(mol: Mol) -> tuple[str, str, str] | None:
    """单片段盐角色 → ("metal"|"halide"|"hx", en, zh)；非盐离子或有机片段返回 None。

    金属/多原子阳离子与卤素/多原子阴离子取词表名（P-71.2），卤化氢取加合物名（P-71.3）。
    """
    smi = Chem.MolToSmiles(mol)
    if smi in POLY_CATION:  # 多原子阳离子（NH4⁺ 等）
        en, zh = POLY_CATION[smi]
        return "metal", en, zh
    if smi in POLY_ANION:  # 多原子阴离子（OH⁻/ClO4⁻/NO3⁻/N3⁻ 等）
        en, zh = POLY_ANION[smi]
        return "halide", en, zh
    if mol.GetNumAtoms() != 1:
        return None
    a = mol.GetAtomWithIdx(0)
    z, q = a.GetAtomicNum(), a.GetFormalCharge()
    if q >= 1 and z in METAL_ION_EN:
        en = METAL_ION_EN[z]
        return "metal", en, METAL_ION_ZH[en]
    if z not in HALO_Z:
        return None
    if q == -1 and a.GetTotalNumHs() == 0:
        return "halide", HALIDE_EN[z], HALIDE_ZH[z]
    if q == 0 and a.GetTotalNumHs() == 1:
        return "hx", HALIDE_HX_EN[z], HALIDE_HX_ZH[z]
    return None


def _mult_word(base: str | None, n: int, mult: dict) -> str | None:
    """按份数给基名加数量前缀（n=1 不加）；缺词表返回 None。"""
    if not base:
        return None
    if n == 1:
        return base
    m = mult.get(n)
    return f"{m}{base}" if m else None


def _net_charge(frag: Mol) -> int:
    """片段的净形式电荷。"""
    return sum(a.GetFormalCharge() for a in frag.GetAtoms())


def pair_unique_ions(frags: tuple[Mol, ...]) -> tuple[Mol, Mol, int, int] | None:
    """唯一可配对的有机阴阳离子 → (阳离子片段, 阴离子片段, n_阳, n_阴)；否则 None。

    要求两侧物种各唯一且整体电荷守恒，供 P-77 二元盐名拼装。
    """
    cats = [f for f in frags if _net_charge(f) > 0]
    ans = [f for f in frags if _net_charge(f) < 0]
    if not cats or not ans:
        return None
    if len({Chem.MolToSmiles(f, isomericSmiles=False) for f in cats}) != 1:
        return None
    if len({Chem.MolToSmiles(f, isomericSmiles=False) for f in ans}) != 1:
        return None
    if sum(_net_charge(f) for f in frags) != 0:
        return None
    return cats[0], ans[0], len(cats), len(ans)


def _from_frags(frags: tuple[Mol, ...]) -> tuple[Mol, dict] | None:
    """从片段中提取单一有机分子与盐元数据；不满足则返回 None。"""
    metals: list[tuple[str, str]] = []   # (en, zh) 阳离子
    anions: list[tuple[str, str]] = []   # (en, zh) 阴离子
    hx: list[tuple[str, str]] = []       # (en, zh) 中性卤化氢 HX
    organics: list[Mol] = []
    for f in frags:
        role = _frag_role(f)
        if role is None:
            organics.append(f)
        elif role[0] == "metal":
            metals.append((role[1], role[2]))
        elif role[0] == "halide":
            anions.append((role[1], role[2]))
        else:
            hx.append((role[1], role[2]))
    if not organics or len({Chem.MolToSmiles(f, isomericSmiles=False) for f in organics}) != 1:
        return None  # 无有机片段，或存在多种有机片段：非简单盐
    organic = organics[0]
    n_org = len(organics)
    charge = sum(a.GetFormalCharge() for a in organic.GetAtoms())
    if metals:  # 金属阳离子：金属种类须唯一，且不与卤素反离子共存
        if anions or hx or len(set(metals)) != 1 or charge >= 0:
            return None
        meta = {"metal": metals[0][0], "metal_zh": metals[0][1],
                "n_metal": len(metals), "n_org": n_org}
    elif anions:  # 阴离子 X⁻：有机阳离子/中性有机碱/有机阴离子均可与之成盐（P-71.2）
        if hx or len(set(anions)) != 1:
            return None
        hx_pair = _HX_BY_HALIDE_EN.get(anions[0][0])
        if (len(anions) == 1 and anions[0][0] == HALIDE_EN[Cl]  # 仅单一 Cl⁻（参考数据里多卤/X≠Cl 一律写卤化物名）
                and charge == 0 and hx_pair is not None
                and not any(a.GetFormalCharge() for a in organic.GetAtoms())):
            # 中性有机碱 + Cl⁻（无内盐）：按盐酸盐加合物命名，与 [Cl-]+HX 输入一致（P-71.3）
            hx_en = _mult_word(hx_pair[0], len(anions), MULT_EN)
            hx_zh = _mult_word(hx_pair[1], len(anions), MULT_ZH)
            if hx_en is None or hx_zh is None:
                return None
            return organic, {"acid_salt": hx_en, "acid_salt_zh": hx_zh, "n_org": n_org}
        en = _mult_word(anions[0][0], len(anions), MULT_EN)
        zh = _mult_word(anions[0][1], len(anions), MULT_ZH)
        if en is None or zh is None:
            return None
        meta = {"halide": en, "halide_zh": zh, "n_org": n_org}
    else:  # 中性卤化氢 HX：有机物须为中性碱（氢卤酸盐）
        if len(set(hx)) != 1 or charge != 0:
            return None
        en = _mult_word(hx[0][0], len(hx), MULT_EN)
        zh = _mult_word(hx[0][1], len(hx), MULT_ZH)
        if en is None or zh is None:
            return None
        meta = {"acid_salt": en, "acid_salt_zh": zh, "n_org": n_org}
    return organic, meta


def dissociate_salt(mol: Mol) -> tuple[Mol, dict]:
    """返回有机 mol 与盐元数据；非简单盐时 meta 为空。"""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return mol, {}
    hit = _from_frags(frags)
    return hit if hit is not None else (mol, {})
