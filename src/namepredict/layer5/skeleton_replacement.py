"""P-23.3.1 骨架置换（'a'）前缀：骨架碳被杂原子置换后的取代前缀。

组内位次升序逗号连接并加数量前缀（`3,14-dioxa`），组间按 P-23.3.1 的
引用顺序用连字符连接（`4-thia-1-aza`）。环外取代基由上层在其前拼接。
"""
from __future__ import annotations

from namepredict.constants import C, MULT_EN, MULT_ZH, P145_SENIOR
from namepredict.tools.lambda_notation import bonding_number, is_nonstandard, lambda_mark


A_PREFIX_EN = {5: "bora", 7: "aza", 8: "oxa", 14: "sila",
               15: "phospha", 16: "thia", 33: "arsa", 34: "selena", 51: "stiba"}
A_PREFIX_ZH = {5: "硼杂", 7: "氮杂", 8: "氧杂", 14: "硅杂",
               15: "磷杂", 16: "硫杂", 33: "砷杂", 34: "硒杂", 51: "锑杂"}


def prefix_from_chain(mol, chain: list[int], *, lambda_ok: bool = True) -> tuple[str, str] | None:
    """由骨架原子序表产出 'a' 前缀；无杂原子返回两个空串，词表外返回 None。

    lambda_ok 时对键数非标准的骨架杂原子在位次后加 λn（P-23.6.1 / P-15.4.1.3）；
    稠环组分名不标 λ（P-25.6「只在完整环系中标注」），调用方传 False。
    """
    by_z: dict[int, list[int]] = {}
    lam: dict[int, int] = {}
    for i, a in enumerate(chain):
        z = mol.GetAtomWithIdx(a).GetAtomicNum()
        if z == C:
            continue
        by_z.setdefault(z, []).append(a)
        atom = mol.GetAtomWithIdx(a)
        if lambda_ok and is_nonstandard(atom):
            lam[i + 1] = bonding_number(atom)
    if not by_z:
        return "", ""
    locants = {z: [chain.index(a) + 1 for a in atoms] for z, atoms in by_z.items()}
    return skeleton_replacement_prefix(locants, lambda_at=lam) or (None, None)


def skeleton_replacement_prefix(locants_by_z: dict[int, list[int]], marks: str = "",
                                lambda_at: dict[int, int] | None = None) -> tuple[str, str] | None:
    """产出 '4-thia-1-aza' / '4-硫杂-1-氮杂' 型前缀串；词表外元素返回 None。

    marks 为位次撇号后缀（组分式螺环第 k 个组分传 k 个撇号，P-24.6）。
    lambda_at 为 {位次: 键数}，命中者在位次后紧贴 λn（P-23.6.1：λ 置于 'a' 前缀之前）。
    """
    marks_at = lambda_at or {}
    parts_en: list[str] = []
    parts_zh: list[str] = []
    for z in P145_SENIOR:  # P-23.3.1 引用顺序（F>Cl>Br>I>O>S>Se>Te>N>P>…）
        locants = sorted(locants_by_z.get(z) or ())
        if not locants:
            continue
        en, zh = A_PREFIX_EN.get(z), A_PREFIX_ZH.get(z)
        mult_en, mult_zh = MULT_EN.get(len(locants)), MULT_ZH.get(len(locants))
        if en is None or zh is None or mult_en is None or mult_zh is None:
            return None  # 词表外元素或数量超出词表：不得臆造前缀
        loc_en = ",".join(f"{x}{marks}" + (lambda_mark(marks_at[x], en=True) if x in marks_at else "")
                          for x in locants)
        loc_zh = ",".join(f"{x}{marks}" + (lambda_mark(marks_at[x]) if x in marks_at else "")
                          for x in locants)
        # 数量词尾 'a' 仅在后接元素名以 'a' 开头时省略（tetraza 对 tetraoxa）
        if mult_en.endswith("a") and en.startswith("a"):
            mult_en = mult_en[:-1]
        parts_en.append(f"{loc_en}-{mult_en}{en}")
        parts_zh.append(f"{loc_zh}-{mult_zh}{zh}")
    return ("-".join(parts_en), "-".join(parts_zh)) if parts_en else None
