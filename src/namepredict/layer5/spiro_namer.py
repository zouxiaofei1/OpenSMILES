"""P-24 螺环母体名：螺词头 + von Baeyer 螺描述符 + 烷词干 + 骨架置换 'a' 前缀。

同时给出两种词干形态：FG 分支要带 "ane/烷" 的完整名（后缀去 e 用），
链式词干引擎要裸词干（自行拼 "ane/烷" 或 "a…-ene"）。L5 不得 import L2，
因此 node 只按鸭子类型读 descriptor / descriptor_superscripts / free_spiro_atoms。
"""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.layer5.skeleton_replacement import prefix_from_chain
from namepredict.layer5.stems import alkane_en, alkane_zh


def spiro_multiplier(n_spiro: int) -> tuple[str, str] | None:
    """螺原子数词头：1 用 spiro/螺，2 起用 dispiro/二螺 型计数词（不是 bispiro）。"""
    if n_spiro == 1:
        return "spiro", "螺"
    en, zh = MULT_EN.get(n_spiro), MULT_ZH.get(n_spiro)
    return (f"{en}spiro", f"{zh}螺") if en and zh else None


def spiro_descriptor_str(descriptor: tuple[int, ...], superscripts: tuple[int, ...]) -> str:
    """螺描述符串：重访螺原子的位次紧接段长之后（判分口径无 ^{}）。"""
    parts: list[str] = []
    for i, d in enumerate(descriptor):
        sup = superscripts[i] if i < len(superscripts) else 0
        parts.append(f"{d}{sup}" if sup else str(d))
    return ".".join(parts)


def spiro_parent_names(mol, node, chain: list[int]):
    """螺环名 → ((完整英, 完整中), (裸词干英, 裸词干中))；不可组装返回 None。"""
    if not chain or not node.descriptor:
        return None
    mult = spiro_multiplier(len(node.free_spiro_atoms))
    if mult is None:
        return None
    a_en, a_zh = prefix_from_chain(mol, chain)
    if a_en is None or a_zh is None:
        return None
    stem_en, stem_zh = alkane_en(len(chain)), alkane_zh(len(chain))
    if not stem_en or not stem_zh:
        return None
    desc = spiro_descriptor_str(node.descriptor, node.descriptor_superscripts)
    body_en, body_zh = f"{mult[0]}[{desc}]", f"{mult[1]}[{desc}]"
    return ((f"{a_en}{body_en}{stem_en}", f"{a_zh}{body_zh}{stem_zh}"),
            (f"{a_en}{body_en}{stem_en[:-3]}", f"{a_zh}{body_zh}{stem_zh[:-1]}"))
