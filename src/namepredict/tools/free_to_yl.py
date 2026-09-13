"""free 母体名 → P-29 -yl 取代基形式的双语转换（含官能团前缀化）。"""

from __future__ import annotations

import re

from namepredict.constants import MONONUCLEAR_YL, zh_bridge_root

_MULT_OL_SUF = re.compile(r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)ol\b")  # 官能团后缀 → 前缀转换

_MONONUCLEAR_NAMES = ("oxidane", "azane", "sulfane", "sulfinyl", "sulfonyl", "imine")  # 本模块转换的单核母体氢化物（P-15.4.1 表 2.1 → 表 1.5 'a' 前缀体系，数据见 constants.MONONUCLEAR_HYDRIDES）；P 酰基（phosphoryl/phosphanyl）走 L5 的 P 专用管线，不经此表。


def _anilino_en(base: str, yl: str) -> str:
    """N-苯基（可带环取代基）的 azane 去氢：phenyl→anilino、4-chlorophenyl→4-chloroanilino（P-62.2.1.1：phenylamino = anilino*）。"""
    return base[: -len("phenyl")] + "anilino" if yl == "amino" and base.endswith("phenyl") else base + yl


def _mononuclear_en(en: str) -> str | None:
    """单核氢化物 free 名 → 去氢取代基名：ethyl-oxidane → ethyloxy、phenyl-azane → anilino。"""
    for en_suf in _MONONUCLEAR_NAMES:
        if en.endswith("-" + en_suf):
            return _anilino_en(en[: -len(en_suf) - 1], MONONUCLEAR_YL[en_suf][1])
    return None


def _mononuclear_zh(zh: str) -> str | None:
    """中文组装名去氢：乙基-氧化烷 → 乙氧基、苯基-氮烷 → 苯胺基、4-氯苯基-氮烷 → 4-氯苯胺基。"""
    for en_suf in _MONONUCLEAR_NAMES:
        zh_suf, _, zy = MONONUCLEAR_YL[en_suf]
        if zh.endswith("-" + zh_suf):
            base = zh[: -len(zh_suf) - 1]
            if zy == "氨基" and base.endswith("苯基"):
                return base[: -len("苯基")] + "苯胺基"
            return zh_bridge_root(base) + zy
    return None

def _try_fg_prefix(en: str, zh: str) -> tuple[str, str] | None:
    """尝试将官能团后缀名转换为取代基前缀形式。"""
    for en_fn, zh_fn in [
        (_mononuclear_en, _mononuclear_zh),
    ]:
        en_out = en_fn(en)
        zh_out = zh_fn(zh) if en_out else None
        if en_out and zh_out:
            return en_out, zh_out
    return None

def free_to_yl(
    en: str, zh: str, attach_locant: int, *, paren: bool = True,
) -> tuple[str, str, bool]:
    """在键合位次处把 free 母体名转 P-29 -yl 双语形式（P-63.2.2 醇→烷氧基、P-63.2.1 硫醇→烷基硫基、P-62.2 1° 胺→烷基氨基）。"""
    fg = _try_fg_prefix(en, zh)
    if fg is not None:
        need_paren = fg[0].endswith("amino") and fg[0] != "amino"  # P-29.3.6：复合前缀（methylamino=CH3-NH-，非普通 amino）需括号与两个独立取代基区分。
        need_paren = need_paren or (fg[0].endswith("anilino") and fg[0] != "anilino")  # 带环取代基的 anilino（4-chloroanilino）与 …phenylamino 同理需括号（gold：(4-chloroanilino)benzoic acid）；裸 anilino 免括。
        return fg[0], fg[1], need_paren
    return None

