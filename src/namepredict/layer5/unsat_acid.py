"""Alkenoic / alkenedioic acid name assembly (E/Z via layer5.stereo)."""
from __future__ import annotations

from namepredict.layer5.stems import (
    ALKANE_EN, ALKANE_ZH, zh_stem,
)
from namepredict.layer5.stereo import _ez_prefix


def _unsat_acid_pair(n, locant, ez, en_sfx, zh_sfx, min_n=2) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or locant is None or n < min_n:
        return None
    return f"{ez}{en[:-3]}-{locant}-{en_sfx}", f"{ez}{zh_stem(zh)}-{locant}-{zh_sfx}"





def alkenamide_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """acrylamide (retained) / (E)-but-2-enamide with E/Z when stereo defined."""
    if n == 3 and numbered.get("ene_locant") == 2 and numbered.get("name_mode") != "pin":
        ez = _ez_prefix(numbered)
        return f"{ez}acrylamide", f"{ez}丙烯酰胺"
    return _unsat_acid_pair(
        n, numbered.get("ene_locant"), _ez_prefix(numbered), "enamide", "烯酰胺",
    )


def _has_ene(numbered: dict) -> bool:
    p = numbered.get("parent") or {}
    return bool(
        numbered.get("ene_locant") or numbered.get("ene_locants")
        or p.get("double_bond") or p.get("double_bonds")
    )


def _has_yne(numbered: dict) -> bool:
    p = numbered.get("parent") or {}
    return bool(numbered.get("yne_locant") or p.get("triple_bond"))


