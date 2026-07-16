"""L5 assembly for supported open-chain polycarboxylic parents."""
from __future__ import annotations

from namepredict.layer5.stems import ALKANE_EN, ALKANE_ZH

MULT_EN = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
MULT_ZH = {2: "二", 3: "三", 4: "四", 5: "五", 6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}


def _alkane_names(n: int) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    return (en, zh) if en and zh else None


def _pair_loc_str(locs: list[int]) -> str:
    return ",".join(str(x) for x in locs)

def polycarboxylic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    locs = numbered.get("cooh_locants") or []
    count = (numbered.get("parent") or {}).get("acid_count")
    return _polycarboxylic_pair(n, locs, count, bool((numbered.get("parent") or {}).get("anion")))
def _polycarboxylic_pair(n: int, locs: list[int], count: int | None, anion: bool) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain or not count or len(locs) != count: return None
    mult = MULT_EN.get(count), MULT_ZH.get(count)
    if not all(mult): return None
    return _polycarboxylic_text(plain, locs, mult, anion)
def _polycarboxylic_text(plain, locs, mult, anion: bool) -> tuple[str, str]:
    en, zh = plain; loc = _pair_loc_str(locs); suffix = f"{mult[0]}carboxylate" if anion else f"{mult[0]}carboxylic acid"
    zsuffix = f"{mult[1]}羧酸根" if anion else f"{mult[1]}羧酸"
    return f"{en}-{loc}-{suffix}", f"{zh}-{loc}-{zsuffix}"
