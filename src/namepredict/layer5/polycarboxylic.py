"""L5 assembly for supported open-chain polycarboxylic parents."""
from __future__ import annotations
from namepredict.constants import MULT_EN, MULT_ZH

def benzene_polycarboxylic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    parent, locs = numbered.get("parent") or {}, numbered.get("cooh_locants") or []
    count = parent.get("acid_count")
    if count not in (2, 3) or len(locs) != count: return None
    en_mult, zh_mult = MULT_EN[count], MULT_ZH[count]
    loc = _pair_loc_str(locs)
    return f"benzene-{loc}-{en_mult}carboxylic acid", f"苯-{loc}-{zh_mult}羧酸"


def _pair_loc_str(locs: list[int]) -> str:
    return ",".join(map(str, locs))

def cycloalkane_polycarboxylic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    parent, locs = numbered.get("parent") or {}, numbered.get("cooh_locants") or []
    count, stem = parent.get("acid_count"), (parent.get("stem_en"), parent.get("stem_zh"))
    if count not in (2, 3) or len(locs) != count or not all(stem): return None
    loc, mult = _pair_loc_str(locs), (MULT_EN[count], MULT_ZH[count])
    en, zh = f"{stem[0]}-{loc}-{mult[0]}carboxylic acid", f"{stem[1]}-{loc}-{mult[1]}甲酸"
    relative = parent.get("relative_stereo_prefix")
    if relative: return f"{relative}-{en}", f"{'顺' if relative == 'cis' else '反'}-{zh}"
    relative = parent.get("relative_stereo_locants")
    return (f"{relative}-{en}", f"{relative}-{zh}") if relative else (en, zh)


def polycarboxylic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """Assemble only the facts materialized by Layer 4."""
    parent, locs = numbered.get("parent") or {}, numbered.get("cooh_locants") or []
    count, stems = parent.get("acid_count"), (numbered.get("stem_en"), numbered.get("stem_zh"))
    if not count or len(locs) != count or not all(stems): return None
    mult = MULT_EN.get(count), MULT_ZH.get(count)
    if not all(mult): return None
    return _polycarboxylic_text(stems, locs, mult, bool(parent.get("anion")))

def _polycarboxylic_text(stems, locs, mult, anion: bool) -> tuple[str, str]:
    en, zh = stems
    suffix = f"{mult[0]}carboxylate" if anion else f"{mult[0]}carboxylic acid"
    zsuffix = f"{mult[1]}羧酸根" if anion else f"{mult[1]}羧酸"
    loc = _pair_loc_str(locs)
    return f"{en}-{loc}-{suffix}", f"{zh}-{loc}-{zsuffix}"
