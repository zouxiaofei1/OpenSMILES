"""Alkenoic / alkenedioic acid name assembly (E/Z via stereo_ez)."""
from __future__ import annotations

from namepredict.layer5.stems import (
    ALKANE_EN, ALKANE_ZH, zh_stem,
)
from namepredict.layer5.stereo_ez import (
    _ez_for_alkenol,
    _ez_prefix,
    ez_for_parent,
)


def _unsat_acid_pair(n, locant, ez, en_sfx, zh_sfx, min_n=2) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or locant is None or n < min_n:
        return None
    return f"{ez}{en[:-3]}-{locant}-{en_sfx}", f"{ez}{zh_stem(zh)}-{locant}-{zh_sfx}"


def _alkenedioic_names(n: int, locant: int | None, ez: str) -> tuple[str, str] | None:
    return _unsat_acid_pair(n, locant, ez, "enedioic acid", "烯二酸", 3)


def _ene_mult_diacid(k: int) -> tuple[str, str]:
    en = {
        2: "dienedioic acid", 3: "trienedioic acid", 4: "tetraenedioic acid",
    }.get(k, "")
    zh = {2: "二烯二酸", 3: "三烯二酸", 4: "四烯二酸"}.get(k, "")
    return en, zh


def _polyalkenedioic_names(n: int, locs, ez: str) -> tuple[str, str] | None:
    """hexa-2,4-dienedioic acid / 己-2,4-二烯二酸."""
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    me, mz = _ene_mult_diacid(len(locs or []))
    if not en or not zh or not me or not locs or len(locs) < 2 or n < 4:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"{ez}{en[:-3]}a-{loc}-{me}", f"{ez}{zh_stem(zh)}-{loc}-{mz}"


def alkenedioic_names(n: int, numbered: dict) -> tuple[str, str] | None:
    locs = numbered.get("ene_locants")
    if locs and len(locs) >= 2:
        return _polyalkenedioic_names(n, locs, _ez_for_alkenol(numbered))
    return _alkenedioic_names(n, numbered.get("ene_locant"), _ez_prefix(numbered))


def _alkenoic_acid_names(n: int, locant: int | None, ez: str) -> tuple[str, str] | None:
    return _unsat_acid_pair(n, locant, ez, "enoic acid", "烯酸")


def _ene_mult_acid(k: int) -> tuple[str, str]:
    en = {
        2: "dienoic acid", 3: "trienoic acid", 4: "tetraenoic acid",
        5: "pentaenoic acid", 6: "hexaenoic acid",
    }.get(k, "")
    zh = {
        2: "二烯酸", 3: "三烯酸", 4: "四烯酸", 5: "五烯酸", 6: "六烯酸",
    }.get(k, "")
    return en, zh


def _polyalkenoic_acid_names(n: int, locs, ez: str) -> tuple[str, str] | None:
    """octadeca-9,11-dienoic acid / 十八-9,11-二烯酸."""
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    me, mz = _ene_mult_acid(len(locs or []))
    if not en or not zh or not me or not locs or len(locs) < 2 or n < 4:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"{ez}{en[:-3]}a-{loc}-{me}", f"{ez}{zh_stem(zh)}-{loc}-{mz}"


def alkenoic_acid_names(n: int, numbered: dict) -> tuple[str, str] | None:
    locs = numbered.get("ene_locants")
    if locs and len(locs) >= 2:
        return _polyalkenoic_acid_names(n, locs, _ez_for_alkenol(numbered))
    return _alkenoic_acid_names(
        n, numbered.get("ene_locant"), _ez_prefix(numbered),
    )


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


def _alkynoic_pair(n, locant, en_sfx, zh_sfx, omit=False) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or n < 2:
        return None
    if omit or locant is None:
        return f"{en[:-3]}{en_sfx}", f"{zh_stem(zh)}{zh_sfx}"
    return f"{en[:-3]}-{locant}-{en_sfx}", f"{zh_stem(zh)}-{locant}-{zh_sfx}"


def alkynoic_acid_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """propynoic acid / but-3-ynoic acid (omit locant when n≤3)."""
    return _alkynoic_pair(
        n, numbered.get("yne_locant"), "ynoic acid", "炔酸",
        numbered.get("omit_yne_locant", False),
    )


def alkynamide_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """propynamide / but-3-ynamide (omit locant when n≤3)."""
    return _alkynoic_pair(
        n, numbered.get("yne_locant"), "ynamide", "炔酰胺",
        numbered.get("omit_yne_locant", False),
    )


def _yne_carbonyl(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "acid" and _has_yne(numbered):
        return alkynoic_acid_names(n, numbered)
    if kind == "amide" and _has_yne(numbered):
        return alkynamide_names(n, numbered)
    return None


def _ene_carbonyl(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if not _has_ene(numbered):
        return None
    if kind == "diacid":
        return alkenedioic_names(n, numbered)
    if kind == "acid":
        return alkenoic_acid_names(n, numbered)
    return alkenamide_names(n, numbered) if kind == "amide" else None


def unsat_carbonyl_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Dispatch unsat acid / diacid / amide stems (yne then ene)."""
    return _yne_carbonyl(kind, n, numbered) or _ene_carbonyl(kind, n, numbered)


def _ene_mult_al(k: int) -> tuple[str, str]:
    en = {
        2: "dienal", 3: "trienal", 4: "tetraenal",
        5: "pentaenal", 6: "hexaenal",
    }.get(k, "")
    zh = {
        2: "二烯醛", 3: "三烯醛", 4: "四烯醛", 5: "五烯醛", 6: "六烯醛",
    }.get(k, "")
    return en, zh


def _polyalkenal_names(n: int, locs, ez: str) -> tuple[str, str] | None:
    """penta-2,4-dienal / 戊-2,4-二烯醛."""
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    me, mz = _ene_mult_al(len(locs or []))
    if not en or not zh or not me or not locs or len(locs) < 2 or n < 4:
        return None
    loc = ",".join(str(x) for x in locs)
    return f"{ez}{en[:-3]}a-{loc}-{me}", f"{ez}{zh_stem(zh)}-{loc}-{mz}"


def alkenal_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """non-2-enal / penta-2,4-dienal with E/Z when stereo defined."""
    locs = numbered.get("ene_locants")
    if locs and len(locs) >= 2:
        return _polyalkenal_names(n, locs, ez_for_parent(numbered))
    return _unsat_acid_pair(
        n, numbered.get("ene_locant"), ez_for_parent(numbered), "enal", "烯醛",
    )


def alkynal_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """propynal / but-3-ynal (omit locant when n≤3)."""
    return _alkynoic_pair(
        n, numbered.get("yne_locant"), "ynal", "炔醛",
        numbered.get("omit_yne_locant", False),
    )


def alkenenitrile_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """octadec-9-enenitrile / (Z)-… with E/Z when stereo defined."""
    return _unsat_acid_pair(
        n, numbered.get("ene_locant"), ez_for_parent(numbered), "enenitrile", "烯腈",
    )


def alkynenitrile_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """propynenitrile / but-3-ynenitrile (omit locant when n≤3)."""
    return _alkynoic_pair(
        n, numbered.get("yne_locant"), "ynenitrile", "炔腈",
        numbered.get("omit_yne_locant", False),
    )


def alkenoate_names(n, locant, parent, ez="") -> tuple[str, str] | None:
    from namepredict.layer5.stems import ester_alkoxy_pair
    alkyl = ester_alkoxy_pair(parent)
    stem = _unsat_acid_pair(n, locant, "", "enoate", "烯酸")
    if not alkyl or not stem:
        return None
    return f"{alkyl[0]} {ez}{stem[0]}", f"{ez}{stem[1]}{alkyl[1]}酯"


def alkynoate_names(n, locant, parent, omit=False) -> tuple[str, str] | None:
    from namepredict.layer5.stems import ester_alkoxy_pair
    alkyl = ester_alkoxy_pair(parent)
    stem = _alkynoic_pair(n, locant, "ynoate", "炔酸", omit)
    if not alkyl or not stem:
        return None
    return f"{alkyl[0]} {stem[0]}", f"{stem[1]}{alkyl[1]}酯"


def _ene_mult_ol(k: int) -> tuple[str, str]:
    en = {2: "diene", 3: "triene", 4: "tetraene", 5: "pentaene"}.get(k, "")
    zh = {2: "二烯", 3: "三烯", 4: "四烯", 5: "五烯"}.get(k, "")
    return en, zh


def _alkenol_names(n, ene_loc, oh_loc, ez="") -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or ene_loc is None or oh_loc is None:
        return None
    return (
        f"{ez}{en[:-3]}-{ene_loc}-en-{oh_loc}-ol",
        f"{ez}{zh_stem(zh)}-{ene_loc}-烯-{oh_loc}-醇",
    )


def _polyalkenol_names(n, locs, oh_loc, ez="") -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    me, mz = _ene_mult_ol(len(locs or []))
    if not en or not zh or not me or oh_loc is None or not locs or len(locs) < 2:
        return None
    loc, en_m = ",".join(str(x) for x in locs), me[:-1] if me.endswith("e") else me
    return (
        f"{ez}{en[:-3]}a-{loc}-{en_m}-{oh_loc}-ol",
        f"{ez}{zh_stem(zh)}-{loc}-{mz}-{oh_loc}-醇",
    )


def alkenol_from(n: int, numbered: dict) -> tuple[str, str] | None:
    ez = _ez_for_alkenol(numbered)
    locs = numbered.get("ene_locants")
    if locs and len(locs) >= 2:
        return _polyalkenol_names(n, locs, numbered.get("oh_locant"), ez)
    return _alkenol_names(n, numbered.get("ene_locant"), numbered.get("oh_locant"), ez)


def _alkynol_names(n, yne_loc, oh_loc) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or yne_loc is None or oh_loc is None:
        return None
    return (
        f"{en[:-3]}-{yne_loc}-yn-{oh_loc}-ol",
        f"{zh_stem(zh)}-{yne_loc}-炔-{oh_loc}-醇",
    )


def alkynol_from(n: int, numbered: dict) -> tuple[str, str] | None:
    """but-3-yn-1-ol / prop-2-yn-1-ol (OH + yne locants retained)."""
    return _alkynol_names(n, numbered.get("yne_locant"), numbered.get("oh_locant"))


def _alkenone_names(n, ene_loc, one_loc, ez="") -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or ene_loc is None or one_loc is None:
        return None
    return (
        f"{ez}{en[:-3]}-{ene_loc}-en-{one_loc}-one",
        f"{ez}{zh_stem(zh)}-{ene_loc}-烯-{one_loc}-酮",
    )


def alkenone_from(n: int, numbered: dict) -> tuple[str, str] | None:
    """but-3-en-2-one / 丁-3-烯-2-酮 with E/Z when stereo defined."""
    return _alkenone_names(
        n, numbered.get("ene_locant"), numbered.get("ketone_locant"),
        ez_for_parent(numbered) or "",
    )


def _alkynone_names(n, yne_loc, one_loc) -> tuple[str, str] | None:
    en, zh = ALKANE_EN.get(n), ALKANE_ZH.get(n)
    if not en or not zh or yne_loc is None or one_loc is None:
        return None
    return (
        f"{en[:-3]}-{yne_loc}-yn-{one_loc}-one",
        f"{zh_stem(zh)}-{yne_loc}-炔-{one_loc}-酮",
    )


def alkynone_from(n: int, numbered: dict) -> tuple[str, str] | None:
    """but-3-yn-2-one / 丁-3-炔-2-酮 (ketone + yne locants retained)."""
    return _alkynone_names(n, numbered.get("yne_locant"), numbered.get("ketone_locant"))


def ester_or_alkenoate(n: int, numbered: dict) -> tuple[str, str] | None:
    """Ester stem; unsat when yne or ene fields present."""
    parent = numbered.get("parent") or {}
    if _has_yne(numbered):
        omit = numbered.get("omit_yne_locant", False)
        return alkynoate_names(n, numbered.get("yne_locant"), parent, omit)
    if _has_ene(numbered):
        return alkenoate_names(
            n, numbered.get("ene_locant"), parent, _ez_prefix(numbered),
        )
    return None
