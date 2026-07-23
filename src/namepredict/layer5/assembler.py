from __future__ import annotations
from namepredict.layer5.stems import (
    ACID_EN, ACID_ZH, ALCOHOL_EN, ALCOHOL_ZH,
    ALDEHYDE_EN, ALDEHYDE_ZH, ALKANE_EN, ALKANE_ZH, ALKOXY_EN, ALKOXY_ZH,
    AMIDE_EN, AMIDE_ZH, ESTER_ACYL_EN, ESTER_ALKYL_EN, ESTER_ALKYL_ZH,
    ETHER_SYM_EN, ETHER_SYM_ZH, NITRILE_EN, NITRILE_ZH, SULFIDE_ALKYL_EN,
    SULFIDE_ALKYL_ZH, SULFIDE_SYM_EN, SULFIDE_SYM_ZH, maybe_anion_names, maybe_metal_salt_names, zh_stem,
)
from namepredict.layer5.benzene_names import (
    benzene_parent_names, benzene_prefix,
    benzoate_parent_names, benzenediamine_names, benzenediol_names, hetero5carboxylic_names,
    join_kind_name, pyridine_kind_names, sat_hetero_carboxylic_names,
)
from namepredict.types import NameResult
from namepredict.constants import MULT_EN, MULT_ZH
def _fail(meta: dict | None = None) -> NameResult: return NameResult(en="", zh="", success=False, source="iupac", meta=meta or {})
def _ok(en: str, zh: str, time_ms: float, source: str) -> NameResult: return NameResult(en=en, zh=zh, success=True, source=source, time_ms=time_ms)
def _pair(en_map: dict, zh_map: dict, n: int) -> tuple[str, str] | None:
    en, zh = en_map.get(n), zh_map.get(n); return (en, zh) if en and zh else None
def _alkane_names(n: int) -> tuple[str, str] | None: return _pair(ALKANE_EN, ALKANE_ZH, n)
def _cyclo_from_alkane(n: int, en_fn, zh_fn) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain or n < 3:
        return None
    return en_fn(plain[0]), zh_fn(plain[1])
def _cycloalkane_names(n: int) -> tuple[str, str] | None:
    return _cyclo_from_alkane(n, lambda e: f"cyclo{e}", lambda z: f"环{z}")
def _cycloalkene_names(n: int) -> tuple[str, str] | None:
    return _cyclo_from_alkane(n, lambda e: f"cyclo{e[:-3]}ene", lambda z: f"环{zh_stem(z)}烯")
def _cyclo_fg_with_loc(
    n: int, loc: int | None, omit: bool, en_suf: str, zh_suf: str,
) -> tuple[str, str] | None:
    """Unsub: cyclohexanol; sub: cyclohexan-1-ol (keep FG@1)."""
    plain = _alkane_names(n)
    if not plain or n < 3:
        return None
    en, zh = plain
    stem_en, stem_zh = f"cyclo{en[:-1]}", f"环{zh_stem(zh)}"
    if omit or loc is None:
        return f"{stem_en}{en_suf}", f"{stem_zh}{zh_suf}"
    return f"{stem_en}-{loc}-{en_suf}", f"{stem_zh}-{loc}-{zh_suf}"
def _cyclo_ene_fg_names(
    n: int, ene_loc: int | None, fg_loc: int | None, en_suf: str, zh_suf: str,
) -> tuple[str, str] | None:
    """cyclohex-2-en-1-ol / 环己-2-烯-1-醇 (FG@1 + endocyclic ene)."""
    plain = _alkane_names(n)
    if not plain or n < 3 or ene_loc is None or fg_loc is None:
        return None
    en, zh = plain
    return (
        f"cyclo{en[:-3]}-{ene_loc}-en-{fg_loc}-{en_suf}",
        f"环{zh_stem(zh)}-{ene_loc}-烯-{fg_loc}-{zh_suf}",
    )
def _cycloalcohol_names(
    n: int, loc: int | None = None, omit: bool = True,
    ene_loc: int | None = None,
) -> tuple[str, str] | None:
    if ene_loc is not None:
        return _cyclo_ene_fg_names(n, ene_loc, loc, "ol", "醇")
    return _cyclo_fg_with_loc(n, loc, omit, "ol", "醇")
def _cycloketone_names(
    n: int, loc: int | None = None, omit: bool = True,
    ene_loc: int | None = None,
) -> tuple[str, str] | None:
    if ene_loc is not None:
        return _cyclo_ene_fg_names(n, ene_loc, loc, "one", "酮")
    return _cyclo_fg_with_loc(n, loc, omit, "one", "酮")
def _cycloamine_names(
    n: int, loc: int | None = None, omit: bool = True,
) -> tuple[str, str] | None:
    return _cyclo_fg_with_loc(n, loc, omit, "amine", "胺")
def _cycloalkanecarboxylic_names(n: int) -> tuple[str, str] | None:
    return _cyclo_from_alkane(n, lambda e: f"cyclo{e[:-1]}ecarboxylic acid", lambda z: f"环{z}甲酸")
def _omit_term_locant(n: int, loc: int | None, omit: bool) -> bool:
    return omit or loc is None or (loc == 1 and n <= 2)
def _alcohol_names(n: int, oh_locant: int | None, omit: bool) -> tuple[str, str] | None:
    if _omit_term_locant(n, oh_locant, omit):
        return _pair(ALCOHOL_EN, ALCOHOL_ZH, n)
    plain = _pair(ALCOHOL_EN, ALCOHOL_ZH, n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-2]}-{oh_locant}-ol", f"{zh_stem(zh)}-{oh_locant}-醇"
def _thiol_names(n: int, sh_locant: int | None, omit: bool) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    if _omit_term_locant(n, sh_locant, omit):
        return f"{en}thiol", f"{zh_stem(zh)}硫醇"
    return f"{en}-{sh_locant}-thiol", f"{zh_stem(zh)}-{sh_locant}-硫醇"
def _pair_loc_str(locs: list[int]) -> str:
    return ",".join(str(x) for x in locs)
def _poly_fg_names(
    n: int, locs: list[int] | None, en_suf: str, zh_suf: str, need: int
) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain or not locs or len(locs) != need:
        return None
    en, zh = plain
    loc = _pair_loc_str(locs)
    return f"{en}-{loc}-{en_suf}", f"{zh}-{loc}-{zh_suf}"
def _polyamine_names(n: int, locs: list[int] | None, kind: str) -> tuple[str, str] | None:
    m = {"diamine": ("diamine", "二胺", 2), "triamine": ("triamine", "三胺", 3), "tetraamine": ("tetraamine", "四胺", 4)}.get(kind)
    return _poly_fg_names(n, locs, *m) if m else None
def _polyol_names(n: int, locs: list[int] | None, kind: str) -> tuple[str, str] | None:
    if kind == "alcohol":
        return _alcohol_names(n, locs[0] if locs else None, False)
    m = {"diol": ("diol", "二醇", 2), "triol": ("triol", "三醇", 3)}.get(kind)
    return _poly_fg_names(n, locs, *m) if m else None
def _unsat_polyol_ok(numbered: dict) -> tuple | None:
    locs = numbered.get("oh_locants")
    kind = (numbered.get("parent") or {}).get("kind", "")
    if not locs or kind not in ("diol", "triol"): return None
    ene, yne = numbered.get("ene_locant"), numbered.get("yne_locant")
    enes = numbered.get("ene_locants")
    if ene or yne or (enes and len(enes) >= 2): return (locs, kind, ene, yne, enes)
    return None
def _monoene_polyol_name(plain, ene, loc, suf, omit: bool = False) -> tuple[str, str]:
    en, zh = plain
    if omit: return f"{en[:-3]}ene-{loc}-{suf[0]}", f"{zh}烯-{loc}-{suf[1]}"
    return f"{en[:-3]}-{ene}-ene-{loc}-{suf[0]}", f"{zh}-{ene}-烯-{loc}-{suf[1]}"
def _monoyne_polyol_name(plain, yne, loc, suf, omit: bool = False) -> tuple[str, str]:
    en, zh = plain
    if omit: return f"{en[:-3]}yne-{loc}-{suf[0]}", f"{zh}炔-{loc}-{suf[1]}"
    return f"{en[:-3]}-{yne}-yne-{loc}-{suf[0]}", f"{zh}-{yne}-炔-{loc}-{suf[1]}"
def _polyene_polyol_name(plain, enes, loc, suf) -> tuple[str, str]:
    en, zh, n = plain[0], plain[1], len(enes)
    el = _pair_loc_str(enes)
    me, mz = {2: "diene", 3: "triene"}.get(n, ""), {2: "二烯", 3: "三烯"}.get(n, "")
    return f"{en[:-3]}a-{el}-{me}-{loc}-{suf[0]}", f"{zh}-{el}-{mz}-{loc}-{suf[1]}"
def _unsat_polyol_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """but-2-ene-1,4-diol / hexa-2,4-diene-1,6-diol — ene/yne before polyol suffix."""
    u = _unsat_polyol_ok(numbered)
    if u is None: return None
    locs, kind, ene, yne, enes = u
    plain, suf = _alkane_names(n), {"diol": ("diol", "二醇"), "triol": ("triol", "三醇")}.get(kind)
    if not plain or not suf: return None
    loc = _pair_loc_str(locs)
    if enes and len(enes) >= 2: return _polyene_polyol_name(plain, enes, loc, suf)
    if ene is not None:
        return _monoene_polyol_name(plain, ene, loc, suf, numbered.get("omit_ene_locant", False))
    return _monoyne_polyol_name(plain, yne, loc, suf, numbered.get("omit_yne_locant", False))
def _amine_names(n: int, am_locant: int | None, omit: bool) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    if _omit_term_locant(n, am_locant, omit):
        return f"{en[:-1]}amine", f"{zh_stem(zh)}胺"
    return f"{en[:-1]}-{am_locant}-amine", f"{zh_stem(zh)}-{am_locant}-胺"
def _acid_names(n: int) -> tuple[str, str] | None:
    return _pair(ACID_EN, ACID_ZH, n)
def _anhydride_from_acid(n: int) -> tuple[str, str] | None:
    plain = _acid_names(n)
    if not plain:
        return None
    en, zh = plain
    return en.replace(" acid", " anhydride"), f"{zh}酐"
def _diacid_from_alkane(n: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain or n < 3:
        return None
    en, zh = plain
    return f"{en}dioic acid", f"{zh_stem(zh)}二酸"
def _diacid_names(n: int) -> tuple[str, str] | None:
    if n == 2:
        return "oxalic acid", "草酸"
    return _diacid_from_alkane(n)
from namepredict.layer5.polycarboxylic import (
    benzene_polycarboxylic_names as _benzene_polycarboxylic_names,
    cycloalkane_polycarboxylic_names as _cycloalkane_polycarboxylic_names,
    polycarboxylic_names as _polycarboxylic_names,
)
def _aldehyde_names(n: int) -> tuple[str, str] | None:
    return _pair(ALDEHYDE_EN, ALDEHYDE_ZH, n)
def _amide_names(n: int) -> tuple[str, str] | None:
    return _pair(AMIDE_EN, AMIDE_ZH, n)
def _nitrile_names(n: int) -> tuple[str, str] | None:
    return _pair(NITRILE_EN, NITRILE_ZH, n)
def _ester_acyl_en(n: int) -> str | None:
    return ESTER_ACYL_EN.get(n)
def _ester_alkyl_pair(alkoxy_n: int) -> tuple[str, str] | None:
    en, zh = ESTER_ALKYL_EN.get(alkoxy_n), ESTER_ALKYL_ZH.get(alkoxy_n)
    return (en, zh) if en and zh else None
def _ester_names(acyl_n: int, parent: dict) -> tuple[str, str] | None:
    from namepredict.layer5.stems import ester_alkoxy_pair
    alkyl = ester_alkoxy_pair(parent)
    acyl = _ester_acyl_en(acyl_n)
    acid_zh = ACID_ZH.get(acyl_n)
    if not alkyl or not acyl or not acid_zh:
        return None
    return f"{alkyl[0]} {acyl}", f"{acid_zh}{alkyl[1]}酯"
def _ketone_plain(n: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-1]}one", f"{zh_stem(zh)}酮"
def _ketone_from_alkane(n: int, locant: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-1]}-{locant}-one", f"{zh_stem(zh)}-{locant}-酮"
def _ketone_names(n: int, locant: int | None) -> tuple[str, str] | None:
    if locant is None:
        return None
    if n <= 2 and locant == 1:
        return _ketone_plain(n)
    return _ketone_from_alkane(n, locant)
def _dione_names(n: int, locs: list[int] | None) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain or not locs or len(locs) != 2:
        return None
    en, zh = plain
    loc = _pair_loc_str(locs)
    return f"{en}-{loc}-dione", f"{zh_stem(zh)}-{loc}-二酮"
def _alkene_plain(n: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-3]}ene", f"{zh_stem(zh)}烯"
def _alkene_with_locant(n: int, locant: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-3]}-{locant}-ene", f"{zh_stem(zh)}-{locant}-烯"
def _alkene_names(n: int, locant: int | None, omit: bool) -> tuple[str, str] | None:
    if omit or locant is None or n <= 3:
        return _alkene_plain(n)
    return _alkene_with_locant(n, locant)
def _ene_mult(k: int) -> tuple[str, str]:
    en = {2: "diene", 3: "triene", 4: "tetraene", 5: "pentaene"}.get(k, "")
    zh = {2: "二烯", 3: "三烯", 4: "四烯", 5: "五烯"}.get(k, "")
    return en, zh
def _polyene_stem(n: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-3]}a", zh_stem(zh)
def _polyene_names(n: int, locs: list[int] | None) -> tuple[str, str] | None:
    stem = _polyene_stem(n)
    if not stem or not locs or len(locs) < 2:
        return None
    me, mz = _ene_mult(len(locs))
    if not me or not mz:
        return None
    loc = _pair_loc_str(locs)
    return f"{stem[0]}-{loc}-{me}", f"{stem[1]}-{loc}-{mz}"
def _cyclopolyene_names(n: int, locs: list[int] | None) -> tuple[str, str] | None:
    pair = _polyene_names(n, locs)
    return (f"cyclo{pair[0]}", f"环{pair[1]}") if pair else None
def _alkyne_retained(n: int) -> tuple[str, str] | None:
    if n == 2:
        return "acetylene", "乙炔"
    if n == 3:
        return "propyne", "丙炔"
    return None
def _alkyne_with_locant(n: int, locant: int) -> tuple[str, str] | None:
    plain = _alkane_names(n)
    if not plain:
        return None
    en, zh = plain
    return f"{en[:-3]}-{locant}-yne", f"{zh_stem(zh)}-{locant}-炔"
def _alkyne_names(n: int, locant: int | None, omit: bool) -> tuple[str, str] | None:
    retained = _alkyne_retained(n)
    if retained is not None:
        return retained
    if locant is None:
        return None
    return _alkyne_with_locant(n, locant)
def _acid_table(n: int) -> dict:
    return {
        "acid": lambda: _acid_names(n), "diacid": lambda: _diacid_names(n),
        "aldehyde": lambda: _aldehyde_names(n), "amide": lambda: _amide_names(n),
    }
def _unsat_aldehyde(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.unsat_acid import _has_yne, alkenal_names, alkynal_names
    if _has_yne(numbered):
        return alkynal_names(n, numbered)
    locs = numbered.get("ene_locants") or []
    if numbered.get("ene_locant") is None and len(locs) < 2:
        return None
    return alkenal_names(n, numbered)
def _polyacid_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "polycarboxylic": return _polycarboxylic_names(n, numbered)
    if kind == "benzene_polycarboxylic": return _benzene_polycarboxylic_names(n, numbered)
    return _cycloalkane_polycarboxylic_names(n, numbered)

def _acid_ald_amide(kind: str, n: int, numbered: dict | None = None) -> tuple[str, str] | None:
    if kind in {"polycarboxylic", "benzene_polycarboxylic", "cycloalkane_polycarboxylic"}: return _polyacid_names(kind, n, numbered or {})
    from namepredict.layer5.unsat_acid import unsat_carbonyl_names
    top = unsat_carbonyl_names(kind, n, numbered or {})
    if top is not None: return top
    if kind == "aldehyde" and numbered:
        top = _unsat_aldehyde(n, numbered)
        if top is not None: return top
    if kind == "cycloalkanecarboxylic": return _cycloalkanecarboxylic_names(n)
    fn = _acid_table(n).get(kind)
    return fn() if fn else None
def _unsat_nitrile(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.unsat_acid import (
        _has_yne, alkenenitrile_names, alkynenitrile_names,
    )
    if _has_yne(numbered):
        return alkynenitrile_names(n, numbered)
    if numbered.get("ene_locant") is not None:
        return alkenenitrile_names(n, numbered)
    return None
def _nitrile_or_none(kind: str, n: int, numbered: dict | None = None) -> tuple[str, str] | None:
    if kind == "nitrile":
        return _unsat_nitrile(n, numbered or {}) or _nitrile_names(n)
    if kind == "anhydride":
        return _anhydride_from_acid(n)
    return None
def _ester_or_alkenoate(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind != "ester":
        return None
    from namepredict.layer5.unsat_acid import ester_or_alkenoate
    top = ester_or_alkenoate(n, numbered)
    if top is not None:
        return top
    return _ester_names(n, numbered.get("parent") or {})
def _ketone_or_alkenone(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.unsat_acid import (
        _has_ene, _has_yne, alkenone_from, alkynone_from,
    )
    if _has_yne(numbered):
        top = alkynone_from(n, numbered)
        if top is not None:
            return top
    if _has_ene(numbered):
        top = alkenone_from(n, numbered)
        if top is not None:
            return top
    return _ketone_names(n, numbered.get("ketone_locant"))
def _cycloketone_from(n: int, numbered: dict) -> tuple[str, str] | None:
    return _cycloketone_names(
        n, numbered.get("ketone_locant"), numbered.get("omit_ketone_locant", True),
        ene_loc=_ene_loc_kept(numbered),
    )
def _cyclo_poly_fg_names(
    n: int, locs: list[int] | None, en_suf: str, zh_suf: str, need: int,
) -> tuple[str, str] | None:
    """cyclohexane-1,2-diol / 环己烷-1,2-二醇 style poly FG on cycloalkane."""
    plain = _alkane_names(n)
    if not plain or not locs or len(locs) != need or n < 3:
        return None
    en, zh = plain
    loc = _pair_loc_str(locs)
    return f"cyclo{en}-{loc}-{en_suf}", f"环{zh}-{loc}-{zh_suf}"
def _cycloalkanediol_names(n: int, locs: list[int] | None) -> tuple[str, str] | None:
    return _cyclo_poly_fg_names(n, locs, "diol", "二醇", 2)
def _cycloalkanedione_names(n: int, locs: list[int] | None) -> tuple[str, str] | None:
    return _cyclo_poly_fg_names(n, locs, "dione", "二酮", 2)
def _ester_ketone(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _ester_or_alkenoate(kind, n, numbered)
    if top is not None: return top
    if kind == "dione": return _dione_names(n, numbered.get("ketone_locants"))
    if kind == "cycloalkanedione":
        return _cycloalkanedione_names(n, numbered.get("ketone_locants"))
    if kind == "ketone": return _ketone_or_alkenone(n, numbered)
    return _cycloketone_from(n, numbered) if kind == "cycloketone" else None
def _carbonyl_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _acid_ald_amide(kind, n, numbered) or _nitrile_or_none(kind, n, numbered)
    return top if top is not None else _ester_ketone(kind, n, numbered)
def _ene_loc_kept(numbered: dict) -> int | None:
    return None if numbered.get("omit_ene_locant") else numbered.get("ene_locant")
def _cyclo_hetero_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "cycloalcohol":
        return _cycloalcohol_names(
            n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", True),
            ene_loc=_ene_loc_kept(numbered),
        )
    if kind != "cycloamine":
        return None
    return _cycloamine_names(
        n, numbered.get("amine_locant"), numbered.get("omit_amine_locant", True),
    )
def _alcohol_or_alkenol(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.unsat_acid import (
        _has_ene, _has_yne, alkenol_from, alkynol_from,
    )
    if _has_yne(numbered):
        top = alkynol_from(n, numbered)
        if top is not None:
            return top
    if _has_ene(numbered):
        top = alkenol_from(n, numbered)
        if top is not None:
            return top
    return _alcohol_names(n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", False))
def _oh_kind_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "alcohol": return _alcohol_or_alkenol(n, numbered)
    if kind == "thiol": return _thiol_names(n, numbered.get("sh_locant"), numbered.get("omit_sh_locant", False))
    if kind == "benzenediol": return benzenediol_names(numbered.get("oh_locants"))
    if kind == "cycloalkanediol": return _cycloalkanediol_names(n, numbered.get("oh_locants"))
    if kind in ("diol", "triol"):
        return _unsat_polyol_names(n, numbered) or _polyol_names(n, numbered.get("oh_locants"), kind)
    return None
def _amine_kind_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind in ("diamine", "triamine", "tetraamine"):
        return _polyamine_names(n, numbered.get("amine_locants"), kind)
    if kind == "benzenediamine":
        return benzenediamine_names(numbered.get("amine_locants"))
    if kind not in ("amine", "sec_amine", "tert_amine"):
        return None
    return _amine_names(
        n, numbered.get("amine_locant"), numbered.get("omit_amine_locant", False)
    )
def _sym_ether_names(n: int) -> tuple[str, str] | None:
    return _pair(ETHER_SYM_EN, ETHER_SYM_ZH, n)
def _asym_ether_names(parent_n: int, alkoxy_n: int) -> tuple[str, str] | None:
    alk = _alkane_names(parent_n)
    en_pre, zh_pre = ALKOXY_EN.get(alkoxy_n), ALKOXY_ZH.get(alkoxy_n)
    if not alk or not en_pre or not zh_pre:
        return None
    en_p, zh_p = alk
    if parent_n <= 2:
        return f"{en_pre}{en_p}", f"{zh_pre}{zh_p}"
    return f"1-{en_pre}{en_p}", f"1-{zh_pre}{zh_p}"
# Functional-class ether arms: (tag, n) → (en radical, zh radical without 基)
_ETHER_ARM_EN = {
    ("n", 1): "methyl", ("n", 2): "ethyl", ("n", 3): "propyl", ("n", 4): "butyl",
    ("ipr", 0): "isopropyl", ("hfip", 0): "hexafluoroisopropyl",
}
_ETHER_ARM_ZH = {
    ("n", 1): "甲", ("n", 2): "乙", ("n", 3): "丙", ("n", 4): "丁",
    ("ipr", 0): "异丙", ("hfip", 0): "六氟异丙",
}
def _ether_arm_pair(tag: str, n: int) -> tuple[str, str] | None:
    key = (tag, int(n))
    en, zh = _ETHER_ARM_EN.get(key), _ETHER_ARM_ZH.get(key)
    return (en, zh) if en and zh else None
def _func_ether_names(arms) -> tuple[str, str] | None:
    """hexafluoroisopropyl methyl ether / 六氟异丙基甲醚 (alpha EN)."""
    if not arms or len(arms) != 2:
        return None
    pairs = [_ether_arm_pair(t, n) for t, n in arms]
    if any(p is None for p in pairs):
        return None
    a, b = sorted(pairs, key=lambda x: x[0])
    return f"{a[0]} {b[0]} ether", f"{a[1]}基{b[1]}醚"
def _ether_names(n: int, numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    arms = parent.get("ether_arms")
    if arms is not None:
        return _func_ether_names(arms)
    alkoxy_n = parent.get("alkoxy_n")
    if alkoxy_n is None:
        return None
    if alkoxy_n == n:
        return _sym_ether_names(n)
    return _asym_ether_names(n, alkoxy_n)
def _sym_sulfide_names(n: int) -> tuple[str, str] | None:
    return _pair(SULFIDE_SYM_EN, SULFIDE_SYM_ZH, n)
def _asym_sulfide_names(n1: int, n2: int) -> tuple[str, str] | None:
    en1, en2 = SULFIDE_ALKYL_EN.get(n1), SULFIDE_ALKYL_EN.get(n2)
    zh1, zh2 = SULFIDE_ALKYL_ZH.get(n1), SULFIDE_ALKYL_ZH.get(n2)
    if not en1 or not en2 or not zh1 or not zh2:
        return None
    a, b = sorted([(en1, zh1), (en2, zh2)], key=lambda x: x[0])
    return f"{a[0]} {b[0]} sulfide", f"{a[1]}{b[1]}硫醚"
def _sulfide_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    ns = parent.get("alkyl_ns")
    if not ns or len(ns) != 2:
        return None
    n1, n2 = int(ns[0]), int(ns[1])
    if n1 == n2:
        return _sym_sulfide_names(n1)
    return _asym_sulfide_names(n1, n2)
def _ether_or_sulfide(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "ether":
        return _ether_names(n, numbered)
    return _sulfide_names(numbered) if kind == "sulfide" else None
def _hetero_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _ether_or_sulfide(kind, n, numbered)
    if top is not None:
        return top
    top = _oh_kind_names(kind, n, numbered)
    if top is not None:
        return top
    cyc = _cyclo_hetero_names(kind, n, numbered)
    if cyc is not None:
        return cyc
    return _amine_kind_names(kind, n, numbered)
def _typed_acid_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "acid":
        return kind
    if facts.relation.value == "exocyclic" and kind == "cycloalkane":
        return "cycloalkanecarboxylic" if facts.multiplicity == 1 else "cycloalkane_polycarboxylic"
    if facts.relation.value == "exocyclic":
        return kind
    return "acid" if facts.multiplicity == 1 else "diacid" if facts.multiplicity == 2 else "polycarboxylic"


def _typed_ketone_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "ketone":
        return kind
    cyclic = kind in {"cycloketone", "cycloalkanedione"}
    if cyclic:
        return "cycloketone" if facts.multiplicity == 1 else "cycloalkanedione"
    return "ketone" if facts.multiplicity == 1 else "dione"


def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    kind = _typed_ketone_kind(_typed_acid_kind(kind, numbered), numbered)
    from namepredict.layer5.phosphate_names import p_fg_names
    from namepredict.layer5.special_fg_names import special_fg_names
    top = special_fg_names(kind, n, numbered) or p_fg_names(kind, n, numbered)
    if top is not None: return top
    top = _hetero_names(kind, n, numbered) or _carbonyl_names(kind, n, numbered)
    return top if top is not None else _unsat_or_alkane(kind, n, numbered)
def _with_ez(pair: tuple[str, str] | None, numbered: dict) -> tuple[str, str] | None:
    if pair is None: return None
    from namepredict.layer5.stereo_ez import ez_for_parent
    ez = ez_for_parent(numbered)
    return f"{ez}{pair[0]}", f"{ez}{pair[1]}"
def _alkene_or_poly(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "alkene":
        omit = numbered.get("omit_ene_locant", False)
        return _with_ez(_alkene_names(n, numbered.get("ene_locant"), omit), numbered)
    if kind == "polyene":
        return _with_ez(_polyene_names(n, numbered.get("ene_locants")), numbered)
    if kind == "cyclopolyene":
        return _cyclopolyene_names(n, numbered.get("ene_locants"))
    return None
def _unsat_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _alkene_or_poly(kind, n, numbered)
    if top is not None:
        return top
    if kind == "cycloalkene":
        return _cycloalkene_names(n)
    if kind == "alkyne":
        omit = numbered.get("omit_yne_locant", False)
        return _alkyne_names(n, numbered.get("yne_locant"), omit)
    return None
_H5COOH_KINDS = frozenset({
    "furancarboxylic", "thiophenecarboxylic", "pyrrolecarboxylic",
    "imidazolecarboxylic", "pyrazolecarboxylic",
})
_SHCOOH_KINDS = frozenset({
    "piperidinecarboxylic", "pyrrolidinecarboxylic", "piperazinecarboxylic",
    "morpholinecarboxylic", "oxolanecarboxylic", "oxanecarboxylic",
    "thiolanecarboxylic", "aziridinecarboxylic",
})
def _parent_stem_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None
def _spiro_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """spiro[4.5]decane / 螺[4.5]癸烷."""
    stem = _parent_stem_names(numbered)
    alk = _alkane_names(n)
    if stem and alk:
        return f"{stem[0]}{alk[0]}", f"{stem[1]}{alk[1]}"
    return None
def _bridged_names(n: int, numbered: dict) -> tuple[str, str] | None:
    """bicyclo[2.2.1]heptane / 双环[2.2.1]庚烷."""
    stem = _parent_stem_names(numbered)
    alk = _alkane_names(n)
    if stem and alk:
        return f"{stem[0]}{alk[0]}", f"{stem[1]}{alk[1]}"
    return None
def _ring_or_alkane(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "cycloalkane": return _cycloalkane_names(n)
    if kind == "benzene": return benzene_parent_names(numbered)
    if kind == "bridged": return _bridged_names(n, numbered)
    if kind == "spiro": return _spiro_names(n, numbered)
    if kind == "benzoate": return benzoate_parent_names(numbered, _build_prefix)
    if kind in _H5COOH_KINDS: return hetero5carboxylic_names(numbered)
    if kind in _SHCOOH_KINDS: return sat_hetero_carboxylic_names(numbered)
    top = pyridine_kind_names(kind, numbered, _build_prefix)
    if top is not None: return top
    stem = _parent_stem_names(numbered)
    return stem if stem is not None else _alkane_names(n)
def _unsat_or_alkane(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    unsat = _unsat_names(kind, n, numbered)
    return unsat if unsat is not None else _ring_or_alkane(kind, n, numbered)
def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)
def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})
from namepredict.layer5.assembler_prefixes import _build_prefix, _prefix_for

def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    from namepredict.layer5.stereo_rs import apply_rs_prefix
    kind, n = _parent_n(numbered)
    effective_kind = _typed_ketone_kind(_typed_acid_kind(kind, numbered), numbered)
    names = _names_for(effective_kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    en, zh = join_kind_name(effective_kind, _prefix_for(numbered, effective_kind, n), names)
    en, zh = maybe_anion_names(numbered, en, zh)
    en, zh = apply_rs_prefix(numbered, en, zh)
    en, zh = maybe_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
