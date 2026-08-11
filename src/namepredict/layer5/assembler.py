from __future__ import annotations
from dataclasses import dataclass
from namepredict.layer5.stems import (
    ACID_EN, ACID_ZH, ALCOHOL_EN, ALCOHOL_ZH,
    ALDEHYDE_EN, ALDEHYDE_ZH, ALKANE_EN, ALKANE_ZH, ALKOXY_EN, ALKOXY_ZH,
    AMIDE_EN, AMIDE_ZH, ESTER_ACYL_EN,
    ETHER_SYM_EN, ETHER_SYM_ZH, NITRILE_EN, NITRILE_ZH, SULFIDE_ALKYL_EN,
    SULFIDE_ALKYL_ZH, SULFIDE_SYM_EN, SULFIDE_SYM_ZH, _en_stem,
    alkane_zh, maybe_anion_names, maybe_metal_salt_names, zh_stem,
)
from namepredict.layer5.benzene_names import (
    benzene_parent_names, phenyl_parent_names,
    benzoate_parent_names, benzenediol_names,
    join_kind_name, pyridine_kind_names,
)
from namepredict.layer5.stereo_ez import _ez_prefix, ez_for_parent
from namepredict.layer5.unsat_acid import alkenamide_names
from namepredict.types import NameResult
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
_Q_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_Q_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
_Q_MULT_EN = {2: "di", 3: "tri", 4: "tetra"}
_Q_MULT_ZH = {2: "二", 3: "三", 4: "四"}


def _quaternary_radicals(lengths: list[int]) -> tuple[str, str]:
    counts = {n: lengths.count(n) for n in set(lengths)}
    en = "".join(f"{_Q_MULT_EN.get(c, '')}{_Q_EN[n]}" for n, c in sorted(counts.items(), key=lambda x: _Q_EN[x[0]]))
    zh = "".join(f"{_Q_MULT_ZH.get(c, '')}{_Q_ZH[n]}" for n, c in sorted(counts.items(), reverse=True))
    return en, zh


def _tetraalkylammonium_names(numbered: dict) -> tuple[str, str] | None:
    lengths = (numbered.get("parent") or {}).get("quaternary_arm_lengths")
    if not lengths or any(n not in _Q_EN for n in lengths):
        return None
    en, zh = _quaternary_radicals(lengths)
    return f"{en}ammonium", f"{zh}铵"

def _amine_kind_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "tetraalkylammonium":
        return _tetraalkylammonium_names(numbered)
    if kind in ("diamine", "triamine", "tetraamine"):
        return _polyamine_names(n, numbered.get("amine_locants"), kind)
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
    if kind in {"cycloalkane", "cycloketone", "cycloalkanedione"}:
        return "cycloketone" if facts.multiplicity == 1 else "cycloalkanedione"
    if kind not in {"ketone", "dione"}:
        return kind
    return "ketone" if facts.multiplicity == 1 else "dione"


def _typed_ring_alcohol_kind(kind: str, numbered: dict, facts) -> str | None:
    if kind in {"cycloalkane", "cycloalcohol", "cycloalkanediol"}:
        return "cycloalcohol" if facts.multiplicity == 1 else "cycloalkanediol"
    if kind in {"benzene", "phenol", "benzenediol"}:
        return "phenol" if facts.multiplicity == 1 else "benzenediol"
    parent = numbered.get("parent") or {}
    scaffold = parent.get("scaffold_identity")
    if parent.get("typed_ring_expression_supported") and scaffold and scaffold.id == "naphthalene":
        return "naphthalenol" if facts.multiplicity == 1 else "naphthalenediol"
    return None


def _typed_alcohol_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "alcohol":
        return kind
    ring_kind = _typed_ring_alcohol_kind(kind, numbered, facts)
    if ring_kind:
        return ring_kind
    if kind not in {"alcohol", "diol", "triol"}:
        return kind
    return "alcohol" if facts.multiplicity == 1 else "diol" if facts.multiplicity == 2 else "triol"


def _typed_amine_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "amine":
        return kind
    if kind in {"cycloalkane", "cycloamine"} and facts.multiplicity == 1:
        return "cycloamine"
    if kind not in {"amine", "diamine", "triamine", "tetraamine"}:
        return kind
    return {1: "amine", 2: "diamine", 3: "triamine", 4: "tetraamine"}.get(facts.multiplicity, kind)


def _typed_expression_kind(kind: str, numbered: dict) -> str:
    kind = _typed_acid_kind(kind, numbered)
    kind = _typed_ketone_kind(kind, numbered)
    kind = _typed_alcohol_kind(kind, numbered)
    return _typed_amine_kind(kind, numbered)
def _with_ez(pair: tuple[str, str] | None, numbered: dict) -> tuple[str, str] | None:
    if pair is None: return None
    from namepredict.layer5.stereo_ez import ez_for_parent
    ez = ez_for_parent(numbered)
    return f"{ez}{pair[0]}", f"{ez}{pair[1]}"


# ===== 链式词干引擎: 数词干 + coda + 词缀后缀 + 位次 + 环 (替代 if-kind 枚举) =====
def _chain_unsat(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    """通用不饱和段引擎: 炔段优先, 烯段其次 — 段式(醇/酮)与融合式(酸)均由 spec 数据驱动."""
    from namepredict.layer5.unsat_acid import _has_ene, _has_yne
    if _has_yne(numbered) and (spec.yne_seg is not None or spec.yne_suf is not None):
        top = _chain_yne(spec, n, numbered)
        if top is not None:
            return top
    if _has_ene(numbered) and (spec.ene_seg is not None or spec.ene_suf_table is not None):
        top = _chain_ene(spec, n, numbered)
        if top is not None:
            return top
    return None


def _chain_yne(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    s, zs = _en_stem(n), _chain_zh_base(n)
    if s is None or zs is None:
        return None
    yne = numbered.get("yne_locant")
    if spec.yne_suf is not None:      # 融合式: 炔后缀 (酸), omit 时无位次
        if n < 2:
            return None
        omit = numbered.get("omit_yne_locant", False)
        if omit or yne is None:
            return f"{s}{spec.yne_suf[0]}", f"{zs}{spec.yne_suf[1]}"
        return f"{s}-{yne}-{spec.yne_suf[0]}", f"{zs}-{yne}-{spec.yne_suf[1]}"
    fg = numbered.get(spec.loc) if spec.loc else None   # 段式: 需 FG 位次 (醇/酮)
    if yne is None or fg is None:
        return None
    return (
        f"{s}-{yne}-{spec.yne_seg[0]}-{fg}-{spec.en_suf}",
        f"{zs}-{yne}-{spec.yne_seg[1]}-{fg}-{spec.zh_suf}",
    )


def _chain_ene(spec: "_Chain", n: int, numbered: dict) -> tuple[str, str] | None:
    s, zs = _en_stem(n), _chain_zh_base(n)
    if s is None or zs is None:
        return None
    enes = numbered.get("ene_locants")
    if enes and len(enes) >= 2:       # 多烯: 词干加 "a", 后缀表带倍数
        if spec.ene_suf_table is not None:   # 融合式 (酸)
            me, mz = spec.ene_suf_table.get(len(enes), ("", ""))
            if not me or not mz or n < spec.ene_n_min:
                return None
            ez = spec.ez_ene_multi(numbered) if spec.ez_ene_multi else ""
            loc = ",".join(str(x) for x in enes)
            return f"{ez}{s}a-{loc}-{me}", f"{ez}{zs}-{loc}-{mz}"
        if spec.mult_seg is not None:        # 段式多烯 (醇)
            fg = numbered.get(spec.loc) if spec.loc else None
            if fg is None:
                return None
            me, mz = spec.mult_seg.get(len(enes), ("", ""))
            if not me or not mz:
                return None
            ez = spec.ez_ene_multi(numbered) if spec.ez_ene_multi else ""
            loc = ",".join(str(x) for x in enes)
            return f"{ez}{s}a-{loc}-{me}-{fg}-{spec.en_suf}", f"{ez}{zs}-{loc}-{mz}-{fg}-{spec.zh_suf}"
        return None
    ene = numbered.get("ene_locant")
    if spec.ene_suf_table is not None:       # 融合式单烯 (酸)
        if spec.ene_special is not None:
            sp = spec.ene_special(n, numbered)
            if sp is not None:
                return sp
        if ene is None or n < spec.ene_single_min:
            return None
        me, mz = spec.ene_suf_table.get(1, ("", ""))
        ez = spec.ez_ene(numbered) if spec.ez_ene else ""
        return f"{ez}{s}-{ene}-{me}", f"{ez}{zs}-{ene}-{mz}"
    if spec.ene_omit_aware:
        ene = _ene_loc_kept(numbered)
    fg = numbered.get(spec.loc) if spec.loc else None   # 段式单烯 (醇/酮)
    if ene is None or fg is None:
        return None
    ez = (spec.ez_ene(numbered) if spec.ez_ene else "") or ""
    return (
        f"{ez}{s}-{ene}-{spec.ene_seg[0]}-{fg}-{spec.en_suf}",
        f"{ez}{zs}-{ene}-{spec.ene_seg[1]}-{fg}-{spec.zh_suf}",
    )


@dataclass(frozen=True)
class _Chain:
    kind: str
    en_suf: str            # 词缀后缀: "ol" / "one" / "ene" / "oic acid"
    zh_suf: str            # "醇" / "酮" / "烯" / "酸"
    coda: str              # 饱和词干后接 "an"; 烯/炔词干接 ""
    loc: str | None        # FG 位次字段 (单 locant)
    omit_key: str
    default_omit: bool
    no_loc: str            # "plain"=无位次输出普通名; "none"=返回 None
    omit_rule: object      # (n, loc, omit) -> bool  True=省略
    plain_maps: tuple | None = None     # 俗名表 (en_dict, zh_dict), 查不到时回落生成
    plain_fn: object = None             # 派生命名 (n)->pair — 酸酐从酸派生
    locs: str | None = None             # 多 FG 位次字段 (list) — diol/dione/环二酮
    need: int | None = None             # 多 FG 位次数量要求
    ene_seg: tuple | None = None        # 单烯段 ("en","烯") — 醇/酮
    yne_seg: tuple | None = None        # 炔段 ("yn","炔") — 醇/酮
    mult_seg: dict | None = None        # 多烯段表 {2:("dien","二烯"),...} — 醇
    ene_suf_table: dict | None = None   # 融合式烯后缀 {1:("enoic acid","烯酸"),2:...} — 酸
    ene_special: object = None          # 融合式单烯俗名钩子 (n, numbered)->pair — amide 丙烯酰胺
    yne_suf: tuple | None = None        # 融合式炔后缀 ("ynoic acid","炔酸") — 酸
    ez_ene: object = None               # (numbered)->str  单烯 E/Z
    ez_ene_multi: object = None         # (numbered)->str  多烯 E/Z
    ene_n_min: int = 4                  # 多烯融合式 n 下限 (polyene 类为 0)
    ene_single_min: int = 2             # 融合式单烯 n 下限 (diacid 为 3)
    ene_omit_aware: bool = False        # 烯段受 omit_ene_locant 影响 (环系 FG)
    cyclic: bool = False                # 恒加环前缀 (cycloalkane/cycloalkene)
    cyclic_unsat: bool = False          # 仅烯/炔段时加环 (cyclopolyene: 无烯回落纯烷烃)
    zh_full: bool = False               # 中文词干保留完整烷烃后缀 "烷" (环烷/回落)
    wrap: object = None                 # (pair, numbered)->pair  整体包裹 (E/Z)


def _chain_zh_base(n: int) -> str | None:
    z = alkane_zh(n)
    return zh_stem(z) if z else None


def _chain_plain(spec: _Chain, s: str, zs: str, n: int) -> tuple[str, str]:
    """普通名: 俗名表 → 派生命名 → 词干拼接."""
    pair = _pair(*spec.plain_maps, n) if spec.plain_maps else None
    if pair is None and spec.plain_fn is not None:
        pair = spec.plain_fn(n)
    return pair if pair is not None else (f"{s}{spec.coda}{spec.en_suf}", f"{zs}{spec.zh_suf}")


def _chain_names(spec: _Chain, n: int, numbered: dict) -> tuple[str, str] | None:
    """单链词干引擎: 数词干 + coda + 词缀后缀 + 位次 + 环; 烯/炔段插入由 spec 数据驱动."""
    alk = _alkane_names(n)
    if not alk:
        return None
    top = _chain_unsat(spec, n, numbered)
    if top is not None:
        if spec.cyclic or spec.cyclic_unsat:
            top = (f"cyclo{top[0]}", f"环{top[1]}")
        return spec.wrap(top, numbered) if spec.wrap is not None else top
    s = _en_stem(n)
    if spec.zh_full:
        zf = alkane_zh(n)
        zs = zf if zf else _chain_zh_base(n)
    else:
        zs = _chain_zh_base(n)
    if s is None or zs is None:
        return None
    if spec.locs is not None:      # 多 FG 位次 (diol/dione/环二酮...): 无位次 → None
        locs = numbered.get(spec.locs)
        if not locs or (spec.need is not None and len(locs) != spec.need):
            return None
        loc_s = ",".join(str(x) for x in locs)
        pair = (f"{s}{spec.coda}-{loc_s}-{spec.en_suf}", f"{zs}-{loc_s}-{spec.zh_suf}")
    else:
        loc = numbered.get(spec.loc) if spec.loc else None
        omit = numbered.get(spec.omit_key, spec.default_omit) if spec.loc else False
        if loc is None:
            if spec.no_loc == "none":
                return None
            pair = _chain_plain(spec, s, zs, n)
        elif spec.omit_rule(n, loc, omit):
            pair = _chain_plain(spec, s, zs, n)
        else:
            pair = (f"{s}{spec.coda}-{loc}-{spec.en_suf}", f"{zs}-{loc}-{spec.zh_suf}")
    if spec.cyclic:
        pair = (f"cyclo{pair[0]}", f"环{pair[1]}")
    return spec.wrap(pair, numbered) if spec.wrap is not None else pair


_AMINE_SPEC = _Chain(kind="amine", en_suf="amine", zh_suf="胺", coda="an",
                     loc="amine_locant", omit_key="omit_amine_locant", default_omit=False,
                     no_loc="plain", omit_rule=_omit_term_locant)

_KIND_TABLE = {
    "alcohol": _Chain(kind="alcohol", en_suf="ol", zh_suf="醇", coda="an",
                      loc="oh_locant", omit_key="omit_oh_locant", default_omit=False,
                      no_loc="plain", omit_rule=_omit_term_locant,
                      plain_maps=(ALCOHOL_EN, ALCOHOL_ZH),
                      ene_seg=("en", "烯"), yne_seg=("yn", "炔"),
                      mult_seg={2: ("dien", "二烯"), 3: ("trien", "三烯"),
                                4: ("tetraen", "四烯"), 5: ("pentaen", "五烯")},
                      ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent),
    "ketone": _Chain(kind="ketone", en_suf="one", zh_suf="酮", coda="an",
                     loc="ketone_locant", omit_key="omit_ketone_locant", default_omit=False,
                     no_loc="none",
                     omit_rule=lambda n, loc, omit: n <= 2 and loc == 1,
                     ene_seg=("en", "烯"), yne_seg=("yn", "炔"),
                     ez_ene=ez_for_parent),
    "alkene": _Chain(kind="alkene", en_suf="ene", zh_suf="烯", coda="",
                     loc="ene_locant", omit_key="omit_ene_locant", default_omit=False,
                     no_loc="plain", omit_rule=lambda n, loc, omit: omit or n <= 3,
                     wrap=_with_ez),
    "acid": _Chain(kind="acid", en_suf="oic acid", zh_suf="酸", coda="an",
                   loc=None, omit_key="", default_omit=False,
                   no_loc="plain", omit_rule=lambda n, loc, omit: False,
                   plain_maps=(ACID_EN, ACID_ZH),
                   ene_suf_table={1: ("enoic acid", "烯酸"),
                                  2: ("dienoic acid", "二烯酸"),
                                  3: ("trienoic acid", "三烯酸"),
                                  4: ("tetraenoic acid", "四烯酸"),
                                  5: ("pentaenoic acid", "五烯酸"),
                                  6: ("hexaenoic acid", "六烯酸")},
                   yne_suf=("ynoic acid", "炔酸"),
                   ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent),
    "thiol": _Chain(kind="thiol", en_suf="thiol", zh_suf="硫醇", coda="ane",
                    loc="sh_locant", omit_key="omit_sh_locant", default_omit=False,
                    no_loc="plain", omit_rule=_omit_term_locant),
    "amine": _AMINE_SPEC,
    "sec_amine": _AMINE_SPEC,
    "tert_amine": _AMINE_SPEC,
    "aldehyde": _Chain(kind="aldehyde", en_suf="anal", zh_suf="醛", coda="an",
                       loc=None, omit_key="", default_omit=False,
                       no_loc="plain", omit_rule=lambda n, loc, omit: False,
                       plain_maps=(ALDEHYDE_EN, ALDEHYDE_ZH),
                       ene_suf_table={1: ("enal", "烯醛"), 2: ("dienal", "二烯醛"),
                                      3: ("trienal", "三烯醛"), 4: ("tetraenal", "四烯醛"),
                                      5: ("pentaenal", "五烯醛"), 6: ("hexaenal", "六烯醛")},
                       yne_suf=("ynal", "炔醛"),
                       ez_ene=ez_for_parent, ez_ene_multi=ez_for_parent),
    "nitrile": _Chain(kind="nitrile", en_suf="anenitrile", zh_suf="腈", coda="an",
                      loc=None, omit_key="", default_omit=False,
                      no_loc="plain", omit_rule=lambda n, loc, omit: False,
                      plain_maps=(NITRILE_EN, NITRILE_ZH),
                      ene_suf_table={1: ("enenitrile", "烯腈")},
                      yne_suf=("ynenitrile", "炔腈"),
                      ez_ene=ez_for_parent),
    "polyene": _Chain(kind="polyene", en_suf="diene", zh_suf="二烯", coda="",
                      loc=None, omit_key="", default_omit=False,
                      no_loc="none", omit_rule=lambda n, loc, omit: False,
                      ene_suf_table={2: ("diene", "二烯"), 3: ("triene", "三烯"),
                                     4: ("tetraene", "四烯"), 5: ("pentaene", "五烯")},
                      ene_n_min=0, wrap=_with_ez),
    "cyclopolyene": _Chain(kind="cyclopolyene", en_suf="", zh_suf="", coda="ane",
                           loc=None, omit_key="", default_omit=False,
                           no_loc="plain", omit_rule=lambda n, loc, omit: False,
                           ene_suf_table={2: ("diene", "二烯"), 3: ("triene", "三烯"),
                                          4: ("tetraene", "四烯"), 5: ("pentaene", "五烯")},
                           ene_n_min=0, cyclic_unsat=True, zh_full=True),
    "alkyne": _Chain(kind="alkyne", en_suf="yne", zh_suf="炔", coda="",
                     loc="yne_locant", omit_key="omit_yne_locant", default_omit=False,
                     no_loc="plain", omit_rule=lambda n, loc, omit: n in (2, 3),
                     plain_maps=({2: "acetylene", 3: "propyne"}, {2: "乙炔", 3: "丙炔"})),
    "cycloalkane": _Chain(kind="cycloalkane", en_suf="", zh_suf="", coda="ane",
                          loc=None, omit_key="", default_omit=False,
                          no_loc="plain", omit_rule=lambda n, loc, omit: False,
                          cyclic=True, zh_full=True),
    "cycloalkene": _Chain(kind="cycloalkene", en_suf="ene", zh_suf="烯", coda="",
                          loc=None, omit_key="", default_omit=False,
                          no_loc="plain", omit_rule=lambda n, loc, omit: False,
                          cyclic=True),
    "dione": _Chain(kind="dione", en_suf="dione", zh_suf="二酮", coda="ane",
                    loc=None, omit_key="", default_omit=False,
                    no_loc="none", omit_rule=lambda n, loc, omit: False,
                    locs="ketone_locants", need=2),
    "cycloalkanedione": _Chain(kind="cycloalkanedione", en_suf="dione", zh_suf="二酮", coda="ane",
                               loc=None, omit_key="", default_omit=False,
                               no_loc="none", omit_rule=lambda n, loc, omit: False,
                               locs="ketone_locants", need=2, cyclic=True, zh_full=True),
    "cycloalcohol": _Chain(kind="cycloalcohol", en_suf="ol", zh_suf="醇", coda="an",
                           loc="oh_locant", omit_key="omit_oh_locant", default_omit=True,
                           no_loc="plain", omit_rule=lambda n, loc, omit: omit,
                           ene_seg=("en", "烯"), ene_omit_aware=True, cyclic=True),
    "cycloketone": _Chain(kind="cycloketone", en_suf="one", zh_suf="酮", coda="an",
                          loc="ketone_locant", omit_key="omit_ketone_locant", default_omit=True,
                          no_loc="plain", omit_rule=lambda n, loc, omit: omit,
                          ene_seg=("en", "烯"), ene_omit_aware=True, cyclic=True),
    "cycloamine": _Chain(kind="cycloamine", en_suf="amine", zh_suf="胺", coda="an",
                         loc="amine_locant", omit_key="omit_amine_locant", default_omit=True,
                         no_loc="plain", omit_rule=lambda n, loc, omit: omit,
                         cyclic=True),
    "diacid": _Chain(kind="diacid", en_suf="dioic acid", zh_suf="二酸", coda="ane",
                     loc=None, omit_key="", default_omit=False,
                     no_loc="plain", omit_rule=lambda n, loc, omit: False,
                     plain_maps=({2: "oxalic acid"}, {2: "草酸"}),
                     ene_suf_table={1: ("enedioic acid", "烯二酸"),
                                    2: ("dienedioic acid", "二烯二酸"),
                                    3: ("trienedioic acid", "三烯二酸"),
                                    4: ("tetraenedioic acid", "四烯二酸")},
                     ene_single_min=3,
                     ez_ene=_ez_prefix, ez_ene_multi=ez_for_parent),
    "diamine": _Chain(kind="diamine", en_suf="diamine", zh_suf="二胺", coda="ane",
                      loc=None, omit_key="", default_omit=False,
                      no_loc="none", omit_rule=lambda n, loc, omit: False,
                      locs="amine_locants", need=2, zh_full=True),
    "triamine": _Chain(kind="triamine", en_suf="triamine", zh_suf="三胺", coda="ane",
                       loc=None, omit_key="", default_omit=False,
                       no_loc="none", omit_rule=lambda n, loc, omit: False,
                       locs="amine_locants", need=3, zh_full=True),
    "tetraamine": _Chain(kind="tetraamine", en_suf="tetraamine", zh_suf="四胺", coda="ane",
                         loc=None, omit_key="", default_omit=False,
                         no_loc="none", omit_rule=lambda n, loc, omit: False,
                         locs="amine_locants", need=4, zh_full=True),
    "anhydride": _Chain(kind="anhydride", en_suf="", zh_suf="", coda="ane",
                        loc=None, omit_key="", default_omit=False,
                        no_loc="plain", omit_rule=lambda n, loc, omit: False,
                        plain_fn=_anhydride_from_acid),
    "cycloalkanediol": _Chain(kind="cycloalkanediol", en_suf="diol", zh_suf="二醇", coda="ane",
                              loc=None, omit_key="", default_omit=False,
                              no_loc="none", omit_rule=lambda n, loc, omit: False,
                              locs="oh_locants", need=2, cyclic=True, zh_full=True),
    "amide": _Chain(kind="amide", en_suf="amide", zh_suf="酰胺", coda="an",
                    loc=None, omit_key="", default_omit=False,
                    no_loc="plain", omit_rule=lambda n, loc, omit: False,
                    plain_maps=(AMIDE_EN, AMIDE_ZH),
                    ene_suf_table={1: ("enamide", "烯酰胺")},
                    ene_special=alkenamide_names,
                    yne_suf=("ynamide", "炔酰胺"),
                    ez_ene=_ez_prefix),
}


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

def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    """Chain-engine dispatch for table kinds, then specific workers."""
    entry = _KIND_TABLE.get(kind)
    if entry is not None:
        return _chain_names(entry, n, numbered)
    if kind == "ether":
        return _ether_names(n, numbered)
    if kind == "sulfide":
        return _sulfide_names(numbered)
    if kind == "thiol":
        return _thiol_names(n, numbered.get("sh_locant"), numbered.get("omit_sh_locant", False))
    if kind == "benzenediol":
        return benzenediol_names(numbered.get("oh_locants"))
    if kind == "cycloalkanediol":
        return _cycloalkanediol_names(n, numbered.get("oh_locants"))
    if kind in ("diol", "triol"):
        top = _unsat_polyol_names(n, numbered)
        if top is not None:
            return top
        return _polyol_names(n, numbered.get("oh_locants"), kind)
    if kind == "cycloalcohol":
        return _cycloalcohol_names(
            n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", True),
            ene_loc=_ene_loc_kept(numbered),
        )
    if kind == "cycloamine":
        return _cycloamine_names(
            n, numbered.get("amine_locant"), numbered.get("omit_amine_locant", True),
        )
    if kind == "tetraalkylammonium":
        return _tetraalkylammonium_names(numbered)
    if kind in ("diamine", "triamine", "tetraamine"):
        return _polyamine_names(n, numbered.get("amine_locants"), kind)
    if kind in ("amine", "sec_amine", "tert_amine"):
        return _amine_names(
            n, numbered.get("amine_locant"), numbered.get("omit_amine_locant", False),
        )
    if kind in ("polycarboxylic", "benzene_polycarboxylic", "cycloalkane_polycarboxylic"):
        return _polyacid_names(kind, n, numbered)
    from namepredict.layer5.unsat_acid import unsat_carbonyl_names
    if kind in ("diacid", "amide"):
        top = unsat_carbonyl_names(kind, n, numbered)
        if top is not None:
            return top
    if kind == "diacid":
        return _diacid_names(n)
    if kind == "amide":
        return _amide_names(n)
    if kind == "aldehyde":
        top = _unsat_aldehyde(n, numbered)
        if top is not None:
            return top
        return _aldehyde_names(n)
    if kind == "nitrile":
        top = _unsat_nitrile(n, numbered)
        if top is not None:
            return top
        return _nitrile_names(n)
    if kind == "anhydride":
        return _anhydride_from_acid(n)
    if kind == "ester":
        from namepredict.layer5.unsat_acid import ester_or_alkenoate
        top = ester_or_alkenoate(n, numbered)
        if top is not None:
            return top
        return _ester_names(n, numbered.get("parent") or {})
    if kind == "dione":
        return _dione_names(n, numbered.get("ketone_locants"))
    if kind == "cycloalkanedione":
        return _cycloalkanedione_names(n, numbered.get("ketone_locants"))
    if kind == "cycloketone":
        return _cycloketone_from(n, numbered)
    if kind == "polyene":
        return _with_ez(_polyene_names(n, numbered.get("ene_locants")), numbered)
    if kind == "cyclopolyene":
        return _cyclopolyene_names(n, numbered.get("ene_locants"))
    if kind == "cycloalkene":
        return _cycloalkene_names(n)
    if kind == "alkyne":
        return _alkyne_names(n, numbered.get("yne_locant"), numbered.get("omit_yne_locant", False))
    if kind == "cycloalkane":
        return _cycloalkane_names(n)
    if kind == "phenyl":
        return phenyl_parent_names(numbered)
    if kind == "benzene":
        return benzene_parent_names(numbered)
    if kind == "benzoate":
        return benzoate_parent_names(numbered, _build_prefix)
    top = pyridine_kind_names(kind, numbered, _build_prefix)
    if top is not None:
        return top
    stem = _parent_stem_names(numbered)
    return stem if stem is not None else _alkane_names(n)



def _parent_stem_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    en, zh = parent.get("stem_en"), parent.get("stem_zh")
    return (en, zh) if en and zh else None

def _ring_or_alkane(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "cycloalkane": return _cycloalkane_names(n)
    if kind == "phenyl": return phenyl_parent_names(numbered)
    if kind == "benzene": return benzene_parent_names(numbered)
    if kind == "benzoate": return benzoate_parent_names(numbered, _build_prefix)
    top = pyridine_kind_names(kind, numbered, _build_prefix)
    if top is not None: return top
    stem = _parent_stem_names(numbered)
    return stem if stem is not None else _alkane_names(n)

def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)
def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})
from namepredict.layer5.assembler_prefixes import _build_prefix, _prefix_for

def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    from namepredict.layer5.stereo_rs import apply_rs_prefix
    kind, n = _parent_n(numbered)
    effective_kind = _typed_expression_kind(kind, numbered)
    names = _names_for(effective_kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    en, zh = join_kind_name(effective_kind, _prefix_for(numbered, effective_kind, n), names)
    en, zh = maybe_anion_names(numbered, en, zh)
    en, zh = apply_rs_prefix(numbered, en, zh)
    en, zh = maybe_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
