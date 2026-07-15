from __future__ import annotations
from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.layer5.stems import (
    ACID_EN, ACID_ZH, ACYL_CHLORIDE_EN, ACYL_CHLORIDE_ZH, ALCOHOL_EN, ALCOHOL_ZH,
    ALDEHYDE_EN, ALDEHYDE_ZH, ALKANE_EN, ALKANE_ZH, ALKOXY_EN, ALKOXY_ZH,
    AMIDE_EN, AMIDE_ZH, ESTER_ACYL_EN, ESTER_ALKYL_EN, ESTER_ALKYL_ZH,
    ETHER_SYM_EN, ETHER_SYM_ZH, NITRILE_EN, NITRILE_ZH, SULFIDE_ALKYL_EN,
    SULFIDE_ALKYL_ZH, SULFIDE_SYM_EN, SULFIDE_SYM_ZH, maybe_anion_names, maybe_metal_salt_names, zh_stem,
)
from namepredict.layer5.benzene_names import (
    arene_fg_parent_names, benzene_parent_names, benzene_prefix,
    benzoate_parent_names, benzenediamine_names, benzenediol_names, hetero5carboxylic_names,
    join_kind_name, pyridine_kind_names, sat_hetero_carboxylic_names,
)
from namepredict.types import NameResult
MULT_EN = {2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
MULT_ZH = {2: "二", 3: "三", 4: "四", 5: "五", 6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}
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
def _cycloalcohol_names(n: int) -> tuple[str, str] | None:
    return _cyclo_from_alkane(n, lambda e: f"cyclo{e[:-1]}ol", lambda z: f"环{zh_stem(z)}醇")
def _cycloketone_names(n: int) -> tuple[str, str] | None:
    return _cyclo_from_alkane(n, lambda e: f"cyclo{e[:-1]}one", lambda z: f"环{zh_stem(z)}酮")
def _cycloamine_names(n: int) -> tuple[str, str] | None:
    return _cyclo_from_alkane(n, lambda e: f"cyclo{e[:-1]}amine", lambda z: f"环{zh_stem(z)}胺")
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
def _diamine_names(n: int, locs: list[int] | None) -> tuple[str, str] | None:
    return _poly_fg_names(n, locs, "diamine", "二胺", 2)
def _polyol_names(n: int, locs: list[int] | None, kind: str) -> tuple[str, str] | None:
    m = {"diol": ("diol", "二醇", 2), "triol": ("triol", "三醇", 3)}.get(kind)
    return _poly_fg_names(n, locs, *m) if m else None
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
    """ethanone / propanone form (no locant; used when n≤2)."""
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
    if numbered.get("ene_locant") is None:
        return None
    return alkenal_names(n, numbered)
def _acid_ald_amide(kind: str, n: int, numbered: dict | None = None) -> tuple[str, str] | None:
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
    if kind == "acyl_chloride":
        return _pair(ACYL_CHLORIDE_EN, ACYL_CHLORIDE_ZH, n)
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
def _ester_ketone(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _ester_or_alkenoate(kind, n, numbered)
    if top is not None: return top
    if kind == "dione":
        return _dione_names(n, numbered.get("ketone_locants"))
    if kind == "ketone":
        return _ketone_names(n, numbered.get("ketone_locant"))
    return _cycloketone_names(n) if kind == "cycloketone" else None
def _carbonyl_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _acid_ald_amide(kind, n, numbered) or _nitrile_or_none(kind, n, numbered)
    return top if top is not None else _ester_ketone(kind, n, numbered)
def _cyclo_hetero_names(kind: str, n: int) -> tuple[str, str] | None:
    if kind == "cycloalcohol":
        return _cycloalcohol_names(n)
    return _cycloamine_names(n) if kind == "cycloamine" else None
def _alcohol_or_alkenol(n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.unsat_acid import (
        _has_ene, _has_yne, alkenol_from, alkynol_from,
    )
    if _has_yne(numbered):
        return alkynol_from(n, numbered)
    if _has_ene(numbered):
        return alkenol_from(n, numbered)
    return _alcohol_names(n, numbered.get("oh_locant"), numbered.get("omit_oh_locant", False))
def _oh_kind_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "alcohol":
        return _alcohol_or_alkenol(n, numbered)
    if kind == "thiol":
        return _thiol_names(n, numbered.get("sh_locant"), numbered.get("omit_sh_locant", False))
    if kind == "benzenediol":
        return benzenediol_names(numbered.get("oh_locants"))
    return _polyol_names(n, numbered.get("oh_locants"), kind) if kind in ("diol", "triol") else None
def _amine_kind_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "diamine":
        return _diamine_names(n, numbered.get("amine_locants"))
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
def _ether_names(n: int, numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
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
    cyc = _cyclo_hetero_names(kind, n)
    if cyc is not None:
        return cyc
    return _amine_kind_names(kind, n, numbered)
def _special_fg_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.carbamate_names import carbamate_names
    from namepredict.layer5.diester_names import diester_names
    from namepredict.layer5.isocyanate_names import iso_kind_names
    from namepredict.layer5.sulfoxide_names import sulfoxide_names
    if kind == "carbamate": return carbamate_names(numbered)
    if kind == "diester": return diester_names(n, numbered)
    if kind in ("isocyanate", "isothiocyanate"): return iso_kind_names(kind, n)
    return sulfoxide_names(numbered) if kind == "sulfoxide" else None
def _names_for(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    from namepredict.layer5.phosphate_names import p_fg_names
    top = _special_fg_names(kind, n, numbered) or p_fg_names(kind, n, numbered)
    if top is not None: return top
    top = _hetero_names(kind, n, numbered) or _carbonyl_names(kind, n, numbered)
    return top if top is not None else _unsat_or_alkane(kind, n, numbered)
def _with_ez(pair: tuple[str, str] | None, numbered: dict) -> tuple[str, str] | None:
    if pair is None:
        return None
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
def _ring_or_alkane(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind == "cycloalkane": return _cycloalkane_names(n)
    if kind == "benzene": return benzene_parent_names(numbered)
    if kind == "benzoate": return benzoate_parent_names(numbered, _build_prefix)
    if kind in _H5COOH_KINDS: return hetero5carboxylic_names(numbered)
    if kind in _SHCOOH_KINDS: return sat_hetero_carboxylic_names(numbered)
    stem = _parent_stem_names(numbered)
    if stem is not None: return stem
    top = pyridine_kind_names(kind, numbered, _build_prefix)
    return top if top is not None else (arene_fg_parent_names(kind) or _alkane_names(n))
def _unsat_or_alkane(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    unsat = _unsat_names(kind, n, numbered)
    return unsat if unsat is not None else _ring_or_alkane(kind, n, numbered)
def _parent_n(numbered: dict) -> tuple[str | None, int]:
    parent = numbered.get("parent") or {}
    return parent.get("kind"), int(parent.get("n_carbons") or 0)
def _unsupported(n: int, kind: str | None) -> NameResult:
    return _fail({"reason": "unsupported", "n_carbons": n, "kind": kind})
def _group_by_stem(substituents: list) -> dict[str, list]:
    groups: dict[str, list] = {}
    for s in substituents:
        groups.setdefault(s.get("en") or "", []).append(s)
    return groups
def _locant_str(subs: list) -> str:
    locs = sorted(int(s["locant"]) for s in subs if "locant" in s)
    return ",".join(str(x) for x in locs)
def _mult_en(n: int) -> str: return MULT_EN.get(n, "")
def _mult_zh(n: int) -> str: return MULT_ZH.get(n, "")
_KEEP_LOCANT_KINDS = frozenset({
    "acid",
    "benzoic", "benzaldehyde", "acetophenone", "pyridinecarboxylic",
    "pyridinecarbonitrile", "benzoate", "benzonitrile", "benzoyl_chloride",
    "cycloalkanecarboxylic"}) | _H5COOH_KINDS | _SHCOOH_KINDS
def _omit_sub_locants(n_carbons: int, substituents: list, kind: str | None = None) -> bool:
    if n_carbons <= 1 or (kind in ("cycloalkane", "benzene") and len(substituents) == 1):
        return True
    if kind in ("sec_amine", "tert_amine", "amide"):
        return {s.get("kind") for s in substituents} <= {"n_alkyl", "n_phenyl", "n_benzyl"}
    if kind in _KEEP_LOCANT_KINDS or kind == "ketone":
        return False
    if any(s.get("paren") or (s.get("en") or "")[:1] == "(" for s in substituents):
        return False
    return n_carbons == 2 and len(substituents) == 1

def _stem_needs_paren(stem: str, subs: list, omit: bool) -> bool:
    """Paren: explicit flag, leading-locant stems, or multi CF3 (EN)."""
    if any(s.get("paren") for s in subs):
        return True
    if stem and stem[0].isdigit():
        return True
    return (not omit) and stem == "trifluoromethyl"
def _wrap_stem(stem: str, need: bool) -> str: return stem if not need else (f"[{stem}]" if "(" in stem else f"({stem})")
def _prefix_one_en(stem: str, subs: list, omit: bool) -> str:
    mult = _mult_en(len(subs))
    s = _wrap_stem(stem, _stem_needs_paren(stem, subs, omit))
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"
def _prefix_one_zh(zh_stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> str:
    mult = _mult_zh(len(subs))
    en = subs[0].get("en") or ""
    need = any(s.get("paren") for s in subs) or (en[:1].isdigit() if en else False)
    if paren_cf3 and zh_stem == "三氟甲基":
        need = True
    s = _wrap_stem(zh_stem, need)
    return f"{mult}{s}" if omit else f"{_locant_str(subs)}-{mult}{s}"
def _sorted_stems(groups: dict[str, list]) -> list[str]: return sorted((k for k in groups if k), key=alkyl_alpha_key)
def _parts_for_stem(stem: str, subs: list, omit: bool, paren_cf3: bool = False) -> tuple[str, str]:
    zh_stem = subs[0].get("zh") or ""
    if (subs[0].get("kind") or "") in ("n_alkyl", "n_phenyl", "n_benzyl"):
        omit = True
    return _prefix_one_en(stem, subs, omit), _prefix_one_zh(zh_stem, subs, omit, paren_cf3)
def _collect_parts(groups: dict[str, list], omit: bool, paren_cf3: bool = False) -> tuple[list[str], list[str]]:
    en_parts: list[str] = []
    zh_parts: list[str] = []
    for stem in _sorted_stems(groups):
        en_p, zh_p = _parts_for_stem(stem, groups[stem], omit, paren_cf3)
        en_parts.append(en_p)
        zh_parts.append(zh_p)
    return en_parts, zh_parts
def _build_prefix(substituents: list, n_carbons: int, kind: str | None = None) -> tuple[str, str]:
    if not substituents:
        return "", ""
    omit = _omit_sub_locants(n_carbons, substituents, kind)
    paren = kind == "benzene" and len(substituents) >= 4
    en_parts, zh_parts = _collect_parts(_group_by_stem(substituents), omit, paren)
    return "-".join(en_parts), "-".join(zh_parts)
def _prefix_for(numbered: dict, kind: str | None, n: int) -> tuple[str, str]:
    if kind == "benzene":
        return benzene_prefix(numbered, _build_prefix)
    _skip = ("benzoate", "pyrimidinamine", "benzothiazolamine",
             "benzoxazolamine", "benzimidazolamine")
    if kind in _skip:
        return "", ""
    return _build_prefix(numbered.get("substituents") or [], n, kind)
def assemble(numbered: dict, *, time_ms: float = 0.0, source: str = "iupac") -> NameResult:
    from namepredict.layer5.stereo_rs import apply_rs_prefix
    kind, n = _parent_n(numbered)
    names = _names_for(kind, n, numbered)
    if not names:
        return _unsupported(n, kind)
    en, zh = join_kind_name(kind, _prefix_for(numbered, kind, n), names)
    en, zh = maybe_anion_names(numbered, en, zh)
    en, zh = apply_rs_prefix(numbered, en, zh)
    en, zh = maybe_metal_salt_names(numbered, en, zh)
    return _ok(en, zh, time_ms, source)
