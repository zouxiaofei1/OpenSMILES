from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.aryl_sub import (
    _benzyl_name,
    _benzyloxy_name,
    _phenoxy_name,
    _phenyl_name,
    _ring_benzyls,
    _ring_benzyloxys,
    _ring_phenoxys,
    _ring_phenyls,
)
from namepredict.layer2.side_alkyl import (
    _c_neighbors,
    _is_2_methylbutan_2_yl,
    _is_isobutyl,
    _is_isopentyl,
    _is_isopropyl,
    _is_neopentyl,
    _is_sec_butyl,
    _is_tert_butyl,
    _is_trifluoromethyl,
    _terminal_halo_z,
    _walk_linear,
    _walk_omega_halo,
)
from namepredict.layer3.alkoxy_names import _extract_alkoxys
from namepredict.layer3.cycloalkyl_names import _one_cycloalkyl_side
from namepredict.layer3.amino_side import _extract_aminos as _extract_aminos_impl


def _strip_ital_prefix(stem: str) -> str:
    if stem.startswith("tert-") or stem.startswith("sec-"):
        return stem[stem.index("-") + 1 :]
    return stem


def _strip_n_prefix(stem: str) -> str:
    if stem.startswith("N,"):
        return stem.split("-")[-1] if "-" in stem else stem
    return stem[2:] if stem.startswith("N-") else stem


def _strip_lead_locant(stem: str) -> str:
    i = 0
    while i < len(stem) and stem[i].isdigit():
        i += 1
    return stem[i + 1 :] if i and i < len(stem) and stem[i] == "-" else stem


def alkyl_alpha_key(stem: str) -> str:
    """Alphanumerical-order key: sec-/tert-/N-/leading locants ignored (P-14.5)."""
    return _strip_lead_locant(_strip_n_prefix(_strip_ital_prefix(stem)))


ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}
# ω-halo n-alkyl: locant = chain length (terminal carbon)
_HALOALKYL_STEM_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_HALOALKYL_STEM_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _side_starts(mol: Mol, chain: list[int]) -> list[tuple[int, int]]:
    chain_set = set(chain)
    out: list[tuple[int, int]] = []
    for c in chain:
        for nb in _c_neighbors(mol, c):
            if nb not in chain_set:
                out.append((c, nb))
    return out


def _make_alkyl(attach: int, path: list[int]) -> dict:
    n = len(path)
    return {
        "kind": "alkyl",
        "n_carbons": n,
        "attach_idx": attach,
        "atoms": path,
        "en": ALKYL_EN[n],
        "zh": ALKYL_ZH[n],
    }


def _make_branch(
    attach: int, atoms: list[int], n: int, en: str, zh: str,
) -> dict:
    return {
        "kind": "alkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": atoms, "en": en, "zh": zh,
    }


def _make_cf3(attach: int, atoms: list[int]) -> dict:
    return {
        "kind": "trifluoromethyl",
        "n_carbons": 1,
        "attach_idx": attach,
        "atoms": atoms,
        "en": "trifluoromethyl",
        "zh": "三氟甲基",
    }


def _haloalkyl_names(n: int, z: int) -> tuple[str, str]:
    se, sz = _HALOALKYL_STEM_EN[n], _HALOALKYL_STEM_ZH[n]
    if n == 1:
        return f"{HALO_EN[z]}{se}", f"{HALO_ZH[z]}{sz}"
    return f"{n}-{HALO_EN[z]}{se}", f"{n}-{HALO_ZH[z]}{sz}"


def _make_haloalkyl(attach: int, path: list[int], z: int) -> dict:
    n = len(path)
    en, zh = _haloalkyl_names(n, z)
    return {
        "kind": "haloalkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": path, "en": en, "zh": zh, "paren": True,
    }


def _one_haloalkyl(mol: Mol, attach: int, start: int, chain_set: set[int]) -> dict | None:
    path = _walk_omega_halo(mol, start, chain_set)
    if not path:
        return None
    z = _terminal_halo_z(mol, path[-1])
    return _make_haloalkyl(attach, path, z) if z else None


def _one_alkyl(mol: Mol, attach: int, start: int, chain_set: set[int]) -> dict | None:
    path = _walk_linear(mol, start, chain_set)
    if path:
        return _make_alkyl(attach, path)
    ha = _one_haloalkyl(mol, attach, start, chain_set)
    return ha if ha is not None else _one_branched(mol, attach, start, chain_set)


_BRANCH_CHECKS = (
    (_is_isopropyl, 3, "isopropyl", "异丙基"),
    (_is_tert_butyl, 4, "tert-butyl", "叔丁基"),
    (_is_2_methylbutan_2_yl, 5, "2-methylbutan-2-yl", "2-甲基丁-2-基"),
    (_is_isobutyl, 4, "isobutyl", "异丁基"),
    (_is_sec_butyl, 4, "sec-butyl", "仲丁基"),
    (_is_neopentyl, 5, "neopentyl", "新戊基"),
    (_is_isopentyl, 5, "isopentyl", "异戊基"),
)


def _one_branched(mol: Mol, attach: int, start: int, chain_set: set[int]) -> dict | None:
    for fn, n, en, zh in _BRANCH_CHECKS:
        atoms = fn(mol, start, chain_set)
        if atoms:
            return _make_branch(attach, atoms, n, en, zh)
    cyc = _one_cycloalkyl_side(mol, attach, start, chain_set)
    if cyc is not None:
        return cyc
    cf3 = _is_trifluoromethyl(mol, start, chain_set)
    return _make_cf3(attach, cf3) if cf3 else None


def _make_halo(attach: int, halo_idx: int, z: int) -> dict:
    return {
        "kind": "halo",
        "attach_idx": attach,
        "atoms": [halo_idx],
        "en": HALO_EN[z],
        "zh": HALO_ZH[z],
    }


def _halo_on_carbon(mol: Mol, c_idx: int) -> list[dict]:
    atom = mol.GetAtomWithIdx(c_idx)
    out: list[dict] = []
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z in HALO_EN:
            out.append(_make_halo(c_idx, n.GetIdx(), z))
    return out


def _extract_halos(mol: Mol, chain: list[int]) -> list[dict]:
    out: list[dict] = []
    for c in chain:
        out.extend(_halo_on_carbon(mol, c))
    return out


def _extract_alkyls(mol: Mol, chain: list[int]) -> list[dict]:
    chain_set = set(chain)
    out: list[dict] = []
    for attach, start in _side_starts(mol, chain):
        one = _one_alkyl(mol, attach, start, chain_set)
        if one is not None:
            out.append(one)
    return out


def _filter_fg_halos(halos: list, parent: dict) -> list:
    if parent.get("kind") != "acyl_chloride":
        return halos
    cl = parent.get("cl_idx")
    return [h for h in halos if cl not in (h.get("atoms") or [])]


_PARENT_OH_KINDS = frozenset(
    {
        "alcohol", "alkenol", "diol", "triol", "cycloalcohol",
        "phenol", "benzenediol", "pyridinol", "benzothiophenol", "quinolinol",
    }
)
_PARENT_NH2_KINDS = frozenset(
    {
        "amine", "diamine", "cycloamine", "sec_amine", "tert_amine",
        "aniline", "pyridinamine", "pyrimidinamine", "benzofuranamine",
        "benzothiazolamine", "benzoxazolamine", "benzimidazolamine",
        "benzenediamine",
    }
)
_PARENT_OXO_KINDS = frozenset({"ketone", "dione", "cycloketone"})


def _make_hydroxy(attach: int, o_idx: int) -> dict:
    return {
        "kind": "hydroxy",
        "attach_idx": attach,
        "atoms": [o_idx],
        "en": "hydroxy",
        "zh": "羟基",
    }


def _make_oxo(attach: int) -> dict:
    return {
        "kind": "oxo",
        "attach_idx": attach,
        "atoms": [attach],
        "en": "oxo",
        "zh": "氧代",
    }


def _make_nitro(attach: int, n_idx: int, o_idxs: list[int]) -> dict:
    return {
        "kind": "nitro",
        "attach_idx": attach,
        "atoms": [n_idx] + list(o_idxs),
        "en": "nitro",
        "zh": "硝基",
    }


def _extract_nitros(info: dict, parent: dict) -> list[dict]:
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for n in info.get("nitros") or []:
        if n["c_idx"] in chain_set:
            out.append(_make_nitro(n["c_idx"], n["n_idx"], n.get("o_idxs") or []))
    return out


def _extract_hydroxys(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") in _PARENT_OH_KINDS:
        return []
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for h in info.get("hydroxyls") or []:
        if h["c_idx"] in chain_set:
            out.append(_make_hydroxy(h["c_idx"], h["o_idx"]))
    return out


def _extract_aminos(info: dict, parent: dict) -> list[dict]:
    return _extract_aminos_impl(info, parent, _PARENT_NH2_KINDS)




def _extract_oxos(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") in _PARENT_OXO_KINDS:
        return []
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for k in info.get("ketones") or []:
        if k["c_idx"] in chain_set:
            out.append(_make_oxo(k["c_idx"]))
    return out


_N_ALKYL_EN = {1: "N-methyl", 2: "N-ethyl", 3: "N-propyl", 4: "N-butyl"}
_N_ALKYL_ZH = {1: "N-甲基", 2: "N-乙基", 3: "N-丙基", 4: "N-丁基"}
_N_STEM_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_N_STEM_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _n_alkyl_sub(en: str, zh: str, attach: int, n: int) -> dict:
    return {
        "kind": "n_alkyl", "n_carbons": n, "attach_idx": attach,
        "atoms": [], "en": en, "zh": zh,
    }
def _n_phenyl_sub(attach: int) -> dict:
    return {
        "kind": "n_phenyl", "n_carbons": 6, "attach_idx": attach,
        "atoms": [], "en": "N-phenyl", "zh": "N-苯基",
    }


def _n_benzyl_en_zh(en: str, zh: str, paren: bool) -> tuple[str, str]:
    if paren or en != "benzyl":
        return f"N-({en})", f"N-({zh})"
    return f"N-{en}", f"N-{zh}"


def _n_benzyl_sub(attach: int, en: str, zh: str, paren: bool) -> dict:
    ne, nz = _n_benzyl_en_zh(en, zh, paren)
    return {
        "kind": "n_benzyl", "n_carbons": 7, "attach_idx": attach,
        "atoms": [], "en": ne, "zh": nz,
    }


def _tert_n_prefix(ns: list[int]) -> tuple[str, str] | None:
    if len(ns) != 2 or any(n not in _N_STEM_EN for n in ns):
        return None
    a, b = ns
    if a == b:
        return f"N,N-di{_N_STEM_EN[a]}", f"N,N-二{_N_STEM_ZH[a]}"
    x, y = sorted(ns, key=lambda n: _N_STEM_EN[n])
    return f"N-{_N_STEM_EN[x]}-N-{_N_STEM_EN[y]}", f"N-{_N_STEM_ZH[x]}-N-{_N_STEM_ZH[y]}"
def _n_alkyl_prefix(parent: dict) -> tuple[str, str, int] | None:
    kind, n = parent.get("kind"), parent.get("n_alkyl_n")
    if kind in ("sec_amine", "amide") and n in _N_ALKYL_EN:
        return _N_ALKYL_EN[n], _N_ALKYL_ZH[n], n
    if kind in ("tert_amine", "amide"):
        pref = _tert_n_prefix(list(parent.get("n_alkyl_ns") or []))
        return (*pref, 0) if pref else None
    return None
def _extract_n_alkyl(parent: dict) -> list[dict]:
    key = "amide_c_idx" if parent.get("kind") == "amide" else "amine_c_idx"
    attach, pref = parent.get(key), _n_alkyl_prefix(parent)
    return [_n_alkyl_sub(*pref[:2], attach, pref[2])] if attach is not None and pref else []
def _extract_n_phenyl(parent: dict) -> list[dict]:
    if parent.get("kind") != "amide" or not parent.get("n_phenyl"):
        return []
    attach = parent.get("amide_c_idx")
    return [_n_phenyl_sub(attach)] if attach is not None else []


def _n_benzyl_named(info: dict, ch2: int) -> tuple[str, str, bool] | None:
    from namepredict.layer2.aryl_sub import _benzyl_name, _ch2_ph_at
    ams = info.get("amides") or []
    n_idx = ams[0]["n_idx"] if ams else -1
    ph = _ch2_ph_at(info["mol"], ch2, n_idx)
    return None if ph is None else _benzyl_name(info["mol"], ph, ch2)


def _extract_n_benzyl(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") != "amide" or not parent.get("n_benzyl"):
        return []
    attach, ch2 = parent.get("amide_c_idx"), parent.get("n_benzyl_ch2")
    if attach is None or ch2 is None:
        return []
    named = _n_benzyl_named(info, ch2)
    return [_n_benzyl_sub(attach, *named)] if named else []


def _make_aryl(kind: str, attach: int, atoms: list[int], en: str, zh: str, paren: bool) -> dict:
    n = 7 if kind in ("benzyl", "benzyloxy") else 6
    return {
        "kind": kind, "attach_idx": attach, "atoms": atoms,
        "n_carbons": n, "en": en, "zh": zh, "paren": paren,
    }


def _one_phenoxy(mol: Mol, p: dict) -> dict:
    en, zh, paren = _phenoxy_name(mol, p["ph"], p["outer_c"])
    return _make_aryl("phenoxy", p["ring_c"], p["atoms"], en, zh, paren)


def _one_phenyl(mol: Mol, p: dict) -> dict:
    en, zh, paren = _phenyl_name(mol, p["ph"], p["outer_c"])
    return _make_aryl("phenyl", p["attach"], p["atoms"], en, zh, paren)


def _one_benzyl(mol: Mol, p: dict) -> dict:
    en, zh, paren = _benzyl_name(mol, p["ph"], p["ch2"])
    return _make_aryl("benzyl", p["attach"], p["atoms"], en, zh, paren)


def _one_benzyloxy(mol: Mol, p: dict) -> dict:
    en, zh, paren = _benzyloxy_name(mol, p["ph"], p["ch2"])
    return _make_aryl("benzyloxy", p["ring_c"], p["atoms"], en, zh, paren)


def _extract_phenoxys(info: dict, parent: dict) -> list[dict]:
    mol: Mol = info["mol"]
    chain = set(parent.get("chain") or [])
    return [_one_phenoxy(mol, p) for p in _ring_phenoxys(info, chain)]


def _extract_phenyls(info: dict, parent: dict) -> list[dict]:
    mol: Mol = info["mol"]
    chain = set(parent.get("chain") or [])
    return [_one_phenyl(mol, p) for p in _ring_phenyls(mol, chain)]


def _extract_benzyls(info: dict, parent: dict) -> list[dict]:
    mol: Mol = info["mol"]
    chain = set(parent.get("chain") or [])
    return [_one_benzyl(mol, p) for p in _ring_benzyls(mol, chain)]


def _extract_benzyloxys(info: dict, parent: dict) -> list[dict]:
    mol: Mol = info["mol"]
    chain = set(parent.get("chain") or [])
    return [_one_benzyloxy(mol, p) for p in _ring_benzyloxys(info, chain)]


def _aryl_outer_starts(info: dict, parent: dict) -> set[int]:
    from namepredict.layer2.heteroaryl_sub import ring_pyridinyls
    mol: Mol = info["mol"]
    chain = set(parent.get("chain") or [])
    ph = {p["outer_c"] for p in _ring_phenyls(mol, chain)}
    bn = {p["outer_c"] for p in _ring_benzyls(mol, chain)}
    py = {p["outer_c"] for p in ring_pyridinyls(mol, chain)}
    return ph | bn | py


def _extract_alkyls_no_aryl(mol: Mol, chain: list[int], skip: set[int]) -> list[dict]:
    chain_set = set(chain)
    out: list[dict] = []
    for attach, start in _side_starts(mol, chain):
        if start in skip:
            continue
        one = _one_alkyl(mol, attach, start, chain_set)
        if one is not None:
            out.append(one)
    return out


def _extract_core_subs(info: dict, parent: dict) -> list:
    mol: Mol = info["mol"]
    chain = parent.get("chain") or []
    halo = _filter_fg_halos(_extract_halos(mol, chain), parent)
    return (
        halo + _extract_hydroxys(info, parent) + _extract_aminos(info, parent)
        + _extract_oxos(info, parent) + _extract_nitros(info, parent)
    )


def _extract_pyridinyls(info: dict, parent: dict) -> list[dict]:
    from namepredict.layer2.heteroaryl_sub import ring_pyridinyls
    mol: Mol = info["mol"]
    chain = set(parent.get("chain") or [])
    out: list[dict] = []
    for p in ring_pyridinyls(mol, chain):
        out.append({
            "kind": "pyridinyl", "attach_idx": p["attach"], "atoms": p["atoms"],
            "n_carbons": 5, "en": p["en"], "zh": p["zh"], "paren": p["paren"],
        })
    return out


def _extract_aryls(info: dict, parent: dict) -> list[dict]:
    return (
        _extract_phenoxys(info, parent) + _extract_phenyls(info, parent)
        + _extract_benzyloxys(info, parent) + _extract_benzyls(info, parent)
        + _extract_pyridinyls(info, parent)
    )


def _extract_n_subs(info: dict, parent: dict) -> list[dict]:
    return (
        _extract_n_alkyl(parent) + _extract_n_phenyl(parent)
        + _extract_n_benzyl(info, parent)
    )


def extract_substituents(info: dict, parent: dict) -> list:
    mol: Mol = info["mol"]
    chain = parent.get("chain") or []
    alkyl = _extract_alkyls_no_aryl(mol, chain, _aryl_outer_starts(info, parent))
    return (
        alkyl + _extract_core_subs(info, parent) + _extract_alkoxys(info, parent)
        + _extract_aryls(info, parent) + _extract_n_subs(info, parent)
    )
