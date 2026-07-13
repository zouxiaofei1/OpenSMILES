from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_alkyl import (
    _c_neighbors,
    _is_isobutyl,
    _is_isopentyl,
    _is_isopropyl,
    _is_neopentyl,
    _is_sec_butyl,
    _is_tert_butyl,
    _is_trifluoromethyl,
    _walk_linear,
)


def alkyl_alpha_key(stem: str) -> str:
    """Alphanumerical-order key: italic prefixes sec-/tert- ignored (P-14.5)."""
    if stem.startswith("tert-") or stem.startswith("sec-"):
        return stem[stem.index("-") + 1 :]
    return stem

ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}


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


def _one_alkyl(mol: Mol, attach: int, start: int, chain_set: set[int]) -> dict | None:
    path = _walk_linear(mol, start, chain_set)
    if path:
        return _make_alkyl(attach, path)
    return _one_branched(mol, attach, start, chain_set)


_BRANCH_CHECKS = (
    (_is_isopropyl, 3, "isopropyl", "异丙基"),
    (_is_tert_butyl, 4, "tert-butyl", "叔丁基"),
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
        "phenol", "benzenediol", "pyridinol",
    }
)
_PARENT_NH2_KINDS = frozenset(
    {
        "amine", "diamine", "cycloamine", "sec_amine", "tert_amine",
        "aniline", "pyridinamine",
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


def _make_amino(attach: int, n_idx: int) -> dict:
    return {
        "kind": "amino",
        "attach_idx": attach,
        "atoms": [n_idx],
        "en": "amino",
        "zh": "氨基",
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
    if parent.get("kind") in _PARENT_NH2_KINDS:
        return []
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for a in info.get("amines") or []:
        if "c_idx" in a and a["c_idx"] in chain_set:
            out.append(_make_amino(a["c_idx"], a["n_idx"]))
    return out


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


_ALKOXY_EN = {1: "methoxy", 2: "ethoxy"}
_ALKOXY_ZH = {1: "甲氧基", 2: "乙氧基"}


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _outer_fwd(mol: Mol, cur: int, prev: int) -> list:
    atom = mol.GetAtomWithIdx(cur)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return []
    return [x for x in _heavies(atom) if x.GetIdx() != prev]


def _alkoxy_n(mol: Mol, o_idx: int, outer_c: int) -> int:
    fwd = _outer_fwd(mol, outer_c, o_idx)
    if not fwd:
        return 1
    if len(fwd) == 1 and fwd[0].GetAtomicNum() == 6:
        return 2 if not _outer_fwd(mol, fwd[0].GetIdx(), outer_c) else 0
    return 0


def _make_alkoxy(attach: int, o_idx: int, atoms: list[int], n: int) -> dict:
    return {
        "kind": "alkoxy", "attach_idx": attach, "atoms": [o_idx] + atoms,
        "n_carbons": n, "en": _ALKOXY_EN[n], "zh": _ALKOXY_ZH[n],
    }


def _alkoxy_atoms(mol: Mol, outer: int, o_idx: int, n: int) -> list[int]:
    if n == 1:
        return [outer]
    fwd = _outer_fwd(mol, outer, o_idx)
    return [outer, fwd[0].GetIdx()] if fwd else [outer]


def _one_ring_alkoxy(mol: Mol, e: dict, chain_set: set[int]) -> dict | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in chain_set) == (c2 in chain_set):
        return None
    ring_c, outer = (c1, c2) if c1 in chain_set else (c2, c1)
    if not mol.GetAtomWithIdx(ring_c).IsInRing():
        return None
    n = _alkoxy_n(mol, o, outer)
    if n not in _ALKOXY_EN:
        return None
    return _make_alkoxy(ring_c, o, _alkoxy_atoms(mol, outer, o, n), n)


def _extract_alkoxys(info: dict, parent: dict) -> list[dict]:
    mol: Mol = info["mol"]
    chain_set = set(parent.get("chain") or [])
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _one_ring_alkoxy(mol, e, chain_set)
        if one is not None:
            out.append(one)
    return out


def extract_substituents(info: dict, parent: dict) -> list:
    mol: Mol = info["mol"]
    chain = parent.get("chain") or []
    halo = _filter_fg_halos(_extract_halos(mol, chain), parent)
    oh = _extract_hydroxys(info, parent)
    nh2 = _extract_aminos(info, parent)
    oxo = _extract_oxos(info, parent)
    nitro = _extract_nitros(info, parent)
    alkox = _extract_alkoxys(info, parent)
    n_sub = _extract_n_alkyl(parent) + _extract_n_phenyl(parent)
    return _extract_alkyls(mol, chain) + halo + oh + nh2 + oxo + nitro + alkox + n_sub
