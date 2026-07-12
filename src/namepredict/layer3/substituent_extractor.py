from __future__ import annotations

from rdkit.Chem import Mol

ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}


def _c_neighbors(mol: Mol, idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]


def _is_pure_alkyl_c(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    return all(n.GetAtomicNum() in (1, 6) for n in atom.GetNeighbors())


def _side_starts(mol: Mol, chain: list[int]) -> list[tuple[int, int]]:
    chain_set = set(chain)
    out: list[tuple[int, int]] = []
    for c in chain:
        for nb in _c_neighbors(mol, c):
            if nb not in chain_set:
                out.append((c, nb))
    return out


def _nb_kind(n: int, prev: int | None, chain_set: set[int], start: int, cur: int) -> str:
    if n == prev:
        return "skip"
    if n in chain_set:
        return "ok" if cur == start else "bad"
    return "free"


def _free_neighbors(
    mol: Mol, cur: int, prev: int | None, chain_set: set[int], start: int
) -> list[int] | None:
    free: list[int] = []
    for n in _c_neighbors(mol, cur):
        kind = _nb_kind(n, prev, chain_set, start, cur)
        if kind == "bad":
            return None
        if kind == "free":
            free.append(n)
    return free


def _next_atom(
    mol: Mol, cur: int, prev: int | None, chain_set: set[int], start: int
) -> int | None | bool:
    """Return next carbon idx, None if end, False if invalid."""
    free = _free_neighbors(mol, cur, prev, chain_set, start)
    if free is None or len(free) > 1:
        return False
    return free[0] if free else None


def _advance(
    mol: Mol, cur: int, prev: int | None, chain_set: set[int], start: int
) -> tuple[int, int | None] | None:
    if not _is_pure_alkyl_c(mol, cur):
        return None
    nxt = _next_atom(mol, cur, prev, chain_set, start)
    if nxt is False:
        return None
    return cur, nxt  # type: ignore[return-value]


def _walk_linear(mol: Mol, start: int, chain_set: set[int]) -> list[int] | None:
    path: list[int] = []
    prev: int | None = None
    cur: int | None = start
    while cur is not None and len(path) < 5:
        step = _advance(mol, cur, prev, chain_set, start)
        if step is None:
            return None
        path.append(step[0])
        prev, cur = step[0], step[1]
    return path if 1 <= len(path) <= 4 else None


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
        path = _walk_linear(mol, start, chain_set)
        if path:
            out.append(_make_alkyl(attach, path))
    return out


def _filter_fg_halos(halos: list, parent: dict) -> list:
    if parent.get("kind") != "acyl_chloride":
        return halos
    cl = parent.get("cl_idx")
    return [h for h in halos if cl not in (h.get("atoms") or [])]


_PARENT_OH_KINDS = frozenset(
    {"alcohol", "alkenol", "diol", "triol", "cycloalcohol", "phenol"}
)
_PARENT_NH2_KINDS = frozenset(
    {"amine", "diamine", "cycloamine", "sec_amine", "tert_amine", "aniline"}
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
def extract_substituents(info: dict, parent: dict) -> list:
    mol: Mol = info["mol"]
    chain = parent.get("chain") or []
    halo = _filter_fg_halos(_extract_halos(mol, chain), parent)
    oh = _extract_hydroxys(info, parent)
    nh2 = _extract_aminos(info, parent)
    oxo = _extract_oxos(info, parent)
    return _extract_alkyls(mol, chain) + halo + oh + nh2 + oxo + _extract_n_alkyl(parent)
