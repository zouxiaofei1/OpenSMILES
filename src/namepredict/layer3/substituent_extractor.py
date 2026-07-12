from __future__ import annotations

from rdkit.Chem import Mol

ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


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


def extract_substituents(info: dict, parent: dict) -> list:
    mol: Mol = info["mol"]
    chain = parent.get("chain") or []
    chain_set = set(chain)
    out: list = []
    for attach, start in _side_starts(mol, chain):
        path = _walk_linear(mol, start, chain_set)
        if path:
            out.append(_make_alkyl(attach, path))
    return out
