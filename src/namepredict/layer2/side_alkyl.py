"""Mono alkyl side-chain topology: linear C1–C4 + retained branched (L2/L3)."""
from __future__ import annotations

from rdkit.Chem import Mol

_HALO_Z = frozenset({9, 17, 35, 53})


def _c_neighbors(mol: Mol, idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]


def _f_neighbors(mol: Mol, idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 9]


def _free_c(mol: Mol, idx: int, blocked: set[int]) -> list[int]:
    return [n for n in _c_neighbors(mol, idx) if n not in blocked]


def _is_pure_alkyl_c(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    return all(n.GetAtomicNum() in (1, 6) for n in atom.GetNeighbors())


def _is_cf3_carbon(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    nbs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(nbs) == 4 and sum(n.GetAtomicNum() == 9 for n in nbs) == 3


def _is_cf3_fluoro(atom) -> bool:
    if atom.GetAtomicNum() != 9:
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    if len(heavies) != 1 or heavies[0].GetAtomicNum() != 6:
        return False
    c = heavies[0]
    if c.IsInRing():
        return False
    nbs = [n for n in c.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(nbs) == 4 and sum(n.GetAtomicNum() == 9 for n in nbs) == 3


def _is_side_halo(atom) -> bool:
    """Terminal halogen on a non-ring carbon (ω-haloalkyl)."""
    if atom.GetAtomicNum() not in _HALO_Z:
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    if len(heavies) != 1 or heavies[0].GetAtomicNum() != 6:
        return False
    return not heavies[0].IsInRing()


def _terminal_halo_z(mol: Mol, idx: int) -> int | None:
    """Atomic number of the single terminal halo on carbon, else None."""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return None
    halos = [n for n in atom.GetNeighbors() if n.GetAtomicNum() in _HALO_Z]
    if len(halos) != 1:
        return None
    ok = all(n.GetAtomicNum() in (1, 6) or n.GetAtomicNum() in _HALO_Z for n in atom.GetNeighbors())
    return halos[0].GetAtomicNum() if ok else None


def _is_omega_halo_c(mol: Mol, idx: int) -> bool:
    return _terminal_halo_z(mol, idx) is not None


def _is_trifluoromethyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    free = _free_c(mol, start, chain)
    if free or not _is_cf3_carbon(mol, start):
        return None
    return [start] if len(_f_neighbors(mol, start)) == 3 else None


def _nb_kind(n: int, prev: int | None, chain: set[int], start: int, cur: int) -> str:
    if n == prev:
        return "skip"
    if n in chain:
        return "ok" if cur == start else "bad"
    return "free"


def _free_neighbors(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> list[int] | None:
    free: list[int] = []
    for n in _c_neighbors(mol, cur):
        kind = _nb_kind(n, prev, chain, start, cur)
        if kind == "bad":
            return None
        if kind == "free":
            free.append(n)
    return free


def _next_atom(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> int | None | bool:
    free = _free_neighbors(mol, cur, prev, chain, start)
    if free is None or len(free) > 1:
        return False
    return free[0] if free else None


def _advance(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> tuple[int, int | None] | None:
    if not _is_pure_alkyl_c(mol, cur):
        return None
    nxt = _next_atom(mol, cur, prev, chain, start)
    if nxt is False:
        return None
    return cur, nxt  # type: ignore[return-value]


def _omega_step(
    mol: Mol, cur: int, prev: int | None, chain: set[int], start: int,
) -> tuple[int, int | None] | None:
    """One step along ω-halo n-alkyl: pure C or terminal mono-halo C."""
    free = _free_neighbors(mol, cur, prev, chain, start)
    if free is None or len(free) > 1:
        return None
    if free:
        return (cur, free[0]) if _is_pure_alkyl_c(mol, cur) else None
    return (cur, None) if _is_omega_halo_c(mol, cur) else None


def _walk_with(
    mol: Mol, start: int, chain: set[int], step_fn, max_len: int = 5,
) -> list[int] | None:
    path: list[int] = []
    prev: int | None = None
    cur: int | None = start
    while cur is not None and len(path) < max_len:
        step = step_fn(mol, cur, prev, chain, start)
        if step is None:
            return None
        path.append(step[0])
        prev, cur = step[0], step[1]
    return path


def _walk_omega_halo(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    path = _walk_with(mol, start, chain, _omega_step)
    if path is None or not (1 <= len(path) <= 4):
        return None
    return path if _is_omega_halo_c(mol, path[-1]) else None


def _walk_linear(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    path = _walk_with(mol, start, chain, _advance)
    return path if path is not None and 1 <= len(path) <= 4 else None


def _walk_linear_n(
    mol: Mol, start: int, chain: set[int], max_n: int = 12,
) -> list[int] | None:
    """Linear pure n-alkyl path of length 1..max_n (C1–C12 for BQ sides)."""
    path = _walk_with(mol, start, chain, _advance, max_len=max_n + 1)
    return path if path is not None and 1 <= len(path) <= max_n else None


def _is_terminal_methyl(mol: Mol, idx: int, parent: int) -> bool:
    if not _is_pure_alkyl_c(mol, idx):
        return False
    return _c_neighbors(mol, idx) == [parent]


def _two_methyls(mol: Mol, idx: int, blocked: set[int]) -> list[int] | None:
    free = _free_c(mol, idx, blocked)
    if len(free) != 2:
        return None
    if not all(_is_terminal_methyl(mol, m, idx) for m in free):
        return None
    return free


def _three_methyls(mol: Mol, idx: int, blocked: set[int]) -> list[int] | None:
    free = _free_c(mol, idx, blocked)
    if len(free) != 3:
        return None
    if not all(_is_terminal_methyl(mol, m, idx) for m in free):
        return None
    return free


def _is_ethyl_arm(mol: Mol, idx: int, parent: int) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, idx):
        return None
    free = _free_c(mol, idx, {parent})
    if len(free) != 1 or not _is_terminal_methyl(mol, free[0], idx):
        return None
    return [idx, free[0]]


def _is_isopropyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    mids = _two_methyls(mol, start, chain)
    return [start, *mids] if mids else None


def _is_tert_butyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    mids = _three_methyls(mol, start, chain)
    return [start, *mids] if mids else None


def _class_tert_arm(mol: Mol, idx: int, parent: int) -> str | None:
    """Classify free C off quaternary attach: 'Me' or 'Et'."""
    if _is_terminal_methyl(mol, idx, parent):
        return "Me"
    if _is_ethyl_arm(mol, idx, parent):
        return "Et"
    return None


def _tert_pentyl_arms(mol: Mol, start: int, free: list[int]) -> list[int] | None:
    """Two Me + one Et off start → atom list without start."""
    kinds = [_class_tert_arm(mol, c, start) for c in free]
    if kinds.count("Me") != 2 or kinds.count("Et") != 1:
        return None
    eth_i = free[kinds.index("Et")]
    eth = _is_ethyl_arm(mol, eth_i, start)
    mes = [c for c, k in zip(free, kinds) if k == "Me"]
    return [*mes, *eth] if eth else None


def _is_2_methylbutan_2_yl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    """parent–C(Me)(Me)–CH2–CH3 (2-methylbutan-2-yl / tert-pentyl)."""
    if not _is_pure_alkyl_c(mol, start):
        return None
    free = _free_c(mol, start, chain)
    if len(free) != 3:
        return None
    arms = _tert_pentyl_arms(mol, start, free)
    return [start, *arms] if arms else None


def _is_tert_pentyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    """Alias for retained-name call sites; same topology as 2-methylbutan-2-yl."""
    return _is_2_methylbutan_2_yl(mol, start, chain)


def _is_isobutyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    mid = _one_free_c(mol, start, chain)
    if mid is None:
        return None
    mids = _two_methyls(mol, mid, {start})
    return [start, mid, *mids] if mids else None


def _sec_from_arms(mol: Mol, start: int, a: int, b: int) -> list[int] | None:
    if _is_terminal_methyl(mol, a, start):
        eth = _is_ethyl_arm(mol, b, start)
        return [start, a, *eth] if eth else None
    if _is_terminal_methyl(mol, b, start):
        eth = _is_ethyl_arm(mol, a, start)
        return [start, b, *eth] if eth else None
    return None


def _is_sec_butyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    free = _free_c(mol, start, chain)
    if len(free) != 2:
        return None
    return _sec_from_arms(mol, start, free[0], free[1])


def _is_neopentyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    mid = _one_free_c(mol, start, chain)
    if mid is None:
        return None
    mids = _three_methyls(mol, mid, {start})
    return [start, mid, *mids] if mids else None


def _one_free_c(mol: Mol, idx: int, blocked: set[int]) -> int | None:
    free = _free_c(mol, idx, blocked)
    if len(free) != 1 or not _is_pure_alkyl_c(mol, free[0]):
        return None
    return free[0]


def _is_isopentyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    if not _is_pure_alkyl_c(mol, start):
        return None
    mid = _one_free_c(mol, start, chain)
    if mid is None:
        return None
    outer = _one_free_c(mol, mid, {start})
    if outer is None:
        return None
    mids = _two_methyls(mol, outer, {mid})
    return [start, mid, outer, *mids] if mids else None


def _has_double(mol: Mol, a: int, b: int) -> bool:
    bond = mol.GetBondBetweenAtoms(a, b)
    return bond is not None and bond.GetBondType().name == "DOUBLE"


def _one_pure_free(mol: Mol, idx: int, blocked: set[int]) -> int | None:
    free = _free_c(mol, idx, blocked)
    if len(free) != 1 or not _is_pure_alkyl_c(mol, free[0]):
        return None
    return free[0]


def _prenyl_outer(mol: Mol, mid: int, start: int) -> int | None:
    free = _free_c(mol, mid, {start})
    if len(free) != 1 or not _has_double(mol, mid, free[0]):
        return None
    return free[0]


def _is_prenyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    """Ring–CH2–CH=C(CH3)2 (3-methylbut-2-enyl)."""
    if not _is_pure_alkyl_c(mol, start):
        return None
    mid = _one_pure_free(mol, start, chain)
    outer = _prenyl_outer(mol, mid, start) if mid is not None else None
    if mid is None or outer is None:
        return None
    mids = _two_methyls(mol, outer, {mid})
    return [start, mid, outer, *mids] if mids else None


def _probe_monocycloalkyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    from namepredict.layer2.side_cycloalkyl import _is_monocycloalkyl
    return _is_monocycloalkyl(mol, start, chain)


def _probe_1_cycloalkylethyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    from namepredict.layer2.side_cycloalkyl import _is_1_cycloalkylethyl
    return _is_1_cycloalkylethyl(mol, start, chain)


# Topology-only side recognition for parent gates (not a naming-mode table).
_TOPOLOGY_SIDE_PROBES = (
    _walk_linear,
    _walk_omega_halo,
    _is_isopropyl,
    _is_tert_butyl,
    _is_2_methylbutan_2_yl,
    _is_isobutyl,
    _is_sec_butyl,
    _is_neopentyl,
    _is_prenyl,
    _is_isopentyl,
    _is_trifluoromethyl,
    _probe_1_cycloalkylethyl,
    _probe_monocycloalkyl,
)


def _first_side(mol: Mol, start: int, chain: set[int], probes) -> list[int] | None:
    for probe in probes:
        got = probe(mol, start, chain)
        if got is not None:
            return got
    return None


def _side_atoms(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    return _first_side(mol, start, chain, _TOPOLOGY_SIDE_PROBES)


def _side_covers(
    mol: Mol, start: int, chain: set[int], outside: set[int],
) -> bool:
    atoms = _side_atoms(mol, start, chain)
    return atoms is not None and set(atoms) == outside


def _side_sets(mol: Mol, chain: set[int], starts: list[int]) -> list[set[int]] | None:
    sets: list[set[int]] = []
    for s in starts:
        atoms = _side_atoms(mol, s, chain)
        if atoms is None:
            return None
        sets.append(set(atoms))
    return sets


def _disjoint_cover(sets: list[set[int]], outside: set[int]) -> bool:
    seen: set[int] = set()
    for s in sets:
        if seen & s:
            return False
        seen |= s
    return seen == outside


def _outside_c_atoms(mol: Mol, ring: set[int]) -> set[int]:
    return {
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    }


def _linear_path_ok(mol: Mol, start: int, ring: set[int], max_n: int) -> list[int] | None:
    """Linear n-alkyl path from start of length 1..max_n, or None."""
    return _walk_linear_n(mol, start, ring, max_n=max_n)


def _linear_n_alkyl_sides_ok(
    mol: Mol, ring: set[int], starts: list[int], max_n: int,
) -> bool:
    """True when every ring side start is linear n-alkyl C1–max_n covering outside C."""
    if not starts:
        return not _outside_c_atoms(mol, ring)
    paths: list[set[int]] = []
    for s in starts:
        path = _linear_path_ok(mol, s, ring, max_n)
        if path is None:
            return False
        paths.append(set(path))
    return _disjoint_cover(paths, _outside_c_atoms(mol, ring))


# Ring outer alkoxy topology lives in side_alkoxy (shared by ring_parent / L3).
from namepredict.layer2.side_alkoxy import (  # noqa: E402
    _outer_alkoxy_n,
    _outer_atoms,
    _parse_outer_alkoxy,
)
