"""Unsubstituted monocyloalkyl side topology C3–C8 (L2; P-29.6)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.side_alkyl import _c_neighbors, _is_pure_alkyl_c, _is_terminal_methyl

_CYCLOALKYL_EN = {
    3: "cyclopropyl", 4: "cyclobutyl", 5: "cyclopentyl",
    6: "cyclohexyl", 7: "cycloheptyl", 8: "cyclooctyl",
}
_CYCLOALKYL_ZH = {
    3: "环丙基", 4: "环丁基", 5: "环戊基",
    6: "环己基", 7: "环庚基", 8: "环辛基",
}


def _bond_single(mol: Mol, a: int, b: int) -> bool:
    bond = mol.GetBondBetweenAtoms(a, b)
    return bond is not None and bond.GetBondType().name == "SINGLE"


def _ring_all_c_single(mol: Mol, atoms: tuple | list) -> bool:
    ids = list(atoms)
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in ids):
        return False
    for i, a in enumerate(ids):
        if not _bond_single(mol, a, ids[(i + 1) % len(ids)]):
            return False
    return True


def _sat_carbocycles(mol: Mol) -> list[set[int]]:
    out: list[set[int]] = []
    for r in mol.GetRingInfo().AtomRings():
        if 3 <= len(r) <= 8 and _ring_all_c_single(mol, r):
            out.append(set(r))
    return out


def _ring_unfused(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _exo_nbs(mol: Mol, i: int, ring: set[int]) -> list[int]:
    return [
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in ring
    ]


def _cyclo_exo_ok(mol: Mol, ring: set[int], attach: int, parent: set[int]) -> bool:
    for i in ring:
        for j in _exo_nbs(mol, i, ring):
            if not (i == attach and j in parent):
                return False
    return True


def _try_side_ring(
    mol: Mol, start: int, parent: set[int], ring: set[int],
) -> set[int] | None:
    if start not in ring or ring & parent:
        return None
    if not _ring_unfused(mol, ring):
        return None
    return ring if _cyclo_exo_ok(mol, ring, start, parent) else None


def _cycloalkyl_ring_at(mol: Mol, start: int, parent: set[int]) -> set[int] | None:
    if start in parent:
        return None
    for ring in _sat_carbocycles(mol):
        hit = _try_side_ring(mol, start, parent, ring)
        if hit is not None:
            return hit
    return None


def _is_monocycloalkyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    ring = _cycloalkyl_ring_at(mol, start, chain)
    return sorted(ring) if ring is not None else None


def _cycloalkyl_names(n: int) -> tuple[str, str] | None:
    en, zh = _CYCLOALKYL_EN.get(n), _CYCLOALKYL_ZH.get(n)
    return (en, zh) if en and zh else None


def _split_me_and_cyc(mol: Mol, start: int, a: int, b: int) -> tuple[int, int] | None:
    if _is_terminal_methyl(mol, a, start):
        return a, b
    if _is_terminal_methyl(mol, b, start):
        return b, a
    return None


def _bridge_free_pair(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    free = [n for n in _c_neighbors(mol, start) if n not in chain]
    return free if len(free) == 2 else None


def _ethyl_bridge_parts(
    mol: Mol, start: int, chain: set[int],
) -> tuple[int, int, set[int]] | None:
    free = _bridge_free_pair(mol, start, chain)
    if free is None:
        return None
    pair = _split_me_and_cyc(mol, start, free[0], free[1])
    if pair is None:
        return None
    me, cyc_c = pair
    ring = _cycloalkyl_ring_at(mol, cyc_c, chain | {start, me})
    return (me, cyc_c, ring) if ring else None


def _is_1_cycloalkylethyl(mol: Mol, start: int, chain: set[int]) -> list[int] | None:
    """parent–CH(CH3)–cycloalkyl (unsubstituted monocyloalkyl)."""
    if not _is_pure_alkyl_c(mol, start) or start in chain:
        return None
    parts = _ethyl_bridge_parts(mol, start, chain)
    if parts is None:
        return None
    me, _cyc, ring = parts
    return [start, me, *sorted(ring)]
