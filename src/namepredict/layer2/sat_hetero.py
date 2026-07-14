"""Saturated monocyclic hetero parents (IUPAC P-22.2.2 / P-22.2.3 / P-25).

Unsubstituted or ring-C simple alkyl / mono-halo Hantzsch–Widman / retained
saturated monoheterocycles: aziridine, oxirane, oxolane, oxane, pyrrolidine,
piperidine, morpholine, piperazine, 1,3-dioxolane, 1,4-dioxane, optional
thiolane. Ring-C sides: linear n-alkyl C1–C4 (multi) + retained branched
(tert-butyl / isopropyl / …); N-alkyl out of scope.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _disjoint_cover,
    _outside_carbons,
    _outside_ok,
    _ring_bonds_single,
    _ring_halo_n,
    _ring_side_starts,
    _side_sets,
)

# (size, z_counts frozenset of (Z, count)) → kind
# z_counts as frozenset of (atomic_num, count) for non-C ring atoms
_MONO = {
    (3, frozenset([(7, 1)])): "aziridine",
    (3, frozenset([(8, 1)])): "oxirane",
    (5, frozenset([(8, 1)])): "oxolane",
    (5, frozenset([(7, 1)])): "pyrrolidine",
    (5, frozenset([(16, 1)])): "thiolane",
    (6, frozenset([(8, 1)])): "oxane",
    (6, frozenset([(7, 1)])): "piperidine",
    (6, frozenset([(8, 1), (7, 1)])): "morpholine",
    (6, frozenset([(7, 2)])): "piperazine",
    (5, frozenset([(8, 2)])): "dioxolane",
    (6, frozenset([(8, 2)])): "dioxane",
}


def _ring_atoms_if_mono(info: dict) -> list[int] | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    return list(rings[0]["atom_ids"])


def _ring_zs(mol: Mol, atom_ids: list[int]) -> list[int]:
    return [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids]


def _z_counts(zs: list[int]) -> frozenset[tuple[int, int]]:
    from collections import Counter
    c = Counter(z for z in zs if z != 6)
    return frozenset(c.items())


def _all_nonarom(mol: Mol, atom_ids: list[int]) -> bool:
    return all(not mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids)


def _only_c_nos(zs: list[int]) -> bool:
    return all(z in (6, 7, 8, 16) for z in zs) and any(z != 6 for z in zs)


def _is_sat_hetero_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None:
        return False
    mol: Mol = info["mol"]
    if not _all_nonarom(mol, atom_ids):
        return False
    if not _ring_bonds_single(mol, tuple(atom_ids)):
        return False
    return _only_c_nos(_ring_zs(mol, atom_ids))


def _hetero_sig(info: dict) -> tuple[int, frozenset] | None:
    if not _is_sat_hetero_core(info):
        return None
    atom_ids = _ring_atoms_if_mono(info)
    zs = _ring_zs(info["mol"], atom_ids)
    return len(atom_ids), _z_counts(zs)


def _base_kind(info: dict) -> str | None:
    sig = _hetero_sig(info)
    return None if sig is None else _MONO.get(sig)


def _ring_hetero_idxs(mol: Mol, atom_ids: list[int]) -> list[int]:
    return [i for i in atom_ids if mol.GetAtomWithIdx(i).GetAtomicNum() != 6]


def _ring_nn_dist(atom_ids: list[int], a: int, b: int) -> int:
    ia, ib = atom_ids.index(a), atom_ids.index(b)
    d = abs(ia - ib)
    return min(d, len(atom_ids) - d)


def _dihetero_dist_ok(info: dict, kind: str) -> bool:
    """1,3-dioxolane / 1,4-dioxane / morpholine O–N=3 / piperazine N–N=3."""
    atom_ids = _ring_atoms_if_mono(info)
    mol: Mol = info["mol"]
    hs = _ring_hetero_idxs(mol, atom_ids)
    if len(hs) != 2:
        return True
    d = _ring_nn_dist(atom_ids, hs[0], hs[1])
    want = {"dioxolane": 2, "dioxane": 3, "morpholine": 3, "piperazine": 3}
    return d == want.get(kind, d)


def _kind_of(info: dict) -> str | None:
    kind = _base_kind(info)
    if kind is None:
        return None
    return kind if _dihetero_dist_ok(info, kind) else None


def _fg_block(info: dict) -> bool:
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_ester") or info.get("has_amide"):
        return True
    return bool(info.get("has_nitrile") or info.get("has_thiol") or info.get("has_nitro"))


def _ring_n_set(mol: Mol, atom_ids: list[int]) -> set[int]:
    return {i for i in atom_ids if mol.GetAtomWithIdx(i).GetAtomicNum() == 7}


def _extra_amine(info: dict, ring_ns: set[int]) -> bool:
    """True if any detected amine is not a ring hetero (principal amine FG)."""
    for a in info.get("amines") or []:
        if a.get("n_idx") not in ring_ns:
            return True
        if a.get("degree") == 1:
            return True
    return False


def _n_unsub(mol: Mol, ring_ns: set[int]) -> bool:
    """Ring N must be N–H (or bare), not N-alkyl."""
    for i in ring_ns:
        if mol.GetAtomWithIdx(i).GetTotalNumHs() < 1:
            return False
    return True


def _mono_methyl_only(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    """True when every outside C is a mono-methyl attach (used by carboxylic)."""
    if not starts:
        return True
    outside = set(_outside_carbons(mol, ring))
    return outside == set(starts)


def _sides_cover(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    """Outside C fully covered by recognized alkyl side probes (C1–C4 / branched)."""
    outside = set(_outside_carbons(mol, ring))
    if not starts:
        return not outside
    sets = _side_sets(mol, ring, starts)
    return sets is not None and _disjoint_cover(sets, outside)


def _subs_ok(mol: Mol, ring: set[int]) -> bool:
    """Allow ring-C linear/branched alkyl sides + monohalo (P-22.2.2 / P-14.3.4)."""
    if _ring_halo_n(mol, ring) > 1:
        return False
    return _sides_cover(mol, ring, _ring_side_starts(mol, ring))


def _is_simple_sat_hetero(info: dict) -> bool:
    if _kind_of(info) is None or _fg_block(info):
        return False
    mol: Mol = info["mol"]
    atom_ids = _ring_atoms_if_mono(info)
    ring, ring_ns = set(atom_ids), _ring_n_set(mol, atom_ids)
    if _extra_amine(info, ring_ns) or not _n_unsub(mol, ring_ns):
        return False
    if not _outside_ok(mol, ring):
        return False
    return _subs_ok(mol, ring)


def _hetero_pair(info: dict) -> list[int]:
    """Hetero indices for orientation: mono → [h]; di → ordered by priority."""
    mol: Mol = info["mol"]
    atom_ids = _ring_atoms_if_mono(info)
    hs = _ring_hetero_idxs(mol, atom_ids)
    if len(hs) <= 1:
        return hs
    # O before N before S for morpholine O=1; piperazine both N; diox both O
    return sorted(hs, key=lambda i: (mol.GetAtomWithIdx(i).GetAtomicNum() != 8, i))


def _hetero_asym(info: dict, hs: list[int]) -> bool:
    """True when dihetero atoms differ (morpholine O≠N); false for N–N / O–O."""
    if len(hs) < 2:
        return False
    mol: Mol = info["mol"]
    return mol.GetAtomWithIdx(hs[0]).GetAtomicNum() != mol.GetAtomWithIdx(hs[1]).GetAtomicNum()


def _sat_hetero_parent(info: dict) -> dict:
    atom_ids = _ring_atoms_if_mono(info)
    kind = _kind_of(info)
    hs = _hetero_pair(info)
    out = {
        "chain": atom_ids, "n_carbons": len(atom_ids), "kind": kind,
        "hetero_idxs": hs, "hetero_asym": _hetero_asym(info, hs),
    }
    if len(hs) == 1:
        out["hetero_idx"] = hs[0]
    return out


def _try_sat_hetero_parent(info: dict) -> dict | None:
    return _sat_hetero_parent(info) if _is_simple_sat_hetero(info) else None
