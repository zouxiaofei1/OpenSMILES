"""Shared fused aromatic 5+6 ring helpers (IUPAC P-22.2.1 / P-25)."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit.Chem import Mol

from namepredict.layer2.heteroarene5 import _ring_nn_dist
from namepredict.layer2.ring_parent import (
    _is_methyl_on_ring,
    _mono_amine_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _fused_56_pair(rings: list[list[int]]) -> tuple[list[int], list[int]] | None:
    fives = [r for r in rings if len(r) == 5]
    sixes = [r for r in rings if len(r) == 6]
    for five in fives:
        for six in sixes:
            if len(set(five) & set(six)) == 2:
                return five, six
    return None


def _two_rings(info: dict) -> tuple[list[int], list[int]] | None:
    """Adjacent ring pair; multi-ring prefers fused 5+6."""
    rings = [list(r["atom_ids"]) for r in (info.get("rings") or [])]
    if len(rings) < 2:
        return None
    if len(rings) == 2:
        return rings[0], rings[1]
    return _fused_56_pair(rings) or (rings[0], rings[1])


def _size_pair(
    r1: list[int], r2: list[int],
) -> tuple[list[int], list[int]] | None:
    """Return (five, six) if sizes are 5 and 6."""
    if len(r1) == 5 and len(r2) == 6:
        return r1, r2
    if len(r1) == 6 and len(r2) == 5:
        return r2, r1
    return None


def _bridge_pair(r1: list[int], r2: list[int]) -> tuple[int, int] | None:
    inter = set(r1) & set(r2)
    if len(inter) != 2:
        return None
    a, b = tuple(inter)
    return a, b


def _bridge_adjacent(mol: Mol, ba: int, bb: int) -> bool:
    return mol.GetBondBetweenAtoms(ba, bb) is not None


def _all_aromatic(mol: Mol, atoms: set[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms)


def _six_all_c(mol: Mol, six: list[int]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in six)


def _fused_bridge(
    five: list[int], six: list[int], mol: Mol,
) -> tuple[int, int] | None:
    bridge = _bridge_pair(five, six)
    if bridge is None or not _bridge_adjacent(mol, *bridge):
        return None
    return bridge


def _fused56(info: dict) -> tuple[list[int], list[int], tuple[int, int]] | None:
    """Return (five, six, bridge) for adjacent 5+6 fusion, else None."""
    pair = _two_rings(info)
    if pair is None:
        return None
    sized = _size_pair(*pair)
    if sized is None:
        return None
    bridge = _fused_bridge(*sized, info["mol"])
    return None if bridge is None else (sized[0], sized[1], bridge)


def _path_along(ring: list[int], start: int, end: int) -> list[int] | None:
    if start not in ring or end not in ring:
        return None
    i, n = ring.index(start), len(ring)
    path: list[int] = []
    for k in range(1, n):
        atom = ring[(i + k) % n]
        if atom == end:
            return path
        path.append(atom)
    return None


def _path_of_len(
    ring: list[int], start: int, end: int, n: int,
) -> list[int] | None:
    for base in (ring, list(reversed(ring))):
        p = _path_along(base, start, end)
        if p is not None and len(p) == n:
            return p
    return None


def _exterior6(six: list[int], a3a: int, a7a: int) -> list[int] | None:
    """Four exterior carbons on the six-ring from 3a toward 7a."""
    return _path_of_len(six, a3a, a7a, 4)


def _hetero_neighbors_on(mol: Mol, h: int, five: list[int]) -> list[int]:
    fset = set(five)
    return [
        n.GetIdx() for n in mol.GetAtomWithIdx(h).GetNeighbors()
        if n.GetIdx() in fset
    ]


def _split_hetero_nb(
    mol: Mol, h: int, five: list[int], bridge: set[int],
) -> tuple[int, int] | None:
    """Return (pos2, a7a) where a7a is bridge neighbor of heteroatom."""
    nbs = _hetero_neighbors_on(mol, h, five)
    if len(nbs) != 2:
        return None
    a, b = nbs
    if a in bridge and b not in bridge:
        return b, a
    if b in bridge and a not in bridge:
        return a, b
    return None


def _chain_atoms(
    mol: Mol, five: list[int], six: list[int], h: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=hetero, 2, 3, 3a, 4, 5, 6, 7, 7a."""
    split = _split_hetero_nb(mol, h, five, {ba, bb})
    if split is None:
        return None
    pos2, a7a = split
    a3a = bb if a7a == ba else ba
    mid = _path_of_len(five, pos2, a3a, 1)
    ext = _exterior6(six, a3a, a7a) if mid else None
    return None if not mid or not ext else [h, pos2, mid[0], a3a] + ext + [a7a]


# --- data-driven mono-hetero fused 5+6 (benzofuran / benzothiophene) ---

_FG_BLOCK_KEYS = (
    "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
    "has_ester", "has_amide", "has_nitrile", "has_amine",
    "has_thiol", "has_nitro", "has_acyl_chloride", "has_anhydride",
)


@dataclass(frozen=True)
class Fused56MonoSpec:
    """Mono-hetero benzo[b] scaffold: one hetero Z on five-ring, four C."""

    kind: str
    hetero_z: int
    hetero_key: str  # parent dict key, e.g. o_idx / s_idx
    sub_cap: int = 1
    n_atoms: int = 9


def _five_one_z(mol: Mol, five: list[int], z: int) -> int | None:
    """Exactly one atom of atomic num z and four C on the five-ring."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(z) != 1 or zs.count(6) != 4:
        return None
    return next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == z)


def _mono_core_ok(
    mol: Mol, five: list[int], six: list[int], spec: Fused56MonoSpec,
) -> int | None:
    atoms = set(five) | set(six)
    if len(atoms) != spec.n_atoms or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_one_z(mol, five, spec.hetero_z)


def _mono_parts(
    info: dict, spec: Fused56MonoSpec,
) -> tuple[list[int], list[int], int, int, int] | None:
    """Return (five, six, hetero_idx, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    h = _mono_core_ok(info["mol"], five, six, spec)
    return None if h is None else (five, six, h, ba, bb)


def _fg_block(info: dict, keys: tuple[str, ...] = _FG_BLOCK_KEYS) -> bool:
    return any(info.get(k) for k in keys)


def _mono_methyl_only(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    if len(starts) != 1:
        return False
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    return outside == starts


def _subs_ok_cap(mol: Mol, ring: set[int], cap: int) -> bool:
    """≤cap simple ring subs: halo + mono-methyl only (no multi-alkyl)."""
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > cap:
        return False
    if not starts:
        return True
    return _mono_methyl_only(mol, ring, starts)


def _ring_set_parts(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _is_simple_mono(info: dict, spec: Fused56MonoSpec) -> bool:
    if _mono_parts(info, spec) is None or _fg_block(info):
        return False
    parts = _mono_parts(info, spec)
    assert parts is not None
    mol, ring = info["mol"], _ring_set_parts(parts)
    if not _outside_ok(mol, ring):
        return False
    return _subs_ok_cap(mol, ring, spec.sub_cap)


def _mono_parent_dict(info: dict, spec: Fused56MonoSpec, **extra) -> dict:
    parts = _mono_parts(info, spec)
    assert parts is not None
    five, six, h, ba, bb = parts
    chain = _chain_atoms(info["mol"], five, six, h, ba, bb) or []
    return {
        "chain": chain, "n_carbons": spec.n_atoms, "kind": spec.kind,
        spec.hetero_key: h, "bridge": [ba, bb], **extra,
    }


def _try_mono_fused56(info: dict, spec: Fused56MonoSpec) -> dict | None:
    """Simple unsubstituted/mono-sub mono-hetero fused 5+6 parent."""
    return _mono_parent_dict(info, spec) if _is_simple_mono(info, spec) else None


# --- data-driven 1,3-dihetero fused 5+6 (benzothiazole / benzoxazole) ---

_AMINE_BLOCK_KEYS = tuple(k for k in _FG_BLOCK_KEYS if k != "has_amine")


@dataclass(frozen=True)
class Fused56Di13Spec:
    """1,3-benzoazole: five-ring 1×Z + 1×N + 3×C, Z–N dist=2."""

    kind: str
    hetero_z: int  # 8=O, 16=S
    hetero_key: str  # o_idx / s_idx
    amine_kind: str
    sub_cap: int = 2
    amine_sub_cap: int = 1
    n_atoms: int = 9


def _five_zn(mol: Mol, five: list[int], z: int) -> tuple[int, int] | None:
    """Exactly one Z, one N, three C; Z–N ring dist 2. Return (z_idx, n_idx)."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(z) != 1 or zs.count(7) != 1 or zs.count(6) != 3:
        return None
    h = next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == z)
    n = next(i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
    return (h, n) if _ring_nn_dist(five, [h, n]) == 2 else None


def _di13_core_ok(
    mol: Mol, five: list[int], six: list[int], spec: Fused56Di13Spec,
) -> tuple[int, int] | None:
    atoms = set(five) | set(six)
    if len(atoms) != spec.n_atoms or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_zn(mol, five, spec.hetero_z)


def _di13_parts(
    info: dict, spec: Fused56Di13Spec,
) -> tuple[list[int], list[int], int, int, int, int] | None:
    """Return (five, six, hetero_idx, n_idx, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    hn = _di13_core_ok(info["mol"], five, six, spec)
    return None if hn is None else (five, six, hn[0], hn[1], ba, bb)


def _methyl_starts_ok(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    """Outside C atoms equal starts; each is methyl or CF3 on ring."""
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl_on_ring(mol, s, ring) for s in starts)


def _di13_subs_ok(mol: Mol, ring: set[int], cap: int) -> bool:
    """Allow ≤cap simple ring subs: halo + methyl + CF3 (multi-ok)."""
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > cap:
        return False
    return True if not starts else _methyl_starts_ok(mol, ring, starts)


def _is_simple_di13(info: dict, spec: Fused56Di13Spec) -> bool:
    if _di13_parts(info, spec) is None or _fg_block(info):
        return False
    parts = _di13_parts(info, spec)
    assert parts is not None
    mol, ring = info["mol"], _ring_set_parts(parts)
    if not _outside_ok(mol, ring):
        return False
    return _di13_subs_ok(mol, ring, spec.sub_cap)


def _di13_parent_dict(info: dict, spec: Fused56Di13Spec, kind: str, **extra) -> dict:
    parts = _di13_parts(info, spec)
    assert parts is not None
    five, six, h, n, ba, bb = parts
    chain = _chain_atoms(info["mol"], five, six, h, ba, bb) or []
    return {
        "chain": chain, "n_carbons": spec.n_atoms, "kind": kind,
        spec.hetero_key: h, "n_idx": n, "nh_idx": h, "bridge": [ba, bb], **extra,
    }


def _try_di13_fused56(info: dict, spec: Fused56Di13Spec) -> dict | None:
    """Simple 1,3-benzoazole parent (unsub / halo+methyl+CF3 ≤ sub_cap)."""
    if not _is_simple_di13(info, spec):
        return None
    return _di13_parent_dict(info, spec, spec.kind)


def _pos2_from_di13(mol: Mol, parts: tuple) -> int | None:
    """C2 = non-bridge five-ring neighbor of hetero (not N)."""
    five, _six, h, n, ba, bb = parts
    bridge = {ba, bb}
    nbs = [
        x.GetIdx() for x in mol.GetAtomWithIdx(h).GetNeighbors()
        if x.GetIdx() in set(five) and x.GetIdx() not in bridge
    ]
    if len(nbs) != 1:
        return None
    return nbs[0] if nbs[0] != n else None


def _primary_amine_on_c2(info: dict, ring: set[int], c2: int) -> dict | None:
    am = _mono_amine_on_ring(info, ring)
    if am is None or am.get("degree") != 1:
        return None
    return am if am.get("c_idx") == c2 else None


def _di13_amine_ctx(
    info: dict, spec: Fused56Di13Spec,
) -> tuple[Mol, set[int], dict] | None:
    """Core + C2 primary amine context, or None."""
    if _di13_parts(info, spec) is None or _fg_block(info, _AMINE_BLOCK_KEYS):
        return None
    parts = _di13_parts(info, spec)
    assert parts is not None
    mol, ring = info["mol"], _ring_set_parts(parts)
    c2 = _pos2_from_di13(mol, parts)
    am = _primary_amine_on_c2(info, ring, c2) if c2 is not None else None
    return (mol, ring, am) if am is not None else None


def _is_simple_di13_amine(info: dict, spec: Fused56Di13Spec) -> bool:
    ctx = _di13_amine_ctx(info, spec)
    if ctx is None:
        return False
    mol, ring, am = ctx
    if not _outside_ok(mol, ring, {am["n_idx"]}):
        return False
    return _di13_subs_ok(mol, ring, spec.amine_sub_cap)


def _try_di13_amine(info: dict, spec: Fused56Di13Spec) -> dict | None:
    """C2 primary amine 1,3-benzoazole parent (≤ amine_sub_cap halo/CF3)."""
    if not _is_simple_di13_amine(info, spec):
        return None
    return _di13_parent_dict(
        info, spec, spec.amine_kind,
        amine_c_idx=info["amines"][0]["c_idx"],
    )
