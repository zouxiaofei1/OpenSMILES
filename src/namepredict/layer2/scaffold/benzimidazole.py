"""Retained 1H-benzimidazole / benzimidazolamine parent (IUPAC P-22.2.1 / P-25).

Fused aromatic 6+5: benzo[d]imidazole. NH=1, N=3 (1,3 on five-ring);
≤2 halo / methyl / CF3; mono primary amine at C2 → benzimidazolamine
(with ≤1 halo or CF3). Also accepts 2-imino dual-NH tautomer of 2-amine.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.fused56 import (
    _all_aromatic,
    _chain_atoms,
    _fused56,
    _six_all_c,
)
from namepredict.layer2.scaffold.heteroarene5 import _ring_nn_dist
from namepredict.layer3.aryl_sub import _aryl_atoms, _aryl_exclude, _aryl_sub_n
from namepredict.layer2.scaffold.ring_parent import (
    _is_methyl_on_ring,
    _mono_amine_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)


def _five_n_idxs(mol: Mol, five: list[int]) -> list[int] | None:
    """Exactly two N and three C on five-ring; return N indices."""
    zs = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in five]
    if zs.count(7) != 2 or zs.count(6) != 3:
        return None
    return [i for i in five if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]


def _single_nh_pair(mol: Mol, five: list[int], ns: list[int]) -> tuple[int, int] | None:
    """1,3 N–N with exactly one NH; return (nh, n)."""
    if _ring_nn_dist(five, ns) != 2:
        return None
    nhs = [i for i in ns if mol.GetAtomWithIdx(i).GetTotalNumHs() >= 1]
    if len(nhs) != 1:
        return None
    nh = nhs[0]
    return nh, (ns[0] if ns[1] == nh else ns[1])


def _five_nn(mol: Mol, five: list[int]) -> tuple[int, int] | None:
    """Exactly two N, three C; N–N ring dist 2 (1,3). Return (nh, n) or None."""
    ns = _five_n_idxs(mol, five)
    return None if ns is None else _single_nh_pair(mol, five, ns)


def _core_ok(mol: Mol, five: list[int], six: list[int]) -> tuple[int, int] | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_nn(mol, five)


def _bim_parts(info: dict) -> tuple[list[int], list[int], int, int, int, int] | None:
    """Return (five, six, nh, n, ba, bb) or None."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    nn = _core_ok(info["mol"], five, six)
    return None if nn is None else (five, six, nn[0], nn[1], ba, bb)


def _is_bim_core(info: dict) -> bool:
    return _bim_parts(info) is not None


def _ring_set(parts: tuple) -> set[int]:
    return set(parts[0]) | set(parts[1])


def _bim_fg_block(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_amine",
        "has_thiol", "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _starts_ok(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl_on_ring(mol, s, ring) for s in starts)


def _bim_subs_ok(info: dict, mol: Mol, ring: set[int], cap: int = 2) -> bool:
    """Allow <=cap ring subs: halo + methyl/CF3 + aryl."""
    n_aryl = _aryl_sub_n(info, ring)
    excl = _aryl_exclude(info, ring)
    h = _ring_halo_n(mol, ring)
    starts = _ring_side_starts(mol, ring, excl)
    if h + len(starts) + n_aryl > cap:
        return False
    return True if not starts else _starts_ok(mol, ring, starts, excl)


def _starts_ok(mol: Mol, ring: set[int], starts: list[int], skip: set[int] | None = None) -> bool:
    skip = skip or set()
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring and a.GetIdx() not in skip
    ]
    if set(outside) != set(starts):
        return False
    return all(_is_methyl_on_ring(mol, s, ring) for s in starts)


def _is_simple_benzimidazole(info: dict) -> bool:
    if not _is_bim_core(info) or _bim_fg_block(info):
        return False
    mol: Mol = info["mol"]
    parts = _bim_parts(info)
    assert parts is not None
    ring = _ring_set(parts)
    allowed = _aryl_atoms(info, ring)
    return _outside_ok(mol, ring, allowed) and _bim_subs_ok(info, mol, ring)


def _build_chain(
    mol: Mol, five: list[int], six: list[int], nh: int, ba: int, bb: int,
) -> list[int] | None:
    """IUPAC order: 1=NH, 2, 3=N, 3a, 4, 5, 6, 7, 7a."""
    return _chain_atoms(mol, five, six, nh, ba, bb)


def _bim_parent_dict(info: dict, kind: str, **extra) -> dict:
    parts = _bim_parts(info)
    assert parts is not None
    five, six, nh, n, ba, bb = parts
    chain = _build_chain(info["mol"], five, six, nh, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": kind,
        "scaffold_id": kind, "nh_idx": nh, "n_idx": n, "bridge": [ba, bb], **extra,
    }


def _benzimidazole_parent(info: dict) -> dict:
    return _bim_parent_dict(info, "benzimidazole")


def _try_benzimidazole_parent(info: dict) -> dict | None:
    return _benzimidazole_parent(info) if _is_simple_benzimidazole(info) else None


def _bim_amine_conflict(info: dict) -> bool:
    keys = (
        "has_acid", "has_aldehyde", "has_ketone", "has_alcohol",
        "has_ester", "has_amide", "has_nitrile", "has_thiol",
        "has_nitro", "has_acyl_chloride", "has_anhydride",
    )
    return any(info.get(k) for k in keys)


def _pos2_from_parts(mol: Mol, parts: tuple) -> int | None:
    five, _six, nh, n, ba, bb = parts
    bridge = {ba, bb}
    nbs = [
        x.GetIdx() for x in mol.GetAtomWithIdx(nh).GetNeighbors()
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


def _bim_amine_ctx(info: dict) -> tuple[Mol, set[int], dict] | None:
    """Core + C2 primary amine context, or None."""
    if not _is_bim_core(info) or _bim_amine_conflict(info):
        return None
    parts = _bim_parts(info)
    assert parts is not None
    mol, ring = info["mol"], _ring_set(parts)
    c2 = _pos2_from_parts(mol, parts)
    am = _primary_amine_on_c2(info, ring, c2) if c2 is not None else None
    return (mol, ring, am) if am is not None else None


def _is_simple_benzimidazolamine(info: dict) -> bool:
    ctx = _bim_amine_ctx(info)
    if ctx is None:
        return False
    mol, ring, am = ctx
    if not _outside_ok(mol, ring, {am["n_idx"]}):
        return False
    return _bim_subs_ok(info, mol, ring, cap=1)


def _benzimidazolamine_parent(info: dict) -> dict:
    return _bim_parent_dict(
        info, "benzimidazolamine", amine_c_idx=info["amines"][0]["c_idx"],
    )


def _try_benzimidazolamine_parent(info: dict) -> dict | None:
    if _is_simple_benzimidazolamine(info):
        return _benzimidazolamine_parent(info)
    return _try_bim_imino_parent(info)


# --- 2-imino dual-NH tautomer of 1H-benzimidazol-2-amine ---


def _both_nh(mol: Mol, five: list[int], ns: list[int]) -> list[int] | None:
    """Both N carry H and N–N ring dist 2; return ns."""
    if _ring_nn_dist(five, ns) != 2:
        return None
    if any(mol.GetAtomWithIdx(i).GetTotalNumHs() < 1 for i in ns):
        return None
    return ns


def _five_dual_nh(mol: Mol, five: list[int]) -> list[int] | None:
    """Two N + three C; both N carry H; N–N ring dist 2. Return N indices."""
    ns = _five_n_idxs(mol, five)
    return None if ns is None else _both_nh(mol, five, ns)


def _imino_core_ok(mol: Mol, five: list[int], six: list[int]) -> list[int] | None:
    atoms = set(five) | set(six)
    if len(atoms) != 9 or not _all_aromatic(mol, atoms):
        return None
    if not _six_all_c(mol, six):
        return None
    return _five_dual_nh(mol, five)


def _imino_parts(
    info: dict,
) -> tuple[list[int], list[int], list[int], int, int] | None:
    """Return (five, six, ns, ba, bb) for dual-NH benzimidazole core."""
    fused = _fused56(info)
    if fused is None:
        return None
    five, six, (ba, bb) = fused
    ns = _imino_core_ok(info["mol"], five, six)
    return None if ns is None else (five, six, ns, ba, bb)


def _exocyclic_nh_on(mol: Mol, c_idx: int, ring: set[int]) -> int | None:
    """Exocyclic N (=NH or -NH) on ring carbon; not aromatic ring N."""
    hits = [
        n.GetIdx() for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 7 and n.GetIdx() not in ring
        and not n.GetIsAromatic()
    ]
    return hits[0] if len(hits) == 1 else None


def _c2_between_ns(mol: Mol, five: list[int], ns: list[int]) -> int | None:
    """Carbon on five-ring bonded to both N (position 2)."""
    nset = set(ns)
    for i in five:
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 6:
            continue
        nbs = {n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()}
        if nset <= nbs:
            return i
    return None


def _imino_pick_nh(
    mol: Mol, five: list[int], six: list[int], ns: list[int], ba: int, bb: int, c2: int,
) -> int | None:
    """Prefer NH that yields a valid chain with c2 at position 2."""
    for nh in ns:
        chain = _chain_atoms(mol, five, six, nh, ba, bb)
        if chain is not None and len(chain) == 9 and chain[1] == c2:
            return nh
    return None


def _imino_c2_exo(
    mol: Mol, five: list[int], ns: list[int], ring: set[int],
) -> tuple[int, int] | None:
    c2 = _c2_between_ns(mol, five, ns)
    if c2 is None:
        return None
    exo = _exocyclic_nh_on(mol, c2, ring)
    return None if exo is None else (c2, exo)


def _imino_from_parts(
    mol: Mol, five: list[int], six: list[int], ns: list[int], ba: int, bb: int,
) -> tuple[Mol, set[int], int, int, int] | None:
    ring = set(five) | set(six)
    ce = _imino_c2_exo(mol, five, ns, ring)
    if ce is None:
        return None
    c2, exo = ce
    nh = _imino_pick_nh(mol, five, six, ns, ba, bb, c2)
    return (mol, ring, c2, exo, nh) if nh is not None else None


def _imino_ctx(info: dict) -> tuple[Mol, set[int], int, int, int] | None:
    """(mol, ring, c2, exo_n, nh) for 2-imino dual-NH, or None."""
    if _bim_amine_conflict(info) or info.get("has_amine"):
        return None
    parts = _imino_parts(info)
    if parts is None:
        return None
    five, six, ns, ba, bb = parts
    return _imino_from_parts(info["mol"], five, six, ns, ba, bb)


def _is_simple_bim_imino(info: dict) -> bool:
    ctx = _imino_ctx(info)
    if ctx is None:
        return False
    mol, ring, _c2, exo, _nh = ctx
    if not _outside_ok(mol, ring, {exo}):
        return False
    return _bim_subs_ok(info, mol, ring, cap=1)


def _bim_imino_dict(
    mol: Mol, five: list[int], six: list[int], ns: list[int],
    ba: int, bb: int, nh: int, c2: int,
) -> dict:
    n_other = ns[0] if ns[1] == nh else ns[1]
    chain = _build_chain(mol, five, six, nh, ba, bb) or []
    return {
        "chain": chain, "n_carbons": 9, "kind": "benzimidazolamine",
        "scaffold_id": "benzimidazolamine",
        "nh_idx": nh, "n_idx": n_other, "bridge": [ba, bb],
        "amine_c_idx": c2,
    }


def _bim_imino_parent(info: dict) -> dict:
    parts, ctx = _imino_parts(info), _imino_ctx(info)
    assert parts is not None and ctx is not None
    five, six, ns, ba, bb = parts
    mol, _ring, c2, _exo, nh = ctx
    return _bim_imino_dict(mol, five, six, ns, ba, bb, nh, c2)


def _try_bim_imino_parent(info: dict) -> dict | None:
    return _bim_imino_parent(info) if _is_simple_bim_imino(info) else None
