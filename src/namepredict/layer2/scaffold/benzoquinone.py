"""1,4-Benzoquinone parent: cyclohexa-2,5-diene-1,4-dione (IUPAC P-64.2).

Single C6 carbocycle + exactly two para ring ketones + two endocyclic C=C.
PIN is systematic (retained 1,4-benzoquinone is general-only).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import (
    _dbl_o_idx,
    _endocyclic_doubles,
    _is_carbocycle_ring,
    _no_hetero_outside,
    _outside_carbons,
    _ring_alkoxy_ethers,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _disjoint_cover, _walk_linear_n


def _ketone_idxs(info: dict) -> list[int]:
    return [e["c_idx"] for e in (info.get("ketones") or []) if "c_idx" in e]


def _bq_ring(info: dict) -> list[int] | None:
    atom_ids = _is_carbocycle_ring(info)
    return list(atom_ids) if atom_ids is not None and len(atom_ids) == 6 else None


def _ring_dist(ring: list[int], a: int, b: int) -> int:
    i, j = ring.index(a), ring.index(b)
    d = abs(i - j)
    return min(d, len(ring) - d)


def _ketones_para(ring: list[int], ket: list[int]) -> bool:
    if len(ket) != 2 or any(k not in ring for k in ket):
        return False
    return _ring_dist(ring, ket[0], ket[1]) == 3


def _bq_doubles_ok(info: dict, ring: list[int]) -> bool:
    return len(_endocyclic_doubles(info, set(ring))) == 2


def _ketone_o_set(mol: Mol, ket: list[int]) -> set[int] | None:
    out: set[int] = set()
    for c in ket:
        o = _dbl_o_idx(mol, c)
        if o is None:
            return None
        out.add(o)
    return out


def _oh_o_set(info: dict, ring: set[int]) -> set[int]:
    return {
        h["o_idx"] for h in (info.get("hydroxyls") or [])
        if h.get("c_idx") in ring and "o_idx" in h
    }


def _alkoxy_atom_set(info: dict, ring: set[int]) -> set[int]:
    out: set[int] = set()
    for a in _ring_alkoxy_ethers(info, ring):
        out.add(a["o_idx"])
        out.update(a["atoms"])
    return out


def _bq_exclude(info: dict, mol: Mol, ring: set[int], ket: list[int]) -> set[int] | None:
    oset = _ketone_o_set(mol, ket)
    if oset is None:
        return None
    return oset | _oh_o_set(info, ring) | _alkoxy_atom_set(info, ring)


def _bq_alkoxy_ok(info: dict, ring: set[int]) -> bool:
    return all(a.get("n") in (1, 2) for a in _ring_alkoxy_ethers(info, ring))


def _bq_side_paths(
    mol: Mol, ring: set[int], exclude: set[int],
) -> list[set[int]] | None:
    """Linear n-alkyl C1–C12 paths from each ring side start, or None."""
    blocked = ring | exclude
    paths: list[set[int]] = []
    for s in _ring_side_starts(mol, ring, exclude):
        path = _walk_linear_n(mol, s, blocked, max_n=12)
        if path is None:
            return None
        paths.append(set(path))
    return paths


def _bq_alkyls_ok(mol: Mol, ring: set[int], exclude: set[int]) -> bool:
    paths = _bq_side_paths(mol, ring, exclude)
    if paths is None:
        return False
    outside = set(_outside_carbons(mol, ring, exclude))
    return _disjoint_cover(paths, outside)


def _bq_subs_ok(info: dict, mol: Mol, ring: list[int], ket: list[int]) -> bool:
    rset = set(ring)
    exclude = _bq_exclude(info, mol, rset, ket)
    if exclude is None or not _no_hetero_outside(mol, rset, exclude):
        return False
    return _bq_alkoxy_ok(info, rset) and _bq_alkyls_ok(mol, rset, exclude)


_BQ_BLOCK = (
    "has_acid", "has_aldehyde", "has_ester", "has_amide",
    "has_nitrile", "has_amine", "has_thiol", "has_nitro",
    "has_acyl_chloride", "has_anhydride",
)


def _bq_fg_block(info: dict) -> bool:
    return any(info.get(k) for k in _BQ_BLOCK)


def _is_simple_benzoquinone(info: dict) -> bool:
    ring, ket = _bq_ring(info), _ketone_idxs(info)
    if ring is None or _bq_fg_block(info) or not _ketones_para(ring, ket):
        return False
    if not _bq_doubles_ok(info, ring):
        return False
    return _bq_subs_ok(info, info["mol"], ring, ket)


def _benzoquinone_parent(info: dict) -> dict:
    ring = _bq_ring(info) or []
    ket = _ketone_idxs(info)
    bonds = _endocyclic_doubles(info, set(ring))
    return {
        "chain": ring, "n_carbons": 6, "kind": "benzoquinone",
        "scaffold_id": "benzoquinone", "ketone_c_idxs": ket,
        "double_bonds": bonds,
    }


def _try_benzoquinone_parent(info: dict) -> dict | None:
    return _benzoquinone_parent(info) if _is_simple_benzoquinone(info) else None


# ── 1,2-benzoquinone (ortho) — same substitution gates, ring distance only ──

def _ketones_ortho(ring: list[int], ket: list[int]) -> bool:
    if len(ket) != 2 or any(k not in ring for k in ket):
        return False
    return _ring_dist(ring, ket[0], ket[1]) == 1


def _is_simple_ortho_benzoquinone(info: dict) -> bool:
    ring, ket = _bq_ring(info), _ketone_idxs(info)
    if ring is None or _bq_fg_block(info) or not _ketones_ortho(ring, ket):
        return False
    if not _bq_doubles_ok(info, ring):
        return False
    return _bq_subs_ok(info, info["mol"], ring, ket)


def _ortho_benzoquinone_parent(info: dict) -> dict:
    ring = _bq_ring(info) or []
    ket = _ketone_idxs(info)
    bonds = _endocyclic_doubles(info, set(ring))
    return {
        "chain": ring, "n_carbons": 6, "kind": "ortho_benzoquinone",
        "scaffold_id": "ortho_benzoquinone", "ketone_c_idxs": ket,
        "double_bonds": bonds,
    }


def _try_ortho_benzoquinone_parent(info: dict) -> dict | None:
    if not _is_simple_ortho_benzoquinone(info):
        return None
    return _ortho_benzoquinone_parent(info)
