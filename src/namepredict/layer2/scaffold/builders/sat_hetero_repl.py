"""Saturated monohetero replacement parent (table-out → x杂环y烷).

Retained `_MONO` (oxolane, oxane, piperidine, …) is handled by sat_hetero.
This builder only accepts cores that miss the retained table.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import _outside_ok, _ring_bonds_single
from namepredict.layer2.scaffold.sat_hetero import (
    _all_nonarom,
    _extra_amine,
    _fg_block,
    _kind_of,
    _n_unsub,
    _ring_atoms_if_mono,
    _ring_hetero_idxs,
    _ring_n_set,
    _ring_zs,
    _subs_ok,
)
from namepredict.layer2.scaffold.hetero_a_names import A_PREFIX_EN
from namepredict.layer2.scaffold.sat_hetero_stem import (
    sat_hetero_stem_en,
    sat_hetero_stem_zh,
)

# IUPAC a-order rank (lower = better low locant); mirror L4 multi_hetero.
_Z_RANK: dict[int, int] = {8: 0, 16: 1, 34: 2, 52: 3, 7: 4, 15: 5, 5: 6, 14: 7}
_ALLOWED_Z = frozenset(A_PREFIX_EN) | {6}


def _only_allowed(zs: list[int]) -> bool:
    return all(z in _ALLOWED_Z for z in zs) and any(z != 6 for z in zs)


def _size_ok(n: int) -> bool:
    return 3 <= n <= 12


def _is_repl_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or not _size_ok(len(atom_ids)):
        return False
    mol: Mol = info["mol"]
    if not _all_nonarom(mol, atom_ids) or not _ring_bonds_single(mol, tuple(atom_ids)):
        return False
    return _only_allowed(_ring_zs(mol, atom_ids))


def _is_table_out(info: dict) -> bool:
    return _is_repl_core(info) and _kind_of(info) is None


def _simple_ok(info: dict) -> bool:
    if not _is_table_out(info) or _fg_block(info):
        return False
    mol: Mol = info["mol"]
    atom_ids = _ring_atoms_if_mono(info)
    ring, ring_ns = set(atom_ids), _ring_n_set(mol, atom_ids)
    if _extra_amine(info, ring_ns) or not _n_unsub(mol, ring_ns):
        return False
    return _outside_ok(mol, ring) and _subs_ok(mol, ring)


def _zmap(mol: Mol, hs: list[int]) -> dict[int, int]:
    return {i: mol.GetAtomWithIdx(i).GetAtomicNum() for i in hs}


def _rotations(chain: list[int]) -> list[list[int]]:
    n = len(chain)
    out: list[list[int]] = []
    for i in range(n):
        rot = chain[i:] + chain[:i]
        out.append(rot)
        out.append(list(reversed(rot)))
    return out


def _hetero_key(order: list[int], hs: list[int], zmap: dict[int, int]) -> tuple:
    pairs = sorted((order.index(h) + 1, _Z_RANK.get(zmap[h], 99)) for h in hs)
    locs = tuple(loc for loc, _ in pairs)
    ranks = tuple(r for _, r in pairs)
    return (locs, ranks, tuple(order))


def _best_order(chain: list[int], hs: list[int], zmap: dict[int, int]) -> list[int]:
    cands = _rotations(chain)
    return min(cands, key=lambda o: _hetero_key(o, hs, zmap))


def _sites(order: list[int], hs: list[int], zmap: dict[int, int]) -> list[tuple[int, int]]:
    return sorted((order.index(h) + 1, zmap[h]) for h in hs)


def _stems(size: int, sites: list[tuple[int, int]]) -> tuple[str, str] | None:
    en, zh = sat_hetero_stem_en(size, sites), sat_hetero_stem_zh(size, sites)
    if en is None or zh is None:
        return None
    return en, zh


def _orient_sites(info: dict) -> tuple[list[int], list[int], dict, list] | None:
    mol: Mol = info["mol"]
    atom_ids = list(_ring_atoms_if_mono(info))
    hs = _ring_hetero_idxs(mol, atom_ids)
    zmap = _zmap(mol, hs)
    order = _best_order(atom_ids, hs, zmap)
    return order, hs, zmap, _sites(order, hs, zmap)


def _parent_dict(
    order: list[int], hs: list[int], zmap: dict, sites: list, en: str, zh: str,
) -> dict:
    return {
        "kind": "sat_hetero_repl", "chain": order, "n_carbons": len(order),
        "stem_en": en, "stem_zh": zh, "hetero_sites": tuple(sites),
        "ring_size": len(order), "hetero_idxs": hs, "hetero_z": zmap,
    }


def _build_parent(info: dict) -> dict | None:
    got = _orient_sites(info)
    if got is None:
        return None
    order, hs, zmap, sites = got
    stems = _stems(len(order), sites)
    return None if stems is None else _parent_dict(order, hs, zmap, sites, *stems)


def try_sat_hetero_repl(info: dict) -> dict | None:
    """Table-out saturated monohetero → replacement parent; retained miss only."""
    return _build_parent(info) if _simple_ok(info) else None
