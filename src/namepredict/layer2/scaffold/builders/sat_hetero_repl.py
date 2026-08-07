"""Saturated monohetero replacement parent (table-out → x杂环y烷).

Retained `_MONO` (oxolane, oxane, piperidine, …) is handled by sat_hetero.
This builder only accepts cores that miss the retained table.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import _ring_bonds_single
from namepredict.layer2.scaffold.sat_hetero import (
    _all_nonarom,
    _ring_atoms_if_mono,
    _ring_zs,
)
from namepredict.layer2.scaffold.hetero_a_names import A_PREFIX_EN

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


def _rotations(chain: list[int]) -> list[list[int]]:
    n = len(chain)
    out: list[list[int]] = []
    for i in range(n):
        rot = chain[i:] + chain[:i]
        out.append(rot)
        out.append(list(reversed(rot)))
    return out


def _parent_dict(
    order: list[int], hs: list[int], zmap: dict, sites: list, en: str, zh: str,
) -> dict:
    return {
        "kind": "sat_hetero_repl", "chain": order, "n_carbons": len(order),
        "stem_en": en, "stem_zh": zh, "hetero_sites": tuple(sites),
        "ring_size": len(order), "hetero_idxs": hs, "hetero_z": zmap,
    }


