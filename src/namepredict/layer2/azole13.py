"""1,3-oxazole / 1,3-thiazole retained parents (IUPAC P-22.2.1 / P-14.3.4).

Monocyclic fully aromatic C3N+O or C3N+S rings with N–O/S ring distance 2
(standard 1,3). O/S fixed as locant 1 (stored as nh_idx for orient reuse),
N as 3. Simple ring subs: halo + straight n-alkyl C1–C3, total ≤2.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.heteroarene5 import (
    _hetero5_fg_block,
    _hetero5_zs,
    _ring_atoms_if_mono,
    _ring_nn_dist,
)
from namepredict.layer2.ring_parent import (
    _outside_ok,
    _ring_halo_n,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _linear_or_omega_halo_sides_ok

# O/S atomic number → kind
_AZOLE13_KIND = {8: "oxazole", 16: "thiazole"}


def _azole13_hetero_pair(info: dict) -> tuple[int, int] | None:
    """Return (hetero_idx O/S, n_idx) for 1,3-oxazole/thiazole core; else None."""
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or len(atom_ids) != 5:
        return None
    mol: Mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):
        return None
    if not _azole13_z_ok(mol, atom_ids):
        return None
    return _azole13_pair_from_ring(mol, atom_ids)


def _azole13_z_ok(mol: Mol, atom_ids: list[int]) -> bool:
    zs = _hetero5_zs(mol, atom_ids)
    if zs.count(6) != 3 or zs.count(7) != 1:
        return False
    return sum(1 for z in zs if z in _AZOLE13_KIND) == 1


def _azole13_pair_from_ring(mol: Mol, ring: list[int]) -> tuple[int, int] | None:
    h = next(i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() in _AZOLE13_KIND)
    n = next(i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
    if _ring_nn_dist(ring, [h, n]) != 2:
        return None
    return h, n


def _azole13_kind(info: dict) -> str | None:
    pair = _azole13_hetero_pair(info)
    if pair is None:
        return None
    z = info["mol"].GetAtomWithIdx(pair[0]).GetAtomicNum()
    return _AZOLE13_KIND.get(z)


def _azole13_subs_ok(mol: Mol, ring: set[int]) -> bool:
    """Allow ≤2 simple ring subs: monohalo and/or n-alkyl C1–C3."""
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > 2:
        return False
    return _linear_or_omega_halo_sides_ok(mol, ring, starts, 3)


def _is_simple_azole13(info: dict) -> bool:
    if _azole13_kind(info) is None or _hetero5_fg_block(info):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring):
        return False
    return _azole13_subs_ok(mol, ring)


def _azole13_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    h, n = _azole13_hetero_pair(info)
    return {
        "chain": ring, "n_carbons": 5, "kind": _azole13_kind(info),
        "nh_idx": h, "n_idx": n, "hetero_idx": h, "n_idxs": [h, n],
    }


def _try_azole13_parent(info: dict) -> dict | None:
    return _azole13_parent(info) if _is_simple_azole13(info) else None
