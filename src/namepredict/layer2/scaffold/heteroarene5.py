"""Retained monocyclic heteroarene parents (IUPAC P-22.2.1 / P-22.1).

Five-membered mono-hetero: furan / thiophene / 1H-pyrrole (optional monomethyl/monohalo).
Five-membered di-aza: 1H-imidazole / 1H-pyrazole (optional monomethyl/monohalo).
Six-membered diazines: pyrimidine / pyrazine / pyridazine (simple ring subs / pyrimidinamine).
"""
from __future__ import annotations

from rdkit.Chem import Mol

# Z of ring hetero → parent kind (5-membered mono)
_HETERO5_KIND = {8: "furan", 16: "thiophene", 7: "pyrrole"}
# min ring distance between two N → diazine kind
_DIAZINE_KIND = {1: "pyridazine", 2: "pyrimidine", 3: "pyrazine"}

def _ring_atoms_if_mono(info: dict) -> list[int] | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    return list(rings[0]["atom_ids"])

def _ring_nn_dist(atom_ids: list[int], n_idxs: list[int]) -> int:
    ia, ib = atom_ids.index(n_idxs[0]), atom_ids.index(n_idxs[1])
    d = abs(ia - ib)
    return min(d, len(atom_ids) - d)

def _ring_n_idxs(mol: Mol, ring: list[int]) -> list[int]:
    return [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]

