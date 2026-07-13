"""Retained monocyclic heteroarene parents (IUPAC P-22.2.1 / P-22.1).

Five-membered mono-hetero: furan / thiophene / 1H-pyrrole (optional monomethyl/monohalo).
Five-membered di-aza: 1H-imidazole (optional monomethyl/monohalo).
Six-membered diazines: pyrimidine / pyrazine / pyridazine (unsubstituted only).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import _outside_ok, _ring_halo_n, _ring_side_starts

# Z of ring hetero → parent kind (5-membered mono)
_HETERO5_KIND = {8: "furan", 16: "thiophene", 7: "pyrrole"}
# min ring distance between two N → diazine kind
_DIAZINE_KIND = {1: "pyridazine", 2: "pyrimidine", 3: "pyrazine"}


def _ring_atoms_if_mono(info: dict) -> list[int] | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    return list(rings[0]["atom_ids"])


def _hetero5_zs(mol: Mol, atom_ids: list[int]) -> list[int]:
    return [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids]


def _is_hetero5_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or len(atom_ids) != 5:
        return False
    mol: Mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):
        return False
    zs = _hetero5_zs(mol, atom_ids)
    return zs.count(6) == 4 and sum(1 for z in zs if z in _HETERO5_KIND) == 1


def _hetero5_idx(info: dict) -> int | None:
    if not _is_hetero5_core(info):
        return None
    mol: Mol = info["mol"]
    for i in info["rings"][0]["atom_ids"]:
        if mol.GetAtomWithIdx(i).GetAtomicNum() in _HETERO5_KIND:
            return i
    return None


def _hetero5_kind(info: dict) -> str | None:
    idx = _hetero5_idx(info)
    if idx is None:
        return None
    z = info["mol"].GetAtomWithIdx(idx).GetAtomicNum()
    if z == 7 and info["mol"].GetAtomWithIdx(idx).GetTotalNumHs() < 1:
        return None
    return _HETERO5_KIND.get(z)


def _hetero5_fg_block(info: dict) -> bool:
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_ester") or info.get("has_amide"):
        return True
    return bool(info.get("has_nitrile") or info.get("has_thiol") or info.get("has_nitro"))


def _mono_methyl_only(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    if len(starts) != 1:
        return False
    outside = [
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    ]
    return outside == starts


def _hetero5_subs_ok(mol: Mol, ring: set[int]) -> bool:
    h, starts = _ring_halo_n(mol, ring), _ring_side_starts(mol, ring)
    if h + len(starts) > 1:
        return False
    if not starts:
        return True
    return _mono_methyl_only(mol, ring, starts)


def _is_simple_hetero5(info: dict) -> bool:
    if _hetero5_kind(info) is None or _hetero5_fg_block(info):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring):
        return False
    return _hetero5_subs_ok(mol, ring)


def _hetero5_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    kind = _hetero5_kind(info)
    return {
        "chain": ring, "n_carbons": 5, "kind": kind,
        "hetero_idx": _hetero5_idx(info),
    }


def _try_hetero5_parent(info: dict) -> dict | None:
    return _hetero5_parent(info) if _is_simple_hetero5(info) else None


def _is_diazine_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or len(atom_ids) != 6:
        return False
    mol: Mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):
        return False
    zs = _hetero5_zs(mol, atom_ids)
    return zs.count(7) == 2 and zs.count(6) == 4


def _diazine_n_idxs(info: dict) -> list[int] | None:
    if not _is_diazine_core(info):
        return None
    mol: Mol = info["mol"]
    return [i for i in info["rings"][0]["atom_ids"] if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]


def _ring_nn_dist(atom_ids: list[int], n_idxs: list[int]) -> int:
    ia, ib = atom_ids.index(n_idxs[0]), atom_ids.index(n_idxs[1])
    d = abs(ia - ib)
    return min(d, len(atom_ids) - d)


def _diazine_kind(info: dict) -> str | None:
    n_idxs = _diazine_n_idxs(info)
    if n_idxs is None or len(n_idxs) != 2:
        return None
    ring = list(info["rings"][0]["atom_ids"])
    return _DIAZINE_KIND.get(_ring_nn_dist(ring, n_idxs))


def _diazine_unsub(mol: Mol, ring: set[int]) -> bool:
    return _ring_halo_n(mol, ring) == 0 and not _ring_side_starts(mol, ring)


def _is_simple_diazine(info: dict) -> bool:
    if _diazine_kind(info) is None or _hetero5_fg_block(info):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring):
        return False
    return _diazine_unsub(mol, ring)


def _diazine_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    return {
        "chain": ring, "n_carbons": 6, "kind": _diazine_kind(info),
        "n_idxs": _diazine_n_idxs(info),
    }


def _try_diazine_parent(info: dict) -> dict | None:
    return _diazine_parent(info) if _is_simple_diazine(info) else None


def _is_imidazole_core(info: dict) -> bool:
    atom_ids = _ring_atoms_if_mono(info)
    if atom_ids is None or len(atom_ids) != 5:
        return False
    mol: Mol = info["mol"]
    if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):
        return False
    zs = _hetero5_zs(mol, atom_ids)
    return zs.count(7) == 2 and zs.count(6) == 3


def _ring_n_idxs(mol: Mol, ring: list[int]) -> list[int]:
    return [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]


def _nh_among(mol: Mol, ns: list[int]) -> int | None:
    hs = [i for i in ns if mol.GetAtomWithIdx(i).GetTotalNumHs() >= 1]
    return hs[0] if len(hs) == 1 else None


def _imidazole_n_pair(info: dict) -> tuple[int, int] | None:
    """Return (nh_idx, n_idx) for 1H-imidazole (1,3-diazole); else None."""
    if not _is_imidazole_core(info):
        return None
    mol, ring = info["mol"], list(info["rings"][0]["atom_ids"])
    ns = _ring_n_idxs(mol, ring)
    if len(ns) != 2 or _ring_nn_dist(ring, ns) != 2:
        return None
    nh = _nh_among(mol, ns)
    return None if nh is None else (nh, ns[0] if ns[1] == nh else ns[1])


def _is_simple_imidazole(info: dict) -> bool:
    if _imidazole_n_pair(info) is None or _hetero5_fg_block(info):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring):
        return False
    return _hetero5_subs_ok(mol, ring)


def _imidazole_parent(info: dict) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    nh, n = _imidazole_n_pair(info)
    return {
        "chain": ring, "n_carbons": 5, "kind": "imidazole",
        "nh_idx": nh, "n_idx": n, "n_idxs": [nh, n],
    }


def _try_imidazole_parent(info: dict) -> dict | None:
    return _imidazole_parent(info) if _is_simple_imidazole(info) else None
