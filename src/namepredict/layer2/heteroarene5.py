"""Retained monocyclic heteroarene parents (IUPAC P-22.2.1 / P-22.1).

Five-membered mono-hetero: furan / thiophene / 1H-pyrrole (optional monomethyl/monohalo).
Five-membered di-aza: 1H-imidazole / 1H-pyrazole (optional monomethyl/monohalo).
Six-membered diazines: pyrimidine / pyrazine / pyridazine (simple ring subs / pyrimidinamine).
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import (
    _arene_alkoxy,
    _mono_amine_on_ring,
    _outside_ok,
    _ring_halo_n,
    _ring_nitro_atoms,
    _ring_nitro_n,
    _ring_side_starts,
)
from namepredict.layer2.side_alkyl import _linear_or_omega_halo_sides_ok
from namepredict.layer2.aryl_sub import _aryl_atoms, _aryl_exclude, _aryl_sub_n

# Z of ring hetero → parent kind (5-membered mono)
_HETERO5_KIND = {8: "furan", 16: "thiophene", 7: "pyrrole"}
# min ring distance between two N → diazine kind
_DIAZINE_KIND = {1: "pyridazine", 2: "pyrimidine", 3: "pyrazine"}


def _ring_atoms_if_mono(info: dict) -> list[int] | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    return list(rings[0]["atom_ids"])


def _unfused_ring(mol: Mol, ring: set[int]) -> bool:
    for i in ring:
        if sum(1 for r in mol.GetRingInfo().AtomRings() if i in r) != 1:
            return False
    return True


def _arom_ring_lists(mol: Mol, n: int, pred) -> list[list[int]]:
    out: list[list[int]] = []
    for r in mol.GetRingInfo().AtomRings():
        ids = list(r)
        if len(ids) != n or not _unfused_ring(mol, set(ids)):
            continue
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ids):
            continue
        if pred(mol, ids):
            out.append(ids)
    return out


def _is_c4n2(mol: Mol, atoms: list[int]) -> bool:
    zs = _hetero5_zs(mol, atoms)
    return zs.count(7) == 2 and zs.count(6) == 4


def _is_c4x_mono(mol: Mol, atoms: list[int]) -> bool:
    zs = _hetero5_zs(mol, atoms)
    return zs.count(6) == 4 and sum(1 for z in zs if z in _HETERO5_KIND) == 1


def _diazine_rings(mol: Mol) -> list[list[int]]:
    return _arom_ring_lists(mol, 6, _is_c4n2)


def _hetero5_rings(mol: Mol) -> list[list[int]]:
    return _arom_ring_lists(mol, 5, _is_c4x_mono)


def _hetero5_zs(mol: Mol, atom_ids: list[int]) -> list[int]:
    return [mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids]


def _is_hetero5_core(info: dict) -> bool:
    return bool(_hetero5_rings(info["mol"]))


def _hetero5_idx_in(mol: Mol, ring: list[int]) -> int | None:
    for i in ring:
        if mol.GetAtomWithIdx(i).GetAtomicNum() in _HETERO5_KIND:
            return i
    return None


def _hetero5_kind_of(mol: Mol, ring: list[int]) -> str | None:
    idx = _hetero5_idx_in(mol, ring)
    if idx is None:
        return None
    z = mol.GetAtomWithIdx(idx).GetAtomicNum()
    if z == 7 and mol.GetAtomWithIdx(idx).GetTotalNumHs() < 1:
        return None
    return _HETERO5_KIND.get(z)


def _hetero5_idx(info: dict) -> int | None:
    rings = _hetero5_rings(info["mol"])
    return None if not rings else _hetero5_idx_in(info["mol"], rings[0])


def _hetero5_kind(info: dict) -> str | None:
    rings = _hetero5_rings(info["mol"])
    return None if not rings else _hetero5_kind_of(info["mol"], rings[0])


def _hetero5_fg_block(info: dict) -> bool:
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_ester") or info.get("has_amide"):
        return True
    return bool(info.get("has_nitrile") or info.get("has_thiol") or info.get("has_nitro"))


def _mono_methyl_only(mol: Mol, ring: set[int], starts: list[int]) -> bool:
    """True when every ring side start is a pure -CH3 (outside C set == starts)."""
    if not starts:
        return True
    outside = {
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring
    }
    return outside == set(starts)


def _hetero5_subs_ok(info: dict, mol: Mol, ring: set[int]) -> bool:
    """Allow ≤2 ring subs: halo / methyl / aryl (P-14.3.4 / P-29.3)."""
    n_aryl = _aryl_sub_n(info, ring)
    excl = _aryl_exclude(info, ring)
    h = _ring_halo_n(mol, ring)
    starts = _ring_side_starts(mol, ring, excl)
    if h + len(starts) + n_aryl > 2:
        return False
    return True if not starts else _linear_or_omega_halo_sides_ok(mol, ring, starts, 2)


def _h5_ring_ok(info: dict, ring: list[int]) -> bool:
    mol, rs = info["mol"], set(ring)
    if not _outside_ok(mol, rs, _aryl_atoms(info, rs)):
        return False
    return _hetero5_subs_ok(info, mol, rs)


def _pick_hetero5_ring(info: dict) -> list[int] | None:
    cands = [r for r in _hetero5_rings(info["mol"]) if _h5_ring_ok(info, r)]
    return min(cands, key=min) if cands else None


def _is_simple_hetero5(info: dict) -> bool:
    if _hetero5_fg_block(info):
        return False
    return _pick_hetero5_ring(info) is not None


def _hetero5_parent(info: dict) -> dict:
    ring = _pick_hetero5_ring(info) or _hetero5_rings(info["mol"])[0]
    mol = info["mol"]
    return {
        "chain": ring, "n_carbons": 5, "kind": _hetero5_kind_of(mol, ring),
        "hetero_idx": _hetero5_idx_in(mol, ring),
    }


def _try_hetero5_parent(info: dict) -> dict | None:
    return _hetero5_parent(info) if _is_simple_hetero5(info) else None


def _is_diazine_core(info: dict) -> bool:
    return bool(_diazine_rings(info["mol"]))


def _n_idxs_in(mol: Mol, ring: list[int]) -> list[int]:
    return [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7]


def _diazine_n_idxs(info: dict, ring: list[int] | None = None) -> list[int] | None:
    rings = _diazine_rings(info["mol"])
    if ring is None:
        ring = rings[0] if rings else None
    if ring is None:
        return None
    ns = _n_idxs_in(info["mol"], ring)
    return ns if len(ns) == 2 else None


def _ring_nn_dist(atom_ids: list[int], n_idxs: list[int]) -> int:
    ia, ib = atom_ids.index(n_idxs[0]), atom_ids.index(n_idxs[1])
    d = abs(ia - ib)
    return min(d, len(atom_ids) - d)


def _diazine_kind_of(ring: list[int], n_idxs: list[int]) -> str | None:
    return _DIAZINE_KIND.get(_ring_nn_dist(ring, n_idxs))


def _diazine_kind(info: dict, ring: list[int] | None = None) -> str | None:
    rings = _diazine_rings(info["mol"])
    ring = ring or (rings[0] if rings else None)
    ns = _diazine_n_idxs(info, ring)
    if ring is None or ns is None:
        return None
    return _diazine_kind_of(ring, ns)


def _diazine_fg_block(info: dict) -> bool:
    """Block principal FGs for simple diazine (amine routes to pyrimidinamine)."""
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_ester") or info.get("has_amide"):
        return True
    if info.get("has_nitrile") or info.get("has_thiol") or info.get("has_amine"):
        return True
    return False


def _outside_c_set(mol: Mol, ring: set[int], skip: set[int]) -> set[int]:
    return {
        a.GetIdx() for a in mol.GetAtoms()
        if a.GetAtomicNum() == 6 and a.GetIdx() not in ring and a.GetIdx() not in skip
    }


def _diazine_side_methyl_ok(
    mol: Mol, ring: set[int], starts: list[int], skip: set[int],
) -> bool:
    """True when every ring C-side is pure -CH3 (outside C set == starts)."""
    if not starts:
        return not _outside_c_set(mol, ring, skip)
    return _outside_c_set(mol, ring, skip) == set(starts)


def _diazine_subs_ok(info: dict, mol: Mol, ring: set[int]) -> bool:
    """Allow ≤4 ring subs: halo / methyl / alkoxy / nitro / aryl."""
    alk, n_alk = _arene_alkoxy(info, ring)
    excl = alk | _aryl_exclude(info, ring)
    h = _ring_halo_n(mol, ring)
    starts = _ring_side_starts(mol, ring, excl)
    n_sub = h + len(starts) + _ring_nitro_n(info, ring) + n_alk + _aryl_sub_n(info, ring)
    if n_sub > 4:
        return False
    return True if not starts else _diazine_side_methyl_ok(mol, ring, starts, excl)


def _diazine_outside_ok(info: dict, mol: Mol, ring: set[int]) -> bool:
    alk = _arene_alkoxy(info, ring)[0]
    allowed = _ring_nitro_atoms(info, ring) | alk | _aryl_atoms(info, ring)
    return _outside_ok(mol, ring, allowed)


def _dz_ring_ok(info: dict, ring: list[int]) -> bool:
    mol, rs = info["mol"], set(ring)
    return _diazine_outside_ok(info, mol, rs) and _diazine_subs_ok(info, mol, rs)


def _pick_diazine_ring(info: dict) -> list[int] | None:
    cands = [r for r in _diazine_rings(info["mol"]) if _dz_ring_ok(info, r)]
    return min(cands, key=min) if cands else None


def _is_simple_diazine(info: dict) -> bool:
    if _diazine_fg_block(info):
        return False
    return _pick_diazine_ring(info) is not None


def _diazine_parent(info: dict) -> dict:
    ring = _pick_diazine_ring(info) or _diazine_rings(info["mol"])[0]
    ns = _diazine_n_idxs(info, ring) or []
    return {
        "chain": ring, "n_carbons": 6, "kind": _diazine_kind_of(ring, ns),
        "n_idxs": ns,
    }


def _try_diazine_parent(info: dict) -> dict | None:
    return _diazine_parent(info) if _is_simple_diazine(info) else None


def _pyrimidinamine_fg_block(info: dict) -> bool:
    """Block non-amine principal FGs for pyrimidinamine parent."""
    if info.get("has_acid") or info.get("has_aldehyde") or info.get("has_ketone"):
        return True
    if info.get("has_alcohol") or info.get("has_ester") or info.get("has_amide"):
        return True
    return bool(info.get("has_nitrile") or info.get("has_thiol"))


def _pyrimidinamine_subs_ok(info: dict, mol: Mol, ring: set[int], am_n: int) -> bool:
    """Extra ring subs: ≤2 of halo / methyl / methoxy (amine is principal FG)."""
    alk, n_alk = _arene_alkoxy(info, ring)
    h = _ring_halo_n(mol, ring)
    starts = _ring_side_starts(mol, ring, alk)
    if h + len(starts) + n_alk > 2:
        return False
    return _diazine_side_methyl_ok(mol, ring, starts, alk)


def _pyrimidinamine_outside(info: dict, mol: Mol, ring: set[int], am_n: int) -> bool:
    alk = _arene_alkoxy(info, ring)[0]
    allowed = {am_n} | alk
    return _outside_ok(mol, ring, allowed)


def _pyrimidinamine_ring_ok(info: dict, ring: list[int]) -> bool:
    mol, rs = info["mol"], set(ring)
    am = _mono_amine_on_ring(info, rs)
    if am is None or am.get("degree") != 1:
        return False
    if not _pyrimidinamine_outside(info, mol, rs, am["n_idx"]):
        return False
    return _pyrimidinamine_subs_ok(info, mol, rs, am["n_idx"])


def _pick_pyrimidinamine_ring(info: dict) -> list[int] | None:
    if _pyrimidinamine_fg_block(info):
        return None
    cands = [
        r for r in _diazine_rings(info["mol"])
        if _diazine_kind(info, r) == "pyrimidine" and _pyrimidinamine_ring_ok(info, r)
    ]
    return min(cands, key=min) if cands else None


def _is_simple_pyrimidinamine(info: dict) -> bool:
    return _pick_pyrimidinamine_ring(info) is not None


def _pyrimidinamine_parent(info: dict) -> dict:
    ring = _pick_pyrimidinamine_ring(info) or _diazine_rings(info["mol"])[0]
    return {
        "chain": ring, "n_carbons": 6, "kind": "pyrimidinamine",
        "n_idxs": _diazine_n_idxs(info, ring),
        "amine_c_idx": info["amines"][0]["c_idx"],
    }


def _try_pyrimidinamine_parent(info: dict) -> dict | None:
    if not _is_simple_pyrimidinamine(info):
        return None
    return _pyrimidinamine_parent(info)


def _is_diazole_core(info: dict) -> bool:
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


def _diazole_n_pair(info: dict, nn_dist: int) -> tuple[int, int] | None:
    """Return (nh_idx, n_idx) for 1H-diazole with given N–N ring distance."""
    if not _is_diazole_core(info):
        return None
    mol, ring = info["mol"], list(info["rings"][0]["atom_ids"])
    ns = _ring_n_idxs(mol, ring)
    if len(ns) != 2 or _ring_nn_dist(ring, ns) != nn_dist:
        return None
    nh = _nh_among(mol, ns)
    return None if nh is None else (nh, ns[0] if ns[1] == nh else ns[1])


def _imidazole_n_pair(info: dict) -> tuple[int, int] | None:
    """Return (nh_idx, n_idx) for 1H-imidazole (1,3-diazole); else None."""
    return _diazole_n_pair(info, 2)


def _pyrazole_n_pair(info: dict) -> tuple[int, int] | None:
    """Return (nh_idx, n_idx) for 1H-pyrazole (1,2-diazole); else None."""
    return _diazole_n_pair(info, 1)


def _is_simple_diazole(info: dict, pair_fn) -> bool:
    if pair_fn(info) is None or _hetero5_fg_block(info):
        return False
    mol, ring = info["mol"], set(info["rings"][0]["atom_ids"])
    if not _outside_ok(mol, ring, _aryl_atoms(info, ring)):
        return False
    return _hetero5_subs_ok(info, mol, ring)


def _diazole_parent(info: dict, kind: str, pair_fn) -> dict:
    ring = list(info["rings"][0]["atom_ids"])
    nh, n = pair_fn(info)
    return {
        "chain": ring, "n_carbons": 5, "kind": kind,
        "nh_idx": nh, "n_idx": n, "n_idxs": [nh, n],
    }


def _is_simple_imidazole(info: dict) -> bool:
    return _is_simple_diazole(info, _imidazole_n_pair)


def _imidazole_parent(info: dict) -> dict:
    return _diazole_parent(info, "imidazole", _imidazole_n_pair)


def _try_imidazole_parent(info: dict) -> dict | None:
    return _imidazole_parent(info) if _is_simple_imidazole(info) else None


def _is_simple_pyrazole(info: dict) -> bool:
    return _is_simple_diazole(info, _pyrazole_n_pair)


def _pyrazole_parent(info: dict) -> dict:
    return _diazole_parent(info, "pyrazole", _pyrazole_n_pair)


def _try_pyrazole_parent(info: dict) -> dict | None:
    return _pyrazole_parent(info) if _is_simple_pyrazole(info) else None
