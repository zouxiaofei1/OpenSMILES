from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.constants import Br, C, Cl, F, H, I, O

from namepredict.layer2.aryl_sub import (
    _arom_c6_ring_lists, _is_unfused_benzene_ring,
)
from namepredict.layer2.side_alkyl import (
    _disjoint_cover, _is_cf3_carbon, _is_cf3_fluoro, _is_omega_halo_c,
    _is_side_halo, _outer_alkoxy_n, _outer_atoms, _side_covers, _side_sets,
)

def _all_carbons_are_c(mol: Mol, atom_ids: tuple) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == C for i in atom_ids)

def _bond_between(mol: Mol, a: int, b: int):
    return mol.GetBondBetweenAtoms(a, b)

def _ring_bonds_single(mol: Mol, atom_ids: tuple) -> bool:
    ids = list(atom_ids)
    for i, a in enumerate(ids):
        b = ids[(i + 1) % len(ids)]
        bond = _bond_between(mol, a, b)
        if bond is None or bond.GetBondType().name != "SINGLE":
            return False
    return True

def _outside_carbons(
    mol: Mol, ring_set: set[int], exclude: set[int] | None = None,
) -> list[int]:
    skip = exclude or set()
    return [
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == C
        and a.GetIdx() not in ring_set
        and a.GetIdx() not in skip
    ]

def _pure_c_bonds(mol: Mol, idx: int) -> bool:
    atom = mol.GetAtomWithIdx(idx)
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        z = other.GetAtomicNum()
        if z not in (1, 6) or (z != 1 and bond.GetBondType().name != "SINGLE"):
            return False
    return True

def _pure_alkyl_outside(mol: Mol, outside: list[int]) -> bool:
    return all(
        _is_cf3_carbon(mol, i) or _is_omega_halo_c(mol, i) or _pure_c_bonds(mol, i)
        for i in outside
    )

def _ring_side_starts(
    mol: Mol, ring_set: set[int], exclude: set[int] | None = None,
) -> list[int]:
    skip = exclude or set()
    starts: list[int] = []
    for r in ring_set:
        for n in mol.GetAtomWithIdx(r).GetNeighbors():
            if n.GetAtomicNum() == C and n.GetIdx() not in ring_set:
                if n.GetIdx() not in skip:
                    starts.append(n.GetIdx())
    return starts

def _is_ring_halo(atom, ring_set: set[int]) -> bool:
    if atom.GetAtomicNum() not in (F, Cl, Br, I):
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != H]
    return len(heavies) == 1 and heavies[0].GetIdx() in ring_set

def _ring_nitro_n(info: dict, ring_set: set[int]) -> int:
    return sum(1 for n in (info.get("nitros") or []) if n["c_idx"] in ring_set)

def _ring_nitro_atoms(info: dict, ring_set: set[int]) -> set[int]:
    out: set[int] = set()
    for n in info.get("nitros") or []:
        if n["c_idx"] in ring_set:
            out.add(n["n_idx"])
            out.update(n.get("o_idxs") or [])
    return out

def _ring_iso_entries(info: dict, ring_set: set[int]) -> list[dict]:
    keys = ("isocyanates", "isothiocyanates")
    return [e for k in keys for e in (info.get(k) or []) if e["r_c_idx"] in ring_set]

def _ring_iso_n(info: dict, ring_set: set[int]) -> int:
    return len(_ring_iso_entries(info, ring_set))

def _ring_iso_atoms(info: dict, ring_set: set[int]) -> set[int]:
    out: set[int] = set()
    for e in _ring_iso_entries(info, ring_set):
        out.update((e["n_idx"], e["c_idx"], e["x_idx"]))
    return out

def _ring_alkoxy_pair(e: dict, ring_set: set[int]) -> tuple[int, int, int] | None:
    c1, c2, o = e["c1"], e["c2"], e["o_idx"]
    if (c1 in ring_set) == (c2 in ring_set):
        return None
    return (o, c1, c2) if c1 in ring_set else (o, c2, c1)

def _ring_alkoxy_one(mol: Mol, e: dict, ring_set: set[int]) -> dict | None:
    pair = _ring_alkoxy_pair(e, ring_set)
    if pair is None or mol.GetAtomWithIdx(pair[2]).GetIsAromatic():
        return None
    o, ring_c, outer = pair
    n = _outer_alkoxy_n(mol, outer, o)
    if n == 0:
        return None
    return {"o_idx": o, "ring_c": ring_c, "outer_c": outer, "n": n,
            "atoms": _outer_atoms(mol, outer, o, n)}

def _ring_alkoxy_ethers(info: dict, ring_set: set[int]) -> list[dict]:
    mol: Mol = info["mol"]
    out: list[dict] = []
    for e in info.get("ethers") or []:
        one = _ring_alkoxy_one(mol, e, ring_set)
        if one is not None:
            out.append(one)
    return out

def _ring_alkoxy_atoms(info: dict, ring_set: set[int]) -> set[int]:
    out: set[int] = set()
    for a in _ring_alkoxy_ethers(info, ring_set):
        out.add(a["o_idx"])
        out.update(a["atoms"])
    return out

def _outside_hetero_ok(atom, ring_set: set[int], allow: set[int]) -> bool:
    if atom.GetAtomicNum() in (1, 6) or atom.GetIdx() in ring_set:
        return True
    if atom.GetIdx() in allow or _is_ring_halo(atom, ring_set):
        return True
    return _is_cf3_fluoro(atom) or _is_side_halo(atom)

def _no_hetero_outside(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    allow = allowed or set()
    return all(_outside_hetero_ok(a, ring_set, allow) for a in mol.GetAtoms())

def _outside_ok(mol: Mol, ring_set: set[int], allowed: set[int] | None = None) -> bool:
    if not _no_hetero_outside(mol, ring_set, allowed):
        return False
    return _pure_alkyl_outside(mol, _outside_carbons(mol, ring_set, allowed or set()))

def _ring_double_count(mol: Mol, atom_ids: tuple) -> int:
    ids = list(atom_ids)
    n = 0
    for i, a in enumerate(ids):
        b = ids[(i + 1) % len(ids)]
        bond = _bond_between(mol, a, b)
        if bond is not None and bond.GetBondType().name == "DOUBLE":
            n += 1
    return n

def _is_carbocycle_ring(info: dict) -> tuple | None:
    rings = info.get("rings") or []
    if len(rings) != 1:
        return None
    mol: Mol = info["mol"]
    atom_ids = rings[0]["atom_ids"]
    if not _all_carbons_are_c(mol, atom_ids):
        return None
    return atom_ids

def _is_cycloalkene_core(info: dict) -> bool:
    atom_ids = _is_carbocycle_ring(info)
    if atom_ids is None:
        return False
    return _ring_double_count(info["mol"], atom_ids) == 1

def _endocyclic_double(info: dict, ring_set: set[int]) -> tuple[int, int] | None:
    bonds = info.get("double_bonds") or []
    if len(bonds) != 1:
        return None
    c1, c2 = bonds[0]["c1"], bonds[0]["c2"]
    if c1 in ring_set and c2 in ring_set:
        return c1, c2
    return None

def _endocyclic_doubles(info: dict, ring_set: set[int]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for db in info.get("double_bonds") or []:
        c1, c2 = db["c1"], db["c2"]
        if c1 in ring_set and c2 in ring_set:
            out.append((c1, c2))
    return out

def _is_cyclopolyene_core(info: dict) -> bool:
    atom_ids = _is_carbocycle_ring(info)
    if atom_ids is None or info.get("triple_bonds"):
        return False
    return _ring_double_count(info["mol"], atom_ids) >= 2

def _mono_oh_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    hydroxyls = info.get("hydroxyls") or []
    if len(hydroxyls) != 1:
        return None
    oh = hydroxyls[0]
    if oh["c_idx"] not in ring_set:
        return None
    return oh

def _cyclo_fg_sides_ok(mol: Mol, ring_set: set[int]) -> bool:
    """Outside C must be fully covered by claimable alkyl side probes."""
    starts = _ring_side_starts(mol, ring_set)
    outside = set(_outside_carbons(mol, ring_set))
    if not starts:
        return not outside
    sets = _side_sets(mol, ring_set, starts)
    return sets is not None and _disjoint_cover(sets, outside)

def _cyclo_fg_parent_ok(info: dict, allowed: set[int]) -> bool:
    """Sat mono carbocycle FG parent: FG heteros + ring halo; claimable sides."""
    if not _is_cycloalkane_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if not _hetero_or_ring_halo(mol, ring_set, allowed):
        return False
    return _cyclo_fg_sides_ok(mol, ring_set)

def _cyclo_ene_fg_ok(info: dict, allowed: set[int]) -> bool:
    """Mono cycloalkene + FG heteros/halo; claimable alkyl sides."""
    if not _is_cycloalkene_core(info):
        return False
    mol: Mol = info["mol"]
    ring_set = set(info["rings"][0]["atom_ids"])
    if _endocyclic_double(info, ring_set) is None:
        return False
    if not _hetero_or_ring_halo(mol, ring_set, allowed):
        return False
    return _cyclo_fg_sides_ok(mol, ring_set)

def _mono_amine_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    amines = info.get("amines") or []
    if len(amines) != 1:
        return None
    am = amines[0]
    c_idx = am.get("c_idx")
    if c_idx is None or c_idx not in ring_set:
        return None
    return am

def _dbl_o_idx(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != "DOUBLE":
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == O:
            return other.GetIdx()
    return None

def _mono_ketone_on_ring(info: dict, ring_set: set[int]) -> dict | None:
    ketones = info.get("ketones") or []
    if len(ketones) != 1:
        return None
    ket = ketones[0]
    if ket["c_idx"] not in ring_set:
        return None
    return ket

def _is_simple_cycloketone(info: dict) -> bool:
    ring_set = set((info.get("rings") or [{}])[0].get("atom_ids") or [])
    ket = _mono_ketone_on_ring(info, ring_set)
    if ket is None:
        return False
    o_idx = _dbl_o_idx(info["mol"], ket["c_idx"])
    return o_idx is not None and _cyclo_fg_parent_ok(info, {o_idx})

def _is_benzene_core(info: dict) -> bool:
    """True if mol has at least one unfused aromatic C6 carbocycle."""
    mol: Mol = info["mol"]
    for ring in _arom_c6_ring_lists(mol):
        if _is_unfused_benzene_ring(mol, set(ring)):
            return True
    return False

def _ring_halo_n(mol: Mol, ring_set: set[int]) -> int:
    return sum(1 for a in mol.GetAtoms() if _is_ring_halo(a, ring_set))

def _benzene_alkyl_ns(
    mol: Mol, ring_set: set[int], starts: list[int],
    exclude: set[int] | None = None,
) -> list[int]:
    outside = set(_outside_carbons(mol, ring_set, exclude))
    sets = _side_sets(mol, ring_set, starts)
    if sets is None or not _disjoint_cover(sets, outside):
        return []
    return [len(s) for s in sets]

def _mono_benzene_ok(
    mol: Mol, ring_set: set[int], starts: list[int],
    exclude: set[int] | None = None,
) -> bool:
    if not starts:
        return True
    if len(starts) != 1:
        return False
    outside = set(_outside_carbons(mol, ring_set, exclude))
    return _side_covers(mol, starts[0], ring_set, outside)

def _multi_benzene_ok(
    mol: Mol, ring_set: set[int], starts: list[int],
    exclude: set[int] | None = None,
) -> bool:
    ns = _benzene_alkyl_ns(mol, ring_set, starts, exclude)
    return len(ns) == len(starts)

def _benzene_subs_ok(
    mol: Mol, ring_set: set[int], n_nitro: int = 0, n_alkoxy: int = 0,
    exclude: set[int] | None = None, n_aryl: int = 0, n_iso: int = 0,
) -> bool:
    h = _ring_halo_n(mol, ring_set)
    starts = _ring_side_starts(mol, ring_set, exclude)
    n_sub = h + len(starts) + n_nitro + n_alkoxy + n_aryl + n_iso
    if n_sub > 4:
        return False
    if n_sub <= 1:
        return _mono_benzene_ok(mol, ring_set, starts, exclude)
    return _multi_benzene_ok(mol, ring_set, starts, exclude)

def _arene_alkoxy(info: dict, ring_set: set[int]) -> tuple[set[int], int]:
    alk = _ring_alkoxy_atoms(info, ring_set)
    return alk, len(_ring_alkoxy_ethers(info, ring_set))

def _hetero_or_ring_halo(mol: Mol, ring_set: set[int], allowed: set[int]) -> bool:
    return all(_outside_hetero_ok(a, ring_set, allowed) for a in mol.GetAtoms())

def _arene_fg_subs_ok(
    mol: Mol, ring_set: set[int], allowed: set[int], n_nitro: int = 0,
    n_amino: int = 0, n_alkoxy: int = 0, exclude: set[int] | None = None,
    n_aryl: int = 0,
) -> bool:
    if not _hetero_or_ring_halo(mol, ring_set, allowed):
        return False
    starts = _ring_side_starts(mol, ring_set, exclude)
    if len(_benzene_alkyl_ns(mol, ring_set, starts, exclude)) != len(starts):
        return False
    n = _ring_halo_n(mol, ring_set) + len(starts) + n_nitro + n_amino + n_alkoxy + n_aryl
    return n <= 4

def _ring_primary_amines(info: dict, ring_set: set[int]) -> list:
    return [
        a for a in (info.get("amines") or [])
        if a.get("degree") == 1 and a.get("c_idx") in ring_set
    ]

def _phenol_amines_ok(info: dict, ring_set: set[int]) -> list | None:
    ring_ams = _ring_primary_amines(info, ring_set)
    if len(info.get("amines") or []) != len(ring_ams):
        return None
    return ring_ams

def _di_oh_on_ring(info: dict, ring_set: set[int]) -> list[dict] | None:
    hydroxyls = info.get("hydroxyls") or []
    if len(hydroxyls) != 2 or any(h["c_idx"] not in ring_set for h in hydroxyls):
        return None
    return hydroxyls

def _is_methyl_on_ring(mol: Mol, s: int, ring_set: set[int]) -> bool:
    if _is_cf3_carbon(mol, s):
        return True
    atom = mol.GetAtomWithIdx(s)
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z == 1 or (z == 6 and n.GetIdx() in ring_set):
            continue
        return False
    return True

# Multi-ring cyclo parent pick lives in cyclo_pick (late re-export avoids cycles).
from namepredict.layer2.scaffold.cyclo_pick import (  # noqa: E402
    _is_cycloalkane_core,
    _is_simple_cycloalkane,
    _pick_cycloalkane_ring,
)
