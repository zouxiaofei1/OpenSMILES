"""构建诱导切割子分子，连接点用 H 封端以用于 free-name。"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Mol


@dataclass(frozen=True)
class CutSubmol:
    mol: object  # RDKit Mol，连接点用 H 封端以用于 free-name
    atom_map: dict[int, int]  # new_idx -> old_idx（新索引→旧索引）
    inv_map: dict[int, int]  # old_idx -> new_idx（旧索引→新索引）
    attach_new: int  # 子分子中的连接原子索引
    attach_old: int
    atoms_old: frozenset[int]


def _ordered(atoms: frozenset[int]) -> list[int]:
    return sorted(atoms)


def _copy_atoms(em: Chem.RWMol, mol: Mol, order: list[int]) -> dict[int, int]:
    inv: dict[int, int] = {}
    for old in order:
        inv[old] = em.AddAtom(mol.GetAtomWithIdx(old))
    return inv


def _copy_bonds(em: Chem.RWMol, mol: Mol, inv: dict[int, int]) -> None:
    for old_a, new_a in inv.items():
        for bond in mol.GetAtomWithIdx(old_a).GetBonds():
            old_b = bond.GetOtherAtomIdx(old_a)
            if old_b not in inv or old_b < old_a:
                continue
            em.AddBond(new_a, inv[old_b], bond.GetBondType())


def _cap_attach_h(em: Chem.RWMol, attach_new: int) -> None:
    atom = em.GetAtomWithIdx(attach_new)
    atom.SetNoImplicit(False)
    atom.UpdatePropertyCache(strict=False)


def _sanitize(em: Chem.RWMol) -> Mol | None:
    try:
        Chem.SanitizeMol(em)
    except Exception:
        return None
    return em.GetMol()


def _pack(out: Mol, inv: dict[int, int], attach_old: int, atoms: frozenset[int]) -> CutSubmol:
    atom_map = {n: o for o, n in inv.items()}
    return CutSubmol(out, atom_map, inv, inv[attach_old], attach_old, frozenset(atoms))


def build_cut_submol(
    mol: Mol, atoms: frozenset[int], attach_old: int,
) -> CutSubmol | None:
    """在 atoms 上的诱导子分子；连接处的自由价用 H 填充。"""
    if attach_old not in atoms:
        return None
    em = Chem.RWMol()
    inv = _copy_atoms(em, mol, _ordered(atoms))
    _copy_bonds(em, mol, inv)
    _cap_attach_h(em, inv[attach_old])
    out = _sanitize(em)
    return None if out is None else _pack(out, inv, attach_old, atoms)


def _add_anchor(em: Chem.RWMol, attach_new: int) -> None:
    #用 dummy 原子（`*`）标记连接原子
    d = em.AddAtom(Chem.Atom(0))
    em.AddBond(attach_new, d, Chem.BondType.SINGLE)


def build_anchor_submol(mol: Mol, atoms: frozenset[int], attach_old: int) -> Mol | None:
    #在 atoms 上的诱导子分子，连接位点用 dummy 原子标记。
    if attach_old not in atoms:
        return None
    em = Chem.RWMol()
    inv = _copy_atoms(em, mol, _ordered(atoms))
    _copy_bonds(em, mol, inv)
    _add_anchor(em, inv[attach_old])
    return _sanitize(em)
