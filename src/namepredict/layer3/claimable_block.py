"""仅所有权可 claim 的侧链块：拓扑、连接、slot——不命名。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rdkit.Chem import BondType, Mol

from namepredict.tools.block_cut import cut_block, side_roots

class SideSlot(str, Enum):
    """侧链连接原子的角色槽位（链碳/环碳/胺 N/其他）。"""
    CHAIN_C = "chain_c"
    RING_C = "ring_c"
    AMINE_N = "amine_n"
    OTHER = "other"


@dataclass(frozen=True)
class ClaimedBlock:
    """可 claim 的外部重原子组分：槽位、连接原子、根原子与原子集。"""
    slot: SideSlot
    attach_parent: int
    root: int
    atoms: frozenset[int]


def _is_amine_n(mol: Mol, n_idx: int) -> bool:
    """胺 N：非芳香、非环员 N（环 N 用环上位次定位）。"""
    atom = mol.GetAtomWithIdx(n_idx)
    return not (atom.GetAtomicNum() != 7 or atom.GetIsAromatic() or atom.IsInRing())


def derive_slot(mol: Mol, attach_parent: int) -> SideSlot:
    """仅从所属连接原子的角色推导 SideSlot。"""
    if _is_amine_n(mol, attach_parent):
        return SideSlot.AMINE_N
    atom = mol.GetAtomWithIdx(attach_parent)
    if atom.GetAtomicNum() == 6:
        return SideSlot.RING_C if atom.IsInRing() else SideSlot.CHAIN_C
    return SideSlot.OTHER


def _attach_parents_of(mol: Mol, atoms: frozenset[int], owned: frozenset[int]) -> set[int]:
    """返回组分原子在所有权集合中的连接点集合。"""
    out: set[int] = set()
    for a in atoms:
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            if n.GetAtomicNum() != 1 and n.GetIdx() in owned:
                out.add(n.GetIdx())
    return out

def claim_block(
    mol: Mol,
    *,
    owned_atoms: frozenset[int],
    attach_parent: int,
    root: int,
    slot: SideSlot,
) -> ClaimedBlock | None:
    """返回完整的外部重原子组分及其连接点，若无则返回 None。"""
    if attach_parent not in owned_atoms:
        return None

    atoms = cut_block(mol, root, owned_atoms)
    if not atoms or len(_attach_parents_of(mol, atoms, owned_atoms)) != 1:
        return None
    return ClaimedBlock(slot=slot, attach_parent=attach_parent, root=root, atoms=atoms)


def _canonical_edge(
    mol: Mol, atoms: frozenset[int], owned: frozenset[int]
) -> tuple[int, int] | None:
    """所属与该组分之间最小（attach_parent, root）边。"""
    edges = [
        (n.GetIdx(), a)
        for a in atoms
        for n in mol.GetAtomWithIdx(a).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() in owned
    ]
    return min(edges) if edges else None

def _has_dbl_o_edge(mol: Mol, atoms: frozenset[int], owned: frozenset[int]) -> bool:
    """外部组分是否有双键连 owned 重原子的氧（主 FG 成分不 claim）。"""
    for a in atoms:
        if mol.GetAtomWithIdx(a).GetAtomicNum() != 8:
            continue
        for n in mol.GetAtomWithIdx(a).GetNeighbors():
            if n.GetAtomicNum() != 1 and n.GetIdx() in owned:
                bond = mol.GetBondBetweenAtoms(a, n.GetIdx())
                if bond is not None and bond.GetBondType() == BondType.DOUBLE and  n.GetAtomicNum() != 6:
                    return True
    return False


def _try_claim(
    mol: Mol, owned: frozenset[int], atoms: frozenset[int]
) -> ClaimedBlock | None:
    """尝试为单一组分建立 claim 并返回其块。"""
    edge = _canonical_edge(mol, atoms, owned)
    if edge is None:
        return None
    if _has_dbl_o_edge(mol, atoms, owned):
        return None
    attach, root = edge
    slot = derive_slot(mol, attach)
    return claim_block(
        mol, owned_atoms=owned, attach_parent=attach, root=root, slot=slot
    )


def _unique_components(mol: Mol, owned: frozenset[int]) -> list[frozenset[int]]:
    """去重枚举所有权之外的重原子连通组分。"""
    seen: set[frozenset[int]] = set()
    out: list[frozenset[int]] = []
    for root in side_roots(mol, owned):
        atoms = cut_block(mol, root, owned)
        if atoms and atoms not in seen:
            seen.add(atoms)
            out.append(atoms)
    return out

def iter_claims(mol: Mol, owned_atoms: frozenset[int]) -> list[ClaimedBlock]:
    """按 canonical 顺序返回所有外部重原子组分的 claim。"""
    claims = [
        c
        for atoms in _unique_components(mol, owned_atoms)
        if (c := _try_claim(mol, owned_atoms, atoms)) is not None
    ]
    return sorted(claims, key=lambda c: (c.attach_parent, c.root, c.slot.value))
