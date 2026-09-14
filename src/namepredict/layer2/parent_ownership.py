"""末端母体所有权：不可变 owned_atoms 的最终确定。"""
from __future__ import annotations

from rdkit.Chem import Mol


def _chain_atoms(parent: dict) -> set[int]:
    """取母体链原子集合。"""
    return set(parent.get("chain") or ())


def _kind_fg_atoms(parent: dict, mol: Mol) -> set[int]:
    """主官能团所有权原子：骨架内（或邻骨架）锚点及其特征原子。"""
    facts = parent.get("principal_expression_facts")
    occurrences = parent.get("principal_occurrences") or ()
    if facts is None:
        return set()
    chain = _chain_atoms(parent)
    anchors = {a for o in occurrences for a in o.parent_anchors}
    atoms = {a for o in occurrences for a in o.characteristic_atoms} or set(facts.characteristic_atoms)
    seeds = anchors & chain
    if not seeds:  # 锚点全在骨架外：exocyclic 基团（苯甲酸的羧基），改取与骨架相邻的锚点
        linked = {n.GetIdx() for i in chain for n in mol.GetAtomWithIdx(i).GetNeighbors()}
        seeds = anchors & linked
    out = set(seeds)
    for i in tuple(out):
        out |= {n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
                if n.GetIdx() in atoms and n.GetIdx() not in anchors}
    return out


def compute_owned_atoms(parent: dict, mol: Mol) -> frozenset[int]:
    """链与主官能团特征原子的并集（末端所有权集合）。"""
    return frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))


def finalize_parent_ownership(parent: dict, mol: Mol) -> dict:
    """一次性复制候选，生成不可变 owned_atoms frozenset。"""
    if isinstance(parent.get("owned_atoms"), frozenset):
        return parent
    return {**parent, "owned_atoms": compute_owned_atoms(parent, mol)}
