"""L4 指示氢位次计算（P-58.2.1）：饱和环位标 'H'。"""
from __future__ import annotations

from rdkit import Chem

from namepredict.layer1.ring_systems import kekulized
from namepredict.layer4.locant_calc import atom_locant


def saturated_ring_atoms(mol, ring_atoms: set[int], exclude: frozenset[int] = frozenset()) -> list[int]:
    """返回仅单键连邻环原子且带氢的饱和环位；exclude 位由 hydro 表达。"""
    if mol is None or not ring_atoms:
        return []
    ring_atoms = {i for i in ring_atoms if i < mol.GetNumAtoms() and mol.GetAtomWithIdx(i).IsInRing()}
    if not ring_atoms:
        return []
    kek = kekulized(mol)
    if kek is None:
        kek = mol  # Kekulize 失败：芳香键仍是 AROMATIC，该位指示氢静默放弃
    out = []
    for i in sorted(ring_atoms):
        if i in exclude:  # 该位已由 hydro 前缀表达（如 2,3 位），不再标指示氢
            continue
        if mol.GetAtomWithIdx(i).GetTotalNumHs() == 0:  # H 数取原始 mol：Kekulize 会给吡啶型 N 补隐式 H。
            continue
        atom = kek.GetAtomWithIdx(i)
        if all(b.GetBondType() == Chem.BondType.SINGLE
               for b in atom.GetBonds() if b.GetOtherAtomIdx(i) in ring_atoms):
            out.append(i)
    return out


def indicated_hydrogen(mol, chain, labels=None, exclude=frozenset(), extra=frozenset()) -> list[str]:
    """返回指示氢位次列表：位次取 labels，缺失时用链序号（P-58.2.1）。"""
    chain = list(chain or ())
    if not chain:
        return []
    facts = {"labels": labels}
    sats = set(saturated_ring_atoms(mol, set(chain), exclude)) | {
        i for i in extra if i in chain and i not in exclude}
    sats = sorted(sats, key=chain.index)
    if len(sats) > 1 and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in sats):  # 互变异构冗余护栏：饱和位全为氮且多于一个时只留最低位次。
        sats = sats[:1]
    out = []
    for atom in sats:
        out.append(str(atom_locant(chain, atom, facts)))
    return out
