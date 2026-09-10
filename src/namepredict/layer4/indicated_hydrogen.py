"""L4 指示氢位次计算（P-58.2.1）：环系中仅以单键连接相邻环原子的饱和环位标 'H'。"""
from __future__ import annotations

from rdkit import Chem


def saturated_ring_atoms(mol, ring_atoms: set[int], exclude: frozenset[int] = frozenset()) -> list[int]:
    """返回环系内仅以单键连接相邻环原子、且带氢的饱和环位（指示氢候选原子）；部分饱和的 mancude 环系（含芳香位）适用，exclude 中的加氢位已由 hydro 前缀表达故排除。"""
    if mol is None or not ring_atoms:
        return []
    ring_atoms = {i for i in ring_atoms if i < mol.GetNumAtoms() and mol.GetAtomWithIdx(i).IsInRing()}
    if not ring_atoms:
        return []
    if not any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms):  # P-58.2.1: 指示氢只标 mancude（含最大非累积双键）环系；全饱和多环（decalin）由 hydro 前缀承载，否则每个饱和碳都会被标 H（benzoxane 类会产出 8H,7H,6H… 的荒谬串）
        return []
    kek = Chem.Mol(mol)
    try:
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        kek = mol
    out = []
    for i in sorted(ring_atoms):
        if i in exclude:  # 该位已由 hydro 前缀表达（如 2,3 位），不再标指示氢
            continue
        if mol.GetAtomWithIdx(i).GetTotalNumHs() == 0:  # H 数取原始 mol：Kekulize 清芳香会给吡啶型 N 补隐式 H，否则其被误判为饱和位（1H-pyridine）
            continue
        atom = kek.GetAtomWithIdx(i)
        if all(b.GetBondType() == Chem.BondType.SINGLE
               for b in atom.GetBonds() if b.GetOtherAtomIdx(i) in ring_atoms):
            out.append(i)
    return out


def indicated_hydrogen(mol, chain, labels=None, exclude=frozenset()) -> list[str]:
    """返回指示氢标签列表（如 ['1H']、['9H']）：位次取整体编号 labels，缺失时用链序号；exclude 为已用 hydro 表达的加氢位。"""
    chain = list(chain or ())
    if not chain:
        return []
    use_labels = bool(labels) and len(labels) == len(chain)
    sats = sorted(saturated_ring_atoms(mol, set(chain), exclude), key=chain.index)
    if len(sats) > 1 and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in sats):  # 互变异构冗余护栏：饱和位全为氮且多于一个时只保留最低位次——亚胺-胺式 SMILES 会让咪唑环出现两个 [nH]，而标准形式（如 1H-imidazo[4,5-c]pyridine）只标一个
        sats = sats[:1]
    out = []
    for atom in sats:
        locant = labels[chain.index(atom)] if use_labels else chain.index(atom) + 1
        out.append(f"{locant}H")
    return out


def indicated_hydrogen_prefix(mol, chain, labels=None, exclude=frozenset()) -> str:
    """返回可直接置于母体名前的指示氢前缀（'1H-' / '1H,2H-'），无则空串；exclude 为已用 hydro 表达的加氢位。"""
    locants = indicated_hydrogen(mol, chain, labels, exclude)
    return f"{','.join(locants)}-" if locants else ""
