"""L4 指示氢位次计算（P-58.2.1）：饱和环位标 'H'。"""
from __future__ import annotations

from rdkit import Chem

from namepredict.layer1.ring_systems import kekulized
from namepredict.layer4.locant_calc import atom_locant


def _flanked_by_exo_double(mol, idx: int, ring_atoms: set[int]) -> bool:
    """环碳两侧环邻位是否都带环外杂原子多重键（如 1,3-二酮的 C2，P-14.4）。"""
    atom = mol.GetAtomWithIdx(idx)
    if atom.GetAtomicNum() != 6:
        return False
    n_exo = 0
    for b in atom.GetBonds():
        other = b.GetOtherAtomIdx(idx)
        if other not in ring_atoms:
            continue
        oa = mol.GetAtomWithIdx(other)
        if any(ob.GetBondType() != Chem.BondType.SINGLE and ob.GetOtherAtomIdx(other) not in ring_atoms
               and mol.GetAtomWithIdx(ob.GetOtherAtomIdx(other)).GetAtomicNum() != 6
               for ob in oa.GetBonds()):
            n_exo += 1
    return n_exo >= 2


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
        if _flanked_by_exo_double(kek, i, ring_atoms):  # 两侧环位已被 =O 类占满：不饱和度已由后缀确定
            continue
        atom = kek.GetAtomWithIdx(i)
        if all(b.GetBondType() == Chem.BondType.SINGLE
               for b in atom.GetBonds() if b.GetOtherAtomIdx(i) in ring_atoms):
            out.append(i)
    return out


def _ring_double_bonds(mol, ring_atoms: set[int]) -> int:
    """环内多重键计数（kekulé 形）：全饱和环与单双键碳环都不需要指示氢。"""
    if mol is None:
        return 0
    n = 0
    for i in ring_atoms:
        atom = mol.GetAtomWithIdx(i)
        for b in atom.GetBonds():
            if b.GetOtherAtomIdx(i) in ring_atoms and i < b.GetOtherAtomIdx(i) \
                    and b.GetBondType() != Chem.BondType.SINGLE:
                n += 1
    return n


def _is_monocycle(mol, chain: list[int]) -> bool:
    """链上原子是否构成单个环（每个原子在链内恰有两个邻位）。"""
    if mol is None or len(chain) < 3:
        return False
    ring = set(chain)
    for i in ring:
        n_in = sum(1 for b in mol.GetAtomWithIdx(i).GetBonds()
                   if b.GetOtherAtomIdx(i) in ring)
        if n_in != 2:
            return False
    return True


def _is_retained_scaffold(scaffold_id: str | None) -> bool:
    """母体是否命中共有保留名模板（萘/噻吨等 mancude 系统）。"""
    if not scaffold_id:
        return False
    from namepredict.layer2.ring_scaffold import get_spec

    return get_spec(scaffold_id) is not None


def indicated_hydrogen(mol, chain, labels=None, exclude=frozenset(), extra=frozenset(),
                       scaffold_id: str | None = None) -> list[str]:
    """返回指示氢位次列表：位次取 labels，缺失时用链序号（P-58.2.1）。"""
    chain = list(chain or ())
    if not chain:
        return []
    facts = {"labels": labels}
    sats = set(saturated_ring_atoms(mol, set(chain), exclude)) | {
        i for i in extra if i in chain and i not in exclude}
    if not extra and _is_monocycle(mol, chain) \
            and _ring_double_bonds(kekulized(mol) or mol, set(chain)) == 1 \
            and not _is_retained_scaffold(scaffold_id):
        sats = set()  # 单环仅一个环内双键（环己烯/环戊烯）：氢位无歧义
    sats = sorted(sats, key=chain.index)
    if len(sats) > 1 and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 7 for i in sats):  # 互变异构冗余护栏：饱和位全为氮且多于一个时只留最低位次。
        sats = sats[:1]
    out = []
    for atom in sats:
        out.append(str(atom_locant(chain, atom, facts)))
    return out
