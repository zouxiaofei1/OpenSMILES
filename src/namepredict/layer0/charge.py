"""L0 酸性质子重定位/电荷归一化"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from namepredict.constants import ACCEPTOR_Z, ACID_CENTERS, DONOR_KIND, O


def _oxo_neighbor(atom, heavy: int):
    """返回邻接的非芳香成酸中心原子（双键氧 ≥1）；无则 None。"""
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() != heavy or n.GetIsAromatic():
            continue
        oxo = sum(
            1
            for b in n.GetBonds()
            if b.GetOtherAtom(n).GetAtomicNum() == O
            and b.GetBondType() == Chem.BondType.DOUBLE
        )
        if oxo >= 1:
            return n
    return None

def _acid_kind(atom) -> str | None:
    """按 ACID_CENTERS 表返回 O 所连成酸中心对应的酸类名。"""
    for z, kind in ACID_CENTERS.items():
        if _oxo_neighbor(atom, z):
            return kind
    return None

def _acid_kind_of_oh(atom) -> str | None:
    """判中性含 H 的 O 是否为质子化强酸 OH；否则 None。"""
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() < 1:
        return None
    return _acid_kind(atom)

def _is_weak_anion(atom) -> bool:
    """判定位点是否为去质子化的弱酸位（可接受质子）。"""
    if atom.GetFormalCharge() != -1 or atom.GetAtomicNum() not in ACCEPTOR_Z:
        return False
    if any(n.GetFormalCharge() == 1 for n in atom.GetNeighbors()):  # 邻接 +1 电荷 → 硝基/N-氧化物等内平衡写法，不当作弱酸位
        return False
    if atom.GetAtomicNum() == O and _acid_kind(atom) is not None:  # 已是强酸共轭碱（羧酸/磷酸/磺酸根），电荷位置合理
        return False
    return True


def _relocate_proton(mol: Mol, a_idx: int, d_idx: int) -> Mol | None:
    """把质子从强酸供体搬到弱酸受体，只做电荷/H 记账；消毒失败返回 None。"""
    m = RWMol(mol)
    acc = m.GetAtomWithIdx(a_idx)
    don = m.GetAtomWithIdx(d_idx)
    acc.SetFormalCharge(0)
    acc.SetNumExplicitHs(acc.GetNumExplicitHs() + 1)
    acc.SetNoImplicit(True)
    don.SetFormalCharge(-1)
    don.SetNumExplicitHs(max(0, don.GetNumExplicitHs() - 1))
    don.SetNoImplicit(True)
    out = m.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    Chem.AssignStereochemistry(out, force=True, cleanIt=False, flagPossibleStereoCenters=True)
    return out


def normalize_acid_charge(mol: Mol) -> Mol:
    """同片段内质子化强酸与去质子化弱酸位共存时逐次搬质子；无改动返回原 mol。"""
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return mol
    out = mol
    for _ in range(mol.GetNumAtoms()):  # 上界：每次消耗一对 donor/acceptor
        frag_of = {i: fi for fi, tup in enumerate(Chem.GetMolFrags(out)) for i in tup}
        donors: list[int] = []
        acceptors: list[int] = []
        for i, a in enumerate(out.GetAtoms()):
            kind = _acid_kind_of_oh(a)
            if kind in DONOR_KIND:
                donors.append(i)
            elif _is_weak_anion(a):  # 弱受体不会同时是强酸供体
                acceptors.append(i)
        if not donors or not acceptors:
            break
        ranks = list(Chem.CanonicalRankAtoms(out))
        donors.sort(key=lambda i: (-DONOR_KIND.index(_acid_kind_of_oh(out.GetAtomWithIdx(i))), ranks[i]))  # 供体：DONOR_KIND 序靠前优先，再取 rank 最小保确定性
        d_idx = donors[0]
        dfrag = frag_of[d_idx]
        same_frag = [i for i in acceptors if frag_of[i] == dfrag]
        if not same_frag:
            break
        a_idx = min(same_frag, key=lambda i: ranks[i])
        nxt = _relocate_proton(out, a_idx, d_idx)
        if nxt is None:
            break
        if Chem.MolToSmiles(nxt) == Chem.MolToSmiles(out):
            break
        out = nxt
    return out
