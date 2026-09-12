"""L0 酸性质子重定位/电荷归一化"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from namepredict.constants import ACCEPTOR_Z, ACID_CENTERS, ACID_KIND_PRIO, DONOR_KIND, O


def _oxo_neighbor(atom, heavy: int, min_oxo: int = 1):
    """返回 atom 邻接的成酸中心原子（非芳香）若其带 ≥min_oxo 个双键氧（C(=O)/P(=O)/S(=O)n）；否则 None。"""
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() != heavy or n.GetIsAromatic():
            continue
        oxo = sum(
            1
            for b in n.GetBonds()
            if b.GetOtherAtom(n).GetAtomicNum() == O
            and b.GetBondType() == Chem.BondType.DOUBLE
        )
        if oxo >= min_oxo:
            return n
    return None

def _acid_kind(atom) -> str | None:
    """按 ACID_CENTERS 表返回 O 所连成酸中心对应的酸类名，不检查 O 自身 H/电荷。"""
    for z, (min_oxo, kind) in ACID_CENTERS.items():
        if _oxo_neighbor(atom, z, min_oxo):
            return kind
    return None

def _acid_kind_of_oh(atom) -> str | None:
    """判定中性含 H 的 O 是否为质子化强酸 OH：carboxyl 连 C(=O)、phospho 连 P(=O)、sulfo 连 S(=O)n；否则 None。"""
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
    """把质子从强酸供体 d_idx 搬到弱酸受体 a_idx：受体中性 +1H、供体 -1 -1H；只做 FormalCharge/H 记账，SanitizeMol 失败返回 None。"""
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


def _frag_of(mol: Mol) -> dict[int, int]:
    """返回 {atom_idx: 片段号} 映射。"""
    return {i: fi for fi, tup in enumerate(Chem.GetMolFrags(mol)) for i in tup}


def normalize_acid_charge(mol: Mol) -> Mol:
    # return mol
    """同一片段内【质子化羧酸】与【去质子化弱酸位】共存时逐次搬质子使负电荷收敛到最强酸；供体只取强酸 OH、受体只取弱酸阴离子故单调收敛，无改动返回原 mol。"""
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return mol
    out = mol
    for _ in range(mol.GetNumAtoms()):  # 上界：每次消耗一对 donor/acceptor
        frag_of = _frag_of(out)
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
        donors.sort(key=lambda i: (-ACID_KIND_PRIO.get(_acid_kind_of_oh(out.GetAtomWithIdx(i)), 0), ranks[i]))  # 供体：酸更强(更负 prio 取负序)优先，再取 canonical rank 最小保证确定性
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
