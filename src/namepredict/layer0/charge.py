"""L0 酸性质子重定位/电荷归一化"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from namepredict.constants import ACCEPTOR_Z, ACID_CENTERS, C, DONOR_KIND, N, O


def _acid_kind(atom) -> str | None:
    """按 ACID_CENTERS 表返回 O 所连成酸中心对应的酸类名（中心非芳香且带 ≥1 双键氧）。"""
    oxo_z: set[int] = set()
    for n in atom.GetNeighbors():
        if n.GetIsAromatic():
            continue
        if any(b.GetOtherAtom(n).GetAtomicNum() == O
               and b.GetBondType() == Chem.BondType.DOUBLE for b in n.GetBonds()):
            oxo_z.add(n.GetAtomicNum())
    return next((kind for z, kind in ACID_CENTERS.items() if z in oxo_z), None)

def _acid_kind_of_oh(atom) -> str | None:
    """判中性含 H 的 O 是否为质子化强酸 OH；否则 None。"""
    if atom.GetAtomicNum() != O or atom.GetFormalCharge() != 0 or atom.GetTotalNumHs() < 1:
        return None
    return _acid_kind(atom)

def _ring_system_of(mol, idx: int) -> set[int]:
    """返回 idx 所在稠环系（共享原子的环并入）的全部环原子索引。"""
    ri = mol.GetRingInfo()
    out: set[int] = set()
    for ring in ri.AtomRings():
        if idx in ring:
            out.update(ring)
    changed = True
    while changed:
        changed = False
        for ring in ri.AtomRings():
            if out & set(ring) and not set(ring) <= out:
                out.update(ring)
                changed = True
    return out


def _ring_bears_strong_acid(mol, ring_atoms: set[int]) -> bool:
    """环系邻域（环原子及其 2 键内）是否有中性强酸基（-COOH/-SO3H/-PO(OH)2）；
    有此酸时酚氧负离子保留作 -olate 母体（酸降为 carboxy/sulfo 前缀）。"""
    near = set(ring_atoms)
    for _ in range(2):  # 环原子向外扩两键：直连芳环的酸、苄基酸（-CH2COOH）皆覆盖
        near |= {n.GetIdx() for i in list(near) for n in mol.GetAtomWithIdx(i).GetNeighbors()}
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != O or _acid_kind_of_oh(a) is None:
            continue
        if a.GetNeighbors()[0].GetIdx() in near:  # 酸中心（C/S/P）落在环系邻域内
            return True
    return False


def _ring_has_nh(mol, ring_atoms: set[int]) -> bool:
    """环系内是否有带 H 的环氮（内酰胺/亚胺醇互变优先取酮式）。"""
    return any(a.GetAtomicNum() == N and a.GetTotalNumHs() > 0
               for a in mol.GetAtoms() if a.GetIdx() in ring_atoms)


def _is_weak_anion(atom) -> bool:
    """判定位点是否为去质子化的弱酸位（可接受质子）。"""
    if atom.GetFormalCharge() != -1 or atom.GetAtomicNum() not in ACCEPTOR_Z:
        return False
    if any(n.GetFormalCharge() == 1 for n in atom.GetNeighbors()):  # 邻接 +1 电荷 → 硝基/N-氧化物等内平衡写法，不当作弱酸位
        return False
    if atom.GetAtomicNum() == O and _acid_kind(atom) is not None:  # 已是强酸共轭碱（羧酸/磷酸/磺酸根），电荷位置合理
        return False
    if atom.GetAtomicNum() == O:  # 芳环氧负离子且环系邻域带强酸：保留 O⁻ 作 -olate 母体（不搬质子）
        mol = atom.GetOwningMol()
        arom = [n for n in atom.GetNeighbors() if n.GetIsAromatic()]
        if arom:
            ring = _ring_system_of(mol, arom[0].GetIdx())
            if _ring_bears_strong_acid(mol, ring) and not _ring_has_nh(mol, ring):
                return False  # 环氮带 H 者按内酰胺/酮式互变优先搬质子（金标 4-oxo-1H 式）
    return True

def _is_base_n(atom) -> bool:
    """判定中性脂肪胺氮（可接受质子的碱位，P-73.1.2）：有 N-H、非酰胺/亚胺、不连 O。"""
    if atom.GetAtomicNum() != N or atom.GetFormalCharge() != 0 or atom.GetIsAromatic():
        return False
    if atom.GetTotalNumHs() < 1:
        return False
    for b in atom.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE:  # =N- / [N+]= 非碱位
            return False
        nb = b.GetOtherAtom(atom)
        if nb.GetIsAromatic():  # 芳胺（苯胺/嘌呤氨基等）：碱性弱，仍按 amino 命名
            return False
        if nb.GetAtomicNum() == O:  # N-O 键（羟胺/N-氧化物）非碱位
            return False
        if nb.GetAtomicNum() == C and any(  # 酰胺/氨基甲酸酯氮：孤对已离域，不质子化
            b2.GetBondType() == Chem.BondType.DOUBLE
            and b2.GetOtherAtom(nb).GetAtomicNum() == O
            for b2 in nb.GetBonds()
        ):
            return False
    return True


def _relocate_proton(mol: Mol, a_idx: int, d_idx: int) -> Mol | None:
    """把质子从强酸供体搬到受体，只做电荷/H 记账；消毒失败返回 None。"""
    m = RWMol(mol)
    acc = m.GetAtomWithIdx(a_idx)
    don = m.GetAtomWithIdx(d_idx)
    acc.SetFormalCharge(acc.GetFormalCharge() + 1)  # 阴离子受体 -1→0；中性胺氮 0→+1
    acc.SetNumExplicitHs(acc.GetNumExplicitHs() + 1)
    acc.SetNoImplicit(True)
    don.SetFormalCharge(don.GetFormalCharge() - 1)
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
    anionic = Chem.GetFormalCharge(mol) < 0  # 净负离子才把中性胺氮当受体（保净电荷内的两性离子式，P-73.1.2）
    for _ in range(mol.GetNumAtoms()):  # 上界：每次消耗一对 donor/acceptor
        frag_of = {i: fi for fi, tup in enumerate(Chem.GetMolFrags(out)) for i in tup}
        donors: list[int] = []
        anions: list[int] = []
        bases: list[int] = []
        for i, a in enumerate(out.GetAtoms()):
            kind = _acid_kind_of_oh(a)
            if kind in DONOR_KIND:
                donors.append(i)
            elif _is_weak_anion(a):  # 已去质子化的弱酸位：受体首选
                anions.append(i)
            elif anionic and _is_base_n(a):  # 无弱酸位时才把中性胺氮当受体
                bases.append(i)
        acceptors = anions or bases
        if not donors or not acceptors:
            break
        ranks = list(Chem.CanonicalRankAtoms(out))
        kinds = {i: DONOR_KIND.index(_acid_kind_of_oh(out.GetAtomWithIdx(i))) for i in donors}  # 供体种类每原子只算一次
        donors.sort(key=lambda i: (-kinds[i], ranks[i]))  # 供体：DONOR_KIND 序靠前优先，再取 rank 最小保确定性
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
