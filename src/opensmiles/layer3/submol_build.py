"""构建诱导切割子分子，连接点用 H 封端以用于 free-name。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import Mol


def _copy_atoms(em: Chem.RWMol, mol: Mol, order: list[int]) -> dict[int, int]:
    """复制原子到可写分子并记录新旧映射。"""
    inv: dict[int, int] = {}
    for old in order:
        inv[old] = em.AddAtom(mol.GetAtomWithIdx(old))
    return inv


def _copy_bonds(em: Chem.RWMol, mol: Mol, inv: dict[int, int]) -> list:
    """复制诱导子图内的键，返回源键列表供立体迁移扫描。"""
    copied: list = []
    for old_a, new_a in inv.items():
        for bond in mol.GetAtomWithIdx(old_a).GetBonds():
            old_b = bond.GetOtherAtomIdx(old_a)
            if old_b not in inv or old_b < old_a:
                continue
            em.AddBond(new_a, inv[old_b], bond.GetBondType())
            copied.append(bond)
    return copied


def _carry_alkene_stereo(em: Chem.RWMol, inv: dict[int, int], bonds,
                         dummy: int | None = None) -> None:
    """把子图内双键的 E/Z 标签照搬到子分子；无法解析则跳过。"""
    for b in bonds:
        if b.GetBondType() is not Chem.BondType.DOUBLE:
            continue
        st = b.GetStereo()
        if st not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            continue
        ca, cb = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if ca not in inv or cb not in inv:
            continue
        refs = b.GetStereoAtoms()
        if len(refs) < 2:
            continue
        ref_a, ref_b = refs

        def _resolve(carbon: int, ref: int) -> int | None:
            """把源端引用映射到子分子：保内则平移，切掉则仅当 dummy 连在本端时以其顶替。"""
            if ref in inv:
                return inv[ref]
            if dummy is not None and em.GetBondBetweenAtoms(dummy, inv[carbon]) is not None:
                return dummy
            return None

        na = _resolve(ca, ref_a)
        nbb = _resolve(cb, ref_b)
        if na is None or nbb is None:
            continue
        nb = em.GetBondBetweenAtoms(inv[ca], inv[cb])
        if nb is None:
            continue
        bgn, end = nb.GetBeginAtomIdx(), nb.GetEndAtomIdx()  # SetStereoAtoms 要求两引用分连新键 begin/end 端
        first, second = (na, nbb) if bgn == inv[ca] else (nbb, na)
        if em.GetBondBetweenAtoms(bgn, first) is None or em.GetBondBetweenAtoms(end, second) is None:
            continue
        nb.SetStereo(st)
        nb.SetStereoAtoms(first, second)


def _sanitize(em: Chem.RWMol) -> Mol | None:
    """对可写分子做 RDKit 消毒，失败返回 None。"""
    try:
        Chem.SanitizeMol(em)
    except Exception:
        return None
    return em.GetMol()


def _external_bond_type(mol: Mol, attach_old: int, atoms: frozenset[int]):
    """返回连接原子与母体间真实键型，无外邻居时回退单键。"""
    for nb in mol.GetAtomWithIdx(attach_old).GetNeighbors():
        if nb.GetAtomicNum() != 1 and nb.GetIdx() not in atoms:
            b = mol.GetBondBetweenAtoms(attach_old, nb.GetIdx())
            if b is not None:
                return b.GetBondType()
    return Chem.BondType.SINGLE


def _add_anchor(em: Chem.RWMol, attach_new: int,
                bond_type: Chem.BondType = Chem.BondType.SINGLE) -> int:
    """在连接原子处添加 dummy 锚点，返回其索引。"""
    d = em.AddAtom(Chem.Atom(0))  # 用 dummy 原子（`*`）标记连接原子；键型由调用方给出
    em.AddBond(attach_new, d, bond_type)
    return d


def build_anchor_submol(mol: Mol, atoms: frozenset[int], attach_old: int) -> Mol | None:
    """构建以 dummy 原子标记连接位点的诱导子分子。"""
    if attach_old not in atoms:
        return None
    em = Chem.RWMol()
    inv = _copy_atoms(em, mol, sorted(atoms))
    copied = _copy_bonds(em, mol, inv)
    d = _add_anchor(em, inv[attach_old], _external_bond_type(mol, attach_old, atoms))
    _carry_alkene_stereo(em, inv, copied, d)
    return _sanitize(em)