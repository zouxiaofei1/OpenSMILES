"""构建诱导切割子分子，连接点用 H 封端以用于 free-name。"""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Mol


@dataclass(frozen=True)
class CutSubmol:
    """切割子分子及其新旧索引映射（连接点用 H 封端以用于 free-name）。"""
    mol: object  # RDKit Mol，连接点用 H 封端以用于 free-name
    atom_map: dict[int, int]  # new_idx -> old_idx（新索引→旧索引）
    inv_map: dict[int, int]  # old_idx -> new_idx（旧索引→新索引）
    attach_new: int  # 子分子中的连接原子索引
    attach_old: int
    atoms_old: frozenset[int]


def _ordered(atoms: frozenset[int]) -> list[int]:
    """返回排序后的原子索引列表。"""
    return sorted(atoms)


def _copy_atoms(em: Chem.RWMol, mol: Mol, order: list[int]) -> dict[int, int]:
    """复制原子到可写分子并记录新旧映射。"""
    inv: dict[int, int] = {}
    for old in order:
        inv[old] = em.AddAtom(mol.GetAtomWithIdx(old))
    return inv


def _copy_bonds(em: Chem.RWMol, mol: Mol, inv: dict[int, int]) -> list:
    """复制诱导子图内的键，返回被复制的源键列表（供 _carry_alkene_stereo 只扫子图内的键，不必遍历整分子）。"""
    copied: list = []
    for old_a, new_a in inv.items():
        for bond in mol.GetAtomWithIdx(old_a).GetBonds():
            old_b = bond.GetOtherAtomIdx(old_a)
            if old_b not in inv or old_b < old_a:
                continue
            em.AddBond(new_a, inv[old_b], bond.GetBondType())
            copied.append(bond)
    return copied


def _carry_alkene_stereo(em: Chem.RWMol, mol: Mol, inv: dict[int, int], bonds,
                         dummy: int | None = None) -> None:
    """迁移诱导子图内双键的 E/Z 立体到子分子：_copy_bonds 只重建键型会丢奇偶，故把源双键 E/Z 标签照搬（标签取原分子即真实立体，不随配基被切/被 * 顶替而重判）；两端引用被保留则映射到新索引，被切配基由本端 sp2 碳上的 dummy 顶替，无法唯一解析（如 H 封端成非手性 CH2）则跳过留无立体。`bonds` 为 _copy_bonds 返回的子图内源键。"""
    if em is None or mol is None:
        return
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
        bgn, end = nb.GetBeginAtomIdx(), nb.GetEndAtomIdx()  # RDKit SetStereoAtoms 要求两引用分别连在新键 begin/end 端；begin/end 由 AddBond 归一化（较小索引），E/Z 只取决于两引用是否同侧，故按 begin/end 换序传参。
        first, second = (na, nbb) if bgn == inv[ca] else (nbb, na)
        if em.GetBondBetweenAtoms(bgn, first) is None or em.GetBondBetweenAtoms(end, second) is None:
            continue
        nb.SetStereo(st)
        nb.SetStereoAtoms(first, second)


def _cap_attach_h(em: Chem.RWMol, attach_new: int) -> None:
    """让连接原子重新显式计算隐式 H 以允许 H 封端。"""
    atom = em.GetAtomWithIdx(attach_new)
    atom.SetNoImplicit(False)
    atom.UpdatePropertyCache(strict=False)


def _sanitize(em: Chem.RWMol) -> Mol | None:
    """对可写分子做 RDKit 消毒，失败返回 None。"""
    try:
        Chem.SanitizeMol(em)
    except Exception:
        return None
    return em.GetMol()


def _pack(out: Mol, inv: dict[int, int], attach_old: int, atoms: frozenset[int]) -> CutSubmol:
    """把结果封装为 CutSubmol 并补全新旧索引映射。"""
    atom_map = {n: o for o, n in inv.items()}
    return CutSubmol(out, atom_map, inv, inv[attach_old], attach_old, frozenset(atoms))




def _external_bond_type(mol: Mol, attach_old: int, atoms: frozenset[int]):
    """返回连接原子与母体(外)原子间的真实键型，无外部邻居(孤立自由基)回退单键；取代基叶若以双键连母体（外环 =CH2 等 *ylidene）须保留真实键级，否则 canonical 收成 *C 会命成饱和 alkyl（methyl 而非 methylidene，式量丢 H2）。"""
    for nb in mol.GetAtomWithIdx(attach_old).GetNeighbors():
        if nb.GetAtomicNum() != 1 and nb.GetIdx() not in atoms:
            b = mol.GetBondBetweenAtoms(attach_old, nb.GetIdx())
            if b is not None:
                return b.GetBondType()
    return Chem.BondType.SINGLE


def _add_anchor(em: Chem.RWMol, attach_new: int,
                bond_type: Chem.BondType = Chem.BondType.SINGLE) -> int:
    """在连接原子处添加 dummy 原子作锚点（键型随母体-取代基真实键级），返回其索引。"""
    d = em.AddAtom(Chem.Atom(0))  # 用 dummy 原子（`*`）标记连接原子；双键叶(=CH2)需用双键，键型由调用方给出。
    em.AddBond(attach_new, d, bond_type)
    return d


def build_anchor_submol(mol: Mol, atoms: frozenset[int], attach_old: int) -> Mol | None:
    """构建以 dummy 原子标记连接位点的诱导子分子。"""
    if attach_old not in atoms:
        return None
    em = Chem.RWMol()
    inv = _copy_atoms(em, mol, _ordered(atoms))
    copied = _copy_bonds(em, mol, inv)
    d = _add_anchor(em, inv[attach_old], _external_bond_type(mol, attach_old, atoms))
    _carry_alkene_stereo(em, mol, inv, copied, d)
    return _sanitize(em)  # print(Chem.MolToSmiles(em))

