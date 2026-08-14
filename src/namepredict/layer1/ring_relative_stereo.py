"""简单饱和环上相对面的 L1 带类型图事实。"""
from __future__ import annotations

from dataclasses import dataclass
from rdkit.Chem import ChiralType, Mol

from namepredict.constants import H


@dataclass(frozen=True)
class RingRelativeStereoIR:
    """以环原子为键的面朝向符号；``None`` 表示无显式立体化学。"""
    faces: tuple[tuple[int, int], ...] = ()
    invalid: bool = False


def _parity(source: list[int], target: list[int]) -> int:
    positions = {value: i for i, value in enumerate(source)}
    order = [positions[value] for value in target]
    return -1 if sum(order[i] > order[j] for i in range(len(order)) for j in range(i + 1, len(order))) % 2 else 1


def _orientation(ring_order: list[int]) -> int:
    start = ring_order.index(min(ring_order))
    rotated = ring_order[start:] + ring_order[:start]
    return 1 if rotated[1] < rotated[-1] else -1


def _face(mol: Mol, atom_id: int, ring_order: list[int], ligand: int) -> int | None:
    atom = mol.GetAtomWithIdx(atom_id); tag = atom.GetChiralTag()
    if tag not in (ChiralType.CHI_TETRAHEDRAL_CW, ChiralType.CHI_TETRAHEDRAL_CCW): return None
    index = ring_order.index(atom_id); previous, following = ring_order[index - 1], ring_order[(index + 1) % len(ring_order)]
    actual = [-1 if n.GetAtomicNum() == H else n.GetIdx() for n in atom.GetNeighbors()]
    if atom.GetNumExplicitHs() + atom.GetNumImplicitHs() and -1 not in actual: actual = [-1] + actual
    target = [ligand, previous, following, -1]
    if set(actual) != set(target): return None
    handed = 1 if tag == ChiralType.CHI_TETRAHEDRAL_CW else -1
    return handed * _parity(actual, target) * _orientation(ring_order)


def _relative_faces(faces: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    ordered = sorted(faces)
    reference = ordered[0][1] if ordered else 1
    return tuple((atom, face * reference) for atom, face in ordered)


def ring_relative_stereo(mol: Mol, ring_order: list[int], ligands: dict[int, int]) -> RingRelativeStereoIR:
    """推导由图定义的环面朝向符号，无需 CIP 标记或构象。"""
    tagged = [a.GetIdx() for a in mol.GetAtoms() if a.GetChiralTag() != ChiralType.CHI_UNSPECIFIED]
    if any(atom not in ligands for atom in tagged) or len(ligands) > 3: return RingRelativeStereoIR(invalid=bool(tagged))
    faces = []
    for atom, ligand in ligands.items():
        face = _face(mol, atom, ring_order, ligand)
        if face is None and tagged: return RingRelativeStereoIR(invalid=True)
        if face is not None: faces.append((atom, face))
    if tagged and len(faces) != len(ligands): return RingRelativeStereoIR(invalid=True)
    return RingRelativeStereoIR(_relative_faces(faces))
